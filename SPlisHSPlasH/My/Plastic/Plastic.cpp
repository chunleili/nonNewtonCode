#include "Plastic.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"
#include "../Utils/myPrint.h"

using namespace SPH;
using namespace GenParam;

int Plastic::ELASTIC_LIMIT = -1;
int Plastic::PLASTIC_LIMIT = -1;
int Plastic::IS_PLASTIC = -1;


Plastic::Plastic(FluidModel *model) :
	ElasticityBase(model)
{
	numParticles = model->numActiveParticles();
	m_restVolumes.resize(numParticles);
	m_current_to_initial_index.resize(numParticles);
	m_initial_to_current_index.resize(numParticles);
	m_initialNeighbors.resize(numParticles);
	m_rotations.resize(numParticles, Matrix3r::Identity());
	m_stress.resize(numParticles);
	m_F.resize(numParticles);

	elasticLimit = 0.001;
	plasticLimit = 0.486;

	m_plasticStrain.resize(numParticles);
	m_isFracture.resize(numParticles,0);
	m_totalStrain.resize(numParticles);
	m_elasticStrain.resize(numParticles);

	initValues();

	model->addField({ "rest volume", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_restVolumes[i]; }, true });
	model->addField({ "rotation", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_rotations[i](0,0); } });
	model->addField({ "stress", FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_stress[i][0]; } });
	model->addField({ "deformation gradient", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_F[i](0,0); } });

	model->addField({ "plastic strain", FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_plasticStrain[i][0]; }, true });
	model->addField({ "isFracture", FieldType::UInt, [&](const unsigned int i) -> int* { return &m_isFracture[i]; }, true  });
}

Plastic::~Plastic(void)
{
	m_model->removeFieldByName("rest volume");
	m_model->removeFieldByName("rotation");
	m_model->removeFieldByName("stress");
	m_model->removeFieldByName("deformation gradient");

	m_model->removeFieldByName("plastic strain");
	m_model->removeFieldByName("isFracture");

}

void Plastic::initParameters()
{
	ElasticityBase::initParameters();

	//for plasticity
	ELASTIC_LIMIT = createNumericParameter("elasticLimit", "elastic limit", &elasticLimit);
	setGroup(ELASTIC_LIMIT, "Elasticity");
	setDescription(ELASTIC_LIMIT, "elastic limit");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(ELASTIC_LIMIT));
	rparam->setMinValue(0.0);

	PLASTIC_LIMIT = createNumericParameter("plasticLimit", "plastic limit", &plasticLimit);
	setGroup(PLASTIC_LIMIT, "Elasticity");
	setDescription(PLASTIC_LIMIT, "plastic limit");
	rparam = static_cast<RealParameter*>(getParameter(PLASTIC_LIMIT));
	rparam->setMinValue(0.0);

	IS_PLASTIC = createBoolParameter("isPlastic","isPlastic",&isPlastic);
	setGroup(IS_PLASTIC, "Elasticity");
}

void Plastic::initValues()
{
	numParticles = m_model->numActiveParticles();
	Simulation *sim = Simulation::getCurrent();
	sim->getNeighborhoodSearch()->find_neighbors();

	FluidModel *model = m_model;
	const unsigned int fluidModelIndex = model->getPointSetIndex();

	// Store the neighbors in the reference configurations and
	// compute the volume of each particle in rest state
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_current_to_initial_index[i] = i;
			m_initial_to_current_index[i] = i;

			// only neighbors in same phase will influence elasticity
			const unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);
			m_initialNeighbors[i].resize(numNeighbors);
			for (unsigned int j = 0; j < numNeighbors; j++)
				m_initialNeighbors[i][j] = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);

			// compute volume
			Real density = model->getMass(i) * sim->W_zero();
			const Vector3r &xi = model->getPosition(i);
			forall_fluid_neighbors_in_same_phase(
				density += model->getMass(neighborIndex) * sim->W(xi - xj);
			)
			m_restVolumes[i] = model->getMass(i) / density;

			//setzero for plastic strain
			m_plasticStrain[i].setZero();
			m_elasticStrain[i].setZero();
		}
	}
}

void Plastic::step()
{
	computeRotations();
	computeStress();
	computeForces();
	fracture();
	m_step++;
	// std::string fname = "m_initial_to_current_index_" + std::to_string(m_step)+ ".txt";
	// printScalarField(fname,m_initial_to_current_index,0);
	// std::string fname = "m_initialNeighbors_" + std::to_string(m_step)+ ".txt";
	// printVectorField(fname,m_initialNeighbors,0);
}


void Plastic::reset()
{
	initValues();
}

void Plastic::performNeighborhoodSearchSort()
{
	if (numParticles == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_restVolumes[0]);
	d.sort_field(&m_current_to_initial_index[0]);
	d.sort_field(&m_plasticStrain[0]);
	d.sort_field(&m_isFracture[0]);
	d.sort_field(&m_elasticStrain[0]);

	for (unsigned int i = 0; i < numParticles; i++)
		m_initial_to_current_index[m_current_to_initial_index[i]] = i;
}

void Plastic::computeRotations()
{
	Simulation *sim = Simulation::getCurrent();

// #pragma omp parallel default(shared)
	{
// #pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_model->getParticleState(i) == ParticleState::Active)
			{

				const unsigned int i0 = m_current_to_initial_index[i];
				const Vector3r &xi = m_model->getPosition(i);
				const Vector3r &xi0 = m_model->getPosition0(i0);
				Matrix3r Apq;
				Apq.setZero();

				const size_t numNeighbors = m_initialNeighbors[i0].size();

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < numNeighbors; j++)
				{
					const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
					// get initial neighbor index considering the current particle order
					const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

					const Vector3r &xj = m_model->getPosition(neighborIndex);
					const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);
					const Vector3r xj_xi = xj - xi;
					const Vector3r xj_xi_0 = xj0 - xi0;
					Apq += m_model->getMass(neighborIndex) * sim->W(xj_xi_0) * (xj_xi * xj_xi_0.transpose());
				}

				Quaternionr q(m_rotations[i]);
				MathFunctions::extractRotation(Apq, q, 10);
				m_rotations[i] = q.matrix();
			}
		}
	}
}

void Plastic::computeStress()
{
	// Elasticity tensor
	Matrix6r C;
	computeC(C);

	// #pragma omp parallel default(shared)
	{
		// #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_model->getParticleState(i) == ParticleState::Active)
			{
				Matrix3r nablaU;
				computeNablaU(i, nablaU);

				m_F[i] = nablaU + Matrix3r::Identity();
				
				Vector6r totalStrain;
				computeTotalStrain(nablaU, totalStrain);
				m_totalStrain[i] = totalStrain;

				if(isPlastic)
					computePlasticStrain(i, m_elasticStrain[i]); //FIXME: with bug

				m_elasticStrain[i] = m_totalStrain[i] - m_plasticStrain[i];
				m_stress[i] = C * m_elasticStrain[i];
			}
			else
				m_stress[i].setZero();
		}
	}
}


void Plastic::computeC(Matrix6r& C)
{
	C.setZero();
	const Real factor = m_youngsModulus / ((static_cast<Real>(1.0) + m_poissonRatio)*(static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));
	C(0, 0) = C(1, 1) = C(2, 2) = factor * (static_cast<Real>(1.0) - m_poissonRatio);
	C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = factor * (m_poissonRatio);
	C(3, 3) = C(4, 4) = C(5, 5) = factor * static_cast<Real>(0.5)*(static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio);
}


/**
 * @brief compute NablaU
 * 
 * @param i: input, particle index
 * @param nablaU: output, gradient of displacement
 */
void Plastic::computeNablaU(int i, Matrix3r &nablaU)
{
	Simulation *sim = Simulation::getCurrent();

	const unsigned int i0 = m_current_to_initial_index[i];
	const Vector3r& xi = m_model->getPosition(i);
	const Vector3r& xi0 = m_model->getPosition0(i0);

	nablaU.setZero();
	const size_t numNeighbors = m_initialNeighbors[i0].size();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int j = 0; j < numNeighbors; j++)
	{
		const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
		// get initial neighbor index considering the current particle order 
		const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

		const Vector3r& xj = m_model->getPosition(neighborIndex);
		const Vector3r& xj0 = m_model->getPosition0(neighborIndex0);

		const Vector3r xj_xi = xj - xi;
		const Vector3r xj_xi_0 = xj0 - xi0;

		const Vector3r uji = m_rotations[i].transpose() * xj_xi - xj_xi_0;

		// subtract because kernel gradient is taken in direction of xji0 instead of xij0
		// Eq(6) in Becker2009
		nablaU -= (m_restVolumes[neighborIndex] * uji) * sim->gradW(xj_xi_0).transpose();
	}	
}

/**
 * @brief compute totalStrain
 * 
 * @param nablaU: input, gradient of displacement
 * @param totalStrain: output
 */
void Plastic::computeTotalStrain(Matrix3r &nablaU, Vector6r & totalStrain)
{
	totalStrain.setZero();
	// compute Cauchy strain: epsilon = 0.5 (nabla u + nabla u^T)
	totalStrain[0] = nablaU(0, 0);						// \epsilon_{00}
	totalStrain[1] = nablaU(1, 1);						// \epsilon_{11}
	totalStrain[2] = nablaU(2, 2);						// \epsilon_{22}
	totalStrain[3] = static_cast<Real>(0.5) * (nablaU(0, 1) + nablaU(1, 0)); // \epsilon_{01}
	totalStrain[4] = static_cast<Real>(0.5) * (nablaU(0, 2) + nablaU(2, 0)); // \epsilon_{02}
	totalStrain[5] = static_cast<Real>(0.5) * (nablaU(1, 2) + nablaU(2, 1)); // \epsilon_{12}
}

/**
 * @brief Compute the plastic strain based on O'Brien 2002. 
 * 
 * Ref: James F. O’Brien et. al. 2002, "Graphical Modeling and Animation of Ductile Fracture"
 */
void Plastic::computePlasticStrain(int i, Vector6r &strain)
{
	//Eq(2) in O'Brien 2002
	Vector6r deviation =  strain;
	Real trace = strain[0] + strain[1] + strain[2];
	trace /= 3.0;
	deviation[0] -= trace;
	deviation[1] -= trace;
	deviation[2] -= trace;

	//Eq(3),  Consider plasticity only if exceeding elasticLimit
	Real deviationNorm = FNorm(deviation);

	if(deviationNorm <= elasticLimit)
		return;
	
	//Eq(4)
	Vector6r strainIncrement  = ( deviationNorm - elasticLimit ) / deviationNorm * deviation;

	//Eq(5) Calculate the accumulated plastic strain
	Vector6r newPlastic = m_plasticStrain[i] + strainIncrement;
	Real plasticNorm = FNorm(newPlastic);
	Real ratio = (plasticLimit / plasticNorm);
	Real min = 1.0 < ratio ? 1.0 : ratio;
	m_plasticStrain[i] = newPlastic * min;

	if (plasticNorm > plasticLimit)
		m_isFracture[i] = 1;
	else
		m_isFracture[i] = 0;
}



void Plastic::computeForces()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_isFracture[i] == 1)
				continue; //jump over the fractured particles

			if (m_model->getParticleState(i) == ParticleState::Active)
			{
				const unsigned int i0 = m_current_to_initial_index[i];
				const Vector3r& xi0 = m_model->getPosition0(i0);

				const size_t numNeighbors = m_initialNeighbors[i0].size();
				Vector3r fi;
				fi.setZero();

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < numNeighbors; j++)
				{
					const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
					// get initial neighbor index considering the current particle order 
					const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

					const Vector3r& xj0 = m_model->getPosition0(neighborIndex0);

					const Vector3r xj_xi_0 = xj0 - xi0;
					const Vector3r gradW0 = sim->gradW(xj_xi_0);

					const Vector3r dji = m_restVolumes[i] * gradW0;
					const Vector3r dij = -m_restVolumes[neighborIndex] * gradW0;

					Vector3r sdji, sdij;
					symMatTimesVec(m_stress[neighborIndex], dji, sdji); //stress times dji
					symMatTimesVec(m_stress[i], dij, sdij);

					const Vector3r fij = -m_restVolumes[neighborIndex] * sdji;
					const Vector3r fji = -m_restVolumes[i] * sdij;

					fi += m_rotations[neighborIndex] * fij - m_rotations[i] * fji;
				}
				fi = 0.5*fi;

				// elastic acceleration
				Vector3r& ai = m_model->getAcceleration(i);
				ai += fi / m_model->getMass(i);
			}
		}
	}
}

void Plastic::fracture()
{}