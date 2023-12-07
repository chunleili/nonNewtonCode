#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "NonNewton.h"
#include "Utils/myPrint.h"
#include "Utils/mathUtils.h"
#include "extern/my_partio/Partio.h"
#include "extern/my_partio/PartioSingleton.h"

using namespace SPH;
using namespace GenParam;

NonNewton::NonNewton(FluidModel *model) :
SurfaceTensionBase(model)
{
	std::cout<<"constructor\n";

	numParticles = model->numActiveParticles();

	m_strainRateNorm.resize(numParticles, 0.0);
	m_nonNewtonViscosity.resize(numParticles, 0.0);
	m_strainRate.resize(numParticles, Vector6r::Zero());
	m_boundaryViscosity.resize(numParticles, 0.0);

	model->addField({ "strainRate", FieldType::Vector6, [&](const unsigned int i) -> Vector6r* { return &m_strainRate[i]; }, true });
	//nonNewtonViscosity already exists in FluidModel
	model->addField({ "strainRateNorm", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_strainRateNorm[i]; }, true });

}

NonNewton::~NonNewton(void)
{
	m_model->removeFieldByName("strainRate");
	m_strainRate.clear();

	m_model->removeFieldByName("strainRateNorm");
	m_strainRateNorm.clear();
}

void NonNewton::init()
{
	std::cout<<"init\n";
	initParameters();
}

void NonNewton::initParameters()
{
	MAX_VISCOSITY = createNumericParameter("max_viscosity", "max_viscosity", &m_maxViscosity);
	AVG_VISCOSITY = createNumericParameter("average_viscosity", "average_viscosity", &m_avgViscosity);
	MIN_VISCOSITY = createNumericParameter("min_viscosity", "min_viscosity", &m_minViscosity);

	setGroup(MAX_VISCOSITY, "Viscosity");
	setGroup(AVG_VISCOSITY, "Viscosity");
	setGroup(MIN_VISCOSITY, "Viscosity");

	setDescription(MAX_VISCOSITY, "Max viscosity of all fluid particles.");
	setDescription(AVG_VISCOSITY, "Average viscosity of all fluid particles.");
	setDescription(MIN_VISCOSITY, "Min viscosity of all fluid particles.");

	getParameter(MAX_VISCOSITY)->setReadOnly(true);
	getParameter(AVG_VISCOSITY)->setReadOnly(true);
	getParameter(MIN_VISCOSITY)->setReadOnly(true);

	NON_NEWTON_METHOD = createEnumParameter("nonNewtonMethod", "nonNewtonMethod", &(int)m_nonNewtonMethod);
	setGroup(NON_NEWTON_METHOD, "Viscosity");
	setDescription(NON_NEWTON_METHOD, "Method for nonNewton.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(NON_NEWTON_METHOD));
	enumParam->addEnumValue("Newtonian", NEWTONIAN_);
	enumParam->addEnumValue("Power Law", POWER_LAW_);
	enumParam->addEnumValue("Cross Model", CROSS_);
	enumParam->addEnumValue("Casson Model", CASSON_);
	enumParam->addEnumValue("Carreau Model", CARREAU_);
	enumParam->addEnumValue("Bingham Model", BINGHAM_);
	enumParam->addEnumValue("Herschel-Bulkley Model", HERSCHEL_BULKLEY_);


	POWER_INDEX = createNumericParameter("power_index", "power_index", &power_index);
	CONSISTENCY_INDEX = createNumericParameter("consistency_index", "consistency_index", &consistency_index);
	VISCOSITY0 = createNumericParameter("viscosity0", "viscosity0", &m_viscosity0);
	VISCOSITY_INF = createNumericParameter("viscosity_inf", "viscosity_inf", &m_viscosity_inf);
	MU_C = createNumericParameter("muC", "muC", &m_muC);
	CRITICAL_STRAIN_RATE = createNumericParameter("criticalStrainRate", "criticalStrainRate", &m_criticalStrainRate);

	setGroup(POWER_INDEX, "Viscosity");
	setGroup(CONSISTENCY_INDEX, "Viscosity");
	setGroup(VISCOSITY0, "Viscosity");
	setGroup(VISCOSITY_INF, "Viscosity");
	setGroup(MU_C, "Viscosity");
	setGroup(CRITICAL_STRAIN_RATE, "Viscosity");

	setDescription(POWER_INDEX, "Power index for power law.");
	setDescription(CONSISTENCY_INDEX, "Consistency index for power law.");
	setDescription(VISCOSITY0, "Initial viscosity.");
	setDescription(VISCOSITY_INF, "Infinite viscosity for the cross model.");
	setDescription(MU_C, "Critical shear rate for the Casson model.");
	setDescription(CRITICAL_STRAIN_RATE, "Critical strain rate for the Herschel-Bulkley model.");


    // SMOOTH_VELOCITY_FLAG = createBoolParameter("smoothVelocityFlag", "smoothVelocityFlag", &m_smoothVelocityFlag);
    // setGroup(SMOOTH_VELOCITY_FLAG, "Damping");
    // setDescription(SMOOTH_VELOCITY_FLAG, "turn on smoothVelocity.");

	// SMOOTH_VELOCITY_FACTOR = createNumericParameter("smoothVelocityFactor", "smoothVelocityFactor", &m_smoothVelocityFactor);
	// setGroup(SMOOTH_VELOCITY_FACTOR, "Damping");
	// setDescription(SMOOTH_VELOCITY_FACTOR, "smoothVelocityFactor(0 to 1). 0 means no smoothing, 1 means use complete average vel, i.e, sum(mj*(vj)/rhoj*Wij).");
	// static_cast<RealParameter*>(getParameter(SMOOTH_VELOCITY_FACTOR))->setMinValue(0.0);
	// static_cast<RealParameter*>(getParameter(SMOOTH_VELOCITY_FACTOR))->setMaxValue(1.0);

    // DAMP_VELOCITY_FLAG = createBoolParameter("dampVelocityFlag", "dampVelocityFlag", &m_dampVelocityFlag);
    // setGroup(DAMP_VELOCITY_FLAG, "Viscosity");
    // setDescription(DAMP_VELOCITY_FLAG, "turn on dampVelocity.");

	// DAMP_VELOCITY_FACTOR = createNumericParameter("dampVelocityFactor", "dampVelocityFactor", &m_dampVelocityFactor);
	// setGroup(DAMP_VELOCITY_FACTOR, "Viscosity");
	// setDescription(DAMP_VELOCITY_FACTOR, "dampVelocityFactor(0 to 1). 0 means no damping, 1 means use complete last step vel.");
	// static_cast<RealParameter*>(getParameter(DAMP_VELOCITY_FACTOR))->setMinValue(0.0);
	// static_cast<RealParameter*>(getParameter(DAMP_VELOCITY_FACTOR))->setMaxValue(1.0);

	CHANGE_WITH_TIME_FLAG = createBoolParameter("changeWithTimeFlag", "changeWithTimeFlag", &m_changeWithTimeFlag);
	VISCOSITY_WITH_TIME_PART = createNumericParameter("viscosityWithTimePart", "viscosityWithTimePart", &m_viscosityWithTimePart);
	COEFF_A = createNumericParameter("coeffA", "coeffA", &m_coeffA);
	COEFF_B = createNumericParameter("coeffB", "coeffB", &m_coeffB);
	COEFF_C = createNumericParameter("coeffC", "coeffC", &m_coeffC);
	COEFF_D = createNumericParameter("coeffD", "coeffD", &m_coeffD);
	setGroup(CHANGE_WITH_TIME_FLAG, "Viscosity");
	setGroup(VISCOSITY_WITH_TIME_PART, "Viscosity");
	setGroup(COEFF_A, "Viscosity");
	setGroup(COEFF_B, "Viscosity");
	setGroup(COEFF_C, "Viscosity");
	setGroup(COEFF_D, "Viscosity");
	setDescription(CHANGE_WITH_TIME_FLAG, "turn on changeWithTime.");
	setDescription(VISCOSITY_WITH_TIME_PART, "viscosityWithTimePart: viscosity = viscosityOriginal + viscosityWithTimePart.");
	setDescription(COEFF_A, "coeffA for viscosity change with time: viscosityWithTimePart = a * t^3 + b * t^2 + c * t + d");
	setDescription(COEFF_B, "coeffB for viscosity change with time: viscosityWithTimePart = a * t^3 + b * t^2 + c * t + d");
	setDescription(COEFF_C, "coeffC for viscosity change with time: viscosityWithTimePart = a * t^3 + b * t^2 + c * t + d");
	setDescription(COEFF_D, "coeffD for viscosity change with time: viscosityWithTimePart = a * t^3 + b * t^2 + c * t + d");



	CONTROLLED_BY_TEMPERATURE_FLAG = createBoolParameter("controlledByTemperature", "controlledByTemperature", &m_controlledByTemperatureFlag);
	setGroup(CONTROLLED_BY_TEMPERATURE_FLAG, "Viscosity");

	SMOOTH_STRAIN_RATE_FACTOR = createNumericParameter("smoothStrainRateFactor", "smoothStrainRateFactor", &m_smoothStrainRateFactor);
	setGroup(SMOOTH_STRAIN_RATE_FACTOR, "Viscosity");

}



void NonNewton::reset()
{
	std::cout<<"reset!\n";

	numParticles = m_model->numActiveParticles();
	
	m_strainRate.resize(numParticles, Vector6r::Zero());
	m_nonNewtonViscosity.resize(numParticles, 0.0);
	m_boundaryViscosity.resize(numParticles, 0.0);
}

// NOTE: SmoothVelocity and DampVelocity are moved to MyTimeStep.cpp
// void NonNewton::dampVelocity()
// {
// 	// moved to MyTimeStep.cpp

// 	// Simulation *sim = Simulation::getCurrent();
// 	// const unsigned int fluidModelIndex = m_model->getPointSetIndex();

// 	// for (unsigned int i = 0; i < numParticles; ++i)
// 	// {
// 	// 	const Vector3r v_last = m_model->getLastVelocity(i);
// 	// 	const Vector3r &vi = m_model->getVelocity(i);
// 	// 	const Vector3r vnew = vi - m_dampVelocityFactor * (vi - v_last);
// 	// 	if (i==0)
// 	// 	{
// 	// 		std::cout<<"vdiff"<<vi - v_last<<std::endl;
// 	// 	}
// 	// 	m_model->setVelocity(i, vnew);
// 	// }
// }


// void NonNewton::smoothVelocity()
// {
// 	Simulation *sim = Simulation::getCurrent();
// 	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

// 	for (unsigned int i = 0; i < numParticles; ++i)
// 	{
// 		const Vector3r &xi = m_model->getPosition(i);
// 		const Vector3r &vi = m_model->getVelocity(i);
// 		const Real density_i = m_model->getDensity(i);
// 		Vector3r avgVel = Vector3r::Zero(); 
// 		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
// 		{
// 			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
// 			const Vector3r &xj = m_model->getPosition(neighborIndex);
// 			const Vector3r &vj = m_model->getVelocity(neighborIndex);
// 			const Real W = sim->W(xi - xj);
// 			// const Vector3r vji = vj - vi;
// 			const Real mj = m_model->getMass(neighborIndex);
// 			const Real rhoj = m_model->getDensity(neighborIndex);
// 			avgVel += vj * mj/rhoj * W;
// 		}
// 		const Vector3r diff = vi - avgVel;
// 		const Vector3r newV = vi - diff * m_smoothVelocityFactor;
// 		m_model->setVelocity(i, newV);
// 	}
// }

//  a * t^3 + b * t^2 + c * t + d
void NonNewton::viscosityChangeWithTime()
{
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	Real t = TimeManager::getCurrent()->getTime();
	m_viscosityWithTimePart = m_coeffA * pow(t, 3) + m_coeffB * pow(t, 2) + m_coeffC * t + m_coeffD;
}


void NonNewton::calcStrainRate()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	for (unsigned int i = 0; i < numParticles; ++i)
	{
		const Vector3r &xi = m_model->getPosition(i);
		const Vector3r &vi = m_model->getVelocity(i);
		const Real density_i = m_model->getDensity(i);

		Matrix3r velGrad = Matrix3r::Zero();
		Matrix3r strainRate = Matrix3r::Zero();

		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
		{
			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
			const Vector3r &xj = m_model->getPosition(neighborIndex);
			const Vector3r &vj = m_model->getVelocity(neighborIndex);
			const Vector3r gradW = sim->gradW(xi - xj);
			const Vector3r vji = vj - vi;

			const Real m = m_model->getMass(neighborIndex);
			velGrad = vji * gradW.transpose();
			strainRate = velGrad + velGrad.transpose();
			strainRate *= m;
		}
		strainRate = (static_cast<Real>(0.5) / density_i) * strainRate;

		Real norm = 0.0f;
		norm = strainRate.norm();
		m_strainRateNorm[i] = norm;
	}
}

Real NonNewton::FNorm(const Vector6r & vec) const
{
	Real res = 0.0f;
	for (int i = 0; i < 6; i++)
		res += vec[i] * vec[i];
	res += vec[3] * vec[3];
	res += vec[4] * vec[4];
	res += vec[5] * vec[5];
	res = sqrt(res)	;
	return res;
}





void NonNewton::step()
{
	static int steps{0};
	steps++;
	// std::cout<<"\nstep: "<<steps<<"\n";
	numParticles = m_model->numActiveParticles();

	if(!m_controlledByTemperatureFlag)
		computeNonNewtonViscosity();

	// smoothVelocity和dampVelocity都已经搬到MyTimeStep中
	// if (m_smoothVelocityFlag)
	// 	smoothVelocity();

	// if (m_dampVelocityFlag)
	// 	dampVelocity();
	
	if (m_changeWithTimeFlag)
		viscosityChangeWithTime();

	if (m_controlledByTemperatureFlag)
		controlledByTemperature();
	
	// std::cout<<"m_viscosityWithTimePart: "<<m_viscosityWithTimePart<<"\n";
	//通过set函数将计算得到的粘度传递给fluidmodel
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] +=  m_viscosityWithTimePart;
		m_model->setNonNewtonViscosity(i,m_nonNewtonViscosity[i]);
	}

	m_maxViscosity = maxField(m_nonNewtonViscosity, numParticles);
	m_minViscosity = minField(m_nonNewtonViscosity, numParticles);
	m_avgViscosity = avgField(m_nonNewtonViscosity, numParticles);
	// [m_minViscosity, m_maxViscosity] = minMaxField(m_nonNewtonViscosity, numParticles);
	// [m_minViscosity, m_maxViscosity, m_avgViscosity] = minMaxAvgField(m_boundaryViscosity, numParticles);

	// if (m_smoothVelocityFlag)
	// 	smoothVelocity();

	// if (m_dampVelocityFlag)
	// 	dampVelocity();
}

void NonNewton::computeNonNewtonViscosity()
{
	calcStrainRate();

	if(m_nonNewtonMethod == NonNewtonMethod::ENUM_POWER_LAW)
	{
		computeViscosityPowerLaw();
	}
	else if(m_nonNewtonMethod == NonNewtonMethod::ENUM_CROSS_MODEL)
	{
		computeViscosityCrossModel();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_CASSON_MODEL)
	{
		computeViscosityCassonModel();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_CARREAU)
	{
		computeViscosityCarreau();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_BINGHAM)
	{
		computeViscosityBingham();
	}
	else if (m_nonNewtonMethod == NonNewtonMethod::ENUM_HERSCHEL_BULKLEY)
	{
		computeViscosityHerschelBulkley();
	}
	else
	{
		computeViscosityNewtonian();
	}

}



void NonNewton::computeViscosityNewtonian() 
{
	std::cout<<"computeViscosityNewtonian!\n";
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = m_viscosity0;
	}
}

void NonNewton::computeViscosityPowerLaw() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = consistency_index * pow(m_strainRateNorm[i], power_index - 1);
	}
}

void NonNewton::computeViscosityCrossModel() 
{
	assert((m_viscosity0 - m_viscosity_inf >= 0.0) && "the viscosity0 must be larger than viscosity_inf");
	if(m_viscosity0 - m_viscosity_inf < 0.0)
	{
		std::cout << "the viscosity0 must be larger than viscosity_inf" << std::endl;
		throw std::runtime_error("the viscosity0 must be larger than viscosity_inf");
	}
	for (unsigned int i = 0; i < numParticles; ++i)
	{

		m_nonNewtonViscosity[i] = m_viscosity_inf +  (m_viscosity0 - m_viscosity_inf) / (1 +  pow(consistency_index * m_strainRateNorm[i], power_index)) ;
	}
}

void NonNewton::computeViscosityCassonModel() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		float res = sqrt(m_muC) +  sqrt(m_yieldStress / m_strainRateNorm[i]);
		m_nonNewtonViscosity[i] = res * res;
	}
}


void NonNewton::computeViscosityCarreau() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		m_nonNewtonViscosity[i] = m_viscosity_inf +  (m_viscosity0 - m_viscosity_inf) / (1.0f +  pow(consistency_index * m_strainRateNorm[i]*m_strainRateNorm[i], (1.0f - power_index)/2.0f)) ;
	}
}



void NonNewton::computeViscosityBingham() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		if(m_strainRateNorm[i] < m_criticalStrainRate)
			m_nonNewtonViscosity[i] = m_viscosity0;
		else
		{
			float tau0 = m_criticalStrainRate * (m_viscosity0 - m_viscosity_inf);
			m_nonNewtonViscosity[i] = tau0 / m_strainRateNorm[i] + m_viscosity_inf;
		}
	}
}

void NonNewton::computeViscosityHerschelBulkley() 
{
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		if(m_strainRateNorm[i] < m_criticalStrainRate)
			m_nonNewtonViscosity[i] = m_viscosity0;
		else
		{
			float tau0 = m_viscosity0 * m_criticalStrainRate - consistency_index * pow(m_criticalStrainRate, power_index);
			m_nonNewtonViscosity[i] = tau0 / m_strainRateNorm[i] + consistency_index * pow(m_strainRateNorm[i], power_index - 1);
		}
	}
}


// void NonNewton::performNeighborhoodSearchSort()
// {
// 	Simulation *sim = Simulation::getCurrent();
// 	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
// 	d.sort_field(&m_nonNewtonViscosity[0]);
// } 


void NonNewton::controlledByTemperature()
{
	//读取温度
	auto* d = Partio::PartioSingleton::getCurrent();
	auto m_particleData = d->getParticlesData();

	Partio::ParticleAttribute tempAttr;
	m_particleData->attributeInfo("temperature", tempAttr);

	for (unsigned int i = 0; i < numParticles; ++i)
	{
		auto tempRef = m_particleData->data<float>(tempAttr, i);
		auto temp = *tempRef;
		if(temp < m_temperature_threshold)
		{
			m_nonNewtonViscosity[i] = m_viscosity0;
		}
		else
		{
			m_nonNewtonViscosity[i] = m_viscosity_inf;
		}
	}
} 



void NonNewton::smoothStrainRate()
{
	if (m_smoothStrainRateFactor == 0.0f)
		return;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int m = 0; m < nModels; m++)
	{
		FluidModel *fm = sim->getFluidModel(m);
		const unsigned int numParticles = fm->numActiveParticles();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int i = 0; i < numParticles; ++i)
			{
				const Vector3r &xi = fm->getPosition(i);
				// const Vector3r &vi = fm->getVelocity(i);
				const Real density_i = fm->getDensity(i);
				Vector6r avg = Vector6r::Zero();
				for (unsigned int j = 0; j < sim->numberOfNeighbors(m, m, i); j++)
				{
					const unsigned int neighborIndex = sim->getNeighbor(m, m, i, j);
					const Vector3r &xj = fm->getPosition(neighborIndex);
					const Real W = sim->W(xi - xj);
					const Real mj = fm->getMass(neighborIndex);
					const Real rhoj = fm->getDensity(neighborIndex);
					avg += m_strainRate[j] * mj / rhoj * W;
				}
				const Vector6r diff = m_strainRate[i] - avg;
				m_strainRate[i] = m_strainRate[i] - diff * m_smoothStrainRateFactor;
			}
		}
	}
}