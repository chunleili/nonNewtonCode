#include "Melting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../Utils/mathUtils.h"
#include "extern/my_partio/Partio.h"
#include "extern/my_partio/PartioSingleton.h"

using namespace SPH;
using namespace GenParam;

int Melting::THRES_HIGH  = -1;
int Melting::THRES_LOW = -1;
int Melting::DIFFUSIVITY = -1;
int Melting::R_SOURCE = -1;
int Melting::POINT_SRC_IDX = -1;



Melting::Melting(FluidModel* model) :
    NonPressureForceBase(model), 
    m_temp(),
    m_thresHigh(static_cast<Real>(1.0)),
    m_thresLow(static_cast<Real>(0.1)),
    m_diffusivity(static_cast<Real>(50.0)),
    m_rSource(static_cast<Real>(0.1))
{
    m_viscosity.resize(model->numParticles(), 0.0);
    m_temp.resize(model->numParticles(), 0.0);
}



Melting::~Melting(void)
{
}

void Melting::initParameters()
{
    THRES_HIGH = createNumericParameter("thresHigh", "threshold high", &m_thresHigh);
	setGroup(THRES_HIGH, "melting");
	setDescription(THRES_HIGH, "higher treshold for melting");
	GenParam::RealParameter* rparam = static_cast<GenParam::RealParameter*>(getParameter(THRES_HIGH));
	rparam->setMinValue(0.0);

    THRES_LOW = createNumericParameter("thresLow", "thresLow", &m_thresLow);
    setGroup(THRES_LOW, "melting");
    setDescription(THRES_LOW, "lower treshold for melting");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(THRES_LOW));
    rparam->setMinValue(0.0);

    DIFFUSIVITY = createNumericParameter("diffusivity", "diffusivity", &m_diffusivity);
    setGroup(DIFFUSIVITY, "melting");
    setDescription(DIFFUSIVITY, "diffusivity of the convection-diffusion equation");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(DIFFUSIVITY));
    rparam->setMinValue(0.0);

    R_SOURCE = createNumericParameter("rSource", "rSource", &m_rSource);
    setGroup(R_SOURCE, "melting");
    setDescription(R_SOURCE, "source term of the convection-diffusion equation");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(R_SOURCE));
    rparam->setMinValue(0.0);

    POINT_SRC_IDX = createNumericParameter("pointSrcIdx", "pointSrcIdx", &m_pointSrcIdx);
	setDescription(POINT_SRC_IDX, "give a particle ID(uint) of which the initial value is not zero to test the diffusion");
    setGroup(POINT_SRC_IDX, "melting");

    DECAY = createNumericParameter("decay", "decay", &m_decay);
	setGroup(DECAY, "melting");
	setDescription(DECAY, "decay index");


    VISCOSITY0 = createNumericParameter("viscosity0", "viscosity0", &m_viscosity0);
	setGroup(VISCOSITY0, "melting");
	setDescription(VISCOSITY0, "initial viscosity");
}

void Melting::step()
{
    if(m_steps==0)
    {
        printf("Reading temperature from Partio! \n");
        readTemperature();
    }

    m_steps++;

    Simulation* sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const Real density0 = m_model->getDensity0();
    FluidModel* model = m_model;
    const Real h = sim->getSupportRadius();
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();


    // compute temperature
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r &xi = m_model->getPosition(i);
            Real &temp_i = model->getTemperature(i);
            Real density_i = m_model->getDensity(i);
            Real source = 0.0;
            Real temp_sum = 0.0;
            Real temp_old = temp_i;
            

            // 已经燃烧的，就让它燃烧得更旺一些吧！
            if(temp_i > m_thresLow)
                source = m_rSource; //Add fuel

            for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
            {
                const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                const Vector3r &xj = model->getPosition(neighborIndex);
                const Real grawWNorm = sim->gradW(xi - xj).norm();

                Real density_j = m_model->getDensity(neighborIndex);
                Real temp_j = model->getTemperature(neighborIndex);
                temp_sum += (m_diffusivity * m_model->getMass(neighborIndex) 
                / (density_j * density_i) * (temp_j - temp_i) * grawWNorm + source) * dt;
            }
            temp_i = temp_sum + temp_old;
            model->setTemperature(i, temp_i);


            //冰激凌融化：decay越大，粘度掉的越快。从mu0开始掉。
            if (temp_i > m_thresLow)
            {
                m_viscosity[i] = m_viscosity0 * exp(-m_decay * temp_i); 
                model->setNonNewtonViscosity(i, m_viscosity[i]);
            }
        }
    }
}

void Melting::reset()
{
    m_steps=0;
}


void Melting::performNeighborhoodSearchSort()
{

    const unsigned int numPart = m_model->numActiveParticles();
    if (numPart == 0)
        return;

    Simulation* sim = Simulation::getCurrent();
    auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
    d.sort_field(&m_temp[0]);
}

void Melting::deferredInit()
{
    initValues();
}

void Melting::initValues()
{}


void Melting::readTemperature()
{
    // 从初始化的partio文件中读取温度
    m_numParticles = m_model->numActiveParticles();

    //  get data from the partio
    auto* d = Partio::PartioSingleton::getCurrent();
    auto m_particleData = d->getParticlesData();

	Partio::ParticleAttribute tempAttr;
	m_particleData->attributeInfo("temperature", tempAttr);

    // compute temperature
    for (int i = 0; i < (int)m_numParticles; i++)
    {
        auto tempRef = m_particleData->data<float>(tempAttr, i);
        m_model->setTemperature(i, tempRef[0]);
    }
}
