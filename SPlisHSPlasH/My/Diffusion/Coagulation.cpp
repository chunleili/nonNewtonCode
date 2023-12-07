#include "Coagulation.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SurfaceParticles/SurfaceParticles.h"
#include "../Utils/mathUtils.h"

using namespace SPH;
using namespace GenParam;

int Coagulation::THRES_HIGH  = -1;
int Coagulation::THRES_LOW = -1;
int Coagulation::DIFFUSIVITY = -1;
int Coagulation::R_SOURCE = -1;
int Coagulation::BOX_MIN = -1;
int Coagulation::BOX_MAX = -1;
int Coagulation::POINT_SRC_VAL = -1;
int Coagulation::POINT_SRC_POS = -1;
int Coagulation::SURFACE_TEMP = -1;
int Coagulation::MELT_SURFACE = -1;
int Coagulation::SURFACE_SOURCE0 = -1;
int Coagulation::MELT_BOX = -1;
int Coagulation::MAX_VISCOSITY = -1;
int Coagulation::AVG_VISCOSITY = -1;
int Coagulation::AVG_TEMP = -1;


Coagulation::Coagulation(FluidModel* model) :
    NonPressureForceBase(model), 
    m_ccf(),
    m_thresHigh(static_cast<Real>(1.0)),
    m_thresLow(static_cast<Real>(0.1)),
    m_diffusivity(static_cast<Real>(50.0)),
    m_rSource(static_cast<Real>(0.1))
{
    m_ccf.resize(model->numParticles(), 0.0);
    m_viscosity.resize(model->numParticles(), 0.0);

    model->addField({ "ccf field", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_ccf[i]; } });

	m_boxMin.setZero();
	m_boxMax = Vector3r(0.3, 0.3, 0.3);

    // m_ccf[1501] = 100.0;
    m_meltSurface = false;
    m_surfaceTemp = 0.0;
    m_surfaceSource0 = 0.0;
    m_surfaceSource.resize(model->numParticles(), 0.0);

    m_isHotwater.resize(model->numParticles(), false);
    model->addField({ "isHotwater", FieldType::UInt, [&](const unsigned int i) -> int* { return &m_isHotwater[i]; } });
}



Coagulation::~Coagulation(void)
{
    m_model->removeFieldByName("ccf field");
    m_ccf.clear();
    m_model->removeFieldByName("isHotwater");
    m_isHotwater.clear();
}

void Coagulation::initParameters()
{
    THRES_HIGH = createNumericParameter("thresHigh", "threshold high", &m_thresHigh);
	setGroup(THRES_HIGH, "coagualtion");
	setDescription(THRES_HIGH, "higher treshold for coagualtion");
	GenParam::RealParameter* rparam = static_cast<GenParam::RealParameter*>(getParameter(THRES_HIGH));
	rparam->setMinValue(0.0);

    THRES_LOW = createNumericParameter("thresLow", "threshold low", &m_thresLow);
    setGroup(THRES_LOW, "coagualtion");
    setDescription(THRES_LOW, "lower treshold for coagualtion");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(THRES_LOW));
    rparam->setMinValue(0.0);

    DIFFUSIVITY = createNumericParameter("diffusivity", "diffusivity", &m_diffusivity);
    setGroup(DIFFUSIVITY, "coagualtion");
    setDescription(DIFFUSIVITY, "diffusivity of the convection-diffusion equation");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(DIFFUSIVITY));
    rparam->setMinValue(0.0);

    R_SOURCE = createNumericParameter("rSource", "rSource", &m_rSource);
    setGroup(R_SOURCE, "coagualtion");
    setDescription(R_SOURCE, "source term of the convection-diffusion equation");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(R_SOURCE));
    rparam->setMinValue(0.0);

    MELT_BOX = createBoolParameter("meltbox", "meltbox", &m_meltBox);
	setDescription(MELT_BOX, "wheter or not add source term in the box region");
    setGroup(MELT_BOX, "coagualtion");
    
    ParameterBase::GetVecFunc<Real> getFct = [&]()-> Real* { return m_boxMin.data(); };
	ParameterBase::SetVecFunc<Real> setFct = [&](Real* val)	
    {
	    m_boxMin = Vector3r(val[0], val[1], val[2]);
    };
	BOX_MIN = createVectorParameter("boxMin", "box min", 3u, getFct, setFct);
	setGroup(BOX_MIN, "coagualtion");
	setDescription(BOX_MIN, "Minimum point of box of which the rSource is not zero.");

	ParameterBase::GetVecFunc<Real> getFct2 = [&]()-> Real* { return m_boxMax.data(); };
    ParameterBase::SetVecFunc<Real> setFct2 = [&](Real* val)
    {
        m_boxMax = Vector3r(val[0], val[1], val[2]);
    };
	BOX_MAX = createVectorParameter("boxMax", "box max", 3u, getFct2, setFct2);
	setGroup(BOX_MAX, "coagualtion");
	setDescription(BOX_MAX, "Maximum point of box of which the rSource is not zero.");



    POINT_SRC_VAL = createNumericParameter("pointSrcVal", "pointSrcVal", &m_pointSrcVal);
	setDescription(POINT_SRC_VAL, "give a particle initial non-zero value to test the diffusion");
    setGroup(POINT_SRC_VAL, "coagualtion");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(POINT_SRC_VAL));
    rparam->setMinValue(0.0);

    POINT_SRC_POS = createNumericParameter("pointSrcPos", "pointSrcPos", &m_pointSrcPos);
	setDescription(POINT_SRC_POS, "give a particle ID(uint) of which the initial value is not zero to test the diffusion");
    setGroup(POINT_SRC_POS, "coagualtion");

    MELT_SURFACE = createBoolParameter("meltSurface", "meltSurface", &m_meltSurface);
    setGroup(MELT_SURFACE, "coagualtion");
    SURFACE_TEMP = createNumericParameter("surfaceTemp", "surfaceTemp", &m_surfaceTemp);
    setGroup(SURFACE_TEMP, "coagualtion");
    SURFACE_SOURCE0 = createNumericParameter("surfaceSource0", "surfaceSource0", &m_surfaceSource0);
    setGroup(SURFACE_SOURCE0, "coagualtion");


    MAX_VISCOSITY = createNumericParameter("max_viscosity", "max_viscosity", &m_maxViscosity);
	setGroup(MAX_VISCOSITY, "Viscosity");
	setDescription(MAX_VISCOSITY, "Max viscosity of all fluid particles.");
	getParameter(MAX_VISCOSITY)->setReadOnly(true);

	AVG_VISCOSITY = createNumericParameter("average_viscosity", "average_viscosity", &m_avgViscosity);
	setGroup(AVG_VISCOSITY, "Viscosity");
	setDescription(AVG_VISCOSITY, "Average viscosity of all fluid particles.");
	getParameter(AVG_VISCOSITY)->setReadOnly(true);

	AVG_TEMP = createNumericParameter("average_temperature", "average_temperature", &m_avgTemp);
	setGroup(AVG_TEMP, "Viscosity");
	setDescription(AVG_TEMP, "Average temperature of all fluid particles.");
	getParameter(AVG_TEMP)->setReadOnly(true);

    DECAY = createNumericParameter("decay", "decay", &m_decay);
    VISCOSITY0 = createNumericParameter("viscosity0", "viscosity0", &m_viscosity0);
	setGroup(DECAY, "coagualtion");
	setGroup(VISCOSITY0, "coagualtion");
	setDescription(DECAY, "decay index");
	setDescription(VISCOSITY0, "initial viscosity");

    HOT_WATER_TEMP = createNumericParameter("hotWaterTemp", "hotWaterTemp", &m_hotWaterTemp);
	setGroup(HOT_WATER_TEMP, "coagualtion");
	setDescription(HOT_WATER_TEMP, "hot water temp");

    IS_SHOWER_SCENE = createBoolParameter("isShowerScene", "isShowerScene", &m_isShowerScene);
	setGroup(IS_SHOWER_SCENE, "coagualtion");

    IS_ICE_CREAM_SCENE = createBoolParameter("isIceCreamScene", "isIceCreamScene", &m_isIceCreamScene);
    setGroup(IS_ICE_CREAM_SCENE, "coagualtion");

    IS_HOT_CUT_SCENE = createBoolParameter("isHotCutScene", "isHotCutScene", &m_isHotCutScene);
    setGroup(IS_HOT_CUT_SCENE, "coagualtion");
}

void Coagulation::step()
{
    steps++;

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

    if(m_pointSrcPos!=-1)
        m_ccf[m_pointSrcPos] = m_pointSrcVal;

    if (m_meltSurface && steps==1)
    {
        std::vector<int> surf_id; 
        findSurfaceParticles(m_model, surf_id);
        printf("\nSurface particles:\n");
        for (size_t i = 0; i < surf_id.size(); i++)
        {
            // printf("%d\t", surf_id[i]);
            m_ccf[surf_id[i]] = m_surfaceTemp;
        }
        std::cout<<"Total surface particles number: "<<surf_id.size()<<std::endl;
    }

    if (m_meltSurface)
    {
        std::vector<int> surf_id; 
        findSurfaceParticles(m_model, surf_id, 10);
        for (size_t i = 0; i < surf_id.size(); i++)
        {
            m_surfaceSource[surf_id[i]] = m_surfaceSource0;
        }
    }


    // compute ccf
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r &xi = m_model->getPosition(i);
            // Real &ccf_i = getCcf(i);
            Real &ccf_i = model->getTemperature(i);
            Real density_i = m_model->getDensity(i);
            Real source = 0.0;
            Real ccf_sum = 0.0;
            Real ccf_old = ccf_i;
            
            if(m_meltBox)
                if ((xi[0] > m_boxMin[0]) && (xi[1] > m_boxMin[1]) && (xi[2] > m_boxMin[2]) &&
                    (xi[0] < m_boxMax[0]) && (xi[1] < m_boxMax[1]) && (xi[2] < m_boxMax[2]))
                    source += m_rSource;
            source += m_surfaceSource[i];
            for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
            {
                const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                const Vector3r &xj = model->getPosition(neighborIndex);
                const Real grawWNorm = sim->gradW(xi - xj).norm();

                Real density_j = m_model->getDensity(neighborIndex);
                // Real ccf_j = getCcf(neighborIndex);
                Real ccf_j = model->getTemperature(neighborIndex);
                ccf_sum += (m_diffusivity * m_model->getMass(neighborIndex) 
                / (density_j * density_i) * (ccf_j - ccf_i) * grawWNorm + source) * dt;
            }
            ccf_i = ccf_sum + ccf_old;
            model->setTemperature(i, ccf_i);

            if(m_isIceCreamScene || m_isHotCutScene)
            {
                m_viscosity[i] = m_viscosity0 * exp(-m_decay * ccf_i); //冰激凌融化用的模型
                model->setNonNewtonViscosity(i, m_viscosity[i]);
            }

            if(m_isShowerScene)
            {
                // 如果是热水（也就是被emitter发射的粒子），则直接设置温度为100
                if(m_model->getParticleState(i) == ParticleState::AnimatedByEmitter)
                {
                    m_isHotwater[i] = 1;
                    model->setTemperature(i, m_hotWaterTemp);
                    m_viscosity[i] = 0.0;
                    model->setNonNewtonViscosity(i, 0.0);
                }
                if(m_isHotwater[i]==1)
                {   
                    // std::cout<<i<<" is hot water"<<std::endl;
                    model->setTemperature(i, m_hotWaterTemp);
                    m_viscosity[i] = 0.0;
                    model->setNonNewtonViscosity(i, 0.0);
                }
                if(m_isHotwater[i]!=1)
                {
                    // 增加直接控制粘度并传给Weiler
                    // m_viscosity[i] = m_viscosity0 * exp(-m_decay * ccf_i); //冰激凌融化用的模型

                    //热水用的模型 直接0或很大年度
                    if(model->getTemperature(i) > 1.0)
                        m_viscosity[i] = 0.0;
                    else
                        m_viscosity[i] = m_viscosity0;
                    model->setNonNewtonViscosity(i, m_viscosity[i]);
                }
            }

            Vector3r &vel_i = m_model->getVelocity(i);
            // 地板摩擦
            if(xi[1] < 0.01)
            {
                vel_i[0] *= 0.1;
                vel_i[2] *= 0.1;
            }
            
            //clamping
            if(vel_i.norm()>4)
                vel_i = vel_i.normalized()*4;
            
        }
    }

    m_maxViscosity = maxField(m_viscosity, numParticles);
	m_avgViscosity = avgField(m_viscosity, numParticles);
    m_avgTemp = avgField(m_ccf, numParticles);
}

void Coagulation::reset()
{
    steps=0;
}


void Coagulation::performNeighborhoodSearchSort()
{

    const unsigned int numPart = m_model->numActiveParticles();
    if (numPart == 0)
        return;

    Simulation* sim = Simulation::getCurrent();
    auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
    d.sort_field(&m_ccf[0]);

}

void Coagulation::deferredInit()
{
    initValues();
}

void Coagulation::initValues()
{

}
