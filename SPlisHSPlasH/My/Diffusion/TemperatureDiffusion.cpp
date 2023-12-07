#include "TemperatureDiffusion.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/My/SurfaceParticles/SurfaceParticles.h"

using namespace SPH;
using namespace GenParam;
int TemperatureDiffusion::DIFFUSIVITY = -1;
int TemperatureDiffusion::R_SOURCE = -1;
int TemperatureDiffusion::POINT_SRC_VAL = -1;
int TemperatureDiffusion::POINT_SRC_POS = -1;
int TemperatureDiffusion::MELT_SURFACE = -1;
int TemperatureDiffusion::SURFACE_TEMP = -1;

TemperatureDiffusion::TemperatureDiffusion(FluidModel* model) :
    NonPressureForceBase(model),
    m_diffusivity(50.0),
    m_rSource(0.0)
{
    // model->setTemperature(1501, 100.0);


}


TemperatureDiffusion::~TemperatureDiffusion(void)
{
}

void TemperatureDiffusion::initParameters()
{
    POINT_SRC_VAL = createNumericParameter("pointSrcVal", "pointSrcVal", &m_pointSrcVal);
	setDescription(POINT_SRC_VAL, "set the coninuous point source temperature value(by particle ID)");
    setGroup(POINT_SRC_VAL, "diffusion");
    GenParam::RealParameter* rparam = static_cast<GenParam::RealParameter*>(getParameter(POINT_SRC_VAL));
    rparam->setMinValue(0.0);

    POINT_SRC_POS = createNumericParameter("pointSrcPos", "pointSrcPos", &m_pointSrcPos);
	setDescription(POINT_SRC_POS, "set the coninuous point source position(by particle ID)");
    setGroup(POINT_SRC_POS, "diffusion");

    MELT_SURFACE = createBoolParameter("MeltSurface", "MeltSurface", &m_meltSurface);
    setGroup(MELT_SURFACE, "diffusion");

    SURFACE_TEMP = createNumericParameter("surfaceTemp", "surfaceTemp", &m_surfaceTemp);
    setGroup(SURFACE_TEMP, "diffusion");
}

void TemperatureDiffusion::step()
{   
    static unsigned int steps = 0;
    steps++;
    Simulation* sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const Real density0 = m_model->getDensity0();
    FluidModel* model = m_model;
    const Real h = sim->getSupportRadius();
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();

    if(m_pointSrcPos!=-1)
        model->setTemperature(m_pointSrcPos, m_pointSrcVal);


    if (m_meltSurface && steps==1)
    {
        std::vector<int> surf_id; 
        findSurfaceParticles(m_model, surf_id);
        printf("\nSurface particles:\n");
        for (size_t i = 0; i < surf_id.size(); i++)
        {
            printf("%d\t", surf_id[i]);
            m_model->setTemperature(surf_id[i], m_surfaceTemp);
        }
        printf("\nTotal surface particles number: %d\n", surf_id.size());
    }

    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r &xi = m_model->getPosition(i);
            Real density_i = m_model->getDensity(i);
            
            Real &temp = model->getTemperature(i);
            Real temp_sum = 0.0;
            Real temp_old = temp;

            for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
            {
                const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                const Vector3r &xj = model->getPosition(neighborIndex);
                const Real grawWNorm = sim->gradW(xi - xj).norm();

                Real density_j = m_model->getDensity(neighborIndex);

                Real temp_j = model->getTemperature(neighborIndex);
                temp_sum += (m_diffusivity * m_model->getMass(neighborIndex) 
                / (density_j * density_i) * (temp_j - temp) * grawWNorm + m_rSource) * dt;
            }
            temp = temp_sum + temp_old;
        }
    }
}

void TemperatureDiffusion::reset()
{}

void TemperatureDiffusion::deferredInit()
{}

void TemperatureDiffusion::initValues()
{

}
