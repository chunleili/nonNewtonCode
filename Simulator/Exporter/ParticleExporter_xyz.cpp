#include "ParticleExporter_xyz.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
#include <limits>

using namespace SPH;
using namespace Utilities;
using namespace std;

ParticleExporter_xyz::ParticleExporter_xyz(SimulatorBase* base) :
	ExporterBase(base)
{
    m_outfile = nullptr;
}

ParticleExporter_xyz::~ParticleExporter_xyz(void)
{
}

void ParticleExporter_xyz::init(const std::string& outputPath)
{
    // define output path for the data
	m_exportPath = FileSystem::normalizePath(outputPath + "/xyz");
}

void ParticleExporter_xyz::step(const unsigned int frame)
{
    // check if the exporter is active
	if (!m_active)
		return;

	Simulation* sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel* model = sim->getFluidModel(i);

        //define the exportFileName
		std::string fileName = "Fluid";
        fileName = fileName + "_"  + std::to_string(frame);
        std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName + ".in");
        
        // Define the fstream object m_outfile. Open the file
        m_outfile = new std::ofstream(exportFileName, std::ios_base::trunc);
        if (!m_outfile->is_open())
        {
            LOG_WARN << "Cannot open a file to save VTK particles.";
            return;
        }

        //output particle number
        const unsigned int numParticles = model->numActiveParticles();
        *m_outfile << numParticles << std::endl;

        //output particle radius
        *m_outfile << sim->getParticleRadius() << std::endl;

        //output positions
        for(unsigned int i = 0; i < numParticles; i++)
        {
            Vector3f pos = model->getPosition(i);

            *m_outfile << pos.x() <<" " <<pos.y() <<" "<<pos.z() << std::endl;
        }
    }

    m_outfile->close();
    delete m_outfile;
    m_outfile = nullptr;

}

void ParticleExporter_xyz::reset()
{
}

void ParticleExporter_xyz::setActive(const bool active)
{
	ExporterBase::setActive(active);
    // create output folder
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}
