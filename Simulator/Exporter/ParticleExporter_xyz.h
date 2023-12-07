#ifndef __ParticleExporter_xyz_h__
#define __ParticleExporter_xyz_h__

#include "ExporterBase.h"
#include <fstream>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH
{
	/** \brief particle exporter, output the positions(xyz) of the particles
     * line1: particle number
     * line2: particle radius
     * line3: bounding box(minimum positions of all particles)
     * line3: bounding box(maximum positions of all particles)
     * after: particle positions
     * 
	*/
	class ParticleExporter_xyz : public ExporterBase
	{
	protected: 
		std::string m_exportPath;
        std::ofstream *m_outfile;

	public:
		ParticleExporter_xyz(SimulatorBase* base);
		ParticleExporter_xyz(const ParticleExporter_xyz&) = delete;
         ParticleExporter_xyz& operator=(const ParticleExporter_xyz&) = delete;
		virtual ~ParticleExporter_xyz(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}
#endif