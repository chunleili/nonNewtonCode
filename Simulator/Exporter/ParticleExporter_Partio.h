#ifndef __ParticleExporter_Partio_h__
#define __ParticleExporter_Partio_h__

#include "ExporterBase.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "extern/partio/src/lib/Partio.h"
#include <future>

namespace SPH
{
	/** \brief Particle exporter for the partio format.
	*/
	class ParticleExporter_Partio : public ExporterBase
	{
	protected: 
		std::string m_exportPath;
		std::string m_particleFile;
		Partio::ParticlesDataMutable* m_particleData;
		std::future<void> m_handle;

		void writeParticlesPartio(const std::string& fileName, FluidModel* model, const unsigned int objId=0xffffffff);

	public:
		ParticleExporter_Partio(SimulatorBase *base);
		ParticleExporter_Partio(const ParticleExporter_Partio&) = delete;
        ParticleExporter_Partio& operator=(const ParticleExporter_Partio&) = delete;
		virtual ~ParticleExporter_Partio(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
