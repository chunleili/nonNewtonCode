#ifndef __ParticleExporter_MyPartio_h__
#define __ParticleExporter_MyPartio_h__

#include "ExporterBase.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "extern/partio/src/lib/Partio.h"
#include <future>

namespace SPH
{
	/** \brief Particle exporter for the partio format.
	*/
	class ParticleExporter_MyPartio : public ExporterBase
	{
	public:
		Partio::ParticlesDataMutable* m_particleData;
	protected: 
		std::string m_exportPath;
		std::string m_particleFile;
		std::future<void> m_handle;

		void writeParticlesPartio(const std::string& fileName, FluidModel* model, const unsigned int objId=0xffffffff);
		void writeParticlesPartio_carried(const std::string& fileName, FluidModel* model, const unsigned int objId=0xffffffff);

	public:
		ParticleExporter_MyPartio(SimulatorBase *base);
		ParticleExporter_MyPartio(const ParticleExporter_MyPartio&) = delete;
        ParticleExporter_MyPartio& operator=(const ParticleExporter_MyPartio&) = delete;
		virtual ~ParticleExporter_MyPartio(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
