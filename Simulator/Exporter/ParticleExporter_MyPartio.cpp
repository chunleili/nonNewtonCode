#include "ParticleExporter_MyPartio.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
// #include "extern/partio/src/lib/Partio.h"
#include "extern/my_partio/Partio.h"
#include "extern/my_partio/PartioSingleton.h"

using namespace SPH;
using namespace Utilities;

ParticleExporter_MyPartio::ParticleExporter_MyPartio(SimulatorBase *base) :
	ExporterBase(base)
{
}

ParticleExporter_MyPartio::~ParticleExporter_MyPartio(void)
{
}

void ParticleExporter_MyPartio::init(const std::string& outputPath)
{
	m_exportPath = FileSystem::normalizePath(outputPath + "/mypartio");
}

void ParticleExporter_MyPartio::step(const unsigned int frame)
{
	if (!m_active)
		return;

	Simulation* sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel* model = sim->getFluidModel(i);
		std::string fileName = "ParticleData";
		if (!m_base->getValue<bool>(SimulatorBase::EXPORT_OBJECT_SPLITTING))
		{
			fileName = fileName + "_" + model->getId() + "_" + std::to_string(frame);
			std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);
			if (m_base->getValue<bool>(SimulatorBase::USE_CARRIED_PARTIO_DATA))
				writeParticlesPartio_carried(exportFileName + ".bgeo.gz", model);
			else
				writeParticlesPartio(exportFileName + ".bgeo.gz", model);
		}
		else
		{
			// object splitting
			for (auto j = 0u; j < m_base->getLastObjectId(); j++)
			{
				std::string fileName2 = fileName + "_" + model->getId() + "_" + std::to_string(j) + "_" + std::to_string(frame);
				std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName2);
				if (m_base->getValue<bool>(SimulatorBase::USE_CARRIED_PARTIO_DATA))
					writeParticlesPartio_carried(exportFileName + ".bgeo.gz", model, j);
				else
					writeParticlesPartio(exportFileName + ".bgeo.gz", model, j);
			}
		}
	}
}

void ParticleExporter_MyPartio::reset()
{
}

void ParticleExporter_MyPartio::setActive(const bool active)
{
	ExporterBase::setActive(active);
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}


void ParticleExporter_MyPartio::writeParticlesPartio_carried(const std::string& fileName, FluidModel* model, const unsigned int objId)
{	
	auto* d = Partio::PartioSingleton::getCurrent();
	m_particleData = d->getParticlesData();

    Partio::ParticleAttribute posAttr;
    m_particleData->attributeInfo("position", posAttr);

	//根据计算结果（存储在model中），更新粒子位置
	//其余的不需要更新了，想更新什么，就从model中取出来，然后更新到单例中
	// BUG FIX: 由于partio中的粒子是从bhclassic读入的，数量是固定的。当model中额外注入了粒子的时候(比如fluid block或者emitter或其他来源)，就会导致partio中的粒子数量不够用，从而导致程序崩溃。
	unsigned int numParInPartio = m_particleData->numParticles();

    for (unsigned int i = 0; i < model->numActiveParticles(); i++)
    {
		int idx = i;

		//当粒子都是partio读入的时候，不需要额外添加粒子
		//当model中的粒子数量大于partio中的粒子数量时，需要额外添加粒子
		if(i >= numParInPartio)
		{
			int idx2 = m_particleData->addParticle();
		}

        float* p = m_particleData->dataWrite<float>(posAttr, idx);

		const Vector3r& x = model->getPosition(i);
        p[0] = x[0];
        p[1] = x[1];
        p[2] = x[2];
    }
    Partio::write(fileName.c_str(), *m_particleData);
}


void ParticleExporter_MyPartio::writeParticlesPartio(const std::string& fileName, FluidModel* model, const unsigned int objId)
{	
	const bool async = m_base->getValue<bool>(SimulatorBase::ASYNC_EXPORT);
	if (async)
	{
		if (m_handle.valid())
			m_handle.wait();
	}

	m_particleData = Partio::create();
	Partio::ParticlesDataMutable& particleData = *m_particleData;

	Partio::ParticleAttribute posAttr = particleData.addAttribute("position", Partio::VECTOR, 3);
	Partio::ParticleAttribute idAttr = particleData.addAttribute("id", Partio::INT, 1);

	// add attributes
	std::vector<std::string> attributes;
	StringTools::tokenize(m_base->getValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES), attributes, ";");

	std::map<unsigned int, int> attrMap;
	std::map<unsigned int, Partio::ParticleAttribute> partioAttrMap;
	for (unsigned int i = 0; i < attributes.size(); i++)
	{
		// position is exported anyway
		if (attributes[i] == "position")
		{
			attrMap[i] = -1;
			continue;
		}

		bool found = false;
		for (unsigned int j = 0; j < model->numberOfFields(); j++)
		{
			const FieldDescription& field = model->getField(j);
			if (field.name == attributes[i])
			{
				found = true;
				if (field.type == Scalar)
				{
					attrMap[i] = j;
					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::FLOAT, 1);
				}
				else if (field.type == UInt)
				{
					attrMap[i] = j;
					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::INT, 1);
				}
				else if (field.type == Vector3)
				{
					attrMap[i] = j;
					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::VECTOR, 3);
				}
				else if (field.type == Matrix3)
				{
					attrMap[i] = j;
					partioAttrMap[i] = particleData.addAttribute(attributes[i].c_str(), Partio::FLOAT, 9);
				}
				else
				{
					attrMap[i] = -1;
					LOG_WARN << "Only scalar and vector and matrix3 fields are currently supported by the partio exporter.";
				}
				break;
			}
		}
		if (!found)
		{
			attrMap[i] = -1;
			LOG_WARN << "Unknown field cannot be exported in partio file: " << attributes[i];
		}
	}

	const unsigned int numParticles = model->numActiveParticles();

	for (unsigned int i = 0; i < numParticles; i++)
	{
		if ((objId != 0xffffffff) && (model->getObjectId(i) != objId))
			continue;
			
		Partio::ParticleIndex index = particleData.addParticle();
		float* pos = particleData.dataWrite<float>(posAttr, index);
		int* id = particleData.dataWrite<int>(idAttr, index);

		const Vector3r& x = model->getPosition(i);
		pos[0] = (float)x[0];
		pos[1] = (float)x[1];
		pos[2] = (float)x[2];

		id[0] = model->getParticleId(i);

		for (unsigned int j = 0; j < attributes.size(); j++)
		{
			const int fieldIndex = attrMap[j];
			if (fieldIndex != -1)
			{
				const FieldDescription& field = model->getField(fieldIndex);
				if (field.type == FieldType::Scalar)
				{
					float* val = particleData.dataWrite<float>(partioAttrMap[j], index);
					*val = (float)*((Real*)field.getFct(i));
				}
				else if (field.type == FieldType::UInt)
				{
					int* val = particleData.dataWrite<int>(partioAttrMap[j], index);
					*val = (int)*((unsigned int*)field.getFct(i));
				}
				else if (field.type == FieldType::Vector3)
				{
					float* val = particleData.dataWrite<float>(partioAttrMap[j], index);
					Eigen::Map<Vector3r> vec((Real*)field.getFct(i));
					val[0] = (float)vec[0];
					val[1] = (float)vec[1];
					val[2] = (float)vec[2];
				}
			}
		}
	}

	m_particleFile = fileName;
	if (async)
		m_handle = std::async(std::launch::async, [&] { Partio::write(m_particleFile.c_str(), particleData, true); particleData.release(); });
	else
	{
		Partio::write(m_particleFile.c_str(), particleData, true);
		particleData.release();
	}
}