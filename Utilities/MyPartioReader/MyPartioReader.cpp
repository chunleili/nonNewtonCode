#include "MyPartioReader.h"
#include "extern/my_partio/Partio.h"
#include "../FileSystem.h"
#include "extern/my_partio/PartioSingleton.h"

using namespace Utilities;
using namespace Partio;

bool MyPartioReader::readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
	std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities)
{
	//大部分拷贝自PartioReaderWriter.cpp
	if (!FileSystem::fileExists(fileName))
		return false;

	// Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());

	//MYADD
	// std::cout<<"Reading Partio and save in the singleton\n";
	auto *d = PartioSingleton::getCurrent();
	d->read(fileName.c_str());
	auto* data = d->getParticlesData();
	// print(data)	;
	// std::cout<<"End reading Partio and save in the singleton\n";


	if (!data)
		return false;

	unsigned int posIndex = 0xffffffff;
	unsigned int velIndex = 0xffffffff;

	for (int i = 0; i < data->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		data->attributeInfo(i, attr);
		if (attr.name == "position")
			posIndex = i;
		else if (attr.name == "velocity")
			velIndex = i;
	}
	

	Partio::ParticleAttribute attr;

	if (posIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int) positions.size();
		positions.resize(fSize + data->numParticles());
		data->attributeInfo(posIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *pos = data->data<float>(attr, i);
			Vector3r x(pos[0], pos[1], pos[2]);
			x = rotation * (x*scale) + translation;
			positions[i + fSize] = x;
		}
	}

	if (velIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int) velocities.size();
		velocities.resize(fSize + data->numParticles());
		data->attributeInfo(velIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *vel = data->data<float>(attr, i);
			Vector3r v(vel[0], vel[1], vel[2]);
			velocities[i + fSize] = v;
		}
	}
	else
	{
		unsigned int fSize = (unsigned int) velocities.size();
		velocities.resize(fSize + data->numParticles());
		for (int i = 0; i < data->numParticles(); i++)
			velocities[i + fSize].setZero();
	}

	// data->release();

	return true;
}
