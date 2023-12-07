#pragma once

#include "SPlisHSPlasH/Common.h"
#include <vector>

namespace Utilities
{
	/** \brief Class for reading and writing partio files.
	 */
	class MyPartioReader
	{
	public:
		static bool MyPartioReader::readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale, std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities);
		static bool MyPartioReader::readParticlesToSingleton(const std::string &fileName);
	};

}
