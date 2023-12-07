#pragma once

#include "SPlisHSPlasH/Common.h"
#include <vector>

namespace Utilities
{
	/** \brief Class for reading ply files.
	 */
	class PlyReader
	{
	public:
		static bool PlyReader::readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale, std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities);
	};

}
