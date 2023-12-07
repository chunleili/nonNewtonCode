#include "PlyReader.h"
#include "extern/happly/happly.h"
#include "FileSystem.h"

using namespace Utilities;

bool PlyReader::readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
	std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities)
{
    if (!FileSystem::fileExists(fileName))
		return false;
    std::cout<<"Reading Ply "<<fileName<<"\n";
    // Construct the data object by reading from file
    happly::PLYData plyIn(fileName);

    // Get mesh-style data from the object
    std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
    // std::vector<std::vector<size_t>> fInd = plyIn.getFaceIndices<size_t>();
    
    //copy vpos to postions
    positions.resize(vPos.size());
    for (int i = 0; i < vPos.size(); i++)
    {
        Vector3r x(vPos[i][0], vPos[i][1], vPos[i][2]);
        x = rotation * (x*scale) + translation;
        positions[i] = x;
    }
    
    velocities.resize(vPos.size());
    for (int i = 0; i < vPos.size(); i++)
    {
        Vector3r x(0, 0, 0);
        velocities[i] = x;
    }

    std::cout<<"End Reading Ply\n";
    return true;
}