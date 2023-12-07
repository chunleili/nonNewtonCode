
#include <Eigen/Eigen>
#include <string>
#include "ExporterBase.h"
#include "SPlisHSPlasH/FluidModel.h"
#include <fstream>

using Vector3R = Vector3f;
using namespace std;
using namespace SPH;

void outputData(Vector3R *input, unsigned int numParticles, std::string filename, long long int fluidNumL, int numL2)
{
    ofstream out;
    //此处需要自行新建model文件夹
    std::string path = "d:/model/";
    std::string file = path + filename;
    out.open(file, ios_base::trunc);
    Vector3R areaMin = base.getSimulationMethod().model.areaBoxMin;

    FluidModel &m = base.getSimulationMethod().model;

    std::cout << "outputing " << filename << "." << std::endl;
    unsigned int numP;
    numP = m.numActiveParticles();
    for (unsigned int i = 0; i < m.numActiveParticles(); i++)
    {
        if (i >= fluidNumL && i < fluidNumL + numL2)
        {
            numP--;
        }
    }
    //粒子数量
    out << numP << std::endl;
    //粒子半径大小
    // out << m.getParticleRadius() << std::endl;
    out << 0.025 << std::endl;

    //六面体场景边界的左下角坐标
    outputVector(m.areaBoxMin, out);
    //六面体场景边界的右上角坐标
    outputVector(m.areaBoxMax, out);
    //各个粒子的坐标值
    for (unsigned int i = 0; i < numParticles; i++)
    {
        if (i < fluidNumL || i >= fluidNumL + numL2)
        {
            outputVector(input[i], out);
        }
    }
    out.close();

    std::cout << "output " << filename << " Done." << std::endl
              << std::endl;
}
