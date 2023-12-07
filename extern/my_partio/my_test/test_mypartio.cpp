#include "Partio.h"
using namespace Partio;
int main()
{
    // 读取文件
    Partio::ParticlesDataMutable* data=Partio::read("E:/codes/try/partio/my_partio/test.bgeo");
    // 计算粒子数目
    std::cout<<"Number of particles "<<data->numParticles()<<std::endl;
}