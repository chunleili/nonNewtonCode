/**
 * @file test_readWriteVector.cpp
 * @author chunleili
 * @brief 测试利用partio读取粒子文件(houdini bgeo)，然后存储到C++的std::vector当中。
 * 以及反过来，从C++的vector写入到partio，从而保存为bgeo文件。
 * @version 
 * @date 2022-11-28
 * 
 * 
 */
#include <string>
#include <vector>
#include "Partio.h"
#include <Eigen/Dense>
using namespace Partio;
using Eigen::Vector3f; //Eigen仅仅使用Vector3f，你也可以换std::array<float,3>之类的。

void readToVector(const std::string filename, std::vector<Vector3f> &positions);
void writeFromVector(const std::string filename, std::vector<Vector3f> &positions);

int main()
{
    //读入文件并构造ParticlesDataMutable对象data
    std::string filename = "E:/Dev/SPlisHSPlasH/extern/my_partio/data/cube.bgeo";
    std::vector<Vector3f> positions;
    readToVector(filename, positions);
}


void readToVector(const std::string filename, std::vector<Vector3f> &positions)
{
    // 先将文件读取到partio内部的数据结构
    Partio::ParticlesDataMutable* data=Partio::read(filename.c_str());
    // 粒子的数目为data->numParticles()
    std::cout<<"Reading partio particles, numParticles: "<<data->numParticles()<<"\n";
    positions.resize(data->numParticles());
    std::cout<<"positions size: "<<positions.size()<<"\n";

    // // 建立一个attribute作为存储粒子位置的attribute
    Partio::ParticleAttribute posAttr;
    posAttr = data->addAttribute("position", Partio::VECTOR, 3);
    // 遍历并拷贝粒子位置到positions(我们要存到的C++ vector)
    for (int i = 0; i < data->numParticles(); i++)
    {
    //     //从partio的数据结构中取出数据，得到的是C数组的形式
            const float* pos=data->data<float>(posAttr,i);
            positions[i][0] = pos[0];
            positions[i][1] = pos[1];
            positions[i][2] = pos[2];
            // std::cout<<"----\n";
            // std::cout<<positions[i]<<std::endl;
    }
    // // 释放内存
	data->release();
}


void writeFromVector(const std::string filename, std::vector<Vector3f> &positions)
{
    Partio::ParticlesDataMutable* data=Partio::read(filename.c_str());
    std::cout<<"Writing partio particles, numParticles: "<<data->numParticles()<<"\n";

    Partio::ParticleAttribute posAttr = data->addAttribute("position", Partio::VECTOR, 3);

    for (int i = 0; i < data->numParticles(); i++)
    {
        int idx = data->addParticle();
        float* p = data->dataWrite<float>(posAttr, idx);

        p[0] = positions[i][0];
        p[1] = positions[i][1];
        p[2] = positions[i][2];
    }
    Partio::write(filename.c_str(), *data);
	data->release();
}