#pragma once
#include <vector>
#include <array>
#include <Eigen/Dense>
#include "partio.h"

using Eigen::Vector3f;

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
            std::cout<<"----\n";
            std::cout<<positions[i]<<std::endl;
    }
    // 释放内存
	data->release();
}