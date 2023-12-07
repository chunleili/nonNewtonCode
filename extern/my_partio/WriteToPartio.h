#pragma once
#include <vector>
#include <array>
#include "partio.h"

void inline writePartio(const std::string& particleFile,
    const std::vector<std::array<float,3>>& particles)
{
    // 创建可读写对象parts
    Partio::ParticlesDataMutable* parts = Partio::create();
    // 创建pos属性
    Partio::ParticleAttribute pos;
    // 设定属性名和格式
    pos = parts->addAttribute("position", Partio::VECTOR, 3);

    for (int i = 0; i < particles.size(); i++)
    {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(pos, idx);

        p[0] = particles[i][0];
        p[1] = particles[i][1];
        p[2] = particles[i][2];
    }

    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}
