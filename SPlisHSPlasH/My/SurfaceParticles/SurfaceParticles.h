#pragma once

#include <vector>
#include "../../FluidModel.h"
//获取哪些粒子是表面粒子
namespace SPH
{
void findSurfaceParticles(FluidModel *model, std::vector<int> &id, unsigned int threshold=10);

} // namespace SPH