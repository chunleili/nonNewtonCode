#include <vector>
#include <algorithm>
#include "mathUtils.h"
#include <tuple>
namespace SPH
{

    //对一个场求和
    float sumField(std::vector<float>& field, unsigned int numParticles)
    {
        float sum = 0.0f;
        for (unsigned int i = 0; i < numParticles; i++)
        {
            sum += field[i];
        }
        return sum;
    }

    float avgField(std::vector<float>& field, unsigned int numParticles)
    {
        return sumField(field, numParticles) / numParticles;
    }

    float maxField(std::vector<float>& field, unsigned int numParticles)
    {
        float max = field[0];
        for (unsigned int i = 1; i < numParticles; i++)
        {
            if (field[i] > max)
                max = field[i];
        }
        return max;
    }

    float minField(std::vector<float>& field, unsigned int numParticles)
    {
        float min = field[0];
        for (unsigned int i = 1; i < numParticles; i++)
        {
            if (field[i] < min)
                min = field[i];
        }
        return min;
    }

    auto minMaxField(std::vector<float>& field, unsigned int numParticles)
    {
        // auto [min, max] = std::minmax_element(field.begin(), field.end());
        // return std::make_tuple(*min, *max);
        float min = field[0];
        float max = field[0];
        for (unsigned int i = 1; i < numParticles; i++)
        {
            if (field[i] < min)
                min = field[i];
            if (field[i] > max)
                max = field[i];
        }
        return std::make_tuple(min, max);
    }

    auto minMaxAvgField(std::vector<float>& field, unsigned int numParticles)
    {
        // auto [min, max] = std::minmax_element(field.begin(), field.end());
        // return std::make_tuple(*min, *max);
        float min = field[0];
        float max = field[0];
        float sum = field[0];
        for (unsigned int i = 1; i < numParticles; i++)
        {
            if (field[i] < min)
                min = field[i];
            if (field[i] > max)
                max = field[i];
            sum += field[i];
        }
        auto avg = sum / numParticles;
        return std::make_tuple(min, max, avg);
    }
    
} // namespace SPH
