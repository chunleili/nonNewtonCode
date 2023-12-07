#pragma once
#include <string>
#include "Partio.h"

namespace Partio
{

class PartioSingleton
{
public:
    static inline PartioSingleton *current = nullptr;
    static inline Partio::ParticlesDataMutable *partio_data = nullptr;
    // static  PartioSingleton *current;
    // static  Partio::ParticlesDataMutable *partio_data;

    static PartioSingleton *getCurrent()
    {
        if (current == nullptr)
        {
            current = new PartioSingleton();
        }
        return current;
    }
    PartioSingleton()
    {
        partio_data = Partio::create();
    }

    void read(const char *c_filename)
    {
        partio_data = Partio::read(c_filename);
    }

    void print()
    {
        Partio::print(partio_data);
    }

    ParticlesDataMutable *getParticlesData()
    {
        return partio_data;
    }
};
// PartioSingleton *PartioSingleton::current = nullptr;
// ParticlesDataMutable *PartioSingleton::partio_data = nullptr;

//usage: see test_heapdata.cpp

} // namespace Partio
