#pragma once

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"
#include "SPlisHSPlasH/Viscosity/ViscosityBase.h"

namespace SPH
{
    class NonNewton_standard : public ViscosityBase
    {
    protected:
        std::vector<Real> m_viscosity_nonNewton;
        std::vector<Real> m_boundaryViscosity_nonNewton;

        virtual void initParameters();

    public:
        static int VISCOSITY_COEFFICIENT_BOUNDARY;

        NonNewton_standard(FluidModel *model);
        virtual ~NonNewton_standard(void);

        virtual void step();
        virtual void reset();

        static NonPressureForceBase* creator(FluidModel* model) { return new NonNewton_standard(model); }
    };
}