#pragma once

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"
#include "SPlisHSPlasH/Viscosity/ViscosityBase.h"

namespace SPH {
    class NonNewton_Bender2017 : public ViscosityBase {
    protected:
        std::vector<Real> m_viscosity_nonNewton;

        std::vector<Vector6r> m_targetStrainRate;
        std::vector<Matrix6r> m_viscosityFactor;
        std::vector<Vector6r> m_viscosityLambda;
        unsigned int m_iterations;
        unsigned int m_maxIter;
        Real m_maxError;

        virtual void initParameters();

    public:
        static int ITERATIONS;
        static int MAX_ITERATIONS;
        static int MAX_ERROR;

        NonNewton_Bender2017(FluidModel *model);
        virtual ~NonNewton_Bender2017(void);

        static NonPressureForceBase* creator(FluidModel* model) { return new NonNewton_Bender2017(model); }

        virtual void step();
        virtual void reset();

        virtual void performNeighborhoodSearchSort();

        void computeTargetStrainRate();
        void computeViscosityFactor();

        /** Matrix product
        */
        FORCE_INLINE void viscoGradientMultTransposeRightOpt(const Eigen::Matrix<Real, 6, 3>& a, const Eigen::Matrix<Real, 6, 3>& b, Matrix6r &c)
        {
            // a(0,1), a(0,2), a(1,0), a(1,2), a(2,0), a(2,1), a(3,2), a(4, 1), a(5, 0) = 0
            c(0,0) = a(0,0) * b(0,0);
            c(0,1) = 0.0;
            c(0,2) = 0.0;
            c(0,3) = a(0,0) * b(3,0);
            c(0,4) = a(0,0) * b(4,0);
            c(0,5) = 0.0;

            c(1,0) = 0.0;
            c(1,1) = a(1,1) * b(1,1);
            c(1,2) = 0.0;
            c(1,3) = a(1,1) * b(3,1);
            c(1,4) = 0.0;
            c(1,5) = a(1,1) * b(5,1);

            c(2,0) = 0.0;
            c(2,1) = 0.0;
            c(2,2) = a(2,2) * b(2,2);
            c(2,3) = 0.0;
            c(2,4) = a(2,2) * b(4,2);
            c(2,5) = a(2,2) * b(5,2);

            c(3,0) = a(3,0) * b(0,0);
            c(3,1) = a(3,1) * b(1,1);
            c(3,2) = 0.0;
            c(3,3) = a(3,0) * b(3,0) + a(3,1) * b(3,1);
            c(3,4) = a(3,0) * b(4,0);
            c(3,5) = a(3,1) * b(5,1);

            c(4,0) = a(4,0) * b(0,0);
            c(4,1) = 0.0;
            c(4, 2) = a(4, 2) * b(2,2);
            c(4, 3) = a(4, 0) * b(3,0);
            c(4, 4) = a(4, 0) * b(4,0) + a(4,2) * b(4,2);
            c(4, 5) = a(4, 2) * b(5,2);

            c(5,0) = 0.0;
            c(5,1) = a(5,1) * b(1,1);
            c(5,2) = a(5,2) * b(2,2);
            c(5,3) = a(5,1) * b(3,1);
            c(5,4) = a(5,2) * b(4,2);
            c(5,5) = a(5,1) * b(5,1) + a(5,2) * b(5,2);
        }

        FORCE_INLINE const Vector6r& getTargetStrainRate(const unsigned int i) const
        {
            return m_targetStrainRate[i];
        }

        FORCE_INLINE Vector6r& getTargetStrainRate(const unsigned int i)
        {
            return m_targetStrainRate[i];
        }

        FORCE_INLINE void setTargetStrainRate(const unsigned int i, const Vector6r &val)
        {
            m_targetStrainRate[i] = val;
        }

        FORCE_INLINE const Matrix6r& getViscosityFactor(const unsigned int i) const
        {
            return m_viscosityFactor[i];
        }

        FORCE_INLINE Matrix6r& getViscosityFactor(const unsigned int i)
        {
            return m_viscosityFactor[i];
        }

        FORCE_INLINE void setViscosityFactor(const unsigned int i, const Matrix6r &val)
        {
            m_viscosityFactor[i] = val;
        }

        FORCE_INLINE const Vector6r& getViscosityLambda(const unsigned int i) const
        {
            return m_viscosityLambda[i];
        }

        FORCE_INLINE Vector6r& getViscosityLambda(const unsigned int i)
        {
            return m_viscosityLambda[i];
        }

        FORCE_INLINE void setViscosityLambda(const unsigned int i, const Vector6r &val)
        {
            m_viscosityLambda[i] = val;
        }
    };
}