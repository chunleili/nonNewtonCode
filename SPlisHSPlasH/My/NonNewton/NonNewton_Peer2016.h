#pragma once

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"
#include "SPlisHSPlasH/Viscosity/ViscosityBase.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"

namespace SPH
{
    class NonNewton_Peer2016 : public ViscosityBase
    {
    protected:
        std::vector<Real> m_viscosity_nonNewton;

        std::vector<Real> m_density;
        std::vector<Matrix3r> m_targetNablaV;
        std::vector<Vector3r> m_omega;
        typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner1D> Solver;
        Solver m_solverV;
        Solver m_solverOmega;
        unsigned int m_iterationsV;
        unsigned int m_iterationsOmega;
        unsigned int m_maxIterV;
        Real m_maxErrorV;
        unsigned int m_maxIterOmega;
        Real m_maxErrorOmega;

        virtual void initParameters();
        void computeDensities();

    public:
        static int ITERATIONS_V;
        static int ITERATIONS_OMEGA;
        static int MAX_ITERATIONS_V;
        static int MAX_ERROR_V;
        static int MAX_ITERATIONS_OMEGA;
        static int MAX_ERROR_OMEGA;

        NonNewton_Peer2016(FluidModel *model);
        virtual ~NonNewton_Peer2016(void);

        static NonPressureForceBase* creator(FluidModel* model) { return new NonNewton_Peer2016(model); }

        virtual void step();
        virtual void reset();

        virtual void performNeighborhoodSearchSort();

        static void matrixVecProdV(const Real* vec, Real *result, void *userData);
        FORCE_INLINE static void diagonalMatrixElementV(const unsigned int row, Real &result, void *userData);

        static void matrixVecProdOmega(const Real* vec, Real *result, void *userData);
        FORCE_INLINE static void diagonalMatrixElementOmega(const unsigned int row, Real &result, void *userData);

        FORCE_INLINE const Matrix3r& getTargetNablaV(const unsigned int i) const
        {
            return m_targetNablaV[i];
        }

        FORCE_INLINE Matrix3r& getTargetNablaV(const unsigned int i)
        {
            return m_targetNablaV[i];
        }

        FORCE_INLINE void setTargetNablaV(const unsigned int i, const Matrix3r &val)
        {
            m_targetNablaV[i] = val;
        }

        FORCE_INLINE const Vector3r& getOmega(const unsigned int i) const
        {
            return m_omega[i];
        }

        FORCE_INLINE Vector3r& getOmega(const unsigned int i)
        {
            return m_omega[i];
        }

        FORCE_INLINE void setOmega(const unsigned int i, const Vector3r &val)
        {
            m_omega[i] = val;
        }
    };
}