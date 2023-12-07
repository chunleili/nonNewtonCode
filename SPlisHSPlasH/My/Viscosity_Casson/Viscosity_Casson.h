#ifndef __Viscosity_Casson_h__
#define __Viscosity_Casson_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Viscosity/ViscosityBase.h"
#include "SPlisHSPlasH/Utilities/MatrixFreeSolver.h"

#define USE_BLOCKDIAGONAL_PRECONDITIONER

namespace SPH
{
	/** \brief viscosity_casson的屎山代码是如何运行的
		mu=0.01 
		mub=0 所以我把boundary相关的全删除了
		cassonViscosity=10000
		threshold=1(也就是温度大于1认为就是融化了)
	*/
	class Viscosity_Casson : public ViscosityBase
	{
	protected:
		Real m_boundaryViscosity;
		unsigned int m_maxIter;
		Real m_maxError;
		unsigned int m_iterations;
		std::vector<Vector3r> m_vDiff;
		Real m_tangentialDistanceFactor;
		static  Vector3r m_boxMin;
		static  Vector3r m_boxMax;

#ifdef USE_BLOCKDIAGONAL_PRECONDITIONER
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, BlockJacobiPreconditioner3D> Solver;
		FORCE_INLINE static void diagonalMatrixElement(const unsigned int row, Matrix3r &result, void *userData);
#else
		typedef Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, JacobiPreconditioner3D> Solver;
		FORCE_INLINE static void diagonalMatrixElement(const unsigned int row, Vector3r &result, void *userData);
#endif	
		Solver m_solver;


		virtual void initParameters();
		
	public:
		static int ITERATIONS;
		static int MAX_ITERATIONS;
		static int MAX_ERROR;
		static int VISCOSITY_COEFFICIENT_BOUNDARY;
		inline static int MAX_VISCOSITY_CASSON = -1;
		inline static int AVG_VISCOSITY_CASSON = -1;
		inline static int THRESHOLD = -1;
		inline static int VISCOSITY0 = -1;

		std::vector<Real> m_viscosity_nonNewton;
		Real m_maxViscosity_casson = 0.0f;
		Real m_avgViscosity_casson = 0.0f;
		Real m_threshold = 1.0f;
		Real m_viscosity0 = 1000.0f;

		Viscosity_Casson(FluidModel *model);
		virtual ~Viscosity_Casson(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Viscosity_Casson(model); }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		static void matrixVecProd(const Real* vec, Real *result, void *userData);

		FORCE_INLINE const Vector3r& getVDiff(const unsigned int i) const
		{
			return m_vDiff[i];
		}

		FORCE_INLINE Vector3r& getVDiff(const unsigned int i)
		{
			return m_vDiff[i];
		}

		FORCE_INLINE void setVDiff(const unsigned int i, const Vector3r& val)
		{
			m_vDiff[i] = val;
		}

        void computeRHS(VectorXr &b, VectorXr &g);

        void applyForces(const VectorXr &x);
    };
}

#endif
