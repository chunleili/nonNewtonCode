#pragma once

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataMyTimeStep.h"
#include "SPlisHSPlasH/SPHKernels.h"

#define USE_WARMSTART
#define USE_WARMSTART_V

namespace SPH
{
	class SimulationDataMyTimeStep;

	/** \brief 更改自DFSPH
	*/
	class MyTimeStep : public TimeStep
	{
	protected:
		SimulationDataMyTimeStep m_simulationData;
		unsigned int m_counter;
		const Real m_eps = static_cast<Real>(1.0e-5);
		bool m_enableDivergenceSolver;
		unsigned int m_iterationsV;
		Real m_maxErrorV;
		unsigned int m_maxIterationsV;

		void computeDFSPHFactor(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void divergenceSolve();
		void divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int index, const int numParticles, const Real h, const Real density0);
		void computeDensityChange(const unsigned int fluidModelIndex, const unsigned int index, const Real h);

#ifdef USE_WARMSTART_V
		void warmstartDivergenceSolve(const unsigned int fluidModelIndex);
#endif
#ifdef USE_WARMSTART
		void warmstartPressureSolve(const unsigned int fluidModelIndex);
#endif

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();
		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		virtual void initParameters();

	public:
		static int SOLVER_ITERATIONS_V;
		static int MAX_ITERATIONS_V;
		static int MAX_ERROR_V;
		static int USE_DIVERGENCE_SOLVER;
		inline static int DAMP_VELOCITY_FLAG = -1;
		inline static int DAMP_VELOCITY_FACTOR = -1;
		inline static int SMOOTH_VELOCITY_FLAG = -1;
		inline static int SMOOTH_VELOCITY_FACTOR = -1;

		MyTimeStep();
		virtual ~MyTimeStep(void);
		virtual void step();
		virtual void reset();
		virtual void resize();
		void dampVelocity();
		void smoothVelocity();

		inline static int m_timeStepNum = 0;
		bool m_dampVelocityFlag = false;
		float m_dampVelocityFactor = 0.0;
		bool m_smoothVelocityFlag = false;
		float m_smoothVelocityFactor = 0.0f;
	};
}
