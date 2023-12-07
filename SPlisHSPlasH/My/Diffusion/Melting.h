#pragma once
#ifndef __Melting_h__
#define __Melting_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	class Melting : public NonPressureForceBase
	{
	protected:
		virtual void initParameters();
		int m_numParticles;

		inline static int m_steps = 0;
		Real m_thresLow;
		Real m_thresHigh;
		Real m_diffusivity;
		Real m_rSource;
		Real m_pointSrcVal=0.0;
		int m_pointSrcIdx=-1;
		Real m_viscosity0 = 1000.0f;
		Real m_decay = 0.1;

		std::vector<Real> m_temp;
		std::vector<Real> m_viscosity;


		virtual void deferredInit();

		void initValues();

	public:
		static int THRES_LOW;
		static int THRES_HIGH;
		static int DIFFUSIVITY;
		static int R_SOURCE;
		static int POINT_SRC_IDX;
		inline static int DECAY = -1;
		inline static int VISCOSITY0 = -1;


		Melting();
		Melting(FluidModel* model);
		virtual ~Melting(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Melting(model); }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		void readTemperature();
	};
}

#endif