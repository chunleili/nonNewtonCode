#pragma once
#ifndef __Coagulation_h__
#define __Coagulation_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	class Coagulation : public NonPressureForceBase
	{
	protected:
		virtual void initParameters();

		Real m_thresLow;
		Real m_thresHigh;
		Real m_diffusivity;
		Real m_rSource;

		Real m_pointSrcVal=0.0;
		int m_pointSrcPos=-1;

		Vector3r m_boxMin;
		Vector3r m_boxMax;
		bool m_meltBox=true;

		std::vector<Real> m_ccf;
		Real m_surfaceTemp;
		bool m_meltSurface;
		Real m_surfaceSource0=0.0; //用户指定的表面源值
		std::vector<Real> m_surfaceSource;
		int steps = 0;
		Real m_viscosity0 = 1000.0f;
		Real m_decay = 0.1;
		std::vector<Real> m_viscosity;
		float m_maxViscosity = 0.0f;
		float m_avgViscosity = 0.0f;
		float m_avgTemp = 0.0f;
		std::vector<int> m_isHotwater; //1代表是热水，其他代表不是
		// std::vector<std::vector<int>> m_neighbors;
		Real m_hotWaterTemp = 0.0;
		bool m_isIceCreamScene = true;
		bool m_isShowerScene = false;
		bool m_isHotCutScene = false;

		virtual void deferredInit();

		void initValues();

	public:
		static int THRES_LOW;
		static int THRES_HIGH;
		static int DIFFUSIVITY;
		static int R_SOURCE;
		static int POINT_SRC_VAL;
		static int POINT_SRC_POS;
		static int BOX_MIN;
		static int BOX_MAX;
		static int MELT_BOX;
		static int SURFACE_TEMP;
		static int MELT_SURFACE;
		static int SURFACE_SOURCE0;
		static int MAX_VISCOSITY;
		static int AVG_VISCOSITY;
		static int AVG_TEMP;
		inline static int DECAY = -1;
		inline static int VISCOSITY0 = -1;
		inline static int HOT_WATER_TEMP = -1;
		inline static int IS_SHOWER_SCENE = -1;
		inline static int IS_ICE_CREAM_SCENE = -1;
		inline static int IS_HOT_CUT_SCENE = -1;

		Coagulation();
		Coagulation(FluidModel* model);
		virtual ~Coagulation(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Coagulation(model); }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

		Real& getCcf(const unsigned int i)
		{
			return m_ccf[i];
		}

	};
}

#endif