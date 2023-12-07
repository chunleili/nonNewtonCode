#ifndef __Plastic_h__
#define __Plastic_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Elasticity/ElasticityBase.h"

namespace SPH
{
	/** \brief 
	 * Implemetation of plasticity
	 * Modified from the Elatticity_Becker2009
	*/
	class Plastic : public ElasticityBase
	{
	protected:
		// initial particle indices, used to access their original positions
		std::vector<unsigned int> m_current_to_initial_index;
		std::vector<unsigned int> m_initial_to_current_index;
		// initial particle neighborhood
		std::vector<std::vector<unsigned int>> m_initialNeighbors;
		// volumes in rest configuration
		std::vector<Real> m_restVolumes;
		std::vector<Matrix3r> m_rotations;
		std::vector<Vector6r> m_stress;
		std::vector<Matrix3r> m_F;

		//added variables
		std::vector<Vector6r> m_plasticStrain; //add Plastic Strain 
		std::vector<int> m_isFracture; // determine wheter a particle should be fracture and release 
		unsigned int numParticles=0;
		std::vector<Vector3r> m_displacement;
		std::vector<Vector6r> m_totalStrain;
		std::vector<Vector6r> m_elasticStrain;
		unsigned int m_step=0;

		void initValues();
		void computeRotations();
		void computeStress();
		void computeForces();

		void computePlasticStrain(int i, Vector6r & totalStrain);
		void computeNablaU(int i, Matrix3r &nablaU);
		void computeTotalStrain(Matrix3r &nablaU, Vector6r & totalStrain);
		void fracture();
		void computeC(Matrix6r &nablaU);

		virtual void initParameters();

		//////////////////////////////////////////////////////////////////////////
		// multiplication of symmetric matrix, represented by a 6D vector, and a 
		// 3D vector
		//////////////////////////////////////////////////////////////////////////
		FORCE_INLINE void symMatTimesVec(const Vector6r & M, const Vector3r & v, Vector3r &res)
		{
			res[0] = M[0] * v[0] + M[3] * v[1] + M[4] * v[2];
			res[1] = M[3] * v[0] + M[1] * v[1] + M[5] * v[2];
			res[2] = M[4] * v[0] + M[5] * v[1] + M[2] * v[2];
		}

		/**
		 * @brief compute the Frobenius norm of a symmetric 3x3matrix 
		 * represented by vector6r
		 * 
		 * @param vec : the symmetric matrix
		 * @return Real : the F norm
		 */
		FORCE_INLINE Real FNorm(Vector6r & vec)
		{
			Real FNorm = 0.0;
			for (int i = 0; i < 6; i++)
				FNorm += vec[i] * vec[i];
			FNorm = sqrt(FNorm)	;

			return FNorm;
		}


	public:
		//user defined parameters
		float elasticLimit;
		float plasticLimit;
		bool isPlastic=1;
		static int ELASTIC_LIMIT;
		static int PLASTIC_LIMIT;
		static int IS_PLASTIC;

		Plastic(FluidModel *model);
		virtual ~Plastic(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Plastic(model); }

		virtual void step();
		virtual void reset();
		virtual void performNeighborhoodSearchSort();

	};
}

#endif
