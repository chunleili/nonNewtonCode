#include <math.h>
#include "ShapeMatching.h"
#include "SPlisHSPlasH/TimeManager.h"
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include "SPlisHSPlasH/Utilities/MathFunctions.h"

#define _USE_MATH_DEFINES

using namespace SPH;
using namespace std;

inline Matrix3r polarDecompose(Matrix3r A_pq);

ShapeMatching::ShapeMatching(FluidModel* model):
NonPressureForceBase(model)
{
    unsigned int numParticles = m_model->numActiveParticles();
    radius_vector.resize(numParticles);
    oldPosition.resize(numParticles);
}

void ShapeMatching::setStates()
{
    unsigned int numParticles = m_model->numActiveParticles();
    for(unsigned int i=0; i< numParticles ; i++)
    {
        m_model->setParticleState(i, ParticleState::ShapeMatching);
    }   
}

ShapeMatching::~ShapeMatching()
{
}

void ShapeMatching::step()
{   
    static int steps = 0;
    if(steps == 0)
        computeBarycenter();
    setStates();
    shapeMatching();
    steps++;
}

void ShapeMatching::shapeMatching()
{
    const int numParticles = (int)m_model->numActiveParticles();
    TimeManager *tm = TimeManager::getCurrent ();
	const Real dt = tm->getTimeStepSize();

    Vector3r g(0.0, -9.8, 0.0);

    for(int i = 0; i < numParticles; i++)
    {
        if (m_model->getParticleState(i) == ParticleState::ShapeMatching)
        {
            Vector3r &pos = m_model->getPosition(i);
            Vector3r &vel = m_model->getVelocity(i);
            Vector3r f ;
            Real &mass = m_model->getMass(i);
            Real mass_inv = 1.0/mass;
            
            oldPosition[i] = pos;
            f = g ;
            vel +=  f * mass_inv * dt;
            pos += vel *dt;
            if(pos[1] < 1e-5)
            {
                pos = oldPosition[i];
                pos[1] = 0.0 + 1e-5;
            }
        }
    }

    //compute the new(matched shape) mass center
    Vector3r c{0.0, 0.0, 0.0};
    for(int i = 0; i < numParticles; i++)
    {
        if (m_model->getParticleState(i) == ParticleState::ShapeMatching)
        {
            Vector3r &pos = m_model->getPosition(i);
            c += pos;
        }
    }
    c /= numParticles;
    // cout<<"c:"<<c<<endl;

    // compute transformation matrix and extract rotation
    Matrix3r A, A_pq;
    A.setZero();
    A_pq.setZero();
    for(int i = 0; i < numParticles; i++)
    {
        Vector3r &pos = m_model->getPosition(i);

        A_pq += (pos - c) * radius_vector[i].transpose();
    }
    // A = A_pq * A_qq;

    //polar decomposition
    Matrix3r R;
    R = polarDecompose(A_pq);
    // R.setIdentity();

    // update vel and pos
    for(int i = 0; i < numParticles; i++)
    {
        if (m_model->getParticleState(i) == ParticleState::ShapeMatching)
        {
            Vector3r &pos = m_model->getPosition(i);
            Vector3r &vel = m_model->getVelocity(i);

            pos = c + R * radius_vector[i];;
            vel = (pos - oldPosition[i]) / dt; 
        }
    }
}


inline Matrix3r polarDecompose(Matrix3r A_pq)
{
    Matrix3r R,S;
    R.setZero();

    // S = (A_pq.transpose() * A_pq).cwiseSqrt();
    // R = A_pq * S.inverse();

    // Quaternionr q(R);
    // MathFunctions::extractRotation(A_pq, q, 10);
    // R = q.matrix();

    Vector3r sigma;
	Matrix3r U, VT;
	MathFunctions::svdWithInversionHandling(A_pq, sigma, U, VT);
	R = U * VT;

    return R;
}


void ShapeMatching::computeBarycenter()
{
    const  int numParticles = (int) m_model->numActiveParticles();
    if (numParticles == 0)
		return;

    Vector3r sum_mass_pos = Vector3r(0.0, 0.0, 0.0);
    Real sum_mass = 0.0;

    for (int i = 0; i < numParticles; i++)
    { 
        const Vector3r &xi = m_model->getPosition(i);
        const Real mass = m_model->getMass(i);
        sum_mass += mass ;
        sum_mass_pos += mass * xi;
    }
    total_mass = sum_mass;
    barycenter = sum_mass_pos / sum_mass;
    
    //calculate  radius vector and A_qq
    Matrix3r q;
    q.setZero();
    for (int i = 0; i < numParticles; i++)
    { 
        const Vector3r &xi = m_model->getPosition(i);
        radius_vector[i] = xi - barycenter;

        q += radius_vector[i]*(radius_vector[i]).transpose();
    }
    A_qq = q.inverse();

}