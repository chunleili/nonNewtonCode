#ifndef __ShapeMatching_h__
#define __ShapeMatching_h__
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"
namespace SPH
{
/**
 * @brief this class is not using the PBD,
 * this class calculate the rigid body translation and rotation for particles directly
 * the merit of this way is that the particles can switch states between rigid body and SPH fluid seamlessly
 * which is useful for phase change simulation(melting and freezing)
 * 
 */
class ShapeMatching : NonPressureForceBase
{
private:
    Quaternionr quaternion{0.0, 0.0, 0.0, 0.0};
    Vector3r barycenter{0.0, 0.0, 0.0};
    Vector3r angular_velocity{0.0, 0.0, 0.0};
    std::vector<Vector3r> radius_vector;
    std::vector<Vector3r> oldPosition;
    Real total_mass{0.0};
    Matrix3r A_qq;


    void setStates();
    void computeBarycenter();
    void collision_response();
    void shapeMatching();

public:
    ShapeMatching(FluidModel *model);
    ~ShapeMatching();
    void step();
    void addForce();
    void addTorque();
    static NonPressureForceBase* creator(FluidModel* model) {return new ShapeMatching(model);}

    Vector3r velocity{0.0, 0.0, 0.0};
    Matrix3r rotationMatrix ;
};

}
#endif