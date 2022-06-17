#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/math/rigid_transform.h"
#include "drake/multibody/fem/mpm-dev/AnalyticLevelSet.h"

namespace drake {
namespace multibody {
namespace mpm {

// A class holding all collision objects' information in a MPM simulation.
class KinematicCollisionObjects {
 public:
    KinematicCollisionObjects() = default;

    // Add a new collision object with the given intial conditions
    void AddCollisionObject(std::unique_ptr<AnalyticLevelSet> level_set,
                            math::RigidTransform<double> pose,
                            math::SpatialVelocity<double> spatial_velocity,
                            double friction_coeff);

    // Advance the state of all collision objects by dt
    void AdvanceOneTimeStep(double dt);

    // Apply the boundary condition prescribed by all collision objects to a
    // point in the space with the given postion and velocity.
    void ApplyBoundaryConditions(const Vector3<double>& position,
                                 Vector3<double>* velocity) const;

 private:
    // A private class represents collision objects in MPM simulations. The
    // collision objects are represented as analytical level sets (instead of
    // particles, so they are not "MPM objects") with initial pose, spatial
    // velocity (translation and angular velocities) and friction coefficient on
    // the boundary. These collision objects are placed in a MPM simulation to
    // simulate the environment that are not coupled with "MPM objects", for
    // example, to enforce the boundary conditions prescribed by this object.
    // We assume the collision object is not deformable.
    class CollisionObject {
     public:
        struct CollisionObjectState {
            math::RigidTransform<double> pose;
            math::SpatialVelocity<double> spatial_velocity;
        };

        CollisionObject(std::unique_ptr<AnalyticLevelSet> level_set,
                        CollisionObjectState initial_state,
                        double friction_coeff);

        // Advance the state (pose) of the collision object by dt
        void AdvanceOneTimeStep(double dt);

        // Apply the boundary condition prescribed by this object to a point in
        // the space with the given postion and velocity. On Boundaries,
        // we apply the Coulumb friction law. Given the velocity v, and denotes
        // its normal and tangential components by vₙ = (v ⋅ n)n, vₜ = v - vₙ,
        // the impulse of colliding with the wall is given by j = -m vₙ. The
        // Coulumb friction law states the amount of friction imposed is at most
        // μ ‖j‖, where \mu is the friction coefficient of the wall.
        // If ‖vₜ‖ <= μ‖vₙ‖,   v_new = 0.0
        // Otherwise    ,   v_new = vₜ - μ‖vₙ‖t, t - tangential direction
        // Then we overwrite the passed in velocity with v_new.
        // For a grid point locates on multiple boundaries, we impose each
        // boundary condition, ordered as in boundaries_, to this grid point. 
        // TODO(yiminlin.tri): May cause unexpected behavior at sharp corners
        void ApplyBoundaryCondition(const Vector3<double>& position,
                                    const Vector3<double>* velocity) const;

     private:
        // Given the normal vector and friction coefficient mu, update the input
        // velocity using Coulumb friction law
        void UpdateVelocityCoulumbFriction(const Vector3<double>& normal,
                                           double mu,
                                           Vector3<double>* velocity) const;

        CollisionObjectState state_;
        std::unique_ptr<AnalyticLevelSet> level_set_;
        double friction_coeff_;
    }

    std::vector<std::unique_ptr<CollisionObject>> collision_objects_;
};  // class BoundaryCondition

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
