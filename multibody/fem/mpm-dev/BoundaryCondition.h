#pragma once

#include <utility>
#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/geometry/proximity/posed_half_space.h"
#include "drake/math/rigid_transform.h"
#include "drake/math/rotation_matrix.h"
#include "drake/multibody/fem/mpm-dev/AnalyticLevelSet.h"

namespace drake {
namespace multibody {
namespace mpm {

// A Boundary condition class holding the boundary conditions' information in
// the domain. The boundaries in the MPM solver are represented as
// PosedHalfSpace. Points with nonpositive distances to the boundary half
// planes are inside the wall. We store the friction coefficients of each
// boundaries.
//            | normal direction             (outside of halfspace, dist > 0)
// ================================================================ wall
//                                           (inside of halfspace, dist < 0)
class BoundaryCondition {
 public:
    struct WallBoundary {
        double friction_coefficient;
        geometry::internal::PosedHalfSpace<double> boundary_space;
    };

    // A moving cylindrical boundary. The geometry in the physical frame is
    // described through the cylinder level set and pose. We assume the cylinder
    // is moving with constant velocity in the physical (canonical) frame.
    struct MovingCylindricalBoundary {
        double friction_coefficient;
        CylinderLevelSet level_set;
        math::RigidTransform<double> pose;
        Vector3<double> velocity;
    };

    BoundaryCondition() = default;
    explicit BoundaryCondition(std::vector<WallBoundary> wall_boundaries,
                               std::vector<MovingCylindricalBoundary>
                                                moving_cylindrical_boundaries);

    int get_num_wall_boundaries() const;
    const std::vector<WallBoundary>& get_wall_boundaries() const;
    const WallBoundary& get_wall_boundary(int index) const;

    void AddWallBoundary(WallBoundary boundary);

    int get_num_moving_cylindrical_boundaries() const;
    const std::vector<MovingCylindricalBoundary>&
        get_moving_cylindrical_boundaries() const;
    const MovingCylindricalBoundary&
        get_moving_cylindrical_boundary(int index) const;

    void AddMovingCylindricalBoundary(MovingCylindricalBoundary boundary);

    // Applies the frictional wall and moving cylindrical boundary conditions to
    // the grid point with the given position at time t, where we apply the
    // Coulumb friction law. Given the velocity v, and denotes its normal and
    // tangential components by vₙ = (v ⋅ n)n, vₜ = v - vₙ, the impulse of
    // colliding with the wall is given by j = -m vₙ. The Coulumb friction law
    // states the amount of friction imposed is at most μ ‖j‖, where \mu is the
    // friction coefficient of the wall.
    // If ‖vₜ‖ <= μ‖vₙ‖,   v_new = 0.0
    // Otherwise    ,   v_new = vₜ - μ‖vₙ‖t, t - tangential direction
    // Then we overwrite the passed in velocity with v_new
    // For a grid point locates on multiple boundaries, we impose each boundary
    // condition, ordered as in boundaries_, to this grid point.
    // TODO(yiminlin.tri): May cause unexpected behavior at sharp corners
    void Apply(double t, const Vector3<double>& position,
                         Vector3<double>* velocity) const;

 private:
    // Given the normal vector and friction coefficient mu, update the input
    // velocity using Coulumb friction law
    void UpdateVelocityCoulumbFriction(const Vector3<double>& normal, double mu,
                                       Vector3<double>* velocity) const;

    std::vector<WallBoundary> wall_boundaries_;
    std::vector<MovingCylindricalBoundary> moving_cylindrical_boundaries_;
};  // class BoundaryCondition

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
