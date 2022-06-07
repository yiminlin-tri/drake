#pragma once

#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/geometry/proximity/posed_half_space.h"

namespace drake {
namespace multibody {
namespace mpm {

// A Boundary condition class holding the boundary conditions' information in
// the domain. The boundaries in the MPM solver are represented as
// PosedHalfSpace, whose normals point from the interior to the exterior of our
// computational domain. Points with nonnegative distances to the boundary half
// planes are outside the computational domain. We store the friction
// coefficients of each boundaries.
// TODO(yiminlin.tri): Only support frictional wall boundary condition
class BoundaryCondition {
 public:
    struct Boundary {
        double friction_coefficient;
        geometry::internal::PosedHalfSpace<double> boundary_space;
    };

    BoundaryCondition() = default;
    explicit BoundaryCondition(const std::vector<Boundary>& boundaries);

    int get_num_boundary() const;
    const std::vector<Boundary>& get_boundaries() const;
    const Boundary& get_boundary(int index) const;

    void AddBoundary(const Boundary& boundary);

    // Applies the frictional wall boundary condition to the grid point with the
    // given position, where we apply the Coulumb friction law. Given the
    // velocity v, and denotes its normal and tangential components by vn =
    // (v\dot n)n, vt = v - vn, the impulse of colliding with the wall is given
    // by j = -m vn. The Coulumb friction law states the amount of friction
    // imposed is at most \mu \|j\|, where \mu is the friction coefficient of
    // the wall.
    // If \|vt\| <= mu \|j\|/m, v_new = 0.0
    // Otherwise              , v_new = vt - mu\|j\|/m * t, t - tangential
    //                                                          direction
    // Note that we use \|j\|/m = \|vn\| in the implementation
    // Then we overwrite the passed in velocity with v_new
    void UpdateFrictionalWallVelocity(const Vector3<double>& position,
                                      Vector3<double>* velocity) const;

 private:
    std::vector<Boundary> boundaries_;
};  // class BoundaryCondition

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
