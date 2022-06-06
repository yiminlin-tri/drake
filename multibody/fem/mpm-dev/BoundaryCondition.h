#pragma once

#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/geometry/proximity/posed_half_space.h"

namespace drake {
namespace multibody {
namespace mpm {

// A Boundary condition class holding the boundary conditions' information
// in the domain. The boundaries in the MPM solver are represented as
// PosedHalfSpace, whose normals point outward from the interior to the exterior
// of our computational domain. And we store the friction coefficients of each
// boundaries.
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

 private:
    std::vector<Boundary> boundaries_;
};  // class BoundaryCondition

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
