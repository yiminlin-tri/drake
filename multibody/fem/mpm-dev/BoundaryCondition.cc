#include "drake/multibody/fem/mpm-dev/BoundaryCondition.h"

namespace drake {
namespace multibody {
namespace mpm {

BoundaryCondition::BoundaryCondition(
                        std::vector<WallBoundary> wall_boundaries,
                        std::vector<MovingCylindricalBoundary>
                                                moving_cylindrical_boundaries):
                wall_boundaries_(std::move(wall_boundaries)),
                moving_cylindrical_boundaries_(
                                 std::move(moving_cylindrical_boundaries)) {}

int BoundaryCondition::get_num_wall_boundaries() const {
    return wall_boundaries_.size();
}

const std::vector<BoundaryCondition::WallBoundary>&
                            BoundaryCondition::get_wall_boundaries() const {
    return wall_boundaries_;
}

const BoundaryCondition::WallBoundary&
                        BoundaryCondition::get_wall_boundary(int index) const {
    DRAKE_ASSERT(index < wall_boundaries_.size());
    return wall_boundaries_[index];
}

void BoundaryCondition::AddWallBoundary(BoundaryCondition::WallBoundary
                                                                    boundary) {
    wall_boundaries_.emplace_back(std::move(boundary));
}

int BoundaryCondition::get_num_moving_cylindrical_boundaries() const {
    return moving_cylindrical_boundaries_.size();
}

const std::vector<BoundaryCondition::MovingCylindricalBoundary>&
                BoundaryCondition::get_moving_cylindrical_boundaries() const {
    return moving_cylindrical_boundaries_;
}

const BoundaryCondition::MovingCylindricalBoundary&
        BoundaryCondition::get_moving_cylindrical_boundary(int index) const {
    DRAKE_ASSERT(index < moving_cylindrical_boundaries_.size());
    return moving_cylindrical_boundaries_[index];
}

void BoundaryCondition::AddMovingCylindricalBoundary(
            BoundaryCondition::MovingCylindricalBoundary boundary) {
    moving_cylindrical_boundaries_.emplace_back(std::move(boundary));
}

void BoundaryCondition::Apply(double t, const Vector3<double>& position,
                                        Vector3<double>* velocity) const {
    // For all wall boundaries
    for (const auto& wall_boundary : wall_boundaries_) {
        // If the grid point is on/in the boundary, enforce the frictional wall
        // boundary condition
        if (wall_boundary.boundary_space.CalcSignedDistance(position) <= 0) {
            // Normal vector of the boundary plane
            const Vector3<double>& n_i = wall_boundary.boundary_space.normal();
            // Normal and tangential components of the velocity
            Vector3<double> vn = velocity->dot(n_i)*n_i;
            Vector3<double> vt = *velocity - vn;
            // Normal and tangential speed
            double vn_norm = vn.norm();
            double vt_norm = vt.norm();

            // Limit the amount of friction to at most eliminating the
            // tangential velocity
            if (vt_norm <= wall_boundary.friction_coefficient*vn_norm) {
                *velocity = Vector3<double>::Zero();
            } else {
                // If the tangential velocity is zero, the updated velocity is
                // zero by above.
                *velocity = vt - wall_boundary.friction_coefficient
                                *vn_norm*vt/vt_norm;
            }
        }
    }

    // For all moving cylindrical boundaries
    for (const auto& moving_cyl_boundary : moving_cylindrical_boundaries_) {
        // particle's position in cylinder's reference frame
        Vector3<double> x_cyl = moving_cyl_boundary.pose.inverse()
                               *(position - t*moving_cyl_boundary.velocity);
        if (moving_cyl_boundary.level_set.InInterior(x_cyl)) {
            Vector3<double> n_i = moving_cyl_boundary.level_set.Normal(x_cyl);
            // particle's velocity in cylinder's reference frame
            Vector3<double> v_ref =
                            moving_cyl_boundary.pose.rotation().inverse()
                           *(*velocity-moving_cyl_boundary.velocity);
            // Normal and tangential components of the reference frame velocity
            Vector3<double> vn = v_ref.dot(n_i)*n_i;
            Vector3<double> vt = v_ref - vn;
            // Normal and tangential speed
            double vn_norm = vn.norm();
            double vt_norm = vt.norm();
            // Limit the amount of friction to at most eliminating the
            // tangential velocity
            Vector3<double> v_ref_new;
            if (vt_norm <=
                moving_cyl_boundary.friction_coefficient*vn_norm) {
                v_ref_new = Vector3<double>::Zero();
            } else {
                // If the tangential velocity is zero, the updated velocity is
                // zero by above.
                v_ref_new = vt -
                            moving_cyl_boundary.friction_coefficient
                           *vn_norm*vt/vt_norm;
            }
            // Transform the reference velocity back to the physical frame
            *velocity = moving_cyl_boundary.pose.rotation()*v_ref_new
                      + moving_cyl_boundary.velocity;
        }
    }
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
