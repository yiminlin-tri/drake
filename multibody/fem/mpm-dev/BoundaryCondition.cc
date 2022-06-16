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
            UpdateVelocityCoulumbFriction(wall_boundary.boundary_space.normal(),
                                          wall_boundary.friction_coefficient,
                                          velocity);
        }
    }

    // For all moving cylindrical boundaries
    for (const auto& moving_cyl_boundary : moving_cylindrical_boundaries_) {
        // particle's position in cylinder's reference frame
        Vector3<double> x_cyl = moving_cyl_boundary.pose.inverse()
                               *(position - t*moving_cyl_boundary.velocity);
        if (moving_cyl_boundary.level_set.InInterior(x_cyl)) {
            // transform particle's velocity to cylinder's reference frame
            *velocity = moving_cyl_boundary.pose.rotation().inverse()
                       *(*velocity-moving_cyl_boundary.velocity);
            UpdateVelocityCoulumbFriction(
                                    moving_cyl_boundary.level_set.Normal(x_cyl),
                                    moving_cyl_boundary.friction_coefficient,
                                    velocity);
            // Transform the reference velocity back to the physical frame
            *velocity = moving_cyl_boundary.pose.rotation()*(*velocity)
                      + moving_cyl_boundary.velocity;
        }
    }
}

void BoundaryCondition::UpdateVelocityCoulumbFriction(const Vector3<double>& n,
                                            double mu,
                                            Vector3<double>* velocity) const {
    // Normal and tangential components of the velocity
    Vector3<double> vn = velocity->dot(n)*n;
    Vector3<double> vt = *velocity - vn;
    // Normal and tangential speed
    double vn_norm = vn.norm();
    double vt_norm = vt.norm();

    // Limit the amount of friction to at most eliminating the
    // tangential velocity
    if (vt_norm <= mu*vn_norm) {
        *velocity = Vector3<double>::Zero();
    } else {
        // If the tangential velocity is zero, the updated velocity is
        // zero by above.
        *velocity = vt - mu*vn_norm*vt/vt_norm;
    }
}
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
