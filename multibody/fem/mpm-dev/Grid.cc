#include "drake/multibody/fem/mpm-dev/Grid.h"

namespace drake {
namespace multibody {
namespace mpm {

Grid::Grid(const Vector3<int>& num_gridpt_1D, double h,
           const Vector3<int>& bottom_corner):
           num_gridpt_(num_gridpt_1D(0)*num_gridpt_1D(1)*num_gridpt_1D(2)),
           num_gridpt_1D_(num_gridpt_1D), h_(h),
           bottom_corner_(bottom_corner) {
    int idx;
    DRAKE_ASSERT(num_gridpt_1D_(0) >= 0);
    DRAKE_ASSERT(num_gridpt_1D_(1) >= 0);
    DRAKE_ASSERT(num_gridpt_1D_(2) >= 0);
    DRAKE_ASSERT(h_ > 0.0);

    indices_ = std::vector<std::pair<int, Vector3<int>>>(num_gridpt_);
    positions_ = std::vector<Vector3<double>>(num_gridpt_);
    velocities_ = std::vector<Vector3<double>>(num_gridpt_);
    masses_ = std::vector<double>(num_gridpt_);
    forces_ = std::vector<Vector3<double>>(num_gridpt_);

    // Initialize the positions of grid points
    for (int k = bottom_corner_(2);
             k < bottom_corner_(2) + num_gridpt_1D_(2); ++k) {
    for (int j = bottom_corner_(1);
             j < bottom_corner_(1) + num_gridpt_1D_(1); ++j) {
    for (int i = bottom_corner_(0);
             i < bottom_corner_(0) + num_gridpt_1D_(0); ++i) {
        idx = Reduce3DIndex(i, j, k);
        indices_[idx] = std::pair<int, Vector3<int>>(idx, {i, j, k});
        positions_[idx] = Vector3<double>{h_*i, h_*j, h_*k};
    }
    }
    }
}

int Grid::get_num_gridpt() const {
    return num_gridpt_;
}

const Vector3<int>& Grid::get_num_gridpt_1D() const {
    return num_gridpt_1D_;
}

double Grid::get_h() const {
    return h_;
}

const Vector3<int>& Grid::get_bottom_corner() const {
    return bottom_corner_;
}

const Vector3<double>& Grid::get_position(int i, int j, int k) const {
    DRAKE_ASSERT(in_index_range(i, j, k));
    return positions_[Reduce3DIndex(i, j, k)];
}

const Vector3<double>& Grid::get_velocity(int i, int j, int k) const {
    DRAKE_ASSERT(in_index_range(i, j, k));
    return velocities_[Reduce3DIndex(i, j, k)];
}

double Grid::get_mass(int i, int j, int k) const {
    DRAKE_ASSERT(in_index_range(i, j, k));
    return masses_[Reduce3DIndex(i, j, k)];
}

const Vector3<double>& Grid::get_force(int i, int j, int k) const {
    DRAKE_ASSERT(in_index_range(i, j, k));
    return forces_[Reduce3DIndex(i, j, k)];
}

void Grid::set_position(int i, int j, int k, const Vector3<double>& position) {
    DRAKE_ASSERT(in_index_range(i, j, k));
    positions_[Reduce3DIndex(i, j, k)] = position;
}

void Grid::set_velocity(int i, int j, int k, const Vector3<double>& velocity) {
    DRAKE_ASSERT(in_index_range(i, j, k));
    velocities_[Reduce3DIndex(i, j, k)] = velocity;
}

void Grid::set_mass(int i, int j, int k, double mass) {
    DRAKE_ASSERT(in_index_range(i, j, k));
    masses_[Reduce3DIndex(i, j, k)] = mass;
}

void Grid::set_force(int i, int j, int k, const Vector3<double>& force) {
    DRAKE_ASSERT(in_index_range(i, j, k));
    forces_[Reduce3DIndex(i, j, k)] = force;
}

void Grid::AccumulateVelocity(int i, int j, int k,
                                        const Vector3<double>& velocity) {
    DRAKE_ASSERT(in_index_range(i, j, k));
    velocities_[Reduce3DIndex(i, j, k)] += velocity;
}

void Grid::AccumulateMass(int i, int j, int k, double mass) {
    DRAKE_ASSERT(in_index_range(i, j, k));
    masses_[Reduce3DIndex(i, j, k)] += mass;
}

void Grid::AccumulateForce(int i, int j, int k, const Vector3<double>& force) {
    DRAKE_ASSERT(in_index_range(i, j, k));
    forces_[Reduce3DIndex(i, j, k)] += force;
}

void Grid::AccumulateVelocity(const Vector3<int>& index_3d,
                                        const Vector3<double>& velocity) {
    AccumulateVelocity(index_3d(0), index_3d(1), index_3d(2), velocity);
}

void Grid::AccumulateMass(const Vector3<int>& index_3d, double mass) {
    AccumulateMass(index_3d(0), index_3d(1), index_3d(2), mass);
}

void Grid::AccumulateForce(const Vector3<int>& index_3d,
                           const Vector3<double>& force) {
    AccumulateForce(index_3d(0), index_3d(1), index_3d(2), force);
}

void Grid::RescaleVelocities() {
    for (int i = 0; i < num_gridpt_; ++i) {
        velocities_[i] = velocities_[i] / masses_[i];
    }
}

void Grid::ResetStates() {
    std::fill(masses_.begin(), masses_.end(), 0.0);
    std::fill(velocities_.begin(), velocities_.end(), Vector3<double>::Zero());
    std::fill(forces_.begin(), forces_.end(), Vector3<double>::Zero());
}

int Grid::Reduce3DIndex(int i, int j, int k) const {
    DRAKE_ASSERT(in_index_range(i, j, k));
    return (k-bottom_corner_(2))*(num_gridpt_1D_(0)*num_gridpt_1D_(1))
         + (j-bottom_corner_(1))*num_gridpt_1D_(0)
         + (i-bottom_corner_(0));
}

int Grid::Reduce3DIndex(const Vector3<int>& index_3d) const {
    return Reduce3DIndex(index_3d(0), index_3d(1), index_3d(2));
}

Vector3<int> Grid::Expand1DIndex(int idx) const {
    return Vector3<int>(
            bottom_corner_(0) + idx % num_gridpt_1D_(0),
            bottom_corner_(1) + (idx / num_gridpt_1D_(0)) % num_gridpt_1D_(1),
            bottom_corner_(2) + idx / (num_gridpt_1D_(0)*num_gridpt_1D_(1)));
}

const std::vector<std::pair<int, Vector3<int>>>& Grid::get_indices() const {
    return indices_;
}

bool Grid::in_index_range(int i, int j, int k) const {
    return ((i < bottom_corner_(0) + num_gridpt_1D_(0)) &&
            (j < bottom_corner_(1) + num_gridpt_1D_(1)) &&
            (k < bottom_corner_(2) + num_gridpt_1D_(2)) &&
            (i >= bottom_corner_(0)) &&
            (j >= bottom_corner_(1)) &&
            (k >= bottom_corner_(2)));
}

bool Grid::in_index_range(const Vector3<int>& index_3d) const {
    return in_index_range(index_3d(0), index_3d(1), index_3d(2));
}

void Grid::UpdateVelocity(double dt) {
    for (int i = 0; i < num_gridpt_; ++i) {
        velocities_[i] += dt*forces_[i]/masses_[i];
    }
}

void Grid::EnforceWallBoundaryCondition(double friction_coefficient,
            const geometry::internal::PosedHalfSpace<double>& boundary_space) {
    // For all grid points
    for (const auto& [index_flat, index_3d] : indices_) {
        // If the grid point is inside the boundary space and its distance to
        // the half plane is at most 2h, enforce BC. If the grid point is
        // outside the boundary space, set the velocity to be 0.0.
        double dist = boundary_space.CalcSignedDistance(positions_[index_flat]);
        if (dist > 0) {
            velocities_[index_flat] = Vector3<double>::Zero();
        } else if (dist >= -2.0*h_) {
            const Vector3<double>& normal = boundary_space.normal();
            const Vector3<double>& v_i = velocities_[index_flat];
            double vn = v_i.dot(normal);           // v \dot normal
            Vector3<double> vt = v_i - vn*normal;  // Tangential velocity
            velocities_[index_flat] = vt
                                    + friction_coefficient*vn*vt.normalized();
        }
    }
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake