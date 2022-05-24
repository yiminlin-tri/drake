#include "drake/multibody/fem/mpm-dev/MPMTransfer.h"

namespace drake {
namespace multibody {
namespace mpm {

MPMTransfer::MPMTransfer(const Grid& grid, const Particles& particles) {
    int idx;
    int num_gridpt = grid.get_num_gridpt();
    Vector3<int> num_gridpt_1D = grid.get_num_gridpt_1D();
    double h = grid.get_h();
    Vector3<int> bottom_corner = grid.get_bottom_corner();
    Vector3<double> pos;
    num_batches_ = num_gridpt;
    num_particles_ = particles.get_num_particles();
    batch_starting_index_.reserve(num_gridpt);
    bases_.reserve(num_gridpt);

    // Initialize the vector of basis based on grids
    for (int k = bottom_corner(2); k < bottom_corner(2)+num_gridpt_1D(2); ++k) {
    for (int j = bottom_corner(1); j < bottom_corner(1)+num_gridpt_1D(1); ++j) {
    for (int i = bottom_corner(0); i < bottom_corner(0)+num_gridpt_1D(0); ++i) {
        idx = grid.Reduce3DIndex(i, j, k);
        bases_[idx] = BSpline(h, grid.get_position(i, j, k));
    }
    }
    }
}

// TODO(yiminlin.tri): inplace sorting? We assume particles all lies in the grid
// Also may be a good idea to remove this routine's dependency on the grid,
// this need future refactoring
void MPMTransfer::SortParticles(const Grid& grid, Particles* particles) {
    int idx_batch_prev = 0;
    int idx_batch_curr = 0;
    double h = grid.get_h();
    int num_particles = particles->get_num_particles();
    std::vector<size_t> sorted_indices(num_particles);
    std::vector<Vector3<double>> positions_sorted(num_particles);
    std::vector<Vector3<double>> velocities_sorted(num_particles);
    std::vector<double> masses_sorted(num_particles);
    std::vector<double> reference_volumes_sorted(num_particles);
    std::vector<Matrix3<double>> deformation_gradients_sorted(num_particles);
    std::vector<Matrix3<double>> kirchhoff_stresses_sorted(num_particles);

    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::stable_sort(sorted_indices.begin(), sorted_indices.end(),
       [&grid, particles, h](size_t i1, size_t i2) {
      return grid.Reduce3DIndex(std::round(particles->get_position_x(i1)/h),
                                std::round(particles->get_position_y(i1)/h),
                                std::round(particles->get_position_z(i1)/h))
           < grid.Reduce3DIndex(std::round(particles->get_position_x(i2)/h),
                                std::round(particles->get_position_y(i2)/h),
                                std::round(particles->get_position_z(i2)/h));});

    batch_starting_index_[0] = 0;
    for (int p = 0; p < num_particles; ++p) {
        positions_sorted[p] =
                        particles->get_position(sorted_indices[p]);
        velocities_sorted[p] =
                        particles->get_velocity(sorted_indices[p]);
        masses_sorted[p] =
                        particles->get_mass(sorted_indices[p]);
        reference_volumes_sorted[p] =
                        particles->get_reference_volume(sorted_indices[p]);
        deformation_gradients_sorted[p] =
                        particles->get_deformation_gradient(sorted_indices[p]);
        kirchhoff_stresses_sorted[p] =
                        particles->get_kirchhoff_stress(sorted_indices[p]);
        idx_batch_prev = idx_batch_curr;
        idx_batch_curr = grid.Reduce3DIndex(
                                        std::round(positions_sorted[p](0)/h),
                                        std::round(positions_sorted[p](1)/h),
                                        std::round(positions_sorted[p](2)/h));
        if (idx_batch_curr != idx_batch_prev) {
            for (int b = idx_batch_curr; b > idx_batch_prev; --b) {
                batch_starting_index_[b] = p;
            }
        }
    }
    if (idx_batch_curr + 1 < num_batches_ - 1) {
        for (int b = idx_batch_curr + 1; b < num_batches_; ++b) {
            batch_starting_index_[b] = num_particles;
        }
    }

    particles->set_positions(positions_sorted);
    particles->set_velocities(velocities_sorted);
    particles->set_masses(masses_sorted);
    particles->set_reference_volumes(reference_volumes_sorted);
    particles->set_deformation_gradients(deformation_gradients_sorted);
    particles->set_kirchhoff_stresses(kirchhoff_stresses_sorted);
}

int MPMTransfer::get_batch_starting_index(int batch_number) const {
    return batch_starting_index_[batch_number];
}

int MPMTransfer::get_num_particles_in_batch(int batch_number) const {
    return (batch_number == num_batches_ - 1) ?
           (num_particles_ - batch_starting_index_[batch_number]) :
           (batch_starting_index_[batch_number+1]
          - batch_starting_index_[batch_number]);
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
