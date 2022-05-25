#include "drake/multibody/fem/mpm-dev/MPMTransfer.h"

namespace drake {
namespace multibody {
namespace mpm {

// TODO(yiminlin.tri): We assume particles all lies in the grid
// Also may be a good idea to remove this routine's dependency on the grid,
// this need future refactoring
void MPMTransfer::SortParticles(const Grid& grid, Particles* particles) {
    int p_new;
    int idx_batch_prev = 0;
    int idx_batch_curr = 0;
    int num_batches = grid.get_num_gridpt();
    double h = grid.get_h();
    int num_particles = particles->get_num_particles();
    std::vector<size_t> sorted_indices(num_particles);
    Vector3<double> xp;
    std::vector<Vector3<int>> batch_indices(num_particles);

    batch_starting_index_.reserve(num_batches);

    // Preallocate the indices of batches
    for (int p = 0; p < num_particles; ++p) {
        batch_indices[p] = Vector3<int>(
                                std::round(particles->get_position_x(p)/h),
                                std::round(particles->get_position_y(p)/h),
                                std::round(particles->get_position_z(p)/h));
    }

    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::stable_sort(sorted_indices.begin(), sorted_indices.end(),
       [&grid, &batch_indices](size_t i1, size_t i2) {
      return grid.Reduce3DIndex(batch_indices[i1](0), batch_indices[i1](1),
                                batch_indices[i1](2))
           < grid.Reduce3DIndex(batch_indices[i2](0), batch_indices[i2](1),
                                batch_indices[i2](2));});

    // Reorder the particles
    particles->Reorder(sorted_indices);

    // Construct starting indices of batches
    batch_starting_index_[0] = 0;
    for (int p = 0; p < num_particles; ++p) {
        p_new = sorted_indices[p];
        idx_batch_prev = idx_batch_curr;
        xp = particles->get_position(p);
        idx_batch_curr = grid.Reduce3DIndex(batch_indices[p_new](0),
                                            batch_indices[p_new](1),
                                            batch_indices[p_new](2));
        if (idx_batch_curr != idx_batch_prev) {
            for (int b = idx_batch_curr; b > idx_batch_prev; --b) {
                batch_starting_index_[b] = p;
            }
        }
    }
    if (idx_batch_curr + 1 < num_batches - 1) {
        for (int b = idx_batch_curr + 1; b < num_batches; ++b) {
            batch_starting_index_[b] = num_particles;
        }
    }
}

int MPMTransfer::get_batch_starting_index(int batch_number) const {
    return batch_starting_index_[batch_number];
}

int MPMTransfer::get_num_particles_in_batch(int batch_number, int num_batches,
                                            int num_particles)
                                                                        const {
    return (batch_number == num_batches - 1) ?
           (num_particles - batch_starting_index_[batch_number]) :
           (batch_starting_index_[batch_number+1]
          - batch_starting_index_[batch_number]);
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
