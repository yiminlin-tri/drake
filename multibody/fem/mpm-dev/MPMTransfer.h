#pragma once

#include <numeric>
#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/mpm-dev/BSpline.h"
#include "drake/multibody/fem/mpm-dev/Grid.h"
#include "drake/multibody/fem/mpm-dev/Particles.h"

namespace drake {
namespace multibody {
namespace mpm {

// A implementation of MPM's Particles to Grid (P2G) and Grid to Particles (G2P)
// operations
class MPMTransfer {
 public:
    explicit MPMTransfer(const Grid& grid, const Particles& particles);

    // TODO(yiminlin.tri): Put it here until final refactor. Currently it is
    // user's responsibility to call SortParticles before every iterations.
    // Sort the particles according to the batch number, in increasing order.
    // As below shown, o denotes the grid points, $ denotes the batch centered
    // around the grid point. # of batch = # of grid points
    // o = = = o = = = o = = = o = = = o
    //||      ||      ||      ||      ||
    //||      ||      ||      ||      ||
    //||      ||      ||      ||      ||
    // o = = = o = = = o = = = o = = = o
    //||      ||      ||      ||      ||
    //||      ||      ||      ||      ||
    //||      ||      ||      ||      ||
    // o = = = o = = = o = = = o = = = o
    //||      ||      ||      ||      ||
    //||      ||   $$$$$$$$   ||      ||
    //||      ||   $  ||  $   ||      ||
    // o = = = o = $ = o =$= = o = = = o
    //||      ||   $  ||  $   ||      ||
    //||      ||   $$$$$$$$   ||      ||
    //||      ||      ||      ||      ||
    // o = = = o = = = o = = = o = = = o
    //||      ||      ||      ||      ||
    //||      ||      ||      ||      ||
    //||      ||      ||      ||      ||
    // o = = = o = = = o = = = o = = = o
    //||      ||      ||      ||      ||
    //||      ||      ||      ||      ||
    //||      ||      ||      ||      ||
    // o = = = o = = = o = = = o = = = o
    // The batches are ordered in a lexiographical ordering, similar to grid
    // points.
    void SortParticles(const Grid& grid, Particles* particles);

    int get_batch_starting_index(int batch_number) const;

    int get_num_particles_in_batch(int batch_number) const;

 private:
    int num_batches_;
    int num_particles_;
    std::vector<BSpline> bases_{};
    // The starting particle index of every batch. For example,
    // batch_starting_index_ = [0, 5, 6, 10] implies batches 0, 1, 2, 3
    // contains 5, 1, 4, #particles-10 particles.
    // The starting indices of empty batches are the starting index of the next
    // non-empty batch. As another example,
    // batch_starting_index_ = [0, 5, 8, 8, 9] implies batches 0, 1, 2, 3, 4
    // contains 5, 3, 0, 1, #particles
    // If the starting index of a batch equals to total number of particles,
    // then the batch has zero elements. For example, if we have 5 total
    // particles, then given the batch_starting_index_ = [0, 2, 4, 5, 5],
    // Each batch holds 2, 2, 1, 0, 0 particles
    // batch_starting_index_ is updated every time we called SortParticles()
    std::vector<int> batch_starting_index_{};
};  // class MPMTransfer

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
