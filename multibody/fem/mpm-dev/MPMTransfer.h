#pragma once

#include <array>
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
    MPMTransfer(const Grid& grid, const Particles& particles);

    // Put it here until final refactor.
    // Sort the particles according to the batch number. As below shown,
    // o denotes the grid points, $ denotes the batch centered around the grid
    // point. # of batch = # of grid points
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

    // Given the particles information, accumulate vectors bases_val_particles
    // and bases_grad_particles.
    void UpdateBasisAndGradientParticles(const Particles& particles);

    // TODO(yiminlin.tri): temporary routine for testing, to be deleted in
    // future PRs
    const std::vector<std::array<double, 27>>> get_bases_val_particles() const;
    const std::vector<std::array<Vector3<double>, 27>>>
                                              get_bases_grad_particles() const;
 
 private:
    std::vector<BSpline> bases_;
    // The starting particle index of every batch. For example,
    // batch_starting_index_ = [0, 5, 6, 10] implies batches 0, 1, 2, 3
    // contains 5, 1, 4, #particles-10 particles. We assume SortParticles()
    // are called.
    std::vector<int> batch_starting_index_;
    // Evaluations and gradient evaluations of BSpline bases on each particle
    // i.e. N_i(x_p), \nabla N_i(x_p)
    // Length of the vector = # of particles.
    // Length of an element in the vector = 27 (max # of affected grid nodes)
    std::vector<std::array<double, 27>>> bases_val_particles_;
    std::vector<std::array<Vector3<double>, 27>>> bases_grad_particles_;
};  // class MPMTransfer

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
