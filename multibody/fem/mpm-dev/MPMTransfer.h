#pragma once

#include <algorithm>
#include <array>
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
    MPMTransfer() {}

 private:
    friend class MPMTransferTest;
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

    // Update the evalutions and gradients of BSpline bases on each particle,
    // and update bases_val_particles_ and bases_grad_particles_
    void UpdateBasisAndGradientParticles(const Grid& grid,
                                         const Particles& particles);

    // Given the position of a particle xp, calculate the index of the batch
    // this particle is in.
    Vector3<int> CalcBatchIndex(const Vector3<double>& xp, double h) const;

    // Evaluations and gradients of BSpline bases on each particle
    // i.e. N_i(x_p), \nabla N_i(x_p)
    // Length of the vector = # of particles.
    // Length of an element in the vector = 27 (max # of affected grid nodes)
    std::vector<std::array<double, 27>> bases_val_particles_{};
    std::vector<std::array<Vector3<double>, 27>> bases_grad_particles_{};
    // A vector holding the number of particles inside each batch
    std::vector<int> batch_sizes_{};
};  // class MPMTransfer

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
