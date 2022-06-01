#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <tuple>
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

    // Sort the particles according to the batch number, and preallocate
    // bases' evaluations for transfer routines. This routine has to be called
    // at the beginning of each time step (before P2G and G2P transfers),
    // otherwise the results will be incorrect.
    void SetUpTransfer(const Grid& grid, Particles* particles);

    // Transfer masses, velocities, and Kirchhoff stresses on the particles
    // to masses, velocities, and forces on the grid
    void TransferParticlesToGrid(const Particles& particles, Grid* grid);

    // Transfer velocities on the grids to velocities and deformation
    // gradients on the particles
    void TransferGridToParticles(const Grid& grid, double dt,
                                 Particles* particles);

 private:
    friend class MPMTransferTest;

    struct GridState {
        double mass;
        Vector3<double> velocity;
        Vector3<double> force;

        void reset() {
            mass = 0.0;
            velocity.setZero();
            force.setZero();
        }
    };

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
    // SortParticles assume particles are within the bound. We define a particle
    // is not in bound if the particle does not lie within the grid, or the
    // particle's position is in the boundary batches. For example, in above
    // case, the batches that correspond to boundary grid points are boundary
    // batches.
    void SortParticles(const Grid& grid, Particles* particles);

    // Update the evalutions and gradients of BSpline bases on each particle,
    // and update bases_val_particles_ and bases_grad_particles_
    void UpdateBasisAndGradientParticles(const Grid& grid,
                                         const Particles& particles);

    // Evaluate (27) bases neighboring to the given batch, at the p-th particle
    // with position xp, and put the results into preallocated vectors
    void EvalBasisOnBatch(int p, const Vector3<double>& xp, const Grid& grid,
                          const Vector3<int>& batch_index_3d,
                          const std::vector<BSpline>& bases);

    // At a particular particle p in batch with batch_index_3d, transfer
    // particle states (m, mv, tau) to (m, mv, f). Note that we temporarily
    // store the momentum into particles' velocities, in TransferParticlesToGrid
    // we will scale the momentum with the updated mass to get the velocities.
    void AccumulateGridStatesOnBatch(int p, double mass_p,
                                     double reference_volume_p,
                                     const Vector3<double>& momentum_p,
                                     const Matrix3<double>& tau_p,
                                     std::array<GridState, 27>* sum_local);

    void WriteBatchStateToGrid(const Vector3<int>& batch_index_3d,
                               const std::array<GridState, 27>& sum_local,
                               Grid* grid);

    // Update particle states F_p^{n+1} and v_p^{n+1}
    void WriteBatchStateToParticles(const std::array<Vector3<double>, 27>&
                                                            batch_velocities,
                                    double dt, int p,
                                    Particles* particles);

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
