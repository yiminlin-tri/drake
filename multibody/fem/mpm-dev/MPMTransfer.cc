#include "drake/multibody/fem/mpm-dev/MPMTransfer.h"

namespace drake {
namespace multibody {
namespace mpm {

void MPMTransfer::UpdateBasisAndGradientParticles(const Grid& grid,
                                                  const Particles& particles) {
    int batch_number, pstart, pend, idx_basis, idx_local;
    int h = grid.get_h();
    Vector3<int> idx_batch, bottom_corner, num_gridpt_1D;
    bottom_corner = grid.get_bottom_corner();
    num_gridpt_1D = grid.get_num_gridpt_1D();

    // O(27*Np), (the outer loop on batches is just for convinience)
    // Loop through all batches, (bi, bj, bk) is the index of the batch's
    // corresponding grid point
    for (int bk = bottom_corner(2); bk < bottom_corner(2) + num_gridpt_1D(2);
                                                                        ++bk) {
    for (int bj = bottom_corner(1); bj < bottom_corner(1) + num_gridpt_1D(1);
                                                                        ++bj) {
    for (int bi = bottom_corner(0); bi < bottom_corner(0) + num_gridpt_1D(0);
                                                                        ++bi) {
        batch_number = grid.Reduce3DIndex(bi, bj, bk);
        pstart = get_batch_starting_index(batch_number);
        pend = pstart + get_num_particles_in_batch(batch_number);
        // Loop through all particles inside the batch
        for (p = pstart; p < pend; ++p) {
            // Loop through all BSplines in the grid domain that contains the
            // particle
            for (int k = bk - 1; k <= bk + 1; ++k) {
            for (int j = bj - 1; j <= bj + 1; ++j) {
            for (int i = bi - 1; i <= bi + 1; ++i) {
                if (grid.in_index_range(i, j, k)) {
                    idx_basis = grid.Reduce3DIndex(i, j, k);
                    idx_local = i + 3*j + 9*k;
                    bases_val_particles_[p][idx_local] = bases_[idx_basis]
                    bases_grad_particles_[p][idx_local] = bases_[idx_basis]
                }
            }
            }
            }
        }
    }
    }
    }
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
