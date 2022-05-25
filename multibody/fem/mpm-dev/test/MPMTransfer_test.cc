#include "drake/multibody/fem/mpm-dev/MPMTransfer.h"

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace mpm {
namespace internal {
namespace {

constexpr double kEps = 4.0 * std::numeric_limits<double>::epsilon();

GTEST_TEST(MPMTransferTest, MPMTransferPreallocateVal) {
    // Construct a grid of 5x5x5 on [-2,2]^3, and place 125 particles
    // on the grid points.
    //                 -2  -1  0   1   2
    //                 o - o - o - o - o
    //                 |   |   |   |   |
    //             o - o - o - o - o - o
    //             |   |   |   |   |   |
    //         o - o - o - o - o - o - o
    //         |   |   |   |   |   |   |
    //     o - o - o - o - o - o - o - o
    //     |   |   |   |   |   |   |   |
    // o - o - o - o - o - o - o - o - o
    // |   |   |   |   |   |   |   |   |
    // o - o - o - o - o - o - o - o - o
    // |   |   |   |   |   |   |   |
    // o - o - o - o - o - o - o - o
    // |   |   |   |   |   |   |
    // o - o - o - o - o - o - o
    // |   |   |   |   |
    // o - o - o - o - o
    int pc;
    std::vector<std::array<double, 27>>> bases_val_particles;
    std::vector<std::array<Vector3<double>, 27>>> bases_grad_particles;
    int h = 1.0;
    Vector3<int> num_gridpt_1D = { 5,  5,  5};
    Vector3<int> bottom_corner = {-2, -2, -2};
    Grid grid = Grid(num_gridpt_1D, h, bottom_corner);
    int num_particles = grid.get_num_gridpt();
    Particles particles = Particles(num_particles);
    MPMTransfer mpm_transfer = MPMTransfer(grid);

    // Set particles' positions to be on grid points
    pc = 0;
    for (int k = bottom_corner(2); k < bottom_corner(2)+num_grid_pt_1D(2);
                                                                        ++k) {
    for (int j = bottom_corner(1); j < bottom_corner(1)+num_grid_pt_1D(1);
                                                                        ++j) {
    for (int i = bottom_corner(0); i < bottom_corner(0)+num_grid_pt_1D(0);
                                                                        ++i) {
        particles.set_position(pc++, grid.get_position(i, j, k));
    }
    }
    }

    // Sort the particles and set up the batches
    mpm_transfer.SortParticles();

    // Preallocate basis evaluations
    mpm_transfer.UpdateBasisAndGradientParticles(particles);
    bases_val_particles  = mpm_transfer.get_bases_val_particles();
    bases_grad_particles = mpm_transfer.get_bases_grad_particles();

    for (int p = 0; p < num_particles; ++p) {

    }
}

}  // namespace
}  // namespace internal
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
