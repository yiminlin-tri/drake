#include "drake/multibody/fem/mpm-dev/MPMTransfer.h"

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace mpm {
namespace internal {
namespace {

constexpr double kEps = 4.0 * std::numeric_limits<double>::epsilon();

GTEST_TEST(MPMTransferTest, MPMTransferSortParticlesTest) {
    // First, construct a 2x2x2 grid centered at (0.0, 0.0, 0.0), with
    // h = 2:
    //         o - o - o
    //         |   |   |
    //     o - o - o - o
    //     |   |   |   |
    // o - o - o - o - o
    // |   |   |   |
    // o - o - o - o
    // |   |   |
    // o - o - o
    // And make 27 particles on the location of 27 grid points. The particles'
    // order are in a reversed lexiographical ordering.
    Vector3<int> num_grid_pt_1D, bottom_corner;
    double h;
    int num_grid_pt, pc;
    Grid grid;
    Particles particles;

    num_grid_pt_1D = Vector3<int>(3, 3, 3);
    h = 2.0;
    bottom_corner = Vector3<int>(-1, -1, -1);

    grid = Grid(num_grid_pt_1D, h, bottom_corner);
    particles = Particles(27);
    EXPECT_EQ(particles.get_num_particles(), 27);
    MPMTransfer mpm_transfer = MPMTransfer(grid, particles);
    num_grid_pt = grid.get_num_gridpt();

    EXPECT_EQ(num_grid_pt, 27);
    EXPECT_EQ(num_grid_pt_1D(0), 3);
    EXPECT_EQ(num_grid_pt_1D(1), 3);
    EXPECT_EQ(num_grid_pt_1D(2), 3);
    pc = num_grid_pt;
    for (int k = bottom_corner(2); k < bottom_corner(2)+num_grid_pt_1D(2);
                                                                        ++k) {
    for (int j = bottom_corner(1); j < bottom_corner(1)+num_grid_pt_1D(1);
                                                                        ++j) {
    for (int i = bottom_corner(0); i < bottom_corner(0)+num_grid_pt_1D(0);
                                                                        ++i) {
        particles.set_position(--pc, grid.get_position(i, j, k));
    }
    }
    }

    // Sanity check
    EXPECT_TRUE(CompareMatrices(particles.get_position(0),
                                Vector3<double>(2.0, 2.0, 2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(1),
                                Vector3<double>(0.0, 2.0, 2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(8),
                                Vector3<double>(-2.0, -2.0, 2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(13),
                                Vector3<double>(0.0, 0.0, 0.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(26),
                                Vector3<double>(-2.0, -2.0, -2.0),
                                std::numeric_limits<double>::epsilon()));

    // Check particles in the correct ordering after sorting
    mpm_transfer.SortParticles(grid, &particles);

    EXPECT_TRUE(CompareMatrices(particles.get_position(0),
                                Vector3<double>(-2.0, -2.0, -2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(1),
                                Vector3<double>(0.0, -2.0, -2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(8),
                                Vector3<double>(2.0, 2.0, -2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(13),
                                Vector3<double>(0.0, 0.0, 0.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(26),
                                Vector3<double>(2.0, 2.0, 2.0),
                                std::numeric_limits<double>::epsilon()));

    // Check batches have the correct starting indices
    for (int i = 0; i < num_grid_pt; ++i) {
        EXPECT_EQ(mpm_transfer.get_batch_starting_index(i), i);
        EXPECT_EQ(mpm_transfer.get_num_particles_in_batch(i), 1);
    }

    // Next, construct a new set of particles, which contains particles
    // at locations (-0.5, 0.5, -0.5), (0.5, 0.5, 0.5), (0.5, -0.5, 0.5)
    // So all particles are at the batch # 13,
    // The batch_starting_index_ shall look like:
    //            i=13
    // (0, 0, ..., 0, 3, 0, ...)
    // and the number of particles in the batch shall look like:
    // (0, 0, ..., 3, 0, 0, ...)
    // this mainly tests the functionality of start indices
    particles = Particles(3);
    mpm_transfer = MPMTransfer(grid, particles);
    num_grid_pt = grid.get_num_gridpt();
    particles.set_position(0, Vector3<double>(-0.5, 0.5, -0.5));
    particles.set_position(1, Vector3<double>(0.5, 0.5, 0.5));
    particles.set_position(2, Vector3<double>(0.5, -0.5, 0.5));
    mpm_transfer.SortParticles(grid, &particles);

    // Check batches have the correct starting indices
    for (int i = 0; i < num_grid_pt; ++i) {
        if (i == 13) {
            EXPECT_EQ(mpm_transfer.get_num_particles_in_batch(i), 3);
            EXPECT_EQ(mpm_transfer.get_batch_starting_index(i), 0);
        } else if (i < 13) {
            EXPECT_EQ(mpm_transfer.get_num_particles_in_batch(i), 0);
            EXPECT_EQ(mpm_transfer.get_batch_starting_index(i), 0);
        } else {
            EXPECT_EQ(mpm_transfer.get_num_particles_in_batch(i), 0);
            EXPECT_EQ(mpm_transfer.get_batch_starting_index(i), 3);
        }
    }

    // A more comprehensive test, where we combine particles in above two
    // cases together:
    // The batch_starting_index_ shall look like:
    //            i=13
    // (0, 1, ..., 13, 17, 18, ...)
    // and the number of particles in the batch shall look like:
    // (1, 1, ..., 4, 1, 1, ...)
    particles = Particles(30);
    mpm_transfer = MPMTransfer(grid, particles);
    num_grid_pt = grid.get_num_gridpt();

    particles.set_position(0, Vector3<double>(-0.5, 0.5, -0.5));
    particles.set_position(1, Vector3<double>(0.5, 0.5, 0.5));
    particles.set_position(2, Vector3<double>(0.5, -0.5, 0.5));
    pc = num_grid_pt;
    for (int k = bottom_corner(2); k < bottom_corner(2)+num_grid_pt_1D(2);
                                                                        ++k) {
    for (int j = bottom_corner(1); j < bottom_corner(1)+num_grid_pt_1D(1);
                                                                        ++j) {
    for (int i = bottom_corner(0); i < bottom_corner(0)+num_grid_pt_1D(0);
                                                                        ++i) {
        particles.set_position(--pc, grid.get_position(i, j, k));
    }
    }
    }

    mpm_transfer.SortParticles(grid, &particles);

    // Check sorting
    EXPECT_TRUE(CompareMatrices(particles.get_position(0),
                                Vector3<double>(-2.0, -2.0, -2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(1),
                                Vector3<double>(0.0, -2.0, -2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(8),
                                Vector3<double>(2.0, 2.0, -2.0),
                                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_position(29),
                                Vector3<double>(2.0, 2.0, 2.0),
                                std::numeric_limits<double>::epsilon()));

    // Check batches have the correct starting indices
    for (int i = 0; i < num_grid_pt; ++i) {
        if (i == 13) {
            EXPECT_EQ(mpm_transfer.get_num_particles_in_batch(i), 4);
            EXPECT_EQ(mpm_transfer.get_batch_starting_index(i), 13);
        } else if (i < 13) {
            EXPECT_EQ(mpm_transfer.get_num_particles_in_batch(i), 1);
            EXPECT_EQ(mpm_transfer.get_batch_starting_index(i), i);
        } else {
            EXPECT_EQ(mpm_transfer.get_num_particles_in_batch(i), 1);
            EXPECT_EQ(mpm_transfer.get_batch_starting_index(i), i+3);
        }
    }
}

}  // namespace
}  // namespace internal
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
