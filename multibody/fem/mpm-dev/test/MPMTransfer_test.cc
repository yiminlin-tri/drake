#include "drake/multibody/fem/mpm-dev/MPMTransfer.h"

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace mpm {

constexpr double kEps = 4.0 * std::numeric_limits<double>::epsilon();
constexpr double TOLERANCE = 1e-10;

class MPMTransferTest : public ::testing::Test {
 protected:
    void SetUp() { }

    void CheckSortOutOfBounds() {
        // First, construct a 3x3x3 grid centered at (0.0, 0.0, 0.0), with
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
        // And make 27 particles on the location of 27 grid points. The
        // particles' order are in a reversed lexiographical ordering.
        Vector3<int> num_grid_pt_1D, bottom_corner;
        double h;
        int num_grid_pt, num_particles, pc;

        num_grid_pt_1D = Vector3<int>(3, 3, 3);
        h = 2.0;
        bottom_corner = Vector3<int>(-1, -1, -1);

        Grid grid = Grid(num_grid_pt_1D, h, bottom_corner);
        num_particles = 27;
        std::unique_ptr<Particles> particles =
                                    std::make_unique<Particles>(num_particles);
        EXPECT_EQ(particles->get_num_particles(), num_particles);
        MPMTransfer mpm_transfer = MPMTransfer();
        num_grid_pt = grid.get_num_gridpt();

        EXPECT_EQ(num_grid_pt, 27);
        EXPECT_EQ(num_grid_pt_1D(0), 3);
        EXPECT_EQ(num_grid_pt_1D(1), 3);
        EXPECT_EQ(num_grid_pt_1D(2), 3);
        pc = particles->get_num_particles();
        for (int k = bottom_corner(2);
                                k < bottom_corner(2)+num_grid_pt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                                j < bottom_corner(1)+num_grid_pt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                                i < bottom_corner(0)+num_grid_pt_1D(0); ++i) {
            particles->set_position(--pc, grid.get_position(i, j, k));
        }
        }
        }

        // Sanity check
        EXPECT_TRUE(CompareMatrices(particles->get_position(0),
                                    Vector3<double>(2.0, 2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(1),
                                    Vector3<double>(0.0, 2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(8),
                                    Vector3<double>(-2.0, -2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(13),
                                    Vector3<double>(0.0, 0.0, 0.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(26),
                                    Vector3<double>(-2.0, -2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));

        // Check particles out of bound error
        EXPECT_THROW(mpm_transfer.SortParticles(grid, particles.get()),
                     std::exception);
    }

    void CheckSort1() {
        // Construct a grid of 5x5x5 on [-4,4]^3, and place 27 particles
        // on the centering 3x3x3 grid points.
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
        Vector3<int> num_grid_pt_1D, bottom_corner;
        double h;
        int num_grid_pt, num_particles, pc;

        num_grid_pt_1D = Vector3<int>(5, 5, 5);
        h = 2.0;
        bottom_corner = Vector3<int>(-2, -2, -2);

        Grid grid = Grid(num_grid_pt_1D, h, bottom_corner);
        num_particles = 27;
        std::unique_ptr<Particles> particles =
                                    std::make_unique<Particles>(num_particles);
        EXPECT_EQ(particles->get_num_particles(), num_particles);
        MPMTransfer mpm_transfer = MPMTransfer();
        num_grid_pt = grid.get_num_gridpt();

        EXPECT_EQ(num_grid_pt, 125);
        EXPECT_EQ(num_grid_pt_1D(0), 5);
        EXPECT_EQ(num_grid_pt_1D(1), 5);
        EXPECT_EQ(num_grid_pt_1D(2), 5);
        pc = particles->get_num_particles();
        for (int k = bottom_corner(2)+1;
                                k < bottom_corner(2)+num_grid_pt_1D(2)-1; ++k) {
        for (int j = bottom_corner(1)+1;
                                j < bottom_corner(1)+num_grid_pt_1D(1)-1; ++j) {
        for (int i = bottom_corner(0)+1;
                                i < bottom_corner(0)+num_grid_pt_1D(0)-1; ++i) {
            particles->set_position(--pc, grid.get_position(i, j, k));
        }
        }
        }

        // Sanity check
        EXPECT_TRUE(CompareMatrices(particles->get_position(0),
                                    Vector3<double>(2.0, 2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(1),
                                    Vector3<double>(0.0, 2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(8),
                                    Vector3<double>(-2.0, -2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(13),
                                    Vector3<double>(0.0, 0.0, 0.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(26),
                                    Vector3<double>(-2.0, -2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));

        // Check particles in the correct ordering after sorting
        mpm_transfer.SortParticles(grid, particles.get());

        EXPECT_TRUE(CompareMatrices(particles->get_position(0),
                                    Vector3<double>(-2.0, -2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(1),
                                    Vector3<double>(0.0, -2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(8),
                                    Vector3<double>(2.0, 2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(13),
                                    Vector3<double>(0.0, 0.0, 0.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(26),
                                    Vector3<double>(2.0, 2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));

        pc = 0;
        for (int k = bottom_corner(2);
                                k < bottom_corner(2)+num_grid_pt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                                j < bottom_corner(1)+num_grid_pt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                                i < bottom_corner(0)+num_grid_pt_1D(0); ++i) {
            if (std::abs(i) <= 1 && std::abs(j) <= 1 && std::abs(k) <= 1) {
                EXPECT_EQ(mpm_transfer.batch_sizes_[pc++], 1);
            } else {
                EXPECT_EQ(mpm_transfer.batch_sizes_[pc++], 0);
            }
        }
        }
        }
    }

    void CheckSort2() {
        // A more comprehensive test, where we add three particles in the batch
        // centered at (0, 0, 0) to the particles constructed in CheckSort1:
        // The number of particles in the batch shall look like:
        // (1, 1, ..., 4, 1, 1, ...)
        Vector3<int> num_grid_pt_1D, bottom_corner;
        double h;
        int num_particles, pc;

        num_grid_pt_1D = Vector3<int>(5, 5, 5);
        h = 2.0;
        bottom_corner = Vector3<int>(-2, -2, -2);

        num_particles = 30;
        Grid grid = Grid(num_grid_pt_1D, h, bottom_corner);
        std::unique_ptr<Particles> particles =
                                    std::make_unique<Particles>(num_particles);
        MPMTransfer mpm_transfer = MPMTransfer();

        particles->set_position(0, Vector3<double>(-0.5, 0.5, -0.5));
        particles->set_position(1, Vector3<double>(0.5, 0.5, 0.5));
        particles->set_position(2, Vector3<double>(0.5, -0.5, 0.5));
        pc = num_particles;
        for (int k = bottom_corner(2)+1;
                                k < bottom_corner(2)+num_grid_pt_1D(2)-1; ++k) {
        for (int j = bottom_corner(1)+1;
                                j < bottom_corner(1)+num_grid_pt_1D(1)-1; ++j) {
        for (int i = bottom_corner(0)+1;
                                i < bottom_corner(0)+num_grid_pt_1D(0)-1; ++i) {
            particles->set_position(--pc, grid.get_position(i, j, k));
        }
        }
        }

        // Sanity check
        EXPECT_TRUE(CompareMatrices(particles->get_position(0),
                                    Vector3<double>(-0.5, 0.5, -0.5),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(1),
                                    Vector3<double>(0.5, 0.5, 0.5),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(2),
                                    Vector3<double>(0.5, -0.5, 0.5),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(3),
                                    Vector3<double>(2.0, 2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(29),
                                    Vector3<double>(-2.0, -2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));

        mpm_transfer.SortParticles(grid, particles.get());

        // Check sorting
        EXPECT_TRUE(CompareMatrices(particles->get_position(0),
                                    Vector3<double>(-2.0, -2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(1),
                                    Vector3<double>(0.0, -2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(8),
                                    Vector3<double>(2.0, 2.0, -2.0),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(29),
                                    Vector3<double>(2.0, 2.0, 2.0),
                                    std::numeric_limits<double>::epsilon()));

        pc = 0;
        for (int k = bottom_corner(2);
                                k < bottom_corner(2)+num_grid_pt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                                j < bottom_corner(1)+num_grid_pt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                                i < bottom_corner(0)+num_grid_pt_1D(0); ++i) {
            if (std::abs(i) <= 1 && std::abs(j) <= 1 && std::abs(k) <= 1) {
                if (i == 0 && j == 0 && k == 0) {
                    EXPECT_EQ(mpm_transfer.batch_sizes_[pc++], 4);
                } else {
                    EXPECT_EQ(mpm_transfer.batch_sizes_[pc++], 1);
                }
            } else {
                EXPECT_EQ(mpm_transfer.batch_sizes_[pc++], 0);
            }
        }
        }
        }
    }

    void CheckSort3() {
        // A test with the same grid case 1 and 2, and we add a particle to
        // batches centered at (0, 0, 0), (1, 1, 1), and (-1, -1, -1).
        // The number of particles in the batch shall look like:
        // i=0           i=13          i=26
        // (1, 0, ..., 0, 1, 0, ..., 0, 1)
        Vector3<int> num_grid_pt_1D, bottom_corner;
        double h;
        int pc, num_particles;

        num_grid_pt_1D = Vector3<int>(5, 5, 5);
        h = 2.0;
        bottom_corner = Vector3<int>(-2, -2, -2);

        num_particles = 3;
        Grid grid = Grid(num_grid_pt_1D, h, bottom_corner);
        std::unique_ptr<Particles> particles =
                                    std::make_unique<Particles>(num_particles);
        MPMTransfer mpm_transfer = MPMTransfer();

        particles->set_position(0, Vector3<double>(1.5, 1.5, 1.5));
        particles->set_position(1, Vector3<double>(-1.5, -1.5, -1.5));
        particles->set_position(2, Vector3<double>(-0.5, 0.5, -0.5));

        // Sanity check
        mpm_transfer.SortParticles(grid, particles.get());

        // Check sorting
        EXPECT_TRUE(CompareMatrices(particles->get_position(0),
                                    Vector3<double>(-1.5, -1.5, -1.5),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(1),
                                    Vector3<double>(-0.5, 0.5, -0.5),
                                    std::numeric_limits<double>::epsilon()));
        EXPECT_TRUE(CompareMatrices(particles->get_position(2),
                                    Vector3<double>(1.5, 1.5, 1.5),
                                    std::numeric_limits<double>::epsilon()));
        pc = 0;
        for (int k = bottom_corner(2);
                                k < bottom_corner(2)+num_grid_pt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                                j < bottom_corner(1)+num_grid_pt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                                i < bottom_corner(0)+num_grid_pt_1D(0); ++i) {
            if ((i == 0 && j == 0 && k == 0) ||
                (i == 1 && j == 1 && k == 1) ||
                (i == -1 && j == -1 && k == -1)) {
                EXPECT_EQ(mpm_transfer.batch_sizes_[pc++], 1);
            } else {
                EXPECT_EQ(mpm_transfer.batch_sizes_[pc++], 0);
            }
        }
        }
        }
    }

    void checkPreallocation() {
        // Construct a grid of 5x5x5 on [-2,2]^3, and place 27 particles
        // on the centering 3x3x3 grid points.
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
        double sum_val;
        Vector3<double> xp, sum_gradient;
        std::vector<std::array<double, 27>> bases_val_particles;
        std::vector<std::array<Vector3<double>, 27>> bases_grad_particles;
        int h = 1.0;
        Vector3<int> num_gridpt_1D = { 5,  5,  5};
        Vector3<int> bottom_corner = {-2, -2, -2};
        Grid grid = Grid(num_gridpt_1D, h, bottom_corner);
        int num_particles = 27;
        std::unique_ptr<Particles> particles =
                                    std::make_unique<Particles>(num_particles);
        MPMTransfer mpm_transfer = MPMTransfer();

        // Set particles' positions to be on grid points
        pc = num_particles;
        for (int k = bottom_corner(2)+1;
                 k < bottom_corner(2)+num_gridpt_1D(2)-1; ++k) {
        for (int j = bottom_corner(1)+1;
                 j < bottom_corner(1)+num_gridpt_1D(1)-1; ++j) {
        for (int i = bottom_corner(0)+1;
                 i < bottom_corner(0)+num_gridpt_1D(0)-1; ++i) {
            particles->set_position(--pc, grid.get_position(i, j, k));
        }
        }
        }

        // Sort the particles and set up the batches and preallocate basis
        // evaluations
        mpm_transfer.SetUpTransfer(grid, particles.get());

        // The particles are sorted, and for all particles, all bases that cover
        // the particle shall have evaluations sum to 1, and gradients sum to
        // 0, by the partition of unity property.
        for (int p = 0; p < num_particles; ++p) {
            EXPECT_EQ(mpm_transfer.bases_val_particles_[p][13], 0.75*0.75*0.75);
            xp = particles->get_position(p);
            sum_val = 0.0;
            sum_gradient = {0.0, 0.0, 0.0};
            for (int i = 0; i < 27; ++i) {
                sum_val += mpm_transfer.bases_val_particles_[p][i];
                sum_gradient += mpm_transfer.bases_grad_particles_[p][i];
            }
            EXPECT_EQ(sum_val, 1.0);
            EXPECT_TRUE(CompareMatrices(sum_gradient,
                                        Vector3<double>::Zero(), kEps));
        }
    }

    void checkP2GForce() {
        // Construct a grid of 3x3x3 on [-2,2]^3, and place 1 particles
        // at the center
        //        -2   0   2
        //         o - o - o
        //         |   |   |
        //     o - o - o - o
        //     |   |   |   |
        // o - o - o - o - o
        // |   |   |   |
        // o - o - o - o
        // |   |   |
        // o - o - o
        double mass_p, reference_volume_p, sum_mass_grid;
        Vector3<double> momentum_p, sum_momentum_grid;
        Matrix3<double> tau_p;
        int h = 2.0;
        Vector3<int> num_gridpt_1D = { 3,  3,  3};
        Vector3<int> bottom_corner = {-1, -1, -1};
        int num_particles = 1;
        grid_ = std::make_unique<Grid>(num_gridpt_1D, h, bottom_corner);
        particles_ = std::make_unique<Particles>(num_particles);
        mpm_transfer_ = std::make_unique<MPMTransfer>();

        // Set particles' positions to be on grid points
        mass_p = 2.0;
        reference_volume_p = 0.1;
        tau_p = 3.0*Matrix3<double>::Identity();
        momentum_p = {0.2, -0.4, 0.6};
        particles_->set_position(0, grid_->get_position(0, 0, 0));
        particles_->set_mass(0, mass_p);
        particles_->set_reference_volume(0, reference_volume_p);
        particles_->set_velocity(0, momentum_p/mass_p);
        particles_->set_kirchhoff_stress(0, tau_p);

        // Sort the particles and set up the batches and preallocate basis
        // evaluations
        mpm_transfer_->SetUpTransfer(*grid_, particles_.get());

        // Transfer particles' information to grid
        mpm_transfer_->TransferParticlesToGrid(*particles_, grid_.get());
        sum_mass_grid = 0.0;
        sum_momentum_grid = {0.0, 0.0, 0.0};
        for (int k = bottom_corner(2);
                 k < bottom_corner(2)+num_gridpt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                 j < bottom_corner(1)+num_gridpt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                 i < bottom_corner(0)+num_gridpt_1D(0); ++i) {
            sum_mass_grid += grid_->get_mass(i, j, k);
            sum_momentum_grid += grid_->get_mass(i, j, k)
                                *grid_->get_velocity(i, j, k);
        }
        }
        }

        // Verify the conservation of mass and momentum
        EXPECT_NEAR(mass_p, sum_mass_grid, TOLERANCE);
        EXPECT_TRUE(CompareMatrices(momentum_p, sum_momentum_grid, TOLERANCE));

        // Test Force
        // At grid point (0, 0, 0), since the particle locates at this grid
        // point, \grad N_(i, j, k)(x_p) = 0 in this case.
        EXPECT_TRUE(CompareMatrices(grid_->get_force(0, 0, 0),
                                    Vector3<double>::Zero(), TOLERANCE));

        // At other grid points, every components of the gradient look like
        // /grad N_(i, j, k)(x_p) = 1/h * \grad N * N * N, the product of
        // a gradient and two evaluations. Evaluations are the same for all
        // grid points, by symmetry of the grid. N = 1/8 in this case.
        // The gradients' compoents have value \pm 1/2, depending on the actual
        // location of the grid point. So /grad N_(i, j, k)(x_p) =
        // [\pm 1/256, \pm 1/256, \pm 1/256].
        // Since we have only one grid point, the forces have values
        // -V0_p*tau_p*[\pm 1/256, \pm 1/256, \pm 1/256]
        Matrix3<double> scale = -1.0/256.0*reference_volume_p*tau_p;
        for (int k = -1; k <= 1; ++k) {
        for (int j = -1; j <= 1; ++j) {
        for (int i = -1; i <= 1; ++i) {
        if ((i != 0) && (j != 0) && (k != 0)) {
            EXPECT_TRUE(CompareMatrices(grid_->get_force(i, j, k),
                                    scale*Vector3<double>(i, j, k), TOLERANCE));
        }
        }
        }
        }
    }

    void checkP2GMassVelocity1() {
        // Construct a grid of 5x5x5 on [-2,2]^3, and place 27 particles
        // on the centering 3x3x3 grid points.
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
        double sum_mass_particles, sum_mass_grid;
        Vector3<double> sum_momentum_particles, sum_momentum_grid;
        int h = 1.0;
        Vector3<int> num_gridpt_1D = { 5,  5,  5};
        Vector3<int> bottom_corner = {-2, -2, -2};
        int num_particles = 27;
        grid_ = std::make_unique<Grid>(num_gridpt_1D, h, bottom_corner);
        particles_ = std::make_unique<Particles>(num_particles);
        mpm_transfer_ = std::make_unique<MPMTransfer>();

        // Set particles' positions to be on grid points
        double dummy_mass, dummy_reference_volume;
        Vector3<double> dummy_velocity, dummy_momentum;
        Matrix3<double> dummy_tau;
        pc = num_particles;
        sum_mass_particles = 0.0;
        sum_momentum_particles = {0.0, 0.0, 0.0};
        for (int k = bottom_corner(2)+1;
                 k < bottom_corner(2)+num_gridpt_1D(2)-1; ++k) {
        for (int j = bottom_corner(1)+1;
                 j < bottom_corner(1)+num_gridpt_1D(1)-1; ++j) {
        for (int i = bottom_corner(0)+1;
                 i < bottom_corner(0)+num_gridpt_1D(0)-1; ++i) {
            dummy_mass = pc+0.1;
            dummy_reference_volume = 2*pc+0.1;
            dummy_velocity = Vector3<double>(1.1*pc, 1.2*pc, 1.3*pc);
            dummy_momentum = dummy_mass*dummy_velocity;
            dummy_tau = pc*Matrix3<double>::Identity();
            --pc;
            particles_->set_position(pc, grid_->get_position(i, j, k));
            particles_->set_mass(pc, dummy_mass);
            sum_mass_particles += dummy_mass;
            particles_->set_reference_volume(pc, dummy_reference_volume);
            particles_->set_velocity(pc, dummy_velocity);
            sum_momentum_particles += dummy_momentum;
            particles_->set_kirchhoff_stress(pc, dummy_tau);
        }
        }
        }

        // Sort the particles and set up the batches and preallocate basis
        // evaluations
        mpm_transfer_->SetUpTransfer(*grid_, particles_.get());

        // Transfer particles' information to grid
        mpm_transfer_->TransferParticlesToGrid(*particles_, grid_.get());
        sum_mass_grid = 0.0;
        sum_momentum_grid = {0.0, 0.0, 0.0};
        for (int k = bottom_corner(2);
                 k < bottom_corner(2)+num_gridpt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                 j < bottom_corner(1)+num_gridpt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                 i < bottom_corner(0)+num_gridpt_1D(0); ++i) {
            sum_mass_grid += grid_->get_mass(i, j, k);
            sum_momentum_grid += grid_->get_mass(i, j, k)
                                *grid_->get_velocity(i, j, k);
        }
        }
        }

        // Verify the conservation of mass and momentum
        EXPECT_NEAR(sum_mass_particles, sum_mass_grid, TOLERANCE);
        EXPECT_TRUE(CompareMatrices(sum_momentum_particles, sum_momentum_grid,
                                    TOLERANCE));
    }

    void checkP2GMassVelocity2() {
        // Construct a grid of 5x5x5 on [-2,2]^3, and place 27 particles
        // on the centering 3x3x3 grid points, then add 3 particles inside
        // the grid
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
        double sum_mass_particles, sum_mass_grid;
        Vector3<double> sum_momentum_particles, sum_momentum_grid;
        int h = 1.0;
        Vector3<int> num_gridpt_1D = { 5,  5,  5};
        Vector3<int> bottom_corner = {-2, -2, -2};
        int num_particles = 27;
        grid_ = std::make_unique<Grid>(num_gridpt_1D, h, bottom_corner);
        particles_ = std::make_unique<Particles>(num_particles);
        mpm_transfer_ = std::make_unique<MPMTransfer>();

        // Set particles' positions to be on grid points
        double dummy_mass, dummy_reference_volume;
        Vector3<double> dummy_velocity, dummy_momentum;
        Matrix3<double> dummy_tau;
        pc = num_particles;
        sum_mass_particles = 0.0;
        sum_momentum_particles = {0.0, 0.0, 0.0};
        for (int k = bottom_corner(2)+1;
                 k < bottom_corner(2)+num_gridpt_1D(2)-1; ++k) {
        for (int j = bottom_corner(1)+1;
                 j < bottom_corner(1)+num_gridpt_1D(1)-1; ++j) {
        for (int i = bottom_corner(0)+1;
                 i < bottom_corner(0)+num_gridpt_1D(0)-1; ++i) {
            dummy_mass = pc+0.1;
            dummy_reference_volume = 2*pc+0.1;
            dummy_velocity = Vector3<double>(1.1*pc, 1.2*pc, 1.3*pc);
            dummy_momentum = dummy_mass*dummy_velocity;
            dummy_tau = pc*Matrix3<double>::Identity();
            --pc;
            particles_->set_position(pc, grid_->get_position(i, j, k));
            particles_->set_mass(pc, dummy_mass);
            sum_mass_particles += dummy_mass;
            particles_->set_reference_volume(pc, dummy_reference_volume);
            particles_->set_velocity(pc, dummy_velocity);
            sum_momentum_particles += dummy_momentum;
            particles_->set_kirchhoff_stress(pc, dummy_tau);
        }
        }
        }

        // Add more particles
        Vector3<double> pos1 = {0.2, 0.1, 0.3};
        Vector3<double> vel1 = {-1.0, -2.0, -3.0};
        double mass1 = 5.0;
        double vol1  = 10.0;
        Matrix3<double> F1 = pos1.asDiagonal();
        Matrix3<double> stress1 = vel1.asDiagonal();

        Vector3<double> pos2 = {0.3, -0.1, 0.6};
        Vector3<double> vel2 = {-9.0, 8.0, -2.0};
        double mass2 = 7.0;
        double vol2  = 3.0;
        Matrix3<double> F2 = pos2.asDiagonal();
        Matrix3<double> stress2 = vel2.asDiagonal();

        Vector3<double> pos3 = {0.2, -0.5, 0.3};
        Vector3<double> vel3 = {2.0, -6.2, 8.0};
        double mass3 = 2.0;
        double vol3  = 12.0;
        Matrix3<double> F3 = pos3.asDiagonal();
        Matrix3<double> stress3 = vel3.asDiagonal();

        particles_->AddParticle(pos1, vel1, mass1, vol1, F1, stress1);
        particles_->AddParticle(pos2, vel2, mass2, vol2, F2, stress2);
        particles_->AddParticle(pos3, vel3, mass3, vol3, F3, stress3);

        num_particles = particles_->get_num_particles();
        sum_mass_particles += (mass1 + mass2 + mass3);
        sum_momentum_particles += (mass1*vel1 + mass2*vel2 + mass3*vel3);

        // Sort the particles and set up the batches and preallocate basis
        // evaluations
        mpm_transfer_->SetUpTransfer(*grid_, particles_.get());

        // Transfer particles' information to grid
        mpm_transfer_->TransferParticlesToGrid(*particles_, grid_.get());
        sum_mass_grid = 0.0;
        sum_momentum_grid = {0.0, 0.0, 0.0};
        for (int k = bottom_corner(2);
                 k < bottom_corner(2)+num_gridpt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                 j < bottom_corner(1)+num_gridpt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                 i < bottom_corner(0)+num_gridpt_1D(0); ++i) {
            sum_mass_grid += grid_->get_mass(i, j, k);
            sum_momentum_grid += grid_->get_mass(i, j, k)
                                *grid_->get_velocity(i, j, k);
        }
        }
        }

        // Verify the conservation of mass and momentum
        EXPECT_NEAR(sum_mass_particles, sum_mass_grid, TOLERANCE);
        EXPECT_TRUE(CompareMatrices(sum_momentum_particles, sum_momentum_grid,
                                    TOLERANCE));
    }

    void checkG2PDeformationGrad() {
        // Construct a grid of 3x3x3 on [-2,2]^3, and place 1 particles
        // at the center
        //        -2   0   2
        //         o - o - o
        //         |   |   |
        //     o - o - o - o
        //     |   |   |   |
        // o - o - o - o - o
        // |   |   |   |
        // o - o - o - o
        // |   |   |
        // o - o - o
        double mass_p;
        Matrix3<double> F_p;
        int h = 2.0;
        Vector3<int> num_gridpt_1D = { 3,  3,  3};
        Vector3<int> bottom_corner = {-1, -1, -1};
        Vector3<double> velocity_i = {0.1, 0.2, 0.3};
        int num_particles = 1;
        double dt = 0.1;
        grid_ = std::make_unique<Grid>(num_gridpt_1D, h, bottom_corner);
        particles_ = std::make_unique<Particles>(num_particles);
        mpm_transfer_ = std::make_unique<MPMTransfer>();

        // Initialize particles' positions to be on grid points
        mass_p = 0.2;
        F_p   = 2.0*Matrix3<double>::Identity();
        particles_->set_position(0, grid_->get_position(0, 0, 0));
        particles_->set_mass(0, mass_p);
        particles_->set_deformation_gradient(0, F_p);

        // Initialize the grid velocity field. We assume it is constant,
        // so the deformation gradient would not change
        for (int k = bottom_corner(2);
                 k < bottom_corner(2)+num_gridpt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                 j < bottom_corner(1)+num_gridpt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                 i < bottom_corner(0)+num_gridpt_1D(0); ++i) {
            grid_->set_velocity(i, j, k, velocity_i);
        }
        }
        }

        // Sort the particles and set up the batches and preallocate basis
        // evaluations
        mpm_transfer_->SetUpTransfer(*grid_, particles_.get());

        // Grid to particles transfer
        mpm_transfer_->TransferGridToParticles(*grid_, dt, particles_.get());

        // Since the velocity field is constant, the deformation gradient
        // doesn't change with respect to time
        EXPECT_TRUE(CompareMatrices(F_p,
                                    particles_->get_deformation_gradient(0),
                                    TOLERANCE));
        EXPECT_TRUE(CompareMatrices(velocity_i, particles_->get_velocity(0),
                                    TOLERANCE));
    }

    void checkG2PMassVelocity1() {
        double sum_mass_particle, sum_mass_grid;
        Vector3<double> sum_momentum_particle, sum_momentum_grid;
        double dt = 0.1;

        // First setup particles and grid using P2G
        checkP2GMassVelocity1();
        Vector3<int> bottom_corner = grid_->get_bottom_corner();
        Vector3<int> num_gridpt_1D = grid_->get_num_gridpt_1D();
        int num_particles = particles_->get_num_particles();

        // Then we do a grid to particle transfer
        mpm_transfer_->TransferGridToParticles(*grid_, dt, particles_.get());

        // Particles' sum of mass and momentum
        sum_mass_particle = 0.0;
        sum_momentum_particle = {0.0, 0.0, 0.0};
        for (int p = 0; p < num_particles; ++p) {
            sum_mass_particle += particles_->get_mass(p);
            sum_momentum_particle += particles_->get_mass(p)
                            *particles_->get_velocity(p);
        }

        // Grid's sum of mass and momentum
        sum_mass_grid = 0.0;
        sum_momentum_grid = {0.0, 0.0, 0.0};
        for (int k = bottom_corner(2);
                 k < bottom_corner(2)+num_gridpt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                 j < bottom_corner(1)+num_gridpt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                 i < bottom_corner(0)+num_gridpt_1D(0); ++i) {
            sum_mass_grid += grid_->get_mass(i, j, k);
            sum_momentum_grid += grid_->get_mass(i, j, k)
                            *grid_->get_velocity(i, j, k);
        }
        }
        }

        // Verify the conservation of mass and momentum
        EXPECT_NEAR(sum_mass_particle, sum_mass_grid, TOLERANCE);
        EXPECT_TRUE(CompareMatrices(sum_momentum_particle, sum_momentum_grid,
                                    TOLERANCE));
    }

    void checkG2PMassVelocity2() {
        double sum_mass_particle, sum_mass_grid;
        Vector3<double> sum_momentum_particle, sum_momentum_grid;
        double dt = 0.1;

        // First setup particles and grid using P2G
        checkP2GMassVelocity2();
        Vector3<int> bottom_corner = grid_->get_bottom_corner();
        Vector3<int> num_gridpt_1D = grid_->get_num_gridpt_1D();
        int num_particles = particles_->get_num_particles();

        // Then we do a grid to particle transfer
        mpm_transfer_->TransferGridToParticles(*grid_, dt, particles_.get());

        // Particles' sum of mass and momentum
        sum_mass_particle = 0.0;
        sum_momentum_particle = {0.0, 0.0, 0.0};
        for (int p = 0; p < num_particles; ++p) {
            sum_mass_particle += particles_->get_mass(p);
            sum_momentum_particle += particles_->get_mass(p)
                            *particles_->get_velocity(p);
        }

        // Grid's sum of mass and momentum
        sum_mass_grid = 0.0;
        sum_momentum_grid = {0.0, 0.0, 0.0};
        for (int k = bottom_corner(2);
                 k < bottom_corner(2)+num_gridpt_1D(2); ++k) {
        for (int j = bottom_corner(1);
                 j < bottom_corner(1)+num_gridpt_1D(1); ++j) {
        for (int i = bottom_corner(0);
                 i < bottom_corner(0)+num_gridpt_1D(0); ++i) {
            sum_mass_grid += grid_->get_mass(i, j, k);
            sum_momentum_grid += grid_->get_mass(i, j, k)
                            *grid_->get_velocity(i, j, k);
        }
        }
        }

        // Verify the conservation of mass and momentum
        EXPECT_NEAR(sum_mass_particle, sum_mass_grid, TOLERANCE);
        EXPECT_TRUE(CompareMatrices(sum_momentum_particle, sum_momentum_grid,
                                    TOLERANCE));
    }

    std::unique_ptr<Particles> particles_;
    std::unique_ptr<Grid> grid_;
    std::unique_ptr<MPMTransfer> mpm_transfer_;
};

namespace {

TEST_F(MPMTransferTest, SortParticlesTest) {
    CheckSortOutOfBounds();
    CheckSort1();
    CheckSort2();
    CheckSort3();
}

TEST_F(MPMTransferTest, SetUpTest) {
    checkPreallocation();
}

TEST_F(MPMTransferTest, P2GTest) {
    checkP2GForce();
    checkP2GMassVelocity1();
    checkP2GMassVelocity2();
}

TEST_F(MPMTransferTest, G2PTest) {
    checkG2PDeformationGrad();
    checkG2PMassVelocity1();
    checkG2PMassVelocity2();
}

}  // namespace
}  // namespace mpm
}  // namespace multibody
}  // namespace drake