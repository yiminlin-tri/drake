#include "drake/multibody/fem/mpm-dev/ExternalForces.h"

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"
#include "drake/multibody/fem/mpm-dev/Grid.h"

namespace drake {
namespace multibody {
namespace mpm {

constexpr double TOLERANCE = 1e-10;

namespace {

GTEST_TEST(ExternalForcesTest, TestApplyGravitationalForces) {
    Vector3<int> num_gridpt_1D = {2, 2, 2};
    double h = 1.0;
    Vector3<int> bottom_corner  = {0, 0, 0};
    Grid grid = Grid(num_gridpt_1D, h, bottom_corner);
    double dt = 0.1;

    for (const auto& [batch_index_flat, batch_index_3d] : grid.get_indices()) {
        Vector3<double> velocity_i =
                            batch_index_flat*Vector3<double>(1.0, 1.0, 1.0);
        double mass_i = batch_index_flat;
        grid.set_velocity(batch_index_3d(0), batch_index_3d(1),
                          batch_index_3d(2), velocity_i);
        grid.set_mass(batch_index_3d(0), batch_index_3d(1),
                      batch_index_3d(2), mass_i);
    }

    // Apply gravitational force
    ExternalForces external_forces = ExternalForces();
    external_forces.ApplyGravitationalForces(dt, &grid);

    for (const auto& [batch_index_flat, batch_index_3d] : grid.get_indices()) {
        Vector3<double> velocity_exact =
                            batch_index_flat*Vector3<double>(1.0, 1.0, 1.0)
                          + dt*Vector3<double>(0.0, 0.0, -9.81);
        double mass_exact = batch_index_flat;
        EXPECT_EQ(mass_exact, grid.get_mass(batch_index_3d(0),
                                            batch_index_3d(1),
                                            batch_index_3d(2)));
        EXPECT_TRUE(CompareMatrices(velocity_exact,
                                    grid.get_velocity(batch_index_3d(0),
                                                      batch_index_3d(1),
                                                      batch_index_3d(2)),
                                    TOLERANCE));
    }
}

}  // namespace
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
