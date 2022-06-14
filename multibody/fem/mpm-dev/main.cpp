#include <array>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <Partio.h>

#include "drake/common/drake_assert.h"
#include "drake/common/eigen_types.h"
#include "drake/common/filesystem.h"
#include "drake/common/temp_directory.h"
#include "drake/multibody/fem/mpm-dev/MPMDriver.h"

namespace drake {
namespace multibody {
namespace mpm {

int DoMain() {
    MPMParameters::PhysicalParameters p_param {
        8e4,                                   // Young's modulus
        0.49,                                  // Poisson's ratio
        0.10,                                  // Friction coefficient
        {0.0, 0.0, -9.81}                      // Gravitational acceleration
    };

    MPMParameters::InitialConditionParameters init_param {
        1272,                                  // Density
        {0.0, 0.0, 0.0},                       // Initial Velocity
    };

    MPMParameters::SolverParameters s_param {
        5e-1,                                  // End time
        2e-4,                                  // Time step size
        1,                                     // Min num of particles per cell
        0.01,                                  // Grid size
        Vector3<int>(23, 23, 23),              // Number of grid points in each
                                               // direction
        Vector3<int>(-1, -1, -1),              // Bottom corner of the grid
    };

    MPMParameters::IOParameters io_param {
        "mpm-test",                            // case name
        "/home/yiminlin/Desktop/output",       // output directory name
        1e-3,                                  // Interval of outputting
    };

    MPMParameters param {p_param, init_param, s_param, io_param};
    auto driver = std::make_unique<MPMDriver>(std::move(param));

    BoundaryCondition::Boundary b0 = {p_param.mu, {{1, 0, 0}, {0, 0, 0}}};
    BoundaryCondition::Boundary b1 = {p_param.mu, {{-1, 0, std::sqrt(3)},
                                                   {0, 0, 0}}};
    std::vector<BoundaryCondition::Boundary> boundaries = {b0, b1};

    Vector3<double> xscale = {0.02, 0.02, 0.02};
    Vector3<double> translation = {0.18, 0.1, 0.18};
    BoxLevelSet level_set = BoxLevelSet(xscale);
    // Translate the level set by:
    math::RigidTransform<double> pose
                                = math::RigidTransform<double>(translation);

    driver->InitializeBoundaryConditions(std::move(boundaries));
    driver->Run(level_set, pose);

    return 0;
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake

int main() {
    return drake::multibody::mpm::DoMain();
}
