#include <iostream>
#include <memory>
#include <random>
#include <string>
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

bool box_indicator(const Vector3<double>& position) {
    if ((position(0) <= 6) && (position(0) >= 4)
     && (position(1) <= 6) && (position(1) >= 4)
     && (position(2) <= 5) && (position(2) >= 3)) {
        return true;
    } else {
        return false;
    }
}

double constant_density(const Vector3<double>& position) {
    return position(0)-position(0)+1.0;
}

Vector3<double> constant_velocity(const Vector3<double>& position) {
    return position-position+Vector3<double>(0.0, 0.0, 1.0);
}

int DoMain() {
    // TODO(yiminlin.tri): ignore level set and use [4, 6]x[4, 6]x[3, 5] box now
    // And we use a grid of size [0,10]^3
    MPMDriver::MPMParameters param {
        "mpm-test",                        // case name
        "/home/yiminlin/Desktop/output",   // output directory name
        1,                                 // Interval of outputting
        5.0,                               // End time
        0.005,                             // Time step size
        box_indicator,                     // Object indicator
        constant_density,                  // Density distribution
        constant_velocity,                 // Velocity distribution
        8.0,                               // Total volume of the object
        0.1,                               // Radius of Poission disk sampling
        1.0,                               // Grid size
        Vector3<int>(11, 11, 11),          // Number of grid points in each
                                           // direction
        Vector3<int>(0, 0, 0),             // Bottom corner of the grid
        1.0,                               // Young's modulus
        0.499,                             // Poisson's ratio
        0.10,                              // Friction coefficient
        {0.0, 0.0, -9.81}                  // Gravitational acceleration
    };

    std::unique_ptr<MPMDriver> driver;
    driver.reset(new MPMDriver(param));

    BoundaryCondition::Boundary b0 = {param.mu, {{-1, 0, 0}, {9, 0, 0}}};
    BoundaryCondition::Boundary b1 = {param.mu, {{1, 0, 0}, {1, 0, 0}}};
    BoundaryCondition::Boundary b2 = {param.mu, {{0, -1, 0}, {0, 9, 0}}};
    BoundaryCondition::Boundary b3 = {param.mu, {{0, 1, 0}, {0, 1, 0}}};
    BoundaryCondition::Boundary b4 = {param.mu, {{0, 0, -1}, {0, 0, 9}}};
    BoundaryCondition::Boundary b5 = {param.mu, {{0, -1, 10}, {0, 0, 1}}};
    std::vector<BoundaryCondition::Boundary> boundaries =
                                                    {b0, b1, b2, b3, b4, b5};

    driver->InitializeBoundaryConditions(boundaries);
    driver->Run();

    return 0;
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake

int main() {
    return drake::multibody::mpm::DoMain();
}
