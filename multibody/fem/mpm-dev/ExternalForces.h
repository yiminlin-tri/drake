#pragma once

#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/mpm-dev/Grid.h"

namespace drake {
namespace multibody {
namespace mpm {

// A class providing the utitlities to apply external forces to the grid
class ExternalForces {
 public:
    ExternalForces() = default;

    // Apply the gravitational forces to grid points. (Equivalently, apply
    // gravitational acceleration to the velocity)
    // v_new = v_prev + dt*(0, 0, -g), where g = 9.81
    void ApplyGravitationalForces(double dt, Grid* grid) const;
};  // class ExternalForces

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
