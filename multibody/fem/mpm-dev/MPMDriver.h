#pragma once

#pragma once

#include <array>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/mpm-dev/CorotatedModel.h"
#include "drake/multibody/fem/mpm-dev/GravitationalForce.h"
#include "drake/multibody/fem/mpm-dev/Grid.h"
#include "drake/multibody/fem/mpm-dev/MPMParameters.h"
#include "drake/multibody/fem/mpm-dev/MPMTransfer.h"
#include "drake/multibody/fem/mpm-dev/Particles.h"
#include "drake/multibody/fem/mpm-dev/particles_to_bgeo.h"
#include "drake/multibody/fem/mpm-dev/poisson_disk_sampling.h"

namespace drake {
namespace multibody {
namespace mpm {

class MPMDriver {
 public:
    explicit MPMDriver(MPMParameters param);

    // Initialize the boundary condition by the boundaries information
    void InitializeBoundaryConditions(std::vector<BoundaryCondition::Boundary>
                                                                    boundaries);

    // Run the MPM simulation with the object represented by the given level set
    // in the reference frame, then transformed into the physical frame by pose.
    void Run(const AnalyticLevelSet& level_set,
             const math::RigidTransform<double>& pose);

 private:
    friend class MPMDriverTest;

    // Initialize particles' positions with Poisson disk sampling. The object's
    // level set in the physical frame is the given level set in the reference
    // frame transformed by pose. We assume every particles have equal reference
    // volumes, then we can initialize particles' masses with the given constant
    // density, Finally, we initialize the velocities of particles with the
    // constant given velocity.
    void InitializeParticles(const AnalyticLevelSet& initial_level_set,
                             const math::RigidTransform<double>& pose);

    // Advance MPM simulation by a single time step. Assuming both grid and
    // particles' state are at time n, a single time step involves a P2G
    // transfer, grid velocities update, a G2P transfer, and particles'
    // velocities update.
    void AdvanceOneTimeStep();

    // For every write_interval, write particles' information to
    // output_directory/case_name($step).bgeo
    void WriteParticlesToBgeo(int step);

    // Initialize particles and grid
    void SetUpDriver(const AnalyticLevelSet& level_set,
                     const math::RigidTransform<double>& pose);

    // Symplectic Euler time stepping till endtime with dt
    void DoTimeStepping();

    MPMParameters param_;
    Particles particles_;
    Grid grid_;
    MPMTransfer mpm_transfer_;
    CorotatedModel corotated_model_;
    GravitationalForce gravitational_force_;
    BoundaryCondition boundary_condition_;
};  // class MPMDriver

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
