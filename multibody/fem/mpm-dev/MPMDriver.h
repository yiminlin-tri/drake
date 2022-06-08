#pragma once

#pragma once

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/mpm-dev/CorotatedModel.h"
#include "drake/multibody/fem/mpm-dev/GravitationalForce.h"
#include "drake/multibody/fem/mpm-dev/Grid.h"
#include "drake/multibody/fem/mpm-dev/MPMTransfer.h"
#include "drake/multibody/fem/mpm-dev/Particles.h"
#include "drake/multibody/fem/mpm-dev/particles_to_bgeo.h"
#include "drake/multibody/fem/mpm-dev/poisson_disk_sampling.h"

namespace drake {
namespace multibody {
namespace mpm {

class MPMDriver {
 public:
    // Indicator function for an object, true,  if position is in the object
    //                                   false, otherwise
    typedef bool (*IndicatorFunction)(const Vector3<double>& position);
    // Return the initial density of the object at position
    typedef double (*DensityField)(const Vector3<double>& position);
    // Return the initial velocity of the object at position
    typedef Vector3<double> (*VelocityField)(const Vector3<double>& position);

    struct MPMParameters {
        // IO parameters. Write output as output_directory/case_name($step).bgeo
        // every write_interval steps.
        std::string case_name;
        std::string output_directory;
        int write_interval;
        // Run the simulation with timestep size dt till endtime
        double endtime;
        double dt;
        // Particles parameters, determines particles' positions, densities, and
        // velocities
        IndicatorFunction object_indicator;
        DensityField density_field;
        VelocityField velocity_field;
        double total_volume;
        double sample_r;
        std::array<double, 3> x_min;                 // Specify the bounding box
        std::array<double, 3> x_max;                 // of initial particles
        // Grid parameters, as documented in Grid Class
        double h;
        Vector3<int> num_gridpt_1D;
        Vector3<int> bottom_corner;
        // Constitutive model Parameters, as documented in CorotatedModel
        double E;                            // Young's modulus
        double nu;                           // Poisson's ratio
        // Boundary condition parameters
        double mu;                           // Friction coefficient
        // Gravitational acceleration
        Vector3<double> g;
    };

    explicit MPMDriver(MPMParameters param);

    // Initialize the boundary condition by the boundaries information
    void InitializeBoundaryConditions(std::vector<BoundaryCondition::Boundary>
                                                                    boundaries);

    void Run();

 private:
    friend class MPMDriverTest;

    // Initialize particles' positions with Poisson disk sampling. Then
    // assuming every particles have equal reference volumes, initialize
    // particles' masses with the density distribution given. Finally,
    // Initialize the initial velocities of particles
    void InitializeParticles();

    // Advance MPM simulation by a single time step. Assuming both grid and
    // particles' state are at time n, a single time step involves a P2G
    // transfer, grid velocities update, a G2P transfer, and particles'
    // velocities update.
    void AdvanceOneTimeStep();

    // For every write_interval, write particles' information to
    // output_directory/case_name($step).bgeo
    void WriteParticlesToBgeo(int step);

    // Initialize particles and grid
    void SetUpDriver();

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
