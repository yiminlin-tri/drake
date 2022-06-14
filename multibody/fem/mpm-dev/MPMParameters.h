#pragma once

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/math/rigid_transform.h"
#include "drake/multibody/fem/mpm-dev/AnalyticLevelSet.h"

namespace drake {
namespace multibody {
namespace mpm {

// A class holding the parameters used in MPMDriver
class MPMParameters {
 public:
    struct PhysicalParameters {
        // Constitutive model Parameters, as documented in CorotatedModel
        double E;                            // Young's modulus
        double nu;                           // Poisson's ratio
        // Boundary condition parameters
        double mu;                           // Friction coefficient
        // Gravitational acceleration
        Vector3<double> g;
    };

    // @pre density is positive
    struct InitialConditionParameters {
        // Density and the initial velocity of the object, we assume the object
        // has uniform density and velocity.
        double density;
        Vector3<double> initial_velocity;
    };

    // @pre endtime nonnegative, dt ∈ (0, endtime]
    // @pre min_num_particles_per_cell >= 1
    struct SolverParameters {
        // Run the simulation with timestep size dt till endtime
        double endtime;
        double dt;
        // User defined parameter to control the minimum number of particles per
        // grid cell.
        int min_num_particles_per_cell;
        // Grid parameters, as documented in Grid Class
        double h;
        Vector3<int> num_gridpt_1D;
        Vector3<int> bottom_corner;
    };

    // IO parameters. Write the output as output_directory/case_name($i).bgeo
    // write_interval is the length of time between outputting, and i is id of
    // the output file, ordered ascendingly in time.
    struct IOParameters {
        std::string case_name;
        std::string output_directory;
        double write_interval;
    };

    MPMParameters(PhysicalParameters p_param,
                  InitialConditionParameters init_param,
                  SolverParameters s_param, IOParameters i_param):
                                physical_param(std::move(p_param)),
                                initial_condition_param(std::move(init_param)),
                                solver_param(std::move(s_param)),
                                io_param(std::move(i_param)) {}

    // Member variables
    PhysicalParameters physical_param;
    InitialConditionParameters initial_condition_param;
    SolverParameters solver_param;
    IOParameters io_param;
};  // class MPMParameters

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
