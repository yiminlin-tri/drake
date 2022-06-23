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
#include "drake/math/roll_pitch_yaw.h"
#include "drake/multibody/fem/mpm-dev/MPMDriver.h"
#include "drake/multibody/math/spatial_velocity.h"

namespace drake {
namespace multibody {
namespace mpm {

int DoMain() {
    MPMParameters::PhysicalParameters p_param {
        {0.0, 0.0, -9.81}                      // Gravitational acceleration
    };

    MPMParameters::SolverParameters s_param {
        1e-0,                                  // End time
        1e-4,                                  // Time step size
        0.01,                                  // Grid size
        Vector3<int>(23, 23, 23),              // Number of grid points in each
                                               // direction
        Vector3<int>(-1, -1, -1),              // Bottom corner of the grid
    };

    MPMParameters::IOParameters io_param {
        "mpm-test",                            // case name
        "/home/yiminlin/Desktop/output",       // output directory name
        0.04,                                  // Interval of outputting
    };

    MPMParameters param {p_param, s_param, io_param};
    auto driver = std::make_unique<MPMDriver>(std::move(param));

    KinematicCollisionObjects objects = KinematicCollisionObjects();

    // Initialize the wall
    multibody::SpatialVelocity<double> wall_velocity;
    wall_velocity.SetZero();
    double wall_mu = 0.8;
    Vector3<double> wall_normal = {0.0, 0.0, 1.0};
    std::unique_ptr<AnalyticLevelSet> wall_level_set =
                            std::make_unique<HalfSpaceLevelSet>(wall_normal);
    Vector3<double> wall_translation = {0.0, 0.0, 0.0};
    math::RigidTransform<double> wall_pose =
                            math::RigidTransform<double>(wall_translation);
    objects.AddCollisionObject(std::move(wall_level_set), wall_pose,
                               wall_velocity, wall_mu);

    // Initialize the cylinder
    multibody::SpatialVelocity<double> cylinder_velocity;
    cylinder_velocity.translational() = Vector3<double>(-0.1, 0.0, 0.0);
    cylinder_velocity.rotational() = Vector3<double>(0.0, 0.0, -M_PI);
    double cylinder_mu = 0.7;
    double cylinder_height = 0.2;
    double cylinder_radius = 0.04;
    std::unique_ptr<AnalyticLevelSet> cylinder_level_set =
                            std::make_unique<CylinderLevelSet>(cylinder_height,
                                                               cylinder_radius);
    Vector3<double> cylinder_translation = {0.25, 0.1, 0.05};
    math::RigidTransform<double> cylinder_pose =
                            math::RigidTransform<double>(cylinder_translation);
    objects.AddCollisionObject(std::move(cylinder_level_set), cylinder_pose,
                               cylinder_velocity, cylinder_mu);

    // Initialize a sphere
    double radius = 0.02;
    SphereLevelSet level_set_sphere = SphereLevelSet(radius);
    Vector3<double> translation_sphere = {0.18, 0.1, 0.03};
    math::RigidTransform<double> pose_sphere =
                            math::RigidTransform<double>(translation_sphere);
    MPMDriver::MaterialParameters m_param_sphere { {8e1, 0.49},
                                                   800,
                                                   {0.0, 0.0, 0.0},
                                                   1
                                                 };

    driver->InitializeKinematicCollisionObjects(std::move(objects));
    driver->InitializeParticles(level_set_sphere, pose_sphere,
                                std::move(m_param_sphere));
    driver->DoTimeStepping();

    return 0;
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake

int main() {
    return drake::multibody::mpm::DoMain();
}
