#include "drake/multibody/fem/mpm-dev/MPMDriver.h"

namespace drake {
namespace multibody {
namespace mpm {

MPMDriver::MPMDriver(MPMParameters param) {
    param_ = param;
}

void MPMDriver::InitializeBoundaryConditions(
                        std::vector<BoundaryCondition::Boundary> boundaries) {
    boundary_condition_ = BoundaryCondition(boundaries);
}

void MPMDriver::Run() {
    SetUpDriver();
    DoTimeStepping();
}

void MPMDriver::SetUpDriver() {
    MPMDriver::InitializeParticles();
    grid_ = Grid(param_.num_gridpt_1D, param_.h, param_.bottom_corner);
}

void MPMDriver::DoTimeStepping() {
    int step = 0;
    // Advance time steps until endtime
    for (double t = 0; t < param_.endtime; t += param_.dt) {
        // If at the final timestep, modify the timestep size to match endtime
        if (param_.endtime - t < param_.dt) {
            param_.dt = param_.endtime - t;
        }
        AdvanceOneTimeStep();
        step++;
        if (step % param_.write_interval == 0) {
            std::cout << "==== MPM Step " << step <<
                         ", at t = " << t << std::endl;
            WriteParticlesToBgeo(step);
        }
    }
}

void MPMDriver::InitializeParticles() {
    double h = param_.h;
    std::array<double, 3> x_min = {{param_.bottom_corner(0)*h,
                                   param_.bottom_corner(1)*h,
                                   param_.bottom_corner(2)*h}};
    std::array<double, 3> x_max =
              {{(param_.bottom_corner(0)+param_.num_gridpt_1D(0)-1)*h,
                (param_.bottom_corner(1)+param_.num_gridpt_1D(1)-1)*h,
                (param_.bottom_corner(2)+param_.num_gridpt_1D(2)-1)*h}};
    std::vector<Vector3<double>> particles_sample_positions =
        thinks::PoissonDiskSampling<double, 3, Vector3<double>>(param_.sample_r,
                                                                x_min, x_max);

    // Pick out sampled particles that are in the object
    int num_samples = particles_sample_positions.size();
    std::vector<Vector3<double>> particles_positions;
    for (int p = 0; p < num_samples; ++p) {
        const Vector3<double>& xp = particles_sample_positions[p];
        if (param_.object_indicator(xp)) {
            particles_positions.emplace_back(xp);
        }
    }

    int num_particles = particles_positions.size();
    // We assume every particle have the same volume
    double reference_volume_p = param_.total_volume/num_particles;

    // Make an empty particles class
    particles_ = Particles();

    // Add particles
    for (int p = 0; p < num_particles; ++p) {
        const Vector3<double>& xp = particles_positions[p];
        particles_.AddParticle(xp, param_.velocity_field(xp),
                                   param_.density_field(xp)*reference_volume_p,
                                   reference_volume_p,
                                   Matrix3<double>::Identity(),
                                   Matrix3<double>::Identity());
    }
}

void MPMDriver::WriteParticlesToBgeo(int step) {
    std::string output_filename = param_.output_directory + "/"
                                + param_.case_name + std::to_string(step)
                                + ".bgeo";
    internal::WriteParticlesToBgeo(output_filename, particles_.get_positions(),
                                                    particles_.get_velocities(),
                                                    particles_.get_masses());
}

void MPMDriver::AdvanceOneTimeStep() {
    CorotatedModel corotated_model = CorotatedModel(param_.E, param_.nu);
    GravitationalForce gravitational_force = GravitationalForce(param_.g);

    // Set up the transfer routines (Preallocations, sort the particles)
    MPMTransfer mpm_transfer = MPMTransfer();
    mpm_transfer.SetUpTransfer(grid_, &particles_);

    // Main Algorithm:
    // P2G
    mpm_transfer.TransferParticlesToGrid(particles_, &grid_);

    // Update grid velocity
    grid_.UpdateVelocity(param_.dt);

    // Apply gravitational force and enforce boundary condition
    gravitational_force.ApplyGravitationalForces(param_.dt, &grid_);
    grid_.EnforceBoundaryCondition(boundary_condition_);

    // G2P
    mpm_transfer.TransferGridToParticles(grid_, param_.dt, &particles_);

    // Advect and update particles
    particles_.UpdateKirchhoffStresses(corotated_model);
    particles_.AdvectParticles(param_.dt);
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
