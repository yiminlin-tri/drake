#include "drake/multibody/fem/mpm-dev/MPMDriver.h"

namespace drake {
namespace multibody {
namespace mpm {

MPMDriver::MPMDriver(MPMParameters param):
                                param_(std::move(param)),
                                grid_(param.solver_param.num_gridpt_1D,
                                      param.solver_param.h,
                                      param.solver_param.bottom_corner),
                                gravitational_force_(param.physical_param.g) {
    DRAKE_DEMAND(param.solver_param.endtime >= 0.0);
    DRAKE_DEMAND(param.solver_param.dt > 0.0);
}

void MPMDriver::InitializeKinematicCollisionObjects(KinematicCollisionObjects
                                                    objects) {
    collision_objects_ = std::move(objects);
}

void MPMDriver::DoTimeStepping() {
    int step = 0;
    int io_step = 0;
    double endtime = param_.solver_param.endtime;
    double dt = param_.solver_param.dt;
    bool is_last_step = false;
    using Clock = std::chrono::steady_clock;
    using Duration = std::chrono::duration<double>;
    std::chrono::time_point<Clock> start_timestepping_time = Clock::now();
    std::chrono::time_point<Clock> start_io_time;

    // Advance time steps until endtime
    for (double t = 0; t < endtime; t += dt) {
        is_last_step = (endtime - t <= dt);
        // If at the final timestep, modify the timestep size to match endtime
        if (endtime - t < dt) {
            dt = endtime - t;
        }
        // Check the CFL condition
        checkCFL(dt);
        AdvanceOneTimeStep(dt);
        step++;
        if ((t >= io_step*param_.io_param.write_interval) || (is_last_step)) {
            std::cout << "==== MPM Step " << step << " iostep " << io_step <<
                         ", at t = " << t+dt << std::endl;
            start_io_time = Clock::now();
            WriteParticlesToBgeo(io_step++);
            const auto elapsed_time
                    = std::chrono::duration_cast<Duration>(Clock::now()
                                                    - start_io_time).count();
            run_time_statistics_.time_io += elapsed_time;
        }
    }

    const auto total_time = std::chrono::duration_cast<Duration>(Clock::now()
                                            - start_timestepping_time).count();
    run_time_statistics_.time_total += total_time;

    PrintRunTimeStatistics();
}

// TODO(yiminlin.tri): Not tested.
void MPMDriver::InitializeParticles(const AnalyticLevelSet& level_set,
                                    const math::RigidTransform<double>& pose,
                                    MaterialParameters m_param) {
    DRAKE_DEMAND(m_param.density > 0.0);
    DRAKE_DEMAND(m_param.min_num_particles_per_cell >= 1);

    double h = param_.solver_param.h;
    const std::array<Vector3<double>, 2> bounding_box =
                                    level_set.get_bounding_box();

    // Distances between generated particles are at at least sample_r apart, and
    // there are at least min_num_particles_per_cell particles per cell. r =
    // argmax(⌊h/r⌋)^3 >= min_num_particles_per_cell, in other words, if we pick
    // particles located at the grid with grid size r, there are at least
    // min_num_particles_per_cell particles in a cell with size h.
    double sample_r =
            h/(std::cbrt(m_param.min_num_particles_per_cell)+1);
    multibody::SpatialVelocity<double> init_v = m_param.initial_velocity;
    std::array<double, 3> xmin = {bounding_box[0][0], bounding_box[0][1],
                                  bounding_box[0][2]};
    std::array<double, 3> xmax = {bounding_box[1][0], bounding_box[1][1],
                                  bounding_box[1][2]};
    // Generate sample particles in the reference frame (centered at the origin
    // with canonical basis e_i)
    std::vector<Vector3<double>> particles_sample_positions =
        thinks::PoissonDiskSampling<double, 3, Vector3<double>>(sample_r,
                                                                xmin, xmax);

    // Pick out sampled particles that are in the object
    int num_samples = particles_sample_positions.size();
    std::vector<Vector3<double>> particles_positions, particles_velocities;
    for (int p = 0; p < num_samples; ++p) {
        // Denote the particle by P, the object by B, the frame centered at the
        // origin of the object by Bo, and the frame of the particle by Bp.

        // The pose and spatial velocity of the object in the world frame
        const math::RigidTransform<double>& X_WB       = pose;
        const multibody::SpatialVelocity<double>& V_WB = init_v;

        // Rotation matrix from the object to the world
        const math::RotationMatrix<double>& Rot_WB     = X_WB.rotation();

        // Sample particles and get the position of the particle with respect to
        // the object in the object's frame
        const Vector3<double>& p_BoBp_B = particles_sample_positions[p];

        // Express the relative position of the particle with respect to the
        // object in the world frame
        const Vector3<double>& p_BoBp_W = Rot_WB*p_BoBp_B;

        // Compute the spatial velocity of the particle
        SpatialVelocity<double> V_WBp   = V_WB.Shift(p_BoBp_W);

        // If the particle is in the level set
        if (level_set.InInterior(p_BoBp_B)) {
            // TODO(yiminlin.tri): Initialize the affine matrix C_p using
            //                     V_WBp.rotational() ?
            particles_velocities.emplace_back(V_WBp.translational());
            // Place the particle's position in world frame
            particles_positions.emplace_back(X_WB*p_BoBp_B);
        }
    }

    int num_particles = particles_positions.size();
    // We assume every particle have the same volume and mass
    double reference_volume_p = level_set.get_volume()/num_particles;
    double init_m = m_param.density*reference_volume_p;

    // Add particles
    for (int p = 0; p < num_particles; ++p) {
        const Vector3<double>& xp = particles_positions[p];
        const Vector3<double>& vp = particles_velocities[p];
        Matrix3<double> elastic_deformation_grad_p
                                                = Matrix3<double>::Identity();
        Matrix3<double> kirchhoff_stress_p = Matrix3<double>::Identity();
        Matrix3<double> B_p                = Matrix3<double>::Zero();
        std::unique_ptr<ElastoPlasticModel> elastoplastic_model_p
                                        = m_param.elastoplastic_model->Clone();
        particles_.AddParticle(xp, vp, init_m, reference_volume_p,
                               elastic_deformation_grad_p,
                               kirchhoff_stress_p,
                               B_p, std::move(elastoplastic_model_p));
    }
}

void MPMDriver::WriteParticlesToBgeo(int io_step) {
    std::string output_filename = param_.io_param.output_directory + "/"
                                + param_.io_param.case_name
                                + std::to_string(io_step) + ".bgeo";
    internal::WriteParticlesToBgeo(output_filename, particles_.get_positions(),
                                                    particles_.get_velocities(),
                                                    particles_.get_masses());
}

void MPMDriver::PrintRunTimeStatistics() const {
    std::cout << "===================" << std::endl;
    std::cout << "==== Total run                   time: "
            << run_time_statistics_.time_total
            << " seconds" << std::endl;
    std::cout << "==== Total IO                    time: "
            << run_time_statistics_.time_io
            << " seconds" << std::endl;
    std::cout << "==== Total stress/plastic update time: "
            << run_time_statistics_.time_update_stress_and_plasticity
            << " seconds" << std::endl;
    std::cout << "==== Total transfer setup        time: "
            << run_time_statistics_.time_setup_transfer
            << " seconds" << std::endl;
    std::cout << "==== Total P2G                   time: "
            << run_time_statistics_.time_P2G
            << " seconds" << std::endl;
    std::cout << "==== Total update velocity       time: "
            << run_time_statistics_.time_update_velocity
            << " seconds" << std::endl;
    std::cout << "==== Total apply forces          time: "
            << run_time_statistics_.time_apply_forces
            << " seconds" << std::endl;
    std::cout << "==== Total enfroce BC            time: "
            << run_time_statistics_.time_enforce_bc
            << " seconds" << std::endl;
    std::cout << "==== Total G2P                   time: "
            << run_time_statistics_.time_G2P
            << " seconds" << std::endl;
    std::cout << "==== Total particles advection   time: "
            << run_time_statistics_.time_advect_particles
            << " seconds" << std::endl;
}

void MPMDriver::checkCFL(double dt) {
    for (const auto& v : particles_.get_velocities()) {
        if (!(std::max({std::abs(dt*v(0)),
                        std::abs(dt*v(1)),
                        std::abs(dt*v(2))}) <= grid_.get_h())) {
            throw std::runtime_error("CFL condition violation");
        }
    }
}

void MPMDriver::AdvanceOneTimeStep(double dt) {
    using Clock = std::chrono::steady_clock;
    using Duration = std::chrono::duration<double>;
    std::chrono::time_point<Clock> start_time;
    Duration::rep elapsed_time;

    // Apply plasticity and update Kirchhoff on particles
    start_time = Clock::now();
    particles_.ApplyPlasticityAndUpdateKirchhoffStresses();
    elapsed_time = std::chrono::duration_cast<Duration>(Clock::now()
                                                      - start_time).count();
    run_time_statistics_.time_update_stress_and_plasticity += elapsed_time;

    // Set up the transfer routines (Preallocations, sort the particles)
    start_time = Clock::now();
    mpm_transfer_.SetUpTransfer(grid_, &particles_);
    elapsed_time = std::chrono::duration_cast<Duration>(Clock::now()
                                                      - start_time).count();
    run_time_statistics_.time_setup_transfer += elapsed_time;

    // Main Algorithm:
    // P2G
    start_time = Clock::now();
    mpm_transfer_.TransferParticlesToGrid(particles_, &grid_);
    elapsed_time = std::chrono::duration_cast<Duration>(Clock::now()
                                                      - start_time).count();
    run_time_statistics_.time_P2G += elapsed_time;

    // Update grid velocity
    start_time = Clock::now();
    grid_.UpdateVelocity(dt);
    elapsed_time = std::chrono::duration_cast<Duration>(Clock::now()
                                                      - start_time).count();
    run_time_statistics_.time_update_velocity += elapsed_time;

    // Apply gravitational force and enforce boundary condition
    start_time = Clock::now();
    gravitational_force_.ApplyGravitationalForces(dt, &grid_);
    elapsed_time = std::chrono::duration_cast<Duration>(Clock::now()
                                                      - start_time).count();
    run_time_statistics_.time_apply_forces += elapsed_time;

    // Update Collision Objects
    start_time = Clock::now();
    collision_objects_.AdvanceOneTimeStep(dt);
    grid_.EnforceBoundaryCondition(collision_objects_);
    elapsed_time = std::chrono::duration_cast<Duration>(Clock::now()
                                                      - start_time).count();
    run_time_statistics_.time_enforce_bc += elapsed_time;

    // G2P
    start_time = Clock::now();
    mpm_transfer_.TransferGridToParticles(grid_, dt, &particles_);
    elapsed_time = std::chrono::duration_cast<Duration>(Clock::now()
                                                      - start_time).count();
    run_time_statistics_.time_G2P += elapsed_time;

    // Advect particles
    start_time = Clock::now();
    particles_.AdvectParticles(dt);
    elapsed_time = std::chrono::duration_cast<Duration>(Clock::now()
                                                      - start_time).count();
    run_time_statistics_.time_advect_particles += elapsed_time;
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
