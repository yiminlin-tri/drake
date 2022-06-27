#include "drake/multibody/fem/mpm-dev/Particles.h"

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"
#include "drake/math/roll_pitch_yaw.h"
#include "drake/math/rotation_matrix.h"
#include "drake/multibody/fem/mpm-dev/CorotatedModel.h"

namespace drake {
namespace multibody {
namespace mpm {
namespace internal {
namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double kEps = 1e-14;
constexpr double TOLERANCE = 1e-10;

GTEST_TEST(ParticlesClassTest, TestAddSetGet) {
    std::vector<Vector3<double>> positions;
    std::vector<Vector3<double>> velocities;
    std::vector<double> masses;
    std::vector<double> reference_volumes;
    std::vector<Matrix3<double>> deformation_gradients;
    std::vector<Matrix3<double>> kirchhoff_stresses;
    std::vector<Matrix3<double>> affine_matrices;
    std::vector<CorotatedModel> corotated_models;

    Vector3<double> pos1 = {1.0, 2.0, 3.0};
    Vector3<double> vel1 = {-1.0, -2.0, -3.0};
    double mass1 = 5.0;
    double vol1  = 10.0;
    Matrix3<double> F1 = pos1.asDiagonal();
    Matrix3<double> stress1 = vel1.asDiagonal();
    Matrix3<double> B1 = 2.0*vel1.asDiagonal();
    CorotatedModel cmodel1 = CorotatedModel(10.0, 0.1);

    Vector3<double> pos2 = {3.0, -1.0, 6.0};
    Vector3<double> vel2 = {-9.0, 8.0, -2.0};
    double mass2 = 7.0;
    double vol2  = 3.0;
    Matrix3<double> F2 = pos2.asDiagonal();
    Matrix3<double> stress2 = vel2.asDiagonal();
    Matrix3<double> B2 = 2.0*vel2.asDiagonal();
    CorotatedModel cmodel2 = CorotatedModel(20.0, 0.2);

    Particles particles = Particles();
    EXPECT_EQ(particles.get_num_particles(), 0);
    particles.AddParticle(pos1, vel1, mass1, vol1, F1, stress1, B1, cmodel1);
    EXPECT_EQ(particles.get_num_particles(), 1);
    particles.AddParticle(pos2, vel2, mass2, vol2, F2, stress2, B2, cmodel2);
    EXPECT_EQ(particles.get_num_particles(), 2);

    // Test get individual element
    EXPECT_TRUE(CompareMatrices(particles.get_position(0), pos1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(0), vel1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass(0), mass1);
    EXPECT_EQ(particles.get_reference_volume(0), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(0), F1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(0), stress1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(0), B1,
                std::numeric_limits<double>::epsilon()));

    EXPECT_TRUE(CompareMatrices(particles.get_position(1), pos2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(1), vel2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass(1), mass2);
    EXPECT_EQ(particles.get_reference_volume(1), vol2);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(1), F2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(1), stress2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(1), B2,
                std::numeric_limits<double>::epsilon()));

    // Test get vectors
    positions             = particles.get_positions();
    velocities            = particles.get_velocities();
    masses                = particles.get_masses();
    reference_volumes     = particles.get_reference_volumes();
    deformation_gradients = particles.get_deformation_gradients();
    kirchhoff_stresses    = particles.get_kirchhoff_stresses();
    affine_matrices       = particles.get_affine_matrices();
    EXPECT_TRUE(CompareMatrices(positions[0], pos1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(velocities[0], vel1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(masses[0], mass1);
    EXPECT_EQ(reference_volumes[0], vol1);
    EXPECT_TRUE(CompareMatrices(deformation_gradients[0], F1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(kirchhoff_stresses[0], stress1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(affine_matrices[0], B1,
                std::numeric_limits<double>::epsilon()));

    EXPECT_TRUE(CompareMatrices(positions[1], pos2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(velocities[1], vel2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(masses[1], mass2);
    EXPECT_EQ(reference_volumes[1], vol2);
    EXPECT_TRUE(CompareMatrices(deformation_gradients[1], F2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(kirchhoff_stresses[1], stress2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(affine_matrices[1], B2,
                std::numeric_limits<double>::epsilon()));

    // Test set individual element
    particles = Particles(2);
    EXPECT_EQ(particles.get_num_particles(), 2);
    particles.set_position(0, pos1);
    particles.set_velocity(0, vel1);
    particles.set_mass(0, mass1);
    particles.set_reference_volume(0, vol1);
    particles.set_deformation_gradient(0, F1);
    particles.set_kirchhoff_stress(0, stress1);
    particles.set_affine_matrix(0, B1);
    particles.set_position(1, pos2);
    particles.set_velocity(1, vel2);
    particles.set_mass(1, mass2);
    particles.set_reference_volume(1, vol2);
    particles.set_deformation_gradient(1, F2);
    particles.set_kirchhoff_stress(1, stress2);
    particles.set_affine_matrix(1, B2);

    EXPECT_TRUE(CompareMatrices(particles.get_position(0), pos1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(0), vel1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass(0), mass1);
    EXPECT_EQ(particles.get_reference_volume(0), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(0), F1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(0), stress1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(0), B1,
                std::numeric_limits<double>::epsilon()));

    EXPECT_TRUE(CompareMatrices(particles.get_position(1), pos2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(1), vel2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass(1), mass2);
    EXPECT_EQ(particles.get_reference_volume(1), vol2);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(1), F2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(1), stress2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(1), B2,
                std::numeric_limits<double>::epsilon()));

    particles.AddParticle(pos1, vel1, mass1, vol1, F1, stress1, B1, cmodel1);
    EXPECT_EQ(particles.get_num_particles(), 3);
    EXPECT_TRUE(CompareMatrices(particles.get_position(2), pos1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(2), vel1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass(2), mass1);
    EXPECT_EQ(particles.get_reference_volume(2), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(2), F1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(2), stress1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(2), B1,
                std::numeric_limits<double>::epsilon()));

    // Test set vector
    particles = Particles(2);
    EXPECT_EQ(particles.get_num_particles(), 2);
    particles.set_positions(positions);
    particles.set_velocities(velocities);
    particles.set_masses(masses);
    particles.set_reference_volumes(reference_volumes);
    particles.set_deformation_gradients(deformation_gradients);
    particles.set_kirchhoff_stresses(kirchhoff_stresses);
    particles.set_affine_matrices(affine_matrices);

    EXPECT_TRUE(CompareMatrices(particles.get_position(0), pos1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(0), vel1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass(0), mass1);
    EXPECT_EQ(particles.get_reference_volume(0), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(0), F1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(0), stress1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(0), B1,
                std::numeric_limits<double>::epsilon()));

    EXPECT_TRUE(CompareMatrices(particles.get_position(1), pos2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(1), vel2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass(1), mass2);
    EXPECT_EQ(particles.get_reference_volume(1), vol2);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(1), F2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(1), stress2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(1), B2,
                std::numeric_limits<double>::epsilon()));

    particles.AddParticle(pos1, vel1, mass1, vol1, F1, stress1, B1, cmodel1);
    EXPECT_EQ(particles.get_num_particles(), 3);
    EXPECT_TRUE(CompareMatrices(particles.get_position(2), pos1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(2), vel1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass(2), mass1);
    EXPECT_EQ(particles.get_reference_volume(2), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(2), F1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(2), stress1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(2), B1,
                std::numeric_limits<double>::epsilon()));
}

GTEST_TEST(ParticlesClassTest, TestReorder) {
    std::vector<Vector3<double>> positions;
    std::vector<Vector3<double>> velocities;
    std::vector<double> masses;
    std::vector<double> reference_volumes;
    std::vector<Matrix3<double>> deformation_gradients;
    std::vector<Matrix3<double>> kirchhoff_stresses;

    Vector3<double> pos1 = {1.0, 2.0, 3.0};
    Vector3<double> vel1 = {-1.0, -2.0, -3.0};
    double mass1 = 5.0;
    double vol1  = 10.0;
    Matrix3<double> F1 = pos1.asDiagonal();
    Matrix3<double> stress1 = vel1.asDiagonal();
    Matrix3<double> B1 = 2.0*vel1.asDiagonal();
    CorotatedModel cmodel1 = CorotatedModel(10.0, 0.1);

    Vector3<double> pos2 = {3.0, -1.0, 6.0};
    Vector3<double> vel2 = {-9.0, 8.0, -2.0};
    double mass2 = 7.0;
    double vol2  = 3.0;
    Matrix3<double> F2 = pos2.asDiagonal();
    Matrix3<double> stress2 = vel2.asDiagonal();
    Matrix3<double> B2 = 2.0*vel2.asDiagonal();
    CorotatedModel cmodel2 = CorotatedModel(20.0, 0.2);

    Vector3<double> pos3 = {3.2, -1.0, 1.0};
    Vector3<double> vel3 = {2.0, -6.2, 8.0};
    double mass3 = 2.0;
    double vol3  = 12.0;
    Matrix3<double> F3 = pos3.asDiagonal();
    Matrix3<double> stress3 = vel3.asDiagonal();
    Matrix3<double> B3 = 2.0*vel3.asDiagonal();
    CorotatedModel cmodel3 = CorotatedModel(30.0, 0.3);

    Particles particles = Particles();
    particles.AddParticle(pos1, vel1, mass1, vol1, F1, stress1, B1, cmodel1);
    particles.AddParticle(pos2, vel2, mass2, vol2, F2, stress2, B2, cmodel2);
    particles.AddParticle(pos3, vel3, mass3, vol3, F3, stress3, B3, cmodel3);

    // Check the original ordering
    EXPECT_TRUE(CompareMatrices(particles.get_position(0), pos1, kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(0), vel1, kEps));
    EXPECT_EQ(particles.get_mass(0), mass1);
    EXPECT_EQ(particles.get_reference_volume(0), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(0), F1,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(0), stress1,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(0), B1,
                kEps));

    EXPECT_TRUE(CompareMatrices(particles.get_position(1), pos2,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(1), vel2,
                kEps));
    EXPECT_EQ(particles.get_mass(1), mass2);
    EXPECT_EQ(particles.get_reference_volume(1), vol2);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(1), F2,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(1), stress2,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(1), B2,
                kEps));

    EXPECT_TRUE(CompareMatrices(particles.get_position(2), pos3,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(2), vel3,
                kEps));
    EXPECT_EQ(particles.get_mass(2), mass3);
    EXPECT_EQ(particles.get_reference_volume(2), vol3);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(2), F3,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(2), stress3,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(2), B3,
                kEps));

    std::vector<size_t> new_order{2, 0, 1};
    particles.Reorder(new_order);
    EXPECT_TRUE(CompareMatrices(particles.get_position(0), pos3,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(0), vel3,
                kEps));
    EXPECT_EQ(particles.get_mass(0), mass3);
    EXPECT_EQ(particles.get_reference_volume(0), vol3);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(0), F3,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(0), stress3,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(0), B3,
                kEps));

    EXPECT_TRUE(CompareMatrices(particles.get_position(1), pos1,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(1), vel1,
                kEps));
    EXPECT_EQ(particles.get_mass(1), mass1);
    EXPECT_EQ(particles.get_reference_volume(1), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(1), F1,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(1), stress1,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(1), B1,
                kEps));

    EXPECT_TRUE(CompareMatrices(particles.get_position(2), pos2,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity(2), vel2,
                kEps));
    EXPECT_EQ(particles.get_mass(2), mass2);
    EXPECT_EQ(particles.get_reference_volume(2), vol2);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient(2), F2,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(2), stress2,
                kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_affine_matrix(2), B2,
                kEps));
}

GTEST_TEST(ParticlesClassTest, TestAdvectAndUpdateKirchhoffStress) {
    std::vector<Vector3<double>> positions;
    std::vector<Vector3<double>> velocities;
    std::vector<double> masses;
    std::vector<double> reference_volumes;
    std::vector<Matrix3<double>> deformation_gradients;
    std::vector<Matrix3<double>> kirchhoff_stresses;
    CorotatedModel coro_model = CorotatedModel(5.0, 0.25);

    const Matrix3<double> R =
        math::RotationMatrix<double>
            (math::RollPitchYaw<double>(pi/2.0, pi, 3.0*pi/2.0)).matrix();
    const Matrix3<double> S = (Matrix3<double>() <<
            6, 1, 2,
            1, 4, 1,
            2, 1, 5).finished();
    const Matrix3<double> tau_exact = (Matrix3<double>() <<
             18724, -84,  40,
            -84,  18764, -44,
             40, -44,  18680).finished();
    const Matrix3<double> F = R * S;

    Vector3<double> pos1 = {1.0, 2.0, 3.0};
    Vector3<double> vel1 = {-1.0, -2.0, -3.0};
    double mass1 = 5.0;
    double vol1  = 10.0;
    Matrix3<double> F1 = F;
    Matrix3<double> stress1 = Matrix3<double>::Zero();
    Matrix3<double> B1 = Matrix3<double>::Zero();

    Vector3<double> pos2 = {3.0, -1.0, 6.0};
    Vector3<double> vel2 = {-9.0, 8.0, -2.0};
    double mass2 = 7.0;
    double vol2  = 3.0;
    Matrix3<double> F2 = R;
    Matrix3<double> stress2 = Matrix3<double>::Ones();
    Matrix3<double> B2 = Matrix3<double>::Ones();

    Particles particles = Particles();
    particles.AddParticle(pos1, vel1, mass1, vol1, F1, stress1, B1, coro_model);
    particles.AddParticle(pos2, vel2, mass2, vol2, F2, stress2, B2, coro_model);

    // Advect particles
    double dt = 0.3;
    particles.AdvectParticles(dt);
    EXPECT_TRUE(CompareMatrices(particles.get_position(0),
                                Vector3<double>(0.7, 1.4, 2.1), kEps));
    EXPECT_TRUE(CompareMatrices(particles.get_position(1),
                                Vector3<double>(0.3, 1.4, 5.4), kEps));

    // Update Kirchhoff stress
    particles.UpdateKirchhoffStresses();
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(0),
                                tau_exact, TOLERANCE));
    EXPECT_TRUE(CompareMatrices(particles.get_kirchhoff_stress(1),
                                Matrix3<double>::Zero(), TOLERANCE));
}

GTEST_TEST(GridClassTest, TestGridSumState) {
    // In this test case, we construct a 2x2x2 grid
    Vector3<double> pos0 = {1.0, 2.0, 3.0};
    Vector3<double> vel0 = {-1.0, -1.0, -1.0};
    double mass0 = 1.0;
    double vol0  = 1.0;
    Matrix3<double> F0 = pos0.asDiagonal();
    Matrix3<double> stress0 = vel0.asDiagonal();
    Matrix3<double> B0 = vel0.asDiagonal();
    CorotatedModel cmodel0 = CorotatedModel(10.0, 0.1);

    Vector3<double> pos1 = {4.0, 5.0, 6.0};
    Vector3<double> vel1 = {-2.0, -2.0, -2.0};
    double mass1 = 2.0;
    double vol1  = 2.0;
    Matrix3<double> F1 = pos1.asDiagonal();
    Matrix3<double> stress1 = vel1.asDiagonal();
    Matrix3<double> B1 = 2.0*vel1.asDiagonal();
    CorotatedModel cmodel1 = CorotatedModel(10.0, 0.1);

    Vector3<double> pos2 = {7.0, 8.0, 9.0};
    Vector3<double> vel2 = {-3.0, -3.0, -3.0};
    double mass2 = 3.0;
    double vol2  = 3.0;
    Matrix3<double> F2 = pos2.asDiagonal();
    Matrix3<double> stress2 = vel2.asDiagonal();
    Matrix3<double> B2 = 3.0*vel2.asDiagonal();
    CorotatedModel cmodel2 = CorotatedModel(20.0, 0.2);

    Particles particles = Particles();
    particles.AddParticle(pos0, vel0, mass0, vol0, F0, stress0, B0, cmodel0);
    particles.AddParticle(pos1, vel1, mass1, vol1, F1, stress1, B1, cmodel1);
    particles.AddParticle(pos2, vel2, mass2, vol2, F2, stress2, B2, cmodel2);

    Particles::ParticlesSumState sum_state = particles.GetParticlesSumState();

    // Check the sums
    EXPECT_EQ(sum_state.sum_mass, 6.0);
    EXPECT_TRUE(CompareMatrices(sum_state.sum_momentum,
                                Vector3<double>(-14.0, -14.0, -14.0), kEps));
    EXPECT_TRUE(CompareMatrices(sum_state.sum_angular_momentum,
                                Vector3<double>(14.0, -28.0, 14.0), kEps));
}

}  // namespace
}  // namespace internal
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
