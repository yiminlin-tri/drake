#include "drake/multibody/fem/mpm-dev/Particles.h"

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace mpm {
namespace internal {
namespace {

constexpr double kEps = 4.0 * std::numeric_limits<double>::epsilon();

GTEST_TEST(ParticlesClassTest, RoundTrip) {
    Vector3<double> pos1 = {1.0, 2.0, 3.0};
    Vector3<double> vel1 = {-1.0, -2.0, -3.0};
    double mass1 = 5.0;
    double vol1  = 10.0;
    Matrix3<double> dg1 = pos1.asDiagonal();
    Matrix3<double> stress1 = vel1.asDiagonal();

    Vector3<double> pos2 = {3.0, -1.0, 6.0};
    Vector3<double> vel2 = {-9.0, 8.0, -2.0};
    double mass2 = 7.0;
    double vol2  = 3.0;
    Matrix3<double> dg2 = pos2.asDiagonal();
    Matrix3<double> stress2 = vel2.asDiagonal();

    Particles particles = Particles();
    particles.add_particle(pos1, vel1, mass1, vol1, dg1, stress1);
    particles.add_particle(pos2, vel2, mass2, vol2, dg2, stress2);

    EXPECT_EQ(particles.get_num_particles(),2);
    EXPECT_TRUE(CompareMatrices(particles.get_position_at(0), pos1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity_at(0), vel1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass_at(0), mass1);
    EXPECT_EQ(particles.get_volume_at(0), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient_at(0), dg1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_stress_at(0), stress1,
                std::numeric_limits<double>::epsilon()));

    EXPECT_TRUE(CompareMatrices(particles.get_position_at(1), pos2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity_at(1), vel2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass_at(1), mass2);
    EXPECT_EQ(particles.get_volume_at(1), vol2);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient_at(1), dg2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_stress_at(1), stress2,
                std::numeric_limits<double>::epsilon()));

    particles = Particles(2);
    particles.set_position_at(0, pos1);
    particles.set_velocity_at(0, vel1);
    particles.set_mass_at(0, mass1);
    particles.set_volume_at(0, vol1);
    particles.set_deformation_gradient_at(0, dg1);
    particles.set_stress_at(0, stress1);
    particles.set_position_at(1, pos2);
    particles.set_velocity_at(1, vel2);
    particles.set_mass_at(1, mass2);
    particles.set_volume_at(1, vol2);
    particles.set_deformation_gradient_at(1, dg2);
    particles.set_stress_at(1, stress2);

    /*
    EXPECT_TRUE(CompareMatrices(particles.get_position_at(0), pos1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity_at(0), vel1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass_at(0), mass1);
    EXPECT_EQ(particles.get_volume_at(0), vol1);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient_at(0), dg1,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_stress_at(0), stress1,
                std::numeric_limits<double>::epsilon()));

    EXPECT_TRUE(CompareMatrices(particles.get_position_at(1), pos2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_velocity_at(1), vel2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_EQ(particles.get_mass_at(1), mass2);
    EXPECT_EQ(particles.get_volume_at(1), vol2);
    EXPECT_TRUE(CompareMatrices(particles.get_deformation_gradient_at(1), dg2,
                std::numeric_limits<double>::epsilon()));
    EXPECT_TRUE(CompareMatrices(particles.get_stress_at(1), stress2,
                std::numeric_limits<double>::epsilon()));
    */
}

}  // namespace
}  // namespace internal
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
