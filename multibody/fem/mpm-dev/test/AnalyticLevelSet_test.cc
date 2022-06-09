#include "drake/multibody/fem/mpm-dev/AnalyticLevelSet.h"

#include <math.h>

#include <gtest/gtest.h>

#include "drake/math/roll_pitch_yaw.h"
#include "drake/math/rotation_matrix.h"

namespace drake {
namespace multibody {
namespace mpm {
namespace internal {
namespace {

constexpr double TOLERANCE = 1e-12;

GTEST_TEST(AnalyticLevelSetTest, BoxTest) {
    const Vector3<double> center = {1.0, 2.0, 3.0};
    const Vector3<double> xscale = {3.0, 2.0, 1.0};
    const math::RotationMatrix<double> rotation0 =
        math::RotationMatrix<double>(math::RollPitchYaw<double>(0.0, 0.0, 0.0));
    const math::RotationMatrix<double> rotation1 =
        math::RotationMatrix<double>
            (math::RollPitchYaw<double>(M_PI/4.0, M_PI/4.0, M_PI/4.0));

    // A box without rotation
    BoxLevelSet box = BoxLevelSet(center, xscale, rotation0);
    EXPECT_EQ(box.get_volume(), 48.0);
    EXPECT_TRUE(box.InInterior(center));
    EXPECT_TRUE(box.InInterior({-1.0, 2.0, 3.0}));
    EXPECT_TRUE(box.InInterior({2.1, 3.9, 3.9}));
    EXPECT_FALSE(box.InInterior({0.0, 5.0, 0.0}));
    EXPECT_FALSE(box.InInterior({-2.0, 1.0, 3.0}));

    // A box with rotation
    box = BoxLevelSet(center, xscale, rotation1);
    EXPECT_EQ(box.get_volume(), 48.0);
    EXPECT_TRUE(box.InInterior(center));
    EXPECT_TRUE(box.InInterior({2.0, 2.0, 2.0}));
    EXPECT_TRUE(box.InInterior({1.0, 4.0, 3.0}));
    EXPECT_FALSE(box.InInterior({0.0, 0.0, 0.0}));
}

GTEST_TEST(AnalyticLevelSetTest, SphereTest) {
    const Vector3<double> center = {0.5, -0.6, 0.7};
    double radius = 1.0;

    SphereLevelSet sphere = SphereLevelSet(center, radius);
    EXPECT_NEAR(sphere.get_volume(), 4.0/3.0*M_PI, TOLERANCE);
    EXPECT_TRUE(sphere.InInterior(center));
    EXPECT_TRUE(sphere.InInterior({0.2, 0.0, 0.2}));
    EXPECT_FALSE(sphere.InInterior({-0.5, -0.6, 0.7}));
    EXPECT_FALSE(sphere.InInterior({0.0, 0.0, 0.0}));
}

}  // namespace
}  // namespace internal
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
