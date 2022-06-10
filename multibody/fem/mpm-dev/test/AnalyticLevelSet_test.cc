#include "drake/multibody/fem/mpm-dev/AnalyticLevelSet.h"

#include <math.h>

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace mpm {
namespace internal {
namespace {

constexpr double TOLERANCE = 1e-12;

GTEST_TEST(AnalyticLevelSetTest, BoxTest) {
    const Vector3<double> xscale = {3.0, 2.0, 1.0};

    // Check InInterior
    BoxLevelSet box = BoxLevelSet(xscale);
    EXPECT_EQ(box.get_volume(), 48.0);
    EXPECT_TRUE(box.InInterior({0.0, 0.0, 0.0}));
    EXPECT_TRUE(box.InInterior({1.0, -1.0, 0.0}));
    EXPECT_TRUE(box.InInterior({-2.5, -1.5, 0.0}));
    EXPECT_FALSE(box.InInterior({1.0, 2.0, 3.0}));
    EXPECT_FALSE(box.InInterior({-3.0, -2.0, -1.0}));

    // Check Normal
    EXPECT_TRUE(CompareMatrices(box.Normal({0.0, 0.0, 0.1}),
                                           Vector3<double>{0.0, 0.0, 1.0},
                                           TOLERANCE));
    EXPECT_TRUE(CompareMatrices(box.Normal({1.0, 1.0, -0.9}),
                                           Vector3<double>{0.0, 0.0, -1.0},
                                           TOLERANCE));
    EXPECT_TRUE(CompareMatrices(box.Normal({1.0, -1.8, 0.0}),
                                           Vector3<double>{0.0, -1.0, 0.0},
                                           TOLERANCE));
    EXPECT_TRUE(CompareMatrices(box.Normal({1.0, 1.7, 0.7}),
                                           Vector3<double>{0.0, 1.0, 0.0},
                                           TOLERANCE));
    EXPECT_TRUE(CompareMatrices(box.Normal({2.5, -1.0, 0.5}),
                                           Vector3<double>{1.0, 0.0, 0.0},
                                           TOLERANCE));
    EXPECT_TRUE(CompareMatrices(box.Normal({-2.0, -1.0, 0.5}),
                                           Vector3<double>{-1.0, 0.0, 0.0},
                                           TOLERANCE));
    EXPECT_TRUE(CompareMatrices(box.Normal({4.0, 3.0, 3.0}),
                                           Vector3<double>{0.0, 0.0, 0.0},
                                           TOLERANCE));

    // Check Bounding Box
    const std::array<Vector3<double>, 2>& bounding_box = box.get_bounding_box();
    EXPECT_TRUE(CompareMatrices(bounding_box[0],
                                Vector3<double>{-3.0, -2.0, -1.0},
                                TOLERANCE));
    EXPECT_TRUE(CompareMatrices(bounding_box[1],
                                Vector3<double>{3.0, 2.0, 1.0},
                                TOLERANCE));
}

GTEST_TEST(AnalyticLevelSetTest, SphereTest) {
    double radius = 2.0;

    SphereLevelSet sphere = SphereLevelSet(radius);

    // Check InInterior
    EXPECT_NEAR(sphere.get_volume(), 32.0/3.0*M_PI, TOLERANCE);
    EXPECT_TRUE(sphere.InInterior({0.0, 0.0, 0.0}));
    EXPECT_TRUE(sphere.InInterior({0.2, 0.0, 0.2}));
    EXPECT_TRUE(sphere.InInterior({-0.5, -0.6, 0.7}));
    EXPECT_FALSE(sphere.InInterior({3.0, 0.0, -1.0}));

    // Check normal
    EXPECT_TRUE(CompareMatrices(sphere.Normal({1.0, 0.0, 0.0}),
                                              Vector3<double>{1.0, 0.0, 0.0},
                                              TOLERANCE));
    EXPECT_TRUE(CompareMatrices(sphere.Normal({0.0, -1.0, 0.0}),
                                              Vector3<double>{0.0, -1.0, 0.0},
                                              TOLERANCE));
    EXPECT_TRUE(CompareMatrices(sphere.Normal({0.0, 0.0, 1.0}),
                                              Vector3<double>{0.0, 0.0, 1.0},
                                              TOLERANCE));

    // Check Bounding Box
    const std::array<Vector3<double>, 2>& bounding_box =
                                                    sphere.get_bounding_box();
    EXPECT_TRUE(CompareMatrices(bounding_box[0],
                                Vector3<double>{-2.0, -2.0, -2.0},
                                TOLERANCE));
    EXPECT_TRUE(CompareMatrices(bounding_box[1],
                                Vector3<double>{2.0, 2.0, 2.0},
                                TOLERANCE));
}

GTEST_TEST(AnalyticLevelSetTest, CylinderTest) {
    double radius = 2.0;
    double height = 0.5;

    CylinderLevelSet cylinder = CylinderLevelSet(height, radius);

    // Check InInterior
    EXPECT_NEAR(cylinder.get_volume(), 4.0*M_PI, TOLERANCE);
    EXPECT_TRUE(cylinder.InInterior({0.0, 0.0, 0.0}));
    EXPECT_TRUE(cylinder.InInterior({0.2, 0.0, 0.2}));
    EXPECT_FALSE(cylinder.InInterior({-0.5, -0.6, 0.7}));
    EXPECT_FALSE(cylinder.InInterior({3.0, 0.0, -1.0}));

    // Check normal
    EXPECT_TRUE(CompareMatrices(cylinder.Normal({1.0, 0.0, 0.0}),
                                                Vector3<double>{1.0, 0.0, 0.0},
                                                TOLERANCE));
    EXPECT_TRUE(CompareMatrices(cylinder.Normal({0.0, -1.0, 0.0}),
                                                Vector3<double>{0.0, -1.0, 0.0},
                                                TOLERANCE));
    EXPECT_TRUE(CompareMatrices(cylinder.Normal({1.0, -1.0, 0.2}),
                                    .5*Vector3<double>{sqrt(2), -sqrt(2), 0.0},
                                    TOLERANCE));

    // Check Bounding Box
    const std::array<Vector3<double>, 2>& bounding_box =
                                                    cylinder.get_bounding_box();
    EXPECT_TRUE(CompareMatrices(bounding_box[0],
                                Vector3<double>{-2.0, -2.0, -0.5},
                                TOLERANCE));
    EXPECT_TRUE(CompareMatrices(bounding_box[1],
                                Vector3<double>{2.0, 2.0, 0.5},
                                TOLERANCE));
}

}  // namespace
}  // namespace internal
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
