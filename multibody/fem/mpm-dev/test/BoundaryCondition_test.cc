#include "drake/multibody/fem/mpm-dev/BoundaryCondition.h"

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace mpm {

constexpr double TOLERANCE = 1e-10;

namespace {

GTEST_TEST(BoundaryConditionTest, TestDefaultConstructor) {
    BoundaryCondition bc = BoundaryCondition();
    EXPECT_EQ(bc.get_num_boundary(), 0);
}

GTEST_TEST(BoundaryConditionTest, TestGetAdd) {
    // Construct a boundary condition class with two boundaries
    double mu0 = 0.1;
    double mu1 = 0.2;
    double mu2 = 0.3;
    Vector3<double> n0 = Vector3<double>(1.0, 0.0, 0.0);
    Vector3<double> p0 = Vector3<double>(2.0, 2.0, 2.0);
    Vector3<double> n1 = Vector3<double>(0.0, 1.0, 0.0);
    Vector3<double> p1 = Vector3<double>(3.0, 3.0, 3.0);
    Vector3<double> n2 = Vector3<double>(0.0, 0.0, 1.0);
    Vector3<double> p2 = Vector3<double>(1.0, 1.0, 1.0);
    BoundaryCondition::Boundary b0 = {mu0, {n0, p0}};
    BoundaryCondition::Boundary b1 = {mu1, {n1, p1}};
    BoundaryCondition::Boundary b2 = {mu2, {n2, p2}};
    std::vector<BoundaryCondition::Boundary> boundaries = {b0, b1};
    BoundaryCondition bc = BoundaryCondition(boundaries);

    EXPECT_EQ(bc.get_num_boundary(), 2);
    const std::vector<BoundaryCondition::Boundary>& boundaries_ret
                                                        = bc.get_boundaries();
    EXPECT_EQ(boundaries_ret[0].friction_coefficient, mu0);
    EXPECT_EQ(bc.get_boundary(0).friction_coefficient, mu0);
    EXPECT_EQ(boundaries_ret[1].friction_coefficient, mu1);
    EXPECT_EQ(bc.get_boundary(1).friction_coefficient, mu1);
    EXPECT_TRUE(CompareMatrices(boundaries_ret[0].boundary_space.normal(), n0));
    EXPECT_TRUE(CompareMatrices(bc.get_boundary(0).boundary_space.normal(),
                                                                           n0));
    EXPECT_TRUE(CompareMatrices(boundaries_ret[1].boundary_space.normal(), n1));
    EXPECT_TRUE(CompareMatrices(bc.get_boundary(1).boundary_space.normal(),
                                                                           n1));

    // Add a boundary space to BoundaryCondition
    bc.AddBoundary(b2);
    EXPECT_EQ(bc.get_num_boundary(), 3);
    EXPECT_EQ(boundaries_ret[0].friction_coefficient, mu0);
    EXPECT_EQ(bc.get_boundary(0).friction_coefficient, mu0);
    EXPECT_EQ(boundaries_ret[1].friction_coefficient, mu1);
    EXPECT_EQ(bc.get_boundary(1).friction_coefficient, mu1);
    EXPECT_EQ(boundaries_ret[2].friction_coefficient, mu2);
    EXPECT_EQ(bc.get_boundary(2).friction_coefficient, mu2);
    EXPECT_TRUE(CompareMatrices(boundaries_ret[0].boundary_space.normal(), n0));
    EXPECT_TRUE(CompareMatrices(bc.get_boundary(0).boundary_space.normal(),
                                                                           n0));
    EXPECT_TRUE(CompareMatrices(boundaries_ret[1].boundary_space.normal(), n1));
    EXPECT_TRUE(CompareMatrices(bc.get_boundary(1).boundary_space.normal(),
                                                                           n1));
    EXPECT_TRUE(CompareMatrices(boundaries_ret[2].boundary_space.normal(), n2));
    EXPECT_TRUE(CompareMatrices(bc.get_boundary(2).boundary_space.normal(),
                                                                           n2));
}

}  // namespace
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
