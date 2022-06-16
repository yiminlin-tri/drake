#include "drake/multibody/fem/mpm-dev/BoundaryCondition.h"

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"
#include "drake/math/roll_pitch_yaw.h"
#include "drake/multibody/fem/mpm-dev/AnalyticLevelSet.h"

namespace drake {
namespace multibody {
namespace mpm {

constexpr double TOLERANCE = 1e-10;

namespace {

GTEST_TEST(BoundaryConditionTest, TestDefaultConstructor) {
    BoundaryCondition bc = BoundaryCondition();
    EXPECT_EQ(bc.get_num_wall_boundaries(), 0);
    EXPECT_EQ(bc.get_num_moving_cylindrical_boundaries(), 0);
}

GTEST_TEST(BoundaryConditionTest, TestGetAdd) {
    // Construct a boundary condition class with two wall boundaries and a
    // moving cylinder boundary
    double mu0 = 0.1;
    double mu1 = 0.2;
    double mu2 = 0.3;
    Vector3<double> n0 = Vector3<double>(1.0, 0.0, 0.0);
    Vector3<double> p0 = Vector3<double>(2.0, 2.0, 2.0);
    Vector3<double> n1 = Vector3<double>(0.0, 1.0, 0.0);
    Vector3<double> p1 = Vector3<double>(3.0, 3.0, 3.0);
    Vector3<double> n2 = Vector3<double>(0.0, 0.0, 1.0);
    Vector3<double> p2 = Vector3<double>(1.0, 1.0, 1.0);
    BoundaryCondition::WallBoundary b0 = {mu0, {n0, p0}};
    BoundaryCondition::WallBoundary b1 = {mu1, {n1, p1}};
    BoundaryCondition::WallBoundary b2 = {mu2, {n2, p2}};

    double f_coeff0 = 0.1;
    double f_coeff1 = 0.2;
    CylinderLevelSet level_set_cyl0 = CylinderLevelSet(0.1, 0.2);
    CylinderLevelSet level_set_cyl1 = CylinderLevelSet(0.3, 0.4);
    Vector3<double> translation_cyl0 = {0.1, 0.2, 0.3};
    Vector3<double> translation_cyl1 = {-0.1, -0.2, -0.3};
    math::RigidTransform<double> pose_cyl0 =
                        math::RigidTransform<double>(translation_cyl0);
    math::RigidTransform<double> pose_cyl1 =
                        math::RigidTransform<double>(translation_cyl1);
    Vector3<double> v_cyl0 = {1.0, 2.0, 3.0};
    Vector3<double> v_cyl1 = {-1.0, -2.0, -3.0};
    BoundaryCondition::MovingCylindricalBoundary c0 =
                                {f_coeff0, level_set_cyl0, pose_cyl0, v_cyl0};
    BoundaryCondition::MovingCylindricalBoundary c1 =
                                {f_coeff1, level_set_cyl1, pose_cyl1, v_cyl1};
    std::vector<BoundaryCondition::WallBoundary> wall_boundaries = {b0, b1};
    std::vector<BoundaryCondition::MovingCylindricalBoundary> cyl_boundaries
                                                                        = {c0};
    BoundaryCondition bc = BoundaryCondition(wall_boundaries, cyl_boundaries);

    EXPECT_EQ(bc.get_num_wall_boundaries(), 2);
    const std::vector<BoundaryCondition::WallBoundary>& wall_boundaries_ret
                                                    = bc.get_wall_boundaries();
    const std::vector<BoundaryCondition::MovingCylindricalBoundary>&
                                                        cyl_boundaries_ret
                                    = bc.get_moving_cylindrical_boundaries();
    EXPECT_EQ(wall_boundaries_ret[0].friction_coefficient, mu0);
    EXPECT_EQ(bc.get_wall_boundary(0).friction_coefficient, mu0);
    EXPECT_EQ(wall_boundaries_ret[1].friction_coefficient, mu1);
    EXPECT_EQ(bc.get_wall_boundary(1).friction_coefficient, mu1);
    EXPECT_TRUE(CompareMatrices(wall_boundaries_ret[0].boundary_space.normal(),
                                                                           n0));
    EXPECT_TRUE(CompareMatrices(bc.get_wall_boundary(0).boundary_space.normal(),
                                                                           n0));
    EXPECT_TRUE(CompareMatrices(wall_boundaries_ret[1].boundary_space.normal(),
                                                                           n1));
    EXPECT_TRUE(CompareMatrices(bc.get_wall_boundary(1).boundary_space.normal(),
                                                                           n1));

    EXPECT_EQ(bc.get_num_moving_cylindrical_boundaries(), 1);
    EXPECT_EQ(cyl_boundaries_ret[0].friction_coefficient, f_coeff0);
    EXPECT_EQ(bc.get_moving_cylindrical_boundary(0).friction_coefficient,
                                                          f_coeff0);
    EXPECT_TRUE(CompareMatrices(cyl_boundaries_ret[0].pose.GetAsMatrix34(),
                                                 pose_cyl0.GetAsMatrix34()));
    EXPECT_TRUE(CompareMatrices(
                bc.get_moving_cylindrical_boundary(0).pose.GetAsMatrix34(),
                                                 pose_cyl0.GetAsMatrix34()));
    EXPECT_TRUE(CompareMatrices(cyl_boundaries_ret[0].velocity, v_cyl0));
    EXPECT_TRUE(CompareMatrices(bc.get_moving_cylindrical_boundary(0).velocity,
                                                                v_cyl0));

    // Add a boundary space and a moving cylinder to BoundaryCondition
    bc.AddWallBoundary(b2);
    bc.AddMovingCylindricalBoundary(c1);
    EXPECT_EQ(bc.get_num_wall_boundaries(), 3);
    EXPECT_EQ(wall_boundaries_ret[0].friction_coefficient, mu0);
    EXPECT_EQ(bc.get_wall_boundary(0).friction_coefficient, mu0);
    EXPECT_EQ(wall_boundaries_ret[1].friction_coefficient, mu1);
    EXPECT_EQ(bc.get_wall_boundary(1).friction_coefficient, mu1);
    EXPECT_EQ(wall_boundaries_ret[2].friction_coefficient, mu2);
    EXPECT_EQ(bc.get_wall_boundary(2).friction_coefficient, mu2);
    EXPECT_TRUE(CompareMatrices(wall_boundaries_ret[0].boundary_space.normal(),
                                                                           n0));
    EXPECT_TRUE(CompareMatrices(bc.get_wall_boundary(0).boundary_space.normal(),
                                                                           n0));
    EXPECT_TRUE(CompareMatrices(wall_boundaries_ret[1].boundary_space.normal(),
                                                                           n1));
    EXPECT_TRUE(CompareMatrices(bc.get_wall_boundary(1).boundary_space.normal(),
                                                                           n1));
    EXPECT_TRUE(CompareMatrices(wall_boundaries_ret[2].boundary_space.normal(),
                                                                           n2));
    EXPECT_TRUE(CompareMatrices(bc.get_wall_boundary(2).boundary_space.normal(),
                                                                           n2));

    EXPECT_EQ(bc.get_num_moving_cylindrical_boundaries(), 2);
    EXPECT_EQ(cyl_boundaries_ret[1].friction_coefficient, f_coeff1);
    EXPECT_EQ(bc.get_moving_cylindrical_boundary(1).friction_coefficient,
                                                          f_coeff1);
    EXPECT_TRUE(CompareMatrices(cyl_boundaries_ret[1].pose.GetAsMatrix34(),
                                                 pose_cyl1.GetAsMatrix34()));
    EXPECT_TRUE(CompareMatrices(
                bc.get_moving_cylindrical_boundary(1).pose.GetAsMatrix34(),
                                                 pose_cyl1.GetAsMatrix34()));
    EXPECT_TRUE(CompareMatrices(cyl_boundaries_ret[1].velocity, v_cyl1));
    EXPECT_TRUE(CompareMatrices(bc.get_moving_cylindrical_boundary(1).velocity,
                                                                v_cyl1));
}

GTEST_TEST(BoundaryConditionTest, TestApplyWall) {
    // Construct a boundary
    double mu0 = 0.1;
    Vector3<double> n0 = Vector3<double>(1.0, 1.0, 1.0);
    Vector3<double> p0 = Vector3<double>(2.0, 2.0, 2.0);
    BoundaryCondition::WallBoundary b0 = {mu0, {n0, p0}};
    std::vector<BoundaryCondition::WallBoundary> wall_boundaries = {b0};
    std::vector<BoundaryCondition::MovingCylindricalBoundary> cyl_boundaries =
                                                                            {};
    BoundaryCondition bc = BoundaryCondition(wall_boundaries, cyl_boundaries);
    double dummy_t = 0.0;

    // First consider a point inside the domain, no boundary condition shall be
    // enforced.
    Vector3<double> pos0 = {3.0, 3.0, 3.0};
    Vector3<double> vel0 = {1.0, -1.0, 0.8};
    bc.Apply(dummy_t, pos0, &vel0);
    ASSERT_TRUE(CompareMatrices(vel0, Vector3<double>{1.0, -1.0, 0.8},
                                                                    TOLERANCE));

    // Next consider a point on the boundary half plane, boundary condition
    // shall be enforced, the velocity is in the normal direction, so the new
    // velocity should be zero.
    Vector3<double> pos1 = {2.0, 2.0, 2.0};
    Vector3<double> vel1 = {10.0, 10.0, 10.0};
    bc.Apply(dummy_t, pos1, &vel1);
    ASSERT_TRUE(CompareMatrices(vel1, Vector3<double>::Zero(), TOLERANCE));

    // Next consider a point on the boundary half plane. The velocity is in not
    // the normal direction, but the normal impulse and the consequent friction
    // is so large that only the tangential velocity is left after hitting the
    // boundary.
    Vector3<double> pos2 = {2.0, 2.0, 2.0};
    Vector3<double> vel2 = {11.0, 10.0, 10.0};
    bc.Apply(dummy_t, pos2, &vel2);
    ASSERT_TRUE(CompareMatrices(vel2, Vector3<double>::Zero(), TOLERANCE));

    // Finally test the frictional wall boundary condition enforcement using
    // a different velocity such that the friction will not be too large
    Vector3<double> pos3 = {2.0, 2.0, 2.0};
    Vector3<double> vel3 = {2.0, 1.0, 1.0};
    bc.Apply(dummy_t, pos3, &vel3);
    ASSERT_TRUE(CompareMatrices(vel3, (1-mu0*4.0/sqrt(2))/3.0
                                *Vector3<double>{2.0, -1.0, -1.0}, TOLERANCE));
}

GTEST_TEST(BoundaryConditionTest, TestApplyMovingCylindricalBC) {
    // Construct a moving cylindrical boundary with central axis parallel to
    // y-axis, centered at (1.0, 2.0, 3.0), and moving with velocity (1, 0, 0)
    math::RollPitchYaw rpw_cylinder = {M_PI/2.0, 0.0, 0.0};
    CylinderLevelSet level_set_cylinder = CylinderLevelSet(2.0, 1.0);
    Vector3<double> translation_cylinder = {1.0, 2.0, 3.0};
    math::RigidTransform<double> pose_cylinder = {rpw_cylinder,
                                                  translation_cylinder};
    double mu = 0.1;
    BoundaryCondition::MovingCylindricalBoundary mc0 = {mu,
                                                        level_set_cylinder,
                                                        pose_cylinder,
                                                        {1.0, 0.0, 0.0}};
    std::vector<BoundaryCondition::WallBoundary> wall_boundaries = {};
    std::vector<BoundaryCondition::MovingCylindricalBoundary> cyl_boundaries =
                                                                        {mc0};
    BoundaryCondition bc = BoundaryCondition(wall_boundaries, cyl_boundaries);
    double t0 = 0.0;

    // First consider a particle with zero initial velocity not in the cylinder,
    // no boundary condition should be enforced.
    Vector3<double> pos0 = {3.0, 2.5, 2.5};
    Vector3<double> vel0 = Vector3<double>::Zero();
    bc.Apply(t0, pos0, &vel0);
    ASSERT_TRUE(CompareMatrices(vel0, Vector3<double>::Zero(), TOLERANCE));

    // After the cylinder move, the point will lie in the cylinder
    // In this case, the relative velocity of the grid point is (-1, 0, 0),
    // vₙ = (v ⋅ n)n = (-1/2, 0, -1/2), vₜ = v - vₙ = (1/2, 0, -1/2)
    // v_new = vₜ - μ‖vₙ‖t = (1/2, 0, -1/2) - 0.1*(1, 0, -1)/2
    // In physical frame, v_new = (1, 0, 0) + (-0.45, 0, -0.45)
    t0 += 1.5;
    bc.Apply(t0, pos0, &vel0);
    ASSERT_TRUE(CompareMatrices(vel0, Vector3<double>{0.55, 0.0, -0.45},
                                TOLERANCE));
}

}  // namespace
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
