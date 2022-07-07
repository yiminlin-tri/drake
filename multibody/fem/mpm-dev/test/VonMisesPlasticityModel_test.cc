#include "drake/multibody/fem/mpm-dev/VonMisesPlasticityModel.h"

#include <memory>
#include <utility>

#include <gtest/gtest.h>

#include "drake/common/test_utilities/eigen_matrix_compare.h"
#include "drake/multibody/fem/mpm-dev/SaintVenantKirchhoffWithHenckyModel.h"

namespace drake {
namespace multibody {
namespace mpm {
namespace internal {
namespace {

constexpr double TOLERANCE = 1e-10;

GTEST_TEST(CorotatedModelTest, CorotatedModelCalculationTest) {
    Matrix3<double> FE = (Matrix3<double>() <<
            6, 1, 2,
            1, 4, 1,
            2, 1, 5).finished();
    Matrix3<double> FP = (Matrix3<double>() <<
            1, 0, 0,
            0, 1, 0,
            0, 0, 1).finished();

    double E      = 100.0;
    double nu     = 0.2;
    std::unique_ptr<SaintVenantKirchhoffWithHenckyModel> hencky_model
                = std::make_unique<SaintVenantKirchhoffWithHenckyModel>(E, nu);
    double mu     = hencky_model->get_mu();
    // We set a tau_c large enough so that the current elastic deformation
    // gradient is in the yield surface
    double tau_c0  = 1000.0;
    // Calculate the Kirchhoff Stress
    Matrix3<double> tau0;
    hencky_model->CalcKirchhoffStress(FE, &tau0);
    // Calculate the deviatoric component of the Kirchhoff stress
    // dev(τ) = τ - 1/3tr(τ)I
    Matrix3<double> tau_dev0 = tau0
                             - 1.0/3.0*tau0.trace()*Matrix3<double>::Identity();
    // The trial stress is in the yield surface if
    // sqrt(3/2)‖ dev(τ) ‖_F ≤ τ_c
    bool in_yield_surface0 = (std::sqrt(3.0/2.0)*tau_dev0.norm()
                           <= tau_c0 + TOLERANCE);
    EXPECT_TRUE(in_yield_surface0);

    VonMisesPlasticityModel plasticity_model0
                                = VonMisesPlasticityModel(tau_c0);
    // Deformation gradient before plasticity
    Matrix3<double> FEprev = FE;
    Matrix3<double> FPprev = FP;
    Matrix3<double> Fprev  = FE*FP;
    plasticity_model0.UpdateDeformationGradients(mu, &FE, &FP);

    // Calculate the new Kirchhoff Stress
    hencky_model->CalcKirchhoffStress(FE, &tau0);
    tau_dev0 = tau0 - 1.0/3.0*tau0.trace()*Matrix3<double>::Identity();
    in_yield_surface0 = (std::sqrt(3.0/2.0)*tau_dev0.norm()
                      <= tau_c0 + TOLERANCE);

    // Nothing shall change in this case
    EXPECT_TRUE(in_yield_surface0);
    EXPECT_TRUE(CompareMatrices(FEprev, FE, TOLERANCE));
    EXPECT_TRUE(CompareMatrices(FPprev, FP, TOLERANCE));
    EXPECT_TRUE(CompareMatrices(Fprev , FE*FP, TOLERANCE));
    EXPECT_NEAR(FP.determinant(), 1.0, TOLERANCE);

    // We set a tau_c small enough so that the current elastic deformation
    // gradient is not in the yield surface
    double tau_c1  = 70.0;
    // Calculate the Kirchhoff Stress
    Matrix3<double> tau1;
    hencky_model->CalcKirchhoffStress(FE, &tau1);
    Matrix3<double> tau_dev1 = tau1
                             - 1.0/3.0*tau1.trace()*Matrix3<double>::Identity();
    bool in_yield_surface1 = (std::sqrt(3.0/2.0)*tau_dev1.norm()
                          <= tau_c1 + TOLERANCE);
    EXPECT_FALSE(in_yield_surface1);
    VonMisesPlasticityModel plasticity_model1
                                = VonMisesPlasticityModel(tau_c1);

    // Deformation gradient before plasticity
    FEprev = FE;
    FPprev = FP;
    Fprev  = FE*FP;
    plasticity_model1.UpdateDeformationGradients(mu, &FE, &FP);

    // Calculate the new Kirchhoff Stress, it should be in the yield surface
    hencky_model->CalcKirchhoffStress(FE, &tau1);
    tau_dev1 = tau1 - 1.0/3.0*tau1.trace()*Matrix3<double>::Identity();
    in_yield_surface1 = (std::sqrt(3.0/2.0)*tau_dev1.norm()
                      <= tau_c1 + TOLERANCE);

    EXPECT_TRUE(in_yield_surface1);
    EXPECT_FALSE(CompareMatrices(FEprev, FE, TOLERANCE));
    EXPECT_FALSE(CompareMatrices(FPprev, FP, TOLERANCE));
    EXPECT_TRUE(CompareMatrices(Fprev , FE*FP, TOLERANCE));
    EXPECT_NEAR(FP.determinant(), 1.0, TOLERANCE);
}

}  // namespace
}  // namespace internal
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
