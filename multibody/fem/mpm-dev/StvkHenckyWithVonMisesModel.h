#pragma once

#include <memory>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/matrix_utilities.h"
#include "drake/multibody/fem/mpm-dev/ElastoPlasticModel.h"

namespace drake {
namespace multibody {
namespace mpm {

// A implementation of Saint-Venant Kirchhoff model, but replace the left
// Cauchy Green strain with the Hencky strain
// https://dl.acm.org/doi/abs/10.1145/2897824.2925906
// The formula of Kirchhoff stress can be found in Section 6.3,
// and von Mises plasticity model, described in Bonet&Wood (BW)
// https://www.klancek.si/sites/default/files/datoteke/files/
// bonet-woodnonlinearcontinuummechanics2ndedition.pdf
// The return mapping algorithm can be found in Box 7.1.
// Assume the constitutive model we are using is the Saint-Venant Kirchhoff
// constitutive model with Hencky strain. A plasticity model composed
// of a yield surface {τ: f(τ) <= 0}, an implicit surface describing the
// elastic/plastic limit. And an associate flow rule, an ODE relating the
// plastic deformation rate and stress, l_p = dγ/dt df(τ)/dτ, where l_p is the
// plastic rate of deformation and dγ/dt is called the plastic multiplier.
// We first outline the general procedure of applying plasticity
// 1. Update the elastic deformation gradient Fᴱ in the G2P update, compute the
//    corresponding Kirchhoff stress τ, and we call it the trial stress.
// 2. Project the trial stress to the yield surface, get proj(τ), if the trial
//    stress is outside the yield surface
// 3. Recover the projected deformation gradient proj(Fᴱ) by proj(τ) by
//    inverting the constitutive relation
//
// In particular, the Von Mises yield function is (BW 7.20a,b)
//      f(τ) = √(3/2) ‖ dev(τ) ‖_F - τ_c,
//                                τ_c is the maximum allowed tensile strength
// The associate flow rule for Von Mises model can be written as (BW 7.45a,b)
//      ln(Σⁿ⁺¹) - ln(Σ) = -Δγ νⁿ⁺¹
// Given the SVD of elastic deformation gradient Fᴱ = U Σ Vᵀ,
// the Hencky strain             ε = 1/2ln(Fᴱ Fᴱ^T) = U ln(Σ) Uᵀ
// the energy density            ψ = U Σ⁻¹(μln(Σ)²+1/2*λtr(ln(Σ))²I) Vᵀ
// the (trial) Kirchhoff stress  τ = U (2μln(Σ)+λtr(ln(Σ))I) Uᵀ
// the deviatoric component of trial stress is
//                  dev(τ) = τ - 1/3*tr(τ)I = U (2μln(Σ)+2/3*μ tr(ln(Σ))I) Uᵀ
// The we can rewrite the associate flow rule with (BW 7.53)
//       dev(τⁿ⁺¹) - dev(τ) = -2μ Δγ νⁿ⁺¹
// τⁿ⁺¹ denotes the projection of the trial stress: proj(τ)
// One key observation is that the flow direction νⁿ⁺¹ is has the same direction
// as the deviatoric components of updated Kirchhoff stress (BW 7.40)
//       νⁿ⁺¹ = df(τⁿ⁺¹)/dτ = dev(τⁿ⁺¹)/(√(3/2) ‖ dev(τⁿ⁺¹) ‖_F)
// As a result, the deviatoric components of updated and trial Kirchhoff stress
// are in the same direction, i.e.
//      dev(τⁿ⁺¹) = k dev(τ)
//                 k: f(τⁿ⁺¹) = 0 on the yield surface
// Solving for k gives us the prjected states. The above derivation leads to the
// same algorithm described in (BW Box 7.1), without hardening. We assume
// no hardening in the plastic model, so the yield surface (function) doesn't
// change with respect to time. This is equivalent to the formulation in
// (BW Box 7.1) with H = 0.
class StvkHenckyWithVonMisesModel: public ElastoPlasticModel {
 public:
    // Yield stress is the minimum stress at which the material undergoes
    // plastic deformation
    // @pre yield_stress > 0
    explicit StvkHenckyWithVonMisesModel(double yield_stress);
    StvkHenckyWithVonMisesModel(double E, double nu, double yield_stress);

    virtual std::unique_ptr<ElastoPlasticModel> Clone() const {
        return std::make_unique<StvkHenckyWithVonMisesModel>(*this);
    }

    void UpdateDeformationGradientAndCalcKirchhoffStress(
            Matrix3<double>* tau,
            Matrix3<double>* trial_elastic_deformation_gradient) const final;

    // Evaluate the yield function f in the plasticity model given the elastic
    // deformation gradient FE. (f <= 0 if the stress is in/on the yield surface
    // , and there is only elastic response. If f > 0, there will be plastic
    // response). Return 0.0 for a purely elastic model.
    double EvalYieldFunction(const Matrix3<double>& FE) const;

 private:
    // Given the Hencky strain ε, assume he SVD of it is U Σ Vᵀ,
    // We store U in `U`, V in `V`, eps_hat in `Σ`, and tr(Σ) in `tr_eps`
    struct StrainData {
        Matrix3<double> U;
        Matrix3<double> V;
        Vector3<double> eps_hat;
        double tr_eps;
    };

    // Given the Kirchhoff stress τ, the deviatoric component of τ is defined as
    // dev(τ) = τ - pJI = τ - 1/3tr(τ)I
    // Store dev(τ)'s singular values into `tau_dev_hat` and the norm of
    // `tau_dev_hat` in `tau_dev_hat_norm`
    struct StressData {
        Vector3<double> tau_dev_hat;
        double tau_dev_hat_norm;
    };

    // Calculate the strain data using the elastic deformation gradient `FE`
    StrainData CalcStrainData(const Matrix3<double>& FE) const;

    // Calculate the stress data using the strain data  `strain_data`
    StressData CalcStressData(const StrainData& strain_data) const;

    // Helper function, calculate Kirchhoff stress given the left singular
    // vectors and the singular values and the trace of the trial Hencky strain
    void CalcKirchhoffStress(const StrainData& trial_strain_data,
                             Matrix3<double>* tau) const;

    // Helper function, evaluate the yield function
    // f(τ) = sqrt(3/2)‖ dev(τ) ‖ - τ_c, where dev(τ) is the deviatoric
    // component of Kirchhoff stress in the princpal frame
    double EvalYieldFunction(const StressData& trialStressData) const;

    // Given the left and right singular vectors of the deformation gradient (U
    // and V), the singular values and the trace of the trial Hencky strain
    // (eps_hat and tr_eps), (stored in trial_strain_data) if the corresponding
    // stress is outside the yield surface, project it back to the yield surface
    // using associative flow rule. Then, invert the constitutive relation and
    // write the elastic deformation gradient corresponding to the projected
    // stress to `FE`, and write the projected Hencky strain in the principal
    // frame to `trial_strain_data->eps_hat`.
    void ProjectDeformationGradientToYieldSurface(
                                        const StressData& trial_stress_data,
                                        StrainData* trial_strain_data,
                                        Matrix3<double>* FE_trial) const;

    // sqrt(3.0/2.0)
    static constexpr double sqrt_32 = 1.224744871391589;
    double yield_stress_;
};  // class StvkHenckyWithVonMisesModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
