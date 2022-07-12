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
// constitutive model with Hencky strain. A Von Mises plasticity model composed
// of a yield surface {τ: f(τ) <= 0}, an implicit surface describing the
// elastic/plastic limit. And an associate flow rule, an ODE relating the
// plastic deformation rate and stress, l_p = dγ/dt df(τ)/dτ, where l_p is the
// plastic rate of deformation and dγ/dt is called the plastic multiplier.
// In particular, the Von Mises yield function is (BW 7.20a,b)
//      f(τ) = √(3/2) ‖ dev(τ) ‖_F - τ_c,
//                                τ_c is the maximum allowed tensile strength
// The associate flow rule for Von Mises model can be written as (BW 7.45a,b)
//      ln(Σⁿ⁺¹) - ln(Σ) = -Δγ νⁿ⁺¹
// We first outline the general procedure of applying plasticity
// 1. Update the elastic stress Fᴱ in the G2P update, we call it the trial
//    stress
// 2. Project the elastic stress to the yield surface, proj(Fᴱ)
// 3. Update the plastic deformation gradient by the identity
//    F = proj(Fᴱ)proj(Fᴾ) = FᴱFᴾ, or proj(Fᴾ) = proj(Fᴱ)⁻¹FᴱFᴾ
// Given the SVD of elastic trial stress Fᴱ = U Σ Vᵀ,
// the Hencky strain     ε = 1/2ln(Fᴱ Fᴱ^T) = U ln(Σ) Uᵀ
// the energy density    ψ = U Σ⁻¹(μln(Σ)²+1/2*λtr(ln(Σ))²I) Vᵀ
// the Kirchhoff stress  τ = U (2μln(Σ)+λtr(ln(Σ))I) Uᵀ
// the deviatoric component of Kirchhoff stress is
//                  dev(τ) = τ - 1/3*tr(τ)I = U (2μln(Σ)+2μ/3tr(ln(Σ))I) Uᵀ
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
// same algorithm described in (BW Box 7.1), without hardening.
class StvkHenckyWithVonMisesModel: public ElastoPlasticModel {
 public:
    // tau_c: maximum allowed tensile strength
    explicit StvkHenckyWithVonMisesModel(double tau_c);
    StvkHenckyWithVonMisesModel(double E, double nu, double tau_c);

    virtual std::unique_ptr<ElastoPlasticModel> Clone() const {
        return std::make_unique<StvkHenckyWithVonMisesModel>(*this);
    }

    void CalcFirstPiolaKirchhoffStress(
            const Matrix3<double>& F, Matrix3<double>* P) const final;

    void CalcKirchhoffStress(const Matrix3<double>& F, Matrix3<double>* tau)
                                                                    const final;
    void CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
                        const Matrix3<double>& F, Matrix3<double>* P,
                        Matrix3<double>* tau) const final;

    void UpdateDeformationGradient(
                                Matrix3<double>* elastic_deformation_gradient)
                                                                    const final;

    void CalcKirchhoffStressAndUpdateDeformationGradient(
                    Matrix3<double>* tau,
                    Matrix3<double>* elastic_deformation_gradient) const final;

 private:
    double tau_c_;
};  // class StvkHenckyWithVonMisesModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
