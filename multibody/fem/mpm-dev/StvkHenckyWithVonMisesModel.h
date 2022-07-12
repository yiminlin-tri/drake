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
// In particular, the Von Mises yield function is (BW 7.20a,b)
//      f(τ) = √(3/2) ‖ dev(τ) ‖_F - τ_c,
//                                τ_c is the maximum allowed tensile strength
// The associate flow rule for Von Mises model can be written as (BW 7.45a,b)
//      ln(Σⁿ⁺¹) - ln(Σ) = -Δγ νⁿ⁺¹
// We first outline the general procedure of applying plasticity
// 1. Update the elastic deformation gradient Fᴱ in the G2P update, compute the
//    corresponding Kirchhoff stress τ, and we call it the trial stress.
// 2. Project the trial stress to the yield surface, proj(τ)
// 3. Recover the projected deformation gradient proj(Fᴱ) by proj(τ) through the
//    constitutive relation
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
// change with respect to time.
class StvkHenckyWithVonMisesModel: public ElastoPlasticModel {
 public:
    // Yield stress is the minimum stress at which the material undergoes
    // plastic deformation
    // @pre yield_stress >= 0
    explicit StvkHenckyWithVonMisesModel(double yield_stress);
    StvkHenckyWithVonMisesModel(double E, double nu, double yield_stress);

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
    // Helper function, calculate Kirchhoff stress given the left singular
    // vectors and the singular values and the trace of the trial Hencky strain
    void CalcKirchhoffStress(const Matrix3<double>& U,
                             const Vector3<double>& eps_hat,
                             double tr_eps, Matrix3<double>* tau) const;

    // Given the left and right singular vectors of the deformation gradient,
    // the singular values and the trace of the trial Hencky strain ε: σ(ε):
    // project the elastic deformation gradient to the yield surface. Update
    // the elastic deformation gradient to be its projection
    void ProjectDeformationGradientToYieldSurface(
                                                const Matrix3<double>& U,
                                                const Matrix3<double>& V,
                                                const Vector3<double>& eps_hat,
                                                double tr_eps,
                                                Matrix3<double>* FE) const;

    double yield_stress_;
};  // class StvkHenckyWithVonMisesModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
