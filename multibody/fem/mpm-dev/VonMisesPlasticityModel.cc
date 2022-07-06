#include "drake/multibody/fem/mpm-dev/VonMisesPlasticityModel.h"

namespace drake {
namespace multibody {
namespace mpm {

VonMisesPlasticityModel::VonMisesPlasticityModel(double tau_c): tau_c_(tau_c) {
    DRAKE_ASSERT(tau_c >= 0);
}

void VonMisesPlasticityModel::UpdateDeformationGradients(double mu,
                                                         double lambda,
                                                         Matrix3<double>* FE,
                                                         Matrix3<double>* FP) {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(*FE, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U = svd.matrixU();
    const Matrix3<double>& V = svd.matrixV();
    Vector3<double> sigma    = svd.singularValues();
    // Singular values of the trial strain ε: σ(ε)
    Vector3<double> eps_hat = sigma.array().log();
    // trace of the trial Hencky strain tr(ε)
    double tr_eps           = eps_hat.sum();
    Vector3<double> tr_eps_vec = tr_eps*Vector3<double>::Ones();
    // Singular values of the trial stress τ: σ(τ)
    // Σᵢ = 2μ log(σᵢ) + λ ∑ᵢ log(σᵢ)
    Vector3<double> tau_hat = lambda*tr_eps_vec + 2*mu*eps_hat;
    // trace of the trial stress tr(τ)
    double tr_tau = tau_hat.sum();
    // Vector with all entries equals to tr(τ), or tr(τ)𝟙
    Vector3<double> tr_tau_vec = tr_tau*Vector3<double>::Ones();

    // If the trial stress τ is in the yield surface f(τ) <= 0, plasticity is
    // not applied.
    // Otherwise, project the trial stress τ in the plastic flow direction
    // df/dτ τ onto the yield surface f(τ) <= 0
    // The trial stress is on the yield surface, f(τ) <= 0, if and only if the
    // singular values of the trial stress satisfies:
    // ‖ σ(τ) - tr(τ)𝟙 ‖ ≤ τ_c
    bool in_yield_surface = (tau_hat - tr_tau_vec).norm() <= tau_c_;
    if (!in_yield_surface) {
        // C is the matrix of singular values of fourth order elasticity tensor
        // C = 2μ𝕀 + λ I⊗I is the fourth order elasticity tensor
        Matrix3<double> C;
        double diag_entry = 2.0*mu+lambda;
        C << diag_entry, lambda, lambda,
             lambda, diag_entry, lambda,
             lambda, lambda, diag_entry;
        // The singular values of trial stress' projection onto yield surface
        // is σ(τ_proj) = p + τ_c d/‖d‖
        // The singular values of trial strain's projection onto yield surface
        // is σ(ε_proj) = C⁻¹σ(τ_proj) = C⁻¹(p + τ_c d/‖d‖)
        // By the definition of Hencky strain,
        // σ(F_proj) = exp(σ(ε_proj))
        Vector3<double> p = 1.0/3.0*tr_tau_vec;
        Vector3<double> d = tau_hat - p;
        Vector3<double> proj_F_hat = (C.householderQr()
                                       .solve(p+tau_c_*d/d.norm()))
                                    .array().exp();
        // Total deformation gradient in the previous time step
        Matrix3<double> Fn  = (*FE)*(*FP);
        // New elastic deformation gradient projected to the yield surface
        // proj(Fᴱ)
        *FE = U*proj_F_hat.asDiagonal()*V.transpose();
        // By F = proj(Fᴱ)proj(Fᴾ) = FᴱFᴾ, calculated the new plastic
        // deformation gradient
        *FP = (*FE).householderQr().solve(Fn);
    }
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake


