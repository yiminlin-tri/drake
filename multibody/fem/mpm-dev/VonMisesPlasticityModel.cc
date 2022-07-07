#include "drake/multibody/fem/mpm-dev/VonMisesPlasticityModel.h"

namespace drake {
namespace multibody {
namespace mpm {

VonMisesPlasticityModel::VonMisesPlasticityModel(double tau_c): tau_c_(tau_c) {
    DRAKE_ASSERT(tau_c >= 0);
}

void VonMisesPlasticityModel::UpdateDeformationGradients(double mu,
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
    // Singular values of the deviatoric component of Kirchhoff stress:
    // dev(τ) = τ - pJI = τ - 1/3tr(τ)I
    // σ(dev(τ)) = 2μlog(σᵢ) - 2μ/3 ∑ᵢ log(σᵢ)
    Vector3<double> tau_dev_hat = 2*mu*(eps_hat-1.0/3.0*tr_eps_vec);
    double tau_dev_hat_norm = tau_dev_hat.norm();

    // If the trial stress τ is in the yield surface f(τ) <= 0, plasticity is
    // not applied.
    // Otherwise, project the trial stress τ in the plastic flow direction
    // df/dτ τ onto the yield surface f(τ) <= 0
    // The trial stress is on the yield surface, f(τ) <= 0, if and only if the
    // singular values of the deviatoric component of trial stress satisfies:
    // f(τ) = sqrt(3/2)‖ σ(dev(τ)) ‖ - τ_c ≤ 0
    double sqrt_32 = sqrt(3.0/2.0);
    double f_tau = sqrt_32*tau_dev_hat_norm - tau_c_;
    bool in_yield_surface =  f_tau <= 0.0;
    if (!in_yield_surface) {
        // The singular values of trial strain's projection onto yield surface
        // σ(ε_proj) =  σ(ε) - Δγ ν, where
        // ν = 1/√(2/3)*σ(dev(τ))/‖σ(dev(τ))‖,
        // Δγ = f(dev(τ))/(3μ)
        // By the definition of Hencky strain,
        // σ(F_proj) = exp(σ(ε_proj))
        Vector3<double> nu = sqrt_32*tau_dev_hat/tau_dev_hat_norm;
        double delta_gamma = f_tau/(3.0*mu);
        Vector3<double> proj_F_hat = (eps_hat - delta_gamma*nu).array().exp();
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


