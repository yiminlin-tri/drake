#include "drake/multibody/fem/mpm-dev/StvkHenckyWithVonMisesModel.h"

namespace drake {
namespace multibody {
namespace mpm {

StvkHenckyWithVonMisesModel::StvkHenckyWithVonMisesModel(double tau_c):
                                    ElastoPlasticModel(), yield_stress_(tau_c) {
    DRAKE_ASSERT(tau_c >= 0);
}

StvkHenckyWithVonMisesModel::StvkHenckyWithVonMisesModel(double E, double nu,
                                                         double tau_c):
                                                    ElastoPlasticModel(E, nu),
                                                    yield_stress_(tau_c) {
    DRAKE_ASSERT(tau_c >= 0);
}

void StvkHenckyWithVonMisesModel::CalcFirstPiolaKirchhoffStress(
        const Matrix3<double>& F, Matrix3<double>* P) const {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U  = svd.matrixU();
    const Matrix3<double>& V  = svd.matrixV();
    Vector3<double> sigma     = svd.singularValues();
    // logs of singular values
    Vector3<double> log_sigma = sigma.array().log();
    // sum of logs of singular values of F: ∑ ln(σᵢ)
    double sum_log_sigma      = log_sigma.sum();
    // Overwrite the the vector sigma [σᵢ] with the singular values of the first
    // Piola Kirchhoff stress:
    // Σᵢ = 1/σᵢ ⋅ [2μ log(σᵢ) + λ ∑ᵢ log(σᵢ)]
    for (int d = 0; d < 3; ++d) {
        sigma(d) = (2*mu_*log_sigma(d) + lambda_*sum_log_sigma)/sigma(d);
    }
    // The first Piola Kirchhoff stress can be then written as
    // P = U Σᵢ Vᵀ, where U, V are left and right singular vectors of the
    //             deformation gradient F
    *P = U * sigma.asDiagonal() * V.transpose();
}

void StvkHenckyWithVonMisesModel::CalcKirchhoffStress(const Matrix3<double>& F,
                                         Matrix3<double>* tau) const {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U = svd.matrixU();
    Vector3<double> sigma    = svd.singularValues();
    // Since the Hencky strain has form ε = 1/2ln(Fᴱ Fᴱ^T)
    // the Hencky strain in the principal frame is ln(Σ), where Fᴱ = U Σ Vᵀ
    Vector3<double> eps_hat = sigma.array().log();
    // trace of the trial Hencky strain tr(ε)
    double tr_eps           = eps_hat.sum();

    CalcKirchhoffStress(U, eps_hat, tr_eps, tau);
}

void StvkHenckyWithVonMisesModel::
                                CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
        const Matrix3<double>& F, Matrix3<double>* P,
        Matrix3<double>* tau) const {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U = svd.matrixU();
    const Matrix3<double>& V = svd.matrixV();
    Vector3<double> sigma    = svd.singularValues();
    // logs of singular values
    Vector3<double> log_sigma = sigma.array().log();
    // sum of logs of singular values of F: ∑ ln(σᵢ)
    double sum_log_sigma      = log_sigma.sum();
    // Overwrite the the vector sigma [σᵢ] with the singular values of the first
    // Piola Kirchhoff stress:
    // Σᵢ = 1/σᵢ ⋅ [2μ log(σᵢ) + λ ∑ᵢ log(σᵢ)]
    for (int d = 0; d < 3; ++d) {
        sigma(d) = (2*mu_*log_sigma(d) + lambda_*sum_log_sigma)/sigma(d);
    }
    // The first Piola Kirchhoff stress can be then written as
    // P = U Σᵢ Vᵀ, where U, V are left and right singular vectors of the
    //             deformation gradient F, and τ = PFᵀ
    *P   = U * sigma.asDiagonal() * V.transpose();
    *tau = (*P)*F.transpose();
}

void StvkHenckyWithVonMisesModel::UpdateDeformationGradient(
                                                    Matrix3<double>* FE) const {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(*FE, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U = svd.matrixU();
    const Matrix3<double>& V = svd.matrixV();
    Vector3<double> sigma    = svd.singularValues();
    // Since the Hencky strain has form ε = 1/2ln(Fᴱ Fᴱ^T)
    // the Hencky strain in the principal frame is ln(Σ), where Fᴱ = U Σ Vᵀ
    Vector3<double> eps_hat = sigma.array().log();
    // trace of the trial Hencky strain tr(ε)
    double tr_eps           = eps_hat.sum();

    ProjectDeformationGradientToYieldSurface(U, V, eps_hat, tr_eps, FE);
}

void StvkHenckyWithVonMisesModel::
                                CalcKirchhoffStressAndUpdateDeformationGradient(
                        Matrix3<double>* tau, Matrix3<double>* FE) const {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(*FE, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U = svd.matrixU();
    const Matrix3<double>& V = svd.matrixV();
    Vector3<double> sigma    = svd.singularValues();
    // Since the Hencky strain has form ε = 1/2ln(Fᴱ Fᴱ^T)
    // the Hencky strain in the principal frame is ln(Σ), where Fᴱ = U Σ Vᵀ
    Vector3<double> eps_hat = sigma.array().log();
    // trace of the trial Hencky strain tr(ε)
    double tr_eps           = eps_hat.sum();

    CalcKirchhoffStress(U, eps_hat, tr_eps, tau);
    ProjectDeformationGradientToYieldSurface(U, V, eps_hat, tr_eps, FE);
}

void StvkHenckyWithVonMisesModel::CalcKirchhoffStress(
                                const Matrix3<double>& U,
                                const Vector3<double>& eps_hat,
                                double tr_eps, Matrix3<double>* tau) const {
    // Calculate the singular values of Kirchhoff stress
    // Σᵢ = 2μ log(σᵢ) + λ ∑ᵢ log(σᵢ)
    Vector3<double> sigma_tau;
    for (int d = 0; d < 3; ++d) {
        sigma_tau(d)  = 2*mu_*eps_hat(d) + lambda_*tr_eps;
    }
    // The Kirchhoff stress can be then written as
    // τ = U Σᵢ Uᵀ, where U are left singular vectors of the deformation
    //             gradient F
    *tau = U * sigma_tau.asDiagonal() * U.transpose();
}

void StvkHenckyWithVonMisesModel::
        ProjectDeformationGradientToYieldSurface(const Matrix3<double>& U,
                                                 const Matrix3<double>& V,
                                                 const Vector3<double>& eps_hat,
                                                 double tr_eps,
                                                 Matrix3<double>* FE) const {
    // Vector of trace of the trial stress
    Vector3<double> tr_eps_vec = tr_eps*Vector3<double>::Ones();
    // The deviatoric component of Kirchhoff stress in the principal frame is:
    // dev(τ) = τ - pJI = τ - 1/3tr(τ)I
    // In the principal frame: dev(τ) = 2μlog(σᵢ) - 2μ/3 ∑ᵢ log(σᵢ)
    Vector3<double> tau_dev_hat = 2*mu_*(eps_hat-1.0/3.0*tr_eps_vec);
    double tau_dev_hat_norm = tau_dev_hat.norm();

    // If the trial stress τ is in the yield surface f(τ) <= 0, plasticity is
    // not applied.
    // Otherwise, project the trial stress τ in the plastic flow direction
    // df/dτ τ onto the yield surface f(τ) <= 0
    // The trial stress is on the yield surface, f(τ) <= 0, if and only if the
    // singular values of the deviatoric component of trial stress satisfies:
    // f(τ) = sqrt(3/2)‖ dev(τ) ‖ - τ_c ≤ 0
    double sqrt_32 = sqrt(3.0/2.0);
    double f_tau = sqrt_32*tau_dev_hat_norm - yield_stress_;
    bool in_yield_surface =  f_tau <= 0.0;
    if (!in_yield_surface) {
        // Trial strain's projection onto yield surface in the principal frame
        // is:
        // ε_proj =  ε - Δγ νⁿ⁺¹,
        //              where ν = νⁿ⁺¹ = 1/√(2/3)*dev(τⁿ⁺¹)/‖dev(τⁿ⁺¹)‖
        //                            = 1/√(2/3)*dev(τ)/‖dev(τ)‖
        //              Since dev(τ) and dev(τⁿ⁺¹) are in the same direction
        // Taking the dot product of the associative flow rule with flow
        // direction ν:
        //     ν [dev(τⁿ⁺¹) - dev(τ)] = -2μ Δγ ‖ν‖
        // ==>         f(τⁿ⁺¹) - f(τ) = -3μ Δγ
        // Since f(τⁿ⁺¹) = ‖ τⁿ⁺¹ ‖ - τ_c = 0, i.e. τⁿ⁺¹ is on the yield surface
        // Δγ = f(τ)/(3μ)
        // By the definition of Hencky strain,
        // F_proj = exp(ε_proj) in the principal frame
        Vector3<double> nu = sqrt_32*tau_dev_hat/tau_dev_hat_norm;
        double delta_gamma = f_tau/(3.0*mu_);
        Vector3<double> proj_F_hat = (eps_hat - delta_gamma*nu).array().exp();
        // New elastic deformation gradient projected to the yield surface
        // proj(Fᴱ)
        *FE = U*proj_F_hat.asDiagonal()*V.transpose();
    }
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
