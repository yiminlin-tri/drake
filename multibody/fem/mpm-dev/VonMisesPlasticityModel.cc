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
    VectorX<double> sigma    = svd.singularValues();
    // Singular values of the trial strain
    Vector3<double> eps_hat = {log(sigma(0)), log(sigma(1)), log(sigma(2))};
    // trace of the trial Hencky strain
    double tr_eps = 0.0;
    for (int d = 0; d < 3; ++d) {
        tr_eps += eps_hat(d);
    }
    Vector3<double> tr_eps_vec = tr_eps*Vector3<double>::Ones();
    // Singular values of the trial stress
    Vector3<double> tau_hat = lambda*tr_eps_vec + 2*mu*eps_hat;
    // trace of the trial stress
    double tr_tau = 0.0;
    for (int d = 0; d < 3; ++d) {
        tr_tau += tau_hat(d);
    }
    Vector3<double> tr_tau_vec = tr_tau*Vector3<double>::Ones();

    // If the trial stress is in the yield surface, nothing is done.
    // Otherwise, project the trial stress in the plastic flow direction onto
    // the yield surface
    Vector3<double> proj_eps_hat;
    bool in_yield_surface = (tau_hat - tr_tau_vec).norm() <= tau_c_;
    if (!in_yield_surface) {
        Matrix3<double> C;
        double diag_entry = 2.0*mu+lambda;
        C << diag_entry, lambda, lambda,
             lambda, diag_entry, lambda,
             lambda, lambda, diag_entry;
        Vector3<double> p = 1.0/3.0*tr_tau_vec;
        Vector3<double> d = tau_hat - p;
        proj_eps_hat = C.householderQr().solve(p+tau_c_*d/d.norm());
        for (int c = 0; c < 3; ++c) {
            proj_eps_hat(c) = exp(proj_eps_hat(c));
        }
        // Total deformation gradient in the previous time step
        Matrix3<double> Fn  = (*FE)*(*FP);
        *FE = U*proj_eps_hat.asDiagonal()*V.transpose();
        *FP = (*FE).householderQr().solve(Fn);
    }
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake


