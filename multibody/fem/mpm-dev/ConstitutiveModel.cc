#include "drake/multibody/fem/mpm-dev/ConstitutiveModel.h"

namespace drake {
namespace multibody {
namespace mpm {

ConstitutiveModel::ConstitutiveModel(): ConstitutiveModel(9e4, 0.49) {}

ConstitutiveModel::ConstitutiveModel(double E, double nu):
                                                mu_(E/(2*(1+nu))),
                                                lambda_(E*nu/(1+nu)/(1-2*nu)) {
    DRAKE_ASSERT(E >= 0);
    DRAKE_ASSERT(nu > -1.0 && nu < 0.5);
}

CorotatedModel::CorotatedModel(): ConstitutiveModel() {}

CorotatedModel::CorotatedModel(double E, double nu): ConstitutiveModel(E, nu) {}

void CorotatedModel::CalcFirstPiolaKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P) const {
    Matrix3<double> R, S, JFinvT;
    double J = F.determinant();
    fem::internal::PolarDecompose<double>(F, &R, &S);
    fem::internal::CalcCofactorMatrix<double>(F, &JFinvT);
    *P = 2.0*mu_*(F-R) + lambda_*(J-1.0)*JFinvT;
}

void CorotatedModel::CalcKirchhoffStress(const Matrix3<double>& F,
                                         EigenPtr<Matrix3<double>> tau) const {
    Matrix3<double> R, S;
    double J = F.determinant();
    fem::internal::PolarDecompose<double>(F, &R, &S);
    *tau = 2.0*mu_*(F-R)*F.transpose()
         + lambda_*(J-1.0)*J*Matrix3<double>::Identity();
}

void CorotatedModel::CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P,
        EigenPtr<Matrix3<double>> tau) const {
    Matrix3<double> R, S, JFinvT;
    double J = F.determinant();
    fem::internal::PolarDecompose<double>(F, &R, &S);
    fem::internal::CalcCofactorMatrix<double>(F, &JFinvT);
    *P = 2.0*mu_*(F-R) + lambda_*(J-1.0)*JFinvT;
    *tau = (*P)*F.transpose();
}

SaintVenantKirchhoffWithHenckyModel::SaintVenantKirchhoffWithHenckyModel():
                                                        ConstitutiveModel() {}

SaintVenantKirchhoffWithHenckyModel::SaintVenantKirchhoffWithHenckyModel(
                                                        double E, double nu):
                                                    ConstitutiveModel(E, nu) {}

void SaintVenantKirchhoffWithHenckyModel::CalcFirstPiolaKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P) const {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U = svd.matrixU();
    const Matrix3<double>& V = svd.matrixV();
    VectorX<double> sigma    = svd.singularValues();
    // sum of logs of singular values of F: ∑ ln(σᵢ)
    double sum_log_sigma = 0.0;
    for (int d = 0; d < 3; ++d) {
        sum_log_sigma += log(sigma(d));
    }
    for (int d = 0; d < 3; ++d) {
        sigma(d) = (2*mu_*log(sigma(d)) + lambda_*sum_log_sigma)/sigma(d);
    }
    *P = U * sigma.asDiagonal() * V.transpose();
}

void SaintVenantKirchhoffWithHenckyModel::CalcKirchhoffStress
                                        (const Matrix3<double>& F,
                                         EigenPtr<Matrix3<double>> tau) const {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U = svd.matrixU();
    VectorX<double> sigma    = svd.singularValues();
    // sum of logs of singular values of F: ∑ ln(σᵢ)
    double sum_log_sigma = 0.0;
    for (int d = 0; d < 3; ++d) {
        sum_log_sigma += log(sigma(d));
    }
    for (int d = 0; d < 3; ++d) {
        sigma(d)  = 2*mu_*log(sigma(d)) + lambda_*sum_log_sigma;
    }
    *tau = U * sigma.asDiagonal() * U.transpose();
}

void SaintVenantKirchhoffWithHenckyModel::
                                CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P,
        EigenPtr<Matrix3<double>> tau) const {
    Eigen::JacobiSVD<Matrix3<double>>
                            svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix3<double>& U = svd.matrixU();
    const Matrix3<double>& V = svd.matrixV();
    VectorX<double> sigma    = svd.singularValues();
    // sum of logs of singular values of F: ∑ ln(σᵢ)
    double sum_log_sigma = 0.0;
    for (int d = 0; d < 3; ++d) {
        sum_log_sigma += log(sigma(d));
    }
    for (int d = 0; d < 3; ++d) {
        sigma(d) = (2*mu_*log(sigma(d)) + lambda_*sum_log_sigma)/sigma(d);
    }
    *P   = U * sigma.asDiagonal() * V.transpose();
    *tau = (*P)*F.transpose();
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake

