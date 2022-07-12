#include "drake/multibody/fem/mpm-dev/CorotatedElasticModel.h"

#include "drake/common/unused.h"

namespace drake {
namespace multibody {
namespace mpm {

CorotatedElasticModel::CorotatedElasticModel(): ElastoPlasticModel() {}

CorotatedElasticModel::CorotatedElasticModel(double E, double nu):
                                                ElastoPlasticModel(E, nu) {}

void CorotatedElasticModel::CalcFirstPiolaKirchhoffStress(
        const Matrix3<double>& F, Matrix3<double>* P) const {
    Matrix3<double> R, S, JFinvT;
    double J = F.determinant();
    fem::internal::PolarDecompose<double>(F, &R, &S);
    fem::internal::CalcCofactorMatrix<double>(F, &JFinvT);
    *P = 2.0*mu_*(F-R) + lambda_*(J-1.0)*JFinvT;
}

void CorotatedElasticModel::CalcKirchhoffStress(const Matrix3<double>& F,
                                         Matrix3<double>* tau) const {
    Matrix3<double> R, S;
    double J = F.determinant();
    fem::internal::PolarDecompose<double>(F, &R, &S);
    *tau = 2.0*mu_*(F-R)*F.transpose()
         + lambda_*(J-1.0)*J*Matrix3<double>::Identity();
}

void CorotatedElasticModel::CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
        const Matrix3<double>& F, Matrix3<double>* P,
        Matrix3<double>* tau) const {
    Matrix3<double> R, S, JFinvT;
    double J = F.determinant();
    fem::internal::PolarDecompose<double>(F, &R, &S);
    fem::internal::CalcCofactorMatrix<double>(F, &JFinvT);
    *P = 2.0*mu_*(F-R) + lambda_*(J-1.0)*JFinvT;
    *tau = (*P)*F.transpose();
}

void CorotatedElasticModel::UpdateDeformationGradient(
                    Matrix3<double>* elastic_deformation_gradient) const {
    unused(elastic_deformation_gradient);
}

void CorotatedElasticModel::CalcKirchhoffStressAndUpdateDeformationGradient(
                        Matrix3<double>* tau, Matrix3<double>* FE) const {
    CalcKirchhoffStress(*FE, tau);
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
