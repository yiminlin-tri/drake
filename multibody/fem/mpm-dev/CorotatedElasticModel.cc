#include "drake/multibody/fem/mpm-dev/CorotatedElasticModel.h"

namespace drake {
namespace multibody {
namespace mpm {

CorotatedElasticModel::CorotatedElasticModel(): ElastoPlasticModel() {}

CorotatedElasticModel::CorotatedElasticModel(double E, double nu):
                                                ElastoPlasticModel(E, nu) {}

void CorotatedElasticModel::CalcKirchhoffStress(const Matrix3<double>& FE,
                                         Matrix3<double>* tau) const {
    Matrix3<double> R, S;
    double J = FE.determinant();
    fem::internal::PolarDecompose<double>(FE, &R, &S);
    *tau = 2.0*mu_*(FE-R)*FE.transpose()
         + lambda_*(J-1.0)*J*Matrix3<double>::Identity();
}

void CorotatedElasticModel::UpdateDeformationGradientAndCalcKirchhoffStress(
                        Matrix3<double>* tau, Matrix3<double>* FE) const {
    CalcKirchhoffStress(*FE, tau);
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
