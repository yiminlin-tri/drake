#pragma once

#include <memory>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/matrix_utilities.h"

namespace drake {
namespace multibody {
namespace mpm {

// A base class providing the interface of constituive and plastic model
class ElastoPlasticModel {
 public:
    // Default material, dough: E = 9e4 Pa, nu = 0.49
    ElastoPlasticModel();

    // Constructor uses Young's modulus E and Poisson's ratio nu
    ElastoPlasticModel(double E, double nu);

    virtual std::unique_ptr<ElastoPlasticModel> Clone() const = 0;

    // Return the first Lamé coefficient
    double get_lambda() const;

    // Return the second Lamé coefficient
    double get_mu() const;

    // First Piola Kirchhoff stress density: P = dpsi/dF
    virtual void CalcFirstPiolaKirchhoffStress(
            const Matrix3<double>& F, Matrix3<double>* P) const = 0;

    // Kirchhoff stress density: tau = P F^T = dpsi/dF F^T
    virtual void CalcKirchhoffStress(const Matrix3<double>& F,
                                     Matrix3<double>* tau) const = 0;
    virtual void CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
                        const Matrix3<double>& F, Matrix3<double>* P,
                        Matrix3<double>* tau) const = 0;

    // Update the elastic deformation gradient according to the plasticity model
    // by projecting the trial elastic stress to the yield surface
    virtual void UpdateDeformationGradient(
                    Matrix3<double>* elastic_deformation_gradient) const = 0;

    virtual void CalcKirchhoffStressAndUpdateDeformationGradient(
                    Matrix3<double>* tau,
                    Matrix3<double>* elastic_deformation_gradient) const = 0;

    virtual ~ElastoPlasticModel() = default;

 protected:
    double mu_;
    double lambda_;
};  // class ElastoPlasticModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
