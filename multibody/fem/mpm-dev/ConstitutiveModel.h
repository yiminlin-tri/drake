#pragma once

#include <iostream>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/matrix_utilities.h"

namespace drake {
namespace multibody {
namespace mpm {

// A base class providing the interface of Constitutive model
class ConstitutiveModel {
 public:
    // Default material, dough: E = 9e4 Pa, nu = 0.49
    ConstitutiveModel();

    // Constructor uses Young's modulus E and Poisson's ratio nu
    ConstitutiveModel(double E, double nu);

    // First Piola Kirchhoff stress density: P = dpsi/dF
    virtual void CalcFirstPiolaKirchhoffStress(
            const Matrix3<double>& F, EigenPtr<Matrix3<double>> P) const = 0;

    // Kirchhoff stress density: tau = P F^T = dpsi/dF F^T
    virtual void CalcKirchhoffStress(const Matrix3<double>& F,
                                     EigenPtr<Matrix3<double>> tau) const = 0;
    virtual void CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
                        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P,
                        EigenPtr<Matrix3<double>> tau) const = 0;

    virtual ~ConstitutiveModel() = default;

 protected:
    double mu_;                         // Parameters in defining
    double lambda_;                     // the energy density function
};  // class ConstitutiveModel

// A implementation of Fixed Corotated Model (Constitutive Model)
class CorotatedModel : public ConstitutiveModel {
 public:
    CorotatedModel();
    CorotatedModel(double E, double nu);

    void CalcFirstPiolaKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P) const final;

    void CalcKirchhoffStress(const Matrix3<double>& F,
                             EigenPtr<Matrix3<double>> tau) const final;
    void CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P,
        EigenPtr<Matrix3<double>> tau) const final;
};  // class CorotatedModel

// A implementation of Saint-Venant Kirchhoff model, but replace the left
// Cauchy Green strain with the Hencky strain
// https://dl.acm.org/doi/abs/10.1145/2897824.2925906
// The formula of Kirchhoff stress can be found in Section 6.3
class SaintVenantKirchhoffWithHenckyModel : public ConstitutiveModel {
 public:
    SaintVenantKirchhoffWithHenckyModel();
    SaintVenantKirchhoffWithHenckyModel(double E, double nu);

    void CalcFirstPiolaKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P) const final;

    void CalcKirchhoffStress(const Matrix3<double>& F,
                             EigenPtr<Matrix3<double>> tau) const final;
    void CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P,
        EigenPtr<Matrix3<double>> tau) const final;
};  // class SaintVenantKirchhoffWithHenckyModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
