#pragma once

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/matrix_utilities.h"

namespace drake {
namespace multibody {
namespace mpm {

// A implementation of 3D quadratic B spline
class CorotatedModel {
 public:
    // Default material: Young's modulus 2.0, Poisson's ratio 0.0
    CorotatedModel();
    // Constructor uses Young's modulus E and Poisson's ratio nu
    CorotatedModel(const double E, const double nu);

    // First Piola Kirchhoff stress: P = dpsi/dF
    void CalcFirstPiolaKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P);

    // Kirchhoff stress: tau = P F^T = dpsi/dF F^T
    void CalcKirchhoffStress(const Matrix3<double>& F,
                             EigenPtr<Matrix3<double>> tau);
    // First Piola Kirchhoff stress: P = dpsi/dF
    void CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
        const Matrix3<double>& F, EigenPtr<Matrix3<double>> P,
        EigenPtr<Matrix3<double>> tau);

 private:
    double mu_;                         // Parameters in defining
    double lambda_;                     // the energy density function
};  // class CorotatedModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
