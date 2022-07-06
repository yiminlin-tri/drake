#pragma once

#include <memory>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/matrix_utilities.h"
#include "drake/multibody/fem/mpm-dev/ConstitutiveModel.h"

namespace drake {
namespace multibody {
namespace mpm {

// A implementation of Fixed Corotated Model (Constitutive Model)
class CorotatedModel : public ConstitutiveModel {
 public:
    CorotatedModel();
    CorotatedModel(double E, double nu);

    std::unique_ptr<ConstitutiveModel> Clone() const {
        return std::make_unique<CorotatedModel>(*this);
    }

    void CalcFirstPiolaKirchhoffStress(
        const Matrix3<double>& F, Matrix3<double>* P) const final;

    void CalcKirchhoffStress(const Matrix3<double>& F,
                             Matrix3<double>* tau) const final;
    void CalcFirstPiolaKirchhoffStressAndKirchhoffStress(
        const Matrix3<double>& F, Matrix3<double>* P,
        Matrix3<double>* tau) const final;
};  // class CorotatedModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
