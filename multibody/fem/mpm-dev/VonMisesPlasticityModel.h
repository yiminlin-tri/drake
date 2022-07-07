#pragma once

#include <memory>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/matrix_utilities.h"

namespace drake {
namespace multibody {
namespace mpm {

// A implementation of von Mises plasticity model, described in
// https://www.proquest.com/docview/2389768700
// The return mapping algorithm can be found in Section 3.3.2.2
class VonMisesPlasticityModel {
 public:
    // tau_c: maximum allowed tensile strength
    explicit VonMisesPlasticityModel(double tau_c);

    std::unique_ptr<VonMisesPlasticityModel> Clone() const {
        return std::make_unique<VonMisesPlasticityModel>(*this);
    }

    void UpdateDeformationGradients(double mu,
                                Matrix3<double>* elastic_deformation_gradient,
                                Matrix3<double>* plastic_deformation_gradient);

 private:
    double tau_c_;
};  // class VonMisesPlasticityModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
