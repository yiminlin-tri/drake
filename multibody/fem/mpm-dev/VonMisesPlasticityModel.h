#pragma once

#include <memory>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/matrix_utilities.h"

namespace drake {
namespace multibody {
namespace mpm {

// A implementation of von Mises plasticity model, described in
// https://www.klancek.si/sites/default/files/datoteke/files/
// bonet-woodnonlinearcontinuummechanics2ndedition.pdf
// The return mapping algorithm can be found in Box 7.1
class VonMisesPlasticityModel {
 public:
    // tau_c: maximum allowed tensile strength
    explicit VonMisesPlasticityModel(double tau_c);

    std::unique_ptr<VonMisesPlasticityModel> Clone() const {
        return std::make_unique<VonMisesPlasticityModel>(*this);
    }

    // Update the elastic and plastic deformation gradient according to the
    // Von Mises plasticity model. mu is the Lame coefficient of the input
    // Hencky constitutive model
    void UpdateDeformationGradients(double mu,
                                Matrix3<double>* elastic_deformation_gradient,
                                Matrix3<double>* plastic_deformation_gradient);

 private:
    double tau_c_;
};  // class VonMisesPlasticityModel

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
