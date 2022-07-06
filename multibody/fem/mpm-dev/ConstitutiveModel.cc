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

double ConstitutiveModel::get_mu() const {  return mu_;  }

double ConstitutiveModel::get_lambda() const {  return lambda_;  }

}  // namespace mpm
}  // namespace multibody
}  // namespace drake

