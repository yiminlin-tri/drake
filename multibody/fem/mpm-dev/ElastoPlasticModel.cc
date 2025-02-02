#include "drake/multibody/fem/mpm-dev/ElastoPlasticModel.h"

namespace drake {
namespace multibody {
namespace mpm {

ElastoPlasticModel::ElastoPlasticModel(): ElastoPlasticModel(9e4, 0.49) {}

ElastoPlasticModel::ElastoPlasticModel(double E, double nu):
                                                lambda_(E*nu/(1+nu)/(1-2*nu)),
                                                mu_(E/(2*(1+nu))) {
    DRAKE_ASSERT(E >= 0);
    DRAKE_ASSERT(nu > -1.0 && nu < 0.5);
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake

