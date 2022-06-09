#include "drake/multibody/fem/mpm-dev/AnalyticLevelSet.h"

namespace drake {
namespace multibody {
namespace mpm {

AnalyticLevelSet::AnalyticLevelSet(const Vector3<double>& center,
                                   double volume): center_(center),
                                                   volume_(volume) {}

double AnalyticLevelSet::get_volume() const {
    return volume_;
}

SphereLevelSet::SphereLevelSet(const Vector3<double>& center, double radius):
                    AnalyticLevelSet(center, 4.0/3.0*M_PI*radius*radius*radius),
                    radius_(radius) {}

bool SphereLevelSet::InInterior(const Vector3<double>& position) const {
    return (position-center_).norm() < radius_;
}

BoxLevelSet::BoxLevelSet(const Vector3<double>& center,
                         const Vector3<double>& xscale,
                         const math::RotationMatrix<double>& rotation):
                    AnalyticLevelSet(center, 8*xscale(0)*xscale(1)*xscale(2)),
                    xscale_(xscale), rotation_(rotation) {
    DRAKE_ASSERT(xscale(0) >= 0);
    DRAKE_ASSERT(xscale(1) >= 0);
    DRAKE_ASSERT(xscale(2) >= 0);
}

bool BoxLevelSet::InInterior(const Vector3<double>& position) const {
    Vector3<double> normalized_coordinate =
                                    rotation_.inverse()*(position-center_);
    return ((std::abs(normalized_coordinate(0)) < xscale_(0))
         && (std::abs(normalized_coordinate(1)) < xscale_(1))
         && (std::abs(normalized_coordinate(2)) < xscale_(2)));
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
