#pragma once

#include <math.h>

#include "drake/common/eigen_types.h"
#include "drake/math/rotation_matrix.h"

namespace drake {
namespace multibody {
namespace mpm {

// A base class providing the interface of primitive geometries' level set.
class AnalyticLevelSet {
 public:
    explicit AnalyticLevelSet(const Vector3<double>& center, double volume);

    // Return true if the position is in the interiror of the level set.
    virtual bool InInterior(const Vector3<double>& position) const = 0;

    // Return the volume enclosed by the level set
    double get_volume() const;

    virtual ~AnalyticLevelSet() = default;

 protected:
    Vector3<double> center_;
    double volume_;
};  // class AnalyticLevelSet

// An analytic level set class for sphere with radius radius_ and center center_
class SphereLevelSet : public AnalyticLevelSet {
 public:
    SphereLevelSet(const Vector3<double>& center, double radius);

    bool InInterior(const Vector3<double>& position) const;

 private:
    double radius_;
};  // class SphereLevelSet

// An analytic level set class for box of size
// [-xscale(0), xscale(0)] X [-xscale(1), xscale(1)] X [xscale(2), xscale(2)],
// then rotated by rotation_ and translated to center at center_
class BoxLevelSet : public AnalyticLevelSet {
 public:
    BoxLevelSet(const Vector3<double>& center, const Vector3<double>& xscale,
                const math::RotationMatrix<double>& rotation);

    bool InInterior(const Vector3<double>& position) const;

 private:
    Vector3<double> xscale_{};
    math::RotationMatrix<double> rotation_{};
};  // class BoxLevelSet

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
