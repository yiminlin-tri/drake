#pragma once

#include <math.h>

#include <algorithm>
#include <array>
#include <limits>

#include "drake/common/eigen_types.h"
#include "drake/math/rotation_matrix.h"

namespace drake {
namespace multibody {
namespace mpm {

// A base class providing the interface of primitive geometries' level set in
// reference configuration
class AnalyticLevelSet {
 public:
    AnalyticLevelSet(double volume,
                     const std::array<Vector3<double>, 2>& bounding_box);

    // Return true if the position is in the interiror of the level set.
    virtual bool InInterior(const Vector3<double>& position) const = 0;

    // Return the outward unit normal of the analytic level set,
    // return 0 if the position is in the interior or outside of the object
    virtual Vector3<double> Normal(const Vector3<double>& position)
                                                                     const = 0;

    // Return the volume enclosed by the level set
    double get_volume() const;

    // Return the bounding bound of the geometry
    // We denote the first component of bounding_box_ as xmin_, second component
    // is xmax_. xmin_, xmax_ represents the bounding box of the geometry
    // i.e. the geometry lies in
    // [xmin_(0), xmax_(0)]X[xmin_(1), xmax_(1)]X[xmin_(2), xmax_(2)]
    const std::array<Vector3<double>, 2>& get_bounding_box() const;

    virtual ~AnalyticLevelSet() = default;

 protected:
    double volume_;
    std::array<Vector3<double>, 2> bounding_box_{};
};  // class AnalyticLevelSet

// An analytic level set class for sphere with radius radius_ and center
// (0, 0, 0)
class SphereLevelSet : public AnalyticLevelSet {
 public:
    explicit SphereLevelSet(double radius);
    bool InInterior(const Vector3<double>& position) const final;
    Vector3<double> Normal(const Vector3<double>& position) const final;

 private:
    double radius_;
};  // class SphereLevelSet

// An analytic level set class for box of size
// [-xscale(0), xscale(0)] X [-xscale(1), xscale(1)] X [xscale(2), xscale(2)]
// centered at (0, 0, 0)
class BoxLevelSet : public AnalyticLevelSet {
 public:
    explicit BoxLevelSet(const Vector3<double>& xscale);
    bool InInterior(const Vector3<double>& position) const final;
    Vector3<double> Normal(const Vector3<double>& position) const final;

 private:
    Vector3<double> xscale_{};
};  // class BoxLevelSet

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
