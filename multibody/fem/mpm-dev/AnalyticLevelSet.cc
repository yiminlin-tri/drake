#include "drake/multibody/fem/mpm-dev/AnalyticLevelSet.h"

namespace drake {
namespace multibody {
namespace mpm {

AnalyticLevelSet::AnalyticLevelSet(double volume,
                            const std::array<Vector3<double>, 2>& bounding_box):
                                volume_(volume), bounding_box_(bounding_box) {}

double AnalyticLevelSet::get_volume() const {
    return volume_;
}

const std::array<Vector3<double>, 2>& AnalyticLevelSet::get_bounding_box()
                                                                        const {
    return bounding_box_;
}

SphereLevelSet::SphereLevelSet(double radius):
                    AnalyticLevelSet(4.0/3.0*M_PI*radius*radius*radius,
                                    {{-radius*Vector3<double>::Ones(),
                                       radius*Vector3<double>::Ones()}}),
                    radius_(radius) {}

bool SphereLevelSet::InInterior(const Vector3<double>& position) const {
    return position.norm() < radius_;
}

Vector3<double> SphereLevelSet::Normal(const Vector3<double>& position) const {
    if (position.norm() == radius_) {
        return position.normalized();
    } else {
        return Vector3<double>::Zero();
    }
}

BoxLevelSet::BoxLevelSet(const Vector3<double>& xscale):
                    AnalyticLevelSet(8*xscale(0)*xscale(1)*xscale(2),
                                    {{-xscale, xscale}}),
                    xscale_(xscale) {
    DRAKE_ASSERT(xscale_(0) >= 0);
    DRAKE_ASSERT(xscale_(1) >= 0);
    DRAKE_ASSERT(xscale_(2) >= 0);
}

bool BoxLevelSet::InInterior(const Vector3<double>& position) const {
    return ((std::abs(position(0)) < xscale_(0))
         && (std::abs(position(1)) < xscale_(1))
         && (std::abs(position(2)) < xscale_(2)));
}

Vector3<double> BoxLevelSet::Normal(const Vector3<double>& position) const {
    // Imagine the box is aligned like:
    //          -------
    //         /     /|
    //         -----  |
    //  left  |     | |   right
    //        |     |/                  y
    //        ------               z | /
    //         bottom                  - x
    bool on_left_plane   = position(0) == -xscale_(0);
    bool on_right_plane  = position(0) ==  xscale_(0);
    bool on_front_plane  = position(1) == -xscale_(1);
    bool on_back_plane   = position(1) ==  xscale_(1);
    bool on_bottom_plane = position(2) == -xscale_(2);
    bool on_top_plane    = position(2) ==  xscale_(2);
    bool on_boundary = ((on_left_plane   || on_right_plane ||
                         on_front_plane  || on_back_plane  ||
                         on_bottom_plane || on_top_plane)
                     && (std::abs(position(0)) <= xscale_(0))
                     && (std::abs(position(1)) <= xscale_(1))
                     && (std::abs(position(2)) <= xscale_(2)));
    if (on_boundary) {
        if (on_left_plane) {
            return Vector3<double>{-1.0, 0.0, 0.0};
        }
        if (on_right_plane) {
            return Vector3<double>{1.0, 0.0, 0.0};
        }
        if (on_front_plane) {
            return Vector3<double>{0.0, -1.0, 0.0};
        }
        if (on_back_plane) {
            return Vector3<double>{0.0, 1.0, 0.0};
        }
        if (on_bottom_plane) {
            return Vector3<double>{0.0, 0.0, -1.0};
        }
        if (on_top_plane) {
            return Vector3<double>{0.0, 0.0, 1.0};
        }
    }
    return Vector3<double>::Zero();
}
}  // namespace mpm
}  // namespace multibody
}  // namespace drake
