#include "drake/multibody/fem/mpm-dev/BoundaryCondition.h"

namespace drake {
namespace multibody {
namespace mpm {

BoundaryCondition::BoundaryCondition(
            const std::vector<Boundary>& boundaries) {
    boundaries_ = boundaries;
}

int BoundaryCondition::get_num_boundary() const {
    return boundaries_.size();
}

const std::vector<BoundaryCondition::Boundary>&
                                    BoundaryCondition::get_boundaries() const {
    return boundaries_;
}

const BoundaryCondition::Boundary&
                            BoundaryCondition::get_boundary(int index) const {
    DRAKE_ASSERT(index < boundaries_.size());
    return boundaries_[index];
}

void BoundaryCondition::AddBoundary(const
                                    BoundaryCondition::Boundary& boundary) {
    boundaries_.emplace_back(boundary);
}



}  // namespace mpm
}  // namespace multibody
}  // namespace drake
