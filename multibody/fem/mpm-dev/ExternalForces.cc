#include "drake/multibody/fem/mpm-dev/ExternalForces.h"

namespace drake {
namespace multibody {
namespace mpm {

void ExternalForces::ApplyGravitationalForces(double dt, Grid* grid) const {
    int i, j, k;
    // Gravitational acceleration
    const Vector3<double> g = {0.0, 0.0, -9.81*dt};
    for (const auto& [batch_index_flat, batch_index_3d] : grid->get_indices()) {
        i = batch_index_3d(0);
        j = batch_index_3d(1);
        k = batch_index_3d(2);
        const Vector3<double>& velocity_i = grid->get_velocity(i, j, k);
        grid->set_velocity(i, j, k, velocity_i + g);
    }
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake

