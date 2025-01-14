#pragma once

#include <limits>
#include <utility>
#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/multibody/fem/mpm-dev/KinematicCollisionObjects.h"
#include "drake/multibody/fem/mpm-dev/TotalMassAndMomentum.h"

namespace drake {
namespace multibody {
namespace mpm {

// A grid class holding vectors of grid points' state. The grid can be
// visually represented as:
//
//                 h
//               o - o - o - o - o - o
//               |   |   |   |   |   |
//           o - o - o - o - o - o - o
//           |   |   |   |   |   |   |
//       o - o - o - o - o - o - o - o
//       |   |   |   |   |   |   |   |
//       o - o - o - o - o - o - o - o
//       |   |   |   |   |   |   |
//       o - o - o - o - o - o - o                  z   y
//       |   |   |   |   |   |                      | /
//       o - o - o - o - o - o                       -- x
//  bottom_corner
//
// The grid can be uniquely represented as a bottom corner and the number of
// grid points in x, y and z directions.
// Since we assume an uniform grid, the distance between neighboring grid points
// are all h. The bottom_corner is the grid node with the smallest (x, y, z)
// values. We assume the bottom corner aligns with the "index space" (i, j, k)
// , which is an 3D integer space (Z \times Z \times Z). The actual physical
// coordinate of the index space is defined through the mapping (ih, jh, kh).
// Positions of grid points will be stored as physical space coordinates.
// The vector num_gridpt_1D stores the number of grid points along x, y, and z
// directions. Figure above is an example with num_gridpt_1D = (6, 3, 4).
// We stores all the data holding states as 1D vectors, with linear indexing
// following the lexiographical order in x, y, z directions
class Grid {
 public:
    Grid() = default;
    Grid(const Vector3<int>& num_gridpt_1D, double h,
         const Vector3<int>& bottom_corner);

    int get_num_gridpt() const;
    const Vector3<int>& get_num_gridpt_1D() const;
    double get_h() const;
    const Vector3<int>& get_bottom_corner() const;

    // For below (i, j, k), we expect users pass in the coordinate in the
    // index space as documented above.
    Vector3<double> get_position(int i, int j, int k) const;
    const Vector3<double>& get_velocity(int i, int j, int k) const;
    double get_mass(int i, int j, int k) const;
    const Vector3<double>& get_force(int i, int j, int k) const;
    const std::vector<std::pair<int, Vector3<int>>>& get_indices() const;

    const Vector3<double>& get_velocity(int idx_flat) const;
    double get_mass(int idx_flat) const;
    const Vector3<double>& get_force(int idx_flat) const;

    void set_velocity(int i, int j, int k, const Vector3<double>& velocity);
    void set_mass(int i, int j, int k, double mass);
    void set_force(int i, int j, int k, const Vector3<double>& force);

    void set_velocity(int idx_flat, const Vector3<double>& velocity);
    void set_mass(int idx_flat, double mass);
    void set_force(int idx_flat, const Vector3<double>& force);

    // Accumulate the state at (i, j, k) with the given value
    void AccumulateVelocity(int i, int j, int k,
                                            const Vector3<double>& velocity);
    void AccumulateMass(int i, int j, int k, double mass);
    void AccumulateForce(int i, int j, int k, const Vector3<double>& force);
    void AccumulateVelocity(const Vector3<int>& index_3d,
                                            const Vector3<double>& velocity);
    void AccumulateMass(const Vector3<int>& index_3d, double mass);
    void AccumulateForce(const Vector3<int>& index_3d,
                         const Vector3<double>& force);

    // Rescale the velocities_ vector by the mass_, used in P2G where we
    // temporarily store momentum mv into velocities
    void RescaleVelocities();

    // Set masses, velocities and forces to be 0
    void ResetStates();

    // Reduce an 3D (i, j, k) index in the index space to a corresponding
    // linear lexiographical ordered index
    int Reduce3DIndex(int i, int j, int k) const;
    int Reduce3DIndex(const Vector3<int>& index_3d) const;

    // Expand a linearly lexiographical ordered index to a 3D index (i, j, k)
    // in the index space
    Vector3<int> Expand1DIndex(int idx) const;

    // Check the passed in (i, j, k) lies within the index range of this Grid
    bool in_index_range(int i, int j, int k) const;
    bool in_index_range(const Vector3<int>& index_3d) const;

    // Assume the explicit grid force at step n f^n is calculated and stored in
    // the member variable, we update the velocity with the formula
    // v^{n+1} = v^n + dt*f^n/m^n
    void UpdateVelocity(double dt);

    // Enforce wall boundary conditions using the given kinematic collision
    // objects. The normal of the given collision object is the outward pointing
    // normal from the interior of the object to the exterior of the object. We
    // strongly impose this Dirichlet boundary conditions.
    void EnforceBoundaryCondition(const KinematicCollisionObjects& objects);

    // Return the sum of mass, momentum and angular momentum of all grid points
    TotalMassAndMomentum GetTotalMassAndMomentum() const;

 private:
    int num_gridpt_;
    Vector3<int> num_gridpt_1D_;              // Number of grid points on the
                                              // grid along x, y, z directions
    double h_;
    Vector3<int> bottom_corner_{};

    // The vector of 1D and 3D indices of grid points, ordered lexiographically
    std::vector<std::pair<int, Vector3<int>>> indices_{};
    std::vector<Vector3<double>> velocities_{};
    std::vector<double> masses_{};
    std::vector<Vector3<double>> forces_{};
};  // class Grid

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
