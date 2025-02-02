#pragma once
#include <utility>
#include <vector>

#include "drake/geometry/geometry_ids.h"
#include "drake/geometry/proximity/deformable_contact_surface.h"
#include "drake/math/rotation_matrix.h"

namespace drake {
namespace geometry {
namespace internal {

// TODO(DamrongGuoy): Change from `struct` to `class`, so we can maintain the
//  invariants that R_CWs correspond to surface normals and prevent someone
//  from changing `contact_surface` without updating `R_CWs`.

/* A wrapper around DeformableContactSurface that provides additional
 information about the geometries/bodies involved in the contact and proximity
 properties of the contacts. */
template <typename T>
struct DeformableRigidContactPair {
  DRAKE_DEFAULT_COPY_AND_MOVE_AND_ASSIGN(DeformableRigidContactPair)
  DeformableRigidContactPair(DeformableContactSurface<T> contact_surface_in,
                             geometry::GeometryId rigid_id_in,
                             int deformable_id_in, const T& k,
                             const T& d, const T& mu)
      : contact_surface(std::move(contact_surface_in)),
        rigid_id(rigid_id_in),
        deformable_id(deformable_id_in),
        stiffness(k),
        dissipation(d),
        friction(mu),
        R_CWs(contact_surface.num_polygons()) {
    for (int ic = 0; ic < contact_surface.num_polygons(); ++ic) {
      const Vector3<T>& nhat_W = contact_surface.polygon_data(ic).unit_normal;
      constexpr int axis = 2;
      auto R_WC = math::RotationMatrix<T>::MakeFromOneUnitVector(nhat_W, axis);
      R_CWs[ic] = R_WC.transpose();
    }
  }

  /* Returns the number of contact points between the rigid and deformable body.
   We have one contact point at the centroid of each contact polygon in the
   contact surface. */
  int num_contact_points() const { return contact_surface.num_polygons(); }

  DeformableContactSurface<T> contact_surface;
  geometry::GeometryId rigid_id;  // The id of the rigid geometry in contact.
  int deformable_id;  // The id of deformable body in contact.
  T stiffness;    // Combined stiffness at the contact point.
  T dissipation;  // Combined dissipation at the contact point.
  T friction;     // Combined friction at the contact point.
  /* The rotation matrix mapping world frame quantities into contact frames at
   each contact point. The i-th contact point has its own contact frame Cᵢ,
   where R_CᵢW = R_CWs[i]. The basis vector Cᵢz is along the contact surface's
   normal at that contact point. */
  std::vector<math::RotationMatrix<T>> R_CWs;
};

}  // namespace internal
}  // namespace geometry
}  // namespace drake
