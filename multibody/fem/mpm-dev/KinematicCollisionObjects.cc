#include "KinematicCollisionObjects.h"

namespace drake {
namespace multibody {
namespace mpm {

void KinematicCollisionObjects::AddCollisionObject(
                        std::unique_ptr<AnalyticLevelSet> level_set,
                        const math::RigidTransform<double>& pose,
                        const math::SpatialVelocity<double>& spatial_velocity,
                        double friction_coeff) {
    collision_objects_.emplace_back({level_set,
                                     {pose, spatial_velocity},
                                     friction_coeff});
}

void KinematicCollisionObjects::AdvanceOneTimeStep() {
    for (const auto& obj : collision_objects_) {
        obj->AdvanceOneTimeStep();
    }
}

void KinematicCollisionObjects::ApplyBoundaryConditions(
                                            const Vector3<double>& position,
                                            Vector3<double>* velocity) const {
    for (const auto& obj : collision_objects_) {
        obj->ApplyBoundaryConditions(position, velocity);
    }
}

KinematicCollisionObjects::CollisionObject(
                                    std::unique_ptr<AnalyticLevelSet> level_set,
                                    CollisionObjectState initial_state,
                                    double friction_coeff):
                                            level_set_(level_set),
                                            state_(initial_state),
                                            friction_coeff_(friction_coeff) {}

KinematicCollisionObjects::CollisionObject::AdvanceOneTimeStep(dt) {
    RigidTransform transform_new = {dt*state_.spatial_velocity.translational(),
                                    };
}

KinematicCollisionObjects::CollisionObject::ApplyBoundaryConditions(
                                        const Vector3<double>& position,
                                        const Vector3<double>* velocity) const {
    // Get the translational component of the relative spatial velocity at point
    // p_WQ (the grid point) between the grid point Q and the collision object R by
    // subtracting the translational component of the spatial velocity of a point
    // (Rq) coincident with p_WQ on the collision object R from the translational
    // velocity of the grid point, v_WQ (`velocity` on input).
    const Vector3<double> p_WQ = position;
    const Vector3<double> v_WQ = *velocity;
    
    // The pose and spatial velocity of the collision object in world frame.
    const RigidTransformd& X_WR = get_world_pose();
    const RigidTransformd& V_WR = get_world_spatial_velocity();

    // Compute the spatial velocity at Rq for the collision object.
    const Vector3<double> p_RoRq_W = p_WQ - X_WR.translation();
    const SpatialVelocity<double> V_WRq = V_WR.Shift(p_RoRq_W);
    
    // Compute the relative (translational) velocity of grid node
    // Q relative to Frame Rq, expressed in the world frame.
    Vector3<double> v_RqQ_W = v_WQ - V_WRq.translational(); 
    
    // Express the relative velocity in the collision object's frame R.
    Vector3<double> v_RqQ_R = X_WR.inverse() * v_RqQ_W;
    
    // Express the relative position in the collision object's frame R.
    const Vector3<double> p_RoRq_R = X_WR.inverse() * p_RoRq_W;
     
    if (!get_level_set().InInterior(p_RoRq_W)) return;
    UpdateVelocityCoulumbFriction(get_level_set().Normal(p_RoRq_R), get_friciton_coeff(), &v_RqQ_R);
    
    // Transform the velocity back to the world frame.
    v_RqQ_W = X_WR * v_RqQ_R;
    v_WQ = v_RqQ_W + V_WRq.translational();
    *velocity = v_WQ;
}

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
