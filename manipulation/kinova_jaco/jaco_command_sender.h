#pragma once

#include "drake/common/drake_copyable.h"
#include "drake/common/drake_deprecated.h"
#include "drake/lcmt_jaco_command.hpp"
#include "drake/manipulation/kinova_jaco/jaco_constants.h"
#include "drake/systems/framework/leaf_system.h"

namespace drake {
namespace manipulation {
namespace kinova_jaco {

/// Creates and outputs lcmt_jaco_command messages.
///
/// Note that this system does not actually send the message to an LCM
/// channel. To send the message, the output of this system should be
/// connected to a
/// systems::lcm::LcmPublisherSystem::Make<lcmt_jaco_command>().
///
/// This system has two vector-valued input ports containing the desired
/// position and velocity.  Finger velocities will be translated to the values
/// used by the Kinova SDK from values appropriate for the finger joints in
/// the Jaco description (see jaco_constants.h).
///
/// This system has one abstract-valued output port of type lcmt_jaco_command.
///
/// @system
/// name: JacoCommandSender
/// input_ports:
/// - position
/// - velocity
/// output_ports:
/// - lcmt_jaco_command
/// @endsystem
///
/// @see `lcmt_jaco_command.lcm` for additional documentation.
class JacoCommandSender : public systems::LeafSystem<double> {
 public:
  DRAKE_NO_COPY_NO_MOVE_NO_ASSIGN(JacoCommandSender)

  JacoCommandSender(int num_joints = kJacoDefaultArmNumJoints,
                    int num_fingers = kJacoDefaultArmNumFingers);

  DRAKE_DEPRECATED("2022-08-01", "Use the other input ports instead.")
  const systems::InputPort<double>& get_input_port() const {
    return *state_input_;
  }

  const systems::InputPort<double>& get_position_input_port() const {
    return *position_input_;
  }
  const systems::InputPort<double>& get_velocity_input_port() const {
    return *velocity_input_;
  }

 private:
  void CalcOutput(const systems::Context<double>&, lcmt_jaco_command*) const;

  const int num_joints_;
  const int num_fingers_;
  const systems::InputPort<double>* state_input_{};
  const systems::InputPort<double>* position_input_{};
  const systems::InputPort<double>* velocity_input_{};
};

}  // namespace kinova_jaco
}  // namespace manipulation
}  // namespace drake
