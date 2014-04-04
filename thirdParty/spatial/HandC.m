function  [H,C] = HandC(manipulator, q, v, f_ext, grav_accn )

% HandC  Calculate coefficients of equation of motion.
% [H,C]=HandC(model,q,v,f_ext,grav_accn) calculates the coefficients of
% the joint-space equation of motion, tau=H(q)vd+C(d,v,f_ext), where q,
% v and vd are the joint position, velocity and acceleration vectors, H
% is the joint-space inertia matrix, C is the vector of gravity,
% external-force and velocity-product terms, and tau is the joint force
% vector.  Algorithm: recursive Newton-Euler for C, and
% Composite-Rigid-Body for H.  f_ext is a cell array specifying external
% forces acting on the bodies.  If f_ext == {} then there are no external
% forces; otherwise, f_ext{i} is a spatial force vector giving the force
% acting on body i, expressed in body i coordinates.  Empty cells in f_ext
% are interpreted as zero forces.  grav_accn is a 3D vector expressing the
% linear acceleration due to gravity.  The arguments f_ext and grav_accn
% are optional, and default to the values {} and [0,0,-9.81], respectively,
% if omitted.
%
% UPDATED by russt:  f_ext is an empty or (sparse) 6 x model.NB matrix

if nargin < 5
  a_grav = [0;0;0;0;0;-9.81];
else
  a_grav = [0;0;0;grav_accn(1);grav_accn(2);grav_accn(3)];
end

if nargin < 4 || isempty(f_ext)
  f_ext = sparse(6, length(manipulator.body));
end

kinsol = doKinematics(manipulator, q, true, false, v);
inertias_world = computeInertiasInWorld(manipulator, kinsol);
H = computeMassMatrix(manipulator, kinsol, inertias_world);
C = computeBiasTerm(manipulator, kinsol, inertias_world, f_ext, a_grav);

end

function ret = computeInertiasInWorld(manipulator, kinsol) % TODO: consider moving to kinsol; is also useful for CMM
NB = length(manipulator.body);
ret = cell(NB, 1);
ret{1} = zeros(6, 6);
for i = 2 : NB
  body = manipulator.body(i);
  bodyToWorld = kinsol.T{i};
  worldToBody = homogTransInv(bodyToWorld);
  AdWorldToBody = transformAdjoint(worldToBody);
  ret{i} = AdWorldToBody' * body.I * AdWorldToBody;
end
end

function H = computeMassMatrix(manipulator, kinsol, inertias_world)
% world frame implementation
NB = length(manipulator.body);
nv = length(kinsol.v);
H = zeros(nv, nv) * kinsol.q(1); % minor adjustment to make TaylorVar work better.
composite_inertias = computeCompositeRigidBodyInertias(manipulator, inertias_world);

for i = 2 : NB
  Ic = composite_inertias{i};
  Si = kinsol.J{i};
  i_indices = manipulator.body(i).velocity_num;
  F = Ic * Si;
  H(i_indices, i_indices) = Si' * F;
  
  j = i;
  while j ~= 2
    j = manipulator.body(j).parent;
    body_j = manipulator.body(j);
    j_indices = body_j.velocity_num;
    Sj = kinsol.J{j};
    Hji = Sj' * F;
    H(j_indices, i_indices) = Hji;
    H(i_indices, j_indices) = Hji';
  end
end
end

function ret = computeCompositeRigidBodyInertias(manipulator, inertias_world)
% computes composite rigid body inertias expressed in world frame
NB = length(inertias_world);
ret = inertias_world;
for i = NB : -1 : 2
  body = manipulator.body(i);
  ret{body.parent} = ret{body.parent} + ret{i};
end
end

function C = computeBiasTerm(manipulator, kinsol, inertias_world, external_forces, gravitational_accel)
nv = length(kinsol.v);
root_accel = -gravitational_accel; % as if we're standing in an elevator that's accelerating upwards
JdotV = kinsol.JdotV;

NB = length(manipulator.body);
net_wrenches = cell(NB, 1);
net_wrenches{1} = zeros(6, 1);
for i = 2 : NB
  twist = kinsol.twists{i};
  spatial_accel = root_accel + JdotV{i};
  external_wrench = external_forces(:, i); % TODO: check if external_forces are expressed in body frame or world frame
  I = inertias_world{i};
  
  net_wrenches{i} = I * spatial_accel - twistAdjoint(twist)' * I * twist - external_wrench;
end

C = zeros(nv, 1);
for i = NB : -1 : 2
  body = manipulator.body(i);
  joint_wrench = net_wrenches{i};
  tau = kinsol.J{i}' * joint_wrench;
  C(body.velocity_num) = tau;
  net_wrenches{body.parent} = net_wrenches{body.parent} + joint_wrench;
end

end