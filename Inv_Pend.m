%% Inverted Pendulum
% RBE502 Group Project
% % % % % % % % % % % % % % % % %
clc
close all
clear all

syms q1 q2 dq1 dq2 ddq1 ddq2 g...
     m l
pi = sym(pi);

sv = {m, l, dq1, dq2, g};
q = [q1 q2].';
dq = [dq1 dq2].';
ddq = [ddq1 ddq2].';

theta = q1;
d = l;
a = q2;
alpha = 0;
[n_dof,rand] = size(q);

%% DH Parameters
%Define transformation matrices
T01 = DH(theta,d,a,alpha);

%% Velocity Kinematics
% 6x3 Jacobian which describes the forward velocity Kinematics
% vx' vy' vz' values
A = [diff(T01(1:3,4), q1), diff(T01(1:3,4), q2)];
% Angular Velocities - wx', wy', wz'. Only rotation about Z axis
B = [0 0; 0 0; 1 1]; % set rotation about z = 1 for all joints.
J = [A;B];

%% Kinetic and Potential Energies
% solve for x_dots
dx1 = [diff(T01(1:3,4),q1)]*dq1;
dx2 = [dx1 + diff(T01(1:3,4),q2)*(dq1+dq2)];

K1 = simplify(.5*m*(transpose(dx1)*dx1));
K2 = simplify(.5*m*(transpose(dx2)*dx2));

% Kinetic
K = K1+K2;

% PE = m*g*h   h = z
% P1 = m1*g*y1
y1 = T01(3,4);

%Set all masses = 0 
P1 = m*g*y1;
% Potential Energy
P = P1;

%% Lagrange Equation
L = simplify(K-P);

%% Derive the dynamics
% Solve for Torque
T1 = simplify(diff(L,dq1)*ddq1 - diff(L,q1)*dq1);
T2 = simplify(diff(L,dq2)*ddq2 - diff(L,q2)*dq2);
Tau = [T1;T2];

% Solve for M
M11 = simplify(expand(T1 - subs(T1,ddq1,0))/ddq1);
M12 = simplify(expand(T1 - subs(T1,ddq2,0))/ddq2);
M21 = simplify(expand(T2 - subs(T2,ddq1,0))/ddq1);
M22 = simplify(expand(T2 - subs(T2,ddq2,0))/ddq2);
M = [M11 M12;
M21 M22];

% solve for G
G1 = subs(T1,{ddq1,ddq2,dq1,dq2},{0,0,0,0});
G2 = subs(T2,{ddq1,ddq2,dq1,dq2},{0,0,0,0});
G = [G1;G2];

% solve for C
C1 = simplify(expand(T1 - M(1,:)*[ddq1 ddq2].' + G(1)));
C2 = simplify(expand(T2 - M(2,:)*[ddq1 ddq2].' + G(2)));
C = [C1;C2];

%     subs(M,[mt,dq1,dq2,dq3,g],[.5,q1_d,q2_d,q3_d,g]);
%     subs(C,[mt,dq3,g],[mt,q3_d,g]);

% state-space format.
syms u1 u2

% 2 DOF robot = 4 states. 
SS = M \ (simplify([u1;u2] - C));
dx = [dq1;dq2;SS(1);SS(2)];

function [Transform] = DH(theta,d,a,alpha)

Transform = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta);
    sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
    0 sin(alpha) cos(alpha) d;
    0 0 0 1];
end
