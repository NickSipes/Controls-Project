%% Inverted Pendulum
% RBE502 Group Project
% % % % % % % % % % % % % % % % %
clc
close all
clear all

%% Kinetic and Potential Energies
% q1 -> body angle
% q2 -> wheel angle
syms m_wheel l_body q1 q2 dq1 dq2 ddq1 ddq2 g r1 r2

I_wheel_rotation = 0.5 * m_wheel * (r1^2 + r2^2);
I_body = m_wheel * l_body^2;

P = m_wheel * g * l_body * sind(q1);
K = 0.5 * I_wheel_rotation * dq2^2 + 0.5*I_body*dq1^2;

%% Lagrange Equation
L = simplify(K-P);

%% Derive the dynamics
% Solve for TorqueDynam
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
syms u1

% 2 DOF robot = 4 states. 
SS = M \ (simplify([u1] - C));
dx = [dq1;dq2;SS(1);SS(2)];
