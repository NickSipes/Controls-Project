%% Inverted Pendulum
% RBE502 Group Project
% % % % % % % % % % % % % % % % %
clc
close all
clear all

%% Kinetic and Potential Energies
% q1 -> body angle
% q2 -> wheel angle
% syms m_body m_wheel l_body q1 q2 dq1 dq2 ddq1 ddq2 g r1 r2...
%     l_originToBodyCM l_originToWheelCM I_Body I_WheelRotation
% 
% %I_WheelRotation = 0.5 * m_wheel * (r1^2 + r2^2);
% %I_Body = m_wheel * l_originToBodyCM^2;
% 
% P = (m_body*l_originToBodyCM + m_wheel*l_originToWheelCM)*g*cos(q1);
% K = 0.5*(m_body*l_originToBodyCM^2 + m_wheel*l_originToWheelCM^2 + I_Body + I_WheelRotation)*dq1^2 + I_WheelRotation*dq1*dq2 + 0.5*I_WheelRotation*dq2^2;

syms m_body m_wheel l_body q1 q2 dq1 dq2 ddq1 ddq2 g r1 r2...
    l_originToBodyCM l_originToWheelCM I_Body I_WheelRotation

%I_WheelRotation = 0.5 * m_wheel * (r1^2 + r2^2);
%I_Body = m_wheel * l_originToBodyCM^2;

P = (m_body*l_originToBodyCM + m_wheel*l_originToWheelCM)*g*cos(q1);
K = 0.5*(m_body*l_originToBodyCM^2 + m_wheel*l_originToWheelCM^2 + I_Body + I_WheelRotation)*dq1^2 + I_WheelRotation*dq1*dq2 + 0.5*I_WheelRotation*dq2^2;

%% Lagrange Equation
L = simplify(K-P);

%% Derive the dynamics
% Solve for TorqueDynamics
% T1 = simplify(diff(L,dq1)*ddq1 - diff(L,q1)*dq1);
% T2 = simplify(diff(L,dq2)*ddq2 - diff(L,q2)*dq2);
[T1,T2] = TorqueDynamics(L, [q1, q2], [dq1, dq2], [ddq1, ddq2]);

% Lineraize the torque equations
T1 = subs(T1, sin(q1), q1);
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

function [T1, T2] = TorqueDynamics(L, q, dq, ddq)

    % Take parital derivatives of q1
    dL_dq = diff(L, dq(1));
    dL_q = diff(L, q(1));

    % Sub variables that are funcitons of time
    syms temp1_q(t) temp1_dq(t) temp1_ddq(t) temp2_q(t) temp2_dq(t) temp2_ddq(t)
    temp = subs(dL_dq, [q(1), dq(1), ddq(1), q(2), dq(2), ddq(2)],...
        [temp1_q(t) temp1_dq(t) temp1_ddq(t) temp2_q(t) temp2_dq(t) temp2_ddq(t)]);
    

    % Take time derivative
    temp = diff(temp, t);
    
    % sub back to orignial variables
    T1 = subs(temp, [temp1_q(t) temp1_dq(t) temp1_ddq(t) temp2_q(t) temp2_dq(t) temp2_ddq(t) diff(temp1_dq(t), t) diff(temp2_dq(t), t)],...
        [q(1), dq(1), ddq(1), q(2), dq(2), ddq(2), ddq(1), ddq(2) ]) - dL_q;
    
    % Take parital derivatives of q1
    dL_dq = diff(L, dq(2));
    dL_q = diff(L, q(2));

    % Sub variables that are funcitons of time
    syms temp1_q(t) temp1_dq(t) temp1_ddq(t) temp2_q(t) temp2_dq(t) temp2_ddq(t)
    temp = subs(dL_dq, [q(1), dq(1), ddq(1), q(2), dq(2), ddq(2)],...
        [temp1_q(t) temp1_dq(t) temp1_ddq(t) temp2_q(t) temp2_dq(t) temp2_ddq(t)]);
    

    % Take time derivative
    temp = diff(temp, t);
    
    % sub back to orignial variables
    T2 = subs(temp, [temp1_q(t) temp1_dq(t) temp1_ddq(t) temp2_q(t) temp2_dq(t) temp2_ddq(t) diff(temp1_dq(t), t) diff(temp2_dq(t), t)],...
        [q(1), dq(1), ddq(1), q(2), dq(2), ddq(2), ddq(1), ddq(2) ]) - dL_q;

end

function motorEquations()
% R


    voltage = Rm*i + Ke*wm

end
