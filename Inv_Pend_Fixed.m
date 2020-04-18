%% Inverted Pendulum
% RBE502 Group Project
% % % % % % % % % % % % % % % % %
clc
close all
clear all

%% Kinetic and Potential Energies
% q1 -> body angle
% q2 -> wheel angle
syms m_body m_wheel q1 q2 dq1 dq2 ddq1 ddq2 g...
    l_originToBodyCM l_originToWheelCM I_Body I_WheelRotation

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

% state-space format.
syms u1

% 2 DOF robot = 4 states. 
SS = M \ (simplify([u1] - C));
dx = [dq1;SS(1);dq2;SS(2)];

dx = getRobotParams(dx);


%% Helper Functions

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

% vizualize the state trajectory
function seeTrajectory(dx)

    

end

% Replace symbolic values with real numbers
function dx = getRobotParams(dx)

    syms m_body m_wheel g...
    l_originToBodyCM l_originToWheelCM I_Body I_WheelRotation
    
    % Parameters with numerical values 
    % (I made these physical values up[N.Sipes-4/18/2020])
    m_body_param            = 10;   % kg
    m_wheel_param           = 1;    % kg
    g_param                 = -9.8; % m/s2
    r1_param                = 0.1;  % m
    r2_param                = 0.15; % m
    l_originToBodyCM_param  = 0.4;  % m
    l_originToWheelCM_param = 0.4;  % m
    I_Body_param            = m_wheel_param * l_originToBodyCM_param^2;
    I_WheelRotation_param   =  0.5 * m_wheel_param * (r1_param^2 + r2_param^2);
    
    params = [m_body_param m_wheel_param g_param...
    l_originToBodyCM_param l_originToWheelCM_param I_Body_param I_WheelRotation_param];
    
    % Symbolic variables to be replaced
    symbolic_vars = [m_body m_wheel g...
    l_originToBodyCM l_originToWheelCM I_Body I_WheelRotation];
    
    for i = 1:length(params)
        dx = subs(dx, symbolic_vars(i), params(i));
    end

end
