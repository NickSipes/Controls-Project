clc
close all
clear all

% Define Globals
global torque;
global ee_pos;
global eepos_error;
global IC;
global tf;
global time;
time = [];
torque=[];
ee_pos = [];
eepos_error = [];
% Parameters
tf = 10.0;

% Get the M C and G matricies
[M,C,G] = getDynamicModel();

% Get Initial Joint Vairables
% q1 q2 q1d q2d
X0 = [deg2rad(-30), 0, 0, 0];
IC = X0;

% ODE45
[T,X] = ode45(@(t,x)RobotModel(t,x,M,C,G),[0 tf],X0);

showRobot(T,X,1)
showStateVariables(T,X)

function title = buildTitleString()
    global setpoint;
    global IC;
    global tf;
    
    setpoint_str = 'Setpoint =';
    IC_str = 'IC =';
    for i = 1:length(setpoint)
        if i <= 2
            setpoint_str = strcat(" ",setpoint_str, num2str(rad2deg(setpoint(i))), ' Deg | '); 
            IC_str = strcat(' ',IC_str, num2str(rad2deg(IC(i))), ' Deg | '); 
        elseif i == 3
            setpoint_str = strcat(' ',setpoint_str, num2str(rad2deg(setpoint(i))), ' Deg/s | '); 
            IC_str = strcat(' ',IC_str, num2str(rad2deg(IC(i))), ' Deg/s | '); 
        elseif i == 4
            setpoint_str = strcat(' ',setpoint_str, num2str(rad2deg(setpoint(i))), ' Deg/s'); 
            IC_str = strcat(' ',IC_str, num2str(rad2deg(IC(i))), ' Deg/s'); 
        end
    end

    title = {['System State Variables for ', num2str(tf), ' seconds'],IC_str, setpoint_str}; 


end

function showStateVariables(T,X)
    global torque
    global time
    
    title_fontsize = 24;
    subplot_fontsize = 18;
    axis_fontsize = 14;
    
    x = T;
    time_label = 'Time [s]';
    linewidth = 4;

    % Figure for the subplots
    sgtitle(buildTitleString(), 'FontSize', title_fontsize); 

    % Plot the robot body
    hold on
    subplot(2,2,1)
    plot(x,rad2deg(sin(X(:,1))), 'LineWidth', linewidth)
    title('Position q1 (body angle)', 'FontSize', subplot_fontsize)
    xlabel(time_label, 'FontSize', axis_fontsize);
    ylabel('Degrees', 'FontSize', axis_fontsize)

    subplot(2,2,2)
    plot(x,rad2deg(sin(X(:,2))), 'LineWidth', linewidth)
    ylim([-360,360]);
    title('Position q2 (reaction wheel angle)', 'FontSize', subplot_fontsize)
    xlabel(time_label, 'FontSize', axis_fontsize);
    ylabel('Degrees', 'FontSize', axis_fontsize)

    subplot(2,2,3)
    plot(x,rad2deg(sin(X(:,3))), 'LineWidth', linewidth)
    title('Velocity q1 (change in body angle)', 'FontSize', subplot_fontsize)
    xlabel(time_label, 'FontSize', axis_fontsize);
    ylabel('Degrees per second', 'FontSize', axis_fontsize)

    subplot(2,2,4)
    plot(x,rad2deg(sin(X(:,4))), 'LineWidth', linewidth)
    title('Velocity q2 (change in reaction wheel angle)', 'FontSize', subplot_fontsize)
    ylim([-360,360]);
    xlabel(time_label, 'FontSize', axis_fontsize);
    ylabel('Degrees per second', 'FontSize', axis_fontsize)
    hold off
    
    % Plot the Joint Torques
    figure();
    plot(time, torque(2,:),'LineWidth', linewidth);
    title('Reaction Wheel Joint Torque', 'FontSize', subplot_fontsize)
    xlabel(time_label, 'FontSize', axis_fontsize);
    ylabel('Torque [Nm]', 'FontSize', axis_fontsize)
end

function showRobot(T,X,record_movie)
reaction_wheel_radius = 0.15;

    for i = 1:length(T)
        % Create the figure and show as image
        figure(1);
        
        % Get the ee position base
        ee = forKin(X(i,1),X(i,2));
        
        % Plot the robot body
        hold on
        title(['Robot at time:', num2str(T(i))]) 
        plot([0, ee(1)], [0, ee(2)]);
        % Plot the reaction wheel position
        plot(ee(1),ee(2),'o', 'MarkerSize', 20);
        plot([ee(1), ee(1) - reaction_wheel_radius*sin(X(i,2))], [ee(2), ee(2) - reaction_wheel_radius*cos(X(i,2))], 'r');
        xlim([-1, 1]);
        ylim([-1, 1]);
        hold off
        
        F(i) = getframe(gcf);
        clf('reset');
    end
    
    if record_movie == 1
        % create the video writer with 1 fps
      writerObj = VideoWriter('Robot_Motion.avi');
      writerObj.FrameRate = 10;
      % set the seconds per image
        % open the video writer
        open(writerObj);
        % write the frames to the video
        for i=1:length(F)-1
            % convert the image to a frame
            frame = F(i) ;    
            writeVideo(writerObj, frame);
        end
        % close the writer object
        close(writerObj)    
    end
    
end

function dx = RobotModel(t,x,M,C,G)

global time;
global torque;
global ee_pos;
global eepos_error;
global setpoint;

time= [time t];

% Desired Joint Values at endpoint
q1_des = 0; %This is the edge of controllability -> deg2rad(17.18);
q2_des = 0; % shouldn't be relevant

% Desired Set-Point Position
theta_d = [q1_des; q2_des];

% Desired Velocity
dtheta_d = [0;0];

% Desired setpoint in state vars
setpoint = [theta_d; dtheta_d];

% Desired Acceleration
ddtheta_d = [0;0];

% Current State
current_thetas = x(1:2,1); 
current_vel = x(3:4,1);

% Desired end effector position
ee_d = forKin(q1_des, q2_des);

% Current position
ee = forKin(current_thetas(1), current_thetas(2));

% Position Error
ee(2) = rad2deg(sin(x(4,1)))/360;
ee_error = ([ee_d(1)-ee(1); ee_d(2)-ee(2)]);

% Record Data
ee_pos = [ee_pos ee];
eepos_error = [eepos_error ee_error];

% Calculate Dynamics
[M,C,G] = dynamicModelResults(x,M,C,G);
invM = inv(M);
invMc = invM*C;
invMg = invM*G;

% Control Law
% PD plus Feed Forward (point to point)
xd = [q1_des; q2_des; 0; 0];
[Md, Cd, Gd] = dynamicModelResults(xd,M,C,G);
tau = PDplusFeedforward(theta_d, dtheta_d, ddtheta_d, current_thetas,...
                        current_vel, Md, Cd);
% Record Data
torque = [torque tau];

% Incremental Change
dx = zeros(4,1);
dx(1) = x(3); % dq1 
dx(2) = x(4); % dq2
dx(3:4) = -invMc .* x(3:4) + invM * tau - invMg; % No gravity compensation
disp(t)
end

% Control law for our robot
function tau = PDplusFeedforward(theta_d, dtheta_d, ddtheta_d, theta, dtheta, Mmatd, Cmatd)
     % Gain Matricies
     Kp = [100 0;
           0 0];
     Kv = [50 0;
           0 0];
     
     % position error
     e = theta_d - theta;
     
     % velocity error
     de = dtheta_d - dtheta;
     
     % Controller output torque
     tau = double((Kp*e + Kv*de) + Cmatd.*dtheta_d + Mmatd*ddtheta_d);
     
     % Take desired q1 torque and calculate what q2 torque is needed to get it
     
     % Output new torque
     % These need to be actually calculated
    tau(2) = -tau(1);
    tau(1) = 0;
end

% Updated
function [q1, q2] = invKin(x, y)

l_originToBodyCM_param  = 0.4;  % m
l_originToWheelCM_param = 0.4;  % m

full_body_length = l_originToBodyCM_param + l_originToWheelCM_param;

q1 = asin(y/full_body_length);
q2 = 0; % Not sure how to do this part
end

% Updated
function ee = forKin(q1, q2)

l_originToBodyCM_param  = 0.4;  % m
l_originToWheelCM_param = 0.4;  % m

full_body_length = l_originToBodyCM_param + l_originToWheelCM_param;

x = full_body_length * sin(q1);
y = full_body_length * cos(q1);

ee = [x;y];
end

function [M, C, G] = dynamicModelResults(x,M,C,G)
    % unpack inputs
    syms q1 q2 dq1 dq2
    joint_vars = [q1 q2 dq1 dq2];
    
    % Sub in the joint values from the current state
%     for i = 1:length(joint_vars)
%         M = subs(M, joint_vars(i), x(i));
%         C = subs(C, joint_vars(i), x(i));
%         G = subs(G, joint_vars(i), x(i));
%     end  

        C = subs(C, joint_vars(1), x(1));
        G = subs(G, joint_vars(1), x(1));
end

function [M, C, G] = getDynamicModel()
    % Get the dynamic model
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
    
    [M, C, G] = getRobotParams(M,C,G);  
end

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

function [M,C,G] = getRobotParams(M,C,G)

    syms m_body m_wheel g...
    l_originToBodyCM l_originToWheelCM I_Body I_WheelRotation
    
    % Parameters with numerical values 
    % (I made these physical values up[N.Sipes-4/18/2020])
    m_body_param            = 10;   % kg
    m_wheel_param           = 1;    % kg
    g_param                 = 9.8; % m/s2
    r1_param                = 0.4;  % m
    r2_param                = 0.5; % m
    l_originToBodyCM_param  = 0.4;  % m
    l_originToWheelCM_param = 0.4;  % m
    I_Body_param            = m_body_param * l_originToBodyCM_param^2 + m_wheel_param * l_originToWheelCM_param^2 ;
    I_WheelRotation_param   =  0.5 * m_wheel_param * (r1_param^2 + r2_param^2);
    
    params = [m_body_param m_wheel_param g_param...
    l_originToBodyCM_param l_originToWheelCM_param I_Body_param I_WheelRotation_param];
    
    % Symbolic variables to be replaced
    symbolic_vars = [m_body m_wheel g...
    l_originToBodyCM l_originToWheelCM I_Body I_WheelRotation];
    
    for i = 1:length(params)
        M = subs(M, symbolic_vars(i), params(i));
        C = subs(C, symbolic_vars(i), params(i));
        G = subs(G, symbolic_vars(i), params(i));
    end

end

