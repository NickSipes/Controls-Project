clear all
%% Simplified 2D robot model

% Robot Fixed Parameters
l1 = 0.5; % meters
l2 = 0.5; % meters

mDriveWheel = 3; % kg
mBatteryPack = 5; % kg
mReactionWheel = 2; % kg

% Moments for rotating wheels
IDriveWheel = 1;

% Moments for points rotation about the robot center of mass (p2)
% Point mass about point is I= mr^2
Ip1 = mDriveWheel*l1^2;
Ip3 = mReactionWheel*l2^2;

% Robot Variables
syms bodyAngle wBodyAngle x1 y1 wDriveWheel wReactionWheel

% Environment constants
g = -9.8;

% Robot CF origin will be at contact point between ground and drive wheel
% Treating p2 as the center of mass for simplicty

% Represent the robot via 3 points masses 
% p1 -> drive wheels
% p2 -> battery pack
% p3 -> reaction wheels

p1 = [x1, y1];
p2 = [x1 + l1*cosd(bodyAngle), y1 + l1*sind(bodyAngle)];
p3 = [x1 + (l1+l2)*cosd(bodyAngle), y1 + (l1+l2)*sind(bodyAngle)];

v1 = calculateVelocity(p1);
v2 = calculateVelocity(p2);
v3 = calculateVelocity(p3);

% Potential Energy
U1 = mDriveWheel*g*p1(2); 
U2 = mBatteryPack*g*p2(2);
U3 = mReactionWheel*g*p3(2);

% Kinetic Energy
T1 = 0.5*mDriveWheel*v1^2 + 0.5*IDriveWheel*wDriveWheel + 0.5*Ip1*wBodyAngle; % Linear + Drive Wheel Rotational + p1 rotation about p2 
T2 = 0.5*mBatteryPack*v2^2; % Linear
T3 = 0.5*mReactionWheel*v3^2; % Linear + Rotation of Reaction Wheel + Rotation of p3 about p2

%% Simulation
% x1, y1, BodyAngle, wBodyAngle
inital_conditions = [0, 0, 90];

SeeRobot(p1, p2, p3, inital_conditions);

%% Helper Funcitons

function velocity = calculateVelocity(point)

    velocity = sqrt(point(1)^2 + point(2)^2);

end


function SeeRobot(p1, p2, p3, init_conditions)
syms x1 y1 bodyAngle

    % Sub x1,y1 for 0,0 since our robot is starting at the origin
    p1 = double(subs(p1, [x1, y1], init_conditions(1:2)));
    p2 = double(subs(p2, [x1, y1, bodyAngle], init_conditions));
    p3 = double(subs(p3, [x1, y1, bodyAngle], init_conditions));

    figure;
    hold on;
    plot([p1(1), p2(1), p3(1)],[p1(2), p2(2), p3(2)], '-*');
    hold off

end