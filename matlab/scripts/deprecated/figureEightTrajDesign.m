% figureEightTrajDesign.m
% Author: Mark Hartigan
% Date  : October 13, 2023
% Description:
%    Generates a figure-eight ground trajectory for a lunar surface user,
%    taking inspiration from a geosynchronous satellite ground track. The
%    data products generated are in-line with Bong-Jun Yang's request.

%% reset
clc, clear, close all;
addpath(genpath(pwd));

%% constants
P = 86400;      % s, 24 hours
R = 1736.0;     % km, moon polar radius
mu = 4.903e3;   % km^3/s^2, gravitational parameter of moon
w = 2*pi / P;   % rad/s, fake rotation rate of moon

a = (mu * (P / (2*pi))^2)^(1/3);    % km, semimajor axis of GSO
i = 60 * pi/180;                    % rad, inclination
e = 0.0;                            % eccentricity
omega =  0 * pi/180;                % rad, argument of periapsis

[r0, v0] = oe2rv(a, e, i, 0, omega, 0, mu); % position and velocity vectors

%% propagate
x0 = [r0; v0];
opts = odeset('AbsTol', 1e-9, 'RelTol', 1e-9);
[T, X] = ode45(@(~,x) twobp(x, mu), 0:1:P, x0, opts);
% plotLunarOrbit(X, 'Moon-centered Inertial');

%% transform to moon-centered moon-fixed frame
X = X';             % rotate to better orientation
n = size(X, 2);     % number of data points
x_MCMF = zeros(size(X));
for i=1:n
    t = T(i);
    L = rotz(w*t);
    r = norm(X(1:3,i)); % orbital radius
    % ignore velocity for now, scale position onto the lunar surface
    x_MCMF(1:3,i) = L * X(1:3,i) * R / r;
end
% plotLunarOrbit(x_MCMF', 'Moon-centered Moon-fixed');

%% rotate to south pole and scale angular resolution to reasonable size
% x0 was placed on the x-axis deliberately s.t. we can just rotate about
% the y-axis to bring these vectors down to the south pole
x_eight = zeros(size(x_MCMF));

% eight placed on south pole and scaled smaller (proportions between angles
% maintained)
L = roty(-pi/2);                % rotation to south pole
Z = [0; 0; -1];                 % z-axis
sc = 0.022;                   % scale factor, ensures max vel ~ 10 km/h

for i=1:n
    r_ = x_MCMF(1:3,i);         % current position vector
    r_ = L * r_;                % rotate to south pole
    
    r = norm(r_);               % current radius
    a = acos(dot(Z, r_) / r);   % current angle w/ z-axis

    k = cross(r_, Z) / norm(cross(r_, Z)); % rotation axis
    da = a * (1 - sc);          % change in angle desired
    % implement rodrigues' rotation formula
    x_eight(1:3,i) = r_ * cos(da) + cross(k, r_) * sin(da) + k * dot(k, r_) * (1 - cos(da));
end
% plotLunarOrbit(x_eight', 'Moon-centered Moon-fixed');

%% transform coordinates into local NED frame (south pole) and compute v/a
% L would just be the identity matrix, so all we need to do is offset the
% origin to the south pole
x_NED = zeros(9, n);
vmag = zeros(n,1);
amag = zeros(n,1);

for i=1:n
    % translate coordinate system to local NED frame
    x_NED(1:3,i) = x_eight(1:3,i) + [0; 0; R];
    % compute velocity just through forward differences
    if i ~= n
        x_NED(4:6,i) = (x_eight(1:3,i+1) - x_eight(1:3,i)) / (T(i+1)-T(i));
    else
        x_NED(4:6,i) = x_NED(4:6,1);
    end

    vmag(i) = norm(x_NED(4:6,i)) * 60 * 60; % convert from km/s to km/h
end

for i=1:n
    % compute acceleration
    if i ~= n
        x_NED(7:9,i) = (x_NED(4:6,i+1) - x_NED(4:6,i)) / (T(i+1)-T(i));
    else
        x_NED(7:9,i) = x_NED(7:9,1);
    end

    amag(i) = norm(x_NED(7:9,i)) * (60 * 60)^2; % convert from km/s^2 to km/h^2
end

% find the maximum velocity to determine if we need to schange sc
fprintf("Maximum velocity: %.2f km/h\n", max(vmag));
% find the maximum acceleration just for fun
fprintf("Maximum acceleration: %.2f km/h^2\n", max(amag));
plotLunarOrbit(x_NED', 'NED Frame', false);

figure();
plot(T, vmag);
grid on; xlabel("Time (s)"); ylabel("Velocity (km/h)");

figure();
plot(T, amag);
grid on; xlabel("Time (s)"); ylabel("Acceleration (km/h^2)");

%% write data to file
data = [T x_NED'];
tab  = array2table(data, 'VariableNames', ...
    {'Time (s)', 'Position.x (km)', 'Position.y (km)', 'Position.z (km)', ...
    'Velocity.x (km/s)', 'Velocity.y (km/s)', 'Velocity.z (km/s)', ...
    'Acceleration.x (km/s^2)', 'Acceleration.y (km/s^2)', ...
    'Acceleration.z (km/s^2)'});
writetable(tab, 'data/figureEightTrajData.txt');

