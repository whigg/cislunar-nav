% haloOrbitFinder.m
% Author: Mark Hartigan
% Date  : October 5, 2023
% Description:
%    Script to find halo orbits based on Dr. Kathleen Howell's 1983 paper
%    (https://adsabs.harvard.edu/full/1984CeMec..32...53H (Howell, 1983)).

%% reset
clc, clear, close all;

%% init
G = 6.67430e-20;    % km^3/(kg*s^2), gravitational constant
m1 = 5.9724e24;     % kg, mass of Earth
m2 = .07346e24;     % kg, mass of Moon
r_12 = 401000;      % km, distance between earth and moon

w = sqrt(G*(m1 + m2)/r_12^3);   % rad/s, rotation rate of system
mu_1 = G*m1;        % km^3/s^2 (earth)
mu_2 = G*m2;        % km^3/s^2 (moon)
mu = m2 / (m1 + m2);
% mu = 0.04;

r1_norm = 6378 / r_12;  % normalized earth radius
r2_norm = 1737 / r_12;  % normalized moon radius
peri = (75139 + 1737) / r_12;   % normalized perilune
x0 = 1 - cosd(45) * peri;       % initial x position
z0 = -sind(45) * peri;           % initial z position

%% integration and optimization settings
phi0 = eye(6);      % initial STM
X0 = [x0 0 z0 0 0.2 0 reshape(phi0, 1, 36)]';
% X0 = [1.057222 0 0.300720 0 -0.238026 0 reshape(phi0, 1, 36)]';
opts = odeset('Events', @periodevent, 'RelTol', 1e-13, 'AbsTol', 1e-13);

tol = 1e-8;         % convergence tolerance
nmax = 10;         % maximum number of iterations

for i=1:nmax
    % integrate for half a period
    [T,X] = ode45(@(t,x) dynamicsSTM_CR3BP(t, x, mu), [0 5], X0, opts);

    if abs(X(end,4)) < tol && abs(X(end,6)) < tol
        break
    elseif i == nmax
        fprintf("Solution failed to converge in %d iterations.\n", i);
        break;
    end

    phi = reshape(X(end, 7:end), 6, 6)'; % final STM
    dXf = dynamicsSTM_CR3BP(0, X(end,:)', mu);       % final state derivative

    % Iteration scheme from (Howell, 1983). See
    % https://adsabs.harvard.edu/full/1984CeMec..32...53H
    W = [phi(4,3) phi(4,5); phi(6,3) phi(6,5)] - 1/dXf(2) * [dXf(4); dXf(6)] * [phi(2,3) phi(2,5)];
    dZYdot = W \ [-dXf(1); -dXf(3)];
    X0(3) = X0(3) + dZYdot(1);
    X0(5) = X0(5) + dZYdot(2);
end


%% plot
% integrate for a full period for visualization
opts2 = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[T,X] = ode45(@(t,x) dynamicsSTM_CR3BP(t, x, mu), [0 T(end)*2], X0, opts2);

% [T,X] = ode45(@(t,x) halodyn(x, 0.04), [0 2.6], [1.220839 0 0.200987 0 -0.310434 0 reshape(phi0, 1, 36)], opts2);

figure();
% trajectory
plot3(X(:,1), X(:,2), X(:,3), 'r', 'LineWidth', 1);
hold on; grid on; axis equal;
X0 = [1.147204488778055 0 -0.151612989935177 0 -0.219996959246244 0 reshape(phi0, 1, 36)]';
[T,X] = ode45(@(t,x) dynamicsSTM_CR3BP(t, x, mu), [0 T(end)*2], X0, opts2);
plot3(X(:,1), X(:,2), X(:,3), 'r', 'LineWidth', 1);
xlabel('x'); ylabel('y'); zlabel('z');
title("L1 and L2 Halo Orbits (CR3BP)");

% planets
% [xx1, yy1, zz1] = ellipsoid(0, 0, 0, r1_norm, r1_norm, r1_norm);    % earth
% earth = surf(xx1, yy1, zz1);
[xx2, yy2, zz2] = ellipsoid(1, 0, 0, r2_norm, r2_norm, r2_norm);    % moon
moon = surf(xx2, yy2, zz2);

%% final values
L1 = 3.275011930135282e+05;
tu = (2*pi / w) / T(end);   % seconds per unit of time
xf = (X0(1)-1) * r_12; zf = X0(3) * r_12; dyf = X0(5) * r_12 / tu;

%% functions
function [value, isterminal, direction] = periodevent(~, x)
%PERIODEVENT Stops ODE45 at crossing of X-Z plane (time = T/2)

value = x(2);
isterminal = 1;
direction = 0;
end