%% reset
clc, clear, close all;

%% init
w = 2*pi / (27.3217 * 24 * 60 * 60);    % rad/s, rotation rate of moon
r = 1737.4;                             % km, radius of moon
g = 1.625e-3;                           % km/s, gravity of moon
v0 = 0.042e-3;          % km/s, max speed of mars rover for reference
v0 = 0;
% 1-3 km, 4-6 km/s, initial position and velocity
ang = 15;
x0 = [r*sind(ang); 0; -r*cosd(ang); cosd(ang)*v0; w*r*sind(ang); sind(ang)*v0];

n = 10;    % number of nodes
t = (flip(chebichev(n)) + 1)/2 * 12 * 86400;  % s, times to eval
% t = (flip(chebichev(n)) + 1)/2 * 86400;  % s, times to eval
pd = random('Normal', 0, 0, 3, n);   % nominally 1e-8 std

%% STM
phi = @(t) [cos(w*t) -sin(w*t) 0        0         0 0;
            sin(w*t)  cos(w*t) 0        0         0 0;
                   0         0 1        0         0 0;
                   0         0 0 cos(w*t) -sin(w*t) 0;
                   0         0 0 sin(w*t)  cos(w*t) 0;
                   0         0 0        0         0 1];

%% eval
% x = zeros(length(x0), length(t));
% x(:,1) = x0;
% for i=2:length(t)
%     x(:,i) = phi(t(i)) * x0;
% end

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
% [t_nom, x_nom] = ode45( ...
%     @(ti,x) surfDyn(x, polyinterp(ti, x(1:3), t, pd), g, r, [0;0;w]), ...
%     [t(1) t(end)], x0, opts);

[t_nom, x_nom] = ode45( ...
    @(ti,x) surfDyn(x, zeros(3,1), g, r, [0;0;w]), ...
    [0 86400*12], x0, opts);

%% postprocess to fix moon intersections
x = x_nom';
% for i = 1:size(x,2) % if position is below surface, put on surface
%     pos = x(1:3,i);
%     if norm(pos) < r, x(1:3,i) = r * pos / norm(pos); end
% end

%% plot
% plot velocity and acceleration
figure(2);
acc = zeros(1,size(x,2));
vel = acc;
for i=1:size(x,2)
    acc(i) = norm(polyinterp(t_nom(i), x(1:3,i), t, pd)) * 1e3;  % convert to m
    vel(i) = norm(x(4:6,i)) * 1e3;              % convert to m
end
yyaxis left;  plot(t_nom, vel);
yyaxis right; plot(t_nom, acc);
xlabel("Time (s)"); grid on;
yyaxis left;  ylabel("v (m/s)");
yyaxis right; ylabel("a (m/s^2)");

% plot trajectory
figure(1);
s = 100;   % scale factor for velocity vectors

% for i=1:size(x,2)
%     quiver3(0, 0, 0, x(1,i), x(2,i), x(3,i), 'b');
%     if i==1, hold on; end
%     quiver3(x(1,i), x(2,i), x(3,i), x(4,i)*s, x(5,i)*s, x(6,i)*s, 'r');
% end
plot3(x(1,:),x(2,:),x(3,:), 'r', 'LineWidth', 2); hold on;
plot3(x0(1), x0(2), x0(3), 'go', 'LineWidth', 2);
% quiver3(x0(1), x0(2), x0(3), x0(4)*1e6, x0(5)*1e6, x0(6)*1e6);

% moon
[I, map] = imread("lroc_color_poles_1k.jpg");
[xx, yy, zz] = ellipsoid(0, 0, 0, r, r, r);
globe = surf(xx, yy, -zz);
set(globe, 'FaceColor', 'texturemap', 'CData', I, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');

hold off; grid on; axis equal;
xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");
legend(["Start", "Path"]);


%% functions
function [poly] = polyinterp(ti, x, t, pd)
%POLYINTERP interpolating polynomial of random points. Adjusts to be
%perpendicular to vector x.
%   Input:
%    - ti; time to interpolate at
%    - x; position (3,)
%    - t; list of times pd is evaluated at (,n)
%    - pd; evaluations of random noise (3,n)

%     poly = zeros(3,1);
    poly = [laginterp(t, pd(1,:), ti);
            laginterp(t, pd(2,:), ti);
            laginterp(t, pd(3,:), ti)];
    
    poly = poly - dot(poly, x) / dot(x, x) * x;
end
