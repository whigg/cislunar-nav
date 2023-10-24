function plotLunarOrbit(x, name, moon)
%PLOTLUNARORBIT Plots the trajectory of a lunar orbit, given the 3- or
%6-state output of ODE45 (or similar propagators).
%   Inputs: (dims),[units]
%    - x   ; (nx6),[km, km/s] time-history of satellite state
%    - name; (1x1), [N/A] string, plot title
%    - moon; optional boolean argument, whether or not to have moon
%            (default true)

if nargin < 3, moon = true; end

figure();
% satellite
plot3(x(:,1), x(:,2), x(:,3), 'b', 'LineWidth', 2);
hold on;
scatter3(x(end,1), x(end,2), x(end,3), 'ro', 'LineWidth', 2);

% moon
if moon
    r = 1737;   % km, lunar equatorial radius
    [I, ~] = imread("lroc_color_poles_1k.jpg");
    [xx, yy, zz] = ellipsoid(0, 0, 0, r, r, r);
    globe = surf(xx, yy, -zz);
    set(globe, 'FaceColor', 'texturemap', 'CData', I, 'FaceAlpha', 1, ...
        'EdgeColor', 'none');
end

hold off; axis equal; grid on;
xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");
legend(["Orbit", "Satellite"], "location", "best");
title(name);
end

