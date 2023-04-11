%% reset
clc, clear, close all;

%% data
true = importdata("true.csv");
est  = importdata("est.csv");
% true = true ./ 1000; est = est ./ 1000; % convert to km
err = true - est;
err = sqrt(sum(err(1:3,:).^2, 1));

%% plot
n = length(true);
% n = 600

figure(1);
scatter3(true(1,1), true(2,1), true(3,1), 'go'); hold on;
scatter3(est(1,1), est(2,1), est(3,1), 'ro');
%scatter3(est(1,end), est(2,end), est(3,end), 'rx');
plot3(est(1,1:n), est(2,1:n), est(3,1:n), 'b');
plot3(true(1,1:n), true(2,1:n), true(3,1:n), 'LineWidth', 2);
grid on; axis equal;

rad = 1737.400;             % km, lunar radius
% moon
[I, map] = imread("lroc_color_poles_1k.jpg");
[xx, yy, zz] = ellipsoid(0, 0, 0, rad, rad, rad);
globe = surf(xx, yy, -zz);
set(globe, 'FaceColor', 'texturemap', 'CData', I, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');

xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");

figure(2);
plot(1:length(err), err);
grid on;
xlabel("count"); ylabel("RMS error (km)");