%% reset
clc, clear, close all;

%% data
true = importdata("true.csv");
est  = importdata("est.csv");
true = true ./ 1000; est = est ./ 1000; % convert to km

%% plot
figure(1);
scatter3(true(1,1), true(2,1), true(3,1), 'go'); hold on;
scatter3(est(1,1), est(2,1), est(3,1), 'ro');
scatter3(est(1,end), est(2,end), est(3,end), 'rx');
plot3(est(1,:), est(2,:), est(3,:), 'b');
plot3(true(1,:), true(2,:), true(3,:));
grid on; %axis equal;

rad = 1737.400;             % km, lunar radius
[xx, yy, zz] = ellipsoid(0, 0, 0, rad, rad, rad);
%globe = surf(xx, yy, -zz, 'FaceColor', 'w', 'EdgeColor', 0.5*[1 1 1]);

xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");