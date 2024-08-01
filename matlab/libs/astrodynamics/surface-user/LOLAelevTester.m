% LOLAelevTester.m
% Author: Mark Hartigan
% Date  : July 17, 2024
% Description:
%    Plot the elevation map provided by LOLA's laser altimeter and GSFC.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% import data
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
START = '2024 May 1 00:00:00';
t0 = cspice_str2et(START);

data = importdata("data/LDEM_80S_80MPP_ADJ.TIF");
dmin = min(data,[],"all");
dmax = max(data,[],"all");
range = 304;            % km, stereographic range
mpp = 80;               % m, meters per pixel
kmpp = mpp * 1e-3;
range = range - kmpp/2;
px2km = -range:kmpp:range;

%% plot
figure();
imagesc([-range range], [range -range], data ./ 1000);
set(gca, 'YDir', 'normal');
xlabel("y_{ME} (km)");ylabel("x_{ME} (km)");
c = colorbar();
c.Label.String = "Elevation from mean radius (km)";
load('roma.mat');
colormap(flipud(roma));

%%
[pts, reftraj, exptraj] = surfacetrajgen(t0,[-pi/2.01 0 1], 4, 1.6, 30*60);
tf = reftraj.breaks(end); 
ts = linspace(t0, tf, 1000);
xs_ref = ppval(reftraj, ts);
xs_exp = ppval(exptraj, ts);

figure();
% R = 1737.4;
% [xx, yy, zz] = ellipsoid(0, 0, 0, R, R, R);
% globe = surf(xx, yy, zz);
% set(globe, 'FaceColor', 'none', 'FaceAlpha', 0, 'EdgeColor', [0 0 0]);
plot3(xs_ref(1,:), xs_ref(2,:), xs_ref(3,:), '--');
hold on;
plot3(xs_exp(1,:), xs_exp(2,:), xs_exp(3,:));
scatter3(pts(1,:), pts(2,:), pts(3,:), "o");
hold off; grid on; axis equal;
xlabel("x_{ME} (km)"); ylabel("y_{ME} (km)"); zlabel("z_{ME} (km)");
legend(["Straight-line path", "True path", "Waypoints"]);
title("Route for surface user near the lunar south pole");

plotformat("IEEE", 0.75, "scaling", 2);

figure();
subplot(2,1,1);
plot((ts - t0)/3600, sqrt(sum(xs_exp(4:6,:).^2, 1))*3600);
ylabel("Speed (km/hr)");
grid on;
title("User velocity and acceleration");

subplot(2,1,2);
plot((ts - t0)/3600, sqrt(sum(xs_exp(7:9,:).^2, 1))*1000);
ylabel("Acceleration (m/s^2)");
xlabel("Time (hrs)");
grid on;
