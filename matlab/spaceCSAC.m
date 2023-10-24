% spaceCSAC.m
% Author: Mark Hartigan
% Date  : October 13, 2023
% Description:
%    Plots statistics about the Microsemi Space CSAC to be used in the Khon
%    satellites.

%% reset
clc, clear, close all;
addpath(genpath(pwd));

%% init
a_hi = 9e-10 / (86400 * 30);
a_lo = 6e-10 / (86400 * 30);
E = @(t, a) 1/2 * a * t.^2;
t = 1:1:3600*24;
c = 299792458;

%% plot
figure();
plot([0 t(end)/3600], [30e-9 30e-9]*c, 'k--', 'LineWidth', 1.5); hold on;
% patch([t flip(t)]/3600, [E(t, a_lo) flip(E(t, a_hi))], 'cyan', 'FaceAlpha', 0.5);
plot(t/3600, E(mod(t,13146), a_hi)*c, 'LineWidth', 1.5);
axis([0 t(end)/3600 0 12]);
grid on;
xlabel("Time (hrs)"); ylabel("Position error (m)");
title("95% RSS Position Error from Clock Model");
legend(["30ns uncertainty threshold", "RSS Error"], 'location', 'best');