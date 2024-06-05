% PlotLSPAngle.m
% Author: Mark Hartigan
% Date  : May 14, 2024
% Description:
%    Plots the angle of the MOON ME frame z-axis w.r.t. the MOON OP frame
%    z-axis. Since the MOON OP x-axis is defined as the cross product of
%    z_ME x z_OP, we can perform a rotation about the x-axis to align the
%    two. Thus, lunar frozen orbits can be represented simply in a frame
%    that accurately captures their visibility to the lunar south pole
%    while still allowing general representations using orbital elements.
%    The purpose of this script is to validate that this angle is
%    relatively constant over long periods of time.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
% SIMULATION START DATE AND TIME
START = '2000 Jan 1 00:00:00';
% DURATION OF SIMULATION
DAYS = 30*365;

%% eval
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
t0 = cspice_str2et(START);
ts = t0:86400:t0 + 86400*DAYS;
tp = (ts - t0) / (86400*365);
ang = zeros(size(ts));
for i=1:length(ts)
    z_OP = cspice_pxform('MOON_ME', 'MOON_OP', ts(i)) * [0;0;1];
    ang(i) = atan2(norm(z_OP(1:2)),z_OP(3)) * 180/pi;
end

%% plot
plotformat("IEEE", 0.75, "scaling", 2);
figure();
plot(tp,ang);
grid on;
xlabel("Time (yrs)");
ylabel("Angle (deg)");
title("Angle between z_{OP} and z_{ME}, starting on J2000");
