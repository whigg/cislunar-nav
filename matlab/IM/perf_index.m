function [P, x2, x3, x4] = perf_index(f, t, x1, ephfun, oe, moon, debug)
%PERF_INDEX Performance index for lunar constellations evaluating PDOP at
%the south pole. To be used in a broader optimizer for true anomalies of
%orbiting satellites.
% Input: (see ta_optimization.m for how structs are defined)
%  - f; state vector of changes to TA for s/c 2-4 (3x1) [rad]
%  - t; time series to evaluate over [s]
%  - x1; data series of s/c 1's state at times t (6xm) [km ... km /s]
%  - ephfun; struct containing f1(t) ... f4(t) to obtain s/c ephemerides
%  - oe; struct containing pertinent orbital elements for s/c
%  - moon; struct containing pertinent information about the moon
%  - debug; boolean, include stats about PDOP?

if nargin < 7, debug = false; end

m = length(t);  % number of data points

% get time-of-flight to arrive at appropriate starting TA
dt2 = f_to_t(oe.f2, oe.f2 + f(1), oe.a2, oe.e2, moon.GM);
dt3 = f_to_t(oe.f3, oe.f3 + f(2), oe.a3, oe.e3, moon.GM);
dt4 = f_to_t(oe.f4, oe.f4 + f(3), oe.a4, oe.e4, moon.GM);
% trajectories of s/c 2-4 based on state f
x2 = ephfun.Khon2(t + dt2);
x3 = ephfun.Khon3(t + dt3);
x4 = ephfun.Khon4(t + dt4);

sats = zeros(m, 3, 4);
sats(:,:,1) = x1(1:3,:)';
sats(:,:,2) = x2(1:3,:)';
sats(:,:,3) = x3(1:3,:)';
sats(:,:,4) = x4(1:3,:)';

% compute the performance index (lower PDOP -> more negative number)
PDOP = zeros(1, m);
for i=1:m
    % cap PDOP at 2000
    PDOP(i) = min(2000, computeDOP_bare(moon.R, sats, i));
end

avail = 1 - length(PDOP(PDOP >= 2000)) / length(PDOP);
P = sum(PDOP) / (m * avail^2);

if debug
    fprintf("Min PDOP : %f\n", min(PDOP));
    fprintf("Mean PDOP: %f\n", mean(PDOP(PDOP < 2000)));
    fprintf("Pct avail: %f%%\n", avail * 100);
end
end

