% constStationkeeping.m
% Author: Mark Hartigan
% Date  : August 1, 2024
% Description:
%    Manage stationkeeping for constellation geometry to keep the optimal
%    geometry.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
load("data/optimization/IM Xopt.mat");
START = '2027 Feb 2 00:00:00';
DAYS = 30;
oes = xopt2oes(xopt);
prop = LunarPropagator(START, oes, 32, 3);

%% propagate
[ts,xs] = prop.run(86400*DAYS, 4032, 'MOON_ME');
[Rcomp, Rexp] = prop.computedriftrates();
prop.plotlastorbits('MOON_OP');

%% roughly compute stationkeeping cost
% Solve the minimax problem where max(Rcomp(i)-x) is minimized. Thus, this
% correlates to minimizing the maximum amount of Delta-V required on any sat.
minopts = optimoptions('fminimax', 'FunctionTolerance', 1e-11, ...
                       'StepTolerance', 1e-12, 'Display', 'off');
Rtarg = fminimax(@(x) abs(Rcomp - x), Rcomp(1), [],[],[],[],[],[],[], minopts);
dvs = zeros(size(Rcomp));

% compute dV/day for each orbit to maintain targeted drift rate
for j=1:length(oes)
    i = oes(j).i;
    w = oes(j).w;
    R1 = oes(j).RAAN;
    R2 = oes(j).RAAN + (Rtarg - Rcomp(j)) * 86400;
    % Keplerian orbits w/ different RAANs (all else equal) intersect at 2
    % diff points. To find, set |r| equal and solve in inertial frame (see
    % 'GT OneNote>Cislunar PNT>Change of right ascension' for full derivation).
    fint = atan2((cos(R2)-cos(R1))*cos(w) - (sin(R2)-sin(R1))*cos(i)*sin(w), ...
                 -((cos(R1)+cos(R2))*sin(w) + (sin(R1)+sin(R2))*cos(i)*cos(w)));
    % change the maneuver point to whatever is closer to apolune
    if abs(fint) < pi/2, fint = fint + pi; end
    % Compute impulse dV required
    [~,v1] = oe2rv(oes(1).a,oes(1).e,i,R1,w,fint,prop.moon.GM);
    [r,v2] = oe2rv(oes(2).a,oes(2).e,i,R2,w,2*pi-fint,prop.moon.GM);
    dvs(j) = norm(v2 - v1) * 1e3;
end

fprintf("Max Delta-V: %.2f m/s/day\n", max(dvs));
fprintf("     ... or: %.3f km/s over 10 years\n", max(dvs)*365*10/1e3);
