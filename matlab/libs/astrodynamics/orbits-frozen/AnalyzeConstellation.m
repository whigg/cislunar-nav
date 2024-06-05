% AnalyzeConstellation.m
% Author: Mark Hartigan
% Date  : May 16, 2024
% Description:
%    Display full constellation and its statistics

% %% reset
% clc, clear, close all;
% addpath(genpath(pwd));
% format long g;          % display long numbers, no scientific notation

%% orbital parameters
% START = '2024 May 1 00:00:00';

% oe1.a = 19907.97;
% oe1.i = 2.2157;
% oe1.e = frozenorbitfinder(oe1.i);
% oe1.w = pi/2;
% oe2 = oe1; oe3 = oe1; oe4 = oe1; oe5 = oe1; oe6 = oe1;
% 
% oe1.RAAN = 2.9335;      oe1.f = 5.0115;
% oe2.RAAN = 4.6715;      oe2.f = 3.3865;
% oe3.RAAN = 4.7302;      oe3.f = 4.9999;
% oe4.RAAN = 5.9584;      oe4.f = 4.2288;
% oe5.RAAN = 6.2824;      oe5.f = 6.0517;
% oe6.RAAN = 4.7612;      oe6.f = 2.7114;

% cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
% t0 = cspice_str2et(START);

% earth information
earth.GM = cspice_bodvrd('EARTH', 'GM', 1);
earth.x = @(tau) cspice_spkpos('EARTH', tau, 'J2000', 'NONE', 'MOON');

% sun information
sun.GM = cspice_bodvrd('SUN', 'GM', 1);
sun.x = @(tau) cspice_spkpos('SUN', tau, 'J2000', 'NONE', 'MOON');

% moon information
moon.GM = cspice_bodvrd('MOON', 'GM', 1);
[R,C,S] = cofloader(userpath + "/astrodynamics/LP165P.cof");
moon.x = @(~) [0;0;0];
moon.R = R * 1e-3;      % convert from m to km
moon.C = C;             % store in moon struct for orbitaldynamics
moon.S = S;             % store in moon struct for orbitaldynamics
moon.frame = 'MOON_ME'; % body-fixed frame of coefficients

%% plot orbits
oes = [oe1 oe2 oe3 oe4 oe5 oe6];
n = length(oes);
frame = 'MOON_OP';
orbitunitsphereplot(oes, t0, moon, earth, frame);

%% propagate orbits
days = 28;
step = 600;
ts = t0:step:t0+86400*days;
m = length(ts);         % number of data points
sats = zeros(m, 6, n);

opts = odeset("RelTol", 1e-8, "AbsTol", 1e-9);

Xs = [];
for i=1:n
    [r,v] = oe2rv(oes(i).a,oes(i).e,oes(i).i,oes(i).RAAN,oes(i).w,oes(i).f,moon.GM);
    x0 = cspice_sxform('MOON_OP', 'J2000', t0) * [r; v];
    [~,X] = ode45(@(t,x) orbitaldynamics(t,x,moon,100,earth), ts, x0, opts);
    X = X';
    for j=1:length(ts)
        X(:,j) = cspice_sxform('J2000', frame, ts(j)) * X(:,j);
    end
    sats(:,:,i) = X';
end

plotLunarOrbit(ts, sats, frame, "LNSS constellation trajectories");

%% get coverage statistics
[pct, minIdx, pts, covered] = coverageLNSP(step, sats, 'SV2', 4);
plotDayCoverage(pts, covered, minIdx, step);
plotDayCoverage(pts, covered, minIdx, step, true);
fprintf("Minimum %% coverage of service volume over any 24h period: %.2f%%\n", pct*100);
fprintf("\t(starting at day %.2f)\n", (ts(minIdx)-t0)/(ts(end)-t0) * days);

%% EVA coverage
[nwin, minPt, minDay] = supportEVA(covered, step);
day = ceil(86400 / step);       % steps in 1 day
p = day*(minDay - 1) + 1;
q = day*minDay;
% coverage of worst point during worst interval
plotGDOP(pts(:,minPt), sats(p:q,:,:), ts(1:day), 10);
% coverage of LSP for entire duration
plotGDOP([0;0;-1738], sats, ts, 10);
