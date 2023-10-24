% TAOptimization.m
% Author: Mark Hartigan
% Date  : October 13, 2023
% Description:
%    Script to optimize Intuitive Machines' Khon constellation based on the
%    specified objective function perf_index().


%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
t0 = convertTo(datetime('9-Feb-2027 00:00:00'), 'juliandate');
t0 = (t0 - 2451545) * 86400;    % JD to seconds past J2000

imdata = cell(1, 12);
for i=0:5
    imdata{2*i + 1} = sprintf('data/NSNS-DRM1.0-Khon%d.bsp', i+1);
    imdata{2*i + 2} = sprintf('data/NSNS-DRM1.0-Khon%d.tpc.pc', i+1);
end

gendata = {'data/naif0012.tls', 'data/de430.bsp', 'data/gm_de431.tpc'};
cspice_furnsh([imdata, gendata]);

% define parameters for moon
moon.R  = 1737.4;       % km, equatorial radius of moon
moon.GM = cspice_bodvrd('MOON','GM',1);
moon.x  = @(t) zeros(3,1);

% define functions to retrieve spacecraft ephemerides
ephfun.Khon1 = @(t) cspice_spkezr('Khon-1', t, 'J2000', 'NONE', 'MOON');
ephfun.Khon2 = @(t) cspice_spkezr('Khon-2', t, 'J2000', 'NONE', 'MOON');
ephfun.Khon3 = @(t) cspice_spkezr('Khon-3', t, 'J2000', 'NONE', 'MOON');
ephfun.Khon4 = @(t) cspice_spkezr('Khon-4', t, 'J2000', 'NONE', 'MOON');

xk1 = ephfun.Khon1(t0); xk2 = ephfun.Khon2(t0);
xk3 = ephfun.Khon3(t0); xk4 = ephfun.Khon4(t0);

[oe.a1, oe.e1, ~, ~, ~, oe.f1] = rv2oe(xk1(1:3), xk1(4:6), moon.GM);
[oe.a2, oe.e2, ~, ~, ~, oe.f2] = rv2oe(xk2(1:3), xk2(4:6), moon.GM);
[oe.a3, oe.e3, ~, ~, ~, oe.f3] = rv2oe(xk3(1:3), xk3(4:6), moon.GM);
[oe.a4, oe.e4, ~, ~, ~, oe.f4] = rv2oe(xk4(1:3), xk4(4:6), moon.GM);

t = t0:60:t0+86400;     % 1 day, time steps every minute
x1 = ephfun.Khon1(t);

%% optimization
warning('off', 'MATLAB:nearlySingularMatrix');
% function to optimize
F = @(z) perf_index(z, t, x1, ephfun, oe, moon);

m = 100;     % number of designs
X = lhsdesign(m, 3);
X = X' * 2*pi;
Xopt = zeros(size(X));

for i=1:m
    % options = optimset('PlotFcns',@optimplotfval);
    Xopt(:,i) = fminsearch(F, X(:,i)); % , options);
    disp(F(Xopt(:,i)));
end

%% minimizer selection
Fopt = zeros(1, size(Xopt, 2));
for i=1:m
    Fopt(i) = F(Xopt(:,i));
end

[~, j] = min(Fopt);
x0 = [0; 0; 0];
xf = mod(Xopt(:,j), 2*pi);

fprintf("\n=== Initial Design ===\n");
perf_index(x0, t, x1, ephfun, oe, moon, true);
fprintf("\n=== Optimized Design ===\n");
perf_index(xf, t, x1, ephfun, oe, moon, true);

%% constellation_metrics
[~, x2, x3, x4] = F(x0);
data = [t' [x1(1:3,:); x2(1:3,:); x3(1:3,:); x4(1:3,:)]'];
constellation_metrics(data);

[~, x2, x3, x4] = F(xf);
data = [t' [x1(1:3,:); x2(1:3,:); x3(1:3,:); x4(1:3,:)]'];
constellation_metrics(data);





