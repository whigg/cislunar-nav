% conoptTester.m
% Author: Mark Hartigan
% Date  : May 28, 2024
% Description:
%    Used for development of ConOpt class to validate functionality.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
START = '2024 May 1 00:00:00';
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
t0 = cspice_str2et(START);

%% ConOpt
conopt = ConOpt(t0, 9, "n", 6, "dt", 900);

a = 11092;
i = 1.094;
e = frozenorbitfinder(i);
RAANs = [0.7945 2.3469 0.8260 0.7787 2.3661 2.3579];
TAs = [2.7586 2.7586 2.8758 1.8872 3.1030 2.8813];

% R = conopt.T_OP2ME;
R = conopt.T_OP2J;
S = conopt.T_J2ME;
% 
% CVG = conopt.coverage([i,a,RAANs,TAs], R, S);
% [sats,fail] = conopt.orbitgenerator(i,a,RAANs,TAs, conopt.T_OP2ME);
% plotLunarOrbit(conopt.ts,sats,'MOON_OP',"test");
% [pct, minIdx, pts, covered] = coverageLNSP(conopt.dt, sats, "SV2", 4);

%% Execution
if false
    % options = optimoptions('patternsearch','PlotFcn',@psplotbestf,'Display','iter', ...
    %     'UseParallel',true,'FunctionTolerance',1e-15,'Algorithm','nups-gps');
    options = optimoptions('ga','PlotFcn',@gaplotbestf,'Display','iter', ...
        'UseParallel',true,'FunctionTolerance',1e-6,'MaxStallGenerations',50, ...
        'PopulationSize', 200);
    
    % xopt = conopt.run(R, options);
    [xopt,fval,exitflag,output,population,scores] = conopt.run(R, S, options);
else
    load('data/RUN7 (28p).mat');
end

%%
conopt.coverage(xopt, R, S)
[sats,fail] = conopt.orbitgenerator(xopt(1),xopt(2),xopt(3:8),xopt(9:14), R, S);
[pct, minIdx, pts, covered] = coverageLNSP(conopt.dt, sats, "SV2", 4);

%% Final evaluation
oe1.i = xopt(1);
oe1.a = xopt(2);
oe1.e = frozenorbitfinder(oe1.i);
oe1.w = pi/2;
oe2 = oe1; oe3 = oe1; oe4 = oe1; oe5 = oe1; oe6 = oe1;

oe1.RAAN = xopt(3); oe1.f = xopt(9);
oe2.RAAN = xopt(4); oe2.f = xopt(10);
oe3.RAAN = xopt(5); oe3.f = xopt(11);
oe4.RAAN = xopt(6); oe4.f = xopt(12);
oe5.RAAN = xopt(7); oe5.f = xopt(13);
oe6.RAAN = xopt(8); oe6.f = xopt(14);

close all;
AnalyzeConstellation
