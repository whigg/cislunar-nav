 % conoptTester.m
% Author: Mark Hartigan
% Date  : May 28, 2024
% Description:
%    Used for development of ConOpt class to validate functionality .

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
% START = '2024 May 1 00:00:00';
START = '2027 Feb 2 00:00:00';
% Which design space to use? Same i's, a's is 0 and different is 1
DSPACE = 1;
% Initialize with previous population? 0 - no, 1 - yes, 2 - use IM as start
PREV = 1;

cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
t0 = cspice_str2et(START);

%% ConOpt execution
if DSPACE 
    conopt = ConOpt2(t0, 13, "n", 6, "dt", 300, "n_sph", 6); 
else
    conopt = ConOpt(t0, 13, "n", 6, "dt", 300, "n_sph", 6); 
end

V = conopt.T_OP2J;
W = conopt.T_J2ME;

if ~PREV
    % options = optimoptions('patternsearch','PlotFcn',@psplotbestf,'Display','iter', ...
    %     'UseParallel',true,'FunctionTolerance',1e-15,'Algorithm','nups-gps');
    options = optimoptions('ga','PlotFcn',@gaplotbestf,'Display','iter', ...
        'UseParallel',true,'FunctionTolerance',1e-6,'MaxStallGenerations',30, ...
        'PopulationSize', 200);
    
    % xopt = conopt.run(R, options);
    [xopt,fval,exitflag,output,population,scores] = conopt.run(V, W, options);
elseif PREV == 1
    load('data/optimization/RUN10 (46p).mat');
    if size(population, 2) == 24
        newpop = population;
    elseif size(population, 2) == 14
        oldpop = population;
        newpop = zeros(200, 24);
        newpop(:,1:6) = repmat(oldpop(:,1), 1, 6);
        newpop(:,7:12) = repmat(oldpop(:,2), 1, 6);
        newpop(:,13:end) = oldpop(:,3:end);
    end

    options = optimoptions('ga','PlotFcn',@gaplotbestf,'Display','iter', ...
        'UseParallel',true,'FunctionTolerance',1e-6,'MaxStallGenerations',30, ...
        'PopulationSize', 200, 'InitialPopulationMatrix', newpop);

    [xopt,fval,exitflag,output,population,scores] = conopt.run(V, W, options);

elseif PREV == 2
    load('data/optimization/IM Xopt.mat');
    newpop = xopt;

    options = optimoptions('ga','PlotFcn',@gaplotbestf,'Display','iter', ...
        'UseParallel',true,'FunctionTolerance',1e-6,'MaxStallGenerations',30, ...
        'PopulationSize', 200, 'InitialPopulationMatrix', newpop);

    [xopt,fval,exitflag,output,population,scores] = conopt.run(V, W, options);
end

conopt.coverage(xopt, V, W)

%% final evaluation

if conopt.flag == 1     % ConOpt used; orbits use same i, a
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
elseif conopt.flag == 2 % ConOpt 2 used; orbits constrained to same RAAN drift
    oe1.w = pi/2;
    oe2 = oe1; oe3 = oe1; oe4 = oe1; oe5 = oe1; oe6 = oe1;
    
    oe1.i = xopt(1); oe1.a = xopt(7);  oe1.e = frozenorbitfinder(oe1.i); oe1.RAAN = xopt(13); oe1.f = xopt(19);
    oe2.i = xopt(2); oe2.a = xopt(8);  oe2.e = frozenorbitfinder(oe2.i); oe2.RAAN = xopt(14); oe2.f = xopt(20);
    oe3.i = xopt(3); oe3.a = xopt(9);  oe3.e = frozenorbitfinder(oe3.i); oe3.RAAN = xopt(15); oe3.f = xopt(21);
    oe4.i = xopt(4); oe4.a = xopt(10); oe4.e = frozenorbitfinder(oe4.i); oe4.RAAN = xopt(16); oe4.f = xopt(22);
    oe5.i = xopt(5); oe5.a = xopt(11); oe5.e = frozenorbitfinder(oe5.i); oe5.RAAN = xopt(17); oe5.f = xopt(23);
    oe6.i = xopt(6); oe6.a = xopt(12); oe6.e = frozenorbitfinder(oe6.i); oe6.RAAN = xopt(18); oe6.f = xopt(24);
else                    % something went wrong
    return;
end

close all;
AnalyzeConstellation
