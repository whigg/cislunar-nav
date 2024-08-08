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
START = '2027 Feb 2 00:00:00';
% Which design space to use? Same i's, a's is 0 and different is 1
DSPACE = 1;
% Initialize with previous population? 0 - no, 1 - yes, 2 - use IM as start
PREV = 2;

cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
t0 = cspice_str2et(START);

%% ConOpt execution
if DSPACE           % same inclination and semimajor axis (ConOpt)
    conopt = ConOpt2(t0, 13, "n", 6, "dt", 300, "n_sph", 6); 
else                % can be different (ConOpt2)
    conopt = ConOpt(t0, 13, "n", 6, "dt", 300, "n_sph", 6); 
end

V = conopt.T_OP2J;
W = conopt.T_J2ME;

if ~PREV            % no previous info
    options = optimoptions('ga','PlotFcn',@gaplotbestf,'Display','iter', ...
        'UseParallel',true,'FunctionTolerance',1e-6,'MaxStallGenerations',30, ...
        'PopulationSize', 200);
    
    [xopt,fval,exitflag,output,population,scores] = conopt.run(V, W, options);
    
elseif PREV == 1    % use old optimization info
    load('data/optimization/RUN29 (56p,1).mat');
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

elseif PREV == 2    % start from IM baseline constellation
    load('data/optimization/IM Xopt.mat');
    newpop = xopt;

    options = optimoptions('ga','PlotFcn',@gaplotbestf,'Display','iter', ...
        'UseParallel',true,'FunctionTolerance',1e-6,'MaxStallGenerations',30, ...
        'PopulationSize', 200, 'InitialPopulationMatrix', newpop);

    [xopt,fval,exitflag,output,population,scores] = conopt.run(V, W, options);
end

%% final evaluation

conopt.coverage(xopt, V, W)
oes = xopt2oes(xopt);

close all;
AnalyzeConstellation
