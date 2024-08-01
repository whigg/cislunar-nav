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
load("data/optimization/RUN29 (56p,1).mat");
START = '2027 Feb 2 00:00:00';
prop = LunarPropagator(START, xopt, 3);

%% propagate
[ts,xs] = prop.run(86400*28, 4032, 'MOON_ME');