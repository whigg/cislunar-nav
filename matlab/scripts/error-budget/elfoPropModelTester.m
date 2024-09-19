% elfoPropModelTester.m
% Author: Mark Hartigan
% Date  : June 26, 2024
% Description:
%    Test the accuracy of propagation model orbitaldynamics in elliptical
%    lunar frozen orbits.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
% starting epoch
START = '2024 May 1 00:00:00';
% frame of data
FRAME = 'J2000';
% number of days in simulation
DAYS = 7;

% load ELFO orbit 
cspice_furnsh('data/gmat-to-spk/ELFO_20240501-20240531.bsp');

ELFO = '-909';                  % SPICE ID of ELFO s/c
t0 = cspice_str2et(START);

%% generate plot
% planetary info
[R,C,S] = cofloader("data/LP165P.cof");
bods = getplanets("MOON", "EARTH", "SUN", "JUPITER");
bods(1).R = R * 1e-3;           % convert from m to km
bods(1).C = C;                  % store in moon struct for orbitaldynamics
bods(1).S = S;                  % store in moon struct for orbitaldynamics
bods(1).frame = 'MOON_ME';      % body-fixed frame of coefficients

% truth data
ts = t0:60:t0 + 86400*DAYS;
xs = cspice_spkezr(ELFO, ts, 'J2000', 'NONE', 'MOON');
plotLunarOrbit(ts, xs', 'J2000', "Truth");
propmodelcomparison(ts, xs, bods);