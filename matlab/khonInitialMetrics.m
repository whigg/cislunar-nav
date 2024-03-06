% khonInitialMetrics.m
% Author: Mark Hartigan
% Date  : October 13, 2023
% Description:
%    Computes initial performance metrics for the Khon PNT constellation.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% SET OPTIONS HERE
% Choose the operating increment here by setting the IOC constant to either
% 'A', 'B', or 'C'. These correspond to the LunaNet Service Providers
% increments. Alternatively, do not set IOC and instead set nsats, links,
% and volume directly.
IOC = 'C';

if strcmp(IOC, 'A')
    nsats = 1;
    links = 1;
    volume = 'SV1';
elseif strcmp(IOC, 'B')
    nsats = 3;
    links = 2;
    volume = 'SV1';
elseif strcmp(IOC, 'C')
    nsats = 5;
    links = 4;
    volume = 'SV2';
end

% OVERWRITE PARAMETERS DOWN HERE IF DESIRED
% nsats = 1;          % number of satellites to consider (min 1, max 5)
% links = 1;          % make sure always <= nsats
% volume = 'SV1';     % options are 'SV1', 'SV2', or 'SV3'

%% init
t0 = convertTo(datetime('2-Feb-2027 00:00:00'), 'juliandate');
t0 = (t0 - 2451545) * 86400;    % JD to seconds past J2000

imdata = cell(1, 12);
for i=0:5
    imdata{2*i + 1} = sprintf('data/NSNS-DRM1.0-Khon%d.bsp', i+1);
    imdata{2*i + 2} = sprintf('data/NSNS-DRM1.0-Khon%d.tpc.pc', i+1);
end

gendata = {'data/naif0012.tls', 'data/de430.bsp', 'data/gm_de431.tpc'};
cspice_furnsh([imdata, gendata]);

moon.R  = 1736;         % km, polar radius of moon
moon.GM = cspice_bodvrd('MOON','GM',1);
moon.x  = @(t) zeros(3,1);

getpos = @(x) x(1:3,:);
Khon1 = @(t) getpos(cspice_spkezr('Khon-1', t, 'J2000', 'NONE', 'MOON'));
Khon2 = @(t) getpos(cspice_spkezr('Khon-2', t, 'J2000', 'NONE', 'MOON'));
Khon3 = @(t) getpos(cspice_spkezr('Khon-3', t, 'J2000', 'NONE', 'MOON'));
Khon4 = @(t) getpos(cspice_spkezr('Khon-4', t, 'J2000', 'NONE', 'MOON'));
Khon5 = @(t) getpos(cspice_spkezr('Khon-5', t, 'J2000', 'NONE', 'MOON'));
Khon6 = @(t) getpos(cspice_spkezr('Khon-6', t, 'J2000', 'NONE', 'MOON'));

%% build trajectories
days = 28;
step = 60*days; % seconds
t = t0:step:t0 + 86400 * days;
xk1 = Khon1(t); xk2 = Khon2(t); xk3 = Khon3(t);
xk4 = Khon4(t); xk5 = Khon5(t); xk6 = Khon6(t);
xk = [xk2; xk3; xk4; xk5; xk6];
xk = xk(1:3*nsats, :);

%% format data for metrics
data = [t' xk'];

% format into usable shape by metric utilities
t = t';                         % convert to column vector
m = length(t);                  % number of data points
n = (size(data,2) - 1) / 3;     % number of satellites

sats = zeros(m, 3, n);
for num = 1:n
    sats(:,:,num) = data(:,2 + 3*(num-1):1 + 3*num);
end

pos = [0; 0; -moon.R];          % position of user

%% coverage
[pct, pts, covered] = coverageLNSP(sats, volume, links);
plotCoverage(pts, covered);
plotCoverage(pts, covered, true);
fprintf("Percent coverage of service volume: %.2f%%\n", pct*100);

%% other metrics
plotSISE(t);                    % signal-in-space errors
plotTrajectories(sats);         % plot relative orbits
% GDOP portion of plot not relevant for nsats < 4
plotGDOP(pos, sats, t);         % geometric dilution of precision
% plot user navigation error only if user can directly compute position
% from measurements (i.e. nsats >= 4)
if nsats >= 4, plotUNE(pos, sats, t); end   
