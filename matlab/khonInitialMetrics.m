% khonInitialMetrics.m
% Author: Mark Hartigan
% Date  : October 13, 2023
% Description:
%    Computes initial performance metrics for the Khon PNT constellation.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

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
% xk = [xk2; xk3; xk4; xk5; xk6];
xk = [xk2; xk3; xk4; xk5; xk6];

%% write data
data = [t' xk'];

% fid = fopen('Khon_4sat_1day.txt', 'w');
% for i=1:size(data,1)
%     fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', data(i,2:end));
% end
% fclose(fid);

%% format data for metrics
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
[pct, pts, covered] = coverageLNSP(sats, 'SV2', 4);
plotCoverage(pts, covered);
plotCoverage(pts, covered, true);
fprintf("Percent coverage of service volume: %.2f\n", pct*100);

%% other metrics
plotGDOP(pos, sats, t);         % geometric dilution of precision
plotSISE(t);                    % signal-in-space errors
plotUNE(pos, sats, t);          % plot user navigation error
plotTrajectories(sats);         % plot relative orbits
