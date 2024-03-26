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
    nsats = 6;
    links = 4;
    volume = 'SV2';
end

% OVERWRITE PARAMETERS DOWN HERE IF DESIRED
% nsats = 1;          % number of satellites to consider (min 1, max 5)
% links = 1;          % make sure always <= nsats
% volume = 'SV1';     % options are 'SV1', 'SV2', or 'SV3'

%% Gather trajectories
t0 = convertTo(datetime('2-Feb-2027 00:00:00'), 'juliandate');
t0 = (t0 - 2451545) * 86400;    % JD to seconds past MOON_ME

imdata = {
    'data/NSNS-DRM1.0-Khon1.bsp' 'data/NSNS-DRM1.0-Khon1.tpc.pc' ...
    'data/NSNS-DRM1.0-Khon2.bsp' 'data/NSNS-DRM1.0-Khon2.tpc.pc' ...
    'data/NSNS-DRM1.0-Khon3.bsp' 'data/NSNS-DRM1.0-Khon3.tpc.pc' ...
    'data/NSNS-DRM1.0-Khon4.bsp' 'data/NSNS-DRM1.0-Khon4.tpc.pc' ...
    'data/NSNS-DRM1.4-Khon5asKhon2.bsp' 'data/NSNS-DRM1.4-Khon5asKhon2.tpc.pc' ...
    'data/NSNS-DRM1.4-Khon6asKhon4.bsp' 'data/NSNS-DRM1.4-Khon6asKhon4.tpc.pc' ...
    'data/NSNS-DRM1.1-Khon7asKhon1.bsp' 'data/NSNS-DRM1.1-Khon7asKhon1.tpc.pc'};



gendata = {'data/naif0012.tls', 'data/de430.bsp', 'data/gm_de431.tpc', ...
    'data/moon_080317.tf', 'data/moon_pa_de421_1900-2050.bpc'};
cspice_furnsh([imdata, gendata]);

moon.R  = 1736;         % km, polar radius of moon
moon.GM = cspice_bodvrd('MOON','GM',1);
moon.x  = @(t) zeros(3,1);

getpos = @(x) x(1:3,:);
Khon1 = @(t) getpos(cspice_spkezr('Khon-1', t, 'MOON_ME', 'NONE', 'MOON'));
Khon2 = @(t) getpos(cspice_spkezr('Khon-2', t, 'MOON_ME', 'NONE', 'MOON'));
Khon3 = @(t) getpos(cspice_spkezr('Khon-3', t, 'MOON_ME', 'NONE', 'MOON'));
Khon4 = @(t) getpos(cspice_spkezr('Khon-4', t, 'MOON_ME', 'NONE', 'MOON'));
Khon5 = @(t) getpos(cspice_spkezr('Khon-5', t, 'MOON_ME', 'NONE', 'MOON'));
Khon6 = @(t) getpos(cspice_spkezr('Khon-6', t, 'MOON_ME', 'NONE', 'MOON'));
Khon7 = @(t) getpos(cspice_spkezr('Khon-7', t, 'MOON_ME', 'NONE', 'MOON'));

days = 28;
step = 60*days; % seconds
t = t0:step:t0 + 86400 * days;
xk1 = Khon1(t); xk2 = Khon2(t); xk3 = Khon3(t);
xk4 = Khon4(t); xk5 = Khon5(t); xk6 = Khon6(t); xk7 = Khon7(t);
xk = [xk2; xk3; xk4; xk5; xk6; xk7];
xk = xk(1:3*nsats, :);

% MANUALLY CHANGE SATELLITE SETUP HERE
% xk = [xk1; xk2; xk3; xk4; xk5; xk6];  % include all 6 sats as nav

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
[pct, minIdx, pts, covered] = coverageLNSP(step, sats, volume, links);
plotDayCoverage(pts, covered, minIdx, step);
plotDayCoverage(pts, covered, minIdx, step, true);
fprintf("Minimum %% coverage of service volume over any 24h period: %.2f%%\n", pct*100);
fprintf("\t(starting at day %.2f)\n", (t(minIdx)-t0)/(t(end)-t0) * days);

%% EVA coverage
[nwin, minPt, minDay] = supportEVA(covered, step);
day = ceil(86400 / step);       % steps in 1 day
p = day*(minDay - 1) + 1;
q = day*minDay;
plotGDOP(pts(:,minPt), sats(p:q,:,:), t(1:day));

%% compute average GDOP at south pole
dop = 0;
m = 0;
tempdop = GDOP(pos, sats);
for i=1:size(tempdop)
    if tempdop(i) < 6
        dop = dop + tempdop(i);
        m = m + 1;
    end
end

fprintf("Average GDOP for south pole: %.3f\n", dop / m);

% %% speed test
% tic
% [pct, minIdx, pts, covered] = coverageLNSP(step, sats, volume, links);
% toc
% 
% % generate evaluation points and trim to SV latitude
% lat =  90 * pi/180;
% l = 120;
% [X,Y,Z] = mySphere(l);
% accept = Z < sin(lat);
% pts = [X; Y; Z];
% pts = pts(:, accept);   % points on surface of unit sphere below lat
% % scale to lunar distances and generate volume
% pts = pts * 1736;       % km, polar radius of moon; only care about surface
% 
% tic
% coverageLNSP_bare(step, sats, pts)
% toc

%% visualization
% visualizeCoverage(pos, sats);

%% other metrics
plotSISE(t);                    % signal-in-space errors
plotTrajectories(sats);         % plot relative orbits
% GDOP portion of plot not relevant for nsats < 4
plotGDOP(pos, sats, t, 10);         % geometric dilution of precision
% plot user navigation error only if user can directly compute position
% from measurements (i.e. nsats >= 4)
if nsats >= 4, plotUNE(pos, sats, t); end   
