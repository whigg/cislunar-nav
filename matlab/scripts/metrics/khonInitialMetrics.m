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
    'spice/nsns/NSNS-DRM1.0-Khon1.bsp' 'spice/nsns/NSNS-DRM1.0-Khon1.tpc.pc' ...
    'spice/nsns/NSNS-DRM1.0-Khon2.bsp' 'spice/nsns/NSNS-DRM1.0-Khon2.tpc.pc' ...
    'spice/nsns/NSNS-DRM1.0-Khon3.bsp' 'spice/nsns/NSNS-DRM1.0-Khon3.tpc.pc' ...
    'spice/nsns/NSNS-DRM1.0-Khon4.bsp' 'spice/nsns/NSNS-DRM1.0-Khon4.tpc.pc' ...
    'spice/nsns/NSNS-DRM1.4-Khon5asKhon2.bsp' 'spice/nsns/NSNS-DRM1.4-Khon5asKhon2.tpc.pc' ...
    'spice/nsns/NSNS-DRM1.4-Khon6asKhon4.bsp' 'spice/nsns/NSNS-DRM1.4-Khon6asKhon4.tpc.pc' ...
    'spice/nsns/NSNS-DRM1.1-Khon7asKhon1.bsp' 'spice/nsns/NSNS-DRM1.1-Khon7asKhon1.tpc.pc'};



gendata = {'spice/nsns/naif0012.tls', 'spice/nsns/de430.bsp', 'spice/nsns/gm_de431.tpc', ...
    'spice/nsns/moon_080317.tf', 'spice/nsns/moon_pa_de421_1900-2050.bpc'};
cspice_furnsh(imdata);
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));

moon.R  = 1736;         % km, polar radius of moon
moon.GM = cspice_bodvrd('MOON','GM',1);
moon.x  = @(t) zeros(3,1);

getpos = @(x) x(1:3,:);
Khon1 = @(t) getpos(cspice_spkezr('Khon-1', t, 'MOON_OP', 'NONE', 'MOON'));
Khon2 = @(t) getpos(cspice_spkezr('Khon-2', t, 'MOON_OP', 'NONE', 'MOON'));
Khon3 = @(t) getpos(cspice_spkezr('Khon-3', t, 'MOON_OP', 'NONE', 'MOON'));
Khon4 = @(t) getpos(cspice_spkezr('Khon-4', t, 'MOON_OP', 'NONE', 'MOON'));
Khon5 = @(t) getpos(cspice_spkezr('Khon-5', t, 'MOON_OP', 'NONE', 'MOON'));
Khon6 = @(t) getpos(cspice_spkezr('Khon-6', t, 'MOON_OP', 'NONE', 'MOON'));
Khon7 = @(t) getpos(cspice_spkezr('Khon-7', t, 'MOON_OP', 'NONE', 'MOON'));

v2 = cspice_spkezr('Khon-2', t0, 'MOON_OP', 'NONE', 'MOON');
[a2,e2,i2,r2,w2,f2] = rv2oe(v2(1:3),v2(4:6),moon.GM);
v3 = cspice_spkezr('Khon-3', t0, 'MOON_OP', 'NONE', 'MOON');
[a3,e3,i3,r3,w3,f3] = rv2oe(v3(1:3),v3(4:6),moon.GM);
v4 = cspice_spkezr('Khon-4', t0, 'MOON_OP', 'NONE', 'MOON');
[a4,e4,i4,r4,w4,f4] = rv2oe(v4(1:3),v4(4:6),moon.GM);
v5 = cspice_spkezr('Khon-5', t0, 'MOON_OP', 'NONE', 'MOON');
[a5,e5,i5,r5,w5,f5] = rv2oe(v5(1:3),v5(4:6),moon.GM);
v6 = cspice_spkezr('Khon-6', t0, 'MOON_OP', 'NONE', 'MOON');
[a6,e6,i6,r6,w6,f6] = rv2oe(v6(1:3),v6(4:6),moon.GM);
v7 = cspice_spkezr('Khon-7', t0, 'MOON_OP', 'NONE', 'MOON');
[a7,e7,i7,r7,w7,f7] = rv2oe(v7(1:3),v7(4:6),moon.GM);
xopt = [[i2 i3 i4 i5 i6 i7], [a2 a3 a4 a5 a6 a7], r2, r3, r4, r5, r6, r7, f2, f3, f4, f5, f6, f7];

days = 28;
step = 600; % seconds
t = t0:step:t0 + 86400 * days;
xk1 = Khon1(t); xk2 = Khon2(t); xk3 = Khon3(t);
xk4 = Khon4(t); xk5 = Khon5(t); xk6 = Khon6(t); xk7 = Khon7(t);
xk = [xk2; xk3; xk4; xk5; xk6; xk7];
xk = xk(1:3*nsats, :);

% for i=1:n
%     R = cspice_pxform('J2000', 'MOON_ME', t(i));
%     x(:,i) = R * x(:,i);
% end

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
 
% [sats,dt] = reshape_report("spice/nsns/8sat_adapted_2NRHO_Report_1month.txt");
% t = (0:size(sats,1)-1)*dt + t0;

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
plotGDOP(pts(:,minPt), sats(p:q,:,:), t(1:day), 10);

%% compute average GDOP at south pole
dop = 0;
m = 0;
tempdop = GDOP(pos, sats);
for i=1:length(tempdop)
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
% plotTrajectories(sats);         % plot relative orbits
plotLunarOrbit(t, sats, 'MOON_OP', "NSNS Constellation");
% GDOP portion of plot not relevant for nsats < 4
plotGDOP(pos, sats, t, 10);         % geometric dilution of precision
% plot user navigation error only if user can directly compute position
% from measurements (i.e. nsats >= 4)
if nsats >= 4, plotUNE(pos, sats, t); end   
