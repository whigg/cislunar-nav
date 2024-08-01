% evaNavSim.m
% Author: Mark Hartigan
% Date  : July 21, 2024
% Description:
%    Simulate the navigation of an astronaut near the lunar south pole
%    during an EVA.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% import data
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
[fLnss, namesLnss] = getlnsshandles();
nLnss = length(fLnss);

%% Generate simulation data
% planetary info
moon = getplanets("MOON");
c = 299792.458;                     % km/s, speed of light

% timing information
START = '15-Oct-2009 00:00:00';     % start date and time
t0 = cspice_str2et(START);
ts = t0:1:t0+4*3600;
tm = t0:10:ts(end);
[~,im] = intersect(ts, tm, "stable");
ns = length(ts);

% truth data
rng(420);
[pts, reftraj, exptraj] = surfacetrajgen(t0,[-pi/2.01 0 1], 4, 1.6, 30*60);
x_user = ppval(exptraj, ts);
% lat = -pi/2.01; lon = 0;
% x0 = [1737.4*[cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)] 0 0 0 0 0 0]';
% exptraj = @(t) repmat(x0, 1, length(t));
% x_user = exptraj(ts);
a_user = x_user(7:9,:);
x_user = x_user(1:6,:);
a = 2e-10 / (86400 * 30);       % Hz/Hz/s, upper end of aging rate a
x0Clk = [-0.8e-7, 0.5e-10, a];  % starting state of clock
[s1,s2,s3] = DiffCoeffCSAC();
xClk = clockStateOverTime(ts, x0Clk, "CSAC");
x_true = [x_user; xClk];
x_true_m = x_true(:,im);
x0 = x_true(:,1);

%% generate measurement model
R = zeros(2*nLnss,2*nLnss);         % measurement covariance matrix
% per LSP SRD Table 3-11, SISE Position (ignore receiver noise)
R(1:nLnss,1:nLnss) = diag(repmat((0.01343/3)^2, 1, nLnss));
% per LSP SRD Table 3-11, SISE Velocity (ignore receiver noise)
R(nLnss+1:end,nLnss+1:end) = diag(repmat((0.0000012/3)^2, 1, nLnss));
% R = diag(repmat((0.01343/3)^2, 1, nLnss));

measurements = LNSSmeasurements("BOTH",tm,fLnss,R,moon,x_true_m,"frame",'MOON_ME');

%% pseudorange-based trilateration
[x_psd, gdop] = measurements.trilaterate();

%% analyze trilateration
ps_pos = sqrt(sum((x_psd(1:3,:) - x_true_m(1:3,:)).^2, 1));

plotformat("IEEE", 0.8, "scaling", 2);

figure();
subplot(2,1,2);
plot((tm - t0)/3600, gdop);
grid on;
xlabel("Time (hrs)");
ylabel("GDOP");
subplot(2,1,1);
scatter((tm - t0)/3600, ps_pos * 1e3, 10, 'filled', 'diamond');
grid on;
ylabel("Error (m)");
title("RSS position error from multilateration");

%% filter
FILTER = "EKF";
N_PART = 200;

getaccel = @(x) x(7:9);
std_a = 80e-9*9.81;
afunc = @(t) getaccel(ppval(exptraj,t)) + mvnrnd([0 0 0]', std_a.^2);
% std_a = 0;
% afunc = @(t) getaccel(exptraj(t)) + mvnrnd([0 0 0]', std_a.^2);
dt = 1;
Qusr = std_a^2 * [1/2*dt^2*eye(3) zeros(3,3); zeros(3,3) dt*eye(3)];
Qclk = [s1^2*dt + s2^2/3*dt^3 + s3^2/20*dt^5 s2^2/2*dt^2 + s3^2/8*dt^4 s3^2/6*dt^3;
        s2^2/2*dt^2 + s3^2/8*dt^4            s2^2*dt + s3^2/3*dt^3     s3^2/2*dt^2;
        s3^2/6*dt^3                          s3^2/2*dt^2               s3^2*dt     ];
Q = [Qusr zeros(6,3); zeros(3,6) Qclk];
clc
opts = odeset("RelTol", 1e-9, "AbsTol", 1e-11);

if strcmpi(FILTER,"EKF")
    filter = EKF("discrete", @(ts,x) surfaceuserdyn(ts,x,afunc), ...
                 @surfaceuserpartials, ...
                 Q, ...
                 @(t,x) measurements.compute(t,x), ...
                 measurements.ymeas, ...
                 @(t,x) measurements.partials(t,x), ...
                 R, tm, "opts", opts, "t_sim", ts);
elseif strcmpi(FILTER,"UKF")
    % UKF
    filter = UKF("discrete", @(ts,x) surfaceuserdyn(ts,x,afunc), ...
                 Q, ...
                 @(t,x) measurements.compute(t,x), ...
                 measurements.ymeas, ...
                 R, tm, "opts", opts, "t_sim", ts);
elseif strcmpi(FILTER,"PARTICLE")
    filter = Particle("discrete", @(ts,x) surfaceuserdyn(ts,x,afunc), ...
                 @surfaceuserpartials, ...
                 Q, ...
                 @(t,x) measurements.compute(t,x), ...
                 measurements.ymeas, ...
                 R, N_PART, ts, "opts", opts);
else
    error("main:filterNotFound", "FILTER must be one of approved types.")
end

% run filter
P0hat = diag([.001 .001 .001 1e-6 1e-6 1e-6 1e-8 1e-10 1e-16].^2);
x0hat = mvnrnd(x0', P0hat)';
filter.run(x0hat, P0hat);

%% performance statistics
x_err = filter.x - x_true;
% load('res/err_100x100.mat');

rerr = sqrt(sum(x_err(1:3,:).^2, 1));
rstd = reshape(sqrt(filter.P(1,1,:)+filter.P(2,2,:)+filter.P(3,3,:))*3, 1, ns);
verr = sqrt(sum(x_err(4:6,:).^2, 1));
vstd = reshape(sqrt(filter.P(4,4,:)+filter.P(5,5,:)+filter.P(6,6,:))*3, 1, ns);
berr = abs(x_err(7,:));
bstd = reshape(sqrt(filter.P(7,7,:))*3, 1, ns);

fprintf("FILTER: " + FILTER + "\n");
fprintf("========================\n");
fprintf("  Mean error: %.3f (%.3f 3-sigma) m\n", mean(rerr)*1e3, mean(rstd)*1e3);
fprintf("              %.3f (%.3f 3-sigma) mm/s\n", mean(verr)*1e6, mean(vstd)*1e6);
fprintf("              %.3f (%.3f 3-sigma) ns\n", mean(berr)*1e9, mean(bstd)*1e9);
fprintf("Median error: %.3f (%.3f 3-sigma) m\n", median(rerr)*1e3, median(rstd)*1e3);
fprintf("              %.3f (%.3f 3-sigma) mm/s\n", median(verr)*1e6, median(vstd)*1e6);
fprintf("              %.3f (%.3f 3-sigma) ns\n", median(berr)*1e9, median(bstd)*1e9);

%% plot performance
plotformat("IEEE", 1, "scaling", 2);
h4 = figure();

tplot = (ts - t0) / 3600;
colors = colororder();
subplot(3,1,1);
plot(tplot, rerr * 1e3);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 rstd * 1e3 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 4 -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Position (m)");
title("State error over time, " + FILTER);

subplot(3,1,2);
plot(tplot, verr * 1e6);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 vstd * 1e6 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 4 -inf inf]);
% axis([0 DAYS*24 0 75]);
ylabel("Velocity (mm/s)");

subplot(3,1,3);
plot(tplot, berr * 1e6);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 bstd * 1e6 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 4 -inf inf]);
% axis([0 DAYS*24 0 0.5]);
legend(["Error", "3\sigma bound"], "Location", "best");
ylabel("Clock phase (\mus)");
xlabel("Time (hrs)");

% % save and move plots to new folder
% h4name = "stateerror_" + FILTER + "_" + int2str(N_SPH);
% saveas(h4, h4name, "epsc");
% saveas(h4, h4name, "png");
% movefile(h4name + ".eps", "../plots");
% movefile(h4name + ".png", "../plots");

% plot satellite visibility
plotformat("IEEE", 0.4, "scaling", 2);
nvis = sum(measurements.visible, 1);
h5 = figure();
plot((tm - t0) / 3600, nvis);
box off;
axis([0 4 -inf inf]);
xlabel("Time (hrs)"); ylabel("# of satellites");
title("LNSS satellites visible to user over simulation");