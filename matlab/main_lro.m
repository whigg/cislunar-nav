% main.m
% Author: Mark Hartigan
% Date  : April 12, 2024
% Description:
%    Main execution for LNSS user navigation simulations; final project for
%    AE 6505: Kalman filtering.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation
tic

%% INPUT FOR SIMULATIONS GO HERE
% LRO Mission Timeline:
% - NOMINAL MISSION                 2009-09-15 to 2010-09-15 (2009-258 to 2010-258)
% - SCIENCE MISSION                 2010-09-16 to 2012-09-15 (2010-259 to 2012-259)
% - EXTENDED SCIENCE MISSION        2012-09-15 to 2014-09-15 (2012-259 to 2014-258)
% - SECOND EXTENDED SCIENCE MISSION 2014-09-15 to 2016-09-15 (2014-258 to 2016-259)
% - THIRD EXTENDED SCIENCE MISSION  2016-09-15 to 2019-09-15 (2016-259 to 2019-258)
% - FOURTH EXTENDED SCIENCE MISSION 2019-09-15 to 2022-10-01 (2019-258 to 2022-274)
% - FIFTH EXTENDED SCIENCE MISSION  2022-10-01 to 2025-10-01 (2022-274 to 2025-274)

START = "15-Oct-2009 00:00:00";     % start date and time
DAYS = 1;                           % days, length of simulation
COMP_FRAME = 'J2000';               % frame to compute in
N_SPH = 10;                         % max degree and order of spherical harmonics to compute
% filter type; options are "EKF", "UKF", "Particle", "LKF", or "LUMVE"
FILTER = "EKF";
N_PART = 10;                        % Particles for particle filter

% Plot settings
PLOT_SCENE = false;                 % boolean, plot all satellite trajectories?
PLOT_DYN = true;                   % boolean, plot different dynamical errors
PLOT_FRAME = 'J2000';               % frame to plot in
LNSS_COLOR = [0 0.4470 0.7410];     % color of LNSS trajectories in plots
LRO_COLOR = [0 0 0]; % color of LRO trajectories in plots

%% load SPICE kernels and spherical harmonic data
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
cspice_furnsh(strcat(userpath,'/kernels/LRO/mk/lro.tm'));

[fLnss, namesLnss] = getlnsshandles();
nLnss = length(fLnss);

[R,C,S] = cofloader("data/LP165P.cof");

%% Generate simulation data
% planetary info
bods = getplanets("MOON", "EARTH", "SUN", "JUPITER");
bods(1).R = R * 1e-3;           % convert from m to km
bods(1).C = C;                  % store in moon struct for orbitaldynamics
bods(1).S = S;                  % store in moon struct for orbitaldynamics
bods(1).frame = 'MOON_ME';      % body-fixed frame of coefficients
moon = bods(1);                 % primary body
sec = bods(2:end);              % secondary bodies

% timing information
t0 = convertTo(datetime(START), "juliandate");
t0 = (t0 - 2451545) * 86400;    % JD to seconds past J2000
ts = t0:60:t0 + 86400*DAYS;     % simulation time steps
ns = length(ts);                % number of simulation time steps

% truth data
x_user = cspice_spkezr('LRO', ts, COMP_FRAME, 'NONE', 'MOON');
a = 2e-10 / (86400 * 30);       % Hz/Hz/s, upper end of aging rate a
x0Clk = [-0.8e-7, 0.5e-10, a];  % starting state of clock
[s1,s2,s3] = DiffCoeffCSAC();
xClk = clockStateOverTime(ts, x0Clk, "CSAC");
x_true = [x_user; xClk];
x0 = x_true(:,1);

%% generate measurement model
R = zeros(2*nLnss,2*nLnss);         % measurement covariance matrix
% per LSP SRD Table 3-11, SISE Position (ignore receiver noise)
R(1:nLnss,1:nLnss) = diag(repmat((0.01343/3)^2, 1, nLnss));
% per LSP SRD Table 3-11, SISE Velocity (ignore receiver noise)
R(nLnss+1:end,nLnss+1:end) = diag(repmat((0.0000012/3)^2, 1, nLnss));

measurements = LNSSmeasurements("BOTH",ts,fLnss,R,moon,x_true);

%% create filter
% process noise
% no uncertainty about velocity, accel ~4e-7 from experimental_process_noise.m
% stds = 3.98642855800357e-07 3.34295767996475e-07 3.81824183418184e-07
Qorb = [0 0 0 1 1 1] * (3e-7 * 100/N_SPH)^2;  
Qclk = [s1 s2 s3].^2;
Q = diag([Qorb Qclk]);

opts = odeset("RelTol", 1e-9, "AbsTol", 1e-11);
if strcmpi(FILTER,"EKF")
    filter = EKF("hybrid", @(t,x) lnss_llouserdyn(t,x,moon,N_SPH,sec), ...
                 @(t,x) lnss_llouserpartials(t,x,moon), ...
                 Q, ...
                 @(t,x) measurements.compute(t,x), ...
                 measurements.ymeas, ...
                 @(t,x) measurements.partials(t,x), ...
                 R, ts, "opts", opts);
elseif strcmpi(FILTER,"UKF")
    % UKF
    filter = UKF("hybrid", @(t,x) lnss_llouserdyn(t,x,moon,N_SPH,sec), ...
                 Q, ...
                 @(t,x) measurements.compute(t,x), ...
                 measurements.ymeas, ...
                 R, ts, "opts", opts);
elseif strcmpi(FILTER,"PARTICLE")
    filter = Particle("hybrid", @(t,x) lnss_llouserdyn(t,x,moon,N_SPH,sec), ...
                 @(t,x) lnss_llouserpartials(t,x,moon), ...
                 Q, ...
                 @(t,x) measurements.compute(t,x), ...
                 measurements.ymeas, ...
                 R, N_PART, ts, "opts", opts);
else
    error("main:filterNotFound", "FILTER must be one of approved types.")
end

% run filter
P0hat = diag([.01 .01 .01 0.00001 0.00001 0.00001 1e-7 1e-10 1e-16].^2);
x0hat = mvnrnd(x0', P0hat)';
filter.run(x0hat, P0hat);

%% performance statistics
clc;
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
close all;

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
axis([0 DAYS*24 -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Position (m)");
title("State error over time, " + FILTER);

subplot(3,1,2);
plot(tplot, verr * 1e6);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 vstd * 1e6 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 DAYS*24 -inf inf]);
% axis([0 DAYS*24 0 75]);
ylabel("Velocity (mm/s)");

subplot(3,1,3);
plot(tplot, berr * 1e6);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 bstd * 1e6 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 DAYS*24 -inf inf]);
% axis([0 DAYS*24 0 0.5]);
legend(["Error", "3\sigma bound"], "Location", "best");
ylabel("Clock phase (\mus)");
xlabel("Time (hrs)");

% save and move plots to new folder
h4name = "stateerror_" + FILTER + "_" + int2str(N_SPH);
saveas(h4, h4name, "epsc");
saveas(h4, h4name, "png");
movefile(h4name + ".eps", "../plots");
movefile(h4name + ".png", "../plots");

% plot satellite visibility
plotformat("IEEE", 0.4, "scaling", 2);
nvis = sum(measurements.visible, 1);
h5 = figure();
plot(tplot, nvis);
box off;
axis([0 DAYS*24 -inf inf]);
xlabel("Time (hrs)"); ylabel("# of satellites");
title("LNSS satellites visible to user over simulation");

%% ancillary plots
plotformat("IEEE", 1, "scaling", 2);
if PLOT_SCENE
    x_lnss = zeros(6,ns,nLnss);     % array of satellite states
    
    for i=1:nLnss           % get LNSS data for plotting
        x_lnss(:,:,i) = fLnss{i}(ts,PLOT_FRAME);
    end

    h2 = figure();
    % Display moon in trajectory plot
    R_me = 1738.1;          % km, moon equatorial radius
    R_mp = 1736;            % km, moon polar radius
    [Imoon, ~] = imread("res/Moon_HermesCelestiaMotherlode.jpg");
    [xx, yy, zz] = ellipsoid(0, 0, 0, R_me, R_me, R_mp);

    % Rotate moon from MOON_ME frame to PLOT_FRAME
    T = cspice_pxform('MOON_ME', PLOT_FRAME, ts(end));
    for i=1:size(xx,1)
        for j=1:size(xx,2)
            tmp = T * [xx(i,j); yy(i,j); zz(i,j)];
            xx(i,j) = tmp(1); yy(i,j) = tmp(2); zz(i,j) = tmp(3);
        end
    end

    globe = surf(xx, yy, -zz);
    set(globe, 'FaceColor', 'texturemap', 'CData', Imoon, 'FaceAlpha', 1, ...
        'EdgeColor', 'none');
    hold on;
    
    for i=1:nLnss           % plot LNSS satellite trajectories
        plot3(x_lnss(1,:,i), x_lnss(2,:,i), x_lnss(3,:,i), ...
            "LineWidth", 1.5, "Color", LNSS_COLOR, "Marker", "diamond", "MarkerFaceColor", LNSS_COLOR, "MarkerIndices", ns);
    end

    % plot user trajectory for same time frame
    plot3(x_user(1,:), x_user(2,:), x_user(3,:), ...
        "LineWidth", 1.5, "Color", LRO_COLOR, "Marker", "o", "MarkerFaceColor", LRO_COLOR, "MarkerIndices", ns);
    
    hold off; grid on; axis equal;
    xlabel("x_{J2000} (km)");
    ylabel("y_{J2000} (km)");
    zlabel("z_{J2000} (km)");
    title(int2str(DAYS) + "-day propagation, starting at " + START);
    legend('', 'LNSS Satellites', '', '', '', '', '', '', '', 'LRO', "Location", "best");

    h3 = figure();
    subplot(2,1,1);
    plot((ts - t0) / 3600, xClk(1,:), "LineWidth", 1.5);
    hold on;
    hold off; grid on;
    ylabel("Clock bias (s)");
    title("LLO CSAC phase offset");

    subplot(2,1,2);
    plot((ts - t0) / 3600, xClk(2,:), "LineWidth", 1.5);
    hold on;
    hold off; grid on;
    xlabel("Time (hrs)");
    ylabel("Clock drift (s/s)");
    title("LLO CSAC frequency offset");
end

if PLOT_DYN         % plotting error for different dynamics
    x_user = cspice_spkezr('LRO', ts, COMP_FRAME, 'NONE', 'MOON');
    propmodelcomparison(ts, x_user, bods);
end

%% end
fprintf("========================\n");
toc