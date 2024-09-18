% keplerModelExperiment.m
% Author: Mark Hartigan
% Date  : August 12, 2024
% Description:
%    Test various theories for a Keplerian model

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% initialize data
% starting epoch
START = '2024 May 1 00:00:00';
% load generic lunar information
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
GM = cspice_bodvrd('MOON', 'GM', 1);
t0 = cspice_str2et(START);

%% initialize filter
% starting state of ELFO, from Cortinovis et al.
% a = 1850;                       % km, semimajor axis
% e = 0.05;                       % eccentricity
% i = 99 * pi/180;                % rad, inclination
% RAAN = 82.5 * pi/180;           % rad, right ascension of ascending node
% w = 294 * pi/180;               % rad, argument of perilune
% a = 6540;                       % km, semimajor axis
% e = 0.6;                        % eccentricity
% i = 56.2 * pi/180;              % rad, inclination
% RAAN = 180 * pi/180;            % rad, right ascension of ascending node
% w = 90 * pi/180;                % rad, argument of perilune
load("data/optimization/RUN29 (56p,1).mat");
oes = xopt2oes(xopt);
a = oes(1).a; e = oes(1).e; i = oes(1).i; RAAN = oes(1).RAAN; w = oes(1).w;
P = 2*pi * sqrt(a^3/GM);        % s, period of orbit
[r,v] = oe2rv(a, e, i, RAAN, w, pi/2, GM);
% initialize propagator (garbage starting state)
ephprop = LunarPropagator(t0, [r; v], 50, 3);

%%
tf = 3600*1;                    % final time of ephemeris propagation
f_eph = ephprop.ephemerisfit("Kepler", tf, 0);

% truth orbit
[ts,x_true] = ephprop.run(tf, 2000, 'J2000');
% ephemeris approximation errors
x_eph = f_eph(ts);
err_x = x_eph - x_true;
perr = sqrt(sum(err_x(1:3,:).^2, 1));
verr = sqrt(sum(err_x(4:6,:).^2, 1));

as_true = zeros(1, length(ts)); as_est = as_true;
es_true = zeros(1, length(ts)); es_est = es_true;
is_true = zeros(1, length(ts)); is_est = is_true;
Rs_true = zeros(1, length(ts)); Rs_est = Rs_true;
ws_true = zeros(1, length(ts)); ws_est = ws_true;
fs_true = zeros(1, length(ts)); fs_est = fs_true;

for i=1:length(ts)
    [at, et, it, Rt, wt, ft] = rv2oe(x_true(1:3,i), x_true(4:6,i), GM);
    [ae, ee, ie, Re, we, fe] = rv2oe(x_eph(1:3,i), x_eph(4:6,i), GM);
    as_true(i) = at; as_est(i) = ae;
    es_true(i) = et; es_est(i) = ee;
    is_true(i) = it; is_est(i) = ie;
    Rs_true(i) = Rt; Rs_est(i) = Re;
    ws_true(i) = wt; ws_est(i) = we;
    fs_true(i) = ft; fs_est(i) = fe;
end


%% plotting
plotformat("IEEE", 0.9, "scaling", 1.5);
figure();
tiledlayout(2,1);
nexttile
plot((ts - t0) / 3600, perr * 1e3);
grid on;
ylabel("Error (m)");
nexttile
plot((ts  -t0) / 3600, verr * 1e6);
grid on;
xlabel("Time (hrs)");
ylabel("Error (mm/s)");


plotformat("IEEE", 0.6, "scaling", 2);
figure();
tiledlayout(2,3);
nexttile
plot((ts - t0) / 3600, as_true - as_est);
ylabel("Semimajor axis error (km)");
nexttile
plot((ts - t0) / 3600, es_true - es_est);
ylabel("Eccentricity error");
nexttile
plot((ts - t0) / 3600, is_true - is_est);
ylabel("Inclination error (rad)");
nexttile
plot((ts - t0) / 3600, Rs_true - Rs_est);
xlabel("Time (hrs)");
ylabel("RAAN error (rad)");
nexttile
plot((ts - t0) / 3600, ws_true - ws_est);
xlabel("Time (hrs)");
ylabel("Arg. of perilune error (rad)");
nexttile
plot((ts - t0) / 3600, fs_true - fs_est);
xlabel("Time (hrs)");
ylabel("True anomaly error (rad)");