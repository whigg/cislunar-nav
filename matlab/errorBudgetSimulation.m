% errorBudgetSimulation.m
% Author: Mark Hartigan
% Date  : June 10, 2024
% Description:
%    Chaining together simulations to develop a demonstration of error
%    budgets for lunar navigation systems.

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
DAYS = 0.5;
% number of spherical harmonics to use in EKF
N_SPH = 16;
% should filter performance statistics and plots be displayed?
FILTER_PERFORMANCE = true;

% load generic lunar information
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
% load ELFO orbit 
cspice_furnsh('data/gmat-to-spk/ELFO_20240501-20240531.bsp');
% load GNSS satellite info
cspice_furnsh(strcat(userpath,'/kernels/GNSS/mk/gnss.tm'));

ELFO = '-909';                  % SPICE ID of ELFO s/c
t0 = cspice_str2et(START);

%% generate truth orbit
ns = 1440;
ts = linspace(t0, t0+86400*DAYS, ns);
sat = zeros(9,ns);

for i=1:ns
    sat(1:6,i) = cspice_spkezr(ELFO, ts(i), FRAME, 'NONE', 'MOON');
end

x0_OP = cspice_sxform(FRAME, 'MOON_OP', t0) * sat(1:6,1);
moon.GM = cspice_bodvrd('MOON', 'GM', 1);
[oe.a,oe.e,oe.i,oe.RAAN,oe.w,oe.f] = rv2oe(x0_OP(1:3), x0_OP(4:6), moon.GM);
% plotLunarOrbit(ts, sat', FRAME, "ELFO");

a = 2e-10 / (86400 * 30);       % Hz/Hz/s, upper end of aging rate a
x0Clk = [-0.8e-7, 0.5e-10, a];  % starting state of clock
[s1,s2,s3] = DiffCoeffCSAC();
xClk = clockStateOverTime(ts, x0Clk, "CSAC");
sat(7:9,:) = xClk;
x0 = sat(:,1);

%% get GNSS satellite setup
[f_GNSS, names_GNSS] = getgnsshandles("BOTH");
n_GNSS = length(f_GNSS);

[R,C,S] = cofloader("data/LP165P.cof");

earth.GM = cspice_bodvrd('EARTH', 'GM', 1);
earth.x  = @(tau) cspice_spkpos('EARTH', tau, 'J2000', 'NONE', 'MOON');
earth.R  = 6378 + 800;   % km, radius of earth + ionosphere

sun.GM = cspice_bodvrd('SUN', 'GM', 1);
sun.x  = @(tau) cspice_spkpos('SUN', tau, 'J2000', 'NONE', 'MOON');
sun.R  = 0;

jupiter.GM = cspice_bodvrd('JUPITER BARYCENTER', 'GM', 1);
jupiter.x  = @(tau) cspice_spkpos('JUPITER BARYCENTER', tau, 'J2000', 'NONE', 'MOON');
jupiter.R  = 0;

bodies = [earth sun jupiter];

moon.GM = cspice_bodvrd('MOON', 'GM', 1);
moon.x = @(~) [0;0;0];
moon.R = R * 1e-3;      % convert from m to km
moon.C = C;             % store in moon struct for orbitaldynamics
moon.S = S;             % store in moon struct for orbitaldynamics
moon.frame = 'MOON_ME'; % body-fixed frame of coefficients

% generate measurement model
R = zeros(2*n_GNSS,2*n_GNSS);       % measurement covariance matrix
% per LSP SRD Table 3-11, SISE Position (ignore receiver noise)
R(1:n_GNSS,1:n_GNSS) = diag(repmat((9.7/1.96*1e-3)^2, 1, n_GNSS));
% per LSP SRD Table 3-11, SISE Velocity (ignore receiver noise)
R(n_GNSS+1:end,n_GNSS+1:end) = diag(repmat((.006/1.96*1e-3)^2, 1, n_GNSS));

measurements = GNSSmeasurements("BOTH",ts,f_GNSS,R,earth,moon,sat);

%% create filter
% process noise
% no uncertainty about velocity, accel ~4e-7 from experimental_process_noise.m
% stds = 3.98642855800357e-07 3.34295767996475e-07 3.81824183418184e-07
Qorb = [0 0 0 1 1 1] * (2e-6)^2;  
Qclk = [s1 s2 s3].^2;
Q = diag([Qorb Qclk]);

opts = odeset("RelTol", 1e-9, "AbsTol", 1e-11);
filter = EKF("hybrid", @(t,x) lnss_elfodyn(t,x,moon,N_SPH,bodies), ...
             @(~,x) lnss_elfopartials(x,moon), ...
             Q, ...
             @(t,x) measurements.compute(t,x), ...
             measurements.ymeas, ...
             @(t,x) measurements.partials(t,x), ...
             R, ts, "opts", opts);

% run filter
P0hat = diag([.03 .03 .03 3e-10 3e-10 3e-10 7e-13 1e-21 1e-32]);
x0hat = mvnrnd(x0', P0hat)';
filter.run(x0hat, P0hat);

%% convert error data to RTN
x_err = filter.x - sat;
x_err_rtn = x_err;
P_rtn = filter.P;
for i=1:ns
    ti = ts(i);
    r_r = sat(1:3,i);
    u_r = r_r / norm(r_r);          % radial direction
    r_n = cross(r_r, sat(4:6,i));
    u_n = r_n / norm(r_n);          % normal direction
    u_t = cross(u_n, u_r);          % tangential direction
    T_J2RTN = [u_r u_t u_n]';
    T = [T_J2RTN zeros(3,3); zeros(3,3) T_J2RTN];
    x_err_rtn(1:6,i) = T * x_err(1:6,i);
    P_rtn(1:6,1:6,i) = T * filter.P(1:6,1:6,i) * T';
end

%% performance statistics
if FILTER_PERFORMANCE
    filterperformance("EKF", ts, x_err, filter.P, ...
                      x_err_rtn, P_rtn, measurements.visible);
end

%% ephemeris fit
t02 = ts(end);              % starting point of ephemeris propagation
x02 = filter.x(:,end);      % starting state of ephemeris propagation

N_APPX = 14;
% span = linspace(-1, 1, N_APPX - 1);
span = chebichev(N_APPX);
tf2 = t02 + 3600 * 4;
t_interp = t02 + (span + 1) * (tf2 - t02) / 2;
t_eval = linspace(t02, tf2, ns);
[tspan, i_interp, ~] = union(t_interp, t_eval, "sorted");
[T,X] = ode45(@(t,x) lnss_elfodyn(t,x,moon,150,bodies), tspan, x02, opts);

basis = @(tau) tau.^(0:N_APPX);
% basis = @(tau) [tau.^(0:N_APPX/3-1) sin(2*pi*tau*(0:N_APPX/3-1)) cos(2*pi*tau*(0:N_APPX/3-1))];
phi = basis(span');
B = pinv(phi);
c_x = B * X(ismember(tspan, t_interp),1);

x_appx = basis(2*(T - t02)/(tf2 - t02) - 1) * c_x;

err_x = x_appx - X(:,1);

plotformat("IEEE", 0.4, "scaling", 2);
figure();
plot((tspan - t02)/3600, err_x * 1e3);
grid on;
xlabel("Time (hrs)"); ylabel("Error (m)");





% TODO:
% - add additional observable (ranging to moon simulating opnav?) for nav
%   solution improvement
% - fit ephemeris polynomials
