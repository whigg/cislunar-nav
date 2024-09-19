% errorPropagationTester.m
% Author: Mark Hartigan
% Date  : June 24, 2024
% Description:
%    Test the accuracy of uncertainty propagation using state partials

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
DAYS = 1/3;
% number of iterations for computing statistics
N_ITER = 100;
% OE # (orbit from RUN16 of optimizer)
OE = 1;
% True anomaly of starting point (deg)
F = 270;

% load generic lunar information
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));

%% Generate truth and estimate data
fname = sprintf("data/error/stats_%d_OE%d_f=%d.mat", N_ITER, OE, F);

try             % try to load pre-existing data
    load(fname);
catch           % if data doesn't exist, run simulations
    % % load ELFO orbit 
    % cspice_furnsh('data/gmat-to-spk/ELFO_20240501-20240531.bsp');
    % 
    % ELFO = '-909';                  % SPICE ID of ELFO s/c
    t0 = cspice_str2et(START);
    
    % planetary info
    bods = getplanets("MOON", "EARTH", "SUN", "JUPITER");
    
    % store spherical harmonic coefficients for the moon from Lunar Prospector
    [R,C,S] = cofloader("data/LP165P.cof");
    bods(1).C = C; bods(1).S = S;
    bods(1).frame = 'MOON_ME';      % body-fixed frame of coefficients
    moon = bods(1);                 % primary body
    sec = bods(2:end);              % secondary bodies
    
    % construct propagator
    % simulation parameters and execution
    ts = t0:60:t0 + DAYS*86400;
    m = length(ts);
    opts = odeset("RelTol", 1e-10, "AbsTol", 1e-12);
    load("data/optimization/RUN28- (56p,1).mat");
    % load("data/optimization/RUN16 (10p,0).mat");
    
    % assembling dynamics
    n = 6;
    f = @(t,x) orbitaldynamics(t, x, moon, 100, bods(2:end));
    dfdx = @(t,x) orbitalpartials(t, x, moon, bods(2));
    Q = zeros(6, 6);
    xopt(18 + OE) = F * pi/180;
    x0s = gaoutputtostates(xopt, 6, t0);
    x0 = x0s(:,OE);
    traj = statetotrajectory(t0, ts(end), x0, f);
    
    p = @(t,v) lyapunov(v, dfdx(t,traj(t)), Q);
    P0 = diag([1.29e-3 1.29e-3 1.29e-3 1.15e-7 1.15e-7 1.15e-7].^2);
    
    % true state propagation (aligns mostly with GMAT but off by a few meters, disagrees in spherical harmonics)
    [~,X] = ode45(f, ts, x0, opts);
    x_true = X';
    
    % covariance propagation
    P0_ = reshape(P0, n*n, 1);
    [~,X] = ode45(p, ts, P0_, opts);
    
    Phist = zeros(n, n, m);
    for i=1:m
        Phist(:,:,i) = reshape(X(i,:), n, n);
    end
    
    % propagate orbits for Monte-Carlo methods
    x_err = zeros(n, m, N_ITER);
    rerr = zeros(N_ITER, m);
    verr = zeros(N_ITER, m);
    
    for i=1:N_ITER
        x0_ = mvnrnd(x0, P0)';
        [~,X] = ode45(f, ts, x0_, opts);
        x_err(:,:,i) = X' - x_true;
        rerr(i,:) = sqrt(sum(x_err(1:3,:,i).^2, 1));
        verr(i,:) = sqrt(sum(x_err(4:6,:,i).^2, 1));
    end

    save(fname);
end

%% Actual statistics computation
rvar_comp = sum(rerr.^2)/(N_ITER-1);
vvar_comp = sum(verr.^2)/(N_ITER-1);
rvar = reshape(Phist(1,1,:)+Phist(2,2,:)+Phist(3,3,:), 1, m);
vvar = reshape(Phist(4,4,:)+Phist(5,5,:)+Phist(6,6,:), 1, m);
Rsqr = 1 - sum((rvar_comp - rvar).^2)/sum((rvar_comp - mean(rvar_comp)).^2);
Rsqv = 1 - sum((vvar_comp - vvar).^2)/sum((vvar_comp - mean(vvar_comp)).^2);
rstd = sqrt(rvar); vstd = sqrt(vvar);
rstd_comp = sqrt(rvar_comp); vstd_comp = sqrt(vvar_comp);

plotformat("IEEE", 1/2, "scaling", 2);
figure();
subplot(1,2,1);
scatter(rstd_comp, rstd);
hold on;
rline = [0 max(rstd)*1.1];
plot(rline, rline, '--');
hold off; grid on;
axis([rline rline]);
xlabel("Experimental \sigma^2"); ylabel("Computed \sigma^2");
title(sprintf("Position (R^2 = %.2f)", Rsqr));

subplot(1,2,2);
scatter(vstd_comp, vstd);
hold on;
vline = [0 max(vstd)*1.1];
plot(vline, vline, '--');
hold off; grid on;
axis([vline vline]);
xlabel("Experimental \sigma^2"); ylabel("Computed \sigma^2");
title(sprintf("Velocity (R^2 = %.2f)", Rsqv));
sgtitle("Fit of computed vs. experimental variance");

tplot = (ts - t0) / 3600;
figure();
plot(tplot, rstd_comp - rstd);

%% plot errors
plotformat("IEEE", 1, "scaling", 2);
h4 = figure();
colors = colororder();

% plot radial position error over time with 3-sigma bound
subplot(2,1,1);
patch([ tplot(1) tplot tplot(end)]', [0 rstd * 3e3 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 rstd_comp * 3e3 0]', colors(4,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);

for i=1:N_ITER      % compute error for individual run
    plot(tplot, rerr(i,:) * 1e3, "Color", colors(1,:), "LineWidth", 1);
end

hold off; box off;
axis([0 DAYS*24 -inf inf]);
ylabel("Position (m)");
title("State error over time");
legend(["Experimental 3\sigma","Computed 3\sigma","MC Run"],"location","northwest");

% plot radial velocity error over time with 3-sigma bound
subplot(2,1,2);
patch([ tplot(1) tplot tplot(end)]', [0 vstd * 3e6 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 vstd_comp * 3e6 0]', colors(4,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);

for i=1:N_ITER      % compute error for individual run
    plot(tplot, verr(i,:) * 1e6, "Color", colors(1,:), "LineWidth", 1);
end

hold off; box off;
axis([0 DAYS*24 -inf inf]);
ylabel("Velocity (mm/s)");
xlabel("Time (hrs)");
