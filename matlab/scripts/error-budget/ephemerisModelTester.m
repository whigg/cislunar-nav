% ephemerisModelTester.m
% Author: Mark Hartigan
% Date  : August 7, 2024
% Description:
%    Test polynomial vs Kepler's equation + polynomial model fits

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
% model to use, "Polynomial" or "Kepler"
models = ["Polynomial", "Kepler"];
N1 = [6,4];         % polynomial coefficients to use, low  (2 less for Kepler)
N2 = [14,12];       % polynomial coefficients to use, high (2 less for Kepler)

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
% initialize propagator (garbage starting state)
ephprop = LunarPropagator(t0, zeros(6,1), 50, 3);

% loop info
tas = (0:12:348) * pi/180;      % consider every f at 12-deg intervals
istep = 0.01;
ints = (istep:istep:1) * P;     % consider fractions of orbital period up to 1
% feasibility bounds
rbnd = 13.43/6;
vbnd = 1.2/6;
bars = zeros(length(tas), 4);

for m=1:2
    r1stds = Inf(length(tas), length(ints));
    r2stds = Inf(length(tas), length(ints));
    v1stds = Inf(length(tas), length(ints));
    v2stds = Inf(length(tas), length(ints));
    
    % get 3-sigma statistics for various starting TAs and intervals
    for j=1:length(tas)
        % set true anomaly of orbit
        [r0,v0] = oe2rv(a, e, i, RAAN, w, tas(j), GM);
        ephprop.x0 = cspice_sxform('MOON_OP', 'J2000', t0) * [r0; v0];
    
        for k=1:length(ints)
            tf = ints(k);           % final time of ephemeris propagation
            f_eph1 = ephprop.ephemerisfit(models(m), tf, N1(m));
            f_eph2 = ephprop.ephemerisfit(models(m), tf, N2(m));
    
            % truth orbit
            [ts,x_true] = ephprop.run(tf, 2000, 'J2000');
            % ephemeris approximation errors
            err1_x = f_eph1(ts) - x_true;
            err2_x = f_eph2(ts) - x_true;
    
            % compute statistics
            r1stds(j,k) = 3*std(sqrt(sum(err1_x(1:3,:).^2, 1)) * 1e3);
            r2stds(j,k) = 3*std(sqrt(sum(err2_x(1:3,:).^2, 1)) * 1e3);
            v1stds(j,k) = 3*std(sqrt(sum(err1_x(4:6,:).^2, 1)) * 1e6);
            v2stds(j,k) = 3*std(sqrt(sum(err2_x(4:6,:).^2, 1)) * 1e6);
    
            if (r1stds(j,k) > rbnd || v1stds(j,k) > vbnd) && ...
                    (r2stds(j,k) > rbnd || v2stds(j,k) > vbnd)
                break;
            end
        end
    end

    % compute and store feasibility
    bars(:,2*(m-1)+1) = sum((r1stds <= rbnd) & (v1stds <= vbnd), 2);
    bars(:,2*(m-1)+2) = sum((r2stds <= rbnd) & (v2stds <= vbnd), 2);
end

%% bar charts
bars2 = bars .* istep;

plotformat("IEEE", 0.4, "scaling", 2, "width", 16, "coloring", "reef");
tiledlayout(1,2, "TileSpacing", "tight", "Padding", "compact");
nexttile
bar(tas * 180/pi, bars2(:,1:2), 1);
xlabel("True Anomaly (\circ)");
ylabel("Max. Feasible Interval (norm.)");
title("Polynomial Model");
axis([-inf inf 0 0.8]);
legend([sprintf("N = %d", N1(1)), sprintf("N = %d", N2(1))], "location", "best");

nexttile
plotformat("IEEE", 0.4, "scaling", 2, "width", 16, "coloring", "dye");
bar(tas * 180/pi, bars2(:,3:4), 1);
xlabel("True Anomaly (\circ)");
title("Keplerian Model");
axis([-inf inf 0 0.8]);
legend([sprintf("N = %d", N1(2)), sprintf("N = %d", N2(2))], "location", "best");
set(gca, "YTickLabel", []);

