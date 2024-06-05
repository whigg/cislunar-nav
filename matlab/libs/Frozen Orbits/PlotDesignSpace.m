% PlotDesignSpace.m
% Author: Mark Hartigan
% Date  : May 14, 2024
% Description:
%    Plots the optimization design space of lunar frozen orbits as a 3D
%    contour plot of inclination/eccentricity, semimajor axis, and
%    orbit precession rate.

% %% reset
% clc, clear, close all;
% addpath(genpath(pwd));
% format long g;          % display long numbers, no scientific notation

%% variables
MINALT = 500;           % km, minimum altitude of spacecraft
STABILITY = true;
START = '2024 May 1 00:00:00';
TOL = 0.1;
res = 2;
is = (39.24:res:140.77) * pi/180;
es = arrayfun(@frozenorbitfinder, is);
as = 1738:res*300:35000;
[IS,AS] = meshgrid(is,as);
[ES,~] = meshgrid(es,as);
RAANS = zeros(size(IS));

%% evaluation
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
t0 = cspice_str2et(START);

% earth information
earth.GM = cspice_bodvrd('EARTH', 'GM', 1);
earth.x = @(tau) cspice_spkpos('EARTH', tau, 'J2000', 'NONE', 'MOON');

% sun information
sun.GM = cspice_bodvrd('SUN', 'GM', 1);
sun.x = @(tau) cspice_spkpos('SUN', tau, 'J2000', 'NONE', 'MOON');

% moon information
moon.GM = cspice_bodvrd('MOON', 'GM', 1);
[R,C,S] = cofloader(userpath + "/astrodynamics/LP165P.cof");
moon.x = @(~) [0;0;0];
moon.R = R * 1e-3;      % convert from m to km
moon.C = C;             % store in moon struct for orbitaldynamics
moon.S = S;             % store in moon struct for orbitaldynamics
moon.frame = 'MOON_ME'; % body-fixed frame of coefficients

opts = odeset("RelTol", 1e-8, "AbsTol", 1e-9);

for i=1:numel(RAANS)
    RAANS(i) = ascendingnodedrift(AS(i),ES(i),IS(i));
    % periapsis radius too small
    if AS(i) < (1738 + MINALT) / (1 - ES(i))
        RAANS(i) = NaN;
        AS(i) = NaN;
        ES(i) = NaN;
        IS(i) = NaN;
    elseif STABILITY    % now check stability
        T = 2*pi*sqrt(AS(i)^3/moon.GM);
        [r0,v0] = oe2rv(AS(i), ES(i), IS(i), 0, pi/4, 0, moon.GM);
        x0 = cspice_sxform('MOON_OP', 'J2000', t0) * [r0; v0];
        [~,X] = ode45(@(t,x) orbitaldynamics(t,x,moon,2,earth,sun), [t0 t0+T], x0, opts);
        xf = cspice_sxform('J2000', 'MOON_OP', t0+T) * X(end,:)';
        [ai, ei, ii, ~, ~, ~] = rv2oe(xf(1:3), xf(4:6), moon.GM);
    
        % after 1 orbital period, semimajor axis, eccentricity, or inclination
        % has changed by more than TOL%
        if abs(ai-AS(i))/AS(i) > TOL || abs(ei-ES(i)) > TOL / 3
            RAANS(i) = NaN;
            AS(i) = NaN;
            ES(i) = NaN;
            IS(i) = NaN;
        end
    end
end

IS = IS * 180/pi;
RAANS = RAANS * 180/pi * 86400;

%% plot 1
h1 = figure();
surf(IS,AS,RAANS,"EdgeColor", "none");
colormap(h1,parula)
c = colorbar;
c.Label.String = "Precession rate (deg/day)";
axis([is(1)*180/pi is(end)*180/pi as(1) as(end) -inf inf]);
xlabel("Inclination (deg)");
ylabel("Semimajor axis (km)");
zlabel("Precession rate (deg/day)");
title("Frozen orbit design space");

% %% plot 2
% close(h1);
% h2 = figure();
% surf(ES,AS,RAANS,"EdgeColor", "none");
% colormap(h2,parula)
% c = colorbar;
% axis([min(es) max(es) as(1) as(end) -inf inf]);
% xlabel("Eccentricity");
% ylabel("Semimajor axis (km)");
% zlabel("Precession rate (deg/day)");