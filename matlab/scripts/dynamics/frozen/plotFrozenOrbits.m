% plotFrozenOrbits.m
% Author: Mark Hartigan
% Date  : May 1, 2024
% Description:
%    Verify functionality of frozenorbitfinder.m by generating and plotting
%    them.

% %% reset
% clc, clear, close all;
% addpath(genpath(pwd));
% format long g;          % display long numbers, no scientific notation

%% SCRIPT INFORMATION
% INCLINATION, DEG
I = 45;
% SIMULATION START DATE AND TIME
START = '2024 May 1 00:00:00';
% DURATION OF SIMULATION
DAYS = 30;
% COLOR OF SPACECRAFT IN PLOT
COLOR = [0 0.4470 0.7410];
% REFERENCE FRAME FOR VIEWING OF ORBIT
PLOT_FRAME = 'MOON_OP';
% PERILUNE ALTITUDE OF SPACECRAFT
ALT = 300;

%% load SPICE kernels
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
t0 = cspice_str2et(START);
mu_m = cspice_bodvrd('MOON', 'GM', 1);

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

%% define orbit
% inclination must be between 39.23 and 140.77 degrees
i = I*pi/180;          % inclination in MOON_ME is i - 6.76 degrees
e = frozenorbitfinder(i);
a = 15000;
RAAN = 270*pi/180;
w = 90*pi/180;
f = 0;
drift = ascendingnodedrift(a,e,i);
[r0,v0] = oe2rv(a,e,i,RAAN,w,f,mu_m);
% z_OP = cspice_pxform('MOON_ME', 'MOON_OP', t0+1) * [0;0;1];
% atan2(norm(z_OP(1:2)),z_OP(3)) * 180/pi
oe1.i = i; oe1.e = e; oe1.a = a; oe1.RAAN = RAAN; oe1.w = w;
x0 = cspice_sxform('MOON_OP', 'J2000', t0) * [r0; v0];

%% simulate orbit
% simulation
ts = t0:60*DAYS:t0 + 86400*DAYS;
ns = length(ts);
opts = odeset("RelTol", 1e-8, "AbsTol", 1e-9);
[~,X] = ode45(@(t,x) orbitaldynamics(t,x,moon,100,earth,sun), ts, x0, opts);
x = X';

%% plot orbit
%  close all;

% transform data to plotting frame if necessary
if ~strcmp(PLOT_FRAME, 'J2000')
    for j=1:size(x,2)
        xp(:,j) = cspice_sxform('J2000', PLOT_FRAME, ts(j)) * x(:,j);
    end
else
    xp = x;
end

h2 = figure();
% Display moon in trajectory plot
R_me = 1738.1;          % km, moon equatorial radius
R_mp = 1736;            % km, moon polar radius
[Imoon, ~] = imread(userpath + "/astrodynamics/Moon_HermesCelestiaMotherlode.jpg");
[xx, yy, zz] = ellipsoid(0, 0, 0, R_me, R_me, R_mp);

% Rotate moon from MOON_ME frame to PLOT_FRAME
T = cspice_pxform('MOON_ME', PLOT_FRAME, ts(end));
for j=1:size(xx,1)
    for k=1:size(xx,2)
        % -z to flip image
        tmp = T * [xx(j,k); yy(j,k); -zz(j,k)];
        xx(j,k) = tmp(1); yy(j,k) = tmp(2); zz(j,k) = tmp(3);
    end
end

globe = surf(xx, yy, zz);
set(globe, 'FaceColor', 'texturemap', 'CData', Imoon, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');
hold on;

% put a little star on the south pole
lsp = T * [0;0;-R_mp];
scatter3(lsp(1), lsp(2), lsp(3), 100, "red", "filled", "pentagram");

% plot user trajectory for same time frame
plot3(xp(1,:), xp(2,:), xp(3,:), "LineWidth", 1.5, "Color", COLOR, ...
    "Marker", "diamond", "MarkerFaceColor", COLOR, "MarkerIndices", ns);

grid on; axis equal;
SUB = strsplit(PLOT_FRAME,"_");
SUB = SUB(end);
xlabel("x_{"+SUB+"} (km)");
ylabel("y_{"+SUB+"} (km)");
zlabel("z_{"+SUB+"} (km)");
title(int2str(DAYS) + "-day propagation, starting at " + START);

%% verify frozen-ness
xf = cspice_sxform(PLOT_FRAME, 'MOON_OP', ts(end)) * x(:,end);
[a2,e2,i2,RAAN2,w2,f2] = rv2oe(xf(1:3), xf(4:6), mu_m);
expdrift = (RAAN2-RAAN) / (ts(end)-ts(1));
fprintf("Computed drift of AN: %.3e\n", drift);
fprintf("Experimental   ''   : %.3e\n", expdrift);