% weakGPSTimeTransfer.m
% Author: Mark Hartigan
% Date  : November 8, 2023
% Description:
%    Simulate time transfer to lunar satellites.

%% reset
clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% var init
max_a = 40 * pi/180;    % rad, maximum off-boresight angle link can close
c = 299792.458;         % km/s, speed of light

% SPICE kernel initialization
t0 = convertTo(datetime('2-Feb-2027 00:00:00'), 'juliandate');
t0 = (t0 - 2451545) * 86400;    % JD to seconds past J2000

imdata = {'data/NSNS-DRM1.0-Khon2.bsp', 'data/NSNS-DRM1.0-Khon2.tpc.pc'};
gendata = {'data/naif0012.tls', 'data/de430.bsp', 'data/gm_de431.tpc'};
cspice_furnsh([imdata, gendata]);

% get Moon, Earth, Sun, and GPS satellite trajectory handles
file = 'G_2023-10-27_300_oc2_060.pos';
[t, moon, earth, sun, sats, satdata] = load_trajectories(file);
% t = t - t(1) + t0;      % realign time to start of Khon timespan
n = length(t);          % number of time steps
M = length(sats);       % number of satellites total

%% Khon-2 trajectory
% realign time to start of Khon timespan for gathering of truth data
x_true = cspice_spkezr('Khon-2', (t-t(1) + t0)', 'J2000', 'NONE', 'MOON');

% Khon 2 nominal trajectory
x0_ = [x_true(:,1); 0];
P0 = [1 1 1 .001 .001 .001 1]';
x0 = random('Normal', x0_, P0, [7,1]); % initial guess
P0 = diag(P0.^2);

%% Precomputed clock trajectory
[Xclk, Qclk] = clockStateOverTime(t - t(1), 'CSAC');

%% trajectory and pseudorange calculations
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
F = @(t,x) moondyn(t, x, moon, earth, sun);
[~, x_nom] = ode45(F, t, x0, opts);   % nominal trajectory
x_nom = x_nom';   % keep with convention of cols = time steps

% compute measured and expected pseudoranges
% UERE with an advanced receiver (Table A.4-1 + B.2-1, SPS PS 2020) plus
% clock error
UERE = @(i) random('Normal',0,0.009736/1.96,1) + Xclk(1,i) * c;
% IM expected OD error ~10m 3-sigma
EPHE = @(~) random('Normal',0,0.01/3,1);
psr_meas = compute_psr(t, x_true, sats, earth, max_a, UERE);
psr_comp = compute_psr(t, x_true, sats, earth, max_a, EPHE);
psr_res = psr_meas - psr_comp;

figure();
plot(t, Xclk(1,:)*c);
hold on;
plot(t, max(psr_res));
hold off; grid on;



%% render scene
figure();

plot3(x_true(1,:), x_true(2,:), x_true(3,:), 'c', 'Linewidth', 2);
hold on;

x_ef = earth.x(t(177)); x_ef = x_ef(1:3);
x_sf = sun.x(t(177)); x_sf = x_sf(1:3);

% show GPS satellites
for k=1:M
    % visible GPS sats
    if psr(k,177) ~= 0
        pt = sats{k}(t(177));
        scatter3(pt(1), pt(2), pt(3), 'go', 'Linewidth', 2);
    end
end

% earth

[Iearth, ~] = imread("data/flat_earth_Largest_still.0330.jpg");
[xx, yy, zz] = ellipsoid(x_ef(1), x_ef(2), x_ef(3), earth.R, earth.R, earth.R);
globe = surf(xx, yy, zz);
set(globe, 'FaceColor', 'texturemap', 'CData', flip(Iearth,1), 'FaceAlpha', 1, ...
    'EdgeColor', 'none');

% moon
[xx, yy, zz] = ellipsoid(0, 0, 0, moon.R, moon.R, moon.R);
globe = surf(xx, yy, -zz);
set(globe, 'FaceColor', 'texturemap', 'CData', Imoon, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');

hold off; grid on; axis equal;
% set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");
legend(["Spacecraft Trajectory", "GNSS Satellite"], 'location', 'best', ...
    'TextColor', 'w');

% %% plot all GNSS satellites and affiliation
% figure();
% subset = [];
% 
% % show GPS satellites
% for i=1:M
%     % GPS trajectories
%     if i == 1
%         p = plot3(satdata(i).x(1,:) + x_ef(1), satdata(i).x(2,:) + x_ef(2), ...
%               satdata(i).x(3,:) + x_ef(3), 'Color', 'g', 'LineWidth', ...
%               1, 'DisplayName', 'GPS');
%         subset = [subset p];
%         hold on;
%     else
%         plot3(satdata(i).x(1,:) + x_ef(1), satdata(i).x(2,:) + x_ef(2), ...
%               satdata(i).x(3,:) + x_ef(3), 'Color', 'g', 'LineWidth', 1);
%     end
% end
% 
% % earth
% [xx, yy, zz] = ellipsoid(x_ef(1), x_ef(2), x_ef(3), earth.R, earth.R, earth.R);
% globe = surf(xx, yy, zz);
% set(globe, 'FaceColor', 'texturemap', 'CData', flip(Iearth,1), 'FaceAlpha', 1, ...
%     'EdgeColor', 'none');
% hold off; axis equal;
% xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");
% % legend(subset);
% title('24-hour flight paths of GPS satellites, 2023-10-27');
% 
% %% cleanup
% cspice_unload(imdata, gendata);

%% functions

function r = G(x, sats)
%G Measurement model of GPS range. Returns vector of computed ranges.
%   Input:
%    - x; state vector (7,1)
%    - sats; vector of satellite positions + clock offset (n,4)

r = zeros(size(sats,2), 1);

for i=1:length(r), r(i) = sqrt(sum((x(1:3) - sats(:,i)).^2)) + x(7); end
end

function H = Ht(x, sats)
r = G(x, sats);      % range measurements
H = zeros(length(r), length(x));

for i=1:length(r)
    % p = r(i) - x(7);    % transform back to pseudorange
    p = r(i);
    H(i,:) = [(x(1)-sats(1,i))/p (x(2)-sats(2,i))/p (x(3)-sats(3,i))/p 0 0 0 1];
end
end