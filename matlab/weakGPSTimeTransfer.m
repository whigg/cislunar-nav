% weakGPSTimeTransfer.m
% Author: Mark Hartigan
% Date  : November 8, 2023
% Description:
%    Simulate time transfer to lunar satellites.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% var init
max_a = 40 * pi/180;    % rad, maximum off-boresight angle link can close
c = 299792458;          % m/s, speed of light

% SPICE kernel initialization
t0 = convertTo(datetime('11-Feb-2027 00:00:00'), 'juliandate');
t0 = (t0 - 2451545) * 86400;    % JD to seconds past J2000

imdata = {'data/NSNS-DRM1.0-Khon2.bsp', 'data/NSNS-DRM1.0-Khon2.tpc.pc'};
gendata = {'data/naif0012.tls', 'data/de430.bsp', 'data/gm_de431.tpc'};
cspice_furnsh([imdata, gendata]);

% get Moon, Earth, Sun, and GPS satellite trajectory handles
file = 'G_2023-10-27_to_2023-11-02_060.pos';
[t, moon, earth, sun, sats, satdata] = load_trajectories(file);
tdays = (t - t(1)) ./ (60 * 60 * 24);   % elapsed days
n = length(t);          % number of time steps
M = length(sats);       % number of satellites total

%% Khon-2 trajectory
% realign time to start of Khon timespan for gathering of truth data
x_true = cspice_spkezr('Khon-2', (t-t(1) + t0)', 'J2000', 'NONE', 'MOON');
R = rotz(-30 * pi/180);
x_true(1:3,:) = R * x_true(1:3,:);
x_true(4:6,:) = R * x_true(4:6,:);

%% Precomputed clock trajectory
[Xclk, Qclk] = clockStateOverTime(t - t(1), 'CSAC');
Xclk = Xclk .* c;                       % convert from time to distance

%% pseudorange calculations
% compute measured and expected pseudoranges
% [UERE; UERRE] with an advanced receiver (Table 3.4-1, 3.4-2, A.4-1, and 
% B.2-1; SPS PS 2020) plus clock error
R_GPS = diag([9.7/1.96; .006/1.96].^2);
UERE = @(i) mvnrnd([0 0], R_GPS*(1.^2))' + Xclk(1:2,i);
% IM expected OD error ~10m, ~1cm/s (?) 3-sigma
R_IM  = diag([4.0; .001/3].^2);
EPHE = @(~) mvnrnd([0 0], R_IM)';

[psr_meas, psrr_meas] = compute_psr(t, x_true, sats, earth, moon, max_a, UERE);
[psr_comp, psrr_comp] = compute_psr(t, x_true, sats, earth, moon, max_a, EPHE);
psr_res = (psr_meas - psr_comp);        % m, pseudorange residuals
psrr_res = (psrr_meas - psrr_comp);     % m/s, ^range-rate residuals

figure();
plot(tdays, Xclk(1,:));
hold on;
plot(tdays, sum(psr_res, 1) ./ sum(psr_res ~= 0, 1));
hold off; grid on;
xlabel("Time (days)"); ylabel("Clock bias (m)");

%% Kalman filter design
%  clock bias and frequency offset is assumed to be zero at start of
%  simulation

% state vector estimate (phase deviation, frequency deviation, frequency drift)
X = zeros(3, n);
X(3,1) = 9e-10 / (86400 * 30) * c;  % m*Hz/Hz, upper end of aging rate a

% covariance matrix over time
Pt = zeros(3, 3, n);

% state transition matrix, tau = t_k - t_k-1
STM = @(tau) [1 tau tau^2/2;
              0 1   tau    ;
              0 0   1       ];

% covariance matrix of error associated w/ Wiener processes
[s1, s2, s3] = DiffCoeffCSAC();
s1 = s1 + 1e-10;
s2 = s2 + 5e-14;
Q = @(tau) ...
     [s1^2*tau + s2^2/3*tau^3 + s3^2/20*tau^5 s2^2/2*tau^2 + s3^2/8*tau^4 s3^2/6*tau^3;
      s2^2/2*tau^2 + s3^2/8^tau^4             s2^2*tau + s3^2/3*tau^3     s3^2/2*tau^2;
      s3^2/6*tau^3                            s3^2/2*tau^2                s3^2*tau     ];

% measurement variances of residuals (GPS + IM)
var_res = diag(R_GPS + R_IM);

for i=2:n
    dt = t(i) - t(i-1);                 % time step
    A = STM(dt);                        % state transition matrix
    Xest = A * X(:,i-1);                % a-priori state estimate
    Pest = A*Pt(:,:,i-1)*A' + Q(dt).*c^2;   % a-priori estimate covariance
    nvis = sum(psr_comp(:,i) ~= 0);     % visible satellites
    % measurements [range (nvis x 1); range-rate (nvis x 1)]
    z = [psr_res(psr_res(:,i) ~= 0, i); psrr_res(psrr_res(:,i) ~= 0, i)];

    % measurement model
    C = [repmat([1 0 0], nvis, 1); repmat([0 1 0], nvis, 1)];
    y = z - C * Xest;                   % measurement pre-fit residual
    % measurement covariance matrix
    R = [eye(nvis) * var_res(1) zeros(nvis, nvis)     ;
         zeros(nvis, nvis)      eye(nvis) * var_res(2)];
    S = C * Pest * C' + R;              % pre-fit residual covariance
    K = Pest * C' / S;                  % optimal Kalman gain
    X(:,i) = Xest + K * y;              % a-posteriori state estimate
    Pt(:,:,i) = (eye(3) - K*C) * Pest;  % a-posteriori estimate covariance
end

%% Plots and statistics
clc, close all;                     % reset since this block is run a lot

% statistics of clock bias
diff = Xclk(1,:) - X(1,:);
fprintf("3-sigma bound of Monte-Carlo run: %.3f m\n", std(diff) * 3);
var_bias = reshape(Pt(1,1,:), [], 1);
bnds = sqrt(var_bias) * 3;
fprintf("Median 3-sigma bound of filter: %.3f m\n", median(bnds));
nview = sum(psr_res ~= 0);

% statistics of clock drift rate
diff2 = Xclk(2,:) - X(2,:);
fprintf("3-sigma bound of Monte-Carlo run: %.3f mm/s\n", std(diff2) * 3e3);
var_drift = reshape(Pt(2,2,:), [], 1);
bnds2 = sqrt(var_drift) * 3e3;
fprintf("Median 3-sigma bound of filter: %.3f mm/s\n", median(bnds2));



% Plot monte-carlo run and performance of clock bias error
figure();
plot(tdays, Xclk(1,:) - X(1,:), 'LineWidth', 1);
hold on;
plot(tdays,  sqrt(var_bias) * 3, 'g--', 'LineWidth', 1);
plot(tdays, -sqrt(var_bias) * 3, 'g--', 'LineWidth', 1, ...
     'HandleVisibility', 'off');
plotViewOutages(tdays, nview, ylim);
% plot(tdays, X(1,:));
% plot(tdays, sum(psr_res) ./ sum(psr_res ~= 0));
hold off; grid on;
xlabel("Time (days)"); ylabel("Clock bias (m)");
title("Onboard clock bias error from Kalman filter");
legend(["MC Run", "3\sigma bound", "View outage"], 'location', 'best');
ylim([-20 20]);

% Plot satellites in view each day
figure();
plot(tdays, nview, 'LineWidth', 1);
grid on;
xlabel("Time (days)"); ylabel("Visible GPS satellites");
title("Number of visible GPS satellites over time");

%% render scene
figure();

plot3(x_true(1,:), x_true(2,:), x_true(3,:), 'c', 'Linewidth', 2);
hold on;

ti = 1;
x_ef = earth.x(t(ti)); x_ef = x_ef(1:3);
x_sf = sun.x(t(ti)); x_sf = x_sf(1:3);

% show GPS satellites
for k=1:M
    % visible GPS sats
    if psr_meas(k,ti) ~= 0
        pt = sats{k}(t(ti));
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
[Imoon, ~] = imread("data/lroc_color_poles_1k.jpg");
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

%% plot all GNSS satellites and affiliation
figure();
subset = [];

% show GPS satellites
for i=1:M
    % GPS trajectories
    if i == 1
        p = plot3(satdata(i).x(1,:) + x_ef(1), satdata(i).x(2,:) + x_ef(2), ...
              satdata(i).x(3,:) + x_ef(3), 'Color', 'g', 'LineWidth', ...
              1, 'DisplayName', 'GPS');
        subset = [subset p];
        hold on;
    else
        plot3(satdata(i).x(1,:) + x_ef(1), satdata(i).x(2,:) + x_ef(2), ...
              satdata(i).x(3,:) + x_ef(3), 'Color', 'g', 'LineWidth', 1);
    end
end

% earth
[xx, yy, zz] = ellipsoid(x_ef(1), x_ef(2), x_ef(3), earth.R, earth.R, earth.R);
globe = surf(xx, yy, zz);
set(globe, 'FaceColor', 'texturemap', 'CData', flip(Iearth,1), 'FaceAlpha', 1, ...
    'EdgeColor', 'none');
hold off; axis equal;
xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");
% legend(subset);
title('24-hour flight paths of GPS satellites, 2023-10-27');

%% cleanup
cspice_unload([imdata, gendata]);
