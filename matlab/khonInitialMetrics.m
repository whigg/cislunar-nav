%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
t0 = convertTo(datetime('9-Feb-2027 00:00:00'), 'juliandate');
t0 = (t0 - 2451545) * 86400;    % JD to seconds past J2000

imdata = cell(1, 12);
for i=0:5
    imdata{2*i + 1} = sprintf('data/NSNS-DRM1.0-Khon%d.bsp', i+1);
    imdata{2*i + 2} = sprintf('data/NSNS-DRM1.0-Khon%d.tpc.pc', i+1);
end

gendata = {'data/naif0012.tls', 'data/de430.bsp', 'data/gm_de431.tpc'};
cspice_furnsh([imdata, gendata]);

moon.R  = 1734;         % km, equatorial radius of moon
moon.GM = cspice_bodvrd('MOON','GM',1);
moon.x  = @(t) zeros(3,1);

getpos = @(x) x(1:3,:);
Khon1 = @(t) getpos(cspice_spkezr('Khon-1', t, 'J2000', 'NONE', 'MOON'));
Khon2 = @(t) getpos(cspice_spkezr('Khon-2', t, 'J2000', 'NONE', 'MOON'));
Khon3 = @(t) getpos(cspice_spkezr('Khon-3', t, 'J2000', 'NONE', 'MOON'));
Khon4 = @(t) getpos(cspice_spkezr('Khon-4', t, 'J2000', 'NONE', 'MOON'));
Khon5 = @(t) getpos(cspice_spkezr('Khon-5', t, 'J2000', 'NONE', 'MOON'));
Khon6 = @(t) getpos(cspice_spkezr('Khon-6', t, 'J2000', 'NONE', 'MOON'));

%% build trajectories
days = 14;
step = 60*14;    % seconds
t = t0:step:t0 + 86400 * days;
xk1 = Khon1(t); xk2 = Khon2(t); xk3 = Khon3(t);
xk4 = Khon4(t); xk5 = Khon5(t); xk6 = Khon6(t);
xk = [xk1; xk2; xk3; xk4; xk5; xk6];
% xk = [xk1; xk2; xk3; xk4];

%% plot relative orbits
figure();

% plot orbits
for i = 0:size(xk,1)/3-1
    plot3(xk(3*i + 1,:), xk(3*i + 2,:), xk(3*i + 3,:), 'LineWidth', 1);
    if i == 0, hold on; end
end

% moon
[Imoon, ~] = imread("data/lroc_color_poles_1k.jpg");
[xx, yy, zz] = ellipsoid(0, 0, 0, moon.R, moon.R, moon.R);
globe = surf(xx, yy, zz);
set(globe, 'FaceColor', 'texturemap', 'CData', Imoon, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');

hold off; grid on; axis equal;
% set(gcf, 'Color', 'k');
% set(gca, 'Color', 'k');
xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");

%% write data
data = [t' xk'];

% fid = fopen('Khon_4sat_1day.txt', 'w');
% for i=1:size(data,1)
%     fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', data(i,2:end));
% end
% fclose(fid);

%% metrics
UNE(data);

constellation_metrics(data);
