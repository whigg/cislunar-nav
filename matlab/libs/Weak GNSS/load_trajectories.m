function [t, moon, earth, sun, sats, satdata] = load_trajectories(file)
%LOAD_TRAJECTORIES Creates callable functions for all the necessary
%planetary and GPS satellite trajectories over the time span specified in
%the provided SP3 or POS file.
%   Input:
%    - file; relative path to the .sp3 or .pos file

moon.R  = 1734;         % km, equatorial radius of moon
earth.R = 6378;         % km, equatorial radius of earth
sun.R   = 695700;       % km, radius of sun

moon.GM  = cspice_bodvrd('MOON','GM',1);
earth.GM = cspice_bodvrd('EARTH','GM',1);
sun.GM   = cspice_bodvrd('SUN','GM',1);

% SPICE planetary data
% [km; km/s] position and velocity of Earth/Sun w.r.t. the moon
earth.x = @(t) cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'MOON');
sun.x   = @(t) cspice_spkezr('SUN'  , t, 'J2000', 'NONE', 'MOON');
moon.x  = @(t) zeros(3,1);

% GPS data
if endsWith(file, '.sp3') || endsWith(file, '.SP3')
    [t, satdata] = read_sp3(file);      % get GPS sp3 data for 2023-04-25
    t = (t - 2451545) * 86400;          % JD to seconds past J2000
elseif endsWith(file, '.pos')
    [t, satdata] = read_pos(file);      % get GPS pos data from JPL GDGPS
    % time already in seconds past j2000
else
    throw(MException('load_trajectories:fileParsing', ...
        'Bad file extension: only .SP3 and .pos supported.'));
end

% figure();
% plot3(satdata(1).x(1,:), satdata(1).x(2,:), satdata(1).x(3,:));
% hold on;
% % earth
% [Iearth, ~] = imread("data/flat_earth_Largest_still.0330.jpg");
% [xx, yy, zz] = ellipsoid(0, 0, 0, earth.R, earth.R, earth.R);
% globe = surf(xx, yy, zz);
% set(globe, 'FaceColor', 'texturemap', 'CData', flip(Iearth,1), 'FaceAlpha', 1, ...
%     'EdgeColor', 'none');
% hold off; axis equal;
% xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");

% generate interpolatory functions for all GPS satellites
sats = cell(length(satdata), 1);
for i=1:length(sats)
    pps = spline(t, satdata(i).x(1:3,:));       % piecewise polynomial structure
    ppsv = spline(t, satdata(i).x(4:6,:));      % piecewise polynomial structure
    % spline interpolation func; add earth position vector at t
    sats{i} = @(t) [ppval(pps, t); ppval(ppsv, t)] + earth.x(t);
end
end

