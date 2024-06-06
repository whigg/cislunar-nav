%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% INPUT
% which GNSS satellites to include; options are "GPS", "GALILEO", or "BOTH"
INCLUDE_GNSS = "GPS";
% boolean, plot GNSS trajectories around Earth in their own figure?
PLOT_GNSS = true;

%% load GNSS trajectories
[fGnss, namesGnss] = getgnsshandles(INCLUDE_GNSS);
nGnss = length(namesGnss);

t0 = convertTo(datetime(START), "juliandate");
t0 = (t0 - 2451545) * 86400;    % JD to seconds past J2000
ts = t0:60:t0 + 86400*DAYS;     % simulation time steps
ns = length(ts);                % number of simulation time steps

sats = zeros(6,ns,nGnss);       % array of satellite states

for i=1:nGnss
    sats(:,:,i) = fGnss{i}(ts);
end

if PLOT_GNSS
    h1 = figure();
    % Display Earth in trajectory plot
    R_ee = 6378;                    % km, equatorial radius of Earth
    R_ep = 6357;                    % km, polar radius of Earth
    [Iearth, ~] = imread("res/ModifiedBlueMarble.jpg");
    [xx, yy, zz] = ellipsoid(0, 0, 0, R_ee, R_ee, R_ep);
    globe = surf(xx, yy, zz);
    set(globe, 'FaceColor', 'texturemap', 'CData', flip(Iearth,1), 'FaceAlpha', 1, ...
        'EdgeColor', 'none');
    
    hold on;
    
    for i=1:nGnss
        % assign colors to each trajectory based on host constellation
        if contains(namesGnss{i}, 'NAVSTAR')
            color = 'b';            % GPS is blue
        elseif contains(namesGnss{i}, 'GSAT')
            color = 'r';            % Galileo is red
        else
            color = 'k';            % unknown is black
        end
        plot3(sats(1,:,i), sats(2,:,i), sats(3,:,i), "LineWidth", 1.5, "Color", color);
    end
    
    hold off; grid on; axis equal;
    xlabel("x_{J2000} (km)", "FontSize", FSZ);
    ylabel("y_{J2000} (km)", "FontSize", FSZ);
    zlabel("z_{J2000} (km)", "FontSize", FSZ);
    title(int2str(DAYS) + "-day GNSS satellite propagation, starting at " + START, "FontSize", FSZ);
    fontname(h1, FONT);
end
