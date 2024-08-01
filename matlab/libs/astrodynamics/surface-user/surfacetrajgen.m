function [waypts,reftraj,exptraj] = surfacetrajgen(t0,loc,npts,vavg,dwell)
%SURFACETRAJGEN Generates a trajectory (position, velocity, acceleration)
%consisting of travel between waypoints for a surface user, given an operating 
%region.
%   Input:
%    - t0; start time in seconds past J2000
%    - loc; latitude, longitude, radius of operations -- in rad, rad, km
%           respectively
%    - npts; number of waypoints
%    - vavg; km/hr, average travel speed between waypoints
%    - dwell; dwell time, in seconds, at waypoints
arguments
    t0      (1,1)   double {mustBePositive}
    loc     (1,3)   double
    npts    (1,1)   {mustBeInteger,mustBePositive}
    vavg    (1,1)   double {mustBePositive}
    dwell   (1,1)   double {mustBePositive}
end

% initialize LOLAMap object
map = LOLAMap("data/LDEM_80S_80MPP_ADJ.TIF", "SOUTH");
% get center point in MOON_ME frame and convert to stereographic
[x_c, y_c] = map.latlon2xy(loc(1), loc(2));
ctr = [x_c y_c];
% get LHS of points in radius of operations
pts = lhsdesign(npts, 2) * diag([loc(3) 2*pi]);
pts = ctr + pts(:,1) .* [cos(pts(:,2)) sin(pts(:,2))];

% create order of waypoints
ordpts = repmat(ctr, npts+2, 1);
i = 1;
while size(pts,1) > 0
    % find closest point
    start = ordpts(i,:);                        % store starting point               
    dists = sqrt(sum((pts - start).^2, 2));     % get distances from start
    [~,j] = min(dists);                         % index of closest point

    % increment index and assign closest point to ordered list
    i = i + 1;
    ordpts(i,:) = pts(j,:);
    pts(j,:) = [];      % remove point from list
end

% create true and estimated piecewise polynomial trajectories
m = 6;                  % number of points to interpolate trajectory
ti = t0;                % track time
breaks = [];
refcoefs = [];
expcoefs = [];

for i=1:size(ordpts,1)-1    % for each interval
    % get reference points in stereographic frame
    x0 = ordpts(i,1);   y0 = ordpts(i,2);
    xf = ordpts(i+1,1); yf = ordpts(i+1,2);
    dist = norm(ordpts(i+1,:) - ordpts(i,:));
    ref = [linspace(x0, xf, m)' linspace(y0, yf, m)'];
    ts = linspace(ti, ti + 3600*dist/vavg, m);
    ti = ti + 3600*dist/vavg;

    % get true points in stereographic frame
    exp = ref;
    exp(2:end-1,:) = mvnrnd(ref(2:end-1,:), eye(2)*(dist/100)^2);
    ref_me = zeros(3,m);
    exp_me = ref_me;

    % convert to MOON_ME frame
    for j=1:m
        ref_me(:,j) = map.xy2me(ref(j,1), ref(j,2), true);
        exp_me(:,j) = map.xy2me(exp(j,1), exp(j,2), true);
    end

    % generate splines for movement and pause at waypoint
    ppi = spline(ts, [[0;0;0] ref_me [0;0;0]]);
    ppj = spline(ts, [[0;0;0] exp_me [0;0;0]]);
    ppause = spline([ti ti+dwell], [[0;0;0] ref_me(:,end) ref_me(:,end) [0;0;0]]);
    ti = ti + dwell;        % update time
    % concatenate list of breaks and coefficients
    breaks = [breaks ppi.breaks(1:end-1) ppause.breaks(1:end-1)];
    refcoefs = [refcoefs; ppi.coefs; ppause.coefs];
    expcoefs = [expcoefs; ppj.coefs; ppause.coefs];
end

% turn trajectory from solely position to pos, vel, and accel
fullref = zeros(size(refcoefs,1)*3,size(refcoefs,2));
fullexp = fullref;
dc = size(refcoefs, 2)-1:-1:1;
ddc = (size(refcoefs, 2)-2:-1:1) .* dc(1:end-1);
for i=1:size(refcoefs,1)/3
    fullref((i-1)*9+1:i*9,:) = [refcoefs((i-1)*3+1:i*3,:)
                                zeros(3,1) refcoefs((i-1)*3+1:i*3,1:end-1) .* dc
                                zeros(3,2) refcoefs((i-1)*3+1:i*3,1:end-2) .* ddc];
    fullexp((i-1)*9+1:i*9,:) = [expcoefs((i-1)*3+1:i*3,:)
                                zeros(3,1) expcoefs((i-1)*3+1:i*3,1:end-1) .* dc
                                zeros(3,2) expcoefs((i-1)*3+1:i*3,1:end-2) .* ddc];
end

breaks = [breaks ppause.breaks(end)];   % cap off break list
reftraj = mkpp(breaks, fullref, 9);
exptraj = mkpp(breaks, fullexp, 9);

% convert waypoints to MOON_ME frame
waypts = zeros(3,size(ordpts,1));
for j=1:size(waypts,2)
    waypts(:,j) = map.xy2me(ordpts(j,1), ordpts(j,2), true);
end

% figure();
% plot(ordpts(:,1), ordpts(:,2), "-o");
% hold on;
% plot(ctr(1), ctr(2), "r+");
% grid on;
end

