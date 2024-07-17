function [outputArg1,outputArg2] = surfacetrajgen(ts,loc,npts,vavg)
%SURFACETRAJGEN Generates a trajectory (position, velocity, acceleration)
%consisting of travel between waypoints for a surface user, given an operating 
%region.
%   Input:
%    - ts; start time, finish time -- in seconds past J2000
%    - loc; latitude, longitude, radius of operations -- in rad, rad, km
%           respectively
%    - npts; number of waypoints
%    - vavg; average travel speed between waypoints
arguments
    ts      (1,2)   double {mustBePositive}
    loc     (1,3)   double
    npts    (1,1)   {mustBeInteger,mustBePositive}
    vavg    (1,1)   double {mustBePositive}
end

% get center point in MOON_ME frame and convert to stereographic
ctr = [sin(loc(1))*cos(loc(2)) sin(loc(1))*sin(loc(2)) cos(loc(1))]';
ctr_st = ctr(1:2) / (1 - ctr(3));   % normalized stereographic
% get LHS of points in radius of operations
pts = lhsdesign(npts, 2);
pts = pts * diag([loc(3) 2*pi]);
pts = pts(:,1) .* [cos(pts(:,2)) sin(pts(:,2))]; 

data = importdata("data/LDEM_80S_80MPP_ADJ.TIF");
dmin = min(data,[],"all");
dmax = max(data,[],"all");
range = 304;            % km, stereographic range
mpp = 80;               % m, meters per pixel
kmpp = mpp * 1e-3;
range = range - kmpp/2;
px2km = -range:kmpp:range;
end

