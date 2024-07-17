function [outputArg1,outputArg2] = surfacetrajgen(bnds,vmax,dwell)
%SURFACETRAJGEN Generates a trajectory (position, velocity, acceleration)
%consisting of travel between waypoints for a surface user, given an operating 
%region.
%   Input:
%    - bnds; altitude from surface (min/max), latitude (min/max), longitude
%            (min/max) -- in (km, rad, rad) respectively
%    - vavg; average travel speed between waypointsl
%    - dwell; dwell time at waypoints (s)
arguments
    bnds
end

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

