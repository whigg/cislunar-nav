function [LAT,LON] = LOLAxy2latlon(X,Y)
%LOLAXY2LATLON Converts X/Y coordinates in stereographic lunar south pole
%frame used by LOLA to latitude and longitude.
%   Input:
%    - X; x coordinate in km (+Y in MOON_ME)
%    - Y; y coordinate in km (+X in MOON_ME)
%   Output:
%    - LAT; lunar latitude in rad (-pi/2 to pi/2)
%    - LON; lunar longitude in rad (0 to 2pi)
arguments
    X (1,1) double
    Y (1,1) double
end

R = sqrt(X^2 + Y^2);
LAT = -pi/2 + 2*atan(0.5 * R/1737.4);  % southern hemisphere
LON = atan2(X,Y);
end

