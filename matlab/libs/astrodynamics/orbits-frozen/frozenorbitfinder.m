function [e,a] = frozenorbitfinder(i, drift)
%FROZENORBITFINDER Given specific target orbital elements, returns a
%complete set of elements for a lunar frozen orbit.
%
%   Frozen orbit solutions are discussed in [1]. To define one appropriately,
%   argument of periapsis must be 90 or 270 deg and inclination must be between 
%   39.23 and 140.77 degrees. Elements are given in the MOON_OP SPICE frame.
%
%   Input:
%    - i; inclination, in radians
%    - drift; optional, drift rate of right ascension of ascending node
%   Output:
%    - e; corresponding orbit eccentricity
%    - a; semimajor axis if drift is provided
arguments
    i (1,1) double
    drift (1,1) double {mustBeNonpositive} = 0
end

% parse and check inputs for validity
if nargin == 1
    drift = 0;
end
if i < 39.23*pi/180 || i > 140.77*pi/180
    error("frozenorbitfinder:iOutOfBounds", ...
        "Inclination must be between 39.23 and 140.77 degrees.")
end

% compute stability conditions; start with fixed inclination, as that is
% most important for reducing GDOP. Then, inclination drives eccentricity
% to set periapsis drift to zero. Finally, if the user desires a specific
% drift rate of the ascending node, a semimajor axis is provided.
e = sqrt(1 - 5/3*cos(i)^2);

G = 6.6743015e-20;      % kg*km^3/s^2, universal gravitational constant
m_M = 0.07346e24;       % kg, mass of moon
m_E = 5.9724e24;        % kg, mass of Earth
a_E = 0.3844e6;         % km, semimajor axis of orbit b/n two
n3sq = (G*(m_M + m_E)) / a_E^3;
a = (G*m_M * (4/3*drift/n3sq*sqrt(1-e^2)/((4*e^2+1)*cos(i)))^2)^(1/3);
end

