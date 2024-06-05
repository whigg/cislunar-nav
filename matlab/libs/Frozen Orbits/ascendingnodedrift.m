function drift = ascendingnodedrift(a,e,i)
%ASCENDINGNODEDRIFT Returns the drift of the ascending node for a given
%lunar frozen orbit.
%   Input:
%    - a; semimajor axis
%    - e; eccentricity
%    - i; inclination
arguments
    a (1,1) double {mustBePositive}
    e (1,1) double {mustBePositive}
    i (1,1) double {mustBePositive}
end

G = 6.6743015e-20;      % km^3/s^2/kg, universal gravitational constant
m_M = 0.07346e24;       % kg, mass of moon
m_E = 5.9724e24;        % kg, mass of Earth
a_E = 0.3844e6;         % km, semimajor axis of orbit b/n two

n3sq = (G*(m_M + m_E)) / a_E^3;
n = sqrt(G*m_M / a^3);

drift = 3/8 * n3sq/n * 1/sqrt(1-e^2) * cos(i) * (-8*e^2-2);
end

