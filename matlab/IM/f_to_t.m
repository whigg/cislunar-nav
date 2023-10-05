function [dt] = f_to_t(f1, f2, a, e, mu)
%F_TO_T Computes the time difference between two points in the same orbit.
% Input:
%  - f1; initial true anomaly [rad]
%  - f2; final true anomaly [rad]
%  - a ; semimajor axis [km]
%  - e ; eccentricity
%  - mu; gravitational parameter [km^3 / s^2]

E1 = 2 * atan2(sqrt(1 - e) * tan(f1 / 2), sqrt(1 + e));
E2 = 2 * atan2(sqrt(1 - e) * tan(f2 / 2), sqrt(1 + e));
M1 = E1 - e*sin(E1);
M2 = E2 - e*sin(E2);
n = sqrt(mu / a^3);

dt = (mod(M2 - M1 + pi, 2*pi) - pi) / n;
end

