function dxdt = lnss_elfodyn(t,x,moon,N,bodies)
%LNSS_ELFODYN Inertial dynamics for a spacecraft in elliptical lunar frozen
%orbits navigating using an LNSS constellation; state is position, velocity, and
%clock info (given in more detail below).
%
%   State: [x_J2000 (km);    y_J2000 (km);    z_J2000 (km);
%           vx_J2000 (km/s); vy_J2000 (km/s); vz_J2000 (km/s);
%           bias (km);       drift (km/s);    drift rate (km/s^2)]
%   Input:
%    - t; simulation time
%    - x; state of spacecraft
%    - moon; struct, {GM: gravitational parameter, x: @(t) position [km], ...}
%           for moon
%    - N; maximum degree and order of harmonics to compute
%    - earth; struct, {GM: gravitational parameter, x: @(t) position [km], ...}
%           for Earth

dxdt = zeros(9,1);
dxdt(1:6) = orbitaldynamics(t,x(1:6),moon,N,bodies);
dxdt(7:9) = clockdynamics(x(7:9));
end

