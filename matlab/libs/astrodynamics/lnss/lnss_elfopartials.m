function A = lnss_elfopartials(x,moon)
%LNSS_ELFOPARTIALS Returns the Jacobian of lnss_elfodyn w.r.t. x.
%
%   State: [x_J2000 (km);    y_J2000 (km);    z_J2000 (km);
%           vx_J2000 (km/s); vy_J2000 (km/s); vz_J2000 (km/s);
%           bias (km);       drift (km/s);    drift rate (km/s^2)]
%   Input:
%   Input:
%    - x; state [pos (km); vel (km/s)] of spacecraft
%    - moon; struct, {GM: gravitational parameter, x: @(t) position [km]}
%           for moon

A = zeros(9,9);
A(1:6,1:6) = orbitalpartials(x(1:6),moon);
A(7:9,7:9) = [0 1 0; 0 0 1; 0 0 0];
end

