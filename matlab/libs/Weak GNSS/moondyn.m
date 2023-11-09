function [dxdt] = moondyn(t, x, moon, earth, sun)
%MOONDYN Force model for a moon-orbiting satellite
%   Input:
%    - t; current time (seconds past J2000)
%    - x; spacecraft state [position (km), velocity (km/s)]
%    - moon ; planet struct {'GM': [km^3/s^2]; 'x': @(t) [km; km/s]}
%    - earth; planet struct {'GM': [km^3/s^2]; 'x': @(t) [km; km/s]}
%    - sun  ; planet struct {'GM': [km^3/s^2]; 'x': @(t) [km; km/s]}

r = x(1:3);         % position of satellite from moon
r_e = earth.x(t);   % position of earth from moon
r_s = sun.x(t);     % position of sun from moon

% point mass of moon, earth, and sun
f = -moon.GM / norm(r)^3 * r;
f = f + earth.GM * ((r_e - r)/norm(r_e - r)^3 - r_e/norm(r_e)^3);
f = f + sun.GM   * ((r_s - r)/norm(r_s - r)^3 - r_s/norm(r_s)^3);

dxdt = [x(4:6); f; 0];
end

