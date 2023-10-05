function [dx] = surfdyn(x, u, g, r, W)
%SURFDYN dynamical equations for lunar surface user
%   Inputs: (dims),[units]
%    - x ; (6x1),[km,km/s] state vector at current time
%    - A ; (6x6),[misc] state matrix (6x6)
%    - u ; (3x1),[km/s^2] input acceleration (walk on planet)
%    - g ; (1x1),[km/s^2] planetary acceleration
%    - r ; (1x1),[km] hard planetary radius
%   Output:
%    - dx; (6x1),[km/s,km/s^2] derivative of state vector

dx = zeros(6,1);
dx(1:3) = x(4:6);

p = x(1:3);
dir = p / norm(p);
vrel = x(4:6) - cross(W, p);    % velocity relative to moon-fixed frame

u = u + 2*cross(W, vrel) + cross(W, cross(W, p)) - dir*g;
if norm(0) <= r, u = u + dir * max(dot(u, -dir), 0); end

dx(4:6) = u;
end