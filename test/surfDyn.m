function [dx] = surfDyn(x, u, g, r, W)
%SURFDYN dynamical equations for lunar surface user
%   Input
%    - x; state vector at current time [pos; vel], (6,)
%    - A; state matrix (6x6)
%    - u; input acceleration (walk on planet)
%    - g; planetary acceleration
%    - r; hard planetary radius

    dx = zeros(6,1);
    dx(1:3) = x(4:6);

    p = x(1:3);
    dir = p / norm(p);
    vrel = x(4:6) - cross(W, p);    % velocity relative to moon-fixed frame

    u = u + 2*cross(W, vrel) + cross(W, cross(W, p)) - dir*g;
    if norm(p) <= r, u = u + dir * max(dot(u, -dir), 0); end

    dx(4:6) = u;
end