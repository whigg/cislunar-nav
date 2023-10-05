function [visible, R] = visibleSats(pos, sats, time, elev)
%VISIBLESATS Compute number of satellites visible to the user (on lunar
%surface).
%   Checks if satellite is below visibility cone -- rotates all points not
%   on the south pole to end up on the south pole.

if nargin < 4, elev = 5 * pi/180; end   % default elevation is 5 deg

r = 1737.4; % km
zcone = @(x,y) -r - tan(elev) * norm([x y]);

rv = [0 0 -r]';
% get Euler axis and angle for rotation
if norm(cross(pos, rv)) ~= 0
    ax = cross(pos, rv) / norm(cross(pos, rv));
    ang = acos(dot(pos, rv) / (norm(pos) * norm(rv)));  % dot product
else
    ax = [1 0 0]';
    if norm(pos - rv) > r, ang = pi; else, ang = 0; end
end

U = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
R = eye(3) + sin(ang) * U + (1 - cos(ang)) * U^2;

visible = [];
for i=1:size(sats,3)
    loc = R * sats(time,:,i)';
    if loc(3) < zcone(loc(1), loc(2))
        visible = [visible i];
    end
end
end

