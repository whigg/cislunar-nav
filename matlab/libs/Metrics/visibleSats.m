function [visible,zcone,R] = visibleSats(pos,sats,time,elev)
%VISIBLESATS Compute number of satellites visible to the user (on lunar
%surface).
%   Checks if satellite is below visibility cone -- rotates all points not
%   on the south pole to end up on the south pole. Then the tangent
%   geometry is given by:
%                      __ /|
%            pos  __ /     |
%            __ /          | r_moon
%       __ / angle        _|
%     /____)_____________|_|
%   The moon is modeled as a perfect sphere with radius 1736 km (polar
%   radius of moon). Then the view cone is computed and given a degree
%   offset specified by elev.
%   
%   Inputs: (dims),[units]
%    - pos ;(3x1),[km] position of user w.r.t. moon centered inertial
%           frame; norm must be greater than 1736 KM (not inside moon)
%    - sats;(nx3xm),[km] positions of m satellites over n time steps
%    - time;(1x1),[N/A] index of time step to evaluate sats at
%    - elev;(1x1),[rad] elevation angle offset of view cone
%   Output:
%    - visible;(1x?),[N/A] indices of visible satellites at time
%    - zcone; function handle to compute in-bounds / out cone
%    - R   ;(3x3),[N/A] rotation matrix to ?-?-Down frame

if nargin < 4, elev = 5 * pi/180; end   % default elevation is 5 deg

r = 1736.0; % km, polar radius of moon (concerned w/ south pole mostly)
ang = asin(r/norm(pos));
zcone = @(x,y) -norm(pos) + tan((pi/2 - ang) - elev) * norm([x y]);

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

