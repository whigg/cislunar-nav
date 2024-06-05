function pts = getpoints(n)
%GETPOINTS Returns points on the lunar surface (dispersed for optimization
%problems, i.e. along the edges of LunaNet SV2.
%   Input:
%    - n; number of points -- 1 at LSP, the rest around the edges
arguments
    n (1,1) {mustBePositive,mustBeInteger}
end

R = 1738;                       % km, radius of moon
pts = zeros(3,n);
pts(:,1) = [0;0;-R];            % south pole
lat = -75 * pi/180;             % highest latitude of SV2
pts(3,2:end) = R * sin(lat);    % z component of remaining points
xy = R * cos(lat);              % x/y component of remaining points
lon = linspace(0,2*pi,n);       % longitudes
pts(1,2:end) = xy * cos(lon(1:end-1));
pts(2,2:end) = xy * sin(lon(1:end-1));
end

