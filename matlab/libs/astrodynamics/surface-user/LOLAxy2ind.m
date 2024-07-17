function [I,J] = LOLAxy2ind(X,Y,N,MAP_SCALE)
%LOLAXY2IND Converts X/Y coordinates in stereographic lunar south pole
%frame used by LOLA to pixel indices for data maps.
%   Input:
%    - X; x coordinate in km (+Y in MOON_ME)
%    - Y; y coordinate in km (+X in MOON_ME)
%    - N; number of rows/cols (images are square)
%    - MAP_SCALE; kilometers per pixel
%   Output:
%    - I; image row
%    - J; image column
arguments
    X (1,1) double
    Y (1,1) double
    N (1,1) {mustBeInteger,mustBePositive}
    MAP_SCALE (1,1) double {mustBePositive}
end

J = round(X/MAP_SCALE + N/2 + 0.5);
I = round(Y/MAP_SCALE + N/2 + 0.5);
end

