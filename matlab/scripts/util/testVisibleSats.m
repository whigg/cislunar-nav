% testVisibleSats.m
% Author: Mark Hartigan
% Date  : October 24, 2023
% Description:
%    Verify that the view cone for different user locations is working
%    correctly.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% test
r = 1736;               % km, polar radius of moon
sz = r + 1000;          % km, size of cone to make
pos = [r; 0; 0];     % km, position of user
sats = zeros(0,0,0);    % km, satellites available (none!)
elev = 5 * pi/180;      % rad, elevation angle

% get rotation matrix and zcone function
[~, zcone, R] = visibleSats(pos, sats, 0, elev);
[X,Y] = meshgrid(-sz:100:sz, -sz:100:sz);   % create grid of points

% fill out z-coordinates
Z = zeros(size(X));                     
for i=1:size(X,1)
    for j=1:size(X,2)
        Z(i,j) = zcone(X(i,j), Y(i,j));
    end
end

figure();

% plot moon
[x, y, z] = ellipsoid(0, 0, 0, r, r, r, 20);
globe = surf(x, y, z);
hold on;

% plot cone
surf(X,Y,Z);
hold off; axis equal; view(15,24);

xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');