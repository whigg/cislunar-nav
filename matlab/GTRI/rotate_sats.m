%% reset
clc, clear, close all;

%% init
r = 1737.4;
x0 = [sin(pi/12)*r, 0., -cos(pi/12)*r];
xf = [0, 0, -r];
ang = acos(dot(x0, xf) / (norm(x0) * norm(xf)));    % dot product formula
ax = cross(x0, xf) / norm(cross(x0, xf));

%% rotate
% axis-angle rotation
U = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
R = cos(ang) * eye(3) + sin(ang) * U + (1 - cos(ang)) * (ax'*ax);
R * x0'