%% reset
clc, clear, close all;

%% init
w = 2*pi / (27.3217 * 24 * 60 * 60);    % rad/s, rotation rate of moon
r = 1737.4;                             % km, radius of moon
% 1-3 km, 4-6 km/s, initial position and velocity
x0 = [sind(15); 0; -cosd(15); 0; w*sind(15)*r; 0] * r;

A = [0 -w 0 0  0 0;
     w  0 0 0  0 0;
     0  0 0 0  0 0;
     0  0 0 0 -w 0;
     0  0 0 w  0 0;
     0  0 0 0  0 0];

phi = @(t) [cos(w*t) -sin(w*t) 0        0         0 0;
            sin(w*t)  cos(w*t) 0        0         0 0;
                   0         0 1        0         0 0;
                   0         0 0 cos(w*t) -sin(w*t) 0;
                   0         0 0 sin(w*t)  cos(w*t) 0;
                   0         0 0        0         0 1];
t = [0 1 3 7] * 24 * 60 * 60;         % s, times to eval

%% eval
x = zeros(length(x0), length(t));
x(:,1) = x0;

for i=2:length(t)
    x(:,i) = phi(t(i)) * x0;
end

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[t_nom, x_nom] = ode45(@(t,x) A*x, [t(1) t(end)], x0, opts);

%% plot
%x = x_nom'; % test ode45 instead

figure(1);
s = 100;   % scale factor for velocity vectors

for i=1:length(t)
    quiver3(0, 0, 0, x(1,i), x(2,i), x(3,i), 'b');
    if i==1, hold on; end
    quiver3(x(1,i), x(2,i), x(3,i), x(4,i)*s, x(5,i)*s, x(6,i)*s, 'r');
end
hold off; grid on; axis equal;
xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");
