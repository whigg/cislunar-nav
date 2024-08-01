function xf = surfaceuserdyn(ts,x0,fa)
%SURFACEUSERDYN Describes the discrete-time dynamics of a lunar surface
%user -- basically just integrating acceleration.
%   Input:
%    - ts; simulation start and end time
%    - x0; state at start of simulation
%    - dt; size of time step
%    - fa; acceleration function
arguments
    ts  (1,2) double {mustBePositive}
    x0  (9,1) double
    fa  (1,1) function_handle
end

t0 = ts(1);
dt = ts(2) - ts(1);
F = [eye(3) eye(3) * dt; zeros(3,3) eye(3)];
G = [eye(3) * 1/2*dt^2; eye(3) * dt; zeros(3,3)];
H = [1 dt 1/2*dt^2; 0 1 dt; 0 0 1];
F = [F zeros(6,3); zeros(3,6) H];
xf = F * x0 + G * fa(t0);
end

