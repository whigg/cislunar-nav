function dfdx = surfaceuserpartials(ts,x0)
%SURFACEUSERPARTIALS Partial derivates of surface user discrete-time
%dynamics.
%   Input:
%    - ts; simulation start and end time
%    - x0; state at start of simulation
arguments
    ts  (1,2) double {mustBePositive}
    x0  (9,1) double
end

dt = ts(2) - ts(1);
dfdx = [eye(3) eye(3) * dt; zeros(3,3) eye(3)];
H = [1 dt 1/2*dt^2; 0 1 dt; 0 0 1];
dfdx = [dfdx zeros(6,3); zeros(3,6) H];
end

