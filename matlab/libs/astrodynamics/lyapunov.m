function dP = lyapunov(P,A,Q)
%LYAPUNOV Describes the continuous-time Lyapunov equations. For use with ODE45 
%or other propagator.
%   Takes the form
%       dP = AP + PA' + Q
%   P must be passed in as a column vector rather than a matrix
%   to play nice with ODE45.
%
%   Inputs: (dims)
%    - P; ((nxn)x1) Vectorized covariance matrix
%    - A; (nxn) linearized state dynamics matrix
%    - Q; (nxn) process noise matrix (may also take form LQL')
arguments
    P   (:,1) double
    A   (:,:) double
    Q   (:,:) double
end

% column vector to square matrix
n = sqrt(size(P, 1));
P = reshape(P, n, n);   
dP = A*P + P*A' + Q;            % Lyapunov equation
% square matrix to column vector
dP = reshape(dP, n*n, 1);   
end
