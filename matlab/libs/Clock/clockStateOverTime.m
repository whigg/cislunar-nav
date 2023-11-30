function [X, V] = clockStateOverTime(t, clk)
%CLOCKSTATEOVERTIME Returns the 3-state history of a specified clock over
%the given time period.
%   Inputs: (dims),[units]
%    - t  ; (nx1),[s] time steps to evaluate clock at, initialized at t(1)
%    - clk; (1x1),[str] clock type -- either 'CSAC', 'RAFS', or 'USO'
%   Output:
%    - X  ; (3xn),[s; 1; 1/s] 3-state clock vector over time interval
%    - V  ; (3x3xn),[?] covariance of states X over time

n = length(t);
% state vector (phase deviation, frequency deviation, frequency drift)
X = zeros(3, n);
% variance of phase deviation
V = zeros(3, 3, n);

if strcmp(clk, 'CSAC')
    [s1, s2, s3] = DiffCoeffCSAC();
    X(3,1) = 9e-10 / (86400 * 30);  % Hz/Hz/s, upper end of aging rate a
    s1 = s1*1.5;
    s2 = s2*1.5;
    s3 = X(3,1)*.5;
elseif strcmp(clk, 'RAFS')
    [s1, s2, s3] = DiffCoeffRAFS();
elseif strcmp(clk, 'USO')
    [s1, s2, s3] = DiffCoeffUSO();
else
    throw(MException('clockStateOverTime:ClockType', ...
        'Invalid clock type of %s. See documentation.', clk));
end

Phi = @(tau) [1 tau tau^2/2;
              0 1   tau    ;
              0 0   1       ];  % state trans. matrix, tau = t_k - t_k-1
% covariance matrix of error associated w/ Wiener processes
Q = @(s1,s2,s3,tau) ...
     [s1^2*tau + s2^2/3*tau^3 + s3^2/20*tau^5 s2^2/2*tau^2 + s3^2/8*tau^4 s3^2/6*tau^3;
      s2^2/2*tau^2 + s3^2/8*tau^4             s2^2*tau + s3^2/3*tau^3     s3^2/2*tau^2;
      s3^2/6*tau^3                            s3^2/2*tau^2                s3^2*tau     ];

for i=2:n
    dt = t(i) - t(i-1);
    STM = Phi(dt);

    % innovation vector, J ~ N(0,Q)
    J = mvnrnd([0 0 0], Q(s1, s2, s3, dt), 1)';
    V(:,:,i) = Q(s1, s2, s3, t(i)); % covariance of X at t(i)
    X(:,i) = STM * X(:,i-1) + J;
end
end

