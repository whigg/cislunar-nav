%% reset
clc, clear, close all;
addpath(genpath(pwd));

%% init
taus = [1 10 100 1000]';            % intervals for allan deviation
stds = [3e-10 1e-10 3e-11 1e-11]';  % deviations for different intervals
a_hi = 9e-10 / (86400 * 30);        % Hz/Hz/s, upper end of aging rate a
a_lo = 6e-10 / (86400 * 30);        % Hz/Hz/s, lower end of aging rate a
c = 299792458;                      % m/s, speed of light

%% solve for diffusion coefficients of noise types
% Assume frequency drift a is constant, c_3 = a. Zucca and Tavella, 2005, 
% give the Allan variance formula as
%    s_y^2 (t) = s_1^2/t + s_2^2/3 * t + c_3^2 / 2 * t^2
% Since c_3 is known, the formula can become
%    s_y^2 (t) - c_3^2 / 2 * t^2 = s_1^2/t + s_2^2/3 * t
% We can then do a polynomial fit to this with the provided data points in
% the MicroSemi CSAC datasheet. The solution to Ax = b takes
% the form (A' * A)^(-1) * A' * b.

A = [1./taus taus];
b = stds.^2 - a_hi^2 / 2 * taus.^2;
coeff = [(A' * A) \ A' * b; a_hi^2 / 2];
s1 = sqrt(coeff(1));                % diffusion coefficient of white noise
s2 = sqrt(3 * coeff(2));            % ^ of random walk frequency noise
s3 = 0;                             % set to zero (aging rate is constant)

% assume constants assoc. w/ Wiener processes mu_1, mu_2, mu_3 = 0

%% Compute uncertainty of clock over time
Phi = @(tau) [1 tau tau^2/2;
              0 1   tau    ;
              0 0   1       ];  % state trans. matrix, tau = t_k - t_k-1
% covariance matrix of error associated w/ Wiener processes
Q = @(tau) [s1^2*tau + s2^2/3*tau^3 + s3^2/20*tau^5 s2^2/2*tau^2 + s3^2/8*tau^4 s3^2/6*tau^3;
            s2^2/2*tau^2 + s3^2/8^tau^4             s2^2*tau + s3^2/3*tau^3     s3^2/2*tau^2;
            s3^2/6*tau^3                            s3^2/2*tau^2                s3^2*tau     ];

dt = 1;                         % time step
t = 0:dt:3600*5;                % measure once per dt for a time period
n = length(t);                  % number of time steps
% state vector (phase deviation, frequency deviation, frequency drift)
X = zeros(3,n);
var = zeros(1,n);               % variance of phase deviation

for i=2:n                       % compute state vector for time span
    J = mvnrnd([0 0 0], Q(dt))';    % innovation vector, J ~ N(0,Q)
    SIG = Q(t(i));                  % covariance of X at t(i)
    var(i) = SIG(1,1);              % variance of X(1,i)
    X(:,i) = Phi(dt) * X(:,i-1) + J;
end

%% plot clock bias over one day
figure();
plot(t ./ 3600, abs(X(1,:)) * c, 'LineWidth', 1.5);
hold on;
plot(t ./ 3600, 3*sqrt(var) * c, 'r--', 'LineWidth', 1.5);
plot(t ./ 3600, 1/2 * a_hi * t.^2 * c, '-.', 'LineWidth', 1.5);
hold off; grid on; xlim([0 t(end)/3600]);
xlabel('Time (hrs)');
ylabel('Phase deviation / clock bias error (m)');
legend(["Sample run", "3\sigma bound", "Old model (only aging)"], ...
    'location', 'northwest');
title('Microsemi Space CSAC Timing Performance');

% %% plot model to confirm it's working
% t = linspace(taus(1), taus(end), 1000)';
% mdl = sqrt([1./t t t.^2] * coeff);
% 
% figure();
% scatter(taus, stds, 'rx');
% hold on;
% plot(t, mdl, 'LineWidth', 1.5);
% hold off; grid on;
% xlabel('\tau (s)'); ylabel('variance');