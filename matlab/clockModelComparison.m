% clockModelComparison.m
% Author: Mark Hartigan
% Date  : November 8, 2023
% Description:
%    Compare the timing performance of a CSAC, USO, and RAFS.

%% reset
clc, clear, close all;
addpath(genpath(pwd));

%% init
a_hi = 9e-10 / (86400 * 30);        % Hz/Hz/s, upper end of aging rate a
c = 299792458;                      % m/s, speed of light
[csac1, csac2, csac3] = DiffCoeffCSAC();    % chip-scale atomic clock
[uso1 , uso2 , uso3 ] = DiffCoeffUSO();     % ultra-stable oscillator
[rafs1, rafs2, rafs3] = DiffCoeffRAFS();    % rubidium atomic frequency std
% assume constants assoc. w/ Wiener processes mu_1, mu_2, mu_3 = 0

%% Compute uncertainty of clock over time
Phi = @(tau) [1 tau tau^2/2;
              0 1   tau    ;
              0 0   1       ];  % state trans. matrix, tau = t_k - t_k-1
% covariance matrix of error associated w/ Wiener processes
Q = @(s1,s2,s3,tau) ...
     [s1^2*tau + s2^2/3*tau^3 + s3^2/20*tau^5 s2^2/2*tau^2 + s3^2/8*tau^4 s3^2/6*tau^3;
      s2^2/2*tau^2 + s3^2/8^tau^4             s2^2*tau + s3^2/3*tau^3     s3^2/2*tau^2;
      s3^2/6*tau^3                            s3^2/2*tau^2                s3^2*tau     ];

dt = 1;                         % time step
t = 0:dt:3600*5;                % measure once per dt for a time period
n = length(t);                  % number of time steps
STM = Phi(dt);                  % time step state transition matrix

% state vector (phase deviation, frequency deviation, frequency drift)
Xcsac = zeros(3,n); Xcsac(3,1) = a_hi;  % add aging rate
Xuso  = zeros(3,n);
Xrafs = zeros(3,n);
varCsac = zeros(1,n);           % variance of phase deviation
varUso  = zeros(1,n);           % variance of phase deviation
varRafs = zeros(1,n);           % variance of phase deviation

for i=2:n                       % compute state vector for time span
    % Chip-scale atomic clock
    % innovation vector, J ~ N(0,Q)
    J = mvnrnd([0 0 0], Q(csac1, csac2, csac3, dt))';
    SIG = Q(csac1, csac2, csac3, t(i)); % covariance of X at t(i)
    varCsac(i) = SIG(1,1);              % variance of X(1,i)
    Xcsac(:,i) = STM * Xcsac(:,i-1) + J;

    % Ultra-stable oscillator
    J = mvnrnd([0 0 0], Q(uso1, uso2, uso3, dt))';
    SIG = Q(uso1, uso2, uso3, t(i));
    varUso(i) = SIG(1,1);
    Xuso(:,i) = STM * Xuso(:,i-1) + J;

    % Rubidium atomic frequency standard
    J = mvnrnd([0 0 0], Q(rafs1, rafs2, rafs3, dt))';
    SIG = Q(rafs1, rafs2, rafs3, t(i));
    varRafs(i) = SIG(1,1);
    Xrafs(:,i) = STM * Xrafs(:,i-1) + J;
end

%% plot CSAC clock bias over one day
figure();
plot(t ./ 3600, abs(Xcsac(1,:)), 'LineWidth', 1.5, 'Color', '#0072BD');
hold on;
plot(t ./ 3600, 3*sqrt(varCsac), '--', 'LineWidth', 1.5, 'Color', '#0072BD');
plot(t ./ 3600, 1/2 * a_hi * t.^2, '-.', 'LineWidth', 1.5);
hold off; grid on; xlim([0 t(end)/3600]);
xlabel('Time (hrs)');
ylabel('Phase deviation / clock bias error (s)');
legend(["Monte-Carlo run", "3\sigma bound", "Old model (only aging)"], ...
    'location', 'northwest');
title('Microsemi Space CSAC Timing Performance');

%% comparison of CSAC, USO, and RAFS
figure();
% individual Monte-Carlo runs
semilogy(t ./ 3600, abs(Xuso(1,:)) * c , 'LineWidth', 1.5, 'Color', '#D95319');
hold on;
semilogy(t ./ 3600, abs(Xcsac(1,:)) * c, 'LineWidth', 1.5, 'Color', '#0072BD');
semilogy(t ./ 3600, abs(Xrafs(1,:)) * c, 'LineWidth', 1.5, 'Color', '#EDB120');
% 3-sigma bounds
semilogy(t ./ 3600, 3*sqrt(varUso) * c , '--', 'LineWidth', 2, 'Color', '#D95319');
semilogy(t ./ 3600, 3*sqrt(varCsac) * c, '--', 'LineWidth', 2, 'Color', '#0072BD');
semilogy(t ./ 3600, 3*sqrt(varRafs) * c, '--', 'LineWidth', 2, 'Color', '#EDB120');
hold off; grid on; xlim([0 t(end)/3600]);
xlabel('Time (hrs)');
ylabel('Phase deviation / clock bias error (m)');
legend(["USO MC run", "CSAC MC run", "RAFS MC run", "USO 3\sigma bound", ...
        "CSAC 3\sigma bound", "RAFS 3\sigma bound"], 'location', 'best', 'NumColumns', 2);
title('Comparison of Different Clock Timing Performances');
