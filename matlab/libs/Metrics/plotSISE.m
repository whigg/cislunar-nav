function plotSISE(t)
%PLOTSISE generates a plot of the signal-in-space errors (SISE) at the user 
%for given data and specified Intuitive Machines ephemeris error levels.
%   Inputs: (dims),[units]
%    - t   ; (nx1),[s] time steps

t = t - t(1);                   % set t0 = 0

sise1 = 3 * SISE(t, 1);         % SISE w/ ephemeris error @ 1m 3sigma
sise2 = 3 * SISE(t, 10);        % SISE w/ ephemeris error @ 10m 3sigma
sise3 = 3 * SISE(t, 100);       % SISE w/ ephemeris error @ 100m 3sigma

figure();
plot(t * 24 / 86400, sise1, 'LineWidth', 1.5);
hold on;
plot(t * 24 / 86400, sise2, '--', 'LineWidth', 1.5);
plot(t * 24 / 86400, sise3, '-.', 'LineWidth', 1.5);
axis([0 (t(end) * 24 / 86400) 0 150]);
grid on;
xlabel("Time (hrs)"); ylabel("3\sigma SISE (m)");
title("Signal-in-Space Errors (SISE) at the User");
set(gcf, 'position', [500, 250, 900, 350]);
legend(["3\sigma OD Uncertainty = 1m", "3\sigma OD Uncertainty = 10m", ...
        "3\sigma OD Uncertainty = 100m"]);
end

