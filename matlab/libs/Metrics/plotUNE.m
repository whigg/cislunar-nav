function plotUNE(pos,sats,t)
%PLOTUNE generates a plot of the user navigation error (UNE) at the user 
%for given data and specified Intuitive Machines ephemeris error levels.
%   Inputs: (dims),[units]
%    - pos ; (3x1),[km] position at which to compute GDOP
%    - sats; (nx3xm),[km] positions of m satellites over n time steps
%    - t   ; (nx1),[s] time steps

t = t - t(1);                   % set t0 = 0

sise1 = SISE(t, 1);             % SISE w/ ephemeris error @ 1m 3sigma
sise2 = SISE(t, 10);            % SISE w/ ephemeris error @ 10m 3sigma
sise3 = SISE(t, 100);           % SISE w/ ephemeris error @ 100m 3sigma

% dilution of precision @ south pole
dop = PDOP(pos, sats);
une1 = dop .* sise1;
une2 = dop .* sise2;
une3 = dop .* sise3;

figure();
plot(t * 24 / 86400, une1, 'LineWidth', 1.5);
hold on;
plot(t * 24 / 86400, une2, '--', 'LineWidth', 1.5);
plot(t * 24 / 86400, une3, '-.', 'LineWidth', 1.5);
axis([0 (t(end) * 24 / 86400) 0 300]);
grid on;
xlabel("Time (hrs)"); ylabel("3\sigma RSS Position Error (m)");
title("User Navigation Error (UNE) at the User");
set(gcf, 'position', [500, 250, 900, 350]);
legend(["3\sigma OD Uncertainty = 1m", "3\sigma OD Uncertainty = 10m", ...
        "3\sigma OD Uncertainty = 100m"]);
end