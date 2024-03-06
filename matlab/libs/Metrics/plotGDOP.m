function plotGDOP(pos,sats,t)
%PLOTGDOP Plots the # of visible satellites and GDOP over time for a 
%specific point.
%   Inputs: (dims),[units]
%    - pos ; (3x1),[km] position at which to compute GDOP
%    - sats; (nx3xm),[km] positions of m satellites over n time steps
%    - t   ; (nx1),[s] time steps

t = t - t(1);                   % set t0 = 0

[dop, nvis] = GDOP(pos, sats);                      % GDOP & # visible sats

% plot number of visible sats and GDOP over time
figure();
subplot(2,1,1);                                     % plot visible sats
plot(t / 86400, nvis, 'LineWidth', 1.5);
axis([0 (t(end) / 86400) 0 5]);
grid on;
ylabel("# of Satellites in View");
title("Satellites in View of User");

subplot(2,1,2);                                     % plot GDOP
plot(t / 86400, dop, 'LineWidth', 1.5);
axis([0 (t(end) / 86400) 1 6]);
grid on;
xlabel("Time (days)"); ylabel("Geometric Dilution of Precision");
title("Geometric Dilution of Precision for User");

set(gcf, 'position', [500, 250, 750, 500]);         % fixed size
end

