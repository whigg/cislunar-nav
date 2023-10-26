function plotPDOP(pos,sats,t)
%PLOTPDOP Plots the # of visible satellites and PDOP over time for a 
%specific point.
%   Inputs: (dims),[units]
%    - pos ; (3x1),[km] position at which to compute PDOP
%    - sats; (nx3xm),[km] positions of m satellites over n time steps
%    - t   ; (nx1),[s] time steps

t = t - t(1);                   % set t0 = 0

[dop, nvis] = PDOP(pos, sats);                      % PDOP & # visible sats

% plot number of visible sats and PDOP over time
figure();
subplot(2,1,1);                                     % plot visible sats
plot(t * 24 / 86400, nvis, 'LineWidth', 1.5);
axis([0 (t(end) * 24 / 86400) 0 8]);
grid on;
ylabel("# of Satellites in View");
title("Satellites in View of User");

subplot(2,1,2);                                     % plot PDOP
plot(t * 24 / 86400, dop, 'LineWidth', 1.5);
axis([0 (t(end) * 24 / 86400) 0 10]);
grid on;
xlabel("Time (hrs)"); ylabel("Geometric Dilution of Precision");
title("Geometric Dilution of Precision for User");

set(gcf, 'position', [500, 250, 750, 500]);         % fixed size
end

