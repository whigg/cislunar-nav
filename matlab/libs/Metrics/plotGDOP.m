function plotGDOP(pos,sats,t,maxGDOP)
%PLOTGDOP Plots the # of visible satellites and GDOP over time for a 
%specific point.
%   Inputs: (dims),[units]
%    - pos ; (3x1),[km] position at which to compute GDOP
%    - sats; (nx3xm),[km] positions of m satellites over n time steps
%    - t   ; (nx1),[s] time steps
%    - maxGDOP; (1x1),[N/A] maximum GDOP to visualize on plot

if nargin < 4, maxGDOP = 6; end

t = t - t(1);                   % set t0 = 0

% Print time in days if > 2 day time span
fmt = 'days';
if t(end) < 86400 * 2, fmt = 'hours'; end
if strcmp(fmt, 'days')
    t = t / 86400;
else
    t = t / 3600;
end

[dop, nvis] = GDOP(pos, sats);                      % GDOP & # visible sats

% plot number of visible sats and GDOP over time
figure();
subplot(2,1,1);                                     % plot visible sats
plot(t, nvis, 'LineWidth', 1.5);
axis([0 t(end) 0 size(sats, 3)]);
grid on;
ylabel("# of Satellites in View");
title("Satellites in View of User");

subplot(2,1,2);                                     % plot GDOP
plot(t, dop, 'LineWidth', 1.5, 'HandleVisibility', 'off');
if maxGDOP > 6
    hold on;
    patch([t(1) t(end) t(end) t(1) t(1)], [6 6 maxGDOP maxGDOP 6], ...
        'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold off;
end
axis([0 t(end) 1 maxGDOP]);
grid on;
xlabel(sprintf("Time (%s)", fmt));
ylabel("Geometric Dilution of Precision");
title("Geometric Dilution of Precision for User");

% set(gcf, 'position', [500, 250, 750, 500]);         % fixed size
end

