function plotViewOutages(t, nview, ybnds)
%PLOTVIEWOUTAGES Overlay shaded regions on current plot that display when
%there are no satellites in view.
%   Inputs: (dims),[units]
%    - t    ; (1xn),[?] time at each of n steps
%    - nview; (1xn),[N/A] # satellites in view for each of n time steps
%    - ybnds; (1x2),[?] [lower,upper] bound of current plot

n = length(nview);
start = -1;
X = [];
Y = [];

for i=1:n
    % start tracking outage if 0 sats seen and outage not already started
    if nview(i) == 0 && start == -1, start = i; end
    % if outage ends or reach end of time period
    if (nview(i) ~= 0 || i == n) && start ~= -1
        % add polygon to be plotted in patch() and end outage
        X = [X [t(start); t(i); t(i); t(start)]];
        Y = [Y [ybnds(2); ybnds(2); ybnds(1); ybnds(1)]];
        start = -1;
    end
end

patch('XData', X, 'YData', Y, 'FaceColor', 'red', 'FaceAlpha', 0.3, ...
      'EdgeColor', 'none');
end

