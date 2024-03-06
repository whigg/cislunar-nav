function plotCoverage(pts,covered,moon)
%PLOTCOVERAGE Plots coverage for a lunar service volume output by
%coverageLNSP.
%   Inputs: (dims),[units]
%    - pts    ; (3xl),[km] points coverage was evaluated at
%    - covered; (nxl) availability (or not) of l points at n time steps
%    - moon   ; (1x1),[N/A] optional (default false) boolean -- plot moon?

if nargin < 3, moon = false; end

n = size(covered, 1);                   % time steps
CVGs = sum(covered, 1) ./ n .* 100;     % coverage for each pt as a %

% put data in table to utilize ColorVariable in scatter3()
tab = array2table([pts' CVGs'], 'VariableNames', {'X', 'Y', 'Z', 'CVG'});

% start plotting
figure();
% eval points
scatter3(tab, 'X', 'Y', 'Z', 'filled', 'ColorVariable', 'CVG', ...
    'HandleVisibility', 'off');
hold on;

% plot minimum point in red
[~, I] = min(CVGs);
scatter3(pts(1,I), pts(2,I), pts(3,I), 'r', 'filled', ...
    'DisplayName', 'Least Coverage');

% optional moon plotting, makes actual data hard to see
if moon
    [I, ~] = imread("lroc_color_poles_1k.jpg");
    r = 1736;                           % km, moon polar radius
    [xx, yy, zz] = ellipsoid(0, 0, 0, r, r, r);
    globe = surf(xx, yy, -zz, 'HandleVisibility', 'off');
    set(globe, 'FaceColor', 'texturemap', 'CData', I, 'FaceAlpha', 1, ...
        'EdgeColor', 'none');
end

hold off; grid on; axis equal;
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
legend('location', 'best');
c = colorbar;
c.Label.String = 'Coverage (% of Earth day)';
clim([0 100]);                          % limit colorbar to 0-100%
end

