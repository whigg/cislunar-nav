function visualizeCoverage(pos,sats)
%VISUALIZECOVERAGE Generates an animation showing spacecraft trajectories
%over time while also displaying GDOP at the lunar south pole.
%   Inputs: (dims),[units]
%    - pos ; (3x1),[km] position at which to compute GDOP
%    - sats; (nx3xm),[km] positions of m satellites over n time steps

dop = GDOP(pos, sats);              % GDOP & # visible sats

r = 1737.4;                         % km, moon radius
m = size(sats, 3);                  % number of satellites
n = size(sats, 1);                  % number of time steps
% colors of up to 5 satellites
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];
visible = visibleSats(pos, sats);   % matrix of visible satellites
tail = 50;                          % max points in traj to visualize
fname = "coverage_vis.gif";

% set up plotting environment
figure();
hold on; grid minor; axis equal;
xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");
% title that also displays GDOP
title(sprintf('Khon 2-6 28-day Trajectories\nGDOP: %0.2f', dop(1)), ...
    'Interpreter', 'Latex');
view(0,0);                    % Setting viewing angle
% plot with no color to set axis limits
plot3(sats(:,1,m), sats(:,2,m), sats(:,3,m), 'Color', 'none');
plot3(sats(:,1,m-1), sats(:,2,m-1), sats(:,3,m-1), 'Color', 'none');

% moon
[I, ~] = imread("lroc_color_poles_1k.jpg");
[xx, yy, zz] = ellipsoid(0, 0, 0, r, r, r);
globe = surf(xx, yy, -zz);
set(globe, 'FaceColor', 'texturemap', 'CData', I, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');

% plot first iteration
p = {}; s = {}; v = {};
for i = 1:m
    p{i} = plot3(sats(1,1,i), sats(1,2,i), sats(1,3,i), ':', ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    s{i} = scatter3(sats(1,1,i), sats(1,2,i), sats(1,3,i), 'filled', ...
        'CData', colors(i,:));
    v{i} = plot3([pos(1), sats(1,1,i)], [pos(2), sats(1,2,i)], ...
        [pos(3), sats(1,3,i)], 'Color', 'none');
    if visible(1,i), v{i}.Color = 'r'; end
end

for j = 1:n                         % iterate over every time step
    for i = 1:m                     % update lines
        % update satellite trajectory (only keep previous 50 points)
        p{i}.XData = sats(max(1,j-tail):j,1,i);
        p{i}.YData = sats(max(1,j-tail):j,2,i);
        p{i}.ZData = sats(max(1,j-tail):j,3,i);
        % update satellite point
        s{i}.XData = sats(j,1,i);
        s{i}.YData = sats(j,2,i);
        s{i}.ZData = sats(j,3,i);
        if visible(j,i)             % update LOS line
            v{i}.XData = [pos(1), sats(j,1,i)];
            v{i}.YData = [pos(1), sats(j,2,i)];
            v{i}.ZData = [pos(1), sats(j,3,i)];
            v{i}.Color = 'r';
        else
            v{i}.Color = 'none';
        end
    end

    % update title
    title(sprintf('Khon 2-6 28-day Trajectories\nGDOP: %0.2f', dop(j)), ...
        'Interpreter', 'Latex');
    pause(0.01)                     % delay
    % save figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if j == 1
        imwrite(imind,cm,fname,'gif','Loopcount',inf,'Delaytime',0.05);
    else
        imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',0.05);
    end
end

hold off; 
end

