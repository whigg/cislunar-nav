function [] = constellation_metrics(data)
%CONSTELLATION_METRICS generates evaluation metrics about a given lunar
%constellation. m data points are provided for k satellites.
% Input:
%  - data; time history of satellite positions (m x (1 + 3*n)) [s, km...]

t = data(:,1);                  % time history
m = length(t);                  % number of data points
n = (size(data,2) - 1) / 3;     % number of satellites
r = 1737.4; % km
sats = zeros(m, 3, n);
for num = 1:n
    sats(:,:,num) = data(:,2 + 3*(num-1):1 + 3*num);
end

% sats = sats(:,:,1:4);

f1 = figure();
set(f1, 'Color', 'w');
hold on; axis equal; view(15,24);
[x, y, z] = ellipsoid(0, 0, 0, r, r, r, 40);
  
% visibility visualization
vis = zeros(size(x));
gslon = zeros(size(x));
gslat = zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        gslon(i,j) = atan2(y(i,j), x(i,j)) * 180/pi;
        gslat(i,j) = atan2(z(i,j), norm([x(i,j) y(i,j)])) * 180/pi;
        for k=1:m
            vis(i,j) = vis(i,j) + length( ...
                       visibleSats([x(i,j) y(i,j) z(i,j)]', sats, k) );
        end
    end
end

vis = vis ./ m;
labels = [];

for i=1:size(sats,3)
    labels = [labels, {sprintf('Khon%d',i)}];
    plot3(sats(:,1,i), sats(:,2,i), sats(:,3,i), 'LineWidth', 1);
end

for i=1:size(sats,3)
    scatter3(sats(end,1,i), sats(end,2,i), sats(end,3,i), 'ro');
end

globe = surf(x, y, z);
set(globe, 'FaceColor', 'texturemap', 'CData', vis, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
c = colorbar;
c.Label.String = "Avg. # of Satellites in Contact";
title("24-hour LNSS Orbit Propagation");
legend(labels, 'location', 'best');

% 2D coverage plot
figure();
gsLatLonAvgview = zeros(length(gslat)^2,3);
track = 1;
for i = 1:length(gslat)
    for j = 1:length(gslat)
        gsLatLonAvgview(track,1) = gslon(i,j);
        gsLatLonAvgview(track,2) = gslat(i,j);
        gsLatLonAvgview(track,3) = vis(i,j);
        track = track + 1;
    end
end

moon = imread('lroc_color_poles_1k.jpg');
moonResize = imresize(moon,1/2.85);
imshow(moonResize, 'XData', [min(gsLatLonAvgview(:,1)) max(gsLatLonAvgview(:,1))], 'YData', [min(gsLatLonAvgview(:,2)) max(gsLatLonAvgview(:,2))]);
hold on
revCM = flip(gsLatLonAvgview(:,3));
revLat = flip(gsLatLonAvgview(:,2));
revLon = flip(gsLatLonAvgview(:,1));
scatter(revLon,gsLatLonAvgview(:,2),15,revCM,'filled')
% xlim([-180,180])
% ylim([-90,90])
ylabel('Latitude')
xlabel('Longitude')
title('Average Number of Satellites in Contact (Khon1-6)')
axis on
ax = gca;
ax.YTickLabel = flipud(ax.YTickLabel);
set(ax,'CLim',[0 4]);
colorbar
% pbaspect([2 1 1])
set(gcf,'position',[750,500,750,400])

% dilution of precision @ south pole
nvis = zeros(size(t));
dop = zeros(size(t));
for k=1:m
    H = diag(computeDOP([0; 0; -r], sats, k));
    nvis(k) = length(visibleSats([0; 0; -r], sats, k));
    dop(k) = sqrt(abs(sum(H(1:3)))); 
end

figure();
subplot(2,1,1);
plot((t-t(1)) * 24 / 86400, nvis, 'LineWidth', 1.5);
hold on;
axis([0 ((t(end)-t(1)) * 24 / 86400) 0 8]);
grid on;
ylabel("# of Satellites in View");
title("Satellites in View of the Lunar South Pole");
subplot(2,1,2);
plot((t-t(1)) * 24 / 86400, dop, 'LineWidth', 1.5);
hold on;
axis([0 ((t(end)-t(1)) * 24 / 86400) 0 10]);
grid on;
xlabel("Time (hrs)"); ylabel("Dilution of Precision");
title("Dilution of Precision at the Lunar South Pole");
set(gcf, 'position', [500, 250, 750, 500]);
end
