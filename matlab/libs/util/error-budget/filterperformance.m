function filterperformance(filter,ts,err,P,err_rtn,P_rtn,visible)
%FILTERPERFORMANCE Produces summary statistics and plots for the
%performance of a joint state and clock navigation filter (both in provided and
%radial, tangential, normal (RTN) frame).
%   Input:
%    - filter; name of filter for display purposes
%    - ts; time stamps
%    - err; state error
%    - P; state covariance
%    - err_rtn; state error in RTN frame
%    - P_rtn; state covariance in RTN frame
%    - visible; visible satellite array

arguments
    filter  (1,1)   string
    ts      (1,:)   double
    err     (9,:)   double
    P       (9,9,:) double
    err_rtn (9,:)   double
    P_rtn   (9,9,:) double
    visible (:,:)   {mustBeInteger}
end

clc;
ns = size(err, 2);
t0 = ts(1);

rerr = sqrt(sum(err(1:3,:).^2, 1));
rstd = reshape(sqrt(P(1,1,:)+P(2,2,:)+P(3,3,:))*3, 1, ns);
verr = sqrt(sum(err(4:6,:).^2, 1));
vstd = reshape(sqrt(P(4,4,:)+P(5,5,:)+P(6,6,:))*3, 1, ns);
berr = abs(err(7,:));
bstd = reshape(sqrt(P(7,7,:))*3, 1, ns);

fprintf("FILTER: " + filter + "\n");
fprintf("========================\n");
fprintf("  Mean error: %.3f (%.3f 3-sigma) m\n", mean(rerr)*1e3, mean(rstd)*1e3);
fprintf("              %.3f (%.3f 3-sigma) mm/s\n", mean(verr)*1e6, mean(vstd)*1e6);
fprintf("              %.3f (%.3f 3-sigma) ns\n", mean(berr)*1e9, mean(bstd)*1e9);
fprintf("Median error: %.3f (%.3f 3-sigma) m\n", median(rerr)*1e3, median(rstd)*1e3);
fprintf("              %.3f (%.3f 3-sigma) mm/s\n", median(verr)*1e6, median(vstd)*1e6);
fprintf("              %.3f (%.3f 3-sigma) ns\n", median(berr)*1e9, median(bstd)*1e9);

close all;

plotformat("IEEE", 1, "scaling", 2);
h4 = figure();

tplot = (ts - t0) / 3600;
colors = colororder();
subplot(3,1,1);
plot(tplot, rerr * 1e3);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 rstd * 1e3 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Position (m)");
title("State error over time, " + filter);

subplot(3,1,2);
plot(tplot, verr * 1e6);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 vstd * 1e6 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 75]);
ylabel("Velocity (mm/s)");

subplot(3,1,3);
plot(tplot, berr * 1e6);
hold on;
patch([ tplot(1) tplot tplot(end)]', [0 bstd * 1e6 0]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 0.5]);
legend(["Error", "3\sigma bound"], "Location", "best");
ylabel("Clock phase (\mus)");
xlabel("Time (hrs)");

% save and move plots to new folder
h4name = "stateerror_" + filter + "_" + string(yyyymmdd(datetime));
saveas(h4, h4name, "epsc");
saveas(h4, h4name, "png");
movefile(h4name + ".eps", "plots");
movefile(h4name + ".png", "plots");

% plot position and velocity error in RTN frame
figure();
subplot(3,1,1);
plot(tplot, err_rtn(1,:) * 1e3);
hold on;
xstd = 3 * sqrt(reshape(P_rtn(1,1,:), 1, ns)) * 1e3;
patch([tplot flip(tplot)]', [xstd -flip(xstd)]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Radial (m)");
title("Position error in RTN frame");

subplot(3,1,2);
plot(tplot, err_rtn(2,:) * 1e3);
hold on;
ystd = 3 * sqrt(reshape(P_rtn(2,2,:), 1, ns)) * 1e3;
patch([tplot flip(tplot)]', [ystd -flip(ystd)]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Tangential (m)");

subplot(3,1,3);
plot(tplot, err_rtn(3,:) * 1e3);
hold on;
zstd = 3 * sqrt(reshape(P_rtn(3,3,:), 1, ns)) * 1e3;
patch([tplot flip(tplot)]', [zstd -flip(zstd)]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Normal (m)");
xlabel("Time (hrs)");
legend(["Error", "3\sigma bound"], "Location", "best");

figure();
subplot(3,1,1);
plot(tplot, err_rtn(4,:) * 1e6);
hold on;
vxstd = 3 * sqrt(reshape(P_rtn(4,4,:), 1, ns)) * 1e6;
patch([tplot flip(tplot)]', [vxstd -flip(vxstd)]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Radial (mm/s)");
title("Velocity error in RTN frame");

subplot(3,1,2);
plot(tplot, err_rtn(5,:) * 1e6);
hold on;
vystd = 3 * sqrt(reshape(P_rtn(5,5,:), 1, ns)) * 1e6;
patch([tplot flip(tplot)]', [vystd -flip(vystd)]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Tangential (mm/s)");

subplot(3,1,3);
plot(tplot, err_rtn(6,:) * 1e6);
hold on;
vzstd = 3 * sqrt(reshape(P_rtn(6,6,:), 1, ns)) * 1e6;
patch([tplot flip(tplot)]', [vzstd -flip(vzstd)]', colors(2,:), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold off; box off;
axis([0 tplot(end) -inf inf]);
% axis([0 DAYS*24 0 50]);
ylabel("Normal (mm/s)");
xlabel("Time (hrs)");
legend(["Error", "3\sigma bound"], "Location", "best");

% plot satellite visibility
plotformat("IEEE", 0.4, "scaling", 2);
nvis = sum(visible, 1);
figure();
plot(tplot, nvis);
box off;
axis([0 tplot(end) -inf inf]);
xlabel("Time (hrs)"); ylabel("# of satellites");
title("GNSS satellites visible to user over simulation");
end

