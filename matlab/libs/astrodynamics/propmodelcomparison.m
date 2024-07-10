function propmodelcomparison(ts,xs,bodies)
%PROPMODELCOMPARISON Generates a plot comparing the efficacy of different
%propagation models for lunar orbits.
%   Inputs:
%    - ts; time history to compare at
%    - xs; state history at times ts (truth)
%    - bodies; info about primary and third bodies in system (first in list
%              is primary, rest are third-body)
arguments
    ts      (1,:)   double
    xs      (6,:)   double
    bodies  (1,:)   struct
end

ns = length(ts);
x0 = xs(:,1);
t0 = ts(1);
opts = odeset("RelTol", 1e-9, "AbsTol", 1e-11);
moon = bodies(1);
sec = bodies(2:end);

[~,X] = ode45(@(t,x) orbitaldynamics(t,x,moon,1,[]), ts, x0, opts);
xKep = X';
[~,X] = ode45(@(t,x) orbitaldynamics(t,x,moon,1,sec), ts, x0, opts);
xThird = X';
[~,X] = ode45(@(t,x) orbitaldynamics(t,x,moon,10,sec), ts, x0, opts);
xSph1 = X';
[~,X] = ode45(@(t,x) orbitaldynamics(t,x,moon,100,sec), ts, x0, opts);
xSph2 = X';

plotLunarOrbit(ts, xKep', 'J2000', "Keplerian");
plotLunarOrbit(ts, xThird', 'J2000', "+ 3rd body");

diffKep = sqrt(sum((xKep(1:3,:) - xs(1:3,:)).^2, 1)) * 1e3;
diffThird = sqrt(sum((xThird(1:3,:) - xs(1:3,:)).^2, 1)) * 1e3;
diffSph1 = sqrt(sum((xSph1(1:3,:) - xs(1:3,:)).^2, 1)) * 1e3;
diffSph2 = sqrt(sum((xSph2(1:3,:) - xs(1:3,:)).^2, 1)) * 1e3;

plotformat("IEEE", 0.5, "scaling", 2);
figure();
semilogy((ts - t0) / 3600, diffKep, "LineWidth", 1.5);
hold on;
semilogy((ts - t0) / 3600, diffThird, "--", "LineWidth", 1.5, "Marker", "x", "MarkerSize", 8, "MarkerIndices", 1:100:ns);
semilogy((ts - t0) / 3600, diffSph1, "-.", "LineWidth", 1.5);
semilogy((ts - t0) / 3600, diffSph2, ":", "LineWidth", 1.5);
hold off; grid on;
axis([0 (ts(end)-t0)/3600 -inf inf]);
xlabel("Time (hrs)");
ylabel("RSS Position Error (m)");
title("Propagation error of different models");
legend(["Keplerian", "+ third body", "+ SH (n,m=10)", "+ SH (n,m=100)"], "Location", "best");
end

