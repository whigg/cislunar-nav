% AnalyzeConstellation.m
% Author: Mark Hartigan
% Date  : May 16, 2024
% Description:
%    Display full constellation and its statistics

DAYS = 90;
step = 86400*DAYS / 7200;
prop = LunarPropagator(START, xopt, 32, 3);
[ts,xs] = prop.run(86400 * DAYS, 7200, 'MOON_ME');
prop.plotlastorbits('MOON_OP');

%% get coverage statistics
sats = permute(xs, [2,1,3]);
[pct, minIdx, pts, covered] = coverageLNSP(step, sats, 'SV2', 4);
plotDayCoverage(pts, covered, minIdx, step);
plotDayCoverage(pts, covered, minIdx, step, true);
fprintf("Minimum %% coverage of service volume over any 24h period: %.2f%%\n", pct*100);
fprintf("\t(starting at day %.2f)\n", (ts(minIdx)-t0)/(ts(end)-t0) * DAYS);

%% EVA coverage
[nwin, minPt, minDay] = supportEVA(covered, step);
day = ceil(86400 / step);       % steps in 1 day
p = day*(minDay - 1) + 1;
q = day*minDay;
% coverage of worst point during worst interval
plotGDOP(pts(:,minPt), sats(p:q,:,:), ts(1:day), 10);
% coverage of LSP for entire duration
plotGDOP([0;0;-1738], sats, ts, 10);
