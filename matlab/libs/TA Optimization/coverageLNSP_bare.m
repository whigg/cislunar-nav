function [CVG] = coverageLNSP_bare(dt,sats,pts)
%COVERAGELNSP_bare Computes the coverage a constellation achieves in a 
%given service volume as a proportion of the total evaluation time.
%   Per LunaNet Service Provider specifications (ESC-LCRNS-REQ-0090), the 
%   volumes are either SV1, SV2, or SV3 and the evaluation period should be
%   over 1 Earth month. The service volume is discretized into evaluation 
%   points. The number of links in view that define coverage and min GDOP 
%   for links >= 4 can be specified.
%
%   Inputs: (dims),[units]
%    - dt    ; (1x1),[s] seconds between steps; must be equal
%    - sats  ; (nx3xm),[km] positions of m satellites over n time steps --
%              time span covered should be 1 Earth month
%    - pts   ; (3xp),[km] p evaluation points
%   Output:
%    - CVG   ; (1x1),[N/A] proportion of time SV is covered, from 0 to 1

links = 4;
maxDOP = 6;

l = size(pts, 2);       % final number of points
n = size(sats, 1);      % number of time steps
covered = zeros(n, l);  % coverage of points at each time step?


for i=1:l           % iterate over every eval point
    [dop, nvis] = GDOP(pts(:,i), sats);

    covered(:,i) = nvis >= links;           % minimum links met
    if links >= 4                           % GDOP < max, if links 4+
        covered(:,i) = covered(:,i) .* (dop < maxDOP);
    end
end

% % proportion of time SV is covered (all points have coverage)
% CVG = sum(sum(covered, 2) == l) / n; 
% % proportion of time SV is covered (min point) over one month
% CVG = min(sum(covered, 1) / n);
% worst 24h period of coverage (in pct)
steps = floor(86400 / dt) + 1;              % time steps in 1 day
CVG = zeros(n-steps+1,1);
for i = 1:(n-steps+1)
    % coverage over day
    CVG(i) = min(sum(covered(i:(i+steps-1),:), 1) / steps);
end

CVG = min(CVG);
end

