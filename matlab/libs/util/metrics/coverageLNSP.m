function [CVG,minIdx,pts,covered] = coverageLNSP(dt,sats,SV,links,maxDOP)
%COVERAGELNSP Computes the coverage a constellation achieves in a given
%service volume as a proportion of the total evaluation time.
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
%    - SV    ; (1x1),[N/A] string of service volume selected -- options are
%              'SV1', 'SV2', or 'SV3'. see Lunar Relay SRD for details.
%    - links ; (1x1),[N/A] # of links required to close
%    - maxDOP; (1x1),[N/A] optional (default 6), max GDOP if links >=4
%   Output:
%    - CVG   ; (1x1),[N/A] proportion of time SV is covered, from 0 to 1
%    - minIdx; (1x1),[N/A] index of start of worst day
%    - pts   ; (3xl),[km] points coverage was evaluated at
%    - covered;(nxl),[N/A] data product for plotting

if     strcmp(SV, 'SV1')    % service volume 1
    lat = -80 * pi/180;     % rad, -80deg latitude or below
    alt = 125;              % km, max altitude
    l = 15000;              % initial # of points for eval
elseif strcmp(SV, 'SV2')    % service volume 2
    lat = -75 * pi/180;
    alt = 200;
    l = 7500;               
elseif strcmp(SV, 'SV3')    % service volume 3
    lat =  90 * pi/180;
    alt = 200;
    l = 120;
else                        % catch invalid arguments for SV
    throw(MException('coverageLNSP:invalidArgument', ...
          '%s is not a valid argument for SV. See documentation.', SV));
end

if nargin < 5, maxDOP = 6; end  % default maxDOP to 6 if not provided

% generate evaluation points and trim to SV latitude
[X,Y,Z] = mySphere(l);
accept = Z < sin(lat);
pts = [X; Y; Z];
pts = pts(:, accept);   % points on surface of unit sphere below lat

% scale to lunar distances and generate volume
R = 1736;               % km, polar radius of moon
alts = linspace(R, R + alt, 4);
pts = [pts*alts(1) pts*alts(2) pts*alts(3) pts*alts(4)];
% pts = getpoints(9);

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

[CVG, minIdx] = min(CVG);
end

