function [nwin,minPt,minDay] = supportEVA(covered,dt)
%SUPPORTEVA Summary of this function goes here
%   Inputs: (dims),[units]
%    - covered; (nxl) availability (or not) of l points at n time steps
%    - dt     ; (1x1),[s] time between steps
%   Output: (dims),[units]
%    - nwin   ; (1x1),[N/A] number of 5-hour windows within each day in a
%               28-day period
%    - minPt  ; (1x1),[N/A] index of minimum point
%    - minDay ; (1x1),[N/A] start of worst day

n = size(covered, 1);           % number of time steps
l = size(covered, 2);           % number of eval points
day = ceil(86400 / dt);         % steps in 1 day
days = floor(n / day);          % days in provided history
fivehrs = floor(5*3600 / dt);   % steps in 5 hours

nwin = 1000;                    % initialize min number of windows
minPt = -1;
minDay = -1;
for i=1:days                    % iterate over every day
    for j=1:l                   % iterate over every point
        nwinj = 0;
        k = 1;
        while k <= day          % iterate over every step in day
            p = day*(i-1) + k;              % start of 5hr window
            q = min(day*i, p+fivehrs-1);    % end of 5hr window
            if sum(covered(p:q,j))/fivehrs == 1
                % if 5 hours from point is covered, add to nwinj
                nwinj = nwinj + 1;
                % skip to end of window
                k = k + fivehrs;
            else
                k = k + 1;
            end
        end

        if nwinj < nwin
            nwin = nwinj;
            minPt = j;
            minDay = i;
        end
    end
end
end

