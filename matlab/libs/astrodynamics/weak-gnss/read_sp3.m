function [t, sats] = read_sp3(file)
%READ_SP3 Gets satellite position and clock offset data for duration of the
%specified SP3 file, provided by the IAC (ONLY WORKS WITH FILES BEGINNING
%WITH IAC0MGXFIN).
%   Input:
%    - file; relative path to the .sp3 file

lines = readlines(file);
lines = lines(29:end);
N = 120;                                    % number of satellites in file
M = floor(length(lines) / (N+1));           % number of timestamps

t = zeros(1,M);
sats = repmat(struct('dt', zeros(1,M), 'x', zeros(3,M)), N, 1);

i = 0;
k = 0;
for j=1:numel(lines)
    line = lines(j);

    if startsWith(line, "EOF"), return;     % leave if EOF
    elseif startsWith(line, "*")            % timestamp of data
        i = i + 1;                          % increment datapoint counter
        k = 0;                              % reset counter
        data = split(line);
        d    = str2double(data(2:end));
        t(i) = juliandate(datetime(d(1), d(2), d(3), d(4), d(5), d(6)));
        t(i) = t(i) - 18 / 86400;           % convert GPST to UTC
        R = ECEF2ECI(t(i));                 % ECEF to ECI rotation matrix
    else
        k = k + 1;                          % next sat
        data = split(line);
        d    = str2double(data(2:5));
        sats(k).dt(i)  = d(4) * 1e-9;       % convert ns to s
        sats(k).x(:,i) = R * d(1:3);        % rotate to ECI
    end
end

end