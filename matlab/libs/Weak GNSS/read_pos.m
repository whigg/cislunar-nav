function [t, sats] = read_pos(file)
%READ_POS Gets satellite position, velocity, and clock offset data for the
%duration of the specified .POS file, provided by the JPL GDGPS -- no idea
%what the actual file format is, just guessing
%   Input:
%    - file; relative path to the .pos file

% import and separate data, extracting relevant info
content = importdata(file);
data = content.data;
textdata = content.textdata;
names = unique(textdata(:,2));              % names of satellites
N = length(names);                          % number of satellites in file
t = unique(data(:,1));                      % unique timestamps
M = length(t);                              % number of time steps
satidx = dictionary(string(names)', 1:N);   % dict to relate name and idx

sats = repmat(struct('dt', zeros(1,M), 'x', zeros(6,M)), N, 1);
% sats = repmat(struct('dt', zeros(1,M), 'x', zeros(3,M)), N, 1);


for k=1:size(data, 1)                       % iterate over row data
    i = find(t == data(k,1), 1);            % time index
    j = satidx(textdata{k,2});              % sat index from name
    % get ECEF to ECI rotation matrix after converting past J2000 to JD
    R = ECEF2ECI(t(i) / 86400 + 2451545);
    sats(j).x(:,i) = data(k,3:8)';          % assign pos/vel vectors
    % sats(j).x(:,i) = R * data(k,3:5)';
end
end