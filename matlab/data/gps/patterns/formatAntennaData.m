% formatAntennaData.m
% Author: Mark Hartigan
% Date  : September 18, 2024
% Description:
%    Organize the GPS block IIR and IIR-M antenna pattern data into a
%    format more machine-readable. Also average the azimuth patterns,
%    making it only dependent on elevation.

%% reset
clc, clear, close all;
addpath(genpath(pwd));

%% file operations
sheets = dir(fullfile("data/gps/patterns", "*.xlsx"));

for i=1:length(sheets)
    % get sheet name and extract info
    sheet = sheets(i);
    name = sheet.name;
    temp = split(name, {'-','.'});
    SVN = temp(1);                      % name in format 'SVNXX'
    L = temp(2);                        % transmit band, 'L1' or 'L2'
    newname = SVN + "-" + L + ".txt";   % filename for formatted data

    % data formatting; average phi angles and +/- theta
    data = importdata(name);
    data = data.data;
    rows = ~isnan(data(:,1));
    data = data(rows,:);                % remove header

    avgphi = mean(data(:,2:end), 2);
    mid = ceil(length(avgphi)/2);
    avgtht = mean([avgphi(mid:end) flipud(avgphi(1:mid))], 2);
    thetas = data(mid:end,1);

    % open and write to file
    id = fopen(sheet.folder + "/" + newname, 'w');
    fprintf(id, "theta (deg),directivity (dB)\n");
    for j=1:mid
        fprintf(id, "%f,%f\n", thetas(j), avgtht(j));
    end
    fclose(id);
end
    