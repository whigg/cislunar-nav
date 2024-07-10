function [sats,dt] = reshape_report(report_path)

% Read data from the script file
file = importdata(report_path);
data = file.data;

while size(data,1) > 2500

    mask = mod(1:size(data, 1), 2) == 0;
    data = data(mask, :);


end
num_sats = size(data,2)./3;


% Reshape XYZ coordinates into a 3D array
n = size(data,1);
sats = reshape(data, [n, 3, num_sats]);
dt = 86400.*31./size(sats,1); %time step in sec for report file over 1 month



end