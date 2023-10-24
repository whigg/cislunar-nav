function [] = UNE(data)
%UNE generates a plot of the user navigation error at lunar south pole for 
%given data.
% Input:
%  - data; time history of satellite positions (m x (1 + 3*n)) [s, km...]

t = data(:,1);                  % time history
m = length(t);                  % number of data points
n = (size(data,2) - 1) / 3;     % number of satellites
r = 1737.4; % km
sats = zeros(m, 3, n);
for num = 1:n
    sats(:,:,num) = data(:,2 + 3*(num-1):1 + 3*num);
end

t = t - t(1);                   % set t0 = 0

sise1 = SISE(t, 10);
sise2 = SISE(t, 100);

% dilution of precision @ south pole
une1 = zeros(size(t));
une2 = zeros(size(t));
for k=1:m
    H = diag(computeDOP([0; 0; -r], sats, k));
    dop = sqrt(abs(sum(H(1:3))));
    une1(k) = dop * sise1(k);
    une2(k) = dop * sise2(k);
end

figure();
plot(t * 24 / 86400, une1, 'LineWidth', 1.5);
hold on;
plot(t * 24 / 86400, une2, '--', 'LineWidth', 1.5);
axis([0 (t(end) * 24 / 86400) 0 300]);
grid on;
xlabel("Time (hrs)"); ylabel("95% RSS Position Error (m)");
title("User Navigation Error at the Lunar South Pole");
set(gcf, 'position', [500, 250, 900, 350]);
legend(["OD Uncertainty = 65m", "OD Uncertainty = 10m"]);

end