function [] = UNE(data)
%UNE generates a plot of the user navigation error for given data.
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

% CSAC information
a_hi = 9e-10 / (86400 * 30);
E = @(t, a) 1/2 * a * t.^2;
c = 299792458;

% dilution of precision @ south pole
une1 = zeros(size(t));
une2 = zeros(size(t));
for k=1:m
    H = diag(computeDOP([0; 0; -r], sats, k));
    dop = sqrt(abs(sum(H(1:3))));
    uere1 = sqrt((E(mod(t(k),13146), a_hi)*c)^2 + 1.96^2 + 10^2);
    uere2 = sqrt((E(mod(t(k),13146), a_hi)*c)^2 + 1.96^2 + 65^2);
    une1(k) = dop * uere1;
    une2(k) = dop * uere2;
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