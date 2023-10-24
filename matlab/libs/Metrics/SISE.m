function ERR = SISE(T,EPH,CLK)
%SISE Computes and returns the signal-in-space error between the broadcast
%and true PVT for LunaNet Service Provider nodes. It establishes an upper
%bound of error in the user-satellite direction.
%   Inputs: (dims),[units]
%    - T  ; (nx1),[s] time steps to evaluate SISE at (cal. at 0s)
%    - EPH; (1x1),[m] optional, 3-sigma error from ephemeris uncertainty
%    - CLK; (1x1),[m] optional, clock update threshold; assumes using the
%           MicroSemi Space CSAC clock
%   Output:
%    - ERR; (1x1),[m] 1D signal-in-space error given the inputs

c = 299792458;                          % m/s, speed of light
if nargin < 3
    CLK = (30e-9 * c) / 1.96;           % UTCOE < 30ns 95%, from GPS

    if nargin < 2
        EPH = 10 / 3;                   % 10m 3-sigma by default
    end
end

% CSAC information
a_hi = 9e-10 / (86400 * 30);
E = @(t, a) 1/2 * a * t.^2;
tmax = sqrt(2 * CLK / c / a_hi);

var = (E(mod(T, tmax), a_hi) .* c).^2;  % clock update error
var = var + (EPH/3)^2;                  % ephemeris error estimate
var = var + 1^2;                        % multipath estimate

ERR = sqrt(var);                        % SISE is std of variance
end