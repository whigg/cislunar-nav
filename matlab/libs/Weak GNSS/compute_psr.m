function [psr, psrr] = compute_psr(t, x_true, sats, earth, max_a, err)
%COMPUTE_PSR Computes pseudoranges and -range-rates between moon-orbiting 
%spacecraft and visible GPS satellites.
%   Input:
%    - t; time history vector (seconds past J2000)
%    - x_true; true spacecraft trajectory
%    - sats; cell array of trajectory function handles for GPS satellites
%    - earth; planet struct {'GM': [km^3/s^2]; 'x': @(t) [km; km/s]}
%    - max_a; maximum GPS satellite beam width angle (rad)
%    - err; function handle that accepts time index input and outputs error
%           (in km) to apply to range calculation

n = size(x_true, 2);
m = length(sats);
psr = zeros(m, n);
psrr = zeros(m, n);

for i=1:n
    pos = x_true(1:3,i);        % spacecraft position
    for j=1:m
        sat = sats{j}(t(i));    % GPS satellite state
        r_s = sat(1:3) - pos;   % GPS sat position w.r.t. s/c
        r_e = earth.x(t(i));    % earth state vector
        r_e = r_e(1:3) - pos;   % earth position w.r.t. s/c
        % GPS sat velocity w.r.t. s/c
        v_s = sat(4:6) + r_e(4:6) - x_true(4:6,i);
        % convert velocity vector into range-rate
        dr_s = dot(v_s,r_s)/dot(r_s,r_s) * r_s;

        % compute earth-center / earth-tangent angle, earth-center / GPS
        % sat angle, and GPS-to-earth / GPS-to-s/c angle
        ang_e = asin(earth.R / norm(r_e));
        ang_s = acos(dot(r_s,r_e) / (norm(r_s)*norm(r_e)));
        ang_a = acos(dot(-r_s,r_e-r_s) / (norm(-r_s)*norm(r_e-r_s)));

        % satellite is NOT (further than earth AND within its view angle)
        % AND s/c is within the GPS satellite beam width
        if ~(norm(r_s) > norm(r_e) && ang_s < ang_e) && ang_a < max_a
            % km, distance b/n s/c & GPS sat + error
            psr(j,i) = norm(r_s) + err(i);
            psrr(j,i) = norm(dr_s);
        end
    end
end
end

