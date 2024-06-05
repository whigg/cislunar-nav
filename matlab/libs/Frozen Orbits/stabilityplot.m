function stabilityplot(t,x,mu,frame)
%STABILITYPLOT Plots eccentricity vs argument of periapsis on a polar plot,
%given state of a satellite over time.
%   Input:
%    - x; time history of satellite state (position, velocity)
%    - mu; gravitational parameter of moon
%    - frame; frame x is supplied in; will transform if not MOON_OP
arguments
    t (1,:) double
    x (6,:) double
    mu (1,1) double
    frame (1,:) char
end

n = size(x, 2);
es = zeros(1,n);
ws = zeros(1,n);

for i=1:n
    x0 = cspice_sxform(frame, 'MOON_OP', t(i)) * x(:,i);
    [~,e,~,~,w,~] = rv2oe(x0(1:3), x0(4:6), mu);
    es(i) = e;  ws(i) = w;
end

figure();
polarplot(ws, es);
title("Eccentricity vs. argument of periapsis");
end

