function x0 = gaoutputtostates(xopt, n, t0)
%GAOUTPUTTOSTATES Converts the output of the genetic algorithm
%constellation optimization to starting states in J2000
%   Inputs:
%    - xopt; optimized keplerian orbital elements in MOON_OP frame
%    - n; number of satellites
%    - t0; seconds past J2000, starting time of OEs
arguments
    xopt (1,:) double
    n    (1,1) {mustBePositive,mustBeInteger}
    t0   (1,1) double
end

x0 = zeros(6, n);
w = pi/2;
es = zeros(1, n);
mu = cspice_bodvrd('MOON', 'GM', 1);

if length(xopt) ~= 4*n
    is(:) = xopt(1) * ones(1, n);
    as(:) = xopt(2) * ones(1, n);
    RAANs = xopt(3:n+2);
    fs    = xopt(n+3:2*n+2);
    es(:) = frozenorbitfinder(xopt(1));
else
    is    = xopt(1:n);
    as    = xopt(n+1:2*n);
    RAANs = xopt(2*n+1:3*n);
    fs    = xopt(3*n+1:4*n);
    for i=1:n, es(i) = frozenorbitfinder(is(i));  end
end

for i=1:n
    [r,v] = oe2rv(as(i), es(i), is(i), RAANs(i), w, fs(i), mu);
    x0(:,i) = cspice_sxform('MOON_OP', 'J2000', t0) * [r; v];
end
end

