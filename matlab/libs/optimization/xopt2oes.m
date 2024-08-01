function oes = xopt2oes(xopt)
%XOPT2OES converts output from ConOpt or ConOpt2 to an array of oe structs
%   Inputs:
%    - xopt; output from ConOpt or ConOpt2
arguments
    xopt    (1,:)   double
end

if length(xopt) == 14       % ConOpt used; orbits use same i, a
    oe1.i = xopt(1);
    oe1.a = xopt(2);
    oe1.e = frozenorbitfinder(oe1.i);
    oe1.w = pi/2;
    oe2 = oe1; oe3 = oe1; oe4 = oe1; oe5 = oe1; oe6 = oe1;
    
    oe1.RAAN = xopt(3); oe1.f = xopt(9);
    oe2.RAAN = xopt(4); oe2.f = xopt(10);
    oe3.RAAN = xopt(5); oe3.f = xopt(11);
    oe4.RAAN = xopt(6); oe4.f = xopt(12);
    oe5.RAAN = xopt(7); oe5.f = xopt(13);
    oe6.RAAN = xopt(8); oe6.f = xopt(14);
elseif length(xopt) == 24   % ConOpt 2 used; orbits constrained to same RAAN drift
    oe1.w = pi/2;
    oe2 = oe1; oe3 = oe1; oe4 = oe1; oe5 = oe1; oe6 = oe1;
    
    oe1.i = xopt(1); oe1.a = xopt(7);  oe1.e = frozenorbitfinder(oe1.i); oe1.RAAN = xopt(13); oe1.f = xopt(19);
    oe2.i = xopt(2); oe2.a = xopt(8);  oe2.e = frozenorbitfinder(oe2.i); oe2.RAAN = xopt(14); oe2.f = xopt(20);
    oe3.i = xopt(3); oe3.a = xopt(9);  oe3.e = frozenorbitfinder(oe3.i); oe3.RAAN = xopt(15); oe3.f = xopt(21);
    oe4.i = xopt(4); oe4.a = xopt(10); oe4.e = frozenorbitfinder(oe4.i); oe4.RAAN = xopt(16); oe4.f = xopt(22);
    oe5.i = xopt(5); oe5.a = xopt(11); oe5.e = frozenorbitfinder(oe5.i); oe5.RAAN = xopt(17); oe5.f = xopt(23);
    oe6.i = xopt(6); oe6.a = xopt(12); oe6.e = frozenorbitfinder(oe6.i); oe6.RAAN = xopt(18); oe6.f = xopt(24);
else                        % wrong size
    error("xopt2oes:InputError", ...
        "Xopt should either be 14 elements (ConOpt output) or 24 (ConOpt2 output).");
end

oes = [oe1 oe2 oe3 oe4 oe5 oe6];
end

