function [H] = computeDOP_bare(r, sats, time)
%COMPUTEDOP_BARE Compute the positional dilution of precision of the 
%current geometry.
%   Stripped-down version of computeDOP() that only considers a user on the
%   lunar south pole.
%   Inputs:
%    - r;  radius of moon
%    - sats; 3D satellite arrays
%    - time; index of sats to examine

vis = visibleSats_bare(r, sats, time);
num = length(vis);
G = ones(num, 4);
pos = [0; 0; -r];

for i = 1:num
    ik = sats(time,:,vis(i)) - pos';
    G(i,1:3) = ik / norm(ik);
end

H = diag(inv(G'*G));            % take only the diagonal elements
H = sqrt(abs(sum(H(1:3))));     % find the std. dev. magnitude
end

