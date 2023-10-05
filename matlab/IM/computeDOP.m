function [H] = computeDOP(pos, sats, time)
%COMPUTEDOP Compute the dilution of precision of the current geometry.
%   Inputs:
%    - pos;  inertial user position
%    - sats; 3D satellite arrays
%    - time; index of sats to examine

vis = visibleSats(pos, sats, time);
num = length(vis);
G = ones(num, 4);

for i = 1:num
    ik = sats(time,:,vis(i)) - pos';
    G(i,1:3) = ik / norm(ik);
end

H = inv(G'*G);

end

