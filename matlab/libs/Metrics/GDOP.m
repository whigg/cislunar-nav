function [GDOP,nvis] = GDOP(pos,sats)
%GDOP Compute the geometric dilution of precision of the current geometry.
%   Inputs: (dims),[units]
%    - pos ; (3x1),[km] position at which to compute GDOP
%    - sats; (nx3xm),[km] positions of m satellites over n time steps
%   Outputs:
%    - GDOP; (nx1),[N/A] GDOP per time step
%    - nvis; (nxm),[N/A] number of visible satellites at each step

visible = visibleSats(pos, sats);
GDOP = NaN(size(visible, 1), 1);
nvis = sum(visible ~= 0, 2);

% temporarily suppress nearly singular matrix messages (from H)
warning('off', 'MATLAB:nearlySingularMatrix');
for i=1:length(GDOP)
    if nvis(i) < 4
        GDOP(i) = inf;
    else
        vis = visible(i, visible(i,:) ~= 0);
        G = ones(nvis(i), 4);
        
        for j = 1:nvis(i)
            ik = sats(i,:,vis(j)) - pos';
            G(j,1:3) = ik / norm(ik);
        end
    
        H = inv(G'*G);
        % take absolute value because sometimes diagonals are negative
        GDOP(i) = sqrt(sum(abs(diag(H))));
    end
end
warning('on', 'MATLAB:nearlySingularMatrix');
end

