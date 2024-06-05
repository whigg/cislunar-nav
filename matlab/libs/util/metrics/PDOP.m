function [PDOP,visible] = PDOP(pos,sats)
%PDOP Compute the positional dilution of precision of the current geometry.
%   Inputs:
%    - pos ; (3x1),[km] position at which to compute PDOP
%    - sats; (nx3xm),[km] positions of m satellites over n time steps
%   Outputs:
%    - PDOP; (nx1),[N/A] PDOP per time step
%    - visible; (nxm),[N/A] indices of visible satellites are nonzero

visible = visibleSats(pos, sats);
PDOP = NaN(size(visible, 1), 1);
nvis = sum(visible ~= 0, 2);

% temporarily suppress nearly singular matrix messages (from H)
warning('off', 'MATLAB:nearlySingularMatrix');
for i=1:length(PDOP)
    if nvis(i) < 4
        PDOP(i) = inf;
    else
        vis = visible(i, visible(i,:) ~= 0);
        G = ones(nvis(i), 4);
        
        for j = 1:nvis(i)
            ik = sats(i,:,vis(j)) - pos';
            G(j,1:3) = ik / norm(ik);
        end
    
        H = inv(G'*G);
        % take absolute value because sometimes diagonals are negative
        PDOP(i) = sqrt(sum(abs(diag(H))));
    end
end
warning('on', 'MATLAB:nearlySingularMatrix');
end

