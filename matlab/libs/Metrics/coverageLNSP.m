function [CVG] = coverageLNSP(sats,SV,links,maxDOP)
%COVERAGELNSP Computes the coverage a constellation achieves in a given
%service volume as a proportion of the total evaluation time.
%   Per LunaNet Service Provider specifications (ESC-LCRNS-REQ-0090), the 
%   volumes are either SV1, SV2, or SV3 and the evaluation period should be
%   over 1 Earth month. The service volume is discretized into evaluation 
%   points. The number of links in view that define coverage and min GDOP 
%   for links >= 4 can be specified.
%
%   Inputs: (dims),[units]
%    - sats  ; (nx3xm),[km] positions of m satellites over n time steps --
%              time span covered should be 1 Earth month
%    - SV    ; (1x1),[N/A] string of service volume selected -- options are
%              'SV1', 'SV2', or 'SV3'. see Lunar Relay SRD for details.
%    - links ; (1x1),[N/A] # of links required to close
%    - maxDOP; (1x1),[N/A] optional (default 6), max GDOP if links >=4
%   Output:
%    - CVG   ; (1x1),[N/A] proportion of time SV is covered, from 0 to 1

if     strcmp(SV, 'SV1')    % service volume 1

elseif strcmp(SV, 'SV2')    % service volume 2

elseif strcmp(SV, 'SV3')    % service volume 3

else                        % catch invalid arguments for SV
    throw(MException('coverageLNSP:invalidArgument', ...
          '%s is not a valid argument for SV. See documentation.', SV));
end

if nargin < 4, maxDOP = 6; end  % default maxDOP to 6 if not provided

CVG = 0;
end

