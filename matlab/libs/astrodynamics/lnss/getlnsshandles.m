function [fLnss, namesLnss] = getlnsshandles(ref)
%GETLNSSHANDLES Returns function handles for each satellite in generated
%constellation; handle accepts time(s) in seconds past J2000 and returns
%[pos (km); vel (km/s)] state(s).
%   Input:
%    - ref; reference frame for SPICE, e.g. 'J2000' or 'MOON_ME'
%   Output:
%    - fLnss; cell array of fcn handles corresponding to namesLnss sats
%    - namesLnss; SPICE names of LNSS satellites included

% LNSS SPICE kernels generated using GMAT data and MKSPK
%
%     https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/ug/mkspk.html
%
% Data provided from 10 Oct 2009 to 9 Nov 2009
cspice_furnsh(strcat(userpath,'/kernels/LNSS/mk/lnss.tm'));

% load LNSS trajectories; SPICE names of LNSS satellites & # of LNSS satellites to include
namesLnss = importdata(strcat(userpath,'/kernels/LNSS/LNSS_sats.txt'));
nLnss = length(namesLnss);
fLnss = cell(nLnss, 1);
for i=1:nLnss
    fLnss{i} = @(tau,ref) cspice_spkezr(num2str(namesLnss(i)), tau, ref, 'NONE', 'MOON');
end
end


