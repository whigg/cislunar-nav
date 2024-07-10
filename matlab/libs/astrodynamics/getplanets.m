function bod = getplanets(varargin)
%GETPLANETS Returns structs containing information about various planets in
%the order they are input.
%   Inputs:
%    - name(s); name(s) of planet(s) desired -- "MOON", "SUN", "EARTH", and
%               "JUPITER" currently supported

bod = struct([]);
for i=1:length(varargin)
    if     strcmpi(varargin{i}, "MOON")
        bod(i).GM = cspice_bodvrd('MOON', 'GM', 1);
        bod(i).x = @(~) [0;0;0];
        bod(i).R = max(cspice_bodvrd('MOON', 'RADII', 3));

    elseif strcmpi(varargin{i}, "EARTH")
        bod(i).GM = cspice_bodvrd('EARTH', 'GM', 1);
        bod(i).x  = @(tau) cspice_spkpos('EARTH', tau, 'J2000', 'NONE', 'MOON');
        bod(i).R  = max(cspice_bodvrd('EARTH', 'RADII', 3));

    elseif strcmpi(varargin{i}, "SUN")
        bod(i).GM = cspice_bodvrd('SUN', 'GM', 1);
        bod(i).x  = @(tau) cspice_spkpos('SUN', tau, 'J2000', 'NONE', 'MOON');
        bod(i).R  = max(cspice_bodvrd('SUN', 'RADII', 3));

    elseif strcmpi(varargin{i}, "JUPITER")
        bod(i).GM = cspice_bodvrd('JUPITER BARYCENTER', 'GM', 1);
        bod(i).x  = @(tau) cspice_spkpos('JUPITER BARYCENTER', tau, 'J2000', 'NONE', 'MOON');
        bod(i).R  = max(cspice_bodvrd('SUN', 'RADII', 3));

    else
        error("geplanets:planetNotFound", ...
            "Provided planet " + varargin{i} + " not supported.");
    end
end
end

