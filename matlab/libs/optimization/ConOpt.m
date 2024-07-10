classdef ConOpt < handle
    %CONOPT Constellation Optimizer. Creates a wrapper for the optimization
    %of lunar navigation constellations, using MATLAB's underlying genetic
    %algorithm ga().
    
    properties
        ts      (1,:)   double                          % simulation times
        dt      (1,1)   double {mustBePositive} = 3600  % simulation time step
        pts     (3,:)   double                          % evaluation points
        n       (1,1)   {mustBePositive,mustBeInteger} = 5  % number of constellation satellites
        m       (1,1)   {mustBePositive,mustBeInteger} = 1  % number of time steps
        p       (1,1)   {mustBePositive,mustBeInteger} = 1  % number of eval points
        n_sph   (1,1)   {mustBePositive,mustBeInteger} = 2  % number of spherical harmonics
        opts    (1,1)   struct                          % ODE45 propagation options
        moon    (1,1)   struct                          % planetary info for moon
        earth   (1,1)   struct                          % planetary info for earth
        T_OP2J  (6,6)   double                          % frame transformation from MOON_OP to J2000 at t0
        T_J2ME  (3,3,:) double                          % frame transformation from J2000 to MOON_ME at ts
        flag    (1,1)   {mustBeInteger} = 1
    end
    
    methods
        function obj = ConOpt(t0,p,varargin)
            %CONOPT Construct a ConOpt instance.
            %   Inputs:
            %    - t0; starting epoch (in seconds past J2000)
            %    - p; number of grid eval points for computing coverage
            %    - dt; optional name-value pair (default 3600), seconds between 
            %          time steps -- smaller yields greater accuracy, larger runs faster
            %    - n; optional name-value pair (default 5), # of satellites in constellation
            %    - opts; optional name-value pair, ODE45 propagation options

            narg = 2;   % number of required args

            opts = odeset("RelTol", 1e-8, "AbsTol", 1e-9);

            if nargin > narg
                for i=1:2:nargin-narg
                    if strcmp(varargin{i}, "dt")
                        obj.dt = varargin{i+1};
                    elseif strcmp(varargin{i}, "n")
                        obj.n = varargin{i+1};
                    elseif strcmp(varargin{i}, "opts")
                        obj.opts = opts;
                    elseif strcmp(varargin{i}, "n_sph")
                        obj.n_sph = varargin{i+1};
                    end
                end
            elseif nargin < narg
                error("EKF:nargin", "Too few arguments.")
            end

            % get evaluation points
            obj.p = p;
            obj.pts = getpoints(p);
            % construct evaluation timespan (30 days, time step of dt)
            obj.ts = t0:obj.dt:t0 + 86400*30;
            obj.m = length(obj.ts);

            % planetary information from SPICE
            % earth information
            earth.GM = cspice_bodvrd('EARTH', 'GM', 1);
            % convert to spline to avoid SPICE calls in dynamics
            ex = cspice_spkpos('EARTH', obj.ts, 'J2000', 'NONE', 'MOON');
            ppex = spline(obj.ts, ex);
            earth.x = @(tau) ppval(ppex, tau);
            
            % moon information
            moon.GM = cspice_bodvrd('MOON', 'GM', 1);
            [R,C,S] = cofloader(userpath + "/astrodynamics/LP165P.cof");
            moon.x = @(~) [0;0;0];
            moon.R = R * 1e-3;      % convert from m to km
            moon.C = C;             % store in moon struct for orbitaldynamics
            moon.S = S;             % store in moon struct for orbitaldynamics
            moon.frame = 'J2000';   % J2000 means don't use coefficients
            obj.earth = earth;
            obj.moon = moon;

            % necessary storage of transformations since SPICE cannot be used in
            % parallel applications
            obj.T_OP2J = cspice_sxform('MOON_OP', 'J2000', t0);
            obj.T_J2ME = zeros(3,3,obj.m);
            for i=1:obj.m
                obj.T_J2ME(:,:,i) = cspice_pxform('J2000', 'MOON_ME', obj.ts(i));
            end
        end

        function [x,fval,exitflag,output,population,scores] = run(obj,R,S,opts)
            %RUN Executes MATLAB's built-in genetic algorithm on the nonsmooth 
            %objective function with linear and nonlinear constraints.

            nvars = 2*obj.n + 2;
            fun = @(x) obj.coverage(x,R,S);
            
            % input vector x = [i,a,RAANs,TAs]

            % bounds:
            % 39.24 * pi/180 <= i <= 140.76 * pi/180
            % 2138 <= a <= 20000
            % 0 <= RAANs(i) <= 2*pi - 0.0001
            % 0 <= TAs(i) <= 2*pi - 0.0001
            lb = [39.24*pi/180 2138 zeros(1, 2*obj.n)];
            ub = [140.76*pi/180 20000 (2*pi-0.0001)*ones(1, 2*obj.n)];

            % nonlinear constraints:
            % a*e - a <= -2138
            % Incorporate nonlinear constraint into objective instead using
            % exterior penalty methods
            % nonlcon = @ConOpt.perilunecon;

            % x0 = load('data/RUN5 (20p).mat');
            % x0.xopt = [1 15000 0 0 0 0 0 0 0 0 0 0 0 0];
            % x = patternsearch(fun,x0.xopt,[],[],[],[],lb,ub,nonlcon,opts);
            % x = ga(fun,nvars,[],[],[],[],lb,ub,nonlcon,opts);
            % suppress integration errors from ODE45 and run
            warning('off','MATLAB:ode45:IntegrationTolNotMet');
            [x,fval,exitflag,output,population,scores] = ga(fun,nvars,[],[],[],[],lb,ub,[],opts);
            % unsuppress integration errors from ODE45
            warning('on','MATLAB:ode45:IntegrationTolNotMet');
        end
        
        function y = coverage(obj,x,R,S)
            %COVERAGE Computes the coverage a constellation achieves in a 
            %given service volume as a proportion of the total evaluation time.
            %   Per LunaNet Service Provider specifications (ESC-LCRNS-REQ-0090), the 
            %   volumes are either SV1, SV2, or SV3 and the evaluation period should be
            %   over 1 Earth month. The service volume is discretized into evaluation 
            %   points. The number of links in view that define coverage and min GDOP 
            %   for links >= 4 can be specified.
            %
            %   Inputs:
            %    - x; input vector [i,a,RAANs,TAs] (1x2*n+2)
            %    - R; rotation matrices to convert from MOON_OP to MOON_ME
            %   Output:
            %    - y; proportion of time SV is covered, from 0 to 1
            warning('off','MATLAB:ode45:IntegrationTolNotMet');

            i = x(1);
            a = x(2);
            RAANs = x(3:obj.n+2);
            TAs = x(obj.n+3:2*obj.n+2);
            y = 100;

            % handle boundary condition here by using pseudo-objective func
            e = frozenorbitfinder(i);
            rp = a*(1-e);
            if rp < 2138
                y = y + (2138 - rp)^2;
                return
            end
            [sats, fail] = obj.orbitgenerator(i,a,RAANs,TAs,R,S);
            if fail     % at least one integration failed; invalid orbit
                return;
            end
               
            links = 4;
            maxDOP = 6;
            
            % coverage of points at each time step?
            covered = zeros(obj.m, obj.p);
            
            for j=1:obj.p               % iterate over every eval point
                [dop, nvis] = GDOP(obj.pts(:,j), sats);
            
                % minimum links met, meets GDOP requirements
                covered(:,j) = (nvis >= links) .* (dop <= maxDOP);
            end
            
            % make obj combo of total % cvg (worst pt) + % worst-day (worst pt)
            % y = y - min(sum(covered, 1) / obj.m);
            % worst 24h period of coverage (in pct)
            steps = floor(86400 / obj.dt) + 1;      % time steps in 1 day
            % win = ones(1, floor(steps * 5/24));     % time steps in 5-hour window
            % CVG = zeros(obj.m-steps+1,1);
            EVA = steps;

            for j = 1:steps:(obj.m-steps+1)
                % coverage over day
                % [CVG(j), k] = min(sum(covered(j:(j+steps-1),:), 1) / steps);
                for k=1:obj.p
                    temp = sum(maxk(cellfun('length', split(char(covered(j:j+steps-1, k)'), char(0))), 2));
                    if temp < EVA, EVA = temp; end
                end
            end

            % TODO:
            %  - expand design to include orbits of same precession rate
            
            % flip negative so minimizing
            % y = y - 100 * min(CVG) - 100 * EVA / steps;
            y = y - 100 * EVA / steps;
        end

        function [sats,fail] = orbitgenerator(obj,i,a,RAANs,TAs,R,S)
            %ORBITGENERATOR creates position histories of satellites over the 
            %simulation time, given orbital elements.
            %   Inputs:
            %    - i; inclination of orbits (39.23 < i < 140.77 deg)
            %    - a; semimajor axis of orbit (perilune must be > 2138 km)
            %    - RAANs; vector of satellite right ascensions
            %    - TAs; vector of satellite starting true anomalies
            arguments
                obj     (1,1)   ConOpt
                i       (1,1)   double {mustBeInclination}
                a       (1,1)   double {mustBeSemimajorAxis}
                RAANs   (1,:)   double {mustBeAngle,mustBe1xN(obj,RAANs)}
                TAs     (1,:)   double {mustBeAngle,mustBe1xN(obj,TAs)}
                R       (6,6)   double
                S       (3,3,:) double
            end

            % TODO: -change to keplerian orbit (no propagation)
            e = frozenorbitfinder(i);
            sats = zeros(obj.m, 3, obj.n);
            fail = false;

            for j=1:obj.n
                [r0,v0] = oe2rv(a,e,i,RAANs(j),pi/2,TAs(j),obj.moon.GM);
                x0 = R * [r0; v0];
                [~,X] = ode45(@(t,x) orbitaldynamics(t,x,obj.moon,obj.n_sph,obj.earth), ...
                              obj.ts, x0, obj.opts);

                % exit from loop if propagation failed
                if size(X,1) ~= obj.m
                    fail = true;
                    return
                end

                X = X(:,1:3)';
                for k=1:obj.m
                    sats(k,:,j) = S(:,:,k) * X(:,k);
                end
            end
        end
    end

    methods(Static)
        function [c,ceq] = perilunecon(x)
        %PERILUNECON Provides the nonlinear perilune constraints for MATLAB's 
        %ga().
        %   Inputs:
        %    - x; input vector [i,a,RAANs,TAs] (1x2*n+2)

        % nonlinear constraints:
        % a*e - a + 2138 <= 0
        e = frozenorbitfinder(x(1));
        a = x(2);

        c = a*e - a + 2138;
        ceq = 0;
        end
    end
end

% custom validation
function mustBeInclination(i)
%MUSTBEINCLINATION tests that the inclination provided is within ELFO bounds

if i >= 140.77 * pi/180 || i <= 39.23 * pi/180
    error("coverage:inclinationOOB", "Inclination must be 39.23 < i < 140.77 deg.")
end
end

function mustBeAngle(a)
%MUSTBEANGLE tests that the value provided is between 0 and 2pi.

if any(a >= 2*pi) || any(a < 0)
    error("coverage:angleOOB", "Angle must be 0 <= a < 2*pi.")
end
end

function mustBe1xN(obj,vec)
%MUSTBE1XN tests that vector provided is (1,n).

if size(vec,1) ~= 1 || size(vec,2) ~= obj.n
    error("coverage:vecSize", "Vector must be of size (1,n).")
end
end

function mustBeSemimajorAxis(a)
%MUSTBESEMIMAJORAXIS tests that value is a valid semimajor axis for ELFOs.
%   Update: modified, perilune constraint is now nonlinear and included in
%   optimizer.

if a < 2138 || a > 20000
    error("coverage:semimajoraxis", "Semimajor axis must be 2138 <= a <= 20000 km.")
end
end