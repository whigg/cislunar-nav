classdef ConOpt2 < handle
    %ConOpt2 Constellation Optimizer. Creates a wrapper for the optimization
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
        T_OP2J  (6,6,2) double                          % frame transformation from MOON_OP to J2000 at t0
        T_J2ME  (3,3,:) double                          % frame transformation from J2000 to MOON_ME at ts
        flag    (1,1)   {mustBeInteger} = 2
    end
    
    methods
        function obj = ConOpt2(t0,p,varargin)
            %ConOpt2 Construct a ConOpt2 instance.
            %   Inputs:
            %    - t0; starting epoch (in seconds past J2000)
            %    - p; number of grid eval points for computing coverage
            %    - dt; optional name-value pair (default 3600), seconds between 
            %          time steps -- smaller yields greater accuracy, larger runs faster
            %    - n; optional name-value pair (default 5), # of satellites in constellation
            %    - opts; optional name-value pair, ODE45 propagation options

            narg = 2;   % number of required args

            obj.opts = odeset("RelTol", 1e-8, "AbsTol", 1e-9);

            if nargin > narg
                for i=1:2:nargin-narg
                    if strcmp(varargin{i}, "dt")
                        obj.dt = varargin{i+1};
                    elseif strcmp(varargin{i}, "n")
                        obj.n = varargin{i+1};
                    elseif strcmp(varargin{i}, "opts")
                        obj.opts = varargin{i+1};
                    elseif strcmp(varargin{i}, "n_sph")
                        obj.n_sph = varargin{i+1};
                    end
                end
            elseif nargin < narg
                error("ConOpt2:nargin", "Too few arguments.")
            end

            % get evaluation points
            obj.p = p;
            obj.pts = getpoints(p);
            % construct evaluation timespan (7 days, time step of dt)
            obj.ts = t0:obj.dt:t0 + 86400*7;
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
            obj.T_OP2J = zeros(6,6,2);
            obj.T_OP2J(:,:,1) = cspice_sxform('MOON_OP', 'J2000', t0);
            obj.T_OP2J(:,:,2) = cspice_sxform('MOON_OP', 'J2000', obj.ts(end));
            obj.T_J2ME = zeros(3,3,obj.m);
            for i=1:obj.m
                obj.T_J2ME(:,:,i) = cspice_pxform('J2000', 'MOON_ME', obj.ts(i));
            end
        end

        function [x,fval,exitflag,output,population,scores] = run(obj,R,S,opts)
            %RUN Executes MATLAB's built-in genetic algorithm on the nonsmooth 
            %objective function with linear and nonlinear constraints.

            nvars = 4*obj.n;
            fun = @(x) obj.coverage(x,R,S);
            
            % input vector x = [is,as,RAANs,TAs]

            % bounds:
            % 39.24 * pi/180 <= i <= 140.76 * pi/180
            % 2038 <= a <= 20000
            % 0 <= RAANs(i) <= 2*pi - 0.0001
            % 0 <= TAs(i) <= 2*pi - 0.0001
            lb = [39.24*pi/180*ones(1, obj.n)  2038*ones(1, obj.n)  zeros(1, 2*obj.n)];
            ub = [140.76*pi/180*ones(1, obj.n) 20000*ones(1, obj.n) (2*pi-0.0001)*ones(1, 2*obj.n)];

            % nonlinear constraints:
            % a*e - a <= -2038
            % Incorporate nonlinear constraint into objective instead using
            % exterior penalty methods
            % nonlcon = @ConOpt2.perilunecon;

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

            is    = x(1:obj.n);
            as    = x(obj.n+1:2*obj.n);
            RAANs = x(2*obj.n+1:3*obj.n);
            TAs   = x(3*obj.n+1:4*obj.n);
            y = 100;

            % handle boundary condition of periapsis here
            dOs = zeros(1, obj.n);
            for i=1:obj.n
                ei = frozenorbitfinder(is(i));
                rp = as(i)*(1 - ei);
                if rp < 2038, y = y + (2038 - rp)^2; end
                dOs(i) = ascendingnodedrift(as(i), ei, is(i));
            end
            % enfore boundary of drift rates being equal
            y = y + (32 * max(dOs - mean(dOs)) / mean(dOs))^2;

            [sats, fail] = obj.orbitgenerator(is,as,RAANs,TAs,R,S);
            warning('on','MATLAB:ode45:IntegrationTolNotMet');
            if fail     % at least one integration failed; invalid orbit
                return;
            end
               
            links = 4;
            maxDOP = 6;
            
            % coverage of points at each time step?
            covered = zeros(obj.m, obj.p);
            dops = zeros(size(covered));
            
            for j=1:obj.p               % iterate over every eval point
                [dop, nvis] = GDOP(obj.pts(:,j), sats);
            
                % minimum links met, meets GDOP requirements
                covered(:,j) = (nvis >= links) .* (dop <= maxDOP);
                dops(:,j) = (nvis >= links) .* (dop <= 20) .* dop;
            end
            
            % make obj combo of total % cvg (worst pt) + % worst-day (worst pt)
            % y = y - min(sum(covered, 1) / obj.m);
            % worst 24h period of coverage (in pct)
            steps = ceil(86400 / obj.dt);             % time steps in 1 day
            % win = ones(1, floor(steps * 5/24));     % time steps in 5-hour window
            % CVG = zeros(obj.m-steps+1,1);
            EVA = steps;
            % EVADOP = maxDOP;

            for j = 1:steps:(obj.m-steps+1)
                % coverage over day
                % [CVG(j), k] = min(sum(covered(j:(j+steps-1),:), 1) / steps);
                for k=1:obj.p
                    temp = sum(maxk(cellfun('length', split(char(covered(j:j+steps-1, k)'), char(0))), 2));
                    if temp < EVA
                        EVA = temp;
                        % daydop = dops(j:j+steps-1, k);
                        % EVADOP = 2*sum(daydop(daydop > maxDOP)) / steps;
                    end
                end
            end
            
            % TODO:
            %  - expand design to include orbits of same precession rate
            
            % flip negative so minimizing
            % y = y - 100 * min(CVG) - 100 * EVA / steps;
            y = y - 100 * EVA / steps;
            % penalize for having high DOP during the day
            % y = y + EVADOP^2;
        end

        function [sats,fail] = orbitgenerator(obj,is,as,RAANs,TAs,R,S)
            %ORBITGENERATOR creates position histories of satellites over the 
            %simulation time, given orbital elements.
            %   Inputs:
            %    - i; inclination of orbits (39.23 < i < 140.77 deg)
            %    - a; semimajor axis of orbit (perilune must be > 2038 km)
            %    - RAANs; vector of satellite right ascensions
            %    - TAs; vector of satellite starting true anomalies
            arguments
                obj     (1,1)   ConOpt2
                is      (1,:)   double {mustBeInclinations,mustBe1xN(obj,is)}
                as      (1,:)   double {mustBeSemimajorAxes,mustBe1xN(obj,as)}
                RAANs   (1,:)   double {mustBeAngle,mustBe1xN(obj,RAANs)}
                TAs     (1,:)   double {mustBeAngle,mustBe1xN(obj,TAs)}
                R       (6,6,2) double
                S       (3,3,:) double
            end

            % TODO: -change to keplerian orbit (no propagation)
            sats = zeros(obj.m, 3, obj.n);
            fail = false;

            for j=1:obj.n
                e = frozenorbitfinder(is(j));
                [r0,v0] = oe2rv(as(j),e,is(j),RAANs(j),pi/2,TAs(j),obj.moon.GM);
                x0 = R(:,:,1) * [r0; v0];
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
        % a*e - a + 2038 <= 0
        e = frozenorbitfinder(x(1));
        a = x(2);

        c = a*e - a + 2038;
        ceq = 0;
        end
    end
end

% custom validation
function mustBeInclinations(i)
%MUSTBEINCLINATION tests that the inclination provided is within ELFO bounds

for j=1:length(i)
    if i(j) >= 140.77 * pi/180 || i(j) <= 39.23 * pi/180
        error("coverage:inclinationOOB", "Inclination must be 39.23 < i < 140.77 deg.")
    end
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

function mustBeSemimajorAxes(a)
%MUSTBESEMIMAJORAXIS tests that value is a valid semimajor axis for ELFOs.
%   Update: modified, perilune constraint is now nonlinear and included in
%   optimizer.

for i=1:length(a)
    if a(i) < 2038 || a(i) > 20000
        error("coverage:semimajoraxis", "Semimajor axis must be 2038 <= a <= 20000 km.")
    end
end
end