classdef LunarPropagator < handle
    %LUNARPROPAGATOR Generic propagation class for lunar satellites.
    %   Can propagate lunar orbits for various lengths of time and starting
    %   conditions. Mainly created to reduce repetition / verbosity of
    %   scripts.
    
    properties
        % struct containing info about the moon
        moon    (1,1)   struct
        % array of structs containing info about additional planets to consider
        sec     (1,:)   struct
        % starting time of sim, seconds past J2000
        t0      (1,1)   double
        % starting states of satellites @ t0 in J2000 frame
        x0      (6,:)   double
        % max degree/order of gravity model to use
        ord     (1,1)   {mustBeInteger,mustBePositive} = 1
        % number of satellites
        nsats   (1,1)   {mustBePositive,mustBeInteger} = 1
        % ODE45 options
        opts    (1,1)   struct
        % info stored from latest run (if frame == '', no run yet)
        ts      (1,:)   double
        xs      (6,:,:) double
        frame   (1,:)   char = ''
    end
    
    methods
        function obj = LunarPropagator(t0,x0,ord,nbods,varargin)
            %LUNARPROPAGATOR Construct a LunarPropagator instance.
            %   Inputs:
            %    - t0; character string, 'DD-MMM-YYYY XX:XX:XX'
            %    - x0; starting states -- either array of OE structs, (6,n)
            %          array of starting states (MOON_OP frame), or xopt output
            %          from Conopt(2)
            %    - nbods; what secondary bodies to include (1: Earth,
            %            2:+sun, 3:+jupiter)
            %    - opts; optional name-value arg, ODE45 integration tolerances
            arguments
                t0      (1,:)
                x0      (:,:)   
                ord     (1,1)   {mustBeInteger,mustBePositive}
                nbods   (1,1)   {mustBeInteger,mustBeNonnegative}
            end
            arguments (Repeating)
                varargin
            end

            obj.ord = ord;
            obj.opts = odeset("RelTol", 1e-9, "AbsTol", 1e-11);
            default = 4;

            if nargin > default
                for i=1:2:nargin-default
                    if strcmp(varargin{i}, "opts")
                        obj.opts = varargin{i+1};
                    end
                end
            elseif nargin < default
                error("LunarPropagator:nargin", "Too few arguments.");
            end
            
            cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
            [R,C,S] = cofloader("data/LP165P.cof");
            
            % planetary info
            bods = getplanets("MOON", "EARTH", "SUN", "JUPITER");
            bods(1).R = R * 1e-3;           % convert from m to km
            bods(1).C = C;                  % store in moon struct for orbitaldynamics
            bods(1).S = S;                  % store in moon struct for orbitaldynamics
            bods(1).frame = 'MOON_ME';      % body-fixed frame of coefficients
            obj.moon = bods(1);             % primary body
            obj.sec = bods(2:nbods+1);      % secondary bodies

            if isa(t0, 'double')
                obj.t0 = t0;
            else
                obj.t0 = cspice_str2et(t0);
            end

            % Parse x0
            errmsg = "x0 must be either an array of OE structs, (6,n) " + ...
                    "array of starting states, or xopt output from Conopt(2)";
            if isa(x0, 'struct')        % oes provided
                oes = x0;
                obj.nsats = length(x0);
            elseif isa(x0, 'double')    % not oes
                % xopt output or states
                if all(size(x0) == [1 14]) || all(size(x0) == [1 24])
                    oes = xopt2oes(x0);
                    obj.nsats = 6;
                elseif size(x0,1) == 6
                    obj.x0 = cspice_sxform('MOON_OP', 'J2000', obj.t0) * x0;
                    obj.nsats = size(x0,2);

                    return;     % return early to avoid oes2x0
                else                            % invalid input
                    error("LunarPropagator:invalidInput", errmsg);
                end
            else                                % invalid input
                error("LunarPropagator:invalidInput", errmsg);
            end

            % Convert oes to states
            obj.x0 = zeros(6,obj.nsats);
            for i=1:obj.nsats
                [r,v] = oe2rv(oes(i).a,oes(i).e,oes(i).i,oes(i).RAAN,oes(i).w,oes(i).f,obj.moon.GM);
                obj.x0(:,i) = cspice_sxform('MOON_OP', 'J2000', obj.t0) * [r; v];
            end
        end
        
        function [ts,xs] = run(obj,tf,n,frame)
            %RUN Propagate the input states for tf seconds (n steps
            %between). Data returned in provided frame.
            %   Input:
            %    - tf; final time, seconds past t0
            %    - n; number of time steps
            %    - frame; reference frame to return data in
            arguments
                obj     (1,1)   LunarPropagator
                tf      (1,1)   double {mustBePositive}
                n       (1,1)   {mustBeInteger,mustBePositive}
                frame   (1,:)   char
            end

            ts = linspace(obj.t0, obj.t0+tf, n);
            xs = zeros(6,n,obj.nsats);
            
            for i=1:obj.nsats
                [~,X] = ode45(@(t,x) orbitaldynamics(t,x,obj.moon,obj.ord,obj.sec), ...
                              ts, obj.x0(:,i), obj.opts);
                X = X';
                for j=1:length(ts)
                    X(:,j) = cspice_sxform('J2000', frame, ts(j)) * X(:,j);
                end
                xs(:,:,i) = X;
            end

            obj.ts = ts;
            obj.xs = xs;
            obj.frame = frame;
        end

        function [ts,xs] = runat(obj,ts,frame)
            %RUNAT Propagate the input states over the provided time steps. 
            %Data returned in provided frame.
            %   Input:
            %    - ts; eval time steps, seconds past t0
            %    - frame; reference frame to return data in
            arguments
                obj     (1,1)   LunarPropagator
                ts      (1,:)   double {mustBeNonnegative}
                frame   (1,:)   char
            end

            n = length(ts);
            ts = obj.t0 + ts;
            xs = zeros(6,n,obj.nsats);
            
            for i=1:obj.nsats
                [~,X] = ode45(@(t,x) orbitaldynamics(t,x,obj.moon,obj.ord,obj.sec), ...
                              ts, obj.x0(:,i), obj.opts);
                X = X';
                for j=1:length(ts)
                    X(:,j) = cspice_sxform('J2000', frame, ts(j)) * X(:,j);
                end
                xs(:,:,i) = X;
            end

            obj.ts = ts;
            obj.xs = xs;
            obj.frame = frame;
        end

        function [fx,C] = ephemerisfit(obj,type,dt,N)
            %EPHEMERISFIT Fits given surrogate model type to propagated data.
            %   Input:
            %    - type; "Kepler" or "polynomial", type of model fit -- fit
            %            to error from solving Kepler's problem, or entire
            %            trajectory
            %    - tf; time after t0 to propagate to
            %    - N; number of interpolation points
            arguments
                obj     (1,1)   LunarPropagator
                type    (1,:)   {mustBeText}
                dt      (1,1)   double {mustBePositive}
                N       (1,1)   {mustBePositive,mustBeInteger}
            end

            % basis = @(tau) cell2mat(arrayfun(@(x) chebyshevT(0:N_APPX,x), tau, 'UniformOutput', false));

            % use chebichev nodes for interpolation
            span = chebichev(N);
            t_interp = (span + 1) * dt / 2;

            % Parse different model types
            if strcmpi(type, "Kepler")
                x_base = obj.keplertool(t_interp);
                f_base = @(tau) obj.keplertool(tau - obj.t0);
            elseif strcmpi(type, "polynomial")
                x_base = zeros(6,N+1);
                f_base =  @(tau) zeros(6,length(tau));
            else
                error("ephemerisfit:invalidType", ...
                    "Model type must be either 'Kepler' or 'polynomial'.");
            end

            % propagate state at interpolation points
            [~,x] = obj.runat(t_interp, 'J2000');
            dx = x - x_base;

            % compute coefficients for basis and generate model function
            phi = (span').^(0:N);
            B = pinv(phi);
            C = B * dx(1:3,:)';
            basis = @(tau) [tau'.^(0:N) ...
                (0:N).*(tau'.^([0 0:N - 1])) * 2/dt];
            D = [C zeros(size(C)); zeros(size(C)) C];

            % final model function
            fx = @(tau) (basis(2*(tau - obj.t0)/dt - 1) * D)' ...
                        + f_base(tau);
        end

        function x = keplertool(obj,ts)
            %KEPLERTOOL Wrapper for Kepler_universal to handle multiple
            %time requests.
            %   Input:
            %    - ts; states to propagate to by solving Kepler's problem
            arguments
                obj (1,1)   LunarPropagator
                ts  (1,:)   double {mustBeNonnegative}
            end

            x = zeros(6,length(ts));
            r0 = obj.x0(1:3);
            v0 = obj.x0(4:6);

            for i=1:length(ts)
                [rf,vf] = Kepler_universal(r0, v0, ts(i), obj.moon.GM, 1e-10);
                x(:,i) = [rf; vf];
            end
        end

        function plotlastorbits(obj,frame)
            %PLOTLASTORBITS Generates a plot of the most recently created
            %satellite trajectories in the provided frame.
            %   Input:
            %    - frame; reference frame to plot trajectories in
            arguments
                obj     (1,1)   LunarPropagator
                frame   (1,:)   char
            end

            data = obj.xs;

            % convert data to new frame if required
            if isempty(obj.frame)
                error("plotlastorbits:noData", ...
                    "No data has been generated yet!");
            elseif ~strcmpi(frame, obj.frame)
                for i=1:size(data,2)
                    T = cspice_sxform(obj.frame, frame, obj.ts(i));
                    for j=1:obj.nsats
                        data(:,i,j) = T * data(:,i,j);
                    end
                end
            end

            plotformat("IEEE", 1, "scaling", 2, "coloring", "science");
            plotLunarOrbit(obj.ts, permute(data, [2,1,3]), frame, "Satellite trajectories");
        end

        function [comp,exp] = computedriftrates(obj)
            %COMPUTEDRIFTRATES Computes the drift of right ascension for
            %each orbit over the propagation period.
            %   Output:
            %    - comp; (1,nsats) drift rates computed from frozen orbit eqs
            %    - exp; (1,nsats) drift rates calculated from propagation
            arguments
                obj (1,1)   LunarPropagator
            end
            
            % throw error if there hasn't been a propagation yet
            if isempty(obj.frame)
                error("plotlastorbits:noData", ...
                    "No data has been generated yet!");
            end

            comp = zeros(1,obj.nsats);
            exp = zeros(1,obj.nsats);
            for j=1:obj.nsats
                xo = cspice_sxform(obj.frame, 'MOON_OP', obj.ts(1)) * obj.xs(:,1,j);
                xf = cspice_sxform(obj.frame, 'MOON_OP', obj.ts(end)) * obj.xs(:,end,j);
                [a,e,i,r0,~,~] = rv2oe(xo(1:3), xo(4:6), obj.moon.GM);
                [~,~,~,rf,~,~] = rv2oe(xf(1:3), xf(4:6), obj.moon.GM);

                comp(j) = ascendingnodedrift(a,e,i);
                exp(j) = (rf - r0) / (obj.ts(end) - obj.ts(1));
            end
        end
    end
end

