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
                t0      (1,:)   char
                x0      (:,:)   
                ord     (1,1)   {mustBeInteger,mustBePositive}
                nbods   (1,1)   {mustBeInteger,mustBePositive}
            end
            arguments (Repeating)
                varargin
            end

            obj.ord = ord;
            obj.opts = odeset("RelTol", 1e-9, "AbsTol", 1e-11);
            default = 3;

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

            obj.t0 = cspice_str2et(t0);

            % Parse x0
            errmsg = "x0 must be either an array of OE structs, (6,n) " + ...
                    "array of starting states, or xopt output from Conopt(2)";
            if all(class(x0) == 'struct')       % oes provided
                oes = x0;
                obj.nsats = length(x0);
            elseif all(class(x0) == 'double')   % not oes
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

        function plotlastorbits()
            %PLOTLASTORBITS Generates a plot of the most recently created
            %satellite trajectories.

            
        end
    end
end

