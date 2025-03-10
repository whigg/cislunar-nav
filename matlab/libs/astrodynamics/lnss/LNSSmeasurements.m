classdef LNSSmeasurements < handle
    %LNSSMEASUREMENTS Handler for determining satellites visible to user,
    %computing measurements, determining measurement partials, etc.
    
    properties
        tmeas   (1,:) double        % measurement times
        fsats   (:,1) cell          % cell array of handles to satellites
        R       (:,:) double        % measurement covariance
        moon    (1,1) struct        % structure containing info about moon (R)
        visible (:,:) double        % visibility of n sats at m time steps
        ymeas   (:,:) double        % measurement array
        user    (9,:) double        % true state of user
        max_a   (1,1) double        % maximum satellite beam width
        frame   (1,:) {mustBeText} = 'J2000'                    % reference frame to get data in
        range   (1,1) logical = false                           % boolean, is range desired?
        rate    (1,1) logical = false                           % boolean, is range-rate desired?
        nsats   (1,1) {mustBePositive, mustBeInteger} = 1       % number of LNSS satellites
        m       (1,1) {mustBePositive, mustBeInteger} = 1       % # of steps in tmeas

        c = 299792.458              % km/s, speed of light
    end
    
    methods
        function obj = LNSSmeasurements(obs,tmeas,fsats,R,moon,user,varargin)
            %LNSSMEASUREMENTS Construct an LNSS measurement object
            %instance.
            %   Inputs:
            %    - obs; what measurements to use ("RANGE", "RATE", or "BOTH")
            %    - tmeas; times at which to provide measurements
            %    - fsats; cell array of satellite state function handles;
            %             {@(t,frame) x(t,frame) [pos (km); vel (km/s)]; ... }
            %    - R; measurement covariance matrix, (nsats,nsats) or (2*nsats,2*nsats) depending on obs
            %    - moon; struct containing info about moon (radius at minimum)
            %    - user; true state of user at measurement times [pos (km); vel (km/s); clk]
            %    - max_a; optional name-value pair, maximum antenna beam-width 
            %             of LNSS (in radians)
            %    - 
            
            max_a = 30 * pi / 180;

            if nargin > 6
                for i=1:2:nargin-6
                    if strcmp(varargin{i}, "max_a")
                        max_a = varargin{i+1};
                    elseif strcmp(varargin{i}, "frame")
                        obj.frame = varargin{i+1};
                    end
                end
            elseif nargin < 6
                error("LNSSmeasurements:nargin", "Too few arguments.")
            end

            % assign measurements desired to booleans
            mustBeOption(obs);
            if strcmpi(obs,"RANGE") || strcmpi(obs,"BOTH")
                obj.range = true;
            end
            if strcmpi(obs,"RATE") || strcmpi(obs,"BOTH")
                obj.rate = true;
            end

            obj.tmeas = tmeas;
            obj.m = length(tmeas);
            obj.fsats = fsats;
            obj.R = R;
            obj.nsats = length(fsats);
            obj.moon = moon;
            obj.user = user;
            obj.max_a = max_a;

            % initialize visibility and measurement arrays
            obj.visible = zeros(obj.nsats, obj.m);
            obj.ymeas = zeros((obj.range+obj.rate) * obj.nsats, obj.m);

            % compute visibility and measurements
            obj.computevisible();
        end

        function y = compute(obj,t,x)
            %COMPUTE Computes measurements for user at time t.
            %   Input:
            %    - t; time of measurement
            %    - x; state of user
            arguments
                obj (1,1) LNSSmeasurements
                t   (1,1) double
                x   (9,1) double
            end

            if ~ismember(t,obj.tmeas)
                error("LNSSmeasurements:timestep", ...
                    "t is not in measurement times provided upon initialization.")
            end

            % initialize measurement vector
            y = zeros((obj.range+obj.rate) * obj.nsats,1);
            % get visibility
            visi = obj.visible(:,find(obj.tmeas == t,1));
            % user properties
            r_u = x(1:3);       % user position vector
            v_u = x(4:6);       % user velocity vector
            r_m = -r_u;         % user->moon position vector

            for i=1:obj.nsats   % consider every satellite
                if visi(i)          % if satellite is visible
                    x_s = obj.fsats{i}(t,obj.frame);  % satellite state
                    r_s = x_s(1:3);                 % satellite position vector
                    v_s = x_s(4:6);                 % satellite velocity vector
                    r_us = r_m + r_s;               % user->sat position vector

                    dr = (v_s - v_u)' * r_us / norm(r_us);  % range-rate

                    if obj.range        % obs == "RANGE"
                        y(i) = norm(r_us) + obj.c*x(7);
                        if obj.rate     % obs == "BOTH"
                            y(obj.nsats+i) = dr + obj.c*x(8);
                        end
                    elseif obj.rate     % obs == "RATE"
                        y(i) = dr + obj.c*x(8);
                    end
                end
            end
        end

        function H = partials(obj,t,x)
            %PARTIALS Computes the partial derivative of the measurement
            %model w.r.t. x at the current state.
            %   Inputs:
            %    - t; time of measurement
            %    - x; state of user
            arguments
                obj (1,1) LNSSmeasurements
                t   (1,1) double
                x   (9,1) double
            end

            if ~ismember(t,obj.tmeas)
                error("LNSSmeasurements:timestep", ...
                    "t is not in measurement times provided upon initialization.")
            end

            % get visibility
            visi = obj.visible(:,find(obj.tmeas == t,1));
            % initialize partials matrix
            H = zeros((obj.range+obj.rate) * obj.nsats, length(x));
            % user properties
            r_u = x(1:3);       % user position vector
            v_u = x(4:6);       % user velocity vector
            r_m = -r_u;         % user->moon position vector

            for i=1:obj.nsats   % consider every satellite
                if visi(i)          % if satellite is visible
                    x_s = obj.fsats{i}(t,obj.frame);  % satellite state
                    r_s = x_s(1:3);                 % satellite position vector
                    v_s = x_s(4:6);                 % satellite velocity vector

                    dr = r_m + r_s;     % relative user->sat position
                    dv = v_s - v_u;     % relative user->sat velocity
                    rho = norm(dr);     % scalar range
                    dvdr = dv'*dr;      % dot product of dv and dr

                    if obj.range        % obs == "RANGE"
                        H(i,:) = [-dr'/rho 0 0 0 obj.c 0 0];
                        if obj.rate     % obs == "BOTH"
                            H(obj.nsats+i,:) = [(dr'*dvdr/rho^3 - dv'/rho) -dr'/rho 0 obj.c 0];
                        end
                    elseif obj.rate     % obs == "RATE"
                        H(i,:) = [(dr'*dvdr/rho^3 - dv'/rho) -dr'/rho 0 obj.c 0];
                    end
                end
            end
        end

        function [xs, gdop] = trilaterate(obj)
            %TRILATERATE Performs pseudorange-based multilateration on the
            %associated measurements, returning the estimated state and
            %GDOP.
            %   Output:
            %    - xs; (4xm) position (km) and clock bias (s)
            %    - gdop; (1xm) geometric dilution of precision at each point

            x_last = [obj.user(1:3,1); obj.user(7)];
            nvis = sum(obj.visible, 1);
            gdop = zeros(1, obj.m);
            xs = inf(4, obj.m);
            
            for i=1:obj.m
                if nvis(i) >= 4
                    dx = ones(4,1);
                    while norm(dx) > 1e-6
                        dp = zeros(obj.nsats, 1);
                        G = ones(obj.nsats, 4);
                        for j=1:obj.nsats
                            x_s = obj.fsats{j}(obj.tmeas(i), 'MOON_ME');
                            dp(j) = obj.ymeas(j,i) - norm(x_s(1:3)-x_last(1:3)) - obj.c*x_last(4);
                            G(j,1:3) = -(x_s(1:3) - x_last(1:3))' / norm(x_s(1:3)-x_last(1:3));
                        end
                        dp = dp(obj.visible(:,i) == 1);
                        G = G(obj.visible(:,i) == 1, :);
                        dx = pinv(G) * dp;
                        x_last = x_last + [dx(1:3); dx(4)/obj.c];
                    end
                    
                    H = inv(G'*G);
                    gdop(i) = sqrt(trace(H));
                    xs(:,i) = x_last;
                else
                    gdop(i) = Inf;
                end
            end
        end
    end

    methods (Access = private)
        function computevisible(obj)
            %COMPUTEVISIBLE Called at instantiation; creates array of
            %satellite visibilities at each time step and generates
            %pseudorange and -range-rate measurements while it's at it.

            for i=1:obj.m       % iterate over time steps
                ti = obj.tmeas(i);          % current time
                x_u = obj.user(:,i);        % user state
                r_u = x_u(1:3);             % user position vector
                v_u = x_u(4:6);             % user velocity vector
                r_m = -r_u;                 % user->moon position vector

                % compute moon-center / moon-tangent angle, 
                ang_m = asin(obj.moon.R / norm(r_m));

                for j=1:obj.nsats       % iterate over each satellite
                    x_s = obj.fsats{j}(ti,obj.frame); % satellite state
                    r_s = x_s(1:3);                 % satellite position vector
                    v_s = x_s(4:6);                 % satellite velocity vector
                    r_us = r_m + r_s;               % user->sat position vector

                    % compute moon-center / LNSS sat angle and LNSS-to-moon / LNSS-to-user angle
                    ang_s = acos(dot(r_us,r_m) / (norm(r_us)*norm(r_m)));
                    ang_a = acos(dot(-r_s,-r_us) / (norm(-r_s)*norm(-r_us)));

                    dr = (v_s - v_u)' * r_us / norm(r_us);  % range-rate

                    % satellite is NOT (further than moon AND within its view angle)
                    % AND s/c is within the LNSS satellite beam width
                    if ~(norm(r_s) > norm(r_m) && ang_s < ang_m) && ang_a < obj.max_a
                        obj.visible(j,i) = 1;       % spacecraft is visible

                        % go thru skill tree to assign measurements
                        if obj.range
                            obj.ymeas(j,i) = norm(r_us) + obj.c*x_u(7) + mvnrnd(0,obj.R(j,j));
                        
                            if obj.rate
                                obj.ymeas(obj.nsats+j,i) = dr + obj.c*x_u(8) + ...
                                    mvnrnd(0,obj.R(obj.nsats+j,obj.nsats+j));
                            end
                        elseif obj.rate
                            obj.ymeas(j,i) = dr + obj.c*x_u(8) + mvnrnd(0,obj.R(j,j));
                        end
                    end
                end
            end
        end
    end
end

% custom validation
function mustBeOption(txt)
%MUSTBEOPTION Tests that string option is allowable

if ~strcmpi(txt,"RANGE") && ~strcmpi(txt,"RATE") && ~strcmpi(txt,"BOTH")
    error("LNSSmeasurements:obsType", ...
                    "Observation type obs must be 'RANGE', 'RATE', or 'BOTH'.")
end
end