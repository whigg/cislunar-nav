classdef LOLAMap < handle
    %LOLA A wrapper object to provide translational functionality for
    %stereographic maps generated from the LOLA data.
    %   See the links below for more details
    %   https://doi.org/10.3847/PSJ/acf3e1
    %   https://pgda.gsfc.nasa.gov/products/90
    
    properties
        map     (:,:)   double
        n       (1,1)   {mustBePositive,mustBeInteger} = 1
        scale   (1,1)   double {mustBePositive} = 1         % km per pixel
        ns      (1,1)   {mustBeInteger} = -1                % hemisphere sign (-1 south, 1 north)
        r       (1,1)   double = 1737.4                     % km, reference radius of moon
    end
    
    methods
        function obj = LOLAMap(path,hemi)
            %LOLAMap Construct an instance of the LOLAMap class, given the
            %path to a .TIF file.
            %   Input:
            %    - path; path to file
            %    - hemi; "SOUTH" or "NORTH", hemisphere of data
            arguments
                path (1,1) {mustBeText}
                hemi (1,1) {mustBeText}
            end

            obj.map = importdata(path);
            obj.n = size(obj.map, 1);

            temp = strsplit(path, "MPP_");
            temp = strsplit(temp(1), "_");
            obj.scale = str2double(temp(end)) * 1e-3;

            if strcmpi(hemi, "SOUTH")
                obj.ns = -1;
            elseif strcmpi(hemi, "NORTH")
                obj.ns = 1;
            else
                error("LOLAMap:hemisphere", ...
                    "Hemisphere must be SOUTH or NORTH.");
            end
        end

        function z = getfromxy(obj,X,Y)
            %GETFROMXY Returns the map value associated with the X/Y
            %coordinates in the stereographic lunar south pole frame.
            %   Input:
            %    - X; x coordinate in km (+Y in MOON_ME)
            %    - Y; y coordinate in km (-X in MOON_ME)
            arguments
                obj (1,1) LOLAMap
                X (1,1) double
                Y (1,1) double
            end

            [I,J] = obj.xy2ind(X,Y);
            z = obj.map(I,J);
        end

        function [X,Y] = latlon2xy(obj,LAT,LON)
            %LATLON2XY Converts latitude and longitude to X/Y coordinates in 
            %stereographic lunar south pole frame used by LOLA.
            %   Input:
            %    - LAT; lunar latitude in rad (-pi/2 to pi/2)
            %    - LON; lunar longitude in rad (-pi to pi)
            %   Output:
            %    - X; x coordinate in km (+Y in MOON_ME)
            %    - Y; y coordinate in km (-X in MOON_ME)
            arguments
                obj (1,1) LOLAMap
                LAT (1,1) double
                LON (1,1) double
            end
            
            R = 2 * obj.r * tan((pi/2-abs(LAT)) / 2) * sign(LAT);
            X = R*sin(LON);
            Y = R*cos(LON);
        end

        function [LAT,LON] = xy2latlon(obj,X,Y)
            %XY2LATLON Converts X/Y coordinates in stereographic lunar south pole
            %frame used by LOLA to latitude and longitude.
            %   Input:
            %    - X; x coordinate in km (+Y in MOON_ME)
            %    - Y; y coordinate in km (-X in MOON_ME)
            %   Output:
            %    - LAT; lunar latitude in rad (-pi/2 to pi/2)
            %    - LON; lunar longitude in rad (-pi to pi)
            arguments
                obj (1,1) LOLAMap
                X (1,1) double
                Y (1,1) double
            end
            
            R = sqrt(X^2 + Y^2);
            LAT = (pi/2 - 2*atan(0.5 * R/obj.r)) * obj.ns;
            LON = atan2(X,Y);
        end

        function [X,Y] = ind2xy(obj,I,J)
            %XY2IND Converts pixel indices for data maps to X/Y coordinates in 
            %stereographic lunar south pole frame used by LOLA.
            %   Input:
            %    - I; image row
            %    - J; image column
            %   Output:
            %    - X; x coordinate in km (+Y in MOON_ME)
            %    - Y; y coordinate in km (-X in MOON_ME)
            arguments
                obj (1,1)   LOLAMap
                I   (1,1)   {mustBePositive,mustBeInteger}
                J   (1,1)   {mustBePositive,mustBeInteger}
            end
            
            X = (J - obj.n/2 - 0.5) * obj.scale;
            Y = (I - obj.n/2 - 0.5) * obj.scale * obj.ns;
        end
        
        function [I,J] = xy2ind(obj,X,Y)
            %XY2IND Converts X/Y coordinates in stereographic lunar south pole
            %frame used by LOLA to pixel indices for data maps.
            %   Input:
            %    - X; x coordinate in km (+Y in MOON_ME)
            %    - Y; y coordinate in km (-X in MOON_ME)
            %   Output:
            %    - I; image row
            %    - J; image column
            arguments
                obj (1,1)   LOLAMap
                X   (1,1)   double
                Y   (1,1)   double
            end
            
            I = round(Y/obj.scale + obj.n/2 + 0.5);
            J = round(X/obj.scale + obj.n/2 + 0.5);
        end

        function [v] = xy2me(obj,X,Y,elev)
            %XY2ME Converts X/Y coordinates in stereographic lunar south pole
            %frame used by LOLA to the moon ME frame. Also optionally uses
            %the DEM to compute altitude.
            %   Input:
            %    - X; x coordinate in km (+Y in MOON_ME)
            %    - Y; y coordinate in km (-X in MOON_ME)
            %    - elev; optional boolean, use DEM if true (default false)
            arguments
                obj     (1,1)   LOLAMap
                X       (1,1)   double
                Y       (1,1)   double
                elev    (1,1)   logical = false
            end

            % get spherical latitude and longitude
            [LAT,LON] = obj.xy2latlon(X, Y);

            % compute altitude, optionally adding from DEM
            Z = 0;
            if elev, Z = obj.getfromxy(X, Y); end
            R = obj.r + Z*1e-3;

            % convert spherical coordinates to Moon ME frame
            [x_me, y_me, z_me] = sph2cart(LON, LAT, R);
            v = [x_me y_me z_me]';
        end
    end
end

