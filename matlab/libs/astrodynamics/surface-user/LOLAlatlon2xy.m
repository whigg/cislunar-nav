function [X,Y] = LOLAlatlon2xy(LAT,LON)
%LOLALATLON2XY Converts latitude and longitude to X/Y coordinates in 
%stereographic lunar south pole frame used by LOLA.
arguments
    LAT (1,1) double
    LON (1,1) double
end

R = 2*1737.4*tan((pi/2-ABS(LAT)) / 2) * sign(LAT);
X = R*sin(LON);
Y = R*cos(LON);
end



