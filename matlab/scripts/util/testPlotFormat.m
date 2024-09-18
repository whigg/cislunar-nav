% testPlotFormat.m
% Author: Mark Hartigan
% Date  : July 31, 2024
% Description:
%    Test the plotformat function to see coloring taken from Scientific
%    Colour Maps.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format short g;         % display short numbers, no scientific notation

%% random lines
Y = (10:15)'- (0:9);

%% plotformat
plotformat("IEEE", 1, "scaling", 2);

figure();
plot(Y);
linestyleorder("mixedstyles");
grid on;
xlabel("x (units)"); ylabel("y (units)"); title("Test color distribution");
