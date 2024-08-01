% LOLAelevTester.m
% Author: Mark Hartigan
% Date  : July 17, 2024
% Description:
%    Plot the elevation map provided by LOLA's laser altimeter and GSFC.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% import data
data = importdata("data/LDEM_80S_80MPP_ADJ.TIF") ./ 1000;
dmin = min(data,[],"all");
dmax = max(data,[],"all");
range = 304;            % km, stereographic range
mpp = 80;               % m, meters per pixel
kmpp = mpp * 1e-3;
range = range - kmpp/2;
px2km = -range:kmpp:range;

%% plot
figure();
imagesc([-range range], [range -range], data);
set(gca, 'YDir', 'normal');
xlabel("y_{ME} (km)");ylabel("x_{ME} (km)");
c = colorbar();
c.Label.String = "Elevation from mean radius (km)";
load('roma.mat');
colormap(flipud(roma));