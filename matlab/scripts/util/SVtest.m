%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% awaugh
n = 7500;
[X,Y,Z] = mySphere(n);

accept = Z < sind(-80);
pts = [X; Y; Z];
pts = pts(:,accept);
alts = linspace(1736, 1736 + 125, 4);
pts = [pts * alts(1) pts * alts(2) pts * alts(3) pts * alts(4)];

figure();
scatter3(pts(1,:), pts(2,:), pts(3,:));
axis equal;