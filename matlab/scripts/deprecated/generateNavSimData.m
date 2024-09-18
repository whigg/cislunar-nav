% generateNavSimData.m
% Author: Mark Hartigan
% Date  : March 6, 2024
% Description:
%    Generates datafiles for the Python navigation simulations, consisting
%    of time history of satellite positions.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
t0 = convertTo(datetime('6-Feb-2027 00:00:00'), 'juliandate');
t0 = (t0 - 2451545) * 86400;    % JD to seconds past MOON_ME

imdata = {
    'data/NSNS-DRM1.0-Khon1.bsp' 'data/NSNS-DRM1.0-Khon1.tpc.pc' ...
    'data/NSNS-DRM1.0-Khon2.bsp' 'data/NSNS-DRM1.0-Khon2.tpc.pc' ...
    'data/NSNS-DRM1.0-Khon3.bsp' 'data/NSNS-DRM1.0-Khon3.tpc.pc' ...
    'data/NSNS-DRM1.0-Khon4.bsp' 'data/NSNS-DRM1.0-Khon4.tpc.pc' ...
    'data/NSNS-DRM1.4-Khon5asKhon2.bsp' 'data/NSNS-DRM1.4-Khon5asKhon2.tpc.pc' ...
    'data/NSNS-DRM1.4-Khon6asKhon4.bsp' 'data/NSNS-DRM1.4-Khon6asKhon4.tpc.pc' ...
    'data/NSNS-DRM1.1-Khon7asKhon1.bsp' 'data/NSNS-DRM1.1-Khon7asKhon1.tpc.pc'};

gendata = {'data/naif0012.tls', 'data/de430.bsp', 'data/gm_de431.tpc', ...
    'data/moon_080317.tf', 'data/moon_pa_de421_1900-2050.bpc'};
cspice_furnsh([imdata, gendata]);

moon.R  = 1736;         % km, polar radius of moon
moon.GM = cspice_bodvrd('MOON','GM',1);
moon.x  = @(t) zeros(3,1);

getpos = @(x) x(1:3,:);
Khon1 = @(t) getpos(cspice_spkezr('Khon-1', t, 'MOON_ME', 'NONE', 'MOON'));
Khon2 = @(t) getpos(cspice_spkezr('Khon-2', t, 'MOON_ME', 'NONE', 'MOON'));
Khon3 = @(t) getpos(cspice_spkezr('Khon-3', t, 'MOON_ME', 'NONE', 'MOON'));
Khon4 = @(t) getpos(cspice_spkezr('Khon-4', t, 'MOON_ME', 'NONE', 'MOON'));
Khon5 = @(t) getpos(cspice_spkezr('Khon-5', t, 'MOON_ME', 'NONE', 'MOON'));
Khon6 = @(t) getpos(cspice_spkezr('Khon-6', t, 'MOON_ME', 'NONE', 'MOON'));
Khon7 = @(t) getpos(cspice_spkezr('Khon-7', t, 'MOON_ME', 'NONE', 'MOON'));

%% build trajectories
days = 1;
step = 60; % seconds
t = t0:step:t0 + 86400 * days;
xk1 = Khon1(t); xk2 = Khon2(t); xk3 = Khon3(t);
xk4 = Khon4(t); xk5 = Khon5(t); xk6 = Khon6(t); xk7 = Khon7(t);
xk = [xk2; xk3; xk4; xk5; xk6];

%% format data and write to file
data = xk';
writematrix(data, 'D:/Documents/Georgia Tech/_PNT/cislunar-nav/python/data/IM/Khon2-6_MOON_ME.txt', 'Delimiter', 'tab');

