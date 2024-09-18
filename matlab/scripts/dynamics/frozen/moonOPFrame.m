% moonOPFrame.m
% Author: Mark Hartigan
% Date  : May 1, 2024
% Description:
%    Verify functionality of lunar Earth Orbit Frame frame kernel,
%    'moon_op_240501.tf'. Frame defined in (Folta and Quinn, 2006).

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% load SPICE kernels and get rotation
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
t0 = cspice_str2et('2024 May 1 00:00:00') - 0;
T = cspice_pxform('J2000','MOON_OP', t0);

%% generate OP frame manually
x_e = cspice_spkezr('EARTH',t0,'J2000','NONE','MOON');
r_e = x_e(1:3);
v_e = x_e(4:6);
z_me = cspice_pxform('MOON_ME','J2000',t0) * [0;0;1];

z_op = cross(r_e,v_e) / norm(cross(r_e,v_e));
x_op = cross(z_me,z_op) / norm(cross(z_me,z_op));
y_op = cross(z_op,x_op) / norm(cross(z_op,x_op));

T_man = [x_op'; y_op'; z_op'];

fprintf("Max error: %.3e\n", norm(T-T_man));

% RESULTS:
% Worst case output was
%     > Max error: 5.323e-16
% which lets us conclude that the spice frame kernel is accurate.