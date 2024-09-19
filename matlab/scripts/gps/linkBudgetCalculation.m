% linkBudgetCalculation.m
% Author: Mark Hartigan
% Date  : September 18, 2024
% Description:
%    Computes the link budget between GPS satellites and a user (with a
%    theoretical isotropic antenna).

%% reset
clc, clear, close all;

%% link budget

% get user position
% get satellite positions
% get misalignment of user from antenna boresight
% get directivity from file
% get gain correction factor from other file
% Antenna gain (dB) = directivity (dB) + GCF (dB)
% get antenna power from ODTBX work