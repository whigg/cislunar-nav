% init.m
% Author: Mark Hartigan
% Date  : September 18, 2024
% Description:
%    Load workspace directories and format printing.

clc, clear, close all;
% add all subdirectories to path
addpath(genpath(pwd));
% display long numbers, no scientific notation
format long g;
% switch to custom plot formatting
plotformat("IEEE", 1, "scaling", 2, "coloring", "science");