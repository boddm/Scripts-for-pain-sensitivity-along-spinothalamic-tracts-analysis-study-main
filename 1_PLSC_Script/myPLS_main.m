%% Main script to run the myPLS toolbox
% 
% ------------------------------Steps-------------------------------------
%
% This script includes
%   1. Call the script that contains PLS inputs and their descriptions
%   2. Call the function to run PLS and plot results

clc; clear; close all;

plsc_dir = '1_PLSC_Script';

addpath(fullfile(plsc_dir, 'myPLS_functions'))
addpath(fullfile(plsc_dir, 'RotmanBaycrest'))
addpath(fullfile(plsc_dir, 'misc'))

%% Define all input parameters
% Modify this script to set up your PLS analysis

myPLS_inputs

%% Check the validity of all inputs
% !!! Be sure to run this function to check your settings before running PLS !!! 

[input, pls_opts, save_opts] = myPLS_initialize(input, pls_opts, save_opts);

%% Save and plot input data
% myPLS_plot_inputs(input, pls_opts, save_opts)

%% Run PLS analysis (including permutation test and bootstrap)

res = myPLS_analysis(input, pls_opts);

%% Save and plot results
% If you have run multiple PLS analyses, correct for multiple comparisons of the resulting p-values before executing the following function 

myPLS_plot_results(res, save_opts);

rmpath(fullfile(plsc_dir, 'myPLS_functions'))
rmpath(fullfile(plsc_dir, 'RotmanBaycrest'))
rmpath(fullfile(plsc_dir, 'misc'))
