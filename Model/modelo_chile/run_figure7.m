% run_figure7.m
% Simple script to generate Figure 7 from existing model results

clear all; close all; clc;

fprintf('\n=== LOADING MODEL RESULTS ===\n');

% Load Exercise 2 results (manufacturing TFP shock)
if ~exist('model_output_IOSOE_ex2.mat', 'file')
    error('Model output not found. Run main_SOE_exercise2.m first.');
end

load('model_output_IOSOE_ex2.mat');
fprintf('✓ Loaded model_output_IOSOE_ex2.mat\n');

% Number of sectors
nsec = 12;

% Now run the plotting script
run('plot_figure7_manufacturing_shock.m');
