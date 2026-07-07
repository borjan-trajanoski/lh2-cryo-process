% runStaticBoiloff.m
% Main script for LH2 static boil-off simulation.
%
% Default case: 3,300-gallon Dewar, 10 days (sanity check).
% Compare against Petitpas (2018) Figs 9-11 once validated.
%
% Run this file. All other files must be in the same folder.

clc; close all; tic;

%% 1. Load parameters
fprintf('Loading parameters...\n');
inputs_StaticBoiloff;

%% 2. Run simulation
fprintf('Starting ODE integration...\n');
results = LH2StaticSimulate('Static boil-off (validation)');

%% 3. Post-process
Data_extraction_static;

%% 4. Plot
plotStaticBoiloff(results);

%% 5. Save
save('results_static_boiloff.mat', 'results', 'LH2Model');
fprintf('Results saved to results_static_boiloff.mat\n');

elapsed = toc;
fprintf('Total elapsed time: %.1f seconds (%.1f minutes)\n', elapsed, elapsed/60);
