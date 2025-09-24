% 2025.09.23 [cw]

clearvars; close all; clc;

% FILES -------------------------------------------------------------------
save_folder = pwd;
load([save_folder, '\rawdata.mat']); % Loads raw data;

% Add folder with modules & functions (which do the main job):
functions_folder = [save_folder, '\functions'];
addpath(functions_folder);

% SEED FOR RANDOMISATION --------------------------------------------------
rng('default');
random_seed = 20241002;

%% ********************************* MAIN *********************************

%% EXPERIMENT 1
[AGG1a, AGG1b] = mod_exp1_aggregator(EXP1a, EXP1b, COLOURNAMING, MunChips);

%% EXPERIMENT2a
AGG2a = mod_exp2_aggregator(EXP2a, MunChips, AGG1b.dir.HeringRYGB);

%% EXPERIMENT2b
AGG2b = mod_exp2_aggregator(EXP2b, MunChips, AGG1b.dir.HeringRYGB);

%% EXPERIMENT 3
AGG3 = mod_exp3_aggregator(EXP3);

%% ********************************* SAVE *********************************

%% SHIFT ALL RAW DATA TO 1 VARIABLE
% Simply to avoid clutter because the raw data will not be used
% subsequently and is only carried forward for the afficionado.
RAW.exp1a   = EXP1a;
RAW.exp1b   = EXP1b;
RAW.exp2a   = EXP2a;
RAW.exp2b   = EXP2b;
RAW.exp3    = EXP3;

%% SAVE PREPARED VARIABLES
save([save_folder, '\aggdata'],...
    'AGG1a', 'AGG1b', 'AGG2a', 'AGG2b', 'AGG3', 'RAW', 'COLOURNAMING', 'MunChips');