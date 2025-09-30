% ANALYSES & FIGURES 
% This script produces the analyses and figures in the manuscript. When you
% run the whole script, it will stop after loading the data (this section)
% to avoid running through all analyses in one go. To run single analyses,
% load the data, and then go to the respective section and run the section
% only (by pressing ctrl + return). If you do want to run the whole script
% in one go, simply delete "return" at the end of the section and make
% yourself a tea because it will take a while :). 
% 2025.09.23 [cw]

clearvars; close all; clc;

% FILES -------------------------------------------------------------------
src_path = [pwd, '\'];
functions_folder = fullfile(src_path, 'functions');
addpath(functions_folder);
load(fullfile(src_path, 'aggdata.mat'));

% SETTINGS FOR OUTPUT ----------------------------------------------------- 
save_folder = fullfile(src_path, 'output_figures_and_tables');
tab_file = fullfile(save_folder, 'Tables.xlsx');
ani_file = fullfile(save_folder, 'animation');
fsz = 9; fsz2 = 10;
punits = 'centimeters';
%saveformat      = 'all'; % To get pdf, but some take long.
saveformat      = 'png';
saveformat_sup  = 'png';
return

%% **************************** MAIN: Article *****************************

%% GRAPHICAL ABSTRACT
close all; clc;
[~, fh] = mod_exp2_analysor(AGG2a, 'Graphical Abstract');

% REFORMAT ----------------------------------------------------------------
fh.Color = 'k';
fh.InvertHardcopy = 'off';
wd = .5; ht = 1; y = 0;
ax = get(fh, 'Children');
set(ax, 'Position', [.25 0 wd, ht]);

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'graphicalabstract', [16 8], saveformat, punits, 300);

%% Figure 1 | Afterimage Stimuli & Models. 
close all; clc;
[SIM, fh] = mod_afterimage_modeller(MunChips, 'Luv');

% REFORMAT ----------------------------------------------------------------
fh.Color = 'w';
fh.InvertHardcopy = 'off';
axs = get(fh, 'Children');
axs = flipud(axs);
xlabel(axs(1),'Green-Red [u*]');
ylabel(axs(1),'Blue-Yelow [v*]');
axs(2).YTickLabel = [];
xlabel(axs(2:5),'');
ylabel(axs(2:5),'');

% POSITION & SIZE
wd = .18; ht = .72; y = 0.15;
x = [.03 .225 .42 .615 .81];
set(axs(1), 'Position', [x(1) y wd, ht]);
set(axs(2), 'Position', [x(2) y wd, ht]);
set(axs(3), 'Position', [x(3) y wd, ht]);
set(axs(4), 'Position', [x(4) y wd, ht]);
set(axs(5), 'Position', [x(5) y wd, ht]);

panel_identifier(1,5, [.12 .05], .05, fsz2, 'rows', axs, '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'fig1_models', [30 7], saveformat, punits, 300);

%% Table 1 | Differences from Opponency in Experiments 1a.

% PRODUCE -----------------------------------------------------------------
RES1 = mod_exp1_analysor(AGG1a, AGG1b, 'None');
Table1 = RES1.ttab_exp1a;

% SAVE --------------------------------------------------------------------
writetable(Table1, tab_file, 'Sheet', 'Table1', 'WriteRowNames',true);

%% Figure 2 | Experiment 1-2.
close all; clc;

% FIGURE 2b-d: Experiment 1.a-b -------------------------------------------
RES = mod_exp1_analysor(AGG1a, AGG1b, 'MAIN');

% FIGURE 2e-g: Exp2.a Aggregated ------------------------------------------
RES = mod_exp2_analysor(AGG2a, 'MAIN');

% FIGURE 2.i-j: Prediction Errors across Experiments ----------------------
RES = mod_prederror_analysor(AGG1a, AGG1b, AGG2a, AGG2b, 'MAIN');

% REFORMAT ----------------------------------------------------------------
[fh, axs] = subplot_recombiner([NaN, 1,1,1, NaN, 2,2,2], 2, 4, [1, 1:3, 1, 1:3]);
set(fh,'InvertHardCopy', 'off', 'Color', 'w');
xlabel(axs(2:3), '');
title(axs(4:6), '');
wd = .2; ht = .36; y = [.58 .08]; x = [.05 .30 .55 .80];
set(axs(1), 'Position', [x(2) y(1) wd, ht]);
set(axs(2), 'Position', [x(3) y(1) wd, ht]);
set(axs(3), 'Position', [x(4) y(1) wd, ht]);
set(axs(4), 'Position', [x(2) y(2) wd, ht]);
set(axs(5), 'Position', [x(3) y(2) wd, ht]);
set(axs(6), 'Position', [x(4) y(2) wd, ht]);

panel_identifier(1,3, [.18 .15], .05, fsz2, 'rows', axs(1:3), '.)', 2);
panel_identifier(1,3, [.18 .15], .05, fsz2, 'rows', axs(4:6), '.)', 6);

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'fig2b-h_afteri', [27 13], saveformat, punits, 300);

% REFORMAT ----------------------------------------------------------------
fh = figure(3);
axs = get(fh, 'Children');
axs = axs([3 1]);
wd1 = .52; wd2 = .24; ht = .75; y = .12; x = [.05 .75];
set(axs(1), 'Position', [x(1) y wd1, ht]);
set(axs(2), 'Position', [x(2) y wd2, ht]);
ylim([0 3.3]);
panel_identifier(1,2, [.07 .16], .05, fsz2, 'rows', axs, '.)', 9);

fig_saver(save_folder, 'fig2i-l_afteri', [27 5], saveformat, punits, 300);

%% METHOD
METH2a = mod_exp2_analysor(AGG2a, 'Method');
METH2b = mod_exp2_analysor(AGG2b, 'Method');

%% TABLE 2: Correlations in Experiment 2

% PRODUCE -----------------------------------------------------------------
RES1 = mod_exp1_analysor(AGG1a, AGG1b, 'Correlation Tables');
Table2.a = RES1.agg_ctab;

RES2a    = mod_exp2_analysor(AGG2a, 'Correlation Tables');
Table2.b = [RES2a.agg_ctab(1,:); RES2a.IndiM.ctab];
Table2.c = [RES2a.agg_ctab(2,:); RES2a.IndiOppoDev.ctab];
Table2.d = [RES2a.agg_ctab(3,:); RES2a.IndiHueFrq.ctab];

RES2b    = mod_exp2_analysor(AGG2b, 'Correlation Tables');
Table2.e = [RES2b.agg_ctab(1,:); RES2b.IndiM.ctab];
Table2.f = [RES2b.agg_ctab(2,:); RES2b.IndiOppoDev.ctab];
Table2.g = [RES2b.agg_ctab(3,:); RES2b.IndiHueFrq.ctab];

% SAVE --------------------------------------------------------------------
writetable(Table2.a, tab_file, 'Sheet', 'Table2a', 'WriteRowNames',true);
writetable(Table2.b, tab_file, 'Sheet', 'Table2b', 'WriteRowNames',true);
writetable(Table2.c, tab_file, 'Sheet', 'Table2c', 'WriteRowNames',true);
writetable(Table2.d, tab_file, 'Sheet', 'Table2d', 'WriteRowNames',true);
writetable(Table2.e, tab_file, 'Sheet', 'Table2e', 'WriteRowNames',true);
writetable(Table2.f, tab_file, 'Sheet', 'Table2f', 'WriteRowNames',true);
writetable(Table2.g, tab_file, 'Sheet', 'Table2g', 'WriteRowNames',true);

%% Table 3 | Model Comparisons.
close all; clc;

% PRODUCE -----------------------------------------------------------------
RES = mod_prederror_analysor(AGG1a, AGG1b, AGG2a, AGG2b, 'Tables');
Table3.a = RES.ttabs.hue1a_across_pp;
Table3.b = RES.ttabs.hue1b_across_pp;
Table3.c = RES.stabs.hue2a_across_sti;
Table3.d = RES.stabs.hue2b_across_sti;
Table3.e = RES.stabs.chro2a_across_sti;
Table3.f = RES.stabs.chro2b_across_sti;

% SAVE --------------------------------------------------------------------
writetable(Table3.a, tab_file, 'Sheet', 'Table3a', 'WriteRowNames',true);
writetable(Table3.b, tab_file, 'Sheet', 'Table3b', 'WriteRowNames',true);
writetable(Table3.c, tab_file, 'Sheet', 'Table3c', 'WriteRowNames',true);
writetable(Table3.d, tab_file, 'Sheet', 'Table3d', 'WriteRowNames',true);
writetable(Table3.e, tab_file, 'Sheet', 'Table3e', 'WriteRowNames',true);
writetable(Table3.f, tab_file, 'Sheet', 'Table3f', 'WriteRowNames',true);

%% Figure 3 | Change of Afterimage Hue with Chroma (Experiment 3).
close all; clc;

% PRODUCE -----------------------------------------------------------------
RES3 = mod_exp3_analysor(AGG3, 'MAIN');

% REFORMAT ----------------------------------------------------------------
[fh, axs] = subplot_recombiner(ones(6,1), 3, 2, 1:6);

xlabel(axs(2), 'Green-Red [u*]');
ylabel(axs(2), 'Blue-Yellow [v*]');
xlabel(axs(3), 'Green-Red [u*]');
ylabel(axs(3), 'Blue-Yellow [v*]');
ylabel(axs(6), []);

axis(axs(3), [-200 200 -300 100]);
axis(axs(4), [-1 1 -1 1]*100);
axis(axs(5), [-1 1 -1 1]*50);
axis(axs(6), [-1 1 -1 1]*25);

set(axs(5), 'XTick', -50:25:50, 'YTick', -50:25:50);
set(axs(6), 'XTick', -25:25:25, 'YTick', -25:25:25);

wd = .35; ht = .27; y = [.7 .37 .04]; x = [.11 .62];
set(axs(1), 'Position', [.14 .69 wd, ht]);
set(axs(2), 'Position', [x(2) y(1) wd, ht]);
set(axs(3), 'Position', [x(1) y(2) wd, ht]);
set(axs(4), 'Position', [x(2) y(2) wd, ht]);
set(axs(5), 'Position', [x(1) y(3) wd, ht]);
set(axs(6), 'Position', [x(2) y(3) wd, ht]);

panel_identifier(3,2, .18, .08, fsz2, 'rows', axs(2:6), '.)', 2);

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'fig3 Saturation', [10 16], saveformat, punits, 300);

%% Table 4 | Correlations of Experiments 3.

% PRODUCE -----------------------------------------------------------------
RES3 = mod_exp3_analysor(AGG3, 'none');

Table4.a = RES3.ctab_oppodev;
Table4.b = RES3.ctab_Mdev;

% SAVE --------------------------------------------------------------------
writetable(Table4.a, tab_file, 'Sheet', 'Table4a', 'WriteRowNames',true);
writetable(Table4.b, tab_file, 'Sheet', 'Table4b', 'WriteRowNames',true);

%% ************************ SUPPLEMENTARY MATERIAL ************************

%% ANIMATIONS & FIGURE S1
close all; clc;

% 2 EXAMPLES THAT DIFFER FROM OPPONENT ------------------------------------

% ANIMATION 1 (60/80)
iAzi = 60;
iChroma = 80;
mod_animation_maker([ani_file, '1'], iAzi, iChroma);

% ANIMATION 2 (240/80)
iAzi = 240;
iChroma = 80;
mod_animation_maker([ani_file, '2'], iAzi, iChroma);

% 2 EXAMPLES THAT COOCCUR WITH OPPONENT -----------------------------------
% ANIMATION 3 (340/80)
iAzi = 340;
iChroma = 80;
mod_animation_maker([ani_file, '3'], iAzi, iChroma);

% ANIMATION 4 (100/80)
iAzi = 100;
iChroma = 80;
mod_animation_maker([ani_file, '4'], iAzi, iChroma);

% REFORMAT (Figure S1) ----------------------------------------------------

% MERGE IN 1 FIGURE:
[fh, axs] = subplot_recombiner(1:4, 2, 2, [1 1 1 1]);
set(fh,'InvertHardCopy', 'off', 'Color', 'w');
axis(axs, 'off', 'square', 'tight');

% POSITION & SIZE:
wd = .46; ht = .46; 
x = [.03 .53]; y = [.51 .01];
set(axs(1), 'Position', [x(1) y(1) wd, ht]);
set(axs(2), 'Position', [x(2) y(1) wd, ht]);
set(axs(3), 'Position', [x(1) y(2) wd, ht]);
set(axs(4), 'Position', [x(2) y(2) wd, ht]);

% PANELS:
panel_identifier(2, 2, 0, 0, fsz2, 'rows', axs, '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS1_animations', [15 15], saveformat_sup, punits, 300);

%% FIGURE S2 | Afterimage Matches & Control Task in Experiment 1a
close all; clc;
RES = mod_exp1_analysor(AGG1a, AGG1b, 'Supplementary');

% FORMAT ------------------------------------------------------------------
axs = get(gcf, 'Children'); axs = flipud(axs);
% POSITION & SIZE:
wd = .4; ht = .9; 
x = [.08 .58]; y = .05;
set(axs(1), 'Position', [x(1) y wd, ht]);
set(axs(2), 'Position', [x(2) y wd, ht]);

% PANELS:
panel_identifier(1, 2, [.15 .1], 0, fsz2, 'rows', axs, '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS2 Afterimage_Matches_&_Control_Task', [20 10], saveformat_sup, punits, 300);

%% Figure S3 | Hue Histograms for Experiment 1a.
close all; clc;
RES = mod_exp1_analysor(AGG1a, AGG1b, 'Histo');

% FORMATTING --------------------------------------------------------------
axs = get(gcf,'Children');
axs = flipud(axs);

% POSITION & SIZE
wd = .19; ht = .38; 
x = [.05 .3 .55 .8]; y = [.57 .08]; 
set(axs(1),  'Position', [x(1) y(1) wd, ht]);
set(axs(2),  'Position', [x(2) y(1) wd, ht]);
set(axs(3),  'Position', [x(3) y(1) wd, ht]);
set(axs(4),  'Position', [x(4) y(1) wd, ht]);
set(axs(5),  'Position', [x(1) y(2) wd, ht]);
set(axs(6),  'Position', [x(2) y(2) wd, ht]);
set(axs(7),  'Position', [x(3) y(2) wd, ht]);
set(axs(8),  'Position', [x(4) y(2) wd, ht]);

set(axs(2:4), 'YTickLabel', []);
set(axs(6:8), 'YTickLabel', []);

% PANEL IDS
panel_identifier(2, 4, [.15 .04], .04, fsz2, 'rows', axs,  '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS3 Hue_Histograms', [20 11], 'png', punits, 300);

%% Figure S4 | Example Distributions of Afterimage Hue Measurements in Experiment 2a.
close all; clc;
RES = mod_exp2_analysor(AGG2a, 'Distribution Hue');

% FORMATTING --------------------------------------------------------------
axs = get(gcf,'Children');
axs = flipud(axs);
wd = .13; ht = .38; y = [.57 .08]; x = [.05 .21 .37 .53 .69 .85];

set(axs(1),  'Position', [x(1) y(1) wd, ht]);
set(axs(2),  'Position', [x(2) y(1) wd, ht]);
set(axs(3),  'Position', [x(3) y(1) wd, ht]);
set(axs(4),  'Position', [x(4) y(1) wd, ht]);
set(axs(5),  'Position', [x(5) y(1) wd, ht]);
set(axs(6),  'Position', [x(6) y(1) wd, ht]);
set(axs(7),  'Position', [x(1) y(2) wd, ht]);
set(axs(8),  'Position', [x(2) y(2) wd, ht]);
set(axs(9),  'Position', [x(3) y(2) wd, ht]);
set(axs(10), 'Position', [x(4) y(2) wd, ht]);
set(axs(11), 'Position', [x(5) y(2) wd, ht]);
set(axs(12), 'Position', [x(6) y(2) wd, ht]);

xlabel(axs(1:6), []);
ylabel(axs(2:6), []);
ylabel(axs(8:12), []);
set(axs(2:6),  'YTickLabel', []);
set(axs(8:12), 'YTickLabel', []);

panel_identifier(2,6, [.22 .04], .04, fsz2, 'rows', axs,  '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS4 distri_hue', [30 11], 'png', punits, 300);

%% Figure S5 | Example Distribution of Afterimage Chroma Measurements in Experiment 2a.
close all; clc;
RES = mod_exp2_analysor(AGG2a, 'Distribution Chroma');

% FORMATTING --------------------------------------------------------------
axs = get(gcf,'Children');
axs = flipud(axs);
wd = .13; ht = .38; y = [.57 .08]; x = [.05 .21 .37 .53 .69 .85];

set(axs(1),  'Position', [x(1) y(1) wd, ht]);
set(axs(2),  'Position', [x(2) y(1) wd, ht]);
set(axs(3),  'Position', [x(3) y(1) wd, ht]);
set(axs(4),  'Position', [x(4) y(1) wd, ht]);
set(axs(5),  'Position', [x(5) y(1) wd, ht]);
set(axs(6),  'Position', [x(6) y(1) wd, ht]);
set(axs(7),  'Position', [x(1) y(2) wd, ht]);
set(axs(8),  'Position', [x(2) y(2) wd, ht]);
set(axs(9),  'Position', [x(3) y(2) wd, ht]);
set(axs(10), 'Position', [x(4) y(2) wd, ht]);
set(axs(11), 'Position', [x(5) y(2) wd, ht]);
set(axs(12), 'Position', [x(6) y(2) wd, ht]);

xlabel(axs(1:6), []);
ylabel(axs(2:6), []);
ylabel(axs(8:12), []);
set(axs(2:6),  'YTickLabel', []);
set(axs(8:12), 'YTickLabel', []);

panel_identifier(2,6, [.22 .04], .04, fsz2, 'rows', axs,  '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS5 distri_chroma', [30 11], 'png', punits, 300);

%% Figure S6 | Normal Probability Plot for Afterimage Hue Measurements in Experiment 2a.
close all; clc;
RES = mod_exp2_analysor(AGG2a, 'QQ hue');

% FORMATTING --------------------------------------------------------------
axs = get(gcf,'Children');
axs = flipud(axs);
wd = .13; ht = .38; y = [.57 .08]; x = [.05 .21 .37 .53 .69 .85];

set(axs(1),  'Position', [x(1) y(1) wd, ht]);
set(axs(2),  'Position', [x(2) y(1) wd, ht]);
set(axs(3),  'Position', [x(3) y(1) wd, ht]);
set(axs(4),  'Position', [x(4) y(1) wd, ht]);
set(axs(5),  'Position', [x(5) y(1) wd, ht]);
set(axs(6),  'Position', [x(6) y(1) wd, ht]);
set(axs(7),  'Position', [x(1) y(2) wd, ht]);
set(axs(8),  'Position', [x(2) y(2) wd, ht]);
set(axs(9),  'Position', [x(3) y(2) wd, ht]);
set(axs(10), 'Position', [x(4) y(2) wd, ht]);
set(axs(11), 'Position', [x(5) y(2) wd, ht]);
set(axs(12), 'Position', [x(6) y(2) wd, ht]);

xlabel(axs(1:6), []);
ylabel(axs(2:6), []);
ylabel(axs(8:12), []);
set(axs(2:6),  'YTickLabel', []);
set(axs(8:12), 'YTickLabel', []);

panel_identifier(2,6, [.22 .04], .04, fsz2, 'rows', axs,  '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS6 QQ_hue', [30 11], 'png', punits, 300);

%% Figure S7 | Normal Probability Plot for Afterimage Chroma Measurements in Experiment 2a.
close all; clc;
[RES, fh] = mod_exp2_analysor(AGG2a, 'QQ chroma');

% FORMATTING --------------------------------------------------------------
axs = get(gcf,'Children');
axs = flipud(axs);
wd = .13; ht = .38; y = [.57 .08]; x = [.05 .21 .37 .53 .69 .85];

set(axs(1),  'Position', [x(1) y(1) wd, ht]);
set(axs(2),  'Position', [x(2) y(1) wd, ht]);
set(axs(3),  'Position', [x(3) y(1) wd, ht]);
set(axs(4),  'Position', [x(4) y(1) wd, ht]);
set(axs(5),  'Position', [x(5) y(1) wd, ht]);
set(axs(6),  'Position', [x(6) y(1) wd, ht]);
set(axs(7),  'Position', [x(1) y(2) wd, ht]);
set(axs(8),  'Position', [x(2) y(2) wd, ht]);
set(axs(9),  'Position', [x(3) y(2) wd, ht]);
set(axs(10), 'Position', [x(4) y(2) wd, ht]);
set(axs(11), 'Position', [x(5) y(2) wd, ht]);
set(axs(12), 'Position', [x(6) y(2) wd, ht]);

xlabel(axs(1:6), []);
ylabel(axs(2:6), []);
ylabel(axs(8:12), []);
set(axs(2:6),  'YTickLabel', []);
set(axs(8:12), 'YTickLabel', []);

panel_identifier(2,6, [.22 .04], .04, fsz2, 'rows', axs,  '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS7 QQ_chroma', [30 11], 'png', punits, 300);

%% Figure S8 | Averages for Individual Observers in Experiment 2a.
close all; clc;
RES = mod_exp2_analysor(AGG2a, 'Individual Adjustments');

% FORMATTING --------------------------------------------------------------
axs = get(gcf,'Children');
axs = flipud(axs);
wd = .16; ht = .36; y = [.59 .1]; x = [.05 .24 .43 .62 .81];
set(axs(1), 'Position', [x(1) y(1) wd, ht]);
ylabel(axs(1), 'Blue-Yellow [v*]');
set(axs(2), 'Position', [x(2) y(1) wd, ht]);
set(axs(3), 'Position', [x(3) y(1) wd, ht]);
set(axs(4), 'Position', [x(4) y(1) wd, ht]);
set(axs(5), 'Position', [x(5) y(1) wd, ht]);
set(axs(6), 'Position', [x(1) y(2) wd, ht]);
ylabel(axs(6), 'Blue-Yellow [v*]');
set(axs(7), 'Position', [x(2) y(2) wd, ht]);
set(axs(8), 'Position', [x(3) y(2) wd, ht]);
xlabel(axs(8), 'Green-Red [u*]');
set(axs(9), 'Position', [x(4) y(2) wd, ht]);
set(axs(10), 'Position', [x(5) y(2) wd, ht]);
set(axs(1:5), 'XLabel', [], 'XTickLabel', []);
set(axs([6 7 9 10]), 'XLabel', []);
set(axs([2:5, 7:10]),'YLabel', [],'YTickLabel', []);

panel_identifier(2,5, [.2 .07], .04, fsz2, 'rows', axs, '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS8 Indi_Averages', [25 11], 'png', punits, 300);

%% Figure S9 | Deviations from Cone-Opponency for Individual Observers in Experiment 2a.
close all; clc;
RES = mod_exp2_analysor(AGG2a, 'Individual Deviations From Opponency');

% FORMATTING --------------------------------------------------------------
axs = get(gcf,'Children');
axs = flipud(axs);
wd = .16; ht = .36; y = [.59 .1]; x = [.06 .25 .44 .63 .82];
set(axs(1), 'Position', [x(1) y(1) wd, ht]);
set(axs(2), 'Position', [x(2) y(1) wd, ht]);
set(axs(3), 'Position', [x(3) y(1) wd, ht]);
set(axs(4), 'Position', [x(4) y(1) wd, ht]);
set(axs(5), 'Position', [x(5) y(1) wd, ht]);
set(axs(6), 'Position', [x(1) y(2) wd, ht]);
set(axs(7), 'Position', [x(2) y(2) wd, ht]);
set(axs(8), 'Position', [x(3) y(2) wd, ht]);
set(axs(9), 'Position', [x(4) y(2) wd, ht]);
set(axs(10), 'Position', [x(5) y(2) wd, ht]);

ylim(axs(:),[-45 45]);
set(axs(:),'YTick', -45:15:45);
xlabel(axs(:), '');
ylabel(axs(:), '');
xlabel(axs(8), 'Inducer Huer [deg]');
ylabel(axs(1), 'Difference from Opponent [deg]');
ylabel(axs(6), 'Difference from Opponent [deg]');
set(axs(1:5),'XTickLabel', []);
set(axs([6 7 9 10]), 'XLabel', []);
set(axs([2:5, 7:10]),'YLabel', [],'YTickLabel', []);

panel_identifier(2,5, [.27 .07], .04, fsz2, 'rows', axs, '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS9 Indi_OppoDev', [25 11], 'png', punits, 300);

%% Figure S10 | Polar Hue Histogram for Individual Observers in Experiment 2a.
close all; clc;
RES = mod_exp2_analysor(AGG2a, 'Individual Hue Histograms');

% FORMATTING --------------------------------------------------------------
axs = get(gcf,'Children');
axs = flipud(axs);
wd = .16; ht = .36; y = [.59 .1]; x = [.06 .25 .44 .63 .82];
set(axs(1), 'Position', [x(1) y(1) wd, ht]);
set(axs(2), 'Position', [x(2) y(1) wd, ht]);
set(axs(3), 'Position', [x(3) y(1) wd, ht]);
set(axs(4), 'Position', [x(4) y(1) wd, ht]);
set(axs(5), 'Position', [x(5) y(1) wd, ht]);
set(axs(6), 'Position', [x(1) y(2) wd, ht]);
set(axs(7), 'Position', [x(2) y(2) wd, ht]);
set(axs(8), 'Position', [x(3) y(2) wd, ht]);
set(axs(9), 'Position', [x(4) y(2) wd, ht]);
set(axs(10), 'Position', [x(5) y(2) wd, ht]);

set(axs([2:5, 7:10]),'YTickLabel', []);

panel_identifier(2,5, [.27 .07], .04, fsz2, 'rows', axs, '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS10 Indi_HueHist', [25 11], 'png', punits, 300);

%% Figure S11 | Afterimage Stimuli and Models in DKL-Space.
close all; clc;
[SIM, fh] = mod_afterimage_modeller(MunChips, 'dkl');

% REFORMAT ----------------------------------------------------------------
fh.Color = 'w';
fh.InvertHardcopy = 'off';
axs = get(fh, 'Children');
axs = flipud(axs);
axs(2).YTickLabel = [];
xlabel(axs(2:5),'');
ylabel(axs(2:5),'');

% POSITION & SIZE
wd = .18; ht = .72; y = 0.15;
x = [.03 .225 .42 .615 .81];
set(axs(1), 'Position', [x(1) y wd, ht]);
set(axs(2), 'Position', [x(2) y wd, ht]);
set(axs(3), 'Position', [x(3) y wd, ht]);
set(axs(4), 'Position', [x(4) y wd, ht]);
set(axs(5), 'Position', [x(5) y wd, ht]);

panel_identifier(1,5, [.12 .05], .05, fsz2, 'rows', axs, '.)');

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS11 dklmodels', [30 7], 'png', punits, 300);

%% Figure S12 | Results for Experiment 2b.
close all; clc;
RES = mod_exp2_analysor(AGG2b, 'MAIN');
RES = mod_exp2_analysor(AGG2b, 'Individual Adjustments');
RES = mod_exp2_analysor(AGG2b, 'Individual Deviations From Opponency');
RES = mod_exp2_analysor(AGG2b, 'Individual Hue Histograms');

% REFORMAT ----------------------------------------------------------------
axs = get(gcf, 'Children'); axs = flipud(axs);
[fig, axs] = subplot_recombiner([ones(1,3) 2 3 4 2 3 4], 3, 3, [1:3, 1 1 1 2 2 2]);

wd = .25; ht = .25; x = [.06 .4 .74]; y = [.05 .38 .71]; y = fliplr(y);

set(axs(1), 'Position', [x(1) y(1) wd, ht], 'Xlim', [-1.5 1.5], 'Ylim', [-1.3 1.7]);
set(axs(2), 'Position', [x(2) y(1) wd, ht], 'Ylim', [-30 30], 'YTick', -45:15:45);
set(axs(3), 'Position', [x(3) y(1) wd, ht]);

set(axs(4), 'Position', [x(1) y(2) wd, ht], 'Xlim', [-1 1], 'Ylim', [-1 1.5]);
set(axs(5), 'Position', [x(2) y(2) wd, ht], 'Ylim', [-1 1]*30, 'YTick', -45:15:45);
set(axs(6), 'Position', [x(3) y(2) wd, ht], 'Xlim', [-1 1]*2.5, 'Ylim', [-1 1]*2.5);
xlabel(axs(4:6),'');
set(axs(7), 'Position', [x(1) y(3) wd, ht], 'Xlim', [-1.5 1.5], 'Ylim', [-1.2 2.2]);
set(axs(8), 'Position', [x(2) y(3) wd, ht], 'Ylim', [-1 1]*40, 'YTick', -45:15:45);
set(axs(9), 'Position', [x(3) y(3) wd, ht], 'Xlim', [-1 1]*3,  'XTick', -3:3, 'Ylim', [-1 1]*3, 'YTick', -3:3);
title(axs(7:9), 'f1');

panel_identifier(3,3, [.15 .15], .05, fsz2, 'rows', axs, '.)', 1);

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS12_results_in_dkl', [20 20], 'png', punits, 300);

%% Figure S13 | Additional Data on the Effect of Inducer Chroma on the Hue of Afterimages (Experiment 3). 
close all; clc;
mod_exp3_analysor(AGG3, 'supplementary');

% FORMAT ------------------------------------------------------------------
axs = get(gcf, 'Children'); axs = flipud(axs);
ylabel(axs([1 5 9 13]), 'Blue-Yellow [v*]');
ylabel(axs([4,8,12,16]),'');
xlabel(axs(13:14), 'Green-Red [u*]');
xlabel(axs([3 4 7 8 11 12]), '');
title(axs(5:16),'');
title(axs(1), {'INDUCER STIMULI', 'cw'});
title(axs(2), {'ADJUSTMENTS', ' '});
title(axs(3), {'DEVIATION FROM OPPONENT', ' '});
title(axs(4), {'DEVIATION FROM MEAN', ' '});

wd = .19; ht = .2; x = [.05 .3 .55 .8]; y = fliplr([.04 .28 .52 .76]); 
%wd = .195; ht = .16; x = [.04 .29 .54 .79]; y = fliplr([.04 .23 .42 .61 .80]); 
for k = 1:4
    for l = 1:4
        n = (k-1)*4+l;
        set(axs(n), 'Position', [x(l) y(k) wd, ht]);
    end
end
panel_identifier(4, 4, .1, .1, fsz, 'rows', axs, '.)', 1);

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS13 Additional PPs in Exp3', [20 21], 'png', punits, 300);


%% Figure S14 | Determinants of Non-Linearities
close all; clc;
cone_adaptation_illustrator;

% FORMAT ------------------------------------------------------------------
axs = get(gcf, 'Children'); axs = flipud(axs);
wd = .3; ht = .37; 
x = [.06 .39 .72]; y = [.58 .08];
set(axs(1), 'Position', [x(1) y(1) wd, ht]);
set(axs(2), 'Position', [x(2) y(1) wd, ht]);
set(axs(4), 'Position', [x(3) y(1) wd, ht]);
set(axs(5), 'Position', [x(1) y(2) wd, ht]);
set(axs(6), 'Position', [x(2) y(2) wd, ht]);
set(axs(7), 'Position', [x(3) y(2) wd, ht]);
set(axs(3), 'Position', [0.45, 0.5, 0.2, 0.0390]); % legend

panel_identifier(2, 3, .25, 0, fsz, 'rows', axs([1 2 4:7]), '.)', 1);

% SAVE --------------------------------------------------------------------
fig_saver(save_folder, 'figS14 Determinants_of_Non_Linearities', [20 11], 'png', punits, 300);
