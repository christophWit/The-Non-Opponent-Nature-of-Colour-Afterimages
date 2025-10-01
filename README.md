The-Non-Opponent-Nature-of-Colour-Afterimages

This repository (https://github.com/christophWit/The-Non-Opponent-Nature-of-Colour-Afterimages) accompanies the article "The Non-Opponent Nature of Colour Afterimages" (Witzel, 2025, Communications Psychology, https://doi.org/10.1038/s44271-025-00331-5). This repository focuses on data analyses in Matlab. Data is also available in Excel format in the linked zenodo repository: https://zenodo.org/uploads/13328099. There are three main scripts and 2 datasets. All functions required for the scripts are either embedded in subfunctions or saved in the subfolder functions. You should thus be able to run each of the 3 scripts without further ado. Details on the files:

DATASETS: 

(1) rawdata.mat = non-aggregated data; it is not completely raw as it has been slightly pre-processed for usability, e.g., participants are relabelled according to the labels in the article, data has been assigned to the three experiments and sub-experiments; stimuli have been organised and sorted, etc. This is the data uploaded to the script aggregate_rawdata.m, which produces the aggregated data aggdata.mat (see below). 
(2) aggdata.mat = aggregated data; this is the data used by the script analyses_and_figures.m; it is obtained by running aggregate_rawdata.m with the input from rawdata.mat (see above). All necessary data preprocessing has been with this data so that analyses and figure production can be run smoothly.

SCRIPTS:

(1) aggregate_rawdata.m = Script that takes the variables from rawdata.mat as input and produces aggdata.mat as output. This script relies on the function in the subfolder "functions." 
(2) analyses_and_figures.m = Script that produces all analyses and figures featuring in the article. Detailed analyses are output in the command window. Figures and tables are saved in the subfolder "output_figures_and_tables." Tables are saved as Excel files, main figures in pdf and png format, supplementary figures in png format. This script relies on the functions in the subfolder "functions." Note: When you run the whole script, it will stop after loading the data (this section) to avoid running through all analyses in one go. To run single analyses, load the data, and then go to the respective section and run the section only (by pressing ctrl + return). If you do want to run the whole script in one go, simply delete "return" at the end of the section and make yourself a tea because it will take a while 😊.  
(3) dkl_proof_checker.m = Independent, hermetic function that illustrates the code for calculating DKL space and the equivalence of the different ways to calculate the isoluminant plane in DKL space (see section on DKL model in the article and supplementary material). It is hermetic as it involves neither input nor output. It is independent as is does not rely on any external functions because all necessary functions are integrated as subfunctions. Thus, it can be run without the data and without the subfolder "functions."
FUNCTIONS: 
This folder contains two types of functions. 
(1) Key functions, some of which are also part of the supplementary material. 
(2) Modules, which are large, independent functions that complete and important task their own. They rely on important functions in the subfolder "functions" but contain all other functions as subfunctions. Modules start with the prefix "mod_" in the function name.

DATA ORGANISATION

Raw data involves the following main variables in columns:
pp = Participant nr; the participant label can be retrieved from the "pps" table.
ss = Experimental session.
bk = Block within the session (if there are several blocks per session, e.g., in Exp1a); note that the blocks of a participant be distributed across several conditions/experiments.
trial = trial number within the block.
col = inducer colour nr; the coordinates of the inducer may be found in the "stimuli" variable.
adaptLuv or inducer = The colour coordinates of the inducer. They provide lightness, hue and chroma in Experiment 1, and lightness, opponent axes (u*, v*), hue and chroma in Experiments 2-3.
pos = refers to the position of the wedge along the ring of comparison colours in Experiment 1; as the position was randomised, it is not relevant for the analyses.
resp and resp_hue = refers to the comparison colour chosen in Experiment1; resp gives the number of the selected comparison, resp_hue the corresponding hue of the selected comparison.
rt = response time in seconds.
adj_Luv = These are the coordinates of the adjusted colour in Experiments 2-3 instead of the responses resp and resp_hue in Experiment 1. Coordinates are provided as for the Inducer: (1) Lightness, (2) u*, (3) u*, (4) hue, and (5) chroma.
