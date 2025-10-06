function [AGG1a, AGG1b] = mod_exp1_aggregator(EXP1a, EXP1b, CN, MunChips)
%2018.08.16 * as subfunction [cw]
%2025.07.07 polished for publication [cw]

header('AGGREGATION - EXPERIMENT 1',2);

%% ******************************* SETTINGS *******************************
SETTINGS.agg_stat           = 'median'; % Statistics to aggregate WITHIN participants (only EEP1.a).
SETTINGS.huefrq_res         = 15; % 15-degree steps produce 24 bins
SETTINGS.huefrq_smooth_n    = 3; % Number of adjacent hue frequency bins to smooth over
SETTINGS.sim_chroma         = 'halfinducer'; % The chroma use for calculating complementary colours
SETTINGS.sim_noise_n        = 1000; % Number of noise simulations for hue frequencies

%% ******************************** METHOD ********************************
header('METHOD',3);

%% PARTICIPANTS
header('PARTICIPANTS',4);

% EXPERIMENT 1a -----------------------------------------------------------
fprintf('\nEXPERIMENT 1a\n');
AGG1a.pps.n = size(EXP1a.pps,1);
AGG1a.pps.sex = tabulator(cellstr(EXP1a.pps.sex));
AGG1a.pps.age = [nanmean(EXP1a.pps.age), std(EXP1a.pps.age)];
fprintf('N: %d\n', AGG1a.pps.n);
fprintf('SEX: %d women + %d men\n', AGG1a.pps.sex.F(1),AGG1a.pps.sex.F(2));
fprintf('AGE: %.2f+%.2f\n', AGG1a.pps.age(1), AGG1a.pps.age(2));

% EXPERIMENT 1b -----------------------------------------------------------
fprintf('\nEXPERIMENT 1b\n');
AGG1b.pps.n = size(EXP1b.pps,1);
AGG1b.pps.sex = tabulator(cellstr(EXP1b.pps.sex));
AGG1b.pps.age = [nanmean(EXP1b.pps.age), std(EXP1b.pps.age)];
fprintf('N: %d\n', AGG1b.pps.n);
fprintf('SEX: %d women + %d men\n', AGG1b.pps.sex.F(1),AGG1b.pps.sex.F(2));
fprintf('AGE: %.2f+%.2f\n', AGG1b.pps.age(1), AGG1b.pps.age(2));

%% APPARATUS & STIMULI
AGG1a.stimuli = stimulus_polisher(EXP1a.stimuli);
AGG1b.stimuli = stimulus_polisher(EXP1b.stimuli);

header('APPARATUS',4);
fprintf('\nMONITOR IN EXPERIMENT 1a\n');
fprintf('%d bit, xyY for RGB = \n', AGG1a.stimuli.mon.bit);
disp(AGG1a.stimuli.mon.xyY);
fprintf('\nMONITOR IN EXPERIMENT 1b\n');
AGG1b.stimuli.mon.bit = log2(AGG1b.stimuli.mon.ldt(end)+1);
fprintf('%d bit, xyY for RGB = \n', AGG1b.stimuli.mon.bit);
disp(AGG1b.stimuli.mon.xyY);

header('STIMULI',4);
fprintf('\nWHITE-POINT & BACKGROUND IN EXP1a\n');
disp(AGG1a.stimuli.wp.xyY);
disp(AGG1a.stimuli.bg.Luv);

fprintf('\nWHITE-POINT & BACKGROUND IN EXP1b\n');
disp(AGG1b.stimuli.wp.xyY);
disp(AGG1b.stimuli.bg.Luv);


%% ***************************** AGGREGATION ******************************

header('DATA AGGREGATION',3);
Inducers1a = table([AGG1a.stimuli.inducer.Luv(:,1), AGG1a.stimuli.inducer.Luv_pol], 'VariableNames', {'inducer'});
Inducers1b = table([AGG1b.stimuli.inducer.Luv(:,1), AGG1b.stimuli.inducer.Luv_pol], 'VariableNames', {'inducer'});

%% EXPERIMENT 1a

header('Aggregate Experiment 1a',4);

% AFTERIMAGES -------------------------------------------------------------
AGG1a.afteri        = aggregation_helper(EXP1a.main, SETTINGS.agg_stat);
AGG1a.afteri.agg    = across_participants_aggregator(AGG1a.afteri.indi, AGG1a.stimuli);
AGG1a.afteri.agg    = [Inducers1a, AGG1a.afteri.agg];

% TYPICALITY (prototypes of categories) -----------------------------------
AGG1a.typi.cn   = EXP1a.stimuli.typi_lbls';
AGG1a.typi      = aggregation_helper(EXP1a.typi, SETTINGS.agg_stat);
AGG1a.typi.agg  = across_participants_aggregator(AGG1a.typi.indi, AGG1a.stimuli);
AGG1a.typi.agg  = [Inducers1a, AGG1a.typi.agg];

% AFTERIMAGES -------------------------------------------------------------
AGG1a.discri        = aggregation_helper(EXP1a.discri, SETTINGS.agg_stat);
AGG1a.discri.agg    = across_participants_aggregator(AGG1a.discri.indi, AGG1a.stimuli);
AGG1a.discri.agg    = [Inducers1a, AGG1a.discri.agg];

% HERING COLOURS (red, yellow, green, blue) -------------------------------
% defined by highly saturated typical inducer colours
[~, ~, indsH] = intersect({'red', 'yellow', 'green', 'blue'}, AGG1a.stimuli.lbl, 'stable');
AGG1a.dir.HeringRYGB = AGG1a.stimuli.inducer.Luv_pol(indsH,1);

%% EXPERIMENT 1b

header('Aggregate Experiment 1b',4);

% AFTERIMAGES -------------------------------------------------------------
AGG1b.afteri        = aggregation_helper(EXP1b.main, SETTINGS.agg_stat);
AGG1b.afteri.agg    = across_participants_aggregator(AGG1b.afteri.indi, AGG1b.stimuli);
AGG1b.afteri.agg    = [Inducers1b, AGG1b.afteri.agg];

% HUE FREQUENCIES (HISTOGRAM) ---------------------------------------------
AGG1b.afteri.hue_frq = hue_frequency_calculator(AGG1b.afteri.indi', SETTINGS.huefrq_res, SETTINGS.huefrq_smooth_n);
AGG1b.afteri.hue_frq = [Inducers1b, AGG1b.afteri.hue_frq];

% COLOUR NAMING -----------------------------------------------------------
AGG1b.cn = mod_cn_aggregator(EXP1b, CN);

% HERING COLOURS (red, yellow, green, blue) -------------------------------
% defined as prototypes of red, yellow, green, and blue
[~, ~, indsH] = intersect({'red', 'yellow', 'green', 'blue'}, AGG1b.cn.colnames, 'stable');
AGG1b.dir.HeringRYGB = AGG1b.cn.typM(indsH,1);

%% ************************ MODELS & PREDICTIONS **************************

header('MODELS & PREDICTIONS',3);

%% PREDICTIONS: Calculate complementary colours

header('COMPLEMENTARY COLOURS FOR EXP1.a', 4);
[AGG1a.model.Hue, AGG1a.model.Chroma, AGG1a.model.Lum, AGG1a.model.dkl_k, AGG1a.model.Munsell] = complementary_calculator(...
    AGG1a.stimuli, 'Luv', MunChips, AGG1a.dir.HeringRYGB, SETTINGS.sim_chroma);

header('COMPLEMENTARY COLOURS FOR EXP1.b', 4);
[AGG1b.model.Hue, AGG1b.model.Chroma, AGG1b.model.Lum, AGG1b.model.dkl_k, AGG1b.model.Munsell] = complementary_calculator(...
    AGG1b.stimuli, 'Luv', MunChips, AGG1b.dir.HeringRYGB, SETTINGS.sim_chroma);

%% PREDICTION ERRORS: Differences between model predictions and data

header('Prediction Errors', 4);
[AGG1a.prederr.hue, AGG1a.prederr.hueM] = predictionerrorcalculator_hue(AGG1a.afteri.indi, AGG1a.model.Hue);
[AGG1b.prederr.hue, AGG1b.prederr.hueM] = predictionerrorcalculator_hue(AGG1b.afteri.indi, AGG1b.model.Hue);

% DEVIATIONS FROM OPPONENCY -----------------------------------------------
AGG1a.OppoDev               = table;
AGG1a.OppoDev.adj           = AGG1a.prederr.hueM.dkl;
AGG1a.OppoDev.adj_smooth    = smoothdata(AGG1a.OppoDev.adj(:,1),'movmean', SETTINGS.huefrq_smooth_n);
AGG1a.OppoDev.sim           = circularsubtractor(AGG1a.model.Hue.ConeCon, AGG1a.model.Hue.dkl, 360,'signed');
AGG1a.OppoDev.diff          = circularsubtractor(AGG1a.prederr.hueM.ConeCon(:,1), AGG1a.prederr.hueM.dkl(:,1), 360,'signed');
AGG1a.OppoDev               = [Inducers1a, AGG1a.OppoDev];

% DEVIATIONS FROM OPPONENCY -----------------------------------------------
AGG1b.OppoDev               = table;
AGG1b.OppoDev.adj           = AGG1b.prederr.hueM.dkl;
AGG1b.OppoDev.adj_smooth    = smoothdata(AGG1b.OppoDev.adj(:,1),'movmean', SETTINGS.huefrq_smooth_n);
AGG1b.OppoDev.sim           = circularsubtractor(AGG1b.model.Hue.ConeCon, AGG1b.model.Hue.dkl, 360,'signed');
AGG1b.OppoDev.diff          = circularsubtractor(AGG1b.prederr.hueM.ConeCon(:,1), AGG1b.prederr.hueM.dkl(:,1), 360,'signed');
AGG1b.OppoDev               = [Inducers1b, AGG1b.OppoDev];

%% Higher resolution for plot
header('COMPLEMENTARY COLOURS AT HIGHER RESOLUTION (EXP1b)', 4);
stimuli360 = AGG1b.stimuli;
Luv_pol = [(0:360)', ones(361,1)*unique(AGG1b.stimuli.inducer.Luv_pol(:,2))];
% [u,v] = pol2cart(...
%     deg2rad(stimuli360.inducer.Luv_pol(:,1)),...
%     stimuli360.inducer.Luv_pol(:,2));
L =  ones(361,1)*unique(AGG1b.stimuli.inducer.Luv(:,1));
% stimuli360.inducer.Luv = [L u v];
stimuli360.inducer = colourconverter([L Luv_pol], 'Luv_pol', 2,...
    stimuli360.wp.xyY, stimuli360.sens, stimuli360.cmf,...
    stimuli360.mon.xyY, stimuli360.mon.ldt, stimuli360.mon.oog,[],stimuli360.mon.rgb_max,...
    stimuli360.bg.xyY(3),0);
Hue = complementary_calculator(...
    stimuli360, 'Luv', MunChips, AGG1b.dir.HeringRYGB, SETTINGS.sim_chroma);

AGG1b.model360 = table;
AGG1b.model360.inducer  = stimuli360.inducer.Luv_pol;
AGG1b.model360.ConeCon  = Hue.ConeCon;
AGG1b.model360.DKL      = Hue.dkl;
AGG1b.model360.OppoDev  = circularsubtractor(Hue.ConeCon(:,1), Hue.dkl(:,1), 360,'signed');

%% SIMULATE HUE FREQUENCIES WITH NOISE
fprintf('\nSimulating Hue frequencies with noise (Exp1b only)\n');
hueSD   = mean(AGG1b.afteri.agg.hue_std);
noisy_azi = noise_maker([AGG1b.model.Hue.ConeCon, AGG1b.model.Chroma.ConeCon], hueSD, 0, SETTINGS.sim_noise_n);
AGG1b.model.hueF = hue_frequency_calculator(noisy_azi, SETTINGS.huefrq_res, SETTINGS.huefrq_smooth_n);
AGG1b.model.hueF = [Inducers1b, AGG1b.model.hueF];

% FOR FIGURES ---------------------------------------------------------
AGG1a.dir.cones = simple_cone_simulator(AGG1a.stimuli);
AGG1b.dir.cones = simple_cone_simulator(AGG1b.stimuli);

%% ****************************** SUBMODULES ******************************

%% across_participants_aggregator
function AGG = across_participants_aggregator(INDI, stimuli)
%2025.07.08 polished for publication [cw]

wp      = stimuli.wp.xyY;
mon_xyY = stimuli.mon.xyY;
mon_ldt = stimuli.mon.ldt;
mon_oog = stimuli.mon.oog;
rgb_max = stimuli.mon.rgb_max;
comp    = stimuli.comp;
sens    = stimuli.sens;
cmf     = stimuli.cmf;
bgY     = stimuli.bg.xyY(3);
N = sum(~isnan(INDI))';

L       = comp(:,1,1);
chro    = comp(:,1,3);
hue_M   = cstats(INDI, 'nanmean')';
hue_std = cstats(INDI, 'nanstd')';
hue_sem = cstats(INDI, 'ste')';
[x,y] = pol2cart(deg2rad(hue_M), chro);
Luv = [L,x,y];
TMP = colourconverter(Luv, 'Luv', 2, wp, sens, cmf,...
    mon_xyY, mon_ldt, mon_oog,[],rgb_max,...
    bgY, 0);
rgb = TMP.rgb/rgb_max;
AGG = table(hue_M, hue_std, hue_sem, N, Luv, chro, rgb);

%% aggregation_helper
function AGG = aggregation_helper(RAW, agg_stat)
% INPUT -------------------------------------------------------------------
% RAW correponds to raw data from afterimage and typical measurements.
% agg_stat defines aggregation across repeated measurements (within pps).
% OUTPUT ------------------------------------------------------------------
% AGG corresponds to aggregated data for afteri
% 2020.05.28 [cw]
% 2020.06.04 added overall frequency [cw]
% 2025.07.08 polished for publication [cw]

colours = unique(RAW.col);
col_n = numel(colours);
pps = unique(RAW.pp);
pp_n = numel(pps);

% AVERAGES ----------------------------------------------------------------
F = zeros(10, col_n);
rF = zeros(10, col_n);
AGG.resp_freq = table(F, rF);
for col = 1:col_n
    inds1 = RAW.col == col;
    tab = tabulator(RAW.resp(inds1));
    AGG.resp_freq.F(tab.id,col)   = tab.F;
    AGG.resp_freq.rF(tab.id,col)  = tab.rF;
    for pp = 1:pp_n
        inds2 = RAW.pp == pp;
        inds  = inds1 & inds2;
        AGG.skip(pp, col)   = sum(RAW.resp(inds)==10);
        AGG.N(pp, col)      = sum(~isnan(RAW.resp_hue(inds)));
        AGG.missed(pp, col) = sum(isnan(RAW.resp(inds)));
        switch lower(agg_stat)
            case 'median'
                AGG.indi(pp,col)    = cstats(RAW.resp_hue(inds), 'nanmedian');
                AGG.rt(pp,col)      = nanmedian(RAW.rt(inds));
            case {'mean', 'average'}
                AGG.indi(pp,col)    = cstats(RAW.resp_hue(inds), 'nanmean');
                AGG.rt(pp,col)      = nanmean(RAW.rt(inds));
        end
    end
end

%% boundary_aggregator
function [new_low, new_up] = boundary_aggregator(bound_low0, bound_up0)
%2020.05.27 [cw]
%2025.09.24 polished for publication [cw]

% DELETE COMPLETELY MISSING DATA FOR PROCESSING (ADD LATER) ---------------
pp_n0 = size(bound_low0,1);
nans = all(isnan(bound_low0),2);
bound_low = bound_low0(~nans,:);
bound_up = bound_up0(~nans,:);

pp_n = size(bound_low,1);

% DETERMINE ORDER OF BOUNDARIES BASED ON AVERAGE BOUNDARIES ---------------
Mctg_low = cstats(bound_low)';
Mctg_up = cstats(bound_up)';

Mctg = cstats([Mctg_low,Mctg_up]')';
[~, sort_inds] = sortrows(Mctg);

% REPLACE NAN BY ZERO WIDTH CATEGORIES ------------------------------------

% Sort according to averages:
sorted_low  = bound_low(:,sort_inds);
sorted_up   = bound_up(:,sort_inds);

% Identify NaNs:
nan_low = isnan(sorted_low);
nan_up = isnan(sorted_up);

% Shift indices to adjacent categories:
nan_low2 = [false(pp_n,1),nan_low]; % 1 back
nan_up2 = [nan_up,false(pp_n,1)]; % 1 further

% Append fist (low) and last (up):
sorted_low2 = [sorted_low,sorted_low(:,1)]; % FOR UP!
sorted_up2 = [sorted_up(:,end),sorted_up]; % FOR LOW!

% Replace NaNs by adjacent non-nans:
sorted_low(nan_low) = sorted_low2(nan_low2);
sorted_up(nan_up) = sorted_up2(nan_up2);

% UNSORT (Restablish original order) --------------------------------------
[~, unsort_inds] =sortrows(sort_inds);

% new_low0 = sorted_low(:,unsort_inds_low);
% new_up0 = sorted_up(:,unsort_inds_up);
new_low0 = sorted_low(:,unsort_inds);
new_up0 = sorted_up(:,unsort_inds);

if sum(isnan(sorted_low(:)) | isnan(sorted_up(:))) > 0
    [new_low0, new_up0] = boundary_aggregator(new_low0, new_up0);   
end

% ADD MISSING ROWS --------------------------------------------------------
cn_n = size(new_low0,2);
new_low = NaN(pp_n0,cn_n);
new_low(~nans,:) = new_low0;
new_up(~nans,:) = new_up0;

%% mod_cn_aggregator
function CNagg = mod_cn_aggregator(AFT360, CN360)
%2020.05.28

agg = aggregator(CN360.categories, {'pp','col'}, 'response'); 
N = size(agg,1);
[~,frq] = mode(agg); 
rfrq = frq/N;
[frq2(:,2), frq2(:,1)] = cinterpolator(CN360.sti.colours.Luv_pol(:,1), frq', (0:5:355)', 360);
[rfrq2(:,2), rfrq2(:,1)] = cinterpolator(CN360.sti.colours.Luv_pol(:,1), rfrq', (0:5:355)', 360);
CNagg.frqtab = [frq2,rfrq2(:,2)];
CNagg.agg_bounds = catboundfinder(agg','orientationdata',CN360.sti.colours.Luv_pol(:,1), 'cat_ids', 1:size(CN360.colnames,2));


% MATCH PPS OF AFT360 & CN360 ---------------------------------------------
[~,inds1,inds2] = intersect(AFT360.pps.pp,CN360.pps.pp, 'stable');

% PROTOTYPES:
N360 = size(AFT360.pps,1);
CNagg.typi_lbls = CN360.adj_labels;
cn_n = size(CNagg.typi_lbls,2);
indi_agg = NaN(N360,cn_n,2);
TMP = aggregator(CN360.prototypes, {'bk','pp','cn'}, 'azi');
for cn = 1:10
    tmp2 = cdissolver(TMP(:,:,cn));
    Medi = median(tmp2);
    out0 = circularsubtractor(tmp2,Medi,360, 'mindist')>120;
    if any(out0(:))
        tmp2(out0) = NaN;
        fprintf('WARNING: Replaced %d for %d with NaN\n',sum(out0(:)),cn); 
    end
    out(:,:,cn) = out0;
    [indi_agg0(:,cn,1), indi_agg0(:,cn,2)]  = cstats(tmp2,{'nanmean', 'ste'}); 
end
indi_agg(inds1,:,:) = indi_agg0(inds2,:,:); % order according to pps in AFTERI360
CNagg.typi      = indi_agg;
CNagg.outliers  = out;
[cn_typM(:,1), cn_typM(:,2)] = cstats(indi_agg(:,:,1),{'nanmean', 'ste'}); 
CNagg.typM = cn_typM;

% BOUNDARIES:
CNagg.colnames = CN360.colnames(1,:);
cn2_n = size(CNagg.colnames,2);

% match pp ids:
bounds = NaN(N360,cn2_n,2);
bounds0 = cat(3,CN360.boundaries.low,CN360.boundaries.up);
bounds(inds1,:,:) = bounds0(inds2,:,:);
CNagg.bounds = bounds;

% Compensate for missing boundaries to avoid distortion of consensus boundaries:
[new_low, new_up] = boundary_aggregator(bounds0(:,:,1),bounds0(:,:,2));
[Mbound_low(1,:),Mbound_low(2,:)] = cstats(new_low, {'nanmean','ste'});
[Mbound_up(1,:),Mbound_up(2,:)] = cstats(new_up, {'nanmean','ste'});

% DELETE IF AGGREAGED BORDERS ARE NONESENSE (BROWN):
CNagg.Mbound = cat(3,Mbound_low', Mbound_up');

%% hue_frequency_calculator
function hueF = hue_frequency_calculator(azis, step, smoother)
% 2025.07.16 [cw]

if nargin < 3
    smoother = 9;
end

azis = mod(azis,360); % Just in case...

edge1 = step/2; 
edgeEnd = 360-edge1; 
edges = [0,edge1:step:edgeEnd,360]; % Extreme bins area half as large to create a bin at 0
frq0 = histcounts(azis, edges);
frq = frq0(1:end-1)';
frq(1)=frq0(1)+frq0(end);
n = size(azis,2);
rF = frq/n;
N  = ones(size(frq,1),1)*n;
smooth_frq  = smoothdata(frq,'movmean',smoother);
smooth_rF   = smoothdata(rF,'movmean',smoother);
hueF = table(N, frq, rF, smooth_frq, smooth_rF);

%% predictionerrorcalculator_hue
function [pe_hue, pe_hueM] = predictionerrorcalculator_hue(match_hue, OppoHue)
% 2025.07.23 [cw]

pe_hue      = table; 
pe_hueM     = table;

varnames = OppoHue.Properties.VariableNames;
for vn = 1:numel(varnames)
    vrnm = varnames{vn};

    % HUE -----------------------------------------------------------------
    oppo_hue = OppoHue.(vrnm);
    pe_hue.(vrnm) = circularsubtractor(match_hue, ones(size(match_hue,1),1)*oppo_hue', 360, 'signed');
    
    % AGGREGATE:
    Mhue = nanmean(pe_hue.(vrnm), 1);
    semH = ste(pe_hue.(vrnm), 0, 1);
    pe_hueM.(vrnm)(:,1)  = Mhue';
    pe_hueM.(vrnm)(:,2)  = semH';
    

end

%% simple_cone_simulator
function CONES = simple_cone_simulator(stimuli, k)
%2025.07.20

if nargin < 2
    k = 10; % S-Cone factor: produces negative chromaticities if too high
end

colspace    = stimuli.src_space;
wp          = stimuli.wp.xyY;

% CREATE INITIAL LMS EXCITATIONS ------------------------------------------
% Luminance of isoluminant background and stimuli defines the limits
% because S cone direction must have L- and M-cones to reach the luminance
% of background and stimuli. 
bg_lms = XYZ2lms(stimuli.bg.XYZ, stimuli.sens, stimuli.cmf);
mx = sum(bg_lms(1:2));
lms = [...
    mx 0 0;...
    0 mx 0;...
    bg_lms(1) bg_lms(2) k*bg_lms(3)];

% CONVERT TO OTHER COLOURS SPACES -----------------------------------------
%XYZ = lms2XYZ(lms, stimuli.sens, stimuli.cmf);
CONES = colourconverter(lms, 'lms', 2,...
    wp, stimuli.sens, stimuli.cmf,...
    stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max,...
    stimuli.bg.xyY(3), 0);

%% stimulus_polisher
function agg_stimuli = stimulus_polisher(exp_stimuli)
% 2025.07.08 [cw]

agg_stimuli             = exp_stimuli;

% CONVERT OLD DATASET TO TABLE FORMAT -------------------------------------
% agg_stimuli.inducer = dataset2table(agg_stimuli.inducer);
% agg_stimuli.bg = dataset2table(agg_stimuli.bg);
% agg_stimuli.wp = dataset2table(agg_stimuli.wp);
% agg_stimuli.inducer.Properties.RowNames = agg_stimuli.lbl;

% CONE SENSITIVITIES ------------------------------------------------------
agg_stimuli.sens = 'ss';
agg_stimuli.cmf  = '1931';

% MONITOR ----------------------------------------------------------------- 

% FORMAT:
agg_stimuli.mon.rgb_max = agg_stimuli.mon.ldt(end); % Maximum RGB
agg_stimuli.mon.bit     = log2(agg_stimuli.mon.rgb_max+1); % Bits

% CALCULATE ACCURATE MONITOR LMS FOR DKL CONVERSION:
agg_stimuli.mon.lms = XYZ2lms(xyY2XYZ(agg_stimuli.mon.xyY), agg_stimuli.sens, agg_stimuli.cmf);

% RELATIVE RGB (to display on conventional monitor):
agg_stimuli.inducer.rgb = agg_stimuli.inducer.rgb/agg_stimuli.mon.rgb_max;
agg_stimuli.bg.rgb      = agg_stimuli.bg.rgb/agg_stimuli.mon.rgb_max; 

% COMPARISON RGB ----------------------------------------------------------
TMP = colourconverter(agg_stimuli.comp, 'Luv_pol', 3,...
    agg_stimuli.wp.xyY, agg_stimuli.sens, agg_stimuli.cmf,...
    agg_stimuli.mon.xyY, agg_stimuli.mon.ldt, agg_stimuli.mon.oog,[], agg_stimuli.mon.rgb_max,...
    agg_stimuli.bg.xyY(3), 0);
agg_stimuli.comp_rgb = TMP.rgb/agg_stimuli.mon.rgb_max;

% TARGET FOR DISCRIM ------------------------------------------------------
tLuv = shiftdim(agg_stimuli.comp(:,5,:),2)';
agg_stimuli.discri_target = colourconverter(tLuv, 'Luv_pol', 2,...
    agg_stimuli.wp.xyY, agg_stimuli.sens, agg_stimuli.cmf,...
    agg_stimuli.mon.xyY, agg_stimuli.mon.ldt, agg_stimuli.mon.oog,[], agg_stimuli.mon.rgb_max,...
    agg_stimuli.bg.xyY(3), 0);
agg_stimuli.discri_target.rgb = agg_stimuli.discri_target.rgb/agg_stimuli.mon.rgb_max; 

%% ***************************** SUBFUNCTIONS *****************************

%% aggregator
function varargout = aggregator(DATA, grpvar_cell, ourvar_cell)
% NOTE: remaining repetitions are filed under the first dimension (rows).
% See: nr(1) = nr(1) * N/prod(nr); 
% 2018.08.19 * [cw]


if isstr(ourvar_cell)
    ourvar_cell = {ourvar_cell};
end

grpvar_cell = flipud(grpvar_cell(:));

for k = 1:numel(grpvar_cell)
    nr(1,k) = numel(unique(DATA.(grpvar_cell{k})));
end
nr = fliplr(nr);

N = size(DATA,1);
nr(1) = nr(1) * N/prod(nr);

for k = 1:numel(ourvar_cell)
    sorted = sortrows(DATA, grpvar_cell);
    varargout{k}  = reshape(sorted.(ourvar_cell{k}),nr);
end

%% cdissolver
function outM = cdissolver(inM, onecycle)
% 2024.07.14 [cw]

if nargin < 2
    onecycle = 360;
end

%inM = mod(inM, onecycle);
for cl = 1:size(inM,2)
    refM1 = ones(size(inM,1),1) * min(inM(:,cl));
    refM2 = ones(size(inM,1),1) * nanmedian(inM(:,cl));
    dffM1 = circularsubtractor(inM(:,cl), refM1, onecycle, 'signed');
    dffM2 = circularsubtractor(inM(:,cl), refM2, onecycle, 'signed');
    if nanmean(abs(dffM2)) < nanmean(abs(dffM1))
        refM = refM2; dffM = dffM2;
    else
        refM = refM1; dffM = dffM1;
    end
%    refM = refM2; dffM = dffM2;
    outM(:,cl) = refM+dffM;
end

%% cinterpolator
function [newy, newx] = cinterpolator(X1, Y1, X2, onecycle)
% 2012sep24-25 * [cw]

if nargin < 4
    onecycle= 360;    
end

% check input:
if size(X1,1) ~= size(Y1,1) | size(X1,2) ~= size(Y1,2)
    error('cinterpolator: x1 and y1 don''t match!');
end

% Delete NaNs:
inds1 = ~isnan(X1);
inds2 = ~isnan(Y1);
inds = inds1&inds2;
X1 = X1(inds);
Y1 = Y1(inds);

% Map to cycle:
xy(:,1) = mod(X1,onecycle);
xy(:,2) = Y1;
xy = sortrows(xy,1);
x1 = xy(:,1);
y1 = xy(:,2);
x2 = mod(X2, onecycle);
%x2 = sort(x2);

% Eliminate multiple x-values by averaging the corresponding y-values:
[x1 y1] = double_killer(x1,y1);
% [x2 dbl_n] = double_killer(x2);
% if dbl_n > 0
%     fprintf('WARNING - cinterpolator: interpolation involved double target x2; deleted doubles!\n');
% end

% Extend x and y so as to obtain a full cycle:
xex = [...
    x1-onecycle;...
    x1;...
    x1+onecycle];
yex = [y1; y1; y1];
% Actual interpolation:
y2 = interp1(xex, yex, x2, 'linear');
newx = x2;
newy = y2;

%% double_killer
function [newx newy dbl_n] = double_killer(x,y)
%2012sep24 [cw]
newx = x;
qy = 0;
if nargin > 1
    qy = 1;
    newy = y;
end
tab = tabulate(newx);
inds = find(tab(:,2) > 1);
dbls = tab(inds,1);
dbl_n = numel(dbls);
for nr = 1:dbl_n
    inds = newx == dbls(nr);
    newx(inds(2:end)) = [];
    if qy
        newy(inds(1)) = mean(y1(inds));
        newy(inds(2:end)) = [];
    end
end

if ~qy
    newy = dbl_n;
end

%% circularsubtractor
function Difference = circularsubtractor(Minuend, Subtrahend, ONECYCLE, dtype)
% 2012aug23-24 [cw]
        
if nargin < 4
    dtype = 'regular';
    if nargin < 3
        ONECYCLE = 360;
    end
end

% map angles to 0:oncecycle:
Minuend     = mod(Minuend, ONECYCLE);
Subtrahend  = mod(Subtrahend, ONECYCLE);

% Calculate difference(s):
regular     = mod(Minuend - Subtrahend, ONECYCLE);
clockwise   = mod(Subtrahend - Minuend, ONECYCLE);

switch lower(dtype)
    case 'regular' % Regular difference (counterclockwise)
        Difference = mod(regular, ONECYCLE);
    case 'signed'
        Difference = NaN(size(Minuend,1), size(Minuend,2));
        inverse = -clockwise;
        inds = regular <= clockwise; % find the closed
        Difference(inds) = regular(inds);
        Difference(~inds) = inverse(~inds);        
    case 'clockwise' % Clockwise difference (inverse)
        Difference = mod(clockwise, ONECYCLE);
    case {'mindist', 'maxdist'} % Minimal distance
        regular = mod(regular, ONECYCLE);
        clockwise = mod(clockwise, ONECYCLE);
        switch lower(dtype)
            case 'mindist' % Minimal distance
                Difference = min(abs(regular),abs(clockwise));
            case 'maxdist' % Maximal distance
                Difference = max(abs(regular),abs(clockwise));
        end
end

%% cstats
function varargout = cstats(inM, stype, onecycle, remap)
% 2015.10.30 [cw]

if nargin < 4
    remap = 1;
    if nargin < 3
        onecycle = 360;
        if nargin < 2
            stype = {'nanmean', 'nanmin', 'nanmax', 'range', 'nanstd', 'ste'};
        end
    end
end

if ischar(stype)
    stype = {stype};
end
dc = cdissolver(inM, onecycle);
for nr = 1:numel(stype)
        prelim = eval([stype{nr}, '(dc);']);
        if remap % map on 0:onecycle
            varargout{nr} = mod(prelim, onecycle); 
        else
            varargout{nr} = prelim;
        end            
end

%% header
function header(txt, lvl, wdth, stl)
%2015.02.04 [cw]

if nargin < 4
    stl = 'classic';
    if nargin < 3
        wdth = 100;
        if nargin < 2
            lvl = 1;
            if nargin < 1
                txt = mfilename;
            end
        end
    end    
end

switch lower(stl)
    case {'classic'}
        if lvl == 1
            xx = floor((wdth-numel(txt))/2);
            fprintf([linemaker(wdth, '_'), '\n\n']);
            fprintf('%s%s\n', linemaker(xx, ' '), txt);
            fprintf([linemaker(wdth, '_'), '\n\n']);
        elseif lvl == 2
            xx = floor((wdth-numel(txt))/2);
            fprintf([linemaker(wdth, '-'), '\n']);
            fprintf('%s%s\n', linemaker(xx, ' '), txt);
            fprintf([linemaker(wdth, '-'), '\n\n']);
        elseif lvl == 3
            xx = floor((wdth-numel(txt))/2);
            fprintf('\n%s%s\n', linemaker(xx, ' '), txt);
            fprintf([linemaker(wdth, '.'), '\n\n']);
        elseif lvl == 4
            fprintf('\n\n----- %s -----\n\n', txt);
        end
end

%% linemaker
function ln = linemaker(width, style)
% 2011apr06 [cw]

if nargin < 2
    style = '_';
    if nargin < 1
        width = 100;
    end
end
ln(1:width) = style;

%% lms2XYZ
function XYZ = lms2XYZ(lms, fndmtls, cmf)
% 2014dec14 * [cw]
% NOT TESTED YET!!!!!

if nargin < 3
    cmf = 'judd';
    if nargin < 2
        fndmtls = 'smithpokorny';
    end
end

M = coneconversionmatrix(fndmtls, cmf);
XYZ = lms/(M'); %apparently more precise than: XYZ = lms*inv(M')

%% coneconversionmatrix
function M = coneconversionmatrix(msrments, cmf)
M = [];
% Conversion matrices towards CIE1931 XYZ taken from:
% Golz, J., & MacLeod, D. I. A. (2003). Colorimetry for CRT displays. J Opt Soc Am A Opt Image Sci Vis, 20(5), 769-781. 

switch lower(msrments)
    case {'smithpokorny', 'sp'}
        % V. C. Smith and J. Pokorny, ‘‘Spectral sensitivity of the
        % foveal cone photopigments between 400 and 500 nm,’’ Vision
        % Res. 15, 161–171 (1975).
        switch lower(cmf)
            case '1931'
                M = [0.15282 0.54383 -0.02795
                    -0.15254 0.45524 0.03355
                    -0.00045 0.00145 0.95449];
                
            case 'judd'
                % Smith & Pokorny cone fundamentals
                % V. C. Smith & J. Pokorny (1975), Vision Res. 15, 161-172.
                %        X          Y       Z [cw]
                M = [0.15514  0.54312  -0.03286    % L alias R
                    -0.15514  0.45684   0.03286    % M alias G
                    0.0      0.0       0.01608];  % S alias B
        end
    case 'smj2'
        % A. Stockman, D. I. A. MacLeod, and N. E. Johnson, ‘‘Spectral
        % sensitivities of the human cones,’’ J. Opt. Soc. Am. A 10,
        % 2491–2521 (1993).
        switch lower(cmf)
            case '1931'
                M = [ 0.18772 0.60445 -0.02517
                    -0.14014 0.43056 0.03773
                    0.02017 -0.04189 1.08472];
        end
    case 'smj10'
        % A. Stockman, D. I. A. MacLeod, and N. E. Johnson, ‘‘Spectral
        % sensitivities of the human cones,’’ J. Opt. Soc. Am. A 10,
        % 2491–2521 (1993).
        switch lower(cmf)
            case '1931'
                M = [ 0.14460 0.62421 -0.00429
                    -0.14506 0.42265 0.05084
                    0.03105 -0.06416 1.10923];
        end
    case {'stockmansharpe', 'ss'}
        % A. Stockman and L. T. Sharpe, ‘‘Spectral sensitivities of the
        % middle- and long-wavelength sensitive cones derived from
        % measurements in observers of known genotype,’’ Vision Res.
        % 40, 1711–1737 (2000).
        switch lower(cmf)
            case '1931'
                M = [0.17156 0.52901 -0.02199
                    -0.15955 0.48553 0.04298
                    0.01916 -0.03989 1.03993];
                % M = [1.94735469 -1.41445123 0.36476327
                %     0.68990272 0.34832189 0
                %     0 0 1.93485343];
        end
end


if isempty(M)
    error('This transformation matrix is not implemented');
end

%% noise_maker
function [noisy_azi, noisy_rad] = noise_maker(AziRad0, hueSD, chrSD, repetitions)
% 2025.04.17

if nargin < 4
    repetitions = 100;
end

sti_n = size(AziRad0,1);
noisy_azi = NaN(sti_n, repetitions);
noisy_rad = NaN(sti_n, repetitions);
for sti = 1:sti_n
    noisy_azi(sti,:) = normrnd(AziRad0(sti,1),hueSD, 1, repetitions);
    noisy_rad(sti,:) = normrnd(AziRad0(sti,2),chrSD, 1, repetitions);
end
noisy_azi = mod(noisy_azi,360);

%% ste
function y = ste(varargin)
% 2019.04.12 [cw]
N = sum(~isnan(varargin{1}),1);
y = sqrt(nanvar(varargin{:}))./sqrt(N);

%% tabulator
function tab = tabulator(RowsPerColumn, deleteF0, sortpercolumn)
% 2024.05.15 [cw]

if nargin < 3
    sortpercolumn = 0;
    if nargin < 2
        deleteF0 = 1;
    end
end

cl_n = size(RowsPerColumn,2);

tab = table;
for cl = 1:cl_n
    tab1 = table;
    data_cl = RowsPerColumn(:,cl);
    if ~iscell(data_cl)
        data_cl(isnan(data_cl)) = [];
    end
    tab0 = tabulate(data_cl);
    if deleteF0
        if ~iscell(data_cl)
            tab0(tab0(:,2)==0,:)= []; % delete 0 frequency cases
        end
    end
    column = ones(size(tab0,1),1)*cl;
    tab1.cl  = column;
    tab1.id  = tab0(:,1);
    if ~iscell(data_cl)
        tab1.F   = tab0(:,2);
        tab1.rF  = tab0(:,3)/100;
    else
        tab1.F   = cell2mat(tab0(:,2));
        tab1.rF  = cell2mat(tab0(:,3))/100;
    end
    if sortpercolumn
        tab1 = sortrows(tab1,'rF','descend');
    end
    if cl == 1
        tab = tab1;
    else
        tab  = [tab;tab1];
    end
end

%% xyY2XYZ
function XYZ = xyY2XYZ(xyY, dim)
%input: x, y, Y values in a row vector or matrix )
% 2012.01.20 [cw]

if nargin < 2
    dim = 2;
end

if dim == 2
    x = xyY(:,1); y = xyY(:,2); Y = xyY(:,3);
elseif dim == 3
    x = xyY(:,:,1); y = xyY(:,:,2); Y = xyY(:,:,3); 
end

X = (Y./y) .* x;
Y = Y;
Z = (Y./y) .* (1-y-x);

if dim == 2
    XYZ = [X Y Z];
elseif dim == 3
    XYZ = cat(3, X, Y, Z);    
end

%% XYZ2lms
function [lms, M] = XYZ2lms(XYZ, fndmtls, cmf)
% Calculates cone excitations (LMS) based on Tristimulus Values (XYZ).
% Conversion matrices from CIE1931 XYZ taken from:
% Golz, J., & MacLeod, D. I. A. (2003). Colorimetry for CRT displays. J Opt Soc Am A Opt Image Sci Vis, 20(5), 769-781. 
% 2020.07.22 [cw]

if nargin < 3
    cmf = 'judd';
    if nargin < 2
        fndmtls = 'smithpokorny';
    end
end

M = coneconversionmatrix(fndmtls, cmf);
lms = XYZ*M';
