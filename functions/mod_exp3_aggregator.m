function AGG = mod_exp3_aggregator(EXP3)
% 2025.07.24 [cw]

SETTINGS.sim_chroma = 'inducer';
SETTINGS.agg_stat   = 'mean';
SETTINGS.colspace   = 'Luv';
AGG = struct;

%% PARTICIPANTS
AGG.pps = EXP3.pps;

%% APPARATUS & STIMULI
for pp = 1:size(AGG.pps,1)
    header(EXP3.pps.pp_lbl{pp}, 4);
    AGG.stimuli{pp} = stimulus_polisher(EXP3.stimuli{pp});
end

%% AGGREGATION
for pp = 1:size(AGG.pps,1)
    header(EXP3.pps.pp_lbl{pp}, 4);
    pp_id = EXP3.pps.pp(pp);
    indsP = EXP3.main.pp == pp_id;
    [AGG.data{pp}, AGG.model{pp}, AGG.model4illu{pp}] = mod_sat_aggregator(EXP3.main(indsP,:), AGG.stimuli{pp}, SETTINGS);
end

%% ************************ SPECIFIC SUBFUNCTIONS *************************

%% mod_sat_aggregator
function [AGG, SIM, ILLU] = mod_sat_aggregator(MAIN, stimuli, SETTINGS)
% 2025.07.25 polished for publication

% INITIALISE OUTPUT -------------------------------------------------------
AGG     = table;
SIM     = struct;
ILLU    = struct;

% SETTINGS ----------------------------------------------------------------
WP = stimuli.wp;
mon = stimuli.mon;
bgY = stimuli.bg.xyY(3);
switch lower(stimuli.src_space) % Always Luv for Experiment 3
    case 'luv'
        bg      = stimuli.bg.Luv;
        Inducer = MAIN.inducer;
        adjvar  = 'adj_Luv';
    case 'dkl'
        bg      = stimuli.bg.dkl;
        Inducer = ind_dkl; 
        adjvar  = 'adj_dkl';
end
Response    = MAIN.(adjvar);

% IDENTIFY + CHECK INDUCERS -----------------------------------------------

% IDENTIFY OOG
if ~isempty(stimuli.oog)
    n = size(stimuli.oog,1);
    for k = 1:n
        killH = round(stimuli.oog.Luv_pol(k,1));
        killC = round(stimuli.oog.Luv_pol(k,2));
        indsH = round(Inducer(:,4)) == killH;
        indsC = round(Inducer(:,5)) == killC;
        kill = indsH & indsC;
        MAIN(kill,:) = [];
        Inducer(kill,:) = [];
        Response(kill,:) = [];
        fprintf('WARNING: %d data points have been deleted for obsolete oog inducer at H = %d & C = %d\n', sum(kill), killH, killC);        
    end
    fprintf('The remaining data is n = %d\n', size(Response,1));        
end

% LIGHTNESS (should be all at L*=70)
L = unique(Inducer(:,1));
if numel(L) > 1 | L ~= bg(:,1)
    disp(L);
    disp(bg(:,1));
    error('Inducer lightness is not as expected.');
else
    fprintf('Inducer & background lightness is: %.1f\n', L);
end

% HUE (varies across pps)
[xtab,~,~,labels] = crosstab(Inducer(:,4),Inducer(:,5));
inds = sum(xtab > 0,2)>1;
hues0 = cellstr2mat(labels(:,1));
ihues = hues0(inds);
hu_n = numel(ihues);
step = unique(diff(ihues));
if numel(step) == 1
    fprintf('%d inducer hues at step size %d have varying chroma.\n', hu_n, step);
else
    fprintf('%d inducer hues at different step sizes have varying chroma.\n', hu_n);
end
exhues = setdiff(Inducer(:,4),ihues);
if numel(exhues) > 0
    fprintf('%d inducer hues have only 1 chroma level and are excluded.\n', numel(exhues));
end

% CHROMA
chromas = cellstr2mat(labels(:,2));
chromas = chromas(~isnan(chromas));
repN = sum(xtab,1);
inds_maxC = repN<hu_n;
maxChr      = chromas(inds_maxC);
nonmaxChr   = chromas(~inds_maxC);
maxN = sum(repN(inds_maxC));
fprintf('%d trials were ran at max chroma, which corresponds to %.2f rounds (of all hues).\n', maxN, maxN/hu_n);
if maxN ~= sum(MAIN.maxC)
    % If the identification works, this should hold.
    keyboard
end

% CHECK MULTIPLE MAX CHROMA (usually just 1):
inds = logical(MAIN.maxC);
[xtab,~,~,labels] = crosstab(Inducer(inds,4),Inducer(inds,5));
multiples = sum(xtab>0,2)>1;
if any(multiples,1)
    fprintf('WARNING: There are %d cases of multiple max chromas; they will be averaged.\n', sum(multiples));
    disp(ihues(multiples)');
end

% AGGREGATION -------------------------------------------------------------

chr_n = numel(nonmaxChr) + 1; % plus 1 for max chroma
AGG.iHue    = NaN(hu_n,chr_n);
AGG.iChr    = NaN(hu_n,chr_n);
AGG.hueM    = NaN(hu_n,chr_n);
AGG.hueS    = NaN(hu_n,chr_n);
AGG.chrM    = NaN(hu_n,chr_n);
AGG.chrS    = NaN(hu_n,chr_n);

for hu = 1:hu_n
    indsH = MAIN.inducer(:,4) == ihues(hu);
    for chr = 1:chr_n
        if chr < chr_n
            indsC = round(MAIN.inducer(:,5)) == nonmaxChr(chr);
        else
            indsC = logical(MAIN.maxC);
        end
        inds = indsH & indsC;
        if sum(inds) == 0 % When hue rounds are incomplete (level of chroma not measured for all hues)
            continue
        end
        ichroma = unique(MAIN.inducer(inds,5));
        if numel(ichroma) > 1
            ichroma = mean(MAIN.inducer(inds,5));
            fprintf('WARNING: Inducer chroma has been averaged to %.2f for hue %d\n', ichroma, ihues(hu));
            disp(unique(MAIN.inducer(inds,5))');
        end

        AGG.iHue(hu,chr) = ihues(hu);
        AGG.iChr(hu,chr) = ichroma;
        AGG.hueM(hu,chr)  = cstats(Response(inds,4), 'nanmean');
        AGG.hueS(hu,chr)  = cstats(Response(inds,4), 'ste');
        AGG.chrM(hu,chr)  = nanmean(Response(inds,5));
        AGG.chrS(hu,chr)  = ste(Response(inds,5));
        switch lower(SETTINGS.agg_stat) % REPLACE MEAN BY MEDIAN
            case {'median'}
                AGG.hueM(hu,1) = cstats(Response(inds,4), 'nanmedian');
                AGG.chrM(hu,1) = nanmedian(Response(inds,5));
        end
    end
end

% SIMULATION / MODEL ------------------------------------------------------

% ADAPTION STRENGTH:
switch lower(SETTINGS.sim_chroma)
    case 'inducer'
        adapt_chroma = AGG.iChr;
    case 'halfinducer'
        adapt_chroma = AGG.iChr/2;
end

bg_lms = XYZ2lms(stimuli.bg.XYZ, stimuli.sens, stimuli.cmf);
bg_dkl = lms2dkl(bg_lms, bg_lms); % [0 0 0] :)
varNames = {'L','u','v','hue','chr','R','G','B', 'dkl_hue'};
nans = NaN(hu_n, chr_n);
SIM.Adapt   = table(nans, nans, nans, nans, nans, nans, nans, nans, nans, 'VariableNames', varNames);
SIM.ConeCon = SIM.Adapt;
SIM.DKL     = SIM.ConeCon;
einzer = ones(hu_n,1);
for chr = 1:chr_n
    
    % ADAPT (= Inducer at adaptation strength)  ---------------------------
    aLuv_pol = [einzer*L, AGG.iHue(:,chr), adapt_chroma(:,chr)];
%    aCOL = colourconverter(aLuv_pol, 'Luv_pol', 2, WP.xyY, mon.xyY,mon.ldt,mon.oog, [], mon.rgb_max, 1, 0);
    aCOL = colourconverter(aLuv_pol, 'Luv_pol', 2,...
        WP.xyY, stimuli.sens, stimuli.cmf,...
        mon.xyY,mon.ldt,mon.oog, [], mon.rgb_max,...
        bgY, 0);
    alms = XYZ2lms(aCOL.XYZ, stimuli.sens, stimuli.cmf);
    adkl = lms2dkl(alms, bg_lms);
    SIM.Adapt.L(:,chr)    = aCOL.Luv(:,1);
    SIM.Adapt.u(:,chr)    = aCOL.Luv(:,2);
    SIM.Adapt.v(:,chr)    = aCOL.Luv(:,3);
    SIM.Adapt.hue(:,chr)  = aCOL.Luv_pol(:,1);
    SIM.Adapt.chr(:,chr)  = aCOL.Luv_pol(:,2);
    SIM.Adapt.R(:,chr)    = aCOL.rgb(:,1)./mon.rgb_max;
    SIM.Adapt.G(:,chr)    = aCOL.rgb(:,2)./mon.rgb_max;
    SIM.Adapt.B(:,chr)    = aCOL.rgb(:,3)./mon.rgb_max;
    SIM.Adapt.dkl_hue(:,chr) = rad2deg(cart2pol(adkl(:,2),adkl(:,3)));
    
    % CONE ADAPTATION -----------------------------------------------------
    induced_lms = afterimage_simulator(alms, bg_lms);
    induced_dkl = lms2dkl(alms, bg_lms);
    induced_XYZ = lms2XYZ(induced_lms, stimuli.sens, stimuli.cmf);
%    ConeAdapt = colourconverter(induced_XYZ, 'xyz',2, WP.xyY, mon.xyY,mon.ldt,mon.oog, [], mon.rgb_max, 1, 0);
    ConeAdapt = colourconverter(induced_XYZ, 'xyz', 2,...
        WP.xyY, stimuli.sens, stimuli.cmf,...
        mon.xyY,mon.ldt,mon.oog, [], mon.rgb_max,...
        bgY, 0);
    SIM.ConeCon.L(:,chr)    = ConeAdapt.Luv(:,1);
    SIM.ConeCon.u(:,chr)    = ConeAdapt.Luv(:,2);
    SIM.ConeCon.v(:,chr)    = ConeAdapt.Luv(:,3);
    SIM.ConeCon.hue(:,chr)  = ConeAdapt.Luv_pol(:,1);
    SIM.ConeCon.chr(:,chr)  = ConeAdapt.Luv_pol(:,2);
    SIM.ConeCon.R(:,chr)    = ConeAdapt.rgb(:,1)./mon.rgb_max;
    SIM.ConeCon.G(:,chr)    = ConeAdapt.rgb(:,2)./mon.rgb_max;
    SIM.ConeCon.B(:,chr)    = ConeAdapt.rgb(:,3)./mon.rgb_max;
    SIM.ConeCon.dkl_hue(:,chr) = rad2deg(cart2pol(induced_dkl(:,2),induced_dkl(:,3)));
    
    % DKL COMPLEMENTARIES -----------------------------------------------------
    induced_dkl = bg_dkl-adkl;
    induced_lms = dkl2lms(induced_dkl, bg_lms);
    induced_XYZ = lms2XYZ(induced_lms, stimuli.sens, stimuli.cmf);
%    DKL = colourconverter(induced_XYZ, 'XYZ', 2, WP.xyY, mon.xyY, mon.ldt, mon.oog, [], mon.rgb_max,1,0);
    DKL = colourconverter(induced_XYZ, 'XYZ', 2,...
        WP.xyY, stimuli.sens, stimuli.cmf,...
        mon.xyY, mon.ldt, mon.oog, [], mon.rgb_max,...
        bgY, 0);
    SIM.DKL.L(:,chr)    = DKL.Luv(:,1);
    SIM.DKL.u(:,chr)    = DKL.Luv(:,2);
    SIM.DKL.v(:,chr)    = DKL.Luv(:,3);
    SIM.DKL.hue(:,chr)  = DKL.Luv_pol(:,1);
    SIM.DKL.chr(:,chr)  = DKL.Luv_pol(:,2);
    SIM.DKL.R(:,chr)    = DKL.rgb(:,1)./mon.rgb_max;
    SIM.DKL.G(:,chr)    = DKL.rgb(:,2)./mon.rgb_max;
    SIM.DKL.B(:,chr)    = DKL.rgb(:,3)./mon.rgb_max;
    SIM.DKL.dkl_hue(:,chr) = rad2deg(cart2pol(induced_dkl(:,2),induced_dkl(:,3)));
    
end
% KILL ROUNDING MATLAB ERRORS
SIM.DKL.u = round(SIM.DKL.u,10);
SIM.DKL.v = round(SIM.DKL.v,10);
SIM.DKL.hue = round(SIM.DKL.hue,10);

% DIFFERENCES FROM MEAN & OPPONENT ----------------------------------------
einzer = ones(1,chr_n);
einzerC = ones(hu_n,1);
oppo = SIM.DKL.hue;

% DATA
Mhue = cstats(AGG.hueM', 'nanmean')';
Mhue = Mhue*einzer;
AGG.MhueM   = Mhue;
AGG.Mdev    = circularsubtractor(AGG.hueM, AGG.MhueM, 360, 'signed');
AGG.OppoDev = circularsubtractor(AGG.hueM, oppo, 360, 'signed');

Mchr = cstats(AGG.chrM, 'nanmean');
Mchr = einzerC*Mchr;
AGG.MchrM   = Mchr;
AGG.MDev_chro = circularsubtractor(AGG.chrM, AGG.MchrM, 360, 'signed');

% CONE ADAPTATION
Mhue = cstats(SIM.ConeCon.hue', 'nanmean')';
Mhue = Mhue*einzer;
SIM.ConeCon.MhueM = Mhue;
SIM.ConeCon.Mdev = circularsubtractor(SIM.ConeCon.hue, Mhue, 360, 'signed');
SIM.ConeCon.OppoDev = circularsubtractor(SIM.ConeCon.hue, oppo, 360, 'signed');

Mchr = cstats(SIM.ConeCon.chr, 'nanmean');
Mchr = einzerC*Mchr;
SIM.ConeCon.MchrM = Mchr;
SIM.ConeCon.MDev_chro = circularsubtractor(SIM.ConeCon.chr, Mchr, 360, 'signed');

% DKL
M = cstats(SIM.DKL.hue', 'nanmean')';
M = M*einzer;
SIM.DKL.Mdev = circularsubtractor(SIM.DKL.hue, M, 360, 'signed');

MC = cstats(SIM.DKL.chr, 'nanmean');
MC = einzerC*MC;
SIM.DKL.MDev_chro = circularsubtractor(SIM.DKL.chr, MC, 360, 'signed');

mx = max(abs(SIM.DKL.Mdev),[],'all');
if round(mx,10) == 0 
    fprintf('As expected: DKL Mdev is zero (max = %.2f)\n', mx);
else
    fprintf('Warning: DKL Mdev is NOT zero (max = %.2f)\n', mx);
end

% MODEL FOR ILLUSTRATION --------------------------------------------------
% This is used to create a diagram that generally illustrates the effect of
% chroma on afterimage hue.

ILLU.Adapt   = table;
ILLU.ConeCon = table;

ihues = (0:1:360)';
inducer_chroma = 10:10:80;
einzer = ones(numel(ihues),1);
switch lower(SETTINGS.sim_chroma)
    case 'inducer'
        adapt_chroma = inducer_chroma;
    case 'halfinducer'
        adapt_chroma = inducer_chroma/2;
end
chr_n = numel(adapt_chroma);
for chr = 1:chr_n

    % ADAPT (= Inducer at adaptation strength)  ---------------------------
    aLuv_pol = [einzer*L, ihues, einzer*adapt_chroma(chr)];
%    aCOL = colourconverter(aLuv_pol, 'Luv_pol', 2, WP.xyY, mon.xyY,mon.ldt,mon.oog, [], mon.rgb_max, 1, 0);
    aCOL = colourconverter(aLuv_pol, 'Luv_pol', 2,...
        WP.xyY, stimuli.sens, stimuli.cmf,...
        mon.xyY,mon.ldt,mon.oog, [], mon.rgb_max,...
        bgY, 0);
    alms = XYZ2lms(aCOL.XYZ, stimuli.sens, stimuli.cmf);
    ILLU.Adapt.L(:,chr)    = aCOL.Luv(:,1);
    ILLU.Adapt.u(:,chr)    = aCOL.Luv(:,2);
    ILLU.Adapt.v(:,chr)    = aCOL.Luv(:,3);
    ILLU.Adapt.hue(:,chr)  = ihues;
    ILLU.Adapt.chr(:,chr)  = aCOL.Luv_pol(:,2);
    ILLU.Adapt.R(:,chr)    = aCOL.rgb(:,1);
    ILLU.Adapt.G(:,chr)    = aCOL.rgb(:,2);
    ILLU.Adapt.B(:,chr)    = aCOL.rgb(:,3);
    adkl = lms2dkl(alms, bg_lms);
    ILLU.Adapt.dkl_hue(:,chr) = rad2deg(cart2pol(adkl(:,2),adkl(:,3)));

    % CONE ADAPTATION -----------------------------------------------------
    induced_lms = afterimage_simulator(alms, bg_lms);
    induced_dkl = lms2dkl(alms, bg_lms);
    induced_XYZ = lms2XYZ(induced_lms, stimuli.sens, stimuli.cmf);
%    ConeAdapt = colourconverter(induced_XYZ, 'xyz',2, WP.xyY, mon.xyY,mon.ldt,mon.oog, [], mon.rgb_max, 1, 0);
    ConeAdapt = colourconverter(induced_XYZ, 'xyz',2,...
        WP.xyY, stimuli.sens, stimuli.cmf,...
        mon.xyY,mon.ldt,mon.oog, [], mon.rgb_max,...
        bgY, 0);
    ILLU.ConeCon.L(:,chr)    = ConeAdapt.Luv(:,1);
    ILLU.ConeCon.u(:,chr)    = ConeAdapt.Luv(:,2);
    ILLU.ConeCon.v(:,chr)    = ConeAdapt.Luv(:,3);
    ILLU.ConeCon.hue(:,chr)  = ConeAdapt.Luv_pol(:,1);
    ILLU.ConeCon.chr(:,chr)  = ConeAdapt.Luv_pol(:,2);
    ILLU.ConeCon.R(:,chr)    = ConeAdapt.rgb(:,1)./mon.rgb_max;
    ILLU.ConeCon.G(:,chr)    = ConeAdapt.rgb(:,2)./mon.rgb_max;
    ILLU.ConeCon.B(:,chr)    = ConeAdapt.rgb(:,3)./mon.rgb_max;
    ILLU.ConeCon.dkl_hue(:,chr) = rad2deg(cart2pol(induced_dkl(:,2),induced_dkl(:,3)));

    ILLU.cones = simple_cone_simulator(stimuli);

end

%% simple_cone_simulator
function CONES = simple_cone_simulator(stimuli, k)
%2025.07.20

if nargin < 2
    k = 10; % S-Cone factor: produces negative chromaticities if too high
end

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
CONES = colourconverter(lms, 'lms', 2,...
    wp, stimuli.sens, stimuli.cmf,...
    stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max,...
    stimuli.bg.xyY(3), 0);

%% stimulus_polisher
function agg_stimuli = stimulus_polisher(exp_stimuli)
% 2025.10.06 [cw]

agg_stimuli             = exp_stimuli;

% OOG CHECK ---------------------------------------------------------------
oog = any(agg_stimuli.inducer.oog,2);
agg_stimuli.oog = table;
if any(oog)
    fprintf('WARNING: There are %d oogs among the %d inducers.\n', sum(oog), size(agg_stimuli.inducer,1));
    agg_stimuli.oog = agg_stimuli.inducer(oog,:);
    agg_stimuli.inducer(oog,:) = [];
    fprintf('OOGs have been deleted, the remaining induceres are %d. Corresponding data will be deleted.\n', size(agg_stimuli.inducer,1));
end

% CONE SENSITIVITIES ------------------------------------------------------
% For CIELUV: Stockman-Sharpe + CMF 1931
agg_stimuli.sens = 'ss'; 
agg_stimuli.cmf  = '1931';

% CALCULATE ACCURATE MONITOR LMS FOR DKL CONVERSION -----------------------
agg_stimuli.mon.lms = XYZ2lms(xyY2XYZ(agg_stimuli.mon.xyY), agg_stimuli.sens, agg_stimuli.cmf);

% ROUNDING ERRORS --------------------------------------------------------- 
% This version only works up from MATLAB 2023:
%agg_stimuli.bg      = round(agg_stimuli.bg,12);
%agg_stimuli.inducer = round(agg_stimuli.inducer,12);
% For backcompatibility, the loop below has been implemented:
varnames = agg_stimuli.bg.Properties.VariableNames;
for vr = 1:numel(varnames) % The loop is necessary for MATLAB versions below 2023a
    agg_stimuli.bg.(varnames{vr})      = round(agg_stimuli.bg.(varnames{vr}),12);
    agg_stimuli.inducer.(varnames{vr}) = round(agg_stimuli.inducer.(varnames{vr}),12);
end

% RELATIVE RGB (to display on conventional monitor):
agg_stimuli.inducer.rgb = agg_stimuli.inducer.rgb/agg_stimuli.mon.rgb_max;
agg_stimuli.bg.rgb      = agg_stimuli.bg.rgb/agg_stimuli.mon.rgb_max; %

% MONITOR GAMUT -----------------------------------------------------------
Y = agg_stimuli.bg.xyY(3);
gam_XYZ = iso_gamuter(Y, agg_stimuli.mon.xyY);
%agg_stimuli.gamut = colourconverter(gam_XYZ,'XYZ',2,'mon', agg_stimuli.mon.xyY, agg_stimuli.mon.ldt, agg_stimuli.mon.oog, [], agg_stimuli.mon.rgb_max, 1, 0);
agg_stimuli.gamut = colourconverter(gam_XYZ,'XYZ',2,...
    'mon', agg_stimuli.sens, agg_stimuli.cmf,...
    agg_stimuli.mon.xyY, agg_stimuli.mon.ldt, agg_stimuli.mon.oog, [], agg_stimuli.mon.rgb_max,...
    agg_stimuli.bg.xyY(3), 0);
agg_stimuli.gamut = agg_stimuli.gamut([1 2 4 3 1],:);

% REMOVE IRRELEVANT MEASUREMENTS ------------------------------------------
[tab,~,~,labels] = crosstab(agg_stimuli.inducer.Luv_pol(:,1),agg_stimuli.inducer.Luv_pol(:,2));
incl = sum(tab,2) >1;
hues0 = cellstr2mat(labels(:,1));
ihues = hues0(incl);
hu_n = numel(ihues);
step = unique(diff(ihues));
if numel(step) == 1
    fprintf('%d inducer hues at step size %d have varying chroma.\n', hu_n, step);
else
    fprintf('%d inducer hues at different step sizes have varying chroma.\n', hu_n);
end
[exhues, delete] = setdiff(agg_stimuli.inducer.Luv_pol(:,1), ihues);
agg_stimuli.inducer(delete,:) = [];
if numel(exhues) > 0
    fprintf('%d inducer hues have only 1 chroma level and are excluded.\n', numel(exhues));
end

% AVERAGE MULTIPLE MAX CHROMAS (only last pp) -----------------------------
inds = logical(agg_stimuli.inducer.maxC);
[xtab,~,~,labels] = crosstab(agg_stimuli.inducer.Luv_pol(inds,1),agg_stimuli.inducer.Luv_pol(inds,2));
multiples = sum(xtab>0,2)>1;
if any(multiples,1)
    hues0 = cellstr2mat(labels(:,1));
    multiN = sum(multiples);
    fprintf('WARNING: There are %d cases of multiple maximum chromas per hue; they will be averaged.\n', multiN);
    multi_hues = hues0(multiples,1);
    for multi = 1:multiN
        fprintf('For %d, chromas: ', multi_hues(multi));
        indsMu = agg_stimuli.inducer.Luv_pol(:,1) == multi_hues(multi);
        indsMx = agg_stimuli.inducer.maxC;
        indsM = indsMu & indsMx;
        fprintf('%d ', agg_stimuli.inducer.Luv_pol(indsM,2)');
        Mchr = mean(agg_stimuli.inducer.Luv_pol(indsM,2));
        agg_stimuli.inducer.Luv_pol(indsM,2) = Mchr;
        fprintf('have been set to %.f.\n', Mchr);
    end
    % Delete the doubles produced by avaraging:
     [~, inds] = (unique(agg_stimuli.inducer.Luv_pol,'rows'));
    agg_stimuli.inducer = agg_stimuli.inducer(inds,:);
    fprintf('After averaging, the remaining induceres are %d.\n', size(agg_stimuli.inducer,1));
end

%% ***************************** SUBFUNCTIONS *****************************

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

%% cellstr2mat
function A = cellstr2mat(cellstring)
%2022.09.16

A = cellfun(@(x) str2double(x),cellstring);

%% circularsubtractor
function Difference = circularsubtractor(Minuend, Subtrahend, ONECYCLE, dtype)
% 2012.08.23-24 [cw]
        
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
                M = [0.15282  0.54383 -0.02795;...
                    -0.15254  0.45524  0.03355;...
                    -0.00045  0.00145  0.95449];
                
            case 'judd'
                % Smith & Pokorny cone fundamentals
                % V. C. Smith & J. Pokorny (1975), Vision Res. 15, 161-172.
                %        X          Y       Z [cw]
                M = [ 0.15514  0.54312  -0.03286    % L alias R
                     -0.15514  0.45684   0.03286    % M alias G
                      0.0      0.0       0.01608];  % S alias B
        end
    case 'smj2'
        % A. Stockman, D. I. A. MacLeod, and N. E. Johnson, ‘‘Spectral
        % sensitivities of the human cones,’’ J. Opt. Soc. Am. A 10,
        % 2491–2521 (1993).
        switch lower(cmf)
            case '1931'
                M = [ 0.18772  0.60445 -0.02517;...
                     -0.14014  0.43056  0.03773;...
                      0.02017 -0.04189  1.08472];
        end
    case 'smj10'
        % A. Stockman, D. I. A. MacLeod, and N. E. Johnson, ‘‘Spectral
        % sensitivities of the human cones,’’ J. Opt. Soc. Am. A 10,
        % 2491–2521 (1993).
        switch lower(cmf)
            case '1931'
                M = [ 0.14460  0.62421 -0.00429;...
                     -0.14506  0.42265  0.05084;...
                      0.03105 -0.06416  1.10923];
        end
    case {'stockmansharpe', 'ss2000', 'ss'}
        % A. Stockman and L. T. Sharpe, ‘‘Spectral sensitivities of the
        % middle- and long-wavelength sensitive cones derived from
        % measurements in observers of known genotype,’’ Vision Res.
        % 40, 1711–1737 (2000).
        switch lower(cmf)
            case '1931'
                % Taken from Golz & MacLeod (2003).
                M = [ 0.17156  0.52901 -0.02199;...
                     -0.15955  0.48553  0.04298;...
                      0.01916 -0.03989  1.03993];
        end
    case {'hunt-pointer-estevez', 'hpe', 'rlab'}
        % https://en.wikipedia.org/wiki/LMS_color_space
        switch lower(cmf)
            case '1931'
                M = [0.38971,   0.68898,   -0.07868;...
                    -0.22981,   1.18340,    0.04641;... 
                     0.00000,   0.00000,    1.00000];
        end
end

if isempty(M)
    error('This transformation matrix is not implemented');
end

%% cstats
function varargout = cstats(inM, stype, onecycle, remap)
% 2015.10.30 [cw].

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
%2012jul02 * [cw]
%2015feb04 level4 [cw]

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

%% dkl2lms
function [lms, M] = dkl2lms(dkl, lms0, monlms)
% Inverts lms2dkl by inverting the matrix.
% 2025.04.23

L0 = lms0(1);
M0 = lms0(2);
S0 = lms0(3);

% MATRIX THAT SCALES RELATIVE TO WP (L0, M0, S0)
M = [...
    1/(L0+M0),  1/(L0+M0),      0;...
    1/L0,       -1/M0,      0;...
    -1/(L0+M0), -1/(L0+M0), 1/S0];

if nargin > 2
    [maxL, maxS] = map2mongamut(lms0, monlms);
    M(2,:) = M(2,:)/maxL;
    M(3,:) = M(3,:)/maxS;
end

dlms = M \ dkl';
dlms = dlms';

% INCREMENTS TO LMS:
lms = dlms+ones(size(dlms,1),1)*lms0;
lms = round(lms,12);

%% iso_gamuter
function [XYZ, only, full] = iso_gamuter(lum, mon_xyY)
%2023.01.14 * [cw]

if nargin < 2
    mon_xyY = srgb;
    if nargin == 0
        lum = 50;
    end
end

mon_XYZ = xyY2XYZ(mon_xyY);

R = mon_xyY(1,3);
G = mon_xyY(2,3);
B = mon_xyY(3,3);

if lum > R+G+B
    error('lum higher than white');
end

kR = helper(lum, R, G, B);
kG = helper(lum, G, B, R);
kG = circshift(kG,1,2);
kB = helper(lum, B, R, G);
kB = circshift(kB,2,2);

full = [kR; kG; kB];
only = full(any(~isnan(full),2),:);
only = uniquetol(only, 10^-12, 'ByRows', true);
% k = convhull(only(:,1),only(:,2),only(:,3));
% only = only(k,:);

xR = bsxfun(@times,mon_XYZ(1,:),only(:,1));
xG = bsxfun(@times,mon_XYZ(2,:),only(:,2));
xB = bsxfun(@times,mon_XYZ(3,:),only(:,3));

XYZ = xR+xG+xB;

%% helper
function kR = helper(lum, R, G, B)

kR = NaN(1,3);
kRg = NaN(1,3);
kRb = NaN(1,3);
kRgb = NaN(1,3); 
kRbg = NaN(1,3);

if lum < R
    kR = [lum/R 0 0];
else
    if lum < R+G
        kRg = [1 (lum-R)/G 0]; 
    else
        kRgb = [1 1 (lum-R-G)/B];
    end
    if lum < R+B
        kRb = [1 0 (lum-R)/B]; 
    else
        kRbg = [1 (lum-R-B)/G 1];
    end
end
kR = [kR; kRg; kRb; kRgb; kRbg];

%% linemaker
function ln = linemaker(width, style)
% 2010dec17 as a subfunction [cw]
% 2011apr06 as a function [cw]
if nargin < 2
    style = '_';
    if nargin < 1
        width = 100;
    end
end
ln(1:width) = style;

%% lms2dkl
function [dkl, M] = lms2dkl(lms, lms0, monlms)
% Converts cone excitations (lms) to DKL coordinates (dkl) realtive to the
% adapting white-point (lms0). Optionally, DKL space may be rescaled to fit
% into the monitor gamut (monlms).
% INPUT:
% lms    = cone excitations (in columns) of the stimuli (in rows)
% lms0   = cone excitations of the adapting stimulus (e.g., the background)
% monlms = L, M, S (in columns) for the monitor primaries R, G, B (in rows)
% OUTPUT:
% dkl   = DKL axis (in columns) for each lms stimulus (in rows).
% M     = transformation matrix that converts from lms to DKL.
% The transformation matrix M is equivalent to the one from Brainard
% (1996), except for the scaling. Input monlms rescales the axis so that
% the isoluminant hue circle fits into the monitor gamut; this is the same
% scaling as in DKLDemo of Psychtoolbox. 
% 2025.04.11
 
% INCREMENTS:
dlms = lms-ones(size(lms,1),1)*lms0;
 
L0 = lms0(1);
M0 = lms0(2);
S0 = lms0(3);
 
% TRANSFORMATION MATRIX THAT SCALES RELATIVE TO WP (L0, M0, S0)
M = [...
    1/(L0+M0),  1/(L0+M0),      0;...
    1/L0,       -1/M0,      0;...
    -1/(L0+M0), -1/(L0+M0), 1/S0];


% RESCALING TO MONITOR GAMUT
if nargin > 2    
    [maxL, maxS] = map2mongamut(M, lms0, monlms);
    M(2,:) = M(2,:)/maxL;
    M(3,:) = M(3,:)/maxS;
end

% APPLY TRANSFORMATION MATRIX
dkl = dlms * M';

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

%% map2mongamut
function [maxLM, maxS] = map2mongamut(M, lms0, monlms)
% Produces scaling factors for the L-M and S-(L+M) axis. This function is
% based on David Brainard's function MaximizeGamutContrast in Psychtoolbox.
% INPUT:
% M = LMS to DKL conversion matrix, as implemented in function lms2dkl.
% lms0 = LMS signal (cone excitations) of the adapting white-point (i.e.,
% background). 
% monlms = LMS of the monitor primaries (at maximum).
% OUTPUT:
% maxLM, maxS = Maximum within gamut expressed relative to the original
% dkl unit; dividing unit length by these maxima makes units the length of
% those maxima. 
% 2025.07.11 based on PTB MaximizeGamutContrast. [cw]

% DEFINE THE POLES IN DKL AT UNIT LENGTH ----------------------------------
pols = [...
    0  1  0;... L-M pole
    0 -1  0;... M-L pole
    0  0  1;... S-(L+M) pole
    0  0 -1;... (L+M)-S pole
    ];

% CONVERT DKL TO CONE INCREMENTS (dlms) -----------------------------------
% Cone increments = lms signal - lms0;
pole_dlms = pols / M';

% CONVERT LMS TO RGB ------------------------------------------------------
% This can be done directly with the cone-increments: 
rgb0 = lms0/(monlms(1:3,1:3));
pole_drgb = pole_dlms/(monlms(1:3,1:3));

% The below code, going through lms signal first, is equivalent to direct
% conversion of dlms:
% pole_lms = pole_dlms+ones(size(pole_dlms,1),1)*lms0;
% pole_lms = round(pole_lms,12);
% pole_rgb= pole_lms/(monlms(1:3,1:3));
% pole_drgb = pole_rgb-rgb0

% CALCULATE THE RATIO BETWEEN DKL UNIT LENGTH AND MAX RGB -----------------
% RGB gamut is 0 (minimum) and 1 (maximum); subtracting the rgb0 (adapting
% wp) allows determing the length away from rgb0. The result of this
% calculation indicates for each primary how many time the original DKL
% unit fits within the gamut of that primary.

% L-M Mechanism:
maxL = (1-rgb0)./(pole_drgb(1,:)); % L-M direction
maxM = (0-rgb0)./(pole_drgb(2,:)); % M-L direction

% S-(L+M) Mechanism:
maxS1 = (1-rgb0)./(pole_drgb(3,:)); % +S direction
maxS2 = (0-rgb0)./(pole_drgb(4,:)); % -S direction

% IDENTIFY SCALING FACTORS ------------------------------------------------
% The primary that defines the limit in a given direction is the smallest
% ratio because that primatry would be the first to fall out of gamut when
% increasing DKL saturation. As both direction of the axes are scaled by
% the same factor, the scaling factor is the minimum length (absolute min)
% of both directions together. 

maxLM = min(abs([maxL,maxM]));
maxS = min(abs([maxS1,maxS2]));

%% ste
function y = ste(varargin)
% Computes Standard error of mean.
% 2009.08.18 [cw]

y = sqrt(nanvar(varargin{:}))/sqrt(length(varargin{1}));

%% xyY2XYZ
function XYZ = xyY2XYZ(xyY, dim)
%input: x, y, Y values in a row vector or matrix )
% 2012jan20 added dim to enable image input [cw]

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

%% XYZ2Luv
function LUV = XYZ2Luv(XYZ, XYZ_bg, dim)

%___OBJECTIVE______________________________________________________________
%calculates the CIE Luv coordinates from CIE31 XYZ coordinates. Input
%are the XYZ coordinates of the stimulus and those of the background
%(usually white/gray point of the monitor)
%___GENEALOGY______________________________________________________________
% 2011mar08 put exact fraction (instead of rounded decimals), enable
% matrix input, relabeled XYZ2Luv [cw] 
% 2012jan20 added option dim for image conversion; corrected bug in
% if-statement by using indices [cw]
% 2015jun05 corrected NaN for Y=0 [cw]


if nargin < 3
    dim = 2;
end

% Whitepoint:
Xn = XYZ_bg(1);
Yn = XYZ_bg(2);
Zn = XYZ_bg(3);


if dim == 2
    X = XYZ(:,1);
    Y = XYZ(:,2);
    Z = XYZ(:,3);
elseif dim == 3
    X = XYZ(:,:,1);
    Y = XYZ(:,:,2);
    Z = XYZ(:,:,3);
end


uprime = (4 * X)./(X + 15 * Y + 3 * Z);
vprime = (9 * Y)./(X + 15 * Y + 3 * Z);
uprimen = (4 * Xn)./(Xn + 15 * Yn + 3 * Zn);
vprimen = (9 * Yn)./(Xn + 15 * Yn + 3 * Zn); 

L = 116 * (Y/Yn).^(1/3) - 16;
L2 = ((29/3)^3) * (Y/Yn);
inds = Y/Yn > (6/29)^3;
L(~inds) = L2(~inds);

u = (13*L) .* (uprime - uprimen);
v = (13*L) .* (vprime - vprimen);

inds = Y == 0;
L(inds) = 0; u(inds) = 0; v(inds) = 0;

if dim == 2
    LUV = [L u v];
elseif dim == 3
    LUV = cat(3, L, u, v);    
end