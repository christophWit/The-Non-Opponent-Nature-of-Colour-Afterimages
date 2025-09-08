function AGG = mod_exp2_aggregator(EXP2ab, MunChips, HeringRYGB)
%2024.07.02 * as subfunction [cw]
%2025.07.16 polished for publication [cw]

header('AGGREGATION - EXPERIMENT 2', 2);

%% ******************************* SETTINGS *******************************
SETTINGS.trgt_colspace = EXP2ab.stimuli{1}.src_space; % Same as source
SETTINGS.agg_stat      = 'mean';
SETTINGS.huefrq_res    = 5; % 5-degree steps produce 72 bins
SETTINGS.smooth_n      = 9; % Number of adjacent hue frequency bins to smooth over
SETTINGS.sim_n         = 1000; % Number of repetitions in simulation, set 1000 for final

%% ***************************** AGGREGATION ******************************
header('DATA AGGREGATION',2);

%% PARTICIPANTS
AGG.pps = EXP2ab.pps;

%% APPARATUS & STIMULI
for pp = 1:size(AGG.pps,1)
    AGG.stimuli{pp} = stimulus_polisher(EXP2ab.stimuli{pp});
end

%% AGGREGATION WITHIN INDIVIDUALS

header('AGGREGATION', 3);
%ihues = (0:5:355)'; % Assumed inducer hues
pp_n = size(AGG.pps,1);
AGG.pps.repN = NaN(pp_n,2);
AGG.IndiM = table;
AGG.hueF  = table;
for pp = 1:pp_n
    header(EXP2ab.pps.pp_lbl{pp}, 4);
    indsP = EXP2ab.main.pp == pp;

    % CHECK NUMBER OF REPETITIONS ---------------------------------------------
    tab = tabulate(EXP2ab.main.inducer(indsP,4));
    AGG.pps.repN(pp,1) = mode(tab(:,2));
    AGG.pps.repN(pp,2) = max(tab(:,2));
    fprintf('REPETITIONS: mode: %d, max: %d\n', AGG.pps.repN(pp,1), AGG.pps.repN(pp,2));

    % AGGREGATE -----------------------------------------------------------
    ihues = AGG.stimuli{pp}.inducer.(AGG.stimuli{pp}.src_polar)(:,1);
    [agg1pp, hueF, AGG.Indi{pp}]  = aggregation_helper(EXP2ab.main(indsP,:), AGG.pps.repN(pp,2), AGG.stimuli{pp}.bg, ihues, SETTINGS);
    AGG.IndiM = tablecol_adder(AGG.IndiM, agg1pp);
    AGG.hueF  = tablecol_adder(AGG.hueF,hueF);
    
end

%% MODELLING
header('MODELLING', 3);
pp_n = size(AGG.pps,1);
AGG.model.Hue    = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.model.Chroma = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.model.Lum    = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.model.hueF   = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.prederr.hueM = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.prederr.hueS = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.prederr.chrM = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.prederr.chrS = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.prederr.ihue = cell(1,pp_n);
AGG.prederr.ichr = cell(1,pp_n);
noisy_azi        = cell(1,pp_n);
IndiOppoDev      = cell(1,pp_n);
for pp = 1:pp_n
    header(EXP2ab.pps.pp_lbl{pp}, 4);

    AGG.pps.chrMM(pp,1) = mean(AGG.IndiM.chrM(:,pp));
    AGG.pps.chrSD(pp,1) = nanmean(nanstd(AGG.Indi{pp}(:,:,5),0,2));
    AGG.pps.hueSD(pp,1) = nanmean(cstats(AGG.Indi{pp}(:,:,4)', 'std', 360)');
    [OppoHue, OppoChroma, OppoLum, dkl_k, Munsell] = complementary_calculator(AGG.stimuli{pp}, SETTINGS.trgt_colspace, MunChips, HeringRYGB, AGG.pps.chrMM(pp,1));
    AGG.pps.dkl_k(pp,1) = dkl_k; % Strength of subtractive adaptation in the cone-opponent DKL model.
    AGG.model.Hue    = tablecol_adder(AGG.model.Hue, OppoHue);
    AGG.model.Chroma = tablecol_adder(AGG.model.Chroma, OppoChroma);
    AGG.model.Lum    = tablecol_adder(AGG.model.Lum, OppoLum);

    [AGG.prederr.ihue{pp}, AGG.prederr.ichr{pp}, pe_hueM, pe_chrM, pe_hueSEM, pe_chrSEM] =...
        predictionerrorcalculator(AGG.Indi{pp}, OppoHue, OppoChroma);
    IndiOppoDev{pp,1}  = AGG.prederr.ihue{pp}.dkl;
    AGG.prederr.hueM  = tablecol_adder(AGG.prederr.hueM, pe_hueM);
    AGG.prederr.hueS  = tablecol_adder(AGG.prederr.hueS, pe_hueSEM);
    AGG.prederr.chrM  = tablecol_adder(AGG.prederr.chrM, pe_chrM);
    AGG.prederr.chrS  = tablecol_adder(AGG.prederr.chrS, pe_chrSEM);

    % SIMULATE HUE FREQUENCIES WITH NOISE ---------------------------------
    noisy_azi{pp,1} = noise_maker([AGG.model.Hue.ConeCon, AGG.model.Chroma.ConeCon], AGG.pps.hueSD(pp,1), AGG.pps.chrSD(pp,1), SETTINGS.sim_n);
    hueF = hue_frequency_calculator(noisy_azi{pp,1}, SETTINGS.huefrq_res, SETTINGS.smooth_n);
    hueF.Properties.RowNames = AGG.hueF.Properties.RowNames;
    AGG.model.hueF  = tablecol_adder(AGG.model.hueF,hueF);

    % FOR FIGURES ---------------------------------------------------------
    AGG.dir(pp).cones = simple_cone_simulator(AGG.stimuli{pp});
    AGG.dir(pp).HeringRYGB = HeringRYGB;

end

% DEVIATIONS FROM OPPONENCY -----------------------------------------------
AGG.IndiOppoDev               = table('RowNames', AGG.IndiM.Properties.RowNames);
AGG.IndiOppoDev.adj           = AGG.prederr.hueM.dkl;
AGG.IndiOppoDev.adj_sem       = AGG.prederr.hueS.dkl;
AGG.IndiOppoDev.adj_smooth    = smoothdata(AGG.prederr.hueM.dkl,'movmean',SETTINGS.smooth_n);
AGG.IndiOppoDev.sim           = circularsubtractor(AGG.model.Hue.ConeCon, AGG.model.Hue.dkl, 360,'signed');
AGG.IndiOppoDev.diff          = circularsubtractor(AGG.prederr.hueM.ConeCon, AGG.prederr.hueM.dkl, 360,'signed');

%% AGGREGATE ACROSS PARTICIPANTS
alldata   = cat(2,AGG.Indi{:}); 
AGG.chromaM = mean(alldata(:,:,5), 1);
AGG.agg = table('RowNames', AGG.IndiM.Properties.RowNames);
weights = AGG.hueF.N; % Weights for weighted average to account for different repN.

% AVERAGES ----------------------------------------------------------------
AGG.agg.hue           = [cstats(alldata(:,:,4)', 'mean')', cstats(alldata(:,:,4)', 'ste')'];
AGG.agg.chroma        = [mean(alldata(:,:,5),2), ste(alldata(:,:,5), 0, 2)];
AGG.agg.chroma_smooth = smoothdata(AGG.agg.chroma(:,1),'movmean',SETTINGS.smooth_n);

% DEVIATIONS FROM OPPONENCY -----------------------------------------------
alloppodev              = cat(2,IndiOppoDev{:}); 
AGG.agg.oppdev          = [mean(alloppodev,2), ste(alloppodev, 0, 2)];
AGG.agg.oppdev_semM     = sum(weights.*AGG.IndiOppoDev.adj_sem,2)./sum(weights,2);
AGG.agg.oppdev_smooth   = smoothdata(AGG.agg.oppdev(:,1),'movmean',SETTINGS.smooth_n);

% HUE FREQUENCIES ---------------------------------------------------------
hueF        = hue_frequency_calculator(alldata(:,:,4), SETTINGS.huefrq_res, SETTINGS.smooth_n);
AGG.agg     = [AGG.agg, hueF];

% SIMULATIONS -------------------------------------------------------------
% WEIGHTED AVERAGE TO ACCOUNT FOR DIFFERENT repN:

% AVERAGES
AGG.agg.sim_hueM = sum(weights.*AGG.model.Hue.ConeCon,2)./sum(weights,2);
AGG.agg.sim_chrM = sum(weights.*AGG.model.Chroma.ConeCon,2)./sum(weights,2);

% DEVIATIONS FROM OPPONENCY
AGG.agg.sim_oppdev = sum(weights.*AGG.IndiOppoDev.sim,2)./sum(weights,2);
AGG.agg.sim_diff   = AGG.agg.oppdev(:,1)-AGG.agg.sim_oppdev;

% HUE FREQUENCIES
AGG.agg.sim_frq   = sum(weights.*AGG.model.hueF.frq,2)./sum(weights,2);
AGG.agg.sim_rF    = AGG.agg.sim_frq./SETTINGS.sim_n;

%% ****************************** SUBMODULES ******************************

%% aggregation_helper
function [agg, hueF, indi] = aggregation_helper(MAIN, repNmax, BG, ihues, SETTINGS)
% 2025.07.16 [cw]

switch lower(SETTINGS.trgt_colspace)
    case 'luv'
        bg = BG.Luv;
        Inducer     = MAIN.inducer;
        Response    = MAIN.adj_Luv; 
    case 'dkl'
        bg          = BG.dkl;
        Inducer     = MAIN.ind_dkl;
        Response    = MAIN.adj_dkl; 
end
Inducer(:,4)     = round(Inducer(:,4));

% DOUBLE-CHECK THAT INDUCERS ARE AS EXPECTED ------------------------------
% LIGHTNESS
L = unique(Inducer(:,1));
if numel(L) > 1 | L ~= bg(:,1)
    disp(L);
    disp(bg(:,1));
    error('Inducer lightness is not as expected.');
else
    fprintf('Inducer lightness is: %.1f\n', L);
end

% HUE
hu_n = numel(ihues);
[tmp] = intersect(Inducer(:,4), ihues);
if numel(tmp) ~= hu_n
    error('Number of inducer hues does not match');
else
    fprintf('There are %d inducer hues.\n', hu_n);
end

% CHROMA
C = unique(Inducer(:,5));
if numel(C) > 1
    warning('There are several levels of inducer chroma.');
    disp(C);
else
    fprintf('Inducer chroma is: %.1f\n', C);
end

% MAIN: ACTUAL AGGREGATION ------------------------------------------------
indi   = NaN(hu_n,repNmax,5);
repN   = NaN(hu_n,1);
hueM   = NaN(hu_n,1);
hueSEM = NaN(hu_n,1);
chrM   = NaN(hu_n,1);
chrSEM = NaN(hu_n,1);
for hu = 1:hu_n
    indsH = Inducer(:,4) == ihues(hu);
    repN(hu,1) = sum(indsH);

    % INDIVIDUAL DATA
    indi(hu,1:sum(indsH),1) = Response(indsH,1);
    indi(hu,1:sum(indsH),2) = Response(indsH,2);
    indi(hu,1:sum(indsH),3) = Response(indsH,3);
    indi(hu,1:sum(indsH),4) = Response(indsH,4);
    indi(hu,1:sum(indsH),5) = Response(indsH,5);

    % AGGREGATED DATA
    hueM(hu,1)   = cstats(Response(indsH,4), 'nanmean');
    hueSEM(hu,1) = cstats(Response(indsH,4), 'ste');
    chrM(hu,1)   = nanmean(Response(indsH,5));
    chrSEM(hu,1) = ste(Response(indsH,5));
    switch lower(SETTINGS.agg_stat) % REPLACE MEAN BY MEDIAN
        case {'median'}
            hueM(hu,1) = cstats(Response(indsH,4), 'nanmedian');
            chrM(hu,1) = nanmedian(Response(indsH,5));
    end
end
chrM_smooth = smoothdata(chrM,'movmean',SETTINGS.smooth_n);
agg = table(repN, hueM, hueSEM, chrM, chrSEM, chrM_smooth);
agg.Properties.RowNames = cellstr(num2str(ihues));

% HUE FREQUENCIES FOR HUE HISTOGRAM ---------------------------------------

% Discard incomplete rounds for computing hue frequencies
incomplete_rounds = any(isnan(indi(:,:,4)),1); % Any Column containing NaNs
if sum(incomplete_rounds)> 0
    n = sum(~isnan(indi(:,incomplete_rounds,4)));
    fprintf('FULLROUNDS: Deleted round %d with %d data points for hue frequencies and aggregation across pps\n', find(incomplete_rounds), n);
    indi = indi(:, ~incomplete_rounds,:);
end
f_n = size(indi,2);
fprintf('CHECK: Full rounds = %d (Should be the same as the mode of repetitions).\n', f_n);

hueF = hue_frequency_calculator(indi(:,:,4), SETTINGS.huefrq_res, SETTINGS.smooth_n);
hueF.Properties.RowNames = cellstr(num2str(ihues));

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
frq_smooth  = smoothdata(frq,'movmean',smoother);
rF_smooth   = smoothdata(rF,'movmean',smoother);

hueF = table(N, frq, rF, frq_smooth, rF_smooth);

%% predictionerrorcalculator
function [pe_hue, pe_chroma, pe_hueM, pe_chromaM, pe_hueSEM, pe_chromaSEM] = predictionerrorcalculator(INDI, OppoHue, OppoChr)
% 2025.07.19 [cw]

match_hue = INDI(:,:,4); % because already converted to space
match_chr = INDI(:,:,5); % because already converted to space
zmatch    = zscore(match_chr);

pe_hue      = table; 
pe_hueM     = table;
pe_hueSEM   = table;

pe_chroma     = table;
pe_chromaM    = table;
pe_chromaSEM  = table;

varnames = OppoHue.Properties.VariableNames;
for vn = 1:numel(varnames)
    vrnm = varnames{vn};

    % HUE -----------------------------------------------------------------
    oppo_hue = OppoHue.(vrnm);
    pe_hue.(vrnm) = circularsubtractor(match_hue, oppo_hue*ones(1,size(match_hue,2)), 360, 'signed');
    
    % AGGREGATE:
    Mhue = mean(pe_hue.(vrnm), 2);
    semH = ste(pe_hue.(vrnm), 0, 2);
    pe_hueM.(vrnm)    = Mhue;
    pe_hueSEM.(vrnm)  = semH;
    
    % CHROMA --------------------------------------------------------------
    oppo_chr = OppoChr.(vrnm);
    if size(unique(oppo_chr, 'rows'),1)==1
        zoppo = zeros(size(oppo_chr));
    else
        zoppo = zscore(oppo_chr);
    end    
    pe_chroma.(vrnm) = zmatch - zoppo*ones(1,size(match_chr,2));

    % AGGREGATE:
    Mchr   = mean(pe_chroma.(vrnm), 2);
    semC = ste(pe_chroma.(vrnm), 0, 2);
    pe_chromaM.(vrnm)     = Mchr;
    pe_chromaSEM.(vrnm)   = semC;

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
XYZ = lms2XYZ(lms, stimuli.sens, stimuli.cmf);
CONES = colourconverter(XYZ, 'XYZ', 2,...
    wp, stimuli.sens, stimuli.cmf,...
    stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max,...
    stimuli.bg.xyY(3), 0);

%[CONES.dkl, CONES.dkl_pol, CONES.lms] = XYZ2dkl(XYZ, bg_lms, stimuli.sens, stimuli.cmf, stimuli.mon.lms);

%% stimulus_polisher
function agg_stimuli = stimulus_polisher(exp_stimuli)
% 2025.07.16 [cw]

agg_stimuli             = exp_stimuli;

switch lower(agg_stimuli.src_space)
    case 'luv'
        agg_stimuli.src_polar = 'Luv_pol';
        agg_stimuli.sens = 'ss'; 
        agg_stimuli.cmf  = '1931';
    case 'dkl'
        agg_stimuli.src_polar = 'dkl_pol';
        agg_stimuli.sens = 'sp'; 
        agg_stimuli.cmf  = 'judd';
end

% CALCULATE ACCURATE MONITOR LMS FOR DKL CONVERSION -----------------------
agg_stimuli.mon.lms = XYZ2lms(xyY2XYZ(agg_stimuli.mon.xyY), agg_stimuli.sens, agg_stimuli.cmf);

% ROUNDING ERRORS --------------------------------------------------------- 
varnames = agg_stimuli.bg.Properties.VariableNames;
for vr = 1:numel(varnames) % The loop is necessary for MATLAB versions below 2023a
    agg_stimuli.bg.(varnames{vr})      = round(agg_stimuli.bg.(varnames{vr}),12);
    agg_stimuli.inducer.(varnames{vr}) = round(agg_stimuli.inducer.(varnames{vr}),12);
end

% RELATIVE RGB (to display on conventional monitor):
agg_stimuli.inducer.rgb = agg_stimuli.inducer.rgb/agg_stimuli.mon.rgb_max;
agg_stimuli.bg.rgb      = agg_stimuli.bg.rgb/agg_stimuli.mon.rgb_max; %

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
% Computes Standard error of mean.
% 2009.08.18 [cw]

y = sqrt(nanvar(varargin{:}))/sqrt(length(varargin{1}));

%% tablercol_adder
function T0 = tablecol_adder(T0,Tnew)
% 2025.07.16

if isempty(T0)
    T0 = Tnew;
    return;
end

vrnms = T0.Properties.VariableNames;
for vr = 1:numel(vrnms)
    T0.(vrnms{vr}) = [T0.(vrnms{vr}), Tnew.(vrnms{vr})];
end

%% XYZ2dkl
function [dkl, dkl_pol, lms] = XYZ2dkl(XYZ, bg_lms, sens, cmf, mon_lms)
% 2025.07.31 [cw]

% LMS ---------------------------------------------------------------------
lms = XYZ2lms(XYZ, sens, cmf);

% DKL ---------------------------------------------------------------------
dkl = lms2dkl(lms, bg_lms, mon_lms);

% POLAR DKL ---------------------------------------------------------------
[dkl_azi, dkl_rad] = cart2pol(dkl(:,2),dkl(:,3));
dkl_azi = rad2deg(dkl_azi);
dkl_pol = [dkl_azi(:), dkl_rad(:)];

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
        % A. Stockman, D. I. A. MacLeod, and N. E. Johnson, ‘‘Spectrallms2X
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