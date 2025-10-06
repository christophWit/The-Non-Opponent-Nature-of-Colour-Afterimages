function [RES, fh] = mod_exp3_analysor(AGG, which_one)
%2025.07.16 adapted + polished for publication [cw]

%% ******************************* SETTINGS *******************************

AGG.colspace = AGG.stimuli{1}.src_space; % Represent results in the colour space in which data was collected.

% INITIALISE OUTPUT STRUCTURE WITH RESULTS
RES = struct; 

%% DESIGN (cosmetics)
DESIGN.fontsize = 7;

%% CORRELATION TABLES

header('Correlation Tables', 3);

pp_n = size(AGG.pps,1);
for pp = 1:pp_n
    sim = reshape(AGG.model{pp}.ConeCon.Mdev,[],1);
    data = reshape(AGG.data{pp}.Mdev,[],1);
    RES.ctab_Mdev(pp,:) = correlator(data, sim);

    sim = reshape(AGG.model{pp}.ConeCon.OppoDev,[],1);
    data = reshape(AGG.data{pp}.OppoDev,[],1);
    RES.ctab_oppodev(pp,:) = correlator(data, sim);

    data = reshape(AGG.data{pp}.chrM,[],1);
    cone = reshape(AGG.model{pp}.ConeCon.chr,[],1);
    dkl = reshape(AGG.model{pp}.DKL.chr,[],1);
    RES.ctab_chroma_cone(pp,:) = correlator(data, cone);
    RES.ctab_chroma_dkl(pp,:)  = correlator(data, dkl);

    hu_n = size(AGG.data{pp},1);
    [p1(pp,1), z1(pp,1)] = cmpcorr_ztest(RES.ctab_chroma_cone.r(pp), RES.ctab_chroma_dkl.r(pp), hu_n, hu_n);

    data = AGG.data{pp}.MDev_chro(:);
    cone = AGG.model{pp}.ConeCon.MDev_chro(:);
    dkl = AGG.model{pp}.DKL.MDev_chro(:);
    RES.ctab_Mdev_chroma_cone(pp,:) = correlator(data, cone);
    RES.ctab_Mdev_chroma_dkl(pp,:)  = correlator(data, dkl);

    [p2(pp,1), z2(pp,1)] = cmpcorr_ztest(RES.ctab_Mdev_chroma_cone.r(pp), RES.ctab_Mdev_chroma_dkl.r(pp), hu_n, hu_n);
end
RES.ctab_Mdev.Properties.RowNames    = AGG.pps.pp_lbl;
RES.ctab_oppodev.Properties.RowNames = AGG.pps.pp_lbl;

RES.ctab_chroma_cone.Properties.RowNames = AGG.pps.pp_lbl;
RES.ctab_chroma_dkl.Properties.RowNames  = AGG.pps.pp_lbl;

RES.ctab_Mdev_chroma_cone.Properties.RowNames = AGG.pps.pp_lbl;
RES.ctab_Mdev_chroma_dkl.Properties.RowNames  = AGG.pps.pp_lbl;

RES.ctab_chroma_comp        = table(z1, p1, 'RowNames', AGG.pps.pp_lbl);
RES.ctab_Mdev_chroma_comp   = table(z2, p2, 'RowNames', AGG.pps.pp_lbl);

header('Deviations from Opponent', 4);
disp(RES.ctab_oppodev);

header('Deviations from Average', 4);
disp(RES.ctab_Mdev);

header('Chroma & Cone Model', 4);
disp(RES.ctab_chroma_cone);

header('Chroma & DKL Model', 4);
disp(RES.ctab_chroma_dkl);

header('Compare Cone & DKL Chroma correlations', 4);
disp(RES.ctab_chroma_comp);

header('Chroma Deviations from Average - Cone', 4);
disp(RES.ctab_Mdev_chroma_cone);

header('Chroma Deviations from Average - DKL', 4);
disp(RES.ctab_Mdev_chroma_dkl);

header('Compare Cone & DKL MDev Chroma correlations', 4);
disp(RES.ctab_Mdev_chroma_comp);

%% ******************************* FIGURES ********************************

%% MAIN
switch lower(which_one)
    case {'main'}
        header('MAIN FIGURE',3);
        fh = figure('Name', 'Chroma - Main', 'NumberTitle','off');
        pp = 2;
        
        % 1 ILLUSTRATION --------------------------------------------------
        subplot(2,3,1)
        mod_illustration_plotter(AGG.model4illu{pp}, DESIGN);
        
        % 2 STIMULI -------------------------------------------------------
        subplot(2,3,2)
        stimulus_plotter(AGG.stimuli{pp}, DESIGN);

        % 3 MODEL ---------------------------------------------------------
        subplot(2,3,3)
        model_plotter(AGG.model4illu{pp}, AGG.model{pp}, DESIGN);

        % 4 DATA ----------------------------------------------------------
        subplot(2,3,4)
        mod_data_plotter(AGG.data{pp}, AGG.model{pp}, AGG.stimuli{pp}, AGG.model4illu{pp}, DESIGN);

        % 5 DEVIATION FROM OPPONENT ---------------------------------------
        subplot(2,3,5)
        DESIGN.axmax = 50;
        corr_plotter(...
            AGG.data{pp}.OppoDev,...
            AGG.model{pp}.ConeCon.OppoDev,...
            AGG.model{pp}.Adapt,...
            RES.ctab_oppodev(pp,:), DESIGN);
        title('DEVIATION FROM OPPONENT', 'FontWeight', 'bold');

        % 6 DEVIATION FROM MEAN -------------------------------------------
        subplot(2,3,6)
        DESIGN.axmax = 30;
        corr_plotter(...
            AGG.data{pp}.Mdev,...
            AGG.model{pp}.ConeCon.Mdev,...
            AGG.model{pp}.Adapt,...
            RES.ctab_Mdev(pp,:), DESIGN);
        title('DEVIATION FROM MEAN', 'FontWeight', 'bold');

end

%% ADDITIONAL PARTICIPANTS
switch lower(which_one)
    case {'supplementary'}
        header('ADDITIONAL PARTICIPANTS',3);
        pps = [1 3 4 5];
        for nr = 1:numel(pps)
            pp = pps(nr);

            % 1 STIMULI ---------------------------------------------------
            subplot(4,4,(nr-1)*4+1)
            stimulus_plotter(AGG.stimuli{pp}, DESIGN);
            title(AGG.pps.pp_lbl{pp}, 'FontWeight', 'bold');
            ylabel('Blue-Yellow [v*]');

            % 2 DATA ------------------------------------------------------
            subplot(4,4,(nr-1)*4+2)
            mod_data_plotter(AGG.data{pp}, AGG.model{pp}, AGG.stimuli{pp}, AGG.model4illu{pp}, DESIGN);
                ylabel('Blue-Yellow [v*]');
            title('');

            % 3 DEVIATION FROM OPPONENT ---------------------------------------
            subplot(4,4,(nr-1)*4+3)
            DESIGN.axmax = 50;
            corr_plotter(...
                AGG.data{pp}.OppoDev,...
                AGG.model{pp}.ConeCon.OppoDev,...
                AGG.model{pp}.Adapt,...
                RES.ctab_oppodev(pp,:), DESIGN);

            % 4 DEVIATION FROM MEAN -------------------------------------------
            subplot(4,4,(nr-1)*4+4)
            DESIGN.axmax = 30;
            corr_plotter(...
                AGG.data{pp}.Mdev,...
                AGG.model{pp}.ConeCon.Mdev,...
                AGG.model{pp}.Adapt,...
                RES.ctab_Mdev(pp,:), DESIGN);
        end
end

%% ************************* SPECIFIC SUBFUNCTIONS ************************

%% cone_plotter
function cone_plotter(CONES, mx, linestyle, colspace)
% 2022.09.13
% 2024.07.23 added colspace

if nargin < 4
    colspace = 'luv';
    if nargin < 3
        linestyle = '-';
    end
end

% RESCALE -----------------------------------------------------------------
switch lower(colspace)
    case 'luv'
        cones = CONES.Luv;
        cones(:,2) = (cones(:,2)./CONES.Luv_pol(:,2))*mx;
        cones(:,3) = (cones(:,3)./CONES.Luv_pol(:,2))*mx;
    case 'dkl'
        cones = CONES.dkl;
        cones(:,2) = (cones(:,2)./CONES.dkl_pol(:,2))*mx;
        cones(:,3) = (cones(:,3)./CONES.dkl_pol(:,2))*mx;
end

% PLOT --------------------------------------------------------------------
plot([0, cones(1,2)],[0, cones(1,3)],...
    linestyle, 'Color', [1 0 0], 'LineWidth',1);
plot([0, cones(2,2)],[0, cones(2,3)],...
    linestyle, 'Color', [0 0.8 0], 'LineWidth',1);
plot([0, cones(3,2)],[0, cones(3,3)],...
    linestyle, 'Color', [0 0 1], 'LineWidth',1);

%% corr_plotter
function corr_plotter(Data, Model, Adapt, ctab, DESIGN)
% 2025.07.27 [cw]

hold on

plot([-1 1]*DESIGN.axmax, [0 0], '-', 'Color', [.5 .5 .5], 'LineWidth', .5);
scatter(Data(:), Model(:),...
    10, [Adapt.R(:), Adapt.G(:), Adapt.B(:)],'filled','MarkerEdgeColor','k');

hold off
set(gca, 'FontSize', DESIGN.fontsize);
axis([-1 1 -1 1]*DESIGN.axmax);
axis square;
stats_reporter(ctab, 'corr', 'tl',2,'k',DESIGN.fontsize);
xlabel('Predicted Differences [deg]');
ylabel('Measured Differences [deg]');

%% mod_data_plotter
function mod_data_plotter(Data, Model, Stimuli, Illu, DESIGN)
%2025.07.27 [cw]

% DATA PREPARATION --------------------------------------------------------

% MONITOR GAMUT + AXIS LENGTH
gamut = Stimuli.gamut.Luv; % Monitor gamut
gamax = [floor(min(gamut(:,2))), ceil(max(gamut(:,2))), floor(min(gamut(:,3))), ceil(max(gamut(:,3)))];

% DATA:
hueM = Data.hueM; % Average hue data
chrM = Data.chrM; % Average chroma data
[data_u, data_v] = pol2cart(deg2rad(hueM), chrM);
N = sum(~isnan(hueM),'all'); % Number of stimuli

% AVERAGE LINE:
% To identify hue deviations from this average line 
MaziM = Data.MhueM(:,1); 
[Mu, Mv] = pol2cart(deg2rad(MaziM), 300);

% INDUCERS:
ihue = Data.iHue; % Inducer hue
ichr = Data.iChr; % Inducer chroma
[i_u,i_v] = pol2cart(deg2rad(ihue), ichr);
ihueM = nanmean(ihue,2); % Make vector out of reptive columns

% MODEL / SIMULATION: 
ahue = round(nanmean(Model.Adapt.hue,2));
ConeCon = Model.ConeCon; % Make vector out of reptive columns
[~, inds] = intersect(Illu.Adapt.hue(:,1),ahue);
ConeCon2 = Illu.ConeCon(inds,:);
mod_chr = [ConeCon.chr,ConeCon2.chr];
mod_u = [ConeCon.u,ConeCon2.u];
mod_v = [ConeCon.v,ConeCon2.v];

% PLOT --------------------------------------------------------------------

hold on
plot(gamut(:,2),gamut(:,3),':', 'Color', [1 1 1]*0);
%plot(data_u(:), data_v(:),  'k.', 'MarkerSize', 5, 'Color', [1 1 1]*.7);
scatter(data_u(:),data_v(:), 1,'k','filled','MarkerFaceAlpha', 1);

for azi = 0:45:345
    hu = ihueM == azi;
    plot([0 Mu(hu)], [0 Mv(hu)], '-', 'Color', [1 1 1]*.5, 'LineWidth', 0.5);

    mod_uv = sortrows([mod_u(hu,:)',mod_v(hu,:)',mod_chr(hu,:)'],3);
    plot([0; mod_uv(:,1)], [0; mod_uv(:,2)], 'r-', 'LineWidth', 1);

    inds = ~isnan(data_u(hu,:));
    h = plot([0 data_u(hu,inds)],  [0 data_v(hu,inds)], 'k-', 'LineWidth', 2, 'Color', [0 0 0 .4]);
end

hold off
%stats_reporter(N, 'N', 'TR', 2, [.5 .5 .5], DESIGN.fontsize);
set(gca, 'FontSize', DESIGN.fontsize);
axis(gamax);
axis square;
title('DATA', 'FontWeight', 'bold');

%% mod_illustration_plotter
function mod_illustration_plotter(Model4Illu, DESIGN)
% 2025.08.19 [cw]

% SETTINGS ----------------------------------------------------------------
 
% EXAMPLE
ex_azi = 60;
[indsH, ~] = find(Model4Illu.Adapt.hue == ex_azi, 1);
ex_chro = [20 30 40 50 70];
[~, indsC] = intersect(Model4Illu.Adapt.chr(indsH,:), ex_chro);
chr_n = numel(indsC);
L = ones(chr_n,1)*unique(Model4Illu.Adapt.L);

% INDUCER
H = Model4Illu.Adapt.hue(indsH,indsC)';
iC = Model4Illu.Adapt.chr(indsH,indsC)';
INDU = colourconverter([L,H,iC], 'Luv_pol', 2);
INDU.rgb = INDU.rgb/255;

% CONE MODEL
H = Model4Illu.ConeCon.hue(indsH,indsC)';
%C = Model4Illu.ConeCon.chr(indsH,indsC)';
simCOL = colourconverter([L,H,iC], 'Luv_pol', 2);
simCOL.rgb = simCOL.rgb/255;

% AVERAGE FOR COMPARISON
compAzi = ones(chr_n,1)*mean(simCOL.Luv_pol(:,1));
compLuv = [L,compAzi,simCOL.Luv_pol(:,2)];
COMP =  colourconverter(compLuv, 'luv_pol',2);
COMP.rgb = COMP.rgb/255;

% CONE-OPPONENT
oLuv = [L,ones(chr_n,1)*(ex_azi-180),iC];
OPPO = colourconverter(oLuv, 'Luv_pol',2);
OPPO.rgb = OPPO.rgb/255;

header('INDUCER RENDERING',4);
disp(INDU);

header('OPPO RENDERING',4);
disp(OPPO);

header('MODEL RENDERING',4);
disp(simCOL);

header('COMP RENDERING',4);
disp(COMP);


% PRODUCTION --------------------------------------------------------------
rect = [100 500];
bands = rect(2)/chr_n;
hold on
for chr = 1:chr_n
    fill(-rect(2)/2 + [(chr-1)*bands chr*bands chr*bands (chr-1)*bands],[0 0 rect(1) rect(1)]+1.5*rect(1),INDU.rgb(chr,:));
    fill(-rect(2)/2 + [(chr-1)*bands chr*bands chr*bands (chr-1)*bands],[0 0 rect(1) rect(1)]+.5*rect(1),OPPO.rgb(chr,:));
    fill(-rect(2)/2 + [(chr-1)*bands chr*bands chr*bands (chr-1)*bands],[0 0 rect(1) rect(1)]-.5*rect(1),simCOL.rgb(chr,:));
    fill(-rect(2)/2 + [(chr-1)*bands chr*bands chr*bands (chr-1)*bands],[0 0 rect(1) rect(1)]-1.5*rect(1),COMP.rgb(chr,:));
    x = -rect(2)/2 + (chr-.5)*bands;
    text(x,2.5*rect(1),sprintf('%d',round(INDU.Luv_pol(chr,2))),...
        'Units','Data','FontSize', DESIGN.fontsize,...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom');
    text(x,2*rect(1),sprintf('%d',round(INDU.Luv_pol(chr,1))),...
        'Units','Data','FontSize', DESIGN.fontsize,...
        'HorizontalAlignment','Center','VerticalAlignment','Middle');
    text(x,rect(1),sprintf('%d',round(OPPO.Luv_pol(chr,1))),...
        'Units','Data','FontSize', DESIGN.fontsize,...
        'HorizontalAlignment','Center','VerticalAlignment','Middle');
    text(x,0,sprintf('%d',round(simCOL.Luv_pol(chr,1))),...
        'Units','Data','FontSize', DESIGN.fontsize,...
        'HorizontalAlignment','Center','VerticalAlignment','Middle');
    text(x,-rect(1),sprintf('%d',round(COMP.Luv_pol(chr,1))),...
        'Units','Data','FontSize', DESIGN.fontsize,...
        'HorizontalAlignment','Center','VerticalAlignment','Middle');
end
hold off
text(-rect(2)/2,rect(1)*2,'Inducer ',...
    'Units','Data','FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','Right','VerticalAlignment','Middle');
% text(rect(2)/2,rect(1)*2,' Constant',...
%     'Units','Data','FontSize', DESIGN.fontsize,...
%     'HorizontalAlignment','Left','VerticalAlignment','Middle');

text(-rect(2)/2,rect(1),'Opponent ',...
    'Units','Data','FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','Right','VerticalAlignment','Middle');
% text(rect(2)/2,rect(1),' Constant',...
%     'Units','Data','FontSize', DESIGN.fontsize,...
%     'HorizontalAlignment','Left','VerticalAlignment','Middle');


text(-rect(2)/2,0,'Prediction ',...
    'Units','Data','FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','Right','VerticalAlignment','Middle');
% text(rect(2)/2,0,' Changing',...
%     'Units','Data','FontSize', DESIGN.fontsize,...
%     'HorizontalAlignment','Left','VerticalAlignment','Middle');

text(-rect(2)/2,-rect(1),'Average ',...
    'Units','Data','FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','Right','VerticalAlignment','Middle');
% text(rect(2)/2,-rect(1),' Constant',...
%     'Units','Data','FontSize', DESIGN.fontsize,...
%     'HorizontalAlignment','Left','VerticalAlignment','Middle');

axis off;

%% model_plotter
function model_plotter(Illu, Model, DESIGN)
% 2025.07.26 [cw]

N = sum(~isnan(Model.ConeCon.u),'all');
axmx = max([Illu.ConeCon.chr(:); Model.ConeCon.chr(:)]);

hold on
cone_plotter(Illu.cones, axmx+50, '-', 'Luv');

% ILLUSTRATION INDUCERS ---------------------------------------------------
fill(Illu.Adapt.u, Illu.Adapt.v, 'w');
plot(Illu.Adapt.u, Illu.Adapt.v, '-', 'Color', [1 1 1]*.7)
scatter(Model.ConeCon.u(:), Model.ConeCon.v(:), 10, [Model.Adapt.R(:),Model.Adapt.G(:),Model.Adapt.B(:)], 'filled', 'MarkerEdgeColor', 'none');

% ILLUSTRATION MODEL ------------------------------------------------------
plot(Illu.ConeCon.u, Illu.ConeCon.v, 'k-', 'LineWidth', 0.5);
for az = 1:45:355
    plot(Illu.ConeCon.u(az,:),Illu.ConeCon.v(az,:), 'k-', 'LineWidth', 0.5);
end

% MODEL OF STIMULI --------------------------------------------------------
% plot(Model.ConeCon.u, Model.ConeCon.v, '.r');
% 
% scatter(Model.ConeCon.u(:), Model.ConeCon.v(:), 10, [Model.Adapt.R(:),Model.Adapt.G(:),Model.Adapt.B(:)], 'filled', 'MarkerEdgeColor', 'none');

hold off

%stats_reporter(N, 'N', 'TR', 2, 'r', DESIGN.fontsize);
set(gca, 'FontSize', DESIGN.fontsize);
axis([-1 1 -1 1]*axmx);
axis square;
title('MODEL', 'FontWeight', 'bold');

%% radius_plotter
function radius_plotter(radius,rgb)
%2020.06.24
%2022.10.14 multiple radii option
if nargin < 2
    rgb = [.5 .5 .5];
end

azi = deg2rad(0:360)';
for rd = 1:numel(radius)
    [x,y] = pol2cart(azi, ones(numel(azi),1)*radius(rd));
    plot(x,y, 'Color', rgb);
end

%% stimulus_plotter
function stimulus_plotter(Stimuli, DESIGN)
% 2025.07.26 polished for pub [cw]

% PREPARE DATA ------------------------------------------------------------
% INDUCERS
% ihue = Data.iHue; % Inducer hue
% ichr = Data.iChr; % Inducer chroma
% [i_u,i_v] = pol2cart(deg2rad(ihue), ichr);
i_u = Stimuli.inducer.Luv(:,2);
i_v = Stimuli.inducer.Luv(:,3);
rgb = Stimuli.inducer.rgb;

% MONITOR GAMUT
gamut = Stimuli.gamut.Luv;
gamax = [floor(min(gamut(:,2))), ceil(max(gamut(:,2))), floor(min(gamut(:,3))), ceil(max(gamut(:,3)))];

gamax = [-1 1 -1 1] * max(abs([gamut(:,2);gamut(:,3)]));

% PLOT --------------------------------------------------------------------
hold on

% CIRCLES FOR INDUCER CHROMA
chrs = 0:10:80;
for chr = 1:numel(chrs)
    radius_plotter(chrs(chr), [1 1 1]*.7);
end
plot(gamut(:,2),gamut(:,3),':', 'Color', [1 1 1]*0); % GAMUT
scatter(i_u, i_v, 10, rgb, 'filled', 'MarkerEdgeColor', 'none'); % DATA
hold off
set(gca, 'FontSize', DESIGN.fontsize);
N = numel(i_u(:));
stats_reporter(N, 'N', 'TR', 2, 'k', DESIGN.fontsize);
axis(gamax);
axis square;
title('INDUCERS', 'FontWeight', 'bold');

%% ************************* GENERAL SUBFUNCTIONS *************************

%% cmpcorr_ztest
function [p, z] = cmpcorr_ztest(r1, r2, N1, N2)
% 2014oct25 * [cw]

if nargin == 0
    N1 = 320;
    N2 = 320;
    r1 = 0.25;
    r2 = 0.27;
end

z1 = fishertransformer(r1);
z2 = fishertransformer(r2);
    
z = (z1-z2) / sqrt( (1/(N1-3)) + 1/(N2-3) );

p = 2*(1-normcdf(abs(z),0,1));
%% correlator
function tab = correlator(X, Y, ns, onetailor, type, lbls)
%2025.07.16 adapted for publication [cw] 

if nargin < 5
    type = 'Pearson';
    if nargin < 4
        onetailor = 0;
        if nargin < 3
            ns = 'ns';
        end
    end
end
    
y_n = size(Y,2); 

for y = 1:y_n
    N = sum(~isnan(X) & ~isnan(Y(:,y)));
    df(y,1) = N-2;
    [r(y,1), p(y,1)] = corr(X,Y(:,y),'rows','pairwise','Type',type);
end

R2 = r.^2*100;

sig = significancer(p,ns,onetailor);

% Confidence intervals ----------------------------------------------------
[rlo, rup] = ci_calculator(r, df+2);
CI  = [rlo, rup];
tab = table(r, df, CI, p, sig, R2); 

if nargin == 6
    tab.Properties.RowNames = lbls;
end

%% ci_calculator
function [rlo, rup] = ci_calculator(R, N)
% Taken from Matlab's corrcoef
% Double-checked with https://www.statskingdom.com/correlation-confidence-interval-calculator.html
% 2025.07.06 * [cw]

alpha = 0.05;

% Confidence bounds are degenerate if abs(r) = 1, NaN if r = NaN.
z = 0.5 * log((1+R)./(1-R));
zalpha = NaN(size(N));
rlo = NaN(size(N));
rup = NaN(size(N));
if any(N>3)
    zalpha(N>3) = (-erfinv(alpha - 1)) .* sqrt(2) ./ sqrt(N(N>3)-3);
end
rlo = tanh(z-zalpha);
rup = tanh(z+zalpha);

%% fishertransformer
function zz = fishertransformer(rr, inverter)
%INTPUT--------------------------------------------------------------------
% rr        = Typically correlation coefficients (r).
% inverter  = Set to 1 to do the inverse transformation (r to z); default:
%             inverter = 0;
%OUTPUT--------------------------------------------------------------------
% zz = Fisher's z-transform (z).
%GENEALOGY-----------------------------------------------------------------
%2011dec15  * [cw]

if nargin < 2
    inverter = 0;
end

if inverter == 0
    zz = 1/2*log((1+rr)./(1-rr));
elseif inverter == 1
    z = rr;
    r = (exp(2*z)-1)./(exp(2*z)+1);
    zz = r;
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

%% significancer
function [sym, num] = significancer(p, ns, onetailor, alpha)
% 2017.02.23 [cw]

if nargin < 4
    alpha = 0.05;
    if nargin < 3
        onetailor = 0;
        if nargin < 2
            ns = 'ns';
        end
    end
end

if onetailor
    p = p/2;
end

n = size(p);

% NOT SIGNIFICANT --------------------------------------------------
num = zeros(n);
switch lower(ns)
    case 'p'
        sym = cellstr(num2str(p, '%.2f'));
    otherwise
        sym(1:n(1),1:n(2)) = cellstr(ns);
end

% AGAINST HYPOTHESIS (one-tailed) -----------------------------------------
inds = p == 1;
sym(inds) = {''};
num(inds) =  1;

% MARGINALLY SIGNIFICANT --------------------------------------------------
inds = p < alpha * 2; % default 0.1
sym(inds) = {'°'};
num(inds) =  0.5;

% SIGNIFICANT -------------------------------------------------------------
inds = p < alpha; % default 0.05
num(inds) = 1;
sym(inds) = {'*'};

% CLEARLY SIGNIFICANT -----------------------------------------------------
inds = p < alpha/5; % default 0.01
num(inds) =  2;
sym(inds) = {'**'};

% HIGHLY SIGNIFICANT ------------------------------------------------------
inds = p < alpha/50; % default 0.001
num(inds) =  3;
sym(inds) = {'***'};

%% stats_reporter
function stats_reporter(statstab, statstype, location, shift, rgb, fsz)
%2021.07.01 [cw]

if nargin < 6
    fsz = 10;
if nargin < 5
    rgb = 'k';
if nargin < 4
    shift = 1;
    if nargin < 3
        location = 'topleft';
        if nargin < 2
            statstype = 'corr';
            if nargin < 1
                statstab = dataset(1, 100, {'***'}, 'VarNames', {'r', 'df','sig'});
            end
        end
    end
end
end
end

switch lower(location)
    case {'topleft','tl'}
        x = 0;
        y = 1;
        HorizontalAlignment = 'left';
        VerticalAlignment = 'top';
        lshift = blanks(shift);
        rshift = '';
    case {'bottomleft','bl'}
        x = 0;
        y = 0;
        HorizontalAlignment = 'left';
        VerticalAlignment = 'bottom';
        lshift = blanks(shift);
        rshift = '';
    case {'topright','tr'}
        x = 1;
        y = 1;
        HorizontalAlignment = 'right';
        VerticalAlignment = 'top';
        lshift = '';
        rshift = blanks(shift);
    case {'bottomright','br'}
        x = 1;
        y = 0;
        HorizontalAlignment = 'right';
        VerticalAlignment = 'bottom';
        lshift = '';
        rshift = blanks(shift);
    otherwise
        error('No valid location setting');
end
switch lower(statstype)
    case {'n','samplesize'}
        txt = sprintf('%sN=%d%s', lshift, statstab, rshift);
    case {'expvar','explainedvariance'}
        txt = sprintf('%s%.2f%%%s', lshift, statstab, rshift);
    case {'corr','correlation'}
        txt = sprintf('%sr(%d)=%.2f%s%s', lshift, statstab.df(1), statstab.r(1), statstab.sig{1}, rshift);
    case {'reg','regression'}
        txt = sprintf('%sR^2=%.0f%%%s%s', lshift, round(statstab.R2(1)), statstab.sig{1}, rshift);
    case {'ttest'}
        txt = sprintf('%st(%d)=%.1f%s%s', lshift, statstab.df(1), statstab.t(1), statstab.sig{1}, rshift);
end

text(x,y, txt,...
    'Units', 'Normalized',...
    'Color', rgb, 'FontSize', fsz,...
    'HorizontalAlignment', HorizontalAlignment, 'VerticalAlignment', VerticalAlignment);
