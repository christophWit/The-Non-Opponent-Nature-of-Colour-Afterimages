function [RES, fh] = mod_exp2_analysor(AGG, which_one)
%2025.07.16 adapted + polished for publication [cw]

%% ******************************* SETTINGS *******************************

AGG.colspace = AGG.stimuli{1}.src_space; % Represent results in the colour space in which data was collected.
fh = [];

% INITIALISE OUTPUT STRUCTURE WITH RESULTS
RES = struct; 


%% DESIGN (cosmetics)
DESIGN.fontsize  = 9;
DESIGN.fontsize2 = 10;
DESIGN.radius    = 50;
DESIGN.axmax     = 51;
DESIGN.sim_rgb   = 'r';
DESIGN.data_rgb  = [1 1 1]*0.5;
DESIGN.examples  = 0:30:330; % For distributions

%% ************************** GRAPHICAL ABSTRACT **************************
switch lower(which_one)
    case {'graphical abstract'}
        fh = figure('Name', 'Graphical Abstract', 'NumberTitle','off');
        mod_graphicalabstractor(AGG);

end

%% ******************************** METHOD ********************************

switch lower(which_one)
    case {'method'}
        header('METHOD',3);

        % PARTICIPANTS ----------------------------------------------------
        RES = pp_reporter(RES, AGG.pps);
end

%% ************************* AGGEGATED (Main) *****************************
RES.agg_ctab(1,:)        = correlator(AGG.agg.chroma(:,1), AGG.agg.sim_chrM);
RES.agg_ctab_smooth(1,:) = correlator(AGG.agg.chroma_smooth, AGG.agg.sim_chrM);

RES.agg_ctab(2,:)        = correlator(AGG.agg.oppdev(:,1), AGG.agg.sim_oppdev);
RES.agg_ctab_smooth(2,:) = correlator(AGG.agg.oppdev_smooth, AGG.agg.sim_oppdev);

RES.agg_ctab(3,:)        = correlator(AGG.agg.rF, AGG.agg.sim_rF);
RES.agg_ctab_smooth(3,:) = correlator(AGG.agg.rF_smooth, AGG.agg.sim_rF);

RES.agg_ctab.smooth_R2 = RES.agg_ctab_smooth.R2;
RES.agg_ctab.Properties.RowNames = {'Chroma', 'OppoDev', 'HueF'};
RES.agg_ctab_smooth.Properties.RowNames = {'Chroma', 'OppoDev', 'HueF'};

switch lower(which_one)
    case {'main', 'aggregated', 'agg'}
        header('MAIN: Aggregated',3);

        % TABLES ----------------------------------------------------------
        header('Correlations', 4);
        disp(RES.agg_ctab);

        header('Correlations - Smoothed', 4);
        disp(RES.agg_ctab_smooth);
        
        % FIGURE ----------------------------------------------------------
        figure('Name', 'Agg-Main', 'NumberTitle','off');
        mod_agg_plotter(AGG, RES.agg_ctab, DESIGN);
end

%% ************************* INDIVIDUAL RESULTS ***************************

header('INDIVIDUAL RESULTS',3);

%% CORRELATION TABLES
pp_n = size(AGG.pps,1);
for pp = 1:pp_n

    % CHROMA --------------------------------------------------------------
    RES.IndiM.ctab(pp,:)        = correlator(AGG.IndiM.chrM(:,pp), AGG.model.Chroma.ConeCon(:,pp));
    RES.IndiM.ctab_smooth(pp,:) = correlator(AGG.IndiM.chrM_smooth(:,pp), AGG.model.Chroma.ConeCon(:,pp));
    
    % DEVIATIONS FROM OPPONENCY -------------------------------------------
    RES.IndiOppoDev.ctab(pp,:)          = correlator(AGG.IndiOppoDev.sim(:,pp), AGG.IndiOppoDev.adj(:,pp));
    RES.IndiOppoDev.ctab_smooth(pp,:)   = correlator(AGG.IndiOppoDev.sim(:,pp), AGG.IndiOppoDev.adj_smooth(:,pp));

    % HUE FREQUENCIES -----------------------------------------------------
    RES.IndiHueFrq.ctab(pp,:)           = correlator(AGG.hueF.rF(:,pp), AGG.model.hueF.rF(:,pp));
    RES.IndiHueFrq.ctab_smooth(pp,:)    = correlator(AGG.hueF.rF_smooth(:,pp), AGG.model.hueF.rF(:,pp));

end
RES.IndiM.ctab.smooth_R2                    = RES.IndiM.ctab_smooth.R2;
RES.IndiM.ctab.Properties.RowNames          = AGG.pps.pp_lbl;
RES.IndiM.ctab_smooth.Properties.RowNames   = AGG.pps.pp_lbl;

RES.IndiOppoDev.ctab.smooth_R2                  = RES.IndiOppoDev.ctab_smooth.R2;
RES.IndiOppoDev.ctab.Properties.RowNames        = AGG.pps.pp_lbl;
RES.IndiOppoDev.ctab_smooth.Properties.RowNames = AGG.pps.pp_lbl;

RES.IndiHueFrq.ctab.Properties.RowNames         = AGG.pps.pp_lbl;
RES.IndiHueFrq.ctab.smooth_R2                   = RES.IndiHueFrq.ctab_smooth.R2;
RES.IndiHueFrq.ctab_smooth.Properties.RowNames  = AGG.pps.pp_lbl;

%% DISPLAY CORRELATION TABLES

switch lower(which_one)
    case {'correlation tables', 'corrtables'}
        header('CORRELATION TABLES',3);

        header('Aggregated', 4);
        disp(RES.agg_ctab);

        header('Aggregated - Smoothed', 4);
        disp(RES.agg_ctab_smooth);

        % CHROMA ----------------------------------------------------------        
        header('Predicted vs Measured Chroma', 4);
        disp(RES.IndiM.ctab);

        % DEVIATIONS FROM OPPONENCY ---------------------------------------        
        header('Deviations from OPPONENCY', 4);
        disp(RES.IndiM.ctab);
        
        % HUE HISTOGRAM ---------------------------------------------------        
        header('HUE HISTOGRAM', 4);
        disp(RES.IndiHueFrq.ctab);
end

%% INDIVIDUAL DATA AND AVERAGES
switch lower(which_one)
    case {'individual adjustments', 'indiadj'}
        header('Individual Adjustments',3);
        fh = figure('Name', 'Indi-Adj', 'NumberTitle','off');
        mod_indiM_plotter(AGG, RES.IndiM.ctab, DESIGN);
end

%% INDIVIDUAL DEVIATIONS FROM OPPONENCY
switch lower(which_one)
    case {'individual deviations from opponency', 'indioppodev'}
        header('Individual Deviations from Opponency',3);
        fh = figure('Name', 'Indi-OppoDev', 'NumberTitle','off');
        mod_indiOppoDev_plotter(AGG, RES.IndiOppoDev.ctab, DESIGN);
end

%% INDIVIDUAL HUE HISTOGRAMS
switch lower(which_one)
    case {'individual hue histograms', 'indihuehist'}
        header('FIG: Individual Hue Histograms',3);                
        figure('Name', 'Indi-HueHist', 'NumberTitle','off');
        mod_indiHueFreq_plotter(AGG, RES.IndiHueFrq.ctab, DESIGN);
end

%% DISTRIBUTION
switch lower(which_one)
    case {'distribution hue', 'huedist', 'distribution'}
        header('FIG: Hue Distribution',3);                
        figure('Name', 'Hue Distribution', 'NumberTitle','off');
        mod_distribution_plotter_hue(AGG, DESIGN);
end
switch lower(which_one)
    case {'distribution chroma', 'chrodist', 'distribution'}
        header('FIG: Chroma Distribution',3);
        figure('Name', 'Chroma Distribution', 'NumberTitle','off');
        mod_distribution_plotter_chroma(AGG, DESIGN);
end

%% QQ PLOTS
switch lower(which_one)
    case {'qq hue'}
        figure('Name', 'QQ Hue', 'NumberTitle','off');
        mod_qq_plotter_hue(AGG, DESIGN);
end
switch lower(which_one)
    case {'qq chroma'}
        figure('Name', 'QQ Chroma', 'NumberTitle','off');
        mod_qq_plotter_chroma(AGG, DESIGN);
end

%% ****************************** SUBMODULES ******************************

%% mod_distribution_plotter_chroma
function mod_distribution_plotter_chroma(AGG, DESIGN)
% 2025.07.21 [cw]

examples = DESIGN.examples;
alldata = cat(2, AGG.Indi{:});
allchroma = alldata(:,:,5);
oppo_chroma = AGG.model.Chroma.dkl(:,1);

ihue = (0:5:355)';
ex_n = numel(examples);
for ex = 1:ex_n
    ax(ex) = subplot(2,ex_n/2,ex);
    inds = ihue == examples(ex);
    chroF = histcounts(allchroma(inds,:), 0:5:70);
    x = (2.5:5:67.5)';
    hold on
    plot([1 1]*oppo_chroma(inds), [0 20], 'r-', 'LineWidth', 1);
    h = bar(x, chroF, 'FaceColor', 'flat');
    h.CData = x*[1 1 1]/70;
    hold off

    % FORMAT --------------------------------------------------------------
    mx(ex) = max(chroF);
    title(sprintf('%d deg Inducer',examples(ex)),'FontWeight','bold');
    if ex == 1
        N = sum(~isnan(allchroma(inds,:)));
        stats_reporter(N,'N', 'topleft', 1, 'k', DESIGN.fontsize);
    end
    if ex <= ex_n/2
        ax(ex).XTickLabels = [];
    end
    if ex == ceil(1.5*ex_n) % Bottom Centre
        xlabel('Chroma Bins');
    elseif mod(ex, ex_n/2) == 1 % First panel of row
        ylabel('Frequency');
    end    
end

% FORMAT --------------------------------------------------------------
axis(ax, 'square');
set(ax, 'FontSize', DESIGN.fontsize,...
    'XLim', [0 70], 'XTick', 0:10:70,...
    'Ylim', [0 15], 'YTick', 0:5:20);
%mx = max(mx);
%ylim(ax, [0 mx+1]);

%% mod_distribution_plotter_hue
function mod_distribution_plotter_hue(AGG, DESIGN)
% 2025.07.21 [cw]

examples = DESIGN.examples;
alldata = cat(2, AGG.Indi{:});
allhue = alldata(:,:,4);
oppo_hue = AGG.model.Hue.dkl(:,1);
rgb = AGG.stimuli{1}.inducer.rgb;

ihue = (0:5:355)';
ex_n = numel(examples);
for ex = 1:ex_n
    ax(ex) = subplot(2,ex_n/2,ex);
    inds = ihue == examples(ex);
    hueF = hue_frequency_calculator(allhue(inds,:), 5);
    mx(ex) = max(hueF.frq);
    
    hold on
    plot([1 1]*oppo_hue(inds), [0 20], 'r-', 'LineWidth', 1);
    h = bar([ihue; 360], [hueF.frq;hueF.frq(1)], 'FaceColor', 'flat');
    h.CData = [rgb;rgb(1,:)];
    plot(examples(ex), 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', rgb(inds,:));
    hold off

    % FORMAT --------------------------------------------------------------
    title(sprintf('%d deg Inducer',examples(ex)),'FontWeight','bold');
    if ex == 1
        N = sum(~isnan(allhue(inds,:)));
        stats_reporter(N,'N', 'topleft', 1, 'k', DESIGN.fontsize);
    end
    if ex <= ex_n/2
        ax(ex).XTickLabels = [];
    end
    if ex == ceil(1.5*ex_n) % Bottom Centre
        xlabel('Hue Bins [deg]');
    elseif mod(ex, ex_n/2) == 1 % First panel of row
        ylabel('Frequency');
    end    
end

% FORMAT ------------------------------------------------------------------
axis(ax, 'square');
mx = max(mx);
set(ax, 'FontSize', DESIGN.fontsize,...
    'XLim', [0 360],  'XTick', 0:90:360,...
    'Ylim', [0 mx+1], 'YTick', 0:5:20);

%% mod_qq_plotter_chroma
function mod_qq_plotter_chroma(AGG, DESIGN)
% 2025.07.28

examples = DESIGN.examples;
alldata = cat(2, AGG.Indi{:});
allchroma = alldata(:,:,5);
rgb = AGG.stimuli{1}.inducer.rgb;

ihue = (0:5:355)';
ex_n = numel(examples);
for ex = 1:ex_n
    ax(ex) = subplot(2,ex_n/2,ex);
    inds = ihue == examples(ex);
    h = normplot(allchroma(inds,:)');
    h(1).Marker ='o';
    h(1).MarkerEdgeColor ='k';
    h(1).MarkerFaceColor = rgb(inds,:);

    % FORMAT --------------------------------------------------------------
    set(ax(ex), 'FontSize', DESIGN.fontsize);
    title(sprintf('%d deg Inducer',examples(ex)),'FontWeight','bold');
    axis square;
    if ex == 1
        N = sum(~isnan(allchroma(inds,:)));
        stats_reporter(N,'N', 'topleft', 1, 'k', DESIGN.fontsize);
    end
end

%% mod_qq_plotter_hue
function mod_qq_plotter_hue(AGG, DESIGN)
% 2025.07.28

examples = DESIGN.examples;
alldata = cat(2, AGG.Indi{:});
allhue = alldata(:,:,4);
oppo_hue = mod(round(AGG.model.Hue.dkl(:,1),10),360);
rgb = AGG.stimuli{1}.inducer.rgb;

ihue = (0:5:355)';
ex_n = numel(examples);
for ex = 1:ex_n
    ax(ex) = subplot(2,ex_n/2,ex);
    hold on
    inds = ihue == examples(ex);
    h = normplot(cdissolver(allhue(inds,:)'));
    plot([1 1]*oppo_hue(inds), [-3 3], 'r-', 'LineWidth', 1);
    h(1).Marker ='o';
    h(1).MarkerEdgeColor ='k';
    h(1).MarkerFaceColor = rgb(inds,:);
    hold off

    % FORMAT --------------------------------------------------------------
    title(sprintf('%d deg Inducer',examples(ex)),'FontWeight','bold');
    set(ax(ex), 'FontSize', DESIGN.fontsize);
    axis square;
    if ex == 1
        N = sum(~isnan(allhue(inds,:)));
        stats_reporter(N,'N', 'topleft', 1, 'k', DESIGN.fontsize);
    end
end

%% mod_graphicalabstractor(AGG);
function mod_graphicalabstractor(AGG)
% 2025.07.22

chroma      = 50;
radius      = 1.3;
markersize  = 50;
DESIGN.axmx = ceil(max(AGG.agg.rF));
DESIGN.axmx = 1.5;
%DESIGN.fontsize = 9;
L = AGG.stimuli{1}.bg.(AGG.colspace)(1);

% DIAGRAM -----------------------------------------------------------------
colourcircler(AGG.colspace,L,chroma,0:360, 100, radius);
% h = huefrq_plotter(...
%     AGG.agg.rF, ...
%     AGG.agg.rF_smooth, ...
%     AGG.agg.sim_rF,...
%     AGG.dir(1).cones,...
%     AGG.colspace,...
%     DESIGN);

h = [];
ihue = (0:5:355)';
ihue_rad = deg2rad(ihue);

rF          = AGG.agg.rF;
rF_smooth   = AGG.agg.rF_smooth;
sim_rF      = AGG.agg.sim_rF;
cones       = AGG.dir(1).cones;
colspace    = AGG.colspace;


hold on
axes_plotter(DESIGN.axmx+1);
cone_plotter(cones, DESIGN.axmx+1, '-', colspace, 3);

% HUE FREQ DATA -----------------------------------------------------------
[x,y] = pol2cart(ihue_rad,rF);
[x,y] = pol2cart(ihue_rad, rF_smooth);
h_data = fill(x,y,[1 1 1]*.9, 'EdgeColor', [1 1 1]*.5, 'LineWidth', 3);
%        h_fill = huefreq_plotter(ihue, AGG.hueF.rF(:,pp));

% SIMULATION OF ADAPTATION ------------------------------------------------
[x,y] = pol2cart(ihue_rad, sim_rF);
h_sim = fill(x,y,[.8 .8 .8], 'EdgeColor', [1 0 0], 'FaceColor', 'none', 'LineWidth', 3);

% SMOOTH ------------------------------------------------------------------
% [x,y] = pol2cart(ihue_rad, rF_smooth);
% h_smooth = fill(x,y,[.8 .8 .8], 'EdgeColor', 'k', 'FaceColor', 'none', 'LineWidth',1);
hold off
text(0, 1, ' PERCEPTION',...
    'Units', 'normalized',...
    'Color', [1 1 1]*.5, 'FontSize', 14, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
text(1, 1, 'CONE MODEL ',...
    'Units', 'normalized',...
    'Color', [1 0 0], 'FontSize', 14, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top');
axis([-1 1 -1 1] * DESIGN.axmx);
axis square;

set(gca, 'XTick', [], 'YTick', []);
box on;

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
frq_smooth = smoothdata(frq,'movmean',smoother);
rF_smooth   = smoothdata(rF,'movmean',smoother);

hueF = table(N, frq, rF, frq_smooth, rF_smooth);



%% huefrq_plotter
function h = huefrq_plotter(rF, rF_smooth, rf_sim, cones, colspace, DESIGN)
% 2025.07.19 [cw]

h = [];
ihue = (0:5:355)';
ihue_rad = deg2rad(ihue);

hold on
axes_plotter(DESIGN.axmx+1);
cone_plotter(cones, DESIGN.axmx+1, '-', colspace);

% HUE FREQ DATA -----------------------------------------------------------
[x,y] = pol2cart(ihue_rad,rF);
h_data = fill(x,y,[1 1 1]*.9, 'EdgeColor', [1 1 1]*.5, 'LineWidth', 0.5);
%        h_fill = huefreq_plotter(ihue, AGG.hueF.rF(:,pp));

% SIMULATION OF ADAPTATION ------------------------------------------------
[x,y] = pol2cart(ihue_rad, rf_sim);
h_sim = fill(x,y,[.8 .8 .8], 'EdgeColor', [1 0 0], 'FaceColor', 'none', 'LineWidth',2);

% SMOOTH ------------------------------------------------------------------
[x,y] = pol2cart(ihue_rad, rF_smooth);
h_smooth = fill(x,y,[.8 .8 .8], 'EdgeColor', 'k', 'FaceColor', 'none', 'LineWidth',1);
hold off
axis([-1 1 -1 1] * DESIGN.axmx);
axis square;
set(gca, 'FontSize', DESIGN.fontsize, 'XTick', -3:3, 'YTick', -3:3);

% HANDLES -----------------------------------------------------------------
h = [h_sim; h_data; h_smooth];

%% Magg_plotter
function h = Magg_plotter(indi, indiM_HC, indiM_smooth, ConeCon_HC, cones, colspace, DESIGN)
% 2025.07.19 [cw]

hold on

axes_plotter(DESIGN.axmx+1);
cone_plotter(cones, DESIGN.axmx*sqrt(2), '-', colspace);

% AVERAGE -----------------------------------------------------------------
xy_pol = sortrows(indiM_HC,1);
[x,y] = pol2cart(deg2rad(xy_pol(:,1)),xy_pol(:,2));
%        h_Mdata = plot(af(x),af(y),'-', 'Color',data_rgb, LineWidth',2);
h_Mdata = fill(x, y, [1 1 1]*.9, 'EdgeColor', [1 1 1]*0, 'LineWidth', 1);

% DATA --------------------------------------------------------------------
x = indi(:,:,2);
y = indi(:,:,3);
h_data = scatter(x(:),y(:), 1,'k','filled','MarkerFaceAlpha',.3);
%h_data = scatter(x(:),y(:), 1,'k','filled');


% SIMULATION OF CONE ADAPTATION -------------------------------------------
xy_pol = sortrows(ConeCon_HC,1);
[x,y] = pol2cart(deg2rad(xy_pol(:,1)),xy_pol(:,2));
% fill(x, y, 'w', 'EdgeColor', sim_rgb); % To cover the lines of the cone directions
% h_sim = fill(x, y, sim_rgb, 'EdgeColor', sim_rgb);
% h_sim.FaceAlpha = .5;
h_sim = fill(x, y, 'w', 'EdgeColor', DESIGN.sim_rgb, 'FaceColor', 'none', 'LineWidth',2);
h_sim.EdgeAlpha = .7;

% SMOOTH ------------------------------------------------------------------
% xy_pol = sortrows(indiM_smooth,1);
% [x,y] = pol2cart(deg2rad(xy_pol(:,1)),xy_pol(:,2));
% h_smooth = plot(af(x),af(y),'-', 'Color', 'k', 'LineWidth', 1);

hold off

% FORMAT ------------------------------------------------------------------
axis([-1 1 -1 1] * DESIGN.axmx);
axis square;
set(gca, 'FontSize', DESIGN.fontsize);
title('DATA','FontWeight','bold');
xlabel(DESIGN.xlbl);
ylabel(DESIGN.ylbl);

% HANDLES -----------------------------------------------------------------
h = [h_sim; h_data; h_Mdata];

%% Mdata_plotter
function h = Mdata_plotter(indi, indiM_HC, indiM_smooth, ConeCon_HC, cones, colspace, DESIGN)
% 2025.07.19 [cw]

hold on

axes_plotter(DESIGN.axmx+1);
cone_plotter(cones, DESIGN.axmx, '-', colspace);

% AVERAGE -----------------------------------------------------------------
xy_pol = sortrows(indiM_HC,1);
[x,y] = pol2cart(deg2rad(xy_pol(:,1)),xy_pol(:,2));
%        h_Mdata = plot(af(x),af(y),'-', 'Color',data_rgb, LineWidth',2);
h_Mdata = fill(x, y, [1 1 1]*.9, 'EdgeColor', [1 1 1]*.7, 'LineWidth', 0.5);

% SIMULATION OF CONE ADAPTATION -------------------------------------------
xy_pol = sortrows(ConeCon_HC,1);
[x,y] = pol2cart(deg2rad(xy_pol(:,1)),xy_pol(:,2));
h_sim = fill(x, y, 'w', 'EdgeColor', DESIGN.sim_rgb, 'FaceColor', 'none', 'LineWidth',2);

% SMOOTH ------------------------------------------------------------------
xy_pol = sortrows(indiM_smooth,1);
[x,y] = pol2cart(deg2rad(xy_pol(:,1)),xy_pol(:,2));
h_smooth = plot(af(x),af(y),'-', 'Color', 'k', 'LineWidth', 1);

% DATA --------------------------------------------------------------------
x = indi(:,:,2);
y = indi(:,:,3);
h_data = scatter(x(:),y(:), 1,'k','filled', 'MarkerFaceAlpha',.5);

hold off

% FORMAT ------------------------------------------------------------------
axis([-1 1 -1 1] * DESIGN.axmx);
axis square;
set(gca, 'FontSize', DESIGN.fontsize);

% HANDLES -----------------------------------------------------------------
h = [h_sim; h_data; h_Mdata; h_smooth];

%% mod_agg_plotter
function mod_agg_plotter(AGG, agg_ctab, DESIGN)
% 2022.09.13
% 2024.07.18 added colspace

switch lower(AGG.colspace)
    case 'dkl'
        DESIGN.radius = 1;
        DESIGN.axmx = 3;
        axmx2 = 2;
        DESIGN.xlbl = 'L-M';
        DESIGN.ylbl = 'S-(L+M)';
        exp_lbl = 'Exp 2b ';
    otherwise
        DESIGN.radius = 50;
        DESIGN.axmx = 100;
        axmx2 = 100;
        DESIGN.xlbl = 'Green-Red [u*]';
        DESIGN.ylbl = 'Blue-Yellow [v*]';
        exp_lbl = 'Exp 2a ';
end

% 1 AVERAGE DATA ----------------------------------------------------------
ax(1) = subplot(1,3,1);

% DIAGRAM 
L = AGG.stimuli{1}.bg.(AGG.colspace)(1);
colourcircler(AGG.colspace,L,DESIGN.radius,0:360,50,DESIGN.radius);
alldata = cat(2,AGG.Indi{:});
Magg_plotter(...
    alldata,...
    [AGG.agg.hue(:,1), AGG.agg.chroma(:,1)],...
    [AGG.agg.hue(:,1), AGG.agg.chroma_smooth(:,1)],...
    [AGG.agg.sim_hueM, AGG.agg.sim_chrM],...
    AGG.dir(1).cones,...
    AGG.colspace,...
    DESIGN);

% ANNOTATIONS:

% Experiment ID and Sample Size
sti_n = size(alldata,1);
rep_n = size(alldata,2);

switch lower(AGG.colspace)
    case 'dkl'
        text(1,1,sprintf('%2d meas \n%2d sti ', rep_n, sti_n),...
            'FontSize', DESIGN.fontsize,...
            'Units','Normalized',...
            'HorizontalAlignment', 'Right', 'VerticalAlignment','Top');
    otherwise
        text(1,1, exp_lbl,...
            'Units', 'Normalized',...
            'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top',...
            'Color', 'k',...
            'FontSize', DESIGN.fontsize, 'FontWeight','bold');
        text(1,1,sprintf('\n%2d meas / %2d sti ', rep_n, sti_n),...
            'FontSize', DESIGN.fontsize,...
            'Units','Normalized',...
            'HorizontalAlignment', 'Right', 'VerticalAlignment','Top');
end

% ANNOTATION CORRELATION -------------------------------------------------------------
text(0,0, sprintf(' r(%d)=%.2f%s\n smooth R^2=%d%%', agg_ctab.df(1), agg_ctab.r(1), agg_ctab.sig{1}, round(agg_ctab.smooth_R2(1))),...
    'Units', 'Normalized',...
    'Color', 'k', 'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');

% 2 DEVIATIONS FROM OPPONENT HUES ----------------------------------------- 

% DIAGRAM:
ax(2) = subplot(1,3,2);
ihue = (0:5:355)';
hold on

% REFRENCE LINES (GRID)
plot([1;1]*[90 180 270], [-1 1]*100, '-', 'Color', [1 1 1]*.8); % Opponent axes
plot([0 360], [0 0], '-', 'Color', [1 1 1]*.8, 'LineWidth', 1); % Zero deviations

xy  = af(sortrows([ihue, AGG.agg.sim_oppdev],1));
xy(end,1) = 360;
%data = af(sortrows([ihue, AGG.agg.oppdev(:,1), AGG.agg.oppdev_semM],1));
data = af(sortrows([ihue, AGG.agg.oppdev],1));
data(end,1) = 360;
xy3 = af(sortrows([ihue, AGG.agg.oppdev_smooth], 1));
xy3(end,1) = 360;
DESIGN.data_rgb = [1 1 1]*0.5;
h = std_areator(data(:,1), data(:,2), data(:,3), DESIGN.data_rgb, .2);
h_data   = plot(data(:,1), data(:,2),'-', 'Color', DESIGN.data_rgb, 'LineWidth', 0.5);
h_sim    = plot(xy(:,1),  xy(:,2), '-', 'Color', DESIGN.sim_rgb, 'LineWidth',1);
h_smooth = plot(xy3(:,1), xy3(:,2),'-', 'Color', 'k', 'LineWidth', 1);
hold off

% ANNOTATION (Correlation)
% text(0,0, sprintf(' r(%d)=%.2f%s', agg_ctab.df(2), agg_ctab.r(2), agg_ctab.sig{2}),...
%     'Units', 'Normalized',...
%     'Color', 'k', 'FontSize', DESIGN.fontsize,...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');
% text(1,0, sprintf(' smooth R^2=%d%% ', round(agg_ctab.smooth_R2(2))),...
%     'Units', 'Normalized',...
%     'Color', 'k', 'FontSize', DESIGN.fontsize,...
%     'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom');

text(0,0, sprintf(' r(%d)=%.2f%s\n sR^2=%d%%', agg_ctab.df(2), agg_ctab.r(2), agg_ctab.sig{2}, round(agg_ctab.smooth_R2(2))),...
    'Units', 'Normalized',...
    'Color', 'k', 'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');

% FORMATTING --------------------------------------------------------------
axis square;
set(gca, 'FontSize', DESIGN.fontsize,...
    'XLim', [0 360], 'XTick', 0:90:360,...
    'YLim', [-1 1]*22);
xlabel('Inducer Hue [deg]');
ylabel('Difference from Opponent'); 
title('DEVIATION FROM OPPONENCY','FontWeight','bold');

% 3 POLAR HUE DIAGRAM -----------------------------------------------------
subplot(1,3,3)

% DIAGRAMS
DESIGN.axmx = ceil(max(AGG.agg.rF));
colourcircler(AGG.colspace,L,DESIGN.radius,0:360,50, 1.4);
h = huefrq_plotter(...
    AGG.agg.rF, ...
    AGG.agg.rF_smooth, ...
    AGG.agg.sim_rF,...
    AGG.dir(1).cones,...
    AGG.colspace,...
    DESIGN);

% ANNOTATIONS ---------------------------------------------------------
text(0,0, sprintf(' r(%d)=%.2f%s\n sR^2=%d%%', agg_ctab.df(3), agg_ctab.r(3), agg_ctab.sig{3}, round(agg_ctab.smooth_R2(3))),...
    'Units', 'Normalized',...
    'Color', 'k', 'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');

% FORMATING -----------------------------------------------------------
xlabel(DESIGN.xlbl); 
ylabel(DESIGN.ylbl);
title('HUE HISTOGRAM','FontWeight','bold');

%% mod_indiM_plotter
function mod_indiM_plotter(AGG, ctab, DESIGN)
% 2025.07.16 adapted for publication

switch lower(AGG.colspace)
    case 'dkl'
        DESIGN.axmx = 2.5;
        DESIGN.chr_txt = ' Chroma\n %0.2f / %0.2f';
        xlbl = 'L-M';
        ylbl = 'S-(L+M)';
    otherwise
        DESIGN.axmx = 100;
        DESIGN.chr_txt = ' Chroma\n %2.0f / %2.0f';
        xlbl = 'u*';
        ylbl = 'v*';
end

pp_n = size(AGG.IndiM.repN,2);
col_n = ceil(pp_n/2);
row_n = 2;
for pp = 1:pp_n

    subplot(row_n, col_n, pp);

    % DIAGRAM -------------------------------------------------------------
    Mdata_plotter(...
        AGG.Indi{pp},...
        [AGG.IndiM.hueM(:,pp), AGG.IndiM.chrM(:,pp)],...
        [AGG.IndiM.hueM(:,pp), AGG.IndiM.chrM_smooth(:,pp)],...
        [AGG.model.Hue.ConeCon(:,pp), AGG.model.Chroma.ConeCon(:,pp)],...
        AGG.dir(pp).cones,...
        AGG.colspace,...
        DESIGN);

    % ANNOTATIONS ---------------------------------------------------------

    % CHROMA (Average / Inducer)
    text(0,1,sprintf(DESIGN.chr_txt, AGG.pps.chrMM(pp), AGG.pps.Chroma(pp)),...
        'FontSize', DESIGN.fontsize,...
        'Units','Normalized',...
        'HorizontalAlignment', 'Left', 'VerticalAlignment','Top');


    % SAMPLE SIZE
    sti_n = size(AGG.Indi{pp},1);
    rep_n = size(AGG.Indi{pp},2);
    text(1,1,sprintf('%2d meas \n%2d sti ', rep_n, sti_n),...
        'FontSize', DESIGN.fontsize,...
        'Units','Normalized',...
        'HorizontalAlignment', 'Right', 'VerticalAlignment','Top');

    % CORRELATION
    text(0,0, sprintf(' r(%d)=%.2f%s\n sR^2=%d%%', ctab.df(pp), ctab.r(pp), ctab.sig{pp}, round(ctab.smooth_R2(pp))),...
        'Units', 'Normalized',...
        'Color', 'k', 'FontSize', DESIGN.fontsize,...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');

    title(AGG.pps.pp_lbl{pp},'FontWeight','bold');
    if pp == pp_n
        xlabel(xlbl);
        ylabel(ylbl);
    end

end

%% mod_indiHueFreq_plotter
function mod_indiHueFreq_plotter(AGG, ctab, DESIGN)
% 2022.09.13
% 2024.07.23 added colspace input

switch lower(AGG.colspace)
    case 'dkl'
        DESIGN.chr_txt = ' Chroma \n %0.2f / %0.2f';
        xlbl = 'L-M';
        ylbl = 'S-(L+M)';
    otherwise
        DESIGN.chr_txt = ' Chroma \n %2.0f / %2.0f';
        xlbl = 'green-red [u*]';
        ylbl = 'blue-yellow [v*]';
end

pp_n = size(AGG.hueF.rF,2);
col_n = ceil(pp_n/2);
row_n = 2;
DESIGN.axmx = ceil(max(AGG.hueF.rF(:)));
for pp = 1:pp_n

    ax = subplot(row_n, col_n, pp);

    % DIAGRAMS ------------------------------------------------------------
    h = huefrq_plotter(...
        AGG.hueF.rF(:,pp), ...
        AGG.hueF.rF_smooth(:,pp), ...
        AGG.model.hueF.rF(:,pp),...
        AGG.dir(pp).cones,...
        AGG.colspace,...
        DESIGN);

    % ANNOTATIONS ---------------------------------------------------------

    % SAMPLE SIZE
    sti_n = size(AGG.Indi{pp},1);
    rep_n = size(AGG.Indi{pp},2);
    text(1,1,sprintf('%2d meas \n%2d sti ', rep_n, sti_n),...
        'FontSize', DESIGN.fontsize,...
        'Units','Normalized',...
        'HorizontalAlignment', 'Right', 'VerticalAlignment','Top');

    % CORRELATION
    text(0,0, sprintf(' r(%d)=%.2f%s\n sR^2=%d%%', ctab.df(pp), ctab.r(pp), ctab.sig{pp}, round(ctab.smooth_R2(pp))),...
        'Units', 'Normalized',...
        'Color', 'k', 'FontSize', DESIGN.fontsize,...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');
    
    % FORMATING -----------------------------------------------------------
    if pp <= col_n
        ax.XTickLabels = [];
    end
    if pp == ceil(1.5*col_n) % Bottom Centre
        xlabel(xlbl);
    elseif mod(pp, col_n) == 1 % First panel of row
        ylabel(ylbl);
    end    
    title(AGG.pps.pp_lbl{pp},'FontWeight','bold');
end    

%% mod_indiOppoDev_plotter
function mod_indiOppoDev_plotter(AGG, ctab, DESIGN)
% 2025.07.19 polished for publication [cw]

switch lower(AGG.colspace)
    case 'dkl'
        DESIGN.axmx = 50;
        xlbl = 'L-M';
        ylbl = 'S-(L+M)';
        DESIGN.chr_txt = ' Chroma\n %0.2f / %0.2f ';
    otherwise
        DESIGN.axmx = 50;
        xlbl = 'u*';
        ylbl = 'v*';
        DESIGN.chr_txt = ' Chroma\n %2.0f / %2.0f ';
end

ihue = (0:5:355)';
pp_n = size(AGG.hueF.rF,2);
col_n = ceil(pp_n/2);
row_n = 2;
for pp = 1:pp_n
    
    ax = subplot(row_n, col_n, pp);

    hold on
    % REFRENCE LINES (GRID)
    plot([1;1]*[90 180 270], [-1 1]*100, '-', 'Color', [1 1 1]*.8); % Opponent axes
    plot([0 360], [0 0], '-', 'Color', [1 1 1]*.8, 'LineWidth', 1); % Zero deviations

    if ~isnan(AGG.hueF.rF(:,pp))             
        xy  = sortrows([ihue, AGG.IndiOppoDev.sim(:,pp)],1);
        xy2 = sortrows([ihue, AGG.IndiOppoDev.adj(:,pp), AGG.IndiOppoDev.adj_sem(:,pp)],1);
        xy3 = sortrows([ihue, AGG.IndiOppoDev.adj_smooth(:,pp)], 1);
        DESIGN.data_rgb = [1 1 1]*0.5;
        h_data   = plot(xy2(:,1), xy2(:,2),'-', 'Color', DESIGN.data_rgb, 'LineWidth', 0.5);
        h_sim    = plot(xy(:,1),  xy(:,2), '-', 'Color', DESIGN.sim_rgb, 'LineWidth',1);
        h_smooth = plot(xy3(:,1), xy3(:,2),'-', 'Color', 'k', 'LineWidth', 1);
    end
    hold off
    axis square;
    ylim([-1 1]*DESIGN.axmx);
    set(gca, 'FontSize', DESIGN.fontsize,...
        'XLim', [0 360], 'XTick', 0:90:360);

    % ANNOTATIONS ---------------------------------------------------------
    
    % CORRELATION
    text(0,0, sprintf(' r(%d)=%.2f%s\n sR^2=%d%%', ctab.df(pp), ctab.r(pp), ctab.sig{pp}, round(ctab.smooth_R2(pp))),...
        'Units', 'Normalized',...
        'Color', 'k', 'FontSize', DESIGN.fontsize,...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');
        
    % FORMAT --------------------------------------------------------------
    if pp <= col_n
        ax.XTickLabels = [];
    end
    if pp == ceil(1.5*col_n) % Bottom Centre
        xlabel('Inducer Hue [deg]');
    elseif mod(pp, col_n) == 1 % First panel of row
        ylabel('Difference from Opponent');
    end
    title(AGG.pps.pp_lbl{pp},'FontWeight','bold');
end    

%% pp_reporter
function RES = pp_reporter(RES, pps)
% 2025.07.19 

header('PARTICIPANTS',4);
N = size(pps,1);
sex     = tabulate(pps.sex); 
women   = sex{1,2};
men     = sex{2,2};
age     = [nanmean(pps.age), nanstd(pps.age)];
fprintf('N = %d participants, \n', N);
fprintf('SEX: %d women + %d men\n', women, men);
fprintf('AGE: M = %.2f+ SD = %.2f\n', age(1), age(2));
RES.pps = table(N, women, men, age);

%% ***************************** SUBFUNCTIONS *****************************

%% af
function y = af(x, dm)
% 2022.09.13 [cw]

if nargin < 2
    if size(x, 1) == 1 % row vector
        dm = 2;
    elseif size(x, 2) == 1 % column vector
        dm = 1;
    else 
        dm = 1;
    end
end    
if dm == 1 % append row
  y = [x; x(1,:)];
elseif dm == 2 % append column
  y = [x x(:,1)];
elseif dm == 3
  y = cat(3,x, x(:,:,1)); 
  %
end

%% axes_plotter
function axes_plotter(mx, rgb, linestyle)
% 2022.09.13

if nargin < 3
    linestyle = '-';
    if nargin < 2
        rgb = [.8 .8 .8];
    end
end

% AXIS:
plot([-1 1; 0 0]'*mx, [0 0; -1 1]'*mx, linestyle, 'Color', rgb);

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

%% circular_remapper
function data = circular_remapper(data, ref, onecycle)
%2018.08.18 

if nargin < 3
    onecycle = 360;
end

halfcycle = onecycle/2;

inds = data-ref < -halfcycle;
data(inds) = data(inds) + onecycle;

inds = data-ref > halfcycle;
data(inds) = data(inds) - onecycle;

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

%% cone_plotter
function cone_plotter(CONES, mx, linestyle, colspace, linewidth)
% 2022.09.13
% 2024.07.23 added colspace
% 2025.07.28 added linewidht

if nargin < 5
    linewidth = 1;
    if nargin < 4
        colspace = 'luv';
        if nargin < 3
            linestyle = '-';
        end
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
    linestyle, 'Color', [1 0 0], 'LineWidth', linewidth);
plot([0, cones(2,2)],[0, cones(2,3)],...
    linestyle, 'Color', [0 0.8 0], 'LineWidth', linewidth);
plot([0, cones(3,2)],[0, cones(3,3)],...
    linestyle, 'Color', [0 0 1], 'LineWidth', linewidth);

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

%% colourcircler
function colourcircler(space, lum, rad, azi, markersize, rad_scale_to, varargin)
% 2020.06.04 based on colorcircler, but added rescaling input "rad_scale_to".

if nargin < 6
    rad_scale_to = rad; % same as colour radius
    if nargin < 5
        markersize = 10;
        if nargin < 4
            azi = (0:1:360)';
            if nargin < 3
                rad = 50;
                if nargin < 2
                    lum = 60;
                    if nargin < 1
                        space = 'Luv';
                    end
                end
            end
        end
    end
end

azi = deg2rad(azi(:));
if numel(rad) == 1 
    rad = ones(size(azi,1),1)*rad;
end

if numel(lum) == 1
    lum = ones(size(azi,1),1)*lum;
end

[x, y] = pol2cart(azi, rad);
switch space
    case 'hsv'
        col = [azi/2/pi,rad,lum]; 
    otherwise
        col = [lum,x,y];
end
COL = colourconverter(col,space);
rgbs = COL.rgb/255;

h = 0;
if ~ishold
    h = 1;
    hold on
end

% RESCALE -----------------------------------------------------------------
x2 = x./rad*rad_scale_to;
y2 = y./rad*rad_scale_to;
% PLOT --------------------------------------------------------------------
scatter(x2,y2,markersize,rgbs, 'filled', 'MarkerEdgeColor', 'none', varargin{:});

if h
    hold off
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

%% sem_areator
function h = std_areator(x, y, sd, rgb, alpha)
% 2025.07.23 [cw]

if nargin < 5
    alpha = .5;
    if nargin < 4
        rgb = 'k';
    end
end

x = x(:); y = y(:);  % column vectors
upper = y + sd;
lower = y - sd;
x2 = [x; flipud(x)];
sd_shade = [upper; flipud(lower)];
h = fill(x2, sd_shade, rgb, 'FaceAlpha', alpha, 'EdgeColor','none');

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
%2021.07.01  [cw]

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
                        statstab = table(1, 100, {'***'}, 'VarNames', {'r', 'df','sig'});
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
