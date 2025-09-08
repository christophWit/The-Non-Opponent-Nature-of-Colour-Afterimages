function RES = mod_exp1_analysor(AGG1a, AGG1b, which_one)
%2025.07.23 adapted + polished for publication [cw]

%% ******************************* SETTINGS *******************************

% INITIALISE OUTPUT STRUCTURE WITH RESULTS
RES = struct; 

%% DESIGN (cosmetics)
DESIGN.data_rgb     = [1 1 1]*0.5;
DESIGN.sim_rgb      = [1 0 0];
DESIGN.fontsize     = 9;
DESIGN.fontsize2    = 10;
DESIGN.font         = 'Arial';
DESIGN.markersize   = 5;

%% ******************************** MAIN **********************************

%% DIFFERENCES FROM OPPONENT
header('Differences from Opponent', 3);
header('Exp 1a');
RES.ttab_exp1a = ttester(AGG1a.prederr.hue.dkl);
RES.ttab_exp1a.Properties.RowNames = AGG1a.stimuli.lbl;
disp(RES.ttab_exp1a);

header('Exp 1b');
RES.ttab_exp1b = ttester(AGG1b.prederr.hue.dkl);
RES.ttab_exp1b.Properties.RowNames = AGG1b.stimuli.lbl;
disp(RES.ttab_exp1b);

%% MAIN
header('Corrleations for Aggregated Data (Exp 1b)', 3);
RES.agg_ctab(1,:)        = correlator(AGG1a.OppoDev.sim, AGG1a.OppoDev.adj(:,1));
RES.agg_ctab(2,:)        = correlator(AGG1b.OppoDev.sim, AGG1b.OppoDev.adj(:,1));
RES.agg_ctab(3,:)        = correlator(AGG1b.model.hueF.rF, AGG1b.afteri.hue_frq.rF);
RES.agg_ctab.Properties.RowNames = {'OppoDev1a', 'OppoDev1b', 'HueF1b'};

disp(RES.agg_ctab);

switch lower(which_one)
    case {'main'}

        % FIGURE ----------------------------------------------------------
        figure('Name', 'Agg-Main', 'NumberTitle','off', 'InvertHardcopy', 'off', 'Color', 'w');
        DESIGN.examples     = [2 4 6 7]; % Example data shown in Figure 2.b
        mod_main_plotter(AGG1a, AGG1b, RES.agg_ctab(2:3,:), DESIGN);

end

%% SUPPLEMENTARY
switch lower(which_one)
    case {'supplementary'}
        figure('Name', 'Supplementary', 'NumberTitle','off', 'InvertHardcopy', 'off', 'Color', 'w');
        DESIGN.examples     = 1:8; % All
        mod_supplementary_plotter(AGG1a, DESIGN);
end

%% HISTOGRAMS
switch lower(which_one)
    case {'histo', 'histograms', 'all'}
        figure('Name', 'histo', 'NumberTitle', 'off', 'InvertHardCopy', 'off', 'Color', 'w');
        histo_plotter(AGG1a, DESIGN);
end


%% ************************* SPECIFIC SUBFUNCTIONS ************************

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

%% histo_plotter
function histo_plotter(AGG, DESIGN)
% 2025.09.07 polished for publication [cw]

Inducer = AGG.stimuli.inducer;
relFreq = AGG.afteri.resp_freq.rF;
COMP    = AGG.stimuli.comp;
CompRGB = AGG.stimuli.comp_rgb;
OPPO    = AGG.model.Hue; 
TypHue  = AGG.typi.agg.hue_M;
%TypRGB  = AGG.typi.agg.rgb;
TypRGB  = [...
    0  .3 1;... blue
    0  .3 1;... blue
    .7 0 1;... purple
    1 .5 .7;... pink
    1 .5 .7;... pink
    .5 .3 .2;... brown
    0 .7 0;... green
    0 .7 0]; % green

mSize = DESIGN.markersize*5;
N = size(AGG.afteri.N,1);

sti_n = 8;
for col = 1:sti_n

    % MAIN DATA -----------------------------------------------------------
    i_lbl = AGG.stimuli.lbl{col};
    i_hue = Inducer.Luv_pol(col,1);
    i_rgb = Inducer.rgb(col,:);
    c_hue = COMP(col,:,2)';
    c_hue = cdissolver(mod(c_hue,360));
    c_rgb = shiftdim(CompRGB(col,:,:),1);
    rF  = relFreq(:,col)*100;
    xlm = [(floor(min(c_hue)-5)/10)*10, (ceil(max(c_hue)+5)/10)*10];

    % REFERENCE HUES (MODELS + AVERAGE) FOR COMPARISON -------------------- 

    % AVERAGE
    M   = AGG.afteri.agg.hue_M(col);
    sem = AGG.afteri.agg.hue_sem(col);

    % MODELS OF COMPLEMENTARITY
    oppo_cone = OPPO.ConeCon(col);
    oppo_lab = OPPO.Lab(col);
    oppo_mun = OPPO.mun(col);
    typi = TypHue(col);
    trgb = TypRGB(col,:);

    % MAPPING THOSE VALUES TO THE RANGE OF THE COMPARISONS/BARS ----------- 
    M           = circular_remapper(M, c_hue(1));
    oppo_lab    = circular_remapper(oppo_lab, c_hue(1));
    oppo_mun    = circular_remapper(oppo_mun, c_hue(1));
    typi        = circular_remapper(typi, c_hue(1));
    
    subplot(2,4,col)        
    hold on        

    % LINES FOR CIELAB & MUNSELL PREDICTIONS ------------------------------
    line([1 1] * oppo_lab, [0 60], 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth',1);
    line([1 1] * oppo_mun, [0 60], 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth',1);
        
    % BARS ----------------------------------------------------------------
    for k=1:9
        bar(c_hue(k),rF(k),'FaceColor', c_rgb(k,:), 'BarWidth', 10, 'LineWidth', 0.5);
    end
    bar(c_hue(5),rF(5),'FaceColor', c_rgb(5,:), 'BarWidth', 10, 'LineWidth', 1.5);
            
    % AVERAGE DATA --------------------------------------------------------
    x = [M-sem, M-sem, M+sem,M+sem];
    y = [0, 100,100, 0];
    h_area = area(x,y);
    set(h_area, 'FaceColor', [0.5 0.5 0.5],'EdgeColor', 'none');
    set(h_area, 'FaceAlpha', 0.5);
    line([1 1] * M, [0 100], 'Color', [0 0 0], 'LineStyle', '-', 'LineWidth',1);
    
    % PREDICTION BY CONE ADAPTATION ---------------------------------------
    line([1 1]*oppo_cone, [0 90],  'LineStyle', ':', 'Color', [1 0 0], 'LineWidth', 2);

    % COLOUR OF INDUCER ---------------------------------------------------
    plot(max(c_hue), 90, 'ko','MarkerFaceColor', i_rgb, 'MarkerSize', DESIGN.markersize); 
    
    % TYPICAL COLOURS -----------------------------------------------------
    scatter(typi, 0, mSize, trgb, 'filled', 'Marker', '^', 'MarkerEdgeColor', 'k');
    
    % TEXT LABELS FOR LAB & MUNSELL ---------------------------------------
    % Only if the are within x-axis limits
    if oppo_lab >= xlm(1) & oppo_lab <= xlm(2)
        text(oppo_lab, 60, 'Lab',...
            'FontSize', DESIGN.fontsize,...
            'Rotation', 90,...
            'BackgroundColor', 'w', 'Margin', 0.1,'LineStyle', 'none',...
            'HorizontalAlignment', 'left', 'VerticalAlignment','Middle');
    end
    if oppo_mun >= xlm(1) & oppo_mun <= xlm(2)
        text(oppo_mun, 60, 'Mun',...
            'FontSize', DESIGN.fontsize,...
            'Rotation', 90,...
            'BackgroundColor', 'None', 'Margin', 0.1,'LineStyle', 'none',...
            'HorizontalAlignment', 'left', 'VerticalAlignment','Middle');
    end

    hold off
    
    % FORMAT --------------------------------------------------------------
    title(i_lbl, 'FontWeight','bold');
    if col == 1
        text(0,1, sprintf('  N = %d', N),...
            'FontSize', DESIGN.fontsize,...
            'Units', 'Normalized',...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'Top');
    end
    if mod(col,4) == 1
        ylabel('Relative frequency [%]');
    elseif col >= 4
        xlabel('CIELUV Hue [deg]');        
    end
    xlim(xlm);
    ylim([0 100]);
    set(gca, 'XTick', -400:20:400, 'YTick', 0:20:100);
    axis square;
end

%% huefrq_plotter
function h = huefrq_plotter(ihue, rF, rF_smooth, rf_sim, cones, colspace, DESIGN)
% 2025.07.19 [cw]

h = [];
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
if  ~isempty(rF_smooth)
    [x,y] = pol2cart(ihue_rad, rF_smooth);
    h_smooth = fill(x,y,[.8 .8 .8], 'EdgeColor', 'k', 'FaceColor', 'none', 'LineWidth',1);
else 
    h_smooth = [];
end
hold off
axis([-1 1 -1 1] * DESIGN.axmx);
axis square;
set(gca, 'FontSize', DESIGN.fontsize);

% HANDLES -----------------------------------------------------------------
h = [h_sim; h_data; h_smooth];

%% mod_main_plotter
function mod_main_plotter(AGG1a, AGG1b, agg_ctab, DESIGN)
% 2025.07.23 polished for publication [cw]

DESIGN.radius = 50;


% 1 EXAMPLE DATA ----------------------------------------------------------

ax(1) = subplot(1,3,1);
polar_plotter(AGG1a, DESIGN);
xlabel(''); 
ylabel('Blue-Yellow [v*]');
axis([-150 150 -150 150]);
title('AVERAGE DATA', 'FontSize', DESIGN.fontsize, 'FontWeight', 'bold');

% ANNOTATION:
pp_n    = size(AGG1a.afteri.indi,1);
txt = sprintf('\n%d pp ', pp_n);
text(ax(1),1,1, txt,...
    'Units', 'Normalized',...
    'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top',...
    'Color', 'k', 'BackgroundColor', get(ax(1),'color'), 'Margin', .001,...
    'FontSize', DESIGN.fontsize);
text(ax(1),1,1, 'Exp 1a ',...
    'Units', 'Normalized',...
    'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top',...
    'Color', 'k', 'BackgroundColor', get(ax(1),'color'), 'Margin', .001,...
    'FontSize', DESIGN.fontsize, 'FontWeight','bold');


% 2 DEVIATIONS FROM OPPONENT HUES ----------------------------------------- 

% DIAGRAM:
ax(2)   = subplot(1,3,2);
ihue    = AGG1b.stimuli.inducer.Luv_pol(:,1);
pe      = AGG1b.prederr.hueM.dkl(:,1);
pe_sd   = AGG1b.prederr.hueM.dkl(:,2);
data = af(sortrows([ihue, pe, pe_sd],1));
data(end,1) = 360;

sim     = AGG1b.OppoDev.sim;
sim  = af(sortrows([ihue, sim],1));
sim     = [AGG1b.model360.inducer(:,1), AGG1b.model360.OppoDev];
sim  = af(sortrows(sim,1));
sim(end,1) = 360;

smooth  = AGG1b.OppoDev.adj_smooth;

hold on

% REFRENCE LINES (GRID)
plot([1;1]*[90 180 270], [-1 1]*100, '-', 'Color', [1 1 1]*.8); % Opponent axes
plot([0 360], [0 0], '-', 'Color', [1 1 1]*.8, 'LineWidth', 1); % Zero deviations

std_areator(data(:,1), data(:,2), data(:,3), DESIGN.data_rgb, .1);
DESIGN.data_rgb = 'k';
h_sim    = plot(sim(:,1),  sim(:,2), '-', 'Color', DESIGN.sim_rgb, 'LineWidth',1);
h_data   = plot(data(:,1), data(:,2),'.-', 'Color', DESIGN.data_rgb, 'LineWidth', 1);
hold off

% ANNOTATIONS:

% Experiment ID and Sample Size
pp_n    = size(AGG1b.afteri.indi,1);
sti_n   = size(AGG1b.afteri.indi,2);
txt = sprintf('\n%d pp / %d sti ', pp_n, sti_n);
text(ax(2),1,1, txt,...
    'Units', 'Normalized',...
    'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top',...
    'Color', 'k', 'BackgroundColor', get(ax(2),'color'), 'Margin', .001,...
    'FontSize', DESIGN.fontsize);
text(ax(2),1,1, 'Exp 1b ',...
    'Units', 'Normalized',...
    'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top',...
    'Color', 'k', 'BackgroundColor', get(ax(2),'color'), 'Margin', .001,...
    'FontSize', DESIGN.fontsize, 'FontWeight','bold');

% Correlation:
text(0,0, sprintf(' r(%d)=%.2f%s', agg_ctab.df(1), agg_ctab.r(1), agg_ctab.sig{1}),...
    'Units', 'Normalized',...
    'Color', 'k', 'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');

% FORMATTING --------------------------------------------------------------
axis square;
set(gca, 'FontSize', DESIGN.fontsize,...
    'XLim', [0 360], 'XTick', 0:90:360,...
    'YLim', [-1 1]*22);
xlabel('Difference from Opponent'); 
ylabel('Inducer Hue [deg]');
title('DEVIATION FROM OPPONENCY', 'FontSize', DESIGN.fontsize, 'FontWeight','bold');

% 3 POLAR HUE DIAGRAM -----------------------------------------------------
ax(3) = subplot(1,3,3);

% DIAGRAMS
DESIGN.axmx = ceil(max(AGG1b.afteri.hue_frq.rF));
L = AGG1b.stimuli.bg.Luv(1);
colourcircler('Luv',L,DESIGN.radius,0:360,50, 1.4);
h = huefrq_plotter(...
    AGG1b.afteri.agg.inducer(:,2),...
    AGG1b.afteri.hue_frq.rF, ...
    [], ...
    AGG1b.model.hueF.rF,...
    AGG1b.dir.cones,...
    'Luv',...
    DESIGN);

% ANNOTATIONS:

% Correlation:
text(0,0, sprintf(' r(%d)=%.2f%s', agg_ctab.df(2), agg_ctab.r(2), agg_ctab.sig{2}),...
    'Units', 'Normalized',...
    'Color', 'k', 'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom');

% FORMATING -----------------------------------------------------------
xlabel('Green-Red [u*]'); 
ylabel('Blue-Yellow [v*]');
title('HUE HISTOGRAM', 'FontSize',DESIGN.fontsize, 'FontWeight','bold');

%% mod_supplementary_plotter
function mod_supplementary_plotter(AGG1a, DESIGN)
% 2025.09.07 polished for publication [cw]

% 1 AFTERIMAGES -----------------------------------------------------------
ax(1) = subplot(1,2,1);
polar_plotter_simple(AGG1a, DESIGN);
xlabel('Gree-Red [u*]'); 
ylabel('Blue-Yellow [v*]');
title('AFTERIMAGE MATCHES', 'FontSize', DESIGN.fontsize, 'FontWeight', 'bold');

% 2 DISCRIMINATION (Non-Illusory Matches) ---------------------------------
ax(2) = subplot(1,2,2);
polar_plotter4discri(AGG1a, DESIGN);
xlabel('Gree-Red [u*]'); 
ylabel('');
axis([-1 1 -1 1]*60);
title('NON-ILLUSORY MATCHES', 'FontSize', DESIGN.fontsize, 'FontWeight', 'bold');

%% polar_plotter
function polar_plotter(AGG, DESIGN)
% 2025.07.23 [cw]

% SELECT EXAMPLE DATA -----------------------------------------------------
rgb = AGG.stimuli.inducer.rgb(DESIGN.examples,:);
inducer = AGG.stimuli.inducer.Luv(DESIGN.examples,2:3);
[oppo(:,1),oppo(:,2)] = pol2cart(deg2rad(AGG.model.Hue.dkl(DESIGN.examples)), 300);
[cone(:,1),cone(:,2)] = pol2cart(deg2rad(AGG.model.Hue.ConeCon(DESIGN.examples)), 300);
results = AGG.afteri.agg.Luv(DESIGN.examples,2:3);

% PLOT --------------------------------------------------------------------
hold on
sti_n = numel(DESIGN.examples);
for sti = 1:sti_n
    % INDUCERS
    plot([0,inducer(sti,1)],[0,inducer(sti,2)], 'k-', 'LineWidth', 2);
    plot([0,inducer(sti,1)],[0,inducer(sti,2)], '-',  'LineWidth', 1, 'Color', rgb(sti,:));

    % CONE-OPPONENT PREDICTIONS
    plot([0 oppo(sti,1)], [0 oppo(sti,2)], '-', 'Color', rgb(sti,:), 'LineWidth', .5);

    % CONE ADAPATION PREDICTIONS
    plot([0 cone(sti,1)]', [0 cone(sti,2)]', 'k-', 'LineWidth', 1);
    plot([0 cone(sti,1)]', [0 cone(sti,2)]', '--', 'Color', rgb(sti,:), 'LineWidth', 1);

    % RESULTS
    plot([0,results(sti,1)],[0,results(sti,2)],'k-', 'LineWidth', 2);
    plot([0,results(sti,1)],[0,results(sti,2)],'-', 'LineWidth', 1, 'Color', rgb(sti,:));

end

% SYMBOLS ON TOP OF EVERYTHING
scatter(inducer(:,1), inducer(:,2), 40, rgb, 'filled', 'o', 'MarkerEdgeColor', 'k');
scatter(results(:,1), results(:,2), 80, rgb, 'filled', 'p', 'MarkerEdgeColor', 'k');

hold off

% FORMATTING --------------------------------------------------------------
set(gca,...
    'XTick', -300:50:300, 'YTick', -300:50:300,...
    'FontName', DESIGN.font, 'FontSize', DESIGN.fontsize,...
    'Color', AGG.stimuli.bg.rgb);
axis([-150 250 -200 200]);
axis square

%% polar_plotter_simple
function polar_plotter_simple(AGG, DESIGN)
% 2025.07.23 [cw]

% SELECT EXAMPLE DATA -----------------------------------------------------
rgb = AGG.stimuli.inducer.rgb(DESIGN.examples,:);
inducer = AGG.stimuli.inducer.Luv(DESIGN.examples,2:3);
[oppo(:,1),oppo(:,2)] = pol2cart(deg2rad(AGG.model.Hue.dkl(DESIGN.examples)), 300);
[cone(:,1),cone(:,2)] = pol2cart(deg2rad(AGG.model.Hue.ConeCon(DESIGN.examples)), 300);
results = AGG.afteri.agg.Luv(DESIGN.examples,2:3);

% PLOT --------------------------------------------------------------------
hold on
sti_n = numel(DESIGN.examples);
for sti = 1:sti_n
    % INDUCERS
    % plot([0,inducer(sti,1)],[0,inducer(sti,2)], 'k-', 'LineWidth', 2);
    % plot([0,inducer(sti,1)],[0,inducer(sti,2)], '-',  'LineWidth', 1, 'Color', rgb(sti,:));

    % CONE-OPPONENT PREDICTIONS
    plot([0 oppo(sti,1)], [0 oppo(sti,2)], '-', 'Color', rgb(sti,:), 'LineWidth', .5);

    % CONE ADAPATION PREDICTIONS
    plot([0 cone(sti,1)]', [0 cone(sti,2)]', 'k-', 'LineWidth', 1);
    plot([0 cone(sti,1)]', [0 cone(sti,2)]', '--', 'Color', rgb(sti,:), 'LineWidth', 1);

    % RESULTS
    plot([0,results(sti,1)],[0,results(sti,2)],'k-', 'LineWidth', 2);
    plot([0,results(sti,1)],[0,results(sti,2)],'-', 'LineWidth', 1, 'Color', rgb(sti,:));

end

% SYMBOLS ON TOP OF EVERYTHING
scatter(inducer(:,1), inducer(:,2), 40, rgb, 'filled', 'o', 'MarkerEdgeColor', 'k');
scatter(results(:,1), results(:,2), 80, rgb, 'filled', 'p', 'MarkerEdgeColor', 'k');

hold off

% FORMATTING --------------------------------------------------------------
set(gca,...
    'XTick', -300:50:300, 'YTick', -300:50:300,...
    'FontName', DESIGN.font, 'FontSize', DESIGN.fontsize,...
    'Color', AGG.stimuli.bg.rgb);
axis([-150 250 -200 200]);
axis square

%% polar_plotter4discri
function polar_plotter4discri(AGG, DESIGN)
% This is the same as polar_plotter, but uses (non-illusory) discrimination
% data and does not show the prediction by cone adaptation to avoid
% clutter. 
% 2025.09.07 [cw]

% SELECT EXAMPLE DATA -----------------------------------------------------
rgb = AGG.stimuli.discri_target.rgb(DESIGN.examples,:);
target = AGG.stimuli.discri_target.Luv(DESIGN.examples,2:3);
[cone(:,1),cone(:,2)] = pol2cart(deg2rad(AGG.model.Hue.ConeCon(DESIGN.examples)), 300);
[oppo(:,1),oppo(:,2)] = pol2cart(deg2rad(AGG.model.Hue.dkl(DESIGN.examples)), 300);
results = AGG.discri.agg.Luv(DESIGN.examples,2:3);

% PLOT --------------------------------------------------------------------
hold on
sti_n = numel(DESIGN.examples);
for sti = 1:sti_n
    % TARGET
    % plot([0,target(sti,1)],[0,target(sti,2)], 'k-', 'LineWidth', 2);
    % plot([0,target(sti,1)],[0,target(sti,2)], '-',  'LineWidth', 1, 'Color', rgb(sti,:));

    % CONE-OPPONENT PREDICTIONS
    plot([0 oppo(sti,1)], [0 oppo(sti,2)], '-', 'Color', rgb(sti,:), 'LineWidth', .5);

    % CONE ADAPATION PREDICTIONS
    plot([0 cone(sti,1)]', [0 cone(sti,2)]', 'k-', 'LineWidth', 1);
    plot([0 cone(sti,1)]', [0 cone(sti,2)]', '--', 'Color', rgb(sti,:), 'LineWidth', 1);

    % RESULTS
    plot([0,results(sti,1)],[0,results(sti,2)],'k-', 'LineWidth', 2);
    plot([0,results(sti,1)],[0,results(sti,2)],'-', 'LineWidth', 1, 'Color', rgb(sti,:));

end

% SYMBOLS ON TOP OF EVERYTHING
scatter(target(:,1), target(:,2), 40, rgb, 'filled', 'o', 'MarkerEdgeColor', 'k');
scatter(results(:,1), results(:,2), 80, rgb, 'filled', 'p', 'MarkerEdgeColor', 'k');

hold off

% FORMATTING --------------------------------------------------------------
set(gca,...
    'FontName', DESIGN.font, 'FontSize', DESIGN.fontsize,...
    'XTick', -100:20:100, 'YTick', -100:20:100,...
    'Color', AGG.stimuli.bg.rgb);
axis square

%% ************************* GENERAL SUBFUNCTIONS *************************

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

%% ste
function y = ste(varargin)
% 2019.04.12 [cw]
N = sum(~isnan(varargin{1}),1);
y = sqrt(nanvar(varargin{:}))./sqrt(N);

%% ttester
function [tab] = ttester(data, ns, onetailor)
% 2020.07.22 [cw]

if nargin < 3
    onetailor = 0;
    if nargin < 2
        ns = 'ns';
    end
end

M = nanmean(data)';
sem = ste(data)';
sd = nanstd(data)';
d = M./sd; % Cohen's d
[h,p,ci,stats] = ttest(data);

df = stats.df';
t = stats.tstat';
ci = ci';
p = p';
sig = significancer(p,ns,onetailor);

tab = table(M, sem, d, t, df, ci, p, sig); 