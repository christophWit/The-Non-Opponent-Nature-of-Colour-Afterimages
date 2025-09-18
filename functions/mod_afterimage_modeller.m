function [SIM, fh] = mod_afterimage_modeller(MunChips, colspace, STI, HeringRYGB, STI1a)
%2025.07.29 polished for pub [cw]

if nargin < 5
    STI1a = exp1a_stimulus_maker;
    if nargin < 4
        HeringRYGB = [10.6014492753623;72.1533333333333;125.133333333333;227.94];
        if nargin < 3
            if nargin < 2
                colspace = 'Luv';
            end
            STI = stimulusstruc_maker(colspace);
        end
    end
end

header;
header(mfilename);

%% DEFAULT SETTINGS
SETTINGS.sim_n      = 1000; % Number of repetitions in simulation, set 100000 for final
DESIGN     = struct(...
    'addplots',     {{...
    'origin', 'bg', 'day',... polar, powerpolar
    'smooth', 'dkl','luv','lab','munsell', 'hsv', 'typi', 'legend', 'comp', 'adapt', 'conesim',... histo
    'axes', 'cones', 'huecircle'}},... 
    'sim_rgb',  'r',...
    'data_rgb', 'k',...
    'data_alpha', .1, ...
    'font', 'Arial',...
    'fontsize',  9,...
    'fontsize2',  10,...
    'linewidth', 1.5,...
    'fatlinewidth', 3,...
    'markersize', 5);


%% COLSPACE
switch lower(colspace)
    case 'luv'
        SIM.polarspace = 'luv_pol';
        SIM.adapt_chroma    = 30; 
        SETTINGS.sim_hueSD  = 11.3585; % Average across 10 pps in Chaser Experiment
        SETTINGS.sim_chrSD  = 7.1029;  % Average across 10 pps in Chaser Experiment
        DESIGN.axmax = 55;
    case 'dkl'
        SIM.polarspace = 'dkl_pol';
        SIM.adapt_chroma    = 0.6;
        SETTINGS.sim_hueSD  = 6.5420; % Average across 2 pps in DKL-Chaser Experiment
        SETTINGS.sim_chrSD  = 0.1144; % Average across 2 pps in DKL-Chaser Experiment
        DESIGN.axmax = 2;
end

%% ******************************* MODELS *********************************

%% OPPONENT/COMPLEMENTARY COLOURS
[SIM.Oppo.Hue, SIM.Oppo.Chroma, SIM.Oppo.Lum_ctrl, SIM.dkl_k, SIM.Oppo.Munsell] = complementary_calculator(STI, colspace, MunChips, HeringRYGB(:,1), SIM.adapt_chroma);

% NOISE
[SIM.Noisy.ConeCon(:,:,1), SIM.Noisy.ConeCon(:,:,2)] = noise_maker([SIM.Oppo.Hue.ConeCon, SIM.Oppo.Chroma.ConeCon], SETTINGS.sim_hueSD, SETTINGS.sim_chrSD, SETTINGS.sim_n);

%% AXES

% CONES -------------------------------------------------------------------
LMS = [...
    1 0 0;...
    0 1 0;...
    1/2 1/2 1];
%XYZ = lms2XYZ(LMS,'ss','1931');
XYZ = lms2XYZ(LMS, 'smithpokorny', 'judd');
%SIM.dir.cones = colourconverter(XYZ, 'XYZ',2, STI.wp.xyY, STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max, 1, 0);
SIM.dir.cones = colourconverter(XYZ, 'XYZ',2,...
    STI.wp.xyY, STI.sens, STI.cmf,...
    STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max,...
    STI.bg.xyY(3), 0);

% DKL ---------------------------------------------------------------------
dkl = [0 -1 0; 0 1 0; 0 0 -1; 0 0 1];
%SIM.dir.dkl = colourconverter(dkl, 'dkl', 2, STI.wp.xyY, STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max, 1, 0);
SIM.dir.dkl = colourconverter(dkl, 'dkl', 2,...
    STI.wp.xyY, STI.sens, STI.cmf,...
    STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max,...
    STI.bg.xyY(3), 0);

% HERING COLOURS (oriignally from Experiment 1b) --------------------------
R = [70, HeringRYGB(1), 100];
Y = [70, HeringRYGB(2), 100];
G = [70, HeringRYGB(3), 100];
B = [70, HeringRYGB(4), 100];
%SIM.dir.HeringRYGB = colourconverter([R; Y; G; B], 'Luv_pol', 2, STI.wp.xyY, STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max,1,0);
SIM.dir.HeringRYGB = colourconverter([R; Y; G; B], 'Luv_pol', 2,...
    STI.wp.xyY, STI.sens, STI.cmf,...
    STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max,...
    STI.bg.xyY(3), 0);
SIM.dir.HeringRYGB.rgb = SIM.dir.HeringRYGB.rgb / STI.mon.rgb_max;

% RGB ---------------------------------------------------------------------
rgb = [1 0 0; 0 1 0; 0 0 1]*STI.mon.rgb_max;
%SIM.dir.rgb = colourconverter(rgb, 'rgb', 2, STI.wp.xyY, STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max, 1, 0);
SIM.dir.rgb = colourconverter(rgb, 'rgb', 2,...
    STI.wp.xyY, STI.sens, STI.cmf,...
    STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max,...
    STI.bg.xyY(3), 0);

%% ******************************* GRAPHIX ********************************

%% MODEL ILLUSTRATION (model)
fh = figure('Name','Simulation of cone adaptation', 'NumberTitle','off');
mod_model_plotter(SIM, STI, STI1a, colspace, DESIGN);

%% ****************************** SUBMODULES ******************************

%% mod_model_plotter
function axs = mod_model_plotter(SIM, STI, STI1a, colspace, DESIGN)
%2022.09.24
%2024.08.01 added predictions across hue for other models.
%2025.07-08 polished for publication.


% DESIGN ------------------------------------------------------------------
ax_rgb      = [1 1 1]*.8;
axmx        = DESIGN.axmax;
sim_rgb     = DESIGN.sim_rgb;
data_rgb    = DESIGN.data_rgb;

% PANEL 2 (CONE-OPPONENCY)
luv_rgb     = [0 0 0];

% PANEL 3 (APPEARANCE MODELS)
lab_rgb     = 'k';
cam_rgb     = [1 0 0];

% PANEL 4 (COLOUR-OPPONENCY)
mun_rgb     = [1 0 0];
rygb_rgb    = [0 0 0];

% SIMULATION SETTINGS -----------------------------------------------------
switch lower(colspace)
    case 'luv'
        polarspace = 'Luv_pol';
        xlbl = 'green-red [u*]';
        ylbl = 'blue-yellow [v*]';
    case 'dkl'
        polarspace = 'dkl_pol';        
        xlbl = 'L-M';
        ylbl = 'S-(L+M)';
end

MON     = STI.mon;
MON.lms = STI.mon.lms;
WP  = STI.wp;
bgY = STI.bg.xyY(3);

Inducer = [STI.inducer.(colspace)(:,1), round(STI.inducer.(polarspace))];
iChroma = mean(Inducer(:,3)); % This is only 1 constant number, depending on colspace.
L = mean(Inducer(:,1)); % Will be either 70 (Luv) or 0 (dkl)

CONESIM = colourconverter([Inducer(:,1), SIM.Oppo.Hue.ConeCon,SIM.Oppo.Chroma.ConeCon], polarspace, 2,...
    WP.xyY, STI.sens, STI.cmf,...
    MON.xyY, MON.ldt, MON.oog,[], MON.rgb_max,...
    bgY, 0);

DKLSIM = colourconverter([Inducer(:,1), SIM.Oppo.Hue.dkl,SIM.Oppo.Chroma.dkl], polarspace, 2,...
    WP.xyY, STI.sens, STI.cmf,...
    MON.xyY, MON.ldt, MON.oog,[], MON.rgb_max,...
    bgY, 0);

LUVSIM = colourconverter([Inducer(:,1), SIM.Oppo.Hue.Luv,SIM.Oppo.Chroma.Luv], polarspace, 2,...
    WP.xyY, STI.sens, STI.cmf,...
    MON.xyY, MON.ldt, MON.oog,[], MON.rgb_max,...
    bgY,0);

LABSIM = colourconverter([Inducer(:,1), SIM.Oppo.Hue.Lab,SIM.Oppo.Chroma.Lab], polarspace, 2,...
    WP.xyY, STI.sens, STI.cmf,...
    MON.xyY, MON.ldt, MON.oog,[], MON.rgb_max,...
    bgY,0);

CAM02 = colourconverter([Inducer(:,1), SIM.Oppo.Hue.cam02,SIM.Oppo.Chroma.cam02], polarspace, 2,...
    WP.xyY, STI.sens, STI.cmf,...
    MON.xyY, MON.ldt, MON.oog,[], MON.rgb_max,...
    bgY,0);

MUN = colourconverter([Inducer(:,1), SIM.Oppo.Hue.mun,SIM.Oppo.Chroma.mun], polarspace, 2,...
    WP.xyY, STI.sens, STI.cmf,...
    MON.xyY, MON.ldt, MON.oog,[], MON.rgb_max,...
    bgY,0);

RYGB = colourconverter([Inducer(:,1), SIM.Oppo.Hue.RYGB,SIM.Oppo.Chroma.RYGB], polarspace, 2,...
    WP.xyY, STI.sens, STI.cmf,...
    MON.xyY, MON.ldt, MON.oog,[], MON.rgb_max,...
    bgY, 0);

DESIGN.huecircle = huecircle_preparer(colspace, L, iChroma, STI);

% EXAMPLES ----------------------------------------------------------------
ex = 60;
ex1 = find(Inducer(:,2)==ex);
ex2 = find(Inducer(:,2)==mod(ex+180, 360));
EX = Inducer([ex1, ex2],1); 
[EX(:,2), EX(:,3)] = pol2cart(deg2rad(Inducer([ex1, ex2],2)), Inducer([ex1, ex2],3));
ex_mrksz = DESIGN.markersize * 5;
ex_rgb = STI.inducer.rgb([ex1, ex2],:);

% EXAMPLES AFTER ADAPTATION (EXoppo)
% NOTE: Variables are ALL in CIELUV!
EXoppo = table;
EXoppo.lms =   [CONESIM.(colspace)(ex1,:);   CONESIM.(colspace)(ex2,:)];
EXoppo.dkl =   [DKLSIM.(colspace)(ex1,:);    DKLSIM.(colspace)(ex2,:)];
EXoppo.luv =   [LUVSIM.(colspace)(ex1,:);    LUVSIM.(colspace)(ex2,:)];
EXoppo.lab =   [LABSIM.(colspace)(ex1,:);    LABSIM.(colspace)(ex2,:)];
EXoppo.cam =   [CAM02.(colspace)(ex1,:);     CAM02.(colspace)(ex2,:)];
EXoppo.mun =   [MUN.(colspace)(ex1,:);       MUN.(colspace)(ex2,:)];
EXoppo.rygb =  [RYGB.(colspace)(ex1,:);      RYGB.(colspace)(ex2,:)];

% 01 STIMULI --------------------------------------------------------------
axs(1) = subplot(1,5,1);
hold on
switch 'axes'
    case lower(DESIGN.addplots)
        axes_plotter(axmx, ax_rgb, '-', DESIGN.linewidth);
end

% DISK ILLUSTRATING COLOUR SPACE:
colourdisk(...
    300,...
    colspace,...
    iChroma, L, [-1 1 -1 1]*iChroma, 1, DESIGN.markersize);

% STIMULI FROM EXPERIMENT 1a
switch lower(colspace)
    case 'luv'
        plot([zeros(8,1), STI1a.Luv(:,2)]',[zeros(8,1), STI1a.Luv(:,3)]','k-', 'LineWidth', 2);
        for sti = 1:8
            plot([0, STI1a.Luv(sti,2)],[0, STI1a.Luv(sti,3)]','-', 'LineWidth', 1, 'Color', STI1a.rgb(sti,:));
        end
end

% STIMULI FROM EXPERIMENT 2a
x = STI.inducer.(colspace)(:,2); % u* if CIELUV, L-M if DKL
y = STI.inducer.(colspace)(:,3); % v* if CIELUV, S if DKL
rgb = STI.inducer.rgb;
plot(x, y,'k.');

% STIMULI FROM EXPERIMENT 1b
switch lower(colspace)
    case 'luv'
        inds = 1:3:numel(x); % for Experiment 1b
        mrksz = DESIGN.markersize * 5;
        scatter(x(inds), y(inds), mrksz, rgb(inds,:), 'filled', 'MarkerEdgeColor', 'k');
end

hold off

% FORMATTING --------------------------------------------------------------
set(gca, 'FontSize', DESIGN.fontsize,...
    'XTick', [-1 0 1]*iChroma,...
    'YTick', [-1 0 1]*iChroma);
axis([-1 1 -1 1]*axmx);
axis square;
xlabel(xlbl);
ylabel(ylbl);
title('STIMULI','FontWeight','bold');

% CONE ADAPTATION ---------------------------------------------------------
axs(2) = subplot(1,5,2);

hold on

% INDUCER HUE CIRCLE:
scatter(...
    DESIGN.huecircle.Lxy(:,2),...
    DESIGN.huecircle.Lxy(:,3),...
    50,... MarkerSize
    DESIGN.huecircle.rgb, 'filled', 'MarkerEdgeColor', 'none');

switch 'rgb'
    case lower(DESIGN.addplots)
        plot([0, SIM.dir.rgb.(colspace)(1,2)],[0, SIM.dir.rgb.(colspace)(1,3)],'-','Color',SIM.dir.rgb.rgb(1,:)/255);
        plot([0, SIM.dir.rgb.(colspace)(2,2)],[0, SIM.dir.rgb.(colspace)(2,3)],'-','Color',SIM.dir.rgb.rgb(2,:)/255);
        plot([0, SIM.dir.rgb.(colspace)(3,2)],[0, SIM.dir.rgb.(colspace)(3,3)],'-','Color',SIM.dir.rgb.rgb(3,:)/255);
end

% NOISY SIMULATION ---------------------------------------------------------
[x, y] = pol2cart(deg2rad(SIM.Noisy.ConeCon(:,:,1)), SIM.Noisy.ConeCon(:,:,2));
scatter(x(:), y(:), .5, 'k', 'filled', 'MarkerFaceColor', data_rgb, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', DESIGN.data_alpha);

switch 'cones'
    case lower(DESIGN.addplots)
        cone_plotter(SIM.dir.cones, colspace, axmx*1.5, ax_rgb, '-', DESIGN.linewidth);
end

% PURE SIMULATION ---------------------------------------------------------
plot(af(CONESIM.(colspace)(:,2),1),af(CONESIM.(colspace)(:,3),1),...
    '-','Color',sim_rgb,'LineWidth',DESIGN.fatlinewidth);

% EXAMPLES ----------------------------------------------------------------
plot([EX(1,2), 0, EXoppo.lms(1,2)], [EX(1,3), 0, EXoppo.lms(1,3)],'-',...
    'LineWidth', DESIGN.fatlinewidth, 'Color', ex_rgb(1,:), 'MarkerFaceColor','k');
plot([EX(2,2), 0, EXoppo.lms(2,2)], [EX(2,3), 0, EXoppo.lms(2,3)],'-',...
    'LineWidth', DESIGN.fatlinewidth, 'Color', ex_rgb(2,:), 'MarkerFaceColor','k');
scatter(EXoppo.lms(:,2), EXoppo.lms(:,3), ex_mrksz, ex_rgb,'o', 'filled', 'MarkerEdgeColor', 'k');

hold off

% FORMATTING --------------------------------------------------------------
set(gca, 'FontSize', DESIGN.fontsize,...
    'XTick', [-1 0 1]*iChroma,...
    'YTick', [-1 0 1]*iChroma, 'YTickLabel', []);
axis([-1 1 -1 1]*axmx);
axis square;
xlabel(xlbl);
ylabel(ylbl);
title('CONE-ADAPTATION','FontWeight','bold');

% CONE OPPONENCY ----------------------------------------------------------
axs(3) = subplot(1,5,3);
hold on

% INDUCER HUE CIRCLE:
scatter(...
    DESIGN.huecircle.Lxy(:,2),...
    DESIGN.huecircle.Lxy(:,3),...
    50,... MarkerSize
    DESIGN.huecircle.rgb, 'filled', 'MarkerEdgeColor', 'none');

% DKL AXES:
dkl_plotter(SIM.dir.dkl, colspace, ax_rgb, axmx, DESIGN.linewidth, 0);

% LUV SIM:
plot(af(LUVSIM.(colspace)(:,2)), af(LUVSIM.(colspace)(:,3)),...
    '-','Color',luv_rgb,'LineWidth',DESIGN.linewidth);

% DKL SIM:
plot(af(DKLSIM.(colspace)(:,2)), af(DKLSIM.(colspace)(:,3)),...
    '-','Color',sim_rgb,'LineWidth',DESIGN.fatlinewidth);

%plot(EX([1 3],2), EX([1 3],3),'ko-', 'LineWidth', 1, 'MarkerFaceColor','k');
plot([EX(1,2), 0, EXoppo.dkl(1,2)], [EX(1,3), 0, EXoppo.dkl(1,3)],...
    '-', 'LineWidth', DESIGN.fatlinewidth, 'Color', ex_rgb(1,:), 'MarkerFaceColor','k');
plot([EX(2,2), 0, EXoppo.dkl(2,2)], [EX(2,3), 0, EXoppo.dkl(2,3)],...
    '-', 'LineWidth', DESIGN.fatlinewidth, 'Color', ex_rgb(2,:), 'MarkerFaceColor','k');
plot([EXoppo.dkl(1,2), 0, EXoppo.dkl(2,2)], [EXoppo.dkl(1,3), 0, EXoppo.dkl(2,3)],...
    ':', 'LineWidth', DESIGN.fatlinewidth, 'Color', ex_rgb(1,:), 'MarkerFaceColor','k');
%scatter(EX([1 3],2), EX([1 3],3), 50, ex_rgb([1 3],:),'ko', 'filled');
%scatter(EX([2 4],2), EX([2 4],3), 200, rgb,'p', 'filled', 'MarkerEdgeColor', 'k');
scatter(EXoppo.dkl(:,2), EXoppo.dkl(:,3), ex_mrksz, ex_rgb,'o', 'filled', 'MarkerEdgeColor', 'k');
scatter(EXoppo.luv(:,2), EXoppo.luv(:,3), ex_mrksz, ex_rgb,'o', 'filled', 'MarkerEdgeColor', 'k');

hold off

% FORMATTING --------------------------------------------------------------
set(gca, 'FontSize', DESIGN.fontsize,...
        'XTick', [-1 0 1]*iChroma,...
        'YTick', [-1 0 1]*iChroma, 'YTickLabel', []);
axis square;
axis([-1 1 -1 1]*axmx);
axis square;
xlabel(xlbl);
ylabel(ylbl);
title('CONE-OPPONENCY','FontWeight','bold');

% APPEARANCE MODELS 1 -----------------------------------------------------
axs(4) = subplot(1,5,4);
hold on
switch 'axes'
    case lower(DESIGN.addplots)
        axes_plotter(axmx, ax_rgb, '-', DESIGN.linewidth);
end

% INDUCER HUE CIRCLE:
scatter(...
    DESIGN.huecircle.Lxy(:,2),...
    DESIGN.huecircle.Lxy(:,3),...
    50,... MarkerSize
    DESIGN.huecircle.rgb, 'filled', 'MarkerEdgeColor', 'none');
plot(af(LABSIM.(colspace)(:,2)), af(LABSIM.(colspace)(:,3)),...
    '-','Color',lab_rgb,'LineWidth',DESIGN.linewidth);
plot(af(CAM02.(colspace)(:,2)), af(CAM02.(colspace)(:,3)),...
    '-','Color',cam_rgb,'LineWidth',DESIGN.fatlinewidth);

% EXAMPLE LINES -----------------------------------------------------------

% LAB
plot([0, EXoppo.lab(1,2)], [0, EXoppo.lab(1,3)],'-',...
    'LineWidth', DESIGN.linewidth, 'Color', ex_rgb(1,:), 'MarkerFaceColor','k');
plot([0, EXoppo.lab(2,2)], [0, EXoppo.lab(2,3)],'-',...
    'LineWidth', DESIGN.linewidth, 'Color', ex_rgb(2,:), 'MarkerFaceColor','k');
scatter(EXoppo.lab(:,2), EXoppo.lab(:,3), ex_mrksz, ex_rgb, 'o', 'filled', 'MarkerEdgeColor', 'k');

% CIECAM
plot([0, EXoppo.cam(1,2)], [0, EXoppo.cam(1,3)],'-',...
    'LineWidth', DESIGN.linewidth, 'Color', ex_rgb(1,:), 'MarkerFaceColor','k');
plot([0, EXoppo.cam(2,2)], [0, EXoppo.cam(2,3)],'-',...
    'LineWidth', DESIGN.linewidth, 'Color', ex_rgb(2,:), 'MarkerFaceColor','k');
scatter(EXoppo.cam(:,2), EXoppo.cam(:,3), ex_mrksz, ex_rgb,'o', 'filled', 'MarkerEdgeColor', 'k');

hold on

% FORMATTING --------------------------------------------------------------
set(gca, 'FontSize', DESIGN.fontsize,...
    'XTick', [-1 0 1]*iChroma,...
    'YTick', [-1 0 1]*iChroma, 'YTickLabel', []);
axis([-1 1 -1 1]*axmx);
axis square;
xlabel(xlbl);
ylabel(ylbl);
title('APPEARANCE MODELS 1','FontWeight','bold');

% APPEARANCE MODELS 2 -----------------------------------------------------
axs(5) = subplot(1,5,5);
hold on

% INDUCER HUE CIRCLE:
scatter(...
    DESIGN.huecircle.Lxy(:,2),...
    DESIGN.huecircle.Lxy(:,3),...
    50,... MarkerSize
    DESIGN.huecircle.rgb, 'filled', 'MarkerEdgeColor', 'none');

typ_x = SIM.dir.HeringRYGB.(colspace)(:,2)*1.5*axmx;
typ_y = SIM.dir.HeringRYGB.(colspace)(:,3)*1.5*axmx;
rgb2 = [1	0	0;...
        1	1	0;...
        0	0.6	0;...
        0	0.2	1];
rgb2 = SIM.dir.HeringRYGB.rgb;
for cn = 1:4
    plot([0,typ_x(cn)],[0,typ_y(cn)], '-', 'Color', 'k','LineWidth',2);
    plot([0,typ_x(cn)],[0,typ_y(cn)], 'Color', rgb2(cn,:),'LineWidth',1);
%    plot([0,typ_x(cn)],[0,typ_y(cn)], ':', 'Color', 'k','LineWidth',2);
end

plot(af(MUN.(colspace)(:,2)), af(MUN.(colspace)(:,3)),...
    '-','Color', mun_rgb,'LineWidth',DESIGN.fatlinewidth);

% EXAMPLES:
scatter(EXoppo.mun(:,2), EXoppo.mun(:,3), ex_mrksz, ex_rgb,'o', 'filled', 'MarkerEdgeColor', 'k');
plot(af(RYGB.(colspace)(:,2)), af(RYGB.(colspace)(:,3)),...
    '.','Color', rygb_rgb, 'LineWidth',DESIGN.linewidth);
scatter(EXoppo.rygb(:,2), EXoppo.rygb(:,3), ex_mrksz, ex_rgb,'o', 'filled', 'MarkerEdgeColor', 'k');
hold on

% FORMATTING --------------------------------------------------------------
set(gca, 'FontSize', DESIGN.fontsize,...
    'XTick', [-1 0 1]*iChroma,...
    'YTick', [-1 0 1]*iChroma, 'YTickLabel', []);
axis([-1 1 -1 1]*axmx);
axis square;
xlabel(xlbl);
ylabel(ylbl);
title('APPEARANCE MODELS 2','FontWeight','bold');

%% ***************************** SUBPLOTTERS ******************************

%% axes_plotter
function axes_plotter(mx, rgb, linestyle, linewidth)
% 2022.09.13
% 2025.04.17

if nargin < 3
    linestyle = '-';
    if nargin < 2
        rgb = [.8 .8 .8];
    end
end

% AXIS:
plot([-1 1; 0 0]'*mx, [0 0; -1 1]'*mx, linestyle, 'Color', rgb, 'LineWidth', linewidth);

%% circle
function c = circle(radius, imradius)
% 2015.07.03

if nargin < 2
  imradius = radius;
end

x = (-imradius: imradius);
y = x;
[X, Y] = meshgrid(x, y);

z = sqrt(X.^2 + Y.^2);

s = sqrt(2*imradius^2);
% c = 1 - double(im2bw(z/s, radius/s));

% To avoid image processing toolbox
c = z/s;
inds = c <= radius/s;
c(inds) = 1; 
c(~inds) = 0;

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

%% colourdisk
function colourdisk(resol, colspace, chroma, lum, ax, withcircle, markersize)
%2024.07.28 based on dkldisk
%2025.04.17 added markersize 

if nargin < 7
    markersize = 20;
if nargin < 6
%    withcircle = 0:15:345;
    withcircle = 1;     
    if nargin < 5
        ax = [];
        if nargin < 4
            lum = 0;
            if nargin < 3
                chroma = 1;
                if nargin < 2
                    colspace = 'dkl';
                    if nargin == 0
                        clc; close all;
                        resol = 300;
                    end
                end
            end
        end
    end
end
end

[mon_xyY, mon_ldt] = srgb;
sz = resol;
mask = circle((sz-1)/2);
L = ones(sz,sz)*lum;
RG = ones(sz,1)*((-(sz-1):2:(sz-1))/(sz-1))*chroma;
BY = ((-(sz-1):2:(sz-1))/(sz-1))'*ones(1,sz)*chroma;
RG = RG .*mask;
BY = BY.*mask;


if isempty(ax)
    ax = [-sz/2 sz/2, -sz/2 sz/2];
end
wd = abs(ax(1)-ax(2));
ht = abs(ax(3)-ax(4));

colspecs = cat(3, L, RG, BY);
%COL = colourconverter(colspecs, colspace, 3, 'mon', mon_xyY, mon_ldt);
COL = colourconverter(colspecs, colspace, 3,...
    'mon', 'ss', '1931',...
    mon_xyY, mon_ldt, [], [], [],...
    [], 0);
rgb1 = COL.rgb/255;
hold on
image(ax(1:2), ax(3:4), rgb1);
if ~isempty(withcircle) & ~isnan(withcircle) & sum(withcircle)>0
    if numel(withcircle) > 1
        azimuth_list = withcircle(:);
    else
        azimuth_list = 0:0.5:360;
    end
    luminance = 0; radius = 1;
    aziradi = deg2rad(azimuth_list);
    [x,y] = pol2cart(aziradi(:), chroma);
    colspecs1 = [ones(size(azimuth_list(:)))*lum, x, y];
    COL1 = colourconverter(colspecs1, colspace, 2,...
        'mon', 'ss', '1931', mon_xyY, mon_ldt, [], [],...
        [], 0);
    x = cos(aziradi)*wd/2;
    y = sin(aziradi)*ht/2;
%    scatter(x,y,markersize+2, 'k', 'filled')
    scatter(x,y,markersize, COL1.rgb/255, 'filled')
end
set(gca, 'Color', shiftdim(COL.rgb(1,1,:))'/255); % Top left is grey background :)
axis xy;
return

%% cone_plotter
function cone_plotter(COL, colspace, mx, rgb, linestyle, linewidth)
% 2022.09.13
% 2025.04.17
% 2025.05.03 added colspace

if nargin < 6
    linewidth = 1;
    if nargin < 5
        linestyle = '-';
        if nargin < 4
            rgb = [.8 .8 .8];
            if nargin < 3
                mx = 100;
                if nargin < 2
                    colspace = 'Luv';
                end
            end
        end
    end
end


if isempty(rgb)
    rgb = [1 0 0; 0 0.8 0; 0 0 1];
elseif size(rgb,1) == 1
    rgb = ones(3,1)*rgb;
end

% RESCALE -----------------------------------------------------------------
switch lower(colspace)
    case 'luv'
        cones = COL.Luv;
        cones(:,2) = (cones(:,2)./COL.Luv_pol(:,2))*mx;
        cones(:,3) = (cones(:,3)./COL.Luv_pol(:,2))*mx;
    case 'dkl'
        cones = COL.dkl;
        cones(:,2) = (cones(:,2)./COL.dkl_pol(:,2))*mx;
        cones(:,3) = (cones(:,3)./COL.dkl_pol(:,2))*mx;
end

% PLOT --------------------------------------------------------------------
plot([0, cones(1,2)],[0, cones(1,3)],...
    linestyle, 'Color', rgb(1,:), 'LineWidth',linewidth);
plot([0, cones(2,2)],[0, cones(2,3)],...
    linestyle, 'Color', rgb(2,:), 'LineWidth',linewidth);
plot([0, cones(3,2)],[0, cones(3,3)],...
    linestyle, 'Color', rgb(3,:), 'LineWidth',linewidth);

%% dkl_plotter
function h = dkl_plotter(COL, colspace, rgb, axmax, linewidth, rot)
%2020.09.24
% 2025.04.17
% 2025.05.03 made it DKL space compatible.

% COLOUR ------------------------------------------------------------------
if isempty(rgb)
    rgb = COL.rgb/255;
elseif size(rgb,1) == 1 % uniform colours
    rgb = ones(4,1)*rgb;
end

switch lower(colspace)
    case {'luv', 'luv_pol'}
        Lxy = COL.Luv;
        mn = min(COL.Luv_pol(:,2));
    case {'dkl', 'dkl_pol'}
        Lxy = COL.dkl;
        mn = min(COL.dkl_pol(:,2));
end
Lxy(:,2:3) = Lxy(:,2:3)./mn*axmax;
hold on
h = plot([0 Lxy(1,2)], [0 Lxy(1,3)],'-', 'Color', rgb(1,:), 'LineWidth', linewidth);
plot([0 Lxy(2,2)], [0 Lxy(2,3)],  '-', 'Color', rgb(2,:), 'LineWidth', linewidth);
plot([0 Lxy(3,2)],[0 Lxy(3,3)], '-', 'Color', rgb(3,:), 'LineWidth', linewidth);
plot([0 Lxy(4,2)],[0 Lxy(4,3)], '-', 'Color', rgb(4,:), 'LineWidth', linewidth);

return % no text (too much)
if rot
    ang = rad2deg(cart2pol(Luv([2,3],2),Luv([2,3],3)));
else
    ang = [0 0];
end
%txtcol = rgb(1,:)*.8;
txtcol = [0 0 0];
text(-axmax,0,' L-M',...
    'Rotation', ang(1), 'Color', txtcol, ...
    'VerticalAlignment', 'Top', 'HorizontalAlignment','Left');
%txtcol = rgb(3,:)*.8;
txtcol = [0 0 0];
text(0,-axmax,'(L+M)-S',...
    'Rotation', ang(2),'Color', txtcol,...
    'VerticalAlignment', 'Bottom', 'HorizontalAlignment','Right');

%% huecircle_preparer
function huecircle = huecircle_preparer(colspace, lightness, chroma, STI)
% 2025.07.31 [cw]

azi = (0:360)';
switch lower(colspace)
    case 'luv'
        chro = ones(361,1)*chroma;
        L = ones(361,1)*lightness;
        [x,y] = pol2cart(deg2rad(azi), chro);
    case 'dkl'
        chro = ones(361,1)*chroma;
        L = ones(361,1)*lightness;
        [x,y] = pol2cart(deg2rad(azi), chro);
end
%COL = colourconverter([L,x,y], colspace, 2, STI.wp.xyY, STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max,1,0);
COL = colourconverter([L,x,y], colspace, 2,...
    STI.wp.xyY, STI.sens, STI.cmf,...
    STI.mon.xyY, STI.mon.ldt, STI.mon.oog, [], STI.mon.rgb_max,...
    STI.bg.xyY(3),0);

rgb = COL.rgb/STI.mon.rgb_max;
Lxy = COL.(colspace);
[hue, chro] = cart2pol(Lxy(:,2),Lxy(:,3));
hue = rad2deg(hue);

huecircle = table(hue, chro, Lxy, rgb);

%% lab_plotter
function h = lab_plotter(rgb, axmax, rot, linewidth)
%2020.09.24

chroma = (0:1:axmax)';
L = ones(size(chroma))*70;
nuller = L*0;
a1 = [L -chroma nuller];
a2 = [L chroma nuller]; 
b1 = [L nuller -chroma]; 
b2 = [L nuller chroma];
A1 = colourconverter(a1, 'Lab', 2);
A2 = colourconverter(a2, 'Lab', 2);
B1 = colourconverter(b1, 'Lab', 2);
B2 = colourconverter(b2, 'Lab', 2);

% COLOUR ------------------------------------------------------------------
if isempty(rgb)
    rgb(1,:) = A1.rgb(end,:)/255;
    rgb(2,:) = A2.rgb(end,:)/255;
    rgb(3,:) = B1.rgb(end,:)/255;
    rgb(4,:) = B2.rgb(end,:)/255;
else % uniform colours
    rgb = ones(3,1)*rgb; 
end

scaler = 2;
hold on
h = plot([nuller A1.Luv(:,2)]*scaler, [nuller A1.Luv(:,3)]*scaler,'-', 'Color', rgb(1,:), 'LineWidth', linewidth);
plot([nuller A2.Luv(:,2)]*scaler, [nuller A2.Luv(:,3)]*scaler,  '-', 'Color', rgb(2,:), 'LineWidth', linewidth);
plot([nuller B1.Luv(:,2)]*scaler, [nuller B1.Luv(:,3)]*scaler, '-', 'Color', rgb(3,:), 'LineWidth', linewidth);
plot([nuller B2.Luv(:,2)]*scaler, [nuller B2.Luv(:,3)]*scaler, '-', 'Color', rgb(4,:), 'LineWidth', linewidth);

if rot
    ang = rad2deg(cart2pol(Luv([2,3],2),Luv([2,3],3)));
else
    ang = [0 0];
end
text(-axmax,0,' a*',...
    'Rotation', ang(1), 'Color', rgb(1,:), ...
    'VerticalAlignment', 'Top', 'HorizontalAlignment','Left');
text(0,-axmax,'b*',...
    'Rotation', ang(2),'Color', rgb(3,:),...
    'VerticalAlignment', 'Bottom', 'HorizontalAlignment','Right');

%% radius_plotter
function radius_plotter(radius,rgb)
%2020.06.24
if nargin < 2
    rgb = [.5 .5 .5];
end

azi = deg2rad(0:360)';
[x,y] = pol2cart(azi, ones(numel(azi),1)*radius);
plot(x,y, 'Color', rgb);

%% huefreq_plotter
function h_fill = huefreq_plotter(hue, frq)
%2022.09.13 [cw]

[x,y] = pol2cart(deg2rad(hue),frq);
h_fill = fill(x,y,[.8 .8 .8], 'EdgeColor', [.5 .5 .5]);

%% exp1a_stimulus_maker
function STI1a = exp1a_stimulus_maker
% 2025.07.29

% STIMULI IN EXPERIMENT 1.a:
STI1a.rgb = [...
    1	0	0
    1	0.600195503421310	0
    1	1	0
    0	1	0
    0	0.900293255131965	0.799608993157380
    0	0.650048875855327	1
    0.600195503421310	0	1
    1	0	1];

STI1a.Luv = [...
    58.1361334527991	233.868625775788	29.4774827738170	7.18384397141728	235.719019412809
    74.8757290730455	118.600110381531	61.9411405547026	27.5766620061866	133.800938246818
    98.1406558153255	3.77171735686390	101.102179775460	87.8635147423300	101.172509147344
    86.5849429776959	-141.121407716693	108.163536831946	142.531375359917	177.804956106226
    80.4637362867371	-108.283558007562	8.04712913324397	175.749857517590	108.582158857080
    62.5254822946730	-60.6482737092590	-88.7804063787166	235.661970665063	107.518248035779
    42.8233526810712	46.0239722842757	-133.228629270313	289.057526514542	140.954154539943
    62.3693777048737	138.584926045370	-106.219431872120	322.531375359917	174.609133307043
    ];

%% stimulusstruc_maker
function stimuli = stimulusstruc_maker(colspace)
% 2025.07.29

switch lower(colspace)
    case 'luv'
        L0 = 70;
        C0 = 50;
        stimuli.src_space   = 'Luv';
        stimuli.polarspace  = 'luv_pol';
        stimuli.sens        = 'ss'; 
        stimuli.cmf         = '1931';
        bgY = [];
    case 'dkl'
        L0 = 0;
        C0 = 1;
        bgY = 'halfmon';
        stimuli.src_space   = 'dkl';
        stimuli.polarspace  = 'dkl_pol';
        stimuli.sens        = 'sp'; 
        stimuli.cmf         = 'judd';
end

% MONITOR (sRGB) ----------------------------------------------------------
[stimuli.mon.xyY, stimuli.mon.ldt] = srgb;
stimuli.mon.oog = [];
stimuli.mon.rgb_max = stimuli.mon.ldt(end);
% CALCULATE ACCURATE MONITOR LMS FOR DKL CONVERSION
stimuli.mon.lms = XYZ2lms(xyY2XYZ(stimuli.mon.xyY), stimuli.sens, stimuli.cmf);

% WHITE-POINT -------------------------------------------------------------
%stimuli.wp = colourconverter([1 1 1]*stimuli.mon.rgb_max, 'rgb', 2, 'mon', stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max,1,0);
stimuli.wp = colourconverter([1 1 1]*stimuli.mon.rgb_max, 'rgb', 2,...
    'mon', stimuli.sens, stimuli.cmf,...
    stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max,...
    bgY, 0);
stimuli.wp.rgb = stimuli.wp.rgb/stimuli.mon.rgb_max;

% BACKGROUND --------------------------------------------------------------
%stimuli.bg = colourconverter([L0, 0  0], stimuli.polarspace, 2, stimuli.wp.xyY, stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max, 1, 0);
stimuli.bg = colourconverter([L0, 0  0], stimuli.polarspace, 2,...
    stimuli.wp.xyY, stimuli.sens, stimuli.cmf,...
    stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max,...
    bgY, 0);
%stimuli.bg.lms = XYZ2lms(stimuli.bg.XYZ, stimuli.sens, stimuli.cmf);
stimuli.bg.rgb = stimuli.bg.rgb/stimuli.mon.rgb_max;

% INDUCER -----------------------------------------------------------------
hue = (0:5:355)';
col_n = numel(hue);
L = L0*ones(col_n,1);
C = C0*ones(col_n,1);
Inducer = [L, hue, C];
%stimuli.inducer = colourconverter(Inducer, stimuli.polarspace, 2, stimuli.wp.xyY, stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max, 1, 0);
stimuli.inducer = colourconverter(Inducer, stimuli.polarspace, 2,...
    stimuli.wp.xyY, stimuli.sens, stimuli.cmf,...
    stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max,...
    stimuli.bg.xyY(3), 0);
%[stimuli.inducer.dkl, stimuli.inducer.dkl_pol, stimuli.inducer.lms] = XYZ2dkl(stimuli.inducer.XYZ, stimuli.bg.lms, stimuli.sens, stimuli.cmf, stimuli.mon.lms);
stimuli.inducer.rgb = stimuli.inducer.rgb/stimuli.mon.rgb_max;

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

%% cdissolver
function outM = cdissolver(inM, onecycle)
% 2012.08.23 * [cw]
% 2024.07.14 changed min to median [cw]

if nargin < 2
    onecycle = 360;
end

%inM = mod(inM, onecycle);
for cl = 1:size(inM,2)
    refM1 = ones(size(inM,1),1) * min(inM(:,cl));
    refM2 = ones(size(inM,1),1) * nanmedian(inM(:,cl));
    dffM1 = circularsubstractor(inM(:,cl), refM1, onecycle, 'signed');
    dffM2 = circularsubstractor(inM(:,cl), refM2, onecycle, 'signed');
    if nanmean(abs(dffM2)) < nanmean(abs(dffM1))
        refM = refM2; dffM = dffM2;
    else
        refM = refM1; dffM = dffM1;
    end
%    refM = refM2; dffM = dffM2;
    outM(:,cl) = refM+dffM;
end

%% cstats
function varargout = cstats(inM, stype, onecycle, remap)
% 2012aug23 [cw]
% 2012sep25 mapped output to cycle [cw]
% 2012sep27 fixed dimension to 1 [cw]
% 2013aug14 deleted dimension specification to alow for min and max, which
% specify dim as the 3. input instead of 2. [cw]
% 2015.10.30 made modulo of result optional with "remap" [cw].

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

%% dialogue_setter
function [SETTINGS] = dialogue_setter(SETTINGS)
% 2009oct8 [cw]
fldnms = fieldnames(SETTINGS);
counter = 0;
for fld = 1:numel(fldnms)
    switch class(SETTINGS.(fldnms{fld}))
        case 'double'
            counter = counter + 1;
            dnames{counter,:}  = fldnms{fld};
            dvalues{counter,:} = num2str(SETTINGS.(fldnms{fld}));
            dformat{counter,:} = 'double';
        case 'char'
            counter = counter + 1;
            dnames{counter,:}  = fldnms{fld};
            dvalues{counter,:} = SETTINGS.(fldnms{fld});
            dformat{counter,:} = 'char';
    end
end
doutput = inputdlg(...
    dnames,...
    'INPUT',...
    1,...
    dvalues);

for fld = 1:numel(dnames)
    switch dformat{fld}
        case 'double'
            SETTINGS.(dnames{fld}) = str2num(doutput{fld});
        case 'char'
            SETTINGS.(dnames{fld}) = doutput{fld};
    end
end

%% daylighter
function day = daylighter(wp_xyY, mon_xyY, mon_ldt, colortemperature, rectify_rgb)
% 2016.03.16 * based on ana_dress_qad [cw]
% 2018.07.11 added optional input colortemperature
% 2018.07.22 disentangled luminance of daylight locus and RGB > added rectify_rgb [cw]

if nargin < 5
    rectify_rgb = 1;
if nargin < 4
    colortemperature = [1500:10:40000,41000:1000:1000000];
if nargin < 2
    [mon_xyY, mon_ldt] = srgb;
    %     mon_xyY = [...
    %         0.615,0.3506,14.4;...
    %         0.295,0.6,43.5;...
    %         0.144,0.076,5.16;...
    %         0.31133,0.34767,60.4];
    if nargin == 0
        wp_xyY = XYZ2xyY(rgb2XYZ([1 1 1],mon_xyY(1:3,:),2));
    end
end
end
end

mon_xyY = mon_xyY(1:3,:);

load B_cieday;
[~, x, y] = GenerateCIEDay(colortemperature, B_cieday);

xyY = [x; y; ones(1,size(x,2))*wp_xyY(3)]';
XYZ = xyY2XYZ(xyY);
Luv = xyY2Luv(xyY,wp_xyY);

% Rectify for RGB to fit into monitor gamut
if rectify_rgb
    xyYr = [x; y; ones(1,size(x,2))*sum(mon_xyY(1:3,3))/2]';
else
    xyYr = xyY;
end

XYZr = xyY2XYZ(xyYr);
[r, g, b] = XYZ2rgb(XYZr(:,1),XYZr(:,2),XYZr(:,3), mon_xyY);
rgbl = [r,g,b];
rgbl(rgbl<0) = 0;
rgbl(rgbl>255) = 255;
rgb = gammacorrector(rgbl, mon_ldt, 2)/255;

day = table(Luv, xyY, XYZ, rgbl, rgb);

%figure; scatter(day.xyY(:,1),day.xyY(:,2),30, day.rgb, 'filled');

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
% 2011.04.06 [cw]
if nargin < 2
    style = '_';
    if nargin < 1
        width = 100;
    end
end
ln(1:width) = style;

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

%% srgb
function [mon_xyY, mon_ldt] = srgb
%2018.01.05 * [cw]

mon_xyY = [0.64,0.33,17;0.3,0.6,57.2;0.15,0.06,5.8;0.3127,0.329,80];
mon_ldt = [0,0,0;0.077399,0.077399,0.077399;0.154799,0.154799,0.154799;0.232198,0.232198,0.232198;0.309598,0.309598,0.309598;0.386997,0.386997,0.386997;0.464396,0.464396,0.464396;0.541796,0.541796,0.541796;0.619195,0.619195,0.619195;0.696594,0.696594,0.696594;0.773994,0.773994,0.773994;0.853367,0.853367,0.853367;0.937509,0.937509,0.937509;1.026303,1.026303,1.026303;1.119818,1.119818,1.119818;1.218123,1.218123,1.218123;1.321287,1.321287,1.321287;1.429375,1.429375,1.429375;1.542452,1.542452,1.542452;1.660583,1.660583,1.660583;1.78383,1.78383,1.78383;1.912253,1.912253,1.912253;2.045914,2.045914,2.045914;2.184872,2.184872,2.184872;2.329185,2.329185,2.329185;2.47891,2.47891,2.47891;2.634105,2.634105,2.634105;2.794824,2.794824,2.794824;2.961123,2.961123,2.961123;3.133055,3.133055,3.133055;3.310673,3.310673,3.310673;3.494031,3.494031,3.494031;3.68318,3.68318,3.68318;3.878171,3.878171,3.878171;4.079055,4.079055,4.079055;4.285881,4.285881,4.285881;4.498698,4.498698,4.498698;4.717556,4.717556,4.717556;4.942502,4.942502,4.942502;5.173584,5.173584,5.173584;5.410848,5.410848,5.410848;5.654341,5.654341,5.654341;5.904108,5.904108,5.904108;6.160196,6.160196,6.160196;6.422649,6.422649,6.422649;6.691512,6.691512,6.691512;6.966827,6.966827,6.966827;7.24864,7.24864,7.24864;7.536993,7.536993,7.536993;7.831928,7.831928,7.831928;8.133488,8.133488,8.133488;8.441715,8.441715,8.441715;8.756651,8.756651,8.756651;9.078335,9.078335,9.078335;9.40681,9.40681,9.40681;9.742115,9.742115,9.742115;10.08429,10.08429,10.08429;10.4333750000000,10.4333750000000,10.4333750000000;10.78941,10.78941,10.78941;11.1524320000000,11.1524320000000,11.1524320000000;11.5224820000000,11.5224820000000,11.5224820000000;11.8995970000000,11.8995970000000,11.8995970000000;12.2838150000000,12.2838150000000,12.2838150000000;12.6751740000000,12.6751740000000,12.6751740000000;13.0737120000000,13.0737120000000,13.0737120000000;13.4794650000000,13.4794650000000,13.4794650000000;13.89247,13.89247,13.89247;14.3127650000000,14.3127650000000,14.3127650000000;14.7403850000000,14.7403850000000,14.7403850000000;15.1753660000000,15.1753660000000,15.1753660000000;15.6177440000000,15.6177440000000,15.6177440000000;16.0675550000000,16.0675550000000,16.0675550000000;16.5248330000000,16.5248330000000,16.5248330000000;16.9896140000000,16.9896140000000,16.9896140000000;17.4619330000000,17.4619330000000,17.4619330000000;17.9418240000000,17.9418240000000,17.9418240000000;18.4293220000000,18.4293220000000,18.4293220000000;18.92446,18.92446,18.92446;19.4272720000000,19.4272720000000,19.4272720000000;19.9377930000000,19.9377930000000,19.9377930000000;20.4560540000000,20.4560540000000,20.4560540000000;20.98209,20.98209,20.98209;21.5159340000000,21.5159340000000,21.5159340000000;22.0576180000000,22.0576180000000,22.0576180000000;22.6071750000000,22.6071750000000,22.6071750000000;23.1646360000000,23.1646360000000,23.1646360000000;23.7300360000000,23.7300360000000,23.7300360000000;24.3034040000000,24.3034040000000,24.3034040000000;24.8847740000000,24.8847740000000,24.8847740000000;25.4741760000000,25.4741760000000,25.4741760000000;26.0716420000000,26.0716420000000,26.0716420000000;26.6772030000000,26.6772030000000,26.6772030000000;27.2908910000000,27.2908910000000,27.2908910000000;27.9127360000000,27.9127360000000,27.9127360000000;28.5427690000000,28.5427690000000,28.5427690000000;29.18102,29.18102,29.18102;29.82752,29.82752,29.82752;30.4822990000000,30.4822990000000,30.4822990000000;31.1453870000000,31.1453870000000,31.1453870000000;31.8168130000000,31.8168130000000,31.8168130000000;32.4966090000000,32.4966090000000,32.4966090000000;33.1848020000000,33.1848020000000,33.1848020000000;33.8814220000000,33.8814220000000,33.8814220000000;34.5864990000000,34.5864990000000,34.5864990000000;35.3000620000000,35.3000620000000,35.3000620000000;36.0221390000000,36.0221390000000,36.0221390000000;36.75276,36.75276,36.75276;37.4919530000000,37.4919530000000,37.4919530000000;38.2397460000000,38.2397460000000,38.2397460000000;38.9961690000000,38.9961690000000,38.9961690000000;39.7612480000000,39.7612480000000,39.7612480000000;40.5350130000000,40.5350130000000,40.5350130000000;41.3174910000000,41.3174910000000,41.3174910000000;42.10871,42.10871,42.10871;42.9086970000000,42.9086970000000,42.9086970000000;43.7174810000000,43.7174810000000,43.7174810000000;44.5350880000000,44.5350880000000,44.5350880000000;45.3615460000000,45.3615460000000,45.3615460000000;46.1968820000000,46.1968820000000,46.1968820000000;47.0411240000000,47.0411240000000,47.0411240000000;47.8942970000000,47.8942970000000,47.8942970000000;48.7564290000000,48.7564290000000,48.7564290000000;49.6275470000000,49.6275470000000,49.6275470000000;50.5076760000000,50.5076760000000,50.5076760000000;51.3968450000000,51.3968450000000,51.3968450000000;52.2950780000000,52.2950780000000,52.2950780000000;53.2024020000000,53.2024020000000,53.2024020000000;54.1188430000000,54.1188430000000,54.1188430000000;55.0444280000000,55.0444280000000,55.0444280000000;55.9791810000000,55.9791810000000,55.9791810000000;56.9231290000000,56.9231290000000,56.9231290000000;57.8762980000000,57.8762980000000,57.8762980000000;58.8387120000000,58.8387120000000,58.8387120000000;59.8103980000000,59.8103980000000,59.8103980000000;60.7913810000000,60.7913810000000,60.7913810000000;61.7816860000000,61.7816860000000,61.7816860000000;62.7813380000000,62.7813380000000,62.7813380000000;63.7903630000000,63.7903630000000,63.7903630000000;64.8087840000000,64.8087840000000,64.8087840000000;65.8366270000000,65.8366270000000,65.8366270000000;66.8739180000000,66.8739180000000,66.8739180000000;67.9206790000000,67.9206790000000,67.9206790000000;68.9769370000000,68.9769370000000,68.9769370000000;70.0427150000000,70.0427150000000,70.0427150000000;71.1180370000000,71.1180370000000,71.1180370000000;72.2029290000000,72.2029290000000,72.2029290000000;73.2974140000000,73.2974140000000,73.2974140000000;74.4015160000000,74.4015160000000,74.4015160000000;75.5152590000000,75.5152590000000,75.5152590000000;76.6386680000000,76.6386680000000,76.6386680000000;77.7717650000000,77.7717650000000,77.7717650000000;78.9145750000000,78.9145750000000,78.9145750000000;80.0671220000000,80.0671220000000,80.0671220000000;81.2294280000000,81.2294280000000,81.2294280000000;82.4015180000000,82.4015180000000,82.4015180000000;83.5834150000000,83.5834150000000,83.5834150000000;84.7751420000000,84.7751420000000,84.7751420000000;85.9767220000000,85.9767220000000,85.9767220000000;87.1881780000000,87.1881780000000,87.1881780000000;88.4095340000000,88.4095340000000,88.4095340000000;89.6408130000000,89.6408130000000,89.6408130000000;90.8820370000000,90.8820370000000,90.8820370000000;92.1332290000000,92.1332290000000,92.1332290000000;93.3944120000000,93.3944120000000,93.3944120000000;94.6656090000000,94.6656090000000,94.6656090000000;95.9468410000000,95.9468410000000,95.9468410000000;97.2381330000000,97.2381330000000,97.2381330000000;98.5395060000000,98.5395060000000,98.5395060000000;99.8509820000000,99.8509820000000,99.8509820000000;101.172584000000,101.172584000000,101.172584000000;102.504334000000,102.504334000000,102.504334000000;103.846254000000,103.846254000000,103.846254000000;105.198366000000,105.198366000000,105.198366000000;106.560693000000,106.560693000000,106.560693000000;107.933256000000,107.933256000000,107.933256000000;109.316077000000,109.316077000000,109.316077000000;110.709177000000,110.709177000000,110.709177000000;112.112579000000,112.112579000000,112.112579000000;113.526305000000,113.526305000000,113.526305000000;114.950375000000,114.950375000000,114.950375000000;116.384811000000,116.384811000000,116.384811000000;117.829635000000,117.829635000000,117.829635000000;119.284868000000,119.284868000000,119.284868000000;120.750532000000,120.750532000000,120.750532000000;122.226647000000,122.226647000000,122.226647000000;123.713235000000,123.713235000000,123.713235000000;125.210317000000,125.210317000000,125.210317000000;126.717914000000,126.717914000000,126.717914000000;128.236047000000,128.236047000000,128.236047000000;129.764737000000,129.764737000000,129.764737000000;131.304005000000,131.304005000000,131.304005000000;132.853871000000,132.853871000000,132.853871000000;134.414357000000,134.414357000000,134.414357000000;135.985483000000,135.985483000000,135.985483000000;137.567270000000,137.567270000000,137.567270000000;139.159738000000,139.159738000000,139.159738000000;140.762907000000,140.762907000000,140.762907000000;142.376799000000,142.376799000000,142.376799000000;144.001434000000,144.001434000000,144.001434000000;145.636832000000,145.636832000000,145.636832000000;147.283012000000,147.283012000000,147.283012000000;148.939997000000,148.939997000000,148.939997000000;150.607804000000,150.607804000000,150.607804000000;152.286456000000,152.286456000000,152.286456000000;153.975971000000,153.975971000000,153.975971000000;155.676371000000,155.676371000000,155.676371000000;157.387673000000,157.387673000000,157.387673000000;159.1099,159.1099,159.1099;160.843070000000,160.843070000000,160.843070000000;162.587203000000,162.587203000000,162.587203000000;164.342319000000,164.342319000000,164.342319000000;166.108438000000,166.108438000000,166.108438000000;167.885578000000,167.885578000000,167.885578000000;169.673761000000,169.673761000000,169.673761000000;171.473005000000,171.473005000000,171.473005000000;173.283330000000,173.283330000000,173.283330000000;175.104755000000,175.104755000000,175.104755000000;176.937299000000,176.937299000000,176.937299000000;178.780982000000,178.780982000000,178.780982000000;180.635824000000,180.635824000000,180.635824000000;182.501843000000,182.501843000000,182.501843000000;184.379058000000,184.379058000000,184.379058000000;186.267489000000,186.267489000000,186.267489000000;188.167154000000,188.167154000000,188.167154000000;190.078073000000,190.078073000000,190.078073000000;192.000265000000,192.000265000000,192.000265000000;193.933749000000,193.933749000000,193.933749000000;195.878543000000,195.878543000000,195.878543000000;197.834666000000,197.834666000000,197.834666000000;199.802137000000,199.802137000000,199.802137000000;201.780975000000,201.780975000000,201.780975000000;203.771198000000,203.771198000000,203.771198000000;205.772826000000,205.772826000000,205.772826000000;207.785876000000,207.785876000000,207.785876000000;209.810367000000,209.810367000000,209.810367000000;211.846319000000,211.846319000000,211.846319000000;213.893748000000,213.893748000000,213.893748000000;215.952674000000,215.952674000000,215.952674000000;218.023115000000,218.023115000000,218.023115000000;220.105089000000,220.105089000000,220.105089000000;222.198615000000,222.198615000000,222.198615000000;224.303711000000,224.303711000000,224.303711000000;226.420395000000,226.420395000000,226.420395000000;228.548685000000,228.548685000000,228.548685000000;230.688599000000,230.688599000000,230.688599000000;232.840156000000,232.840156000000,232.840156000000;235.003373000000,235.003373000000,235.003373000000;237.178269000000,237.178269000000,237.178269000000;239.364861000000,239.364861000000,239.364861000000;241.563167000000,241.563167000000,241.563167000000;243.773205000000,243.773205000000,243.773205000000;245.994993000000,245.994993000000,245.994993000000;248.228549000000,248.228549000000,248.228549000000;250.473890000000,250.473890000000,250.473890000000;252.731035000000,252.731035000000,252.731035000000;255,255,255];

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
