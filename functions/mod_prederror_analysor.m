function RES = mod_prederror_analysor(AGG1a, AGG1b, AGG2a, AGG2b, which_one)
% 2025.07.23 [cw]

header('PREDICTION ERRORS ACROSS EXPERIMENTS',1);

%% SETTINGS
DESIGN.fontsize     = 9;
DESIGN.markersize   = 10;

%% ********************************* HUE **********************************
header('HUE',2);
space = {'ConeCon', 'dkl','Luv', 'Lab', 'cam02', 'mun', 'RYGB'};

%% PREPARE DATA
RES.dev_hue1a_across_sti    = table;
RES.dev_hue1a_across_pp     = table;
RES.dev_hue1b_across_sti    = table;
RES.dev_hue1b_across_pp     = table;
RES.dev_hue2a_across_sti    = table;
RES.dev_hue2a_across_pp     = table;
RES.dev_hue2b_across_sti    = table;
RES.dev_hue2b_across_pp     = table;
for spc = 1:numel(space)
    RES.dev_hue1a_across_sti.(space{spc}) = abs(AGG1a.prederr.hueM.(space{spc})(:,1));
    RES.dev_hue1a_across_pp.(space{spc})  = nanmean(abs(AGG1a.prederr.hue.(space{spc})),2);

    RES.dev_hue1b_across_sti.(space{spc}) = abs(AGG1b.prederr.hueM.(space{spc})(:,1));
    RES.dev_hue1b_across_pp.(space{spc})  = nanmean(abs(AGG1b.prederr.hue.(space{spc})),2);

    RES.dev_hue2a_across_sti.(space{spc}) = nanmean(abs(AGG2a.prederr.hueM.(space{spc})),2);
    RES.dev_hue2a_across_pp.(space{spc})  = nanmean(abs(AGG2a.prederr.hueM.(space{spc})),1)';

    RES.dev_hue2b_across_sti.(space{spc}) = nanmean(abs(AGG2b.prederr.hueM.(space{spc})),2);
    RES.dev_hue2b_across_pp.(space{spc})  = nanmean(abs(AGG2b.prederr.hueM.(space{spc})),1)';
end

% CALCULATE NOISE ESTIMATIONS IN EXP1 BASED ON INTERINDIVIDUAL DIFFERENCES
% The differences from the group mean are calculated:
RES.dev_hue1a_across_sti.Mdev = mean(circularsubtractor(AGG1a.afteri.indi, ones(size(AGG1a.afteri.indi,1),1)*AGG1a.afteri.agg.hue_M', 360, 'mindist'),1)';
RES.dev_hue1b_across_sti.Mdev = nanmean(circularsubtractor(AGG1b.afteri.indi, ones(size(AGG1b.afteri.indi,1),1)*AGG1b.afteri.agg.hue_M', 360, 'mindist'),1)';
RES.dev_hue1a_across_pp.Mdev  = mean(circularsubtractor(AGG1a.afteri.indi, ones(size(AGG1a.afteri.indi,1),1)*AGG1a.afteri.agg.hue_M', 360, 'mindist'),2);
RES.dev_hue1b_across_pp.Mdev  = nanmean(circularsubtractor(AGG1b.afteri.indi, ones(size(AGG1b.afteri.indi,1),1)*AGG1b.afteri.agg.hue_M', 360, 'mindist'),2);

n = size(AGG2a.IndiM.hueM,2);
RES.dev_hue2a_across_sti.Mdev = mean(circularsubtractor(AGG2a.IndiM.hueM, AGG2a.agg.hue(:,1)*ones(1,n),360, 'mindist'),2);
n = size(AGG2b.IndiM.hueM,2);
RES.dev_hue2b_across_sti.Mdev = mean(circularsubtractor(AGG2b.IndiM.hueM, AGG2b.agg.hue(:,1)*ones(1,n),360, 'mindist'),2);

% CALCULATE NOISE ESTIMATIONS IN EXP2 BASED ON INTRAINDIVIDUAL DIFFERENCES
% The differences from each individual's mean are calculated:
tmp = [];
for pp = 1:numel(AGG2a.Indi)
    matches = AGG2a.Indi{pp}(:,:,4);
    matchesM = AGG2a.IndiM.hueM(:,pp);
    n = size(matches,2);
    tmp(:,pp) = mean(circularsubtractor(matches, matchesM*ones(1,n),360, 'mindist'),2);    
end
RES.dev_hue2a_across_sti.Mdev = mean(tmp,2);
RES.dev_hue2a_across_pp.Mdev  = mean(tmp,1)';

tmp = [];
for pp = 1:numel(AGG2b.Indi)
    matches = AGG2b.Indi{pp}(:,:,4);
    matchesM = AGG2b.IndiM.hueM(:,pp);
    n = size(matches,2);
    tmp(:,pp) = mean(circularsubtractor(matches, matchesM*ones(1,n),360, 'mindist'),2);
end
RES.dev_hue2b_across_sti.Mdev = mean(tmp,2);
RES.dev_hue2b_across_pp.Mdev  = mean(tmp,1)';


%% TESTS(HUE)
space = [space, {'Mdev'}];
for spc = 2:numel(space)

    % EXPERIMENT 1a -------------------------------------------------------
    dev = RES.dev_hue1a_across_sti.(space{spc})-RES.dev_hue1a_across_sti.ConeCon;
    RES.ttabs.hue1a_across_sti(spc-1,:) = ttester(dev);
    RES.wtabs.hue1a_across_sti(spc-1,:) = wilcoxer(dev);
    RES.stabs.hue1a_across_sti(spc-1,:) = signtester(dev);
    
    dev = RES.dev_hue1a_across_pp.(space{spc})-RES.dev_hue1a_across_pp.ConeCon;
    RES.ttabs.hue1a_across_pp(spc-1,:) = ttester(dev);
    RES.wtabs.hue1a_across_pp(spc-1,:) = wilcoxer(dev);
    RES.stabs.hue1a_across_pp(spc-1,:) = signtester(dev);

    % EXPERIMENT 1b -------------------------------------------------------
    dev = RES.dev_hue1b_across_sti.(space{spc})-RES.dev_hue1b_across_sti.ConeCon;
    RES.ttabs.hue1b_across_sti(spc-1,:) = ttester(dev);
    RES.wtabs.hue1b_across_sti(spc-1,:) = wilcoxer(dev);
    RES.stabs.hue1b_across_sti(spc-1,:) = signtester(dev);

    dev = RES.dev_hue1b_across_pp.(space{spc})-RES.dev_hue1b_across_pp.ConeCon;
    RES.ttabs.hue1b_across_pp(spc-1,:) = ttester(dev);
    RES.wtabs.hue1b_across_pp(spc-1,:) = wilcoxer(dev);
    RES.stabs.hue1b_across_pp(spc-1,:) = signtester(dev);
    
    % EXPERIMENT 2a -------------------------------------------------------
    dev = RES.dev_hue2a_across_sti.(space{spc})-RES.dev_hue2a_across_sti.ConeCon;
    RES.ttabs.hue2a_across_sti(spc-1,:) = ttester(dev);
    RES.wtabs.hue2a_across_sti(spc-1,:) = wilcoxer(dev);
    RES.stabs.hue2a_across_sti(spc-1,:) = signtester(dev);
    
    dev = RES.dev_hue2a_across_pp.(space{spc})-RES.dev_hue2a_across_pp.ConeCon;
    RES.ttabs.hue2a_across_pp(spc-1,:) = ttester(dev);
    RES.wtabs.hue2a_across_pp(spc-1,:) = wilcoxer(dev);
    RES.stabs.hue2a_across_pp(spc-1,:) = signtester(dev);

    % EXPERIMENT 2b -------------------------------------------------------
    dev = RES.dev_hue2b_across_sti.(space{spc})-RES.dev_hue2b_across_sti.ConeCon;
    dev = round(dev,10);
    RES.ttabs.hue2b_across_sti(spc-1,:) = ttester(dev);
    RES.wtabs.hue2b_across_sti(spc-1,:) = wilcoxer(dev);
    RES.stabs.hue2b_across_sti(spc-1,:) = signtester(dev);

    dev = RES.dev_hue2b_across_pp.(space{spc})-RES.dev_hue2b_across_pp.ConeCon;
    dev = round(dev,10);
    RES.ttabs.hue2b_across_pp(spc-1,:) = ttester(dev);
    RES.wtabs.hue2b_across_pp(spc-1,:) = wilcoxer(dev);
    RES.stabs.hue2b_across_pp(spc-1,:) = signtester(dev);
    
end
name = {'ttabs', 'wtabs','stabs'};
for k = 1:3
    RES.(name{k}).hue1a_across_sti.Properties.RowNames  = space(2:end);
    RES.(name{k}).hue1a_across_pp.Properties.RowNames   = space(2:end);
    RES.(name{k}).hue1b_across_sti.Properties.RowNames  = space(2:end);
    RES.(name{k}).hue1b_across_pp.Properties.RowNames   = space(2:end);
    RES.(name{k}).hue2a_across_sti.Properties.RowNames  = space(2:end);
    RES.(name{k}).hue2a_across_pp.Properties.RowNames   = space(2:end);
    RES.(name{k}).hue2b_across_sti.Properties.RowNames  = space(2:end);
    RES.(name{k}).hue2b_across_pp.Properties.RowNames   = space(2:end);
end

header('EXP1a', 3);

header('Across Colours', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.hue1a_across_sti);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.hue1a_across_sti);

header('Across PPs', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.hue1a_across_pp);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.hue1a_across_pp);

header('EXP1b', 3);

header('Across Colours', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.hue1b_across_sti);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.hue1b_across_sti);

header('Across PPs', 4);
disp(RES.ttabs.hue1b_across_pp);
disp(RES.wtabs.hue1b_across_pp);

header('EXP2a', 3);

header('Across Colours', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.hue2a_across_sti);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.hue2a_across_sti);

header('Across PPs', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.hue2a_across_pp);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.hue2a_across_pp);

header('EXP2b', 3);

header('Across Colours', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.hue2b_across_sti);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.hue2b_across_sti);

header('Across PPs', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.hue2b_across_pp);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.hue2b_across_pp);

%% ******************************** CHROMA ********************************
% There is no data on chroma for Experiment 1

header('CHROMA',2);
space = {'ConeCon', 'dkl','Luv', 'Lab', 'cam02', 'mun', 'RYGB'};

%% PREPARE DATA
RES.dev_chro2a_across_sti    = table;
RES.dev_chro2a_across_pp     = table;
RES.dev_chro2b_across_sti    = table;
RES.dev_chro2b_across_pp     = table;
for spc = 1:numel(space)

    % EXPERIMENT 2a
    RES.dev_chro2a_across_sti.(space{spc}) = nanmean(abs(AGG2a.prederr.chrM.(space{spc})),2);
    RES.dev_chro2a_across_pp.(space{spc})  = nanmean(abs(AGG2a.prederr.chrM.(space{spc})),1)';

    % EXPERIMENT 2b
    RES.dev_chro2b_across_sti.(space{spc}) = nanmean(abs(AGG2b.prederr.chrM.(space{spc})),2);
    RES.dev_chro2b_across_pp.(space{spc})  = nanmean(abs(AGG2b.prederr.chrM.(space{spc})),1)';

end
tmp = [];
for pp = 1:10
    zmatch = zscore(AGG2a.Indi{pp}(:,:,5));
    zmatchM = mean(zmatch,2); % Average across trials
    tmp(:,pp) = mean(abs(bsxfun(@minus, zmatch, zmatchM)),2);
end
RES.dev_chro2a_across_sti.Mdev = mean(tmp,2);
RES.dev_chro2a_across_pp.Mdev  = mean(tmp,1)';

tmp = [];
for pp = 1:2
    zmatch = zscore(AGG2b.Indi{pp}(:,:,5));
    zmatchM = mean(zmatch,2); % Average across trials
    tmp(:,pp) = mean(abs(bsxfun(@minus, zmatch, zmatchM)),2);
end
RES.dev_chro2b_across_sti.Mdev = mean(tmp,2);
RES.dev_chro2b_across_pp.Mdev  = mean(tmp,1)';

%% TESTS (Chroma)

space = [space, {'Mdev'}];
for spc = 2:numel(space)
    
    % EXPERIMENT 2a
    dev = RES.dev_chro2a_across_sti.(space{spc})-RES.dev_chro2a_across_sti.ConeCon;
    RES.ttabs.chro2a_across_sti(spc-1,:) = ttester(dev);
    RES.wtabs.chro2a_across_sti(spc-1,:) = wilcoxer(dev);
    RES.stabs.chro2a_across_sti(spc-1,:) = signtester(dev);
    
    dev = RES.dev_chro2a_across_pp.(space{spc})-RES.dev_chro2a_across_pp.ConeCon;
    RES.ttabs.chro2a_across_pp(spc-1,:) = ttester(dev);
    RES.wtabs.chro2a_across_pp(spc-1,:) = wilcoxer(dev);
    RES.stabs.chro2a_across_pp(spc-1,:) = signtester(dev);

    % EXPERIMENT 2b
    dev = RES.dev_chro2b_across_sti.(space{spc})-RES.dev_chro2b_across_sti.ConeCon;
    RES.ttabs.chro2b_across_sti(spc-1,:) = ttester(dev);
    RES.wtabs.chro2b_across_sti(spc-1,:) = wilcoxer(dev);
    RES.stabs.chro2b_across_sti(spc-1,:) = signtester(dev);

    dev = RES.dev_chro2b_across_pp.(space{spc})-RES.dev_chro2b_across_pp.ConeCon;
    RES.ttabs.chro2b_across_pp(spc-1,:) = ttester(dev);
    RES.wtabs.chro2b_across_pp(spc-1,:) = wilcoxer(dev);
    RES.stabs.chro2b_across_pp(spc-1,:) = signtester(dev);
    
end
name = {'ttabs', 'wtabs', 'stabs'};
for k = 1:3
    RES.(name{k}).chro2a_across_sti.Properties.RowNames  = space(2:end);
    RES.(name{k}).chro2a_across_pp.Properties.RowNames   = space(2:end);
    RES.(name{k}).chro2b_across_sti.Properties.RowNames  = space(2:end);
    RES.(name{k}).chro2b_across_pp.Properties.RowNames   = space(2:end);
end

header('EXP2a', 3);

header('Across Colours', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.chro2a_across_sti);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.chro2a_across_sti);

header('Across PPs', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.chro2a_across_pp);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.chro2a_across_pp);

header('EXP2b', 3);

header('Across Colours', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.chro2b_across_sti);
fprintf('Wilcoxon:\n');
disp(RES.wtabs.chro2b_across_sti);

header('Across PPs', 4);
fprintf('T-Test:\n');
disp(RES.ttabs.chro2b_across_pp);

fprintf('Wilcoxon:\n');
disp(RES.wtabs.chro2b_across_pp);

fprintf('Sign:\n');
disp(RES.stabs.chro2b_across_pp);

%% ******************************* GRAPHIX ********************************
figure('Name', 'Prediction Errors across Experiments', 'NumberTitle', 'off');

ax(1) = subplot(1,3,[1 2]);
h_hue = plot_huebox(RES, 'sti', DESIGN);

ax(2) = subplot(1,3,3);
h_chroma = plot_chromabox(RES, 'sti', DESIGN);

%% ****************************** SUBMODULES ******************************

%% plot_chromabox
function H = plot_chromabox(RES, which_one, DESIGN)
% 2025.07.24 [cw]

% SET UP ------------------------------------------------------------------
if nargin < 2
    which_one = 'sti';
end

rgb = [...
    .5 .5 .5;...
    0  .7 0;...
    1  .2 .2;...
    1  .5 .2;...
    .2 .2 1;...
    .2 .5 1;...
    .7 .2 1;...
    .2 .7 .7];

switch lower(which_one)
    case 'sti'
        dev2a = RES.dev_chro2a_across_sti;
        dev2b = RES.dev_chro2b_across_sti;
    case 'pp'
        dev2a = RES.dev_chro2a_across_pp;
        dev2b = RES.dev_chro2b_across_pp;
end

% PLOT --------------------------------------------------------------------
n2a = size(dev2a,1);
n2b = size(dev2b,1);
x0 = [0 9 18 27];
space = dev2a.Properties.VariableNames;
space = [space(end), space(1:end-1)];
space_lbl = {'Noise', 'Cones', 'DKL', 'CIELUV', 'CIELAB', 'CIECAM02', 'Munsell', 'Hering'};
spc_n = numel(space);
hold on
H = [];
for sp = 1:spc_n
    x = x0+sp;
    chroma = [...
        dev2a.(space{sp});...
        dev2b.(space{sp})];
    X = [...
        ones(n2a,1)*x(1);...
        ones(n2b,1)*x(2)];
    h = boxchart(X, chroma);
    h.BoxFaceColor      = rgb(sp,:);
    h.WhiskerLineColor  = rgb(sp,:);
    h.MarkerColor       = rgb(sp,:);
    h.JitterOutliers    = 'on';
    h.MarkerStyle       = '.';
    h.MarkerSize        = DESIGN.fontsize;
    H = [H;h];
end
set(gca,...
    'FontSize', DESIGN.fontsize,...
    'Xlim', [0 18], 'XTick', 4.5:9:32, 'XTickLabel', {'Exp 2a', 'Exp 2b (DKL)', });
ylabel('Abs. Differences [zscores]');
title('HUE-SPECIFIC CHROMA', 'FontWeight', 'bold');

%% plot_huebox
function H = plot_huebox(RES, which_one, DESIGN)
% 2025.07.24 [cw]

% SET UP ------------------------------------------------------------------
if nargin < 2
    which_one = 'sti';
end

rgb = [...
    .5 .5 .5
    0  .7 0;...
    1  .2 .2;...
    1  .5 .2;...
    .2 .2 1;...
    .2 .5 1;...
    .7 .2 1;...
    .2 .7 .7];

switch lower(which_one)
    case 'sti'
        dev1a = RES.dev_hue1a_across_sti;
        dev1b = RES.dev_hue1b_across_sti;
        dev2a = RES.dev_hue2a_across_sti;
        dev2b = RES.dev_hue2b_across_sti;
    case 'pp'
        dev1a = RES.dev_hue1a_across_pp;
        dev1b = RES.dev_hue1b_across_pp;
        dev2a = RES.dev_hue2a_across_pp;
        dev2b = RES.dev_hue2b_across_pp;
end

% PLOT --------------------------------------------------------------------
n1a = size(dev1a,1);
n1b = size(dev1b,1);
n2a = size(dev2a,1);
n2b = size(dev2b,1);
x0 = [0 9 18 27];
space = dev1a.Properties.VariableNames;
space = [space(end), space(1:end-1)];
space_lbl = {'Noise', 'Cones', 'DKL', 'CIELUV', 'CIELAB', 'CIECAM02', 'Munsell', 'Hering'};
spc_n = numel(space);
hold on
H = [];
for sp = 1:spc_n
    x = x0+sp;
    hue = [...
        dev1a.(space{sp});...
        dev1b.(space{sp});...
        dev2a.(space{sp});...
        dev2b.(space{sp})];
    X = [...
        ones(n1a,1)*x(1);...
        ones(n1b,1)*x(2);...
        ones(n2a,1)*x(3);...
        ones(n2b,1)*x(4)];
    h = boxchart(X, hue);
    h.BoxFaceColor      = rgb(sp,:);
    h.WhiskerLineColor  = rgb(sp,:);
    h.MarkerColor       = rgb(sp,:);
    h.JitterOutliers    = 'on';
    h.MarkerStyle       = '.';
    h.MarkerSize        = DESIGN.fontsize;
    H = [H;h];
end
set(gca,...
    'FontSize', DESIGN.fontsize,...
    'Xlim', [0 36], 'XTick', 4.5:9:32, 'XTickLabel', {'Exp 1a', 'Exp 1b', 'Exp 2a', 'Exp 2b (DKL)', });
ylabel('Abs. Differences [deg]');
legend(H, space_lbl, 'Location', 'EastOutside');
title('HUE', 'FontWeight', 'bold');

%% swarmhue_plotter
function swarmhue_plotter(RES, which_one)
% 2025.07.24 [cw]

if nargin < 2
    which_one = 'sti';
end

mrksz = 10;

switch lower(which_one)
    case 'sti'
        M1a = table2array(RES.dev_hue1a_across_sti);
        M1b = table2array(RES.dev_hue1b_across_sti);
        M2a = table2array(RES.dev_hue2a_across_sti);
        M2b = table2array(RES.dev_hue2a_across_sti);
    case 'pp'
        M1a = table2array(RES.dev_hue1a_across_pp);
        M1b = table2array(RES.dev_hue1b_across_pp);
        M2a = table2array(RES.dev_hue2a_across_pp);
        M2b = table2array(RES.dev_hue2b_across_pp);
end

hold on
swarmchart(ones(8,1),dev, mrksz,AGG1a.stimuli.inducer.rgb, 'filled', 'MarkerEdgeColor', 'k');
swarmchart(ones(8,1)*2,dev, mrksz,AGG1a.stimuli.inducer.rgb, 'filled', 'MarkerEdgeColor', 'k');
dev = abs(AGG1a.prederr.hueM.Luv(:,1));
boxplot(dev, 'Positions',3);
dev = abs(AGG1a.prederr.hueM.Lab(:,1));
boxplot(dev, 'Positions',4);
%swarmchart(ones(8,1)*3,dev, mrksz,AGG1a.stimuli.inducer.rgb, 'filled', 'MarkerEdgeColor', 'k');
dev = abs(AGG1a.prederr.hueM.cam02(:,1));
dev = mean(abs(AGG1a.prederr.hue.cam02),2);
boxplot(dev, 'Positions',5);
%swarmchart(ones(8,1)*4,dev, mrksz,AGG1a.stimuli.inducer.rgb, 'filled', 'MarkerEdgeColor', 'k');
dev = abs(AGG1a.prederr.hueM.mun(:,1));
boxplot(dev, 'Positions',6);
%swarmchart(ones(8,1)*5,dev, mrksz,AGG1a.stimuli.inducer.rgb, 'filled', 'MarkerEdgeColor', 'k');
dev = abs(AGG1a.prederr.hueM.RYGB(:,1));
h = boxplot(dev, 'Positions',7);
set(h(:), 'MarkerEdgeColor', 'g');
%swarmchart(ones(8,1)*6,dev, mrksz,AGG1a.stimuli.inducer.rgb, 'filled', 'MarkerEdgeColor', 'k');

hold off


%% hue_bars
function hue_bars()
% 2025.07.24 [cw]

fsz = 9;
fsz2 = 10;

x = 1:7;
x1 = x;
x2 = x+8;
x3 = x+16;
x4 = x+24;
rgb = [...
    0 0.7 0;...
    1 0.2 0.2;...
    1 .5 0.2;...
    0.2 0.2 1;...
    0.2 0.5 1;...
    0.7 .2 1;...
    0.2 .7 .7];
hold on
for k = 1:7
    h(k) = bar(x1(k),M1(k),'FaceColor', rgb(k,:), 'BarWidth', 1);
    bar(x2(k),M2(k),'FaceColor', rgb(k,:), 'BarWidth', 1);
    bar(x3(k),M3(k),'FaceColor', rgb(k,:), 'BarWidth', 1);
    bar(x4(k),M4(k),'FaceColor', rgb(k,:), 'BarWidth', 1);
end
errorbar(x1,M1,S1, 'k.');
errorbar(x2,M2,S2, 'k.');
errorbar(x3,M3,S3, 'k.');
errorbar(x4,M4,S4, 'k.');
significance_plotter1(x1([1 4]),M1([1,3]),S1([1 3]),ttab.sig(1));
significance_plotter1(x2([1 4]),M2([1 4]),S2([1 4]),ttab.sig(2));
significance_plotter1(x3([1 4]),M3([1 3]),S3([1 3]),ttab.sig(3));
significance_plotter1(x4([1 4]),M4([1 3]),S4([1 3]),ttab.sig(4));
hold off
significance_plotter0(x1(5),M1(5),S1(5),ttab_cam02.sig(1));
significance_plotter0(x2(5),M2(5),S2(5),ttab_cam02.sig(2));
significance_plotter0(x3(5),M3(5),S3(5),ttab_cam02.sig(3));
significance_plotter0(x4(5),M4(5),S4(5),ttab_cam02.sig(4));
set(gca, 'FontSize', fsz,...
    'XTick', mean([x1(:), x2(:), x3(:), x4(:)]), 'XTickLabel', {'Exp1a', 'Exp1b', 'Exp2a', 'Exp2b (DKL)'});
legend(h, spc_lbls, 'Location', 'EastOutside');
ylabel('Hue Differences [deg]');
title('HUE', 'FontWeight', 'bold');

axs(2) = subplot(1,3,3)
inds = 1:7; % With Hering Chroma, which is meaningless 
inds = 1:6; % Exclude Hering Chroma because it's meaningless
hold on
% plot([x1(1)-.5 x1(end)+.5],[1 1]*cM3(3),'k:');
% plot([x2(1)-.5 x2(end)+.5],[1 1]*cM4(2),'k:');
% % EXCLUDE MUNSELL BECAUSE CHROMA DOES NOT MAKE SENSE:
% cM3(6) = NaN;
% cM4(6) = NaN;
for k = 1:numel(inds)
    h2(k) = bar(x1(k),cM3(k),'FaceColor', rgb(k,:), 'BarWidth', 1);
    bar(x2(k),cM4(k),'FaceColor', rgb(k,:), 'BarWidth', 1);
end
errorbar(x1(inds),cM3(inds),cS3(inds), 'k.');
errorbar(x2(inds),cM4(inds),cS4(inds), 'k.');
significance_plotter1(x1([1 end]),cM3([1 2]),cS3([1 2]),ttab_c.sig(1));
significance_plotter1(x2([1 2]),cM4([1 2]),cS4([1 2]),ttab_c.sig(2));
hold off
set(gca, 'FontSize', fsz,...
    'XTick', mean([x1(:), x2(:)]), 'XTickLabel', {'Exp2a','Exp2b (DKL)'});
%legend(h, spc_lbls);
ylabel('Chroma Differences [zscores]');
title('HUE-SPECIFIC CHROMA', 'FontWeight', 'bold');

% REFORMAT ----------------------------------------------------------------
wd1 = .52; wd2 = .24; ht = .75; y = .12; x = [.05 .75];
set(axs(1), 'Position', [x(1) y wd1, ht]);
set(axs(2), 'Position', [x(2) y wd2, ht]);
% title(axs(1:3), {''});
% set(axs(3), 'Position', [x(3) y wd, ht]);
ylim([0 2]);


%% ************************* GENERAL SUBFUNCTIONS *************************

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

%% signtester
function stab = signtester(diffs)
% 2025.07.28 [cw]

n = sum(~isnan(diffs));
[p, h, stats] = signtest(diffs);
z = stats.zval;
sign = stats.sign;
sig = significancer(p,'ns',0);

stab = table(n, sign, z, p, sig);

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

%% wilcoxer
function tab = wilcoxer(data)
% 2021.04.07

if nargin == 0
    data = rand(20,1)-0.5;
end

inds = data>0;
rnk = tiedrank(abs(data));
W(1) = sum(rnk(inds));
W(2) = sum(rnk(~inds));
n = sum(~isnan(data));

[p,h,stats] = signrank(data);
S = stats.signedrank;
if isfield(stats, 'zval')
    z = stats.zval;
else
    z = NaN(size(stats,1),1);
end
drctn = sign(W(1)-W(2));

sig = significancer(p,'ns',0);

tab = table(S, W, n, z, p, sig, drctn);

if nargout == 0
    fprintf('Wilcoxon sign rank test (two-sided): S = %d (W+ = %d; W- = %d; n = %d), p = %1.3f (direction = %d)\n', S, W(1), W(2), n, p, drctn);
end