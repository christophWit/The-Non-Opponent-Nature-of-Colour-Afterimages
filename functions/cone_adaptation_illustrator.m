function fh = cone_adaptation_illustrator
% 2025.08.07 [cw]

DESIGN.fontsize     = 8;
DESIGN.fontsize2    = 9;
DESIGN.axmax        = 100;

%% ILLUSTRATE AFTERIMAGE EXTREMES (Fig S9)
close all;
BG = colourconverter([1 1 1]*190, 'rgb',2,'mon');
bg = BG.lms;

% EXAMPLE COLOURS ---------------------------------------------------------
C = (0:100)';
einzer = ones(size(C,1),1);
L = einzer*BG.Luv(1);
H1 = einzer*60;
H2 = einzer*240;
EX1 = colourconverter([L, H1, C], 'Luv_pol', 2, 'mon');
EX2 = colourconverter([L, H2, C], 'Luv_pol', 2, 'mon');
[induced_EX1.lms] = afterimage_simulator(EX1.lms, bg);
[induced_EX2.lms] = afterimage_simulator(EX2.lms, bg);
induced_EX1  = colourconverter(induced_EX1.lms, 'lms',2,'mon');
induced_EX2  = colourconverter(induced_EX2.lms, 'lms',2,'mon');
rgb_ex1 = EX1.rgb(end,:)/255;
rgb_ex2 = EX2.rgb(end,:)/255;

% INCREASING SINGLE INDUCER CONE EXCITATIONS ------------------------------
inducer0 = (1:0.01:2)';
einzer = ones(size(inducer0));
inducer_L = [inducer0*bg(1), einzer*bg(2), einzer*bg(3)];
inducer_L  = colourconverter(inducer_L, 'lms',2,'mon');

inducer_M = [einzer*bg(1), inducer0*bg(2), einzer*bg(3)];
inducer_M  = colourconverter(inducer_M, 'lms',2,'mon');

inducer_S = [einzer*bg(1), einzer*bg(2), inducer0*bg(3)];
inducer_S  = colourconverter(inducer_S, 'lms',2,'mon');

induced_L = afterimage_simulator(inducer_L.lms, bg);
induced_L   = colourconverter(induced_L, 'lms',2,'mon');

induced_M = afterimage_simulator(inducer_M.lms, bg);
induced_M   = colourconverter(induced_M, 'lms',2,'mon');

induced_S = afterimage_simulator(inducer_S.lms, bg);
induced_S   = colourconverter(induced_S, 'lms',2,'mon');

rgb = [...
    1 0 0;...
    0 .7 0;...
    0 0 1];

fh = figure;

% INCREASING SINGLE CONE EXCITATIONS --------------------------------------
subplot(2,3,1)
hold on
plot(inducer_L.lms(:,1), induced_L.lms(:,1),'-', 'Color', rgb(1,:));
plot(inducer_M.lms(:,2), induced_M.lms(:,2),'-', 'Color', rgb(2,:));
plot(inducer_S.lms(:,3), induced_S.lms(:,3),'-', 'Color', rgb(3,:));
plot(inducer_L.lms(1,1), induced_L.lms(1,1),'.', 'Color', rgb(1,:));
plot(inducer_M.lms(1,2), induced_M.lms(1,2),'.', 'Color', rgb(2,:));
plot(inducer_S.lms(1,3), induced_S.lms(1,3),'.', 'Color', rgb(3,:));
text(inducer_L.lms(1,1), induced_L.lms(1,1),' BG ',...
    'Color', rgb(1,:),...
    'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','Center', 'VerticalAlignment', 'Bottom');
text(inducer_M.lms(1,2), induced_M.lms(1,2),' BG ',...
    'Color', rgb(2,:),...
    'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','Center', 'VerticalAlignment', 'Bottom');
text(inducer_S.lms(1,3), induced_S.lms(1,3),' BG ',...
    'Color', rgb(3,:),...
    'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','left', 'VerticalAlignment', 'Middle');
hold off
set(gca, 'FontSize', DESIGN.fontsize);
text(inducer_L.lms(end,1), induced_L.lms(end,1), 'L',...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'Top',...
    'FontSize', DESIGN.fontsize,...
    'Color', rgb(1,:));
text(inducer_M.lms(end,2), induced_M.lms(end,2), 'M',...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'Top',...
    'FontSize', DESIGN.fontsize,...
    'Color', rgb(2,:));
text(inducer_S.lms(end,3), induced_S.lms(end,3), 'S',...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'Top',...
    'FontSize', DESIGN.fontsize,...
    'Color', rgb(3,:));
%xlabel('Inducer Cone Excitation');
ylabel('Adapted Cone Excitation');
axis square;
title('CONE EXCITATIONS', 'FontWeight', 'Bold');
text(-.5, .5, 'INCREASING',...
    'Units', 'Normalized', 'Rotation', 90,...
    'FontSize', DESIGN.fontsize2, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

subplot(2,3,2)
hold on
% FOR LEGEND
h0 = plot(inducer_L.dkl(:,2),inducer_L.dkl(:,3),'k--');
h1 = plot(induced_L.dkl(:,2),induced_L.dkl(:,3),'k-');

plot(EX1.dkl(:,2), EX1.dkl(:,3),'--', 'Color', rgb_ex1, 'LineWidth',.5);
plot(EX2.dkl(:,2), EX2.dkl(:,3),'--', 'Color', rgb_ex2, 'LineWidth',.5);
plot(induced_EX1.dkl(:,2), induced_EX1.dkl(:,3),'-', 'Color', rgb_ex1, 'LineWidth',.5);
plot(induced_EX2.dkl(:,2), induced_EX2.dkl(:,3),'-', 'Color', rgb_ex2, 'LineWidth',.5);

plot(inducer_L.dkl(:,2),inducer_L.dkl(:,3),'--', 'Color', rgb(1,:));
plot(inducer_M.dkl(:,2),inducer_M.dkl(:,3),'--', 'Color', rgb(2,:));
plot(inducer_S.dkl(:,2),inducer_S.dkl(:,3),'--', 'Color', rgb(3,:));
plot(induced_L.dkl(:,2),induced_L.dkl(:,3),'-', 'Color', rgb(1,:));
plot(induced_M.dkl(:,2),induced_M.dkl(:,3),'-', 'Color', rgb(2,:));
plot(induced_S.dkl(:,2),induced_S.dkl(:,3),'-', 'Color', rgb(3,:));
hold off
set(gca, 'FontSize', DESIGN.fontsize);
ylabel('(L+M)-S');
axis([-1 1 -1 1]);
axis square;
title('DKL', 'FontWeight', 'Bold');
legend([h0;h1], {'inducer', 'induced'}, 'Location', 'SouthWest', 'Orientation', 'Horizontal');

subplot(2,3,3)
hold on
plot(EX1.Luv(:,2), EX1.Luv(:,3),'--', 'Color', rgb_ex1, 'LineWidth',.5);
plot(EX2.Luv(:,2), EX2.Luv(:,3),'--', 'Color', rgb_ex2, 'LineWidth',.5);
plot(induced_EX1.Luv(:,2), induced_EX1.Luv(:,3),'-', 'Color', rgb_ex1, 'LineWidth',.5);
plot(induced_EX2.Luv(:,2), induced_EX2.Luv(:,3),'-', 'Color', rgb_ex2, 'LineWidth',.5);

plot(inducer_L.Luv(:,2), inducer_L.Luv(:,3),'--', 'Color', rgb(1,:));
plot(inducer_M.Luv(:,2), inducer_M.Luv(:,3),'--', 'Color', rgb(2,:));
plot(inducer_S.Luv(:,2), inducer_S.Luv(:,3),'--', 'Color', rgb(3,:));
plot(induced_L.Luv(:,2), induced_L.Luv(:,3),'-', 'Color', rgb(1,:));
plot(induced_M.Luv(:,2), induced_M.Luv(:,3),'-', 'Color', rgb(2,:));
plot(induced_S.Luv(:,2), induced_S.Luv(:,3),'-', 'Color', rgb(3,:));
hold off
set(gca, 'FontSize', DESIGN.fontsize);
axis([-1 1 -1 1]*DESIGN.axmax);
axis square;
ylabel('v*');
title('CIELUV', 'FontWeight', 'Bold');


%% ADAPTING TWO CONES BY DECREASING EXCITATIONS OF ONE CONE

inducer0 = (.1:0.01:1)';
einzer = ones(size(inducer0));

inducerLM = [inducer0*bg(1), einzer*bg(2),   einzer*bg(3)];
inducerLM = colourconverter(inducerLM, 'lms',2,'mon');
inducerMS = [einzer*bg(1),   inducer0*bg(2), einzer*bg(3)];
inducerMS = colourconverter(inducerMS, 'lms',2,'mon');
inducerSL = [einzer*bg(1),   einzer*bg(2),   inducer0*bg(3)];
inducerSL = colourconverter(inducerSL, 'lms',2,'mon');

induced_LM = afterimage_simulator(inducerLM.lms, bg);
induced_LM = colourconverter(induced_LM, 'lms', 2,'mon');
induced_MS = afterimage_simulator(inducerMS.lms, bg);
induced_MS = colourconverter(induced_MS, 'lms', 2,'mon');
induced_SL = afterimage_simulator(inducerSL.lms, bg);
induced_SL = colourconverter(induced_SL, 'lms', 2,'mon');

subplot(2,3,4)
hold on
plot(inducerLM.lms(:,1), induced_LM.lms(:,1),'-', 'Color', rgb(1,:));
plot(inducerMS.lms(:,2), induced_MS.lms(:,2),'-', 'Color', rgb(2,:));
plot(inducerSL.lms(:,3), induced_SL.lms(:,3),'-', 'Color', rgb(3,:));
plot(inducerLM.lms(end,1), induced_LM.lms(end,1),'.', 'Color', rgb(1,:));
plot(inducerMS.lms(end,2), induced_MS.lms(end,2),'.', 'Color', rgb(2,:));
plot(inducerSL.lms(end,3), induced_SL.lms(end,3),'.', 'Color', rgb(3,:));
hold off
set(gca, 'FontSize', DESIGN.fontsize);
xlabel('Inducer Cone Excitation');
ylabel('Adapted Cone Excitation');
axis square;
text(-.5, .5, 'DECREASING',...
    'Units', 'Normalized', 'Rotation', 90,...
    'FontSize', DESIGN.fontsize2, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'Middle');
text(inducerLM.lms(end,1), induced_LM.lms(end,1),' BG ',...
    'Color', rgb(1,:),...
    'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','center', 'VerticalAlignment', 'Bottom');
text(inducerMS.lms(end,2), induced_MS.lms(end,2),' BG ',...
    'Color', rgb(2,:),...
    'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','center', 'VerticalAlignment', 'Bottom');
text(inducerSL.lms(end,3), induced_SL.lms(end,3),' BG ',...
    'Color', rgb(3,:),...
    'FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','center', 'VerticalAlignment', 'Bottom');

subplot(2,3,5)
hold on
plot(EX1.dkl(:,2), EX1.dkl(:,3),'--', 'Color', rgb_ex1, 'LineWidth',.5);
plot(EX2.dkl(:,2), EX2.dkl(:,3),'--', 'Color', rgb_ex2, 'LineWidth',.5);
plot(induced_EX1.dkl(:,2), induced_EX1.dkl(:,3),'-', 'Color', rgb_ex1, 'LineWidth',.5);
plot(induced_EX2.dkl(:,2), induced_EX2.dkl(:,3),'-', 'Color', rgb_ex2, 'LineWidth',.5);

plot(inducerLM.dkl(:,2), inducerLM.dkl(:,3),'--', 'Color', rgb(1,:));
plot(inducerMS.dkl(:,2), inducerMS.dkl(:,3),'--', 'Color', rgb(2,:));
plot(inducerSL.dkl(:,2), inducerSL.dkl(:,3),'--', 'Color', rgb(3,:));
plot(induced_LM.dkl(:,2), induced_LM.dkl(:,3),'-', 'Color', rgb(1,:));
plot(induced_MS.dkl(:,2), induced_MS.dkl(:,3),'-', 'Color', rgb(2,:));
plot(induced_SL.dkl(:,2), induced_SL.dkl(:,3),'-', 'Color', rgb(3,:));
hold off
set(gca, 'FontSize', DESIGN.fontsize);
xlabel('L-M');
ylabel('(L+M)-S');
axis([-1 1 -1 1]);
axis square;

subplot(2,3,6)
hold on
plot(EX1.Luv(:,2), EX1.Luv(:,3),'--', 'Color', rgb_ex1, 'LineWidth',.5);
plot(EX2.Luv(:,2), EX2.Luv(:,3),'--', 'Color', rgb_ex2, 'LineWidth',.5);
plot(induced_EX1.Luv(:,2), induced_EX1.Luv(:,3),'-', 'Color', rgb_ex1, 'LineWidth',.5);
plot(induced_EX2.Luv(:,2), induced_EX2.Luv(:,3),'-', 'Color', rgb_ex2, 'LineWidth',.5);

plot(inducerLM.Luv(:,2), inducerLM.Luv(:,3),'--', 'Color', rgb(1,:));
plot(inducerMS.Luv(:,2), inducerMS.Luv(:,3),'--', 'Color', rgb(2,:));
plot(inducerSL.Luv(:,2), inducerSL.Luv(:,3),'--', 'Color', rgb(3,:));
plot(induced_LM.Luv(:,2), induced_LM.Luv(:,3),'-', 'Color', rgb(1,:));
plot(induced_MS.Luv(:,2), induced_MS.Luv(:,3),'-', 'Color', rgb(2,:));
plot(induced_SL.Luv(:,2), induced_SL.Luv(:,3),'-', 'Color', rgb(3,:));
hold off
set(gca, 'FontSize', DESIGN.fontsize);
xlabel('u*');
ylabel('v*');
axis([-1 1 -1 1]*DESIGN.axmax);
axis square;

axs = get(gcf,'Children');
axs = flipud(axs);
