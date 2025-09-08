function [induced_lms, cone_increments, cone_contrasts] = afterimage_simulator(inducer_lms, bg_lms, weights, adapt_strength, probe_lms)
% 2020.08.18 * based on "simulate_adaptation" [cw]
% 2022.10.18 added input probe_lms if different from background [cw]
% 2023.01.19 added option to leave weights & adapt_strength empty [cw]

if nargin < 5
    if nargin < 4
        adapt_strength = [];
        if nargin < 3
            weights = [1 1 1];
            if nargin < 2
                close all; clc;
                example;
                return;
            end
        end
    end
    probe_lms = bg_lms;
end

if isempty(weights)
    weights = [1 1 1];
end
if isempty(adapt_strength)
    adapt_strength = 1;
end


% WEIGHT CONES (optional) -------------------------------------------------
% depending on (1) completeness vs partial adaptation; (2) relative
% adaptation of cones compared one to another
inducer_lms2 = bsxfun(@times, inducer_lms, weights);

% CALCULATE COLOUR OF PROBE RELATIVE TO INDUCER --------------------------- 
% i.e. taking INDUCER as the point of adaptation.
cone_increments = bsxfun(@minus,probe_lms,inducer_lms2); % cone increments = [L-L0, M-M0, S-S0] (Kaiser & Boynton, 1996, p.564)
cone_contrasts = bsxfun(@rdivide, cone_increments,inducer_lms2); % cone contrasts = [dL/L0, dM/M0, dS/S0]

% INDUCTION: Calculate the induced colour relative to the background ------
% i.e. taking the BG as the point of adaptation.

% ARBITRARY SCALING (ignore, leave at 1)
adapted_cone_contrasts = cone_contrasts*adapt_strength;

% MAIN (induction):
cone_increments = bsxfun(@times, adapted_cone_contrasts,bg_lms);
induced_lms = bsxfun(@plus, cone_increments,bg_lms);

%% ***************************** SUBFUNCTIONS *****************************

%% example
function [inducer_lms, bg_lms] = example
%2020.08.23

% SRGB
mon_xyY = [0.640000000000000,0.330000000000000,17;0.300000000000000,0.600000000000000,57.2000000000000;0.150000000000000,0.0600000000000000,5.80000000000000;0.312700000000000,0.329000000000000,80];

%mon_xyY = [0.684670000000000,0.311100000000000,26.3810000000000;0.213800000000000,0.726300000000000,69.8960000000000;0.152100000000000,0.0453000000000000,4.78210000000000;0.328700000000000,0.350130000000000,100.440000000000];

% MON FOR CHASER EXP:
mon_xyY = [0.646700000000000,0.335400000000000,42.4750000000000;0.312950000000000,0.601700000000000,190.350000000000;0.147950000000000,0.0503500000000000,17.0100000000000;0.295600000000000,0.319650000000000,251.950000000000];


% INUDCERS = circle in CIELUV ---------------------------------------------
azis = (5:5:360)'; 
RAD = 42;
L = 70;

sti_n = size(azis,1);
light = ones(sti_n,1)*L;
chroma = ones(sti_n,1)*RAD;
INDUCER = colourconverter([light azis chroma], 'Luv_pol',2,'mon',mon_xyY);
inducer_lms = INDUCER.lms;
%inducer_lms = XYZ2lms(INDUCER.XYZ,'ss','1931');


% BACKGROUND ISOLUMINANT WITH INDUCERS ------------------------------------
BG = colourconverter([L 0 0], 'Luv',2,'mon',mon_xyY);
bg_lms = BG.lms;
%bg_lms = XYZ2lms(BG.XYZ,'ss','1931');
disp(bg_lms);

induced_lms = afterimage_simulater(inducer_lms, bg_lms);
XYZ = lms2XYZ(induced_lms,'ss','1931');
INDUCED = colourconverter(XYZ, 'xyz',2,'mon',mon_xyY);

figure('Name','EXAMPLE', 'NumberTitle','off');
subplot(1,2,1)
scatter(INDUCED.Luv(:,2),INDUCED.Luv(:,3),40,INDUCED.rgb/255, 'filled','MarkerEdgeColor','k');
axis([-60 60 -100 60]);
axis square;
title('INDUCED COLOURS', 'FontWeight', 'bold');

subplot(1,2,2)
frq0 = histcounts(INDUCED.Luv_pol(:,1),[0,2.5:5:357.5,360]);
frq = frq0(1:end-1)';
frq(1)=frq0(1)+frq0(end);
[x,y] = pol2cart(deg2rad(azis),frq);
h1 = fill(x,y,[.8 .8 .8], 'EdgeColor', [.5 .5 .5]);
axis square;
title('HUE HISTOGRAM', 'FontWeight', 'bold');
