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
