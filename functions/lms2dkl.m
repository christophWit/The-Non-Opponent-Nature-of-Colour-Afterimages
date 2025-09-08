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
