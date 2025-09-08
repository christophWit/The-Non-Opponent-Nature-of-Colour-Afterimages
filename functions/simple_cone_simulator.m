function CONES = simple_cone_simulator(stimuli, k)
% Calculates the directions of excitations of each single, isolated cone.
% 2025.07.20

if nargin < 2
    k = 10; % S-Cone factor: produces negative chromaticities if too high
end

colspace    = stimuli.src_space;
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
%XYZ = lms2XYZ(lms, stimuli.sens, stimuli.cmf);
CONES = colourconverter(lms, 'lms', 2,...
    wp, stimuli.sens, stimuli.cmf,...
    stimuli.mon.xyY, stimuli.mon.ldt, stimuli.mon.oog, [], stimuli.mon.rgb_max,...
    stimuli.bg.xyY(3), 0);
