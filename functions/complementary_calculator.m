function [OppoHue, OppoChroma, OppoLum, dkl_k, Munsell] = complementary_calculator(STI, trgt_colspace, CHIPS, HeringRYGB, adapt_chroma)
% --------------------------------- INPUT ---------------------------------
% 01 STI
% = The stimulus structure with all information about stimuli; this
% structure is the same as the one called "stimuli" as part of the
% aggregated data AGG.
% 02 trgt_colspace
% = The space in which complementary colours will be represented; this is
% not the same as the src_colspace, which is included in the stimulus
% structure STI (see 01).
% 03 CHIPS
% = Table with Munsell Chips as prepared for and included in these
% experiment data files.
% 04 HeringRYGB 
% = carries the CIELUV hue directions (azimuth in degree) of the Hering
% colours Red, Yellow, Green, and Blue; must be CIELUV (!), typically from
% Experiment 1a or 1b.
% 05 adapt_chroma 
% = determines the strenght of adaptation; can be set to: 
% 'inducer'     = full adaptation to the inducer
% 'halfinducer' = adaptation to half the inducer chroma (as in Experiment 1 & 3).
% a number      = typically for the empirical grand average (as in Experiment 2).
% -------------------------------- OUTPUT ---------------------------------
% 01 OppoHue
% Complementary Hues according to each space. Note that all those hues are
% represented in the trgt_colspace (CIELUV or DKL).
% 02 OppoChroma
% Complementary chroma according to each space. Note that chroma is also
% represented in the trgt_colspace (CIELUV or DKL).
% 03 OppoLum
% Complementary luminance; this is mainly for control and double-checking
% as it should be constant (apart from small differences when different
% colour matching functions came into play). 
% 04 dkl_k
% The factor of subtractive adaptation used to simulate subtractive
% adaptation in DKL space.
% 05 Munsell
% The CIELAB and MUNSEL coordinates of inducers and complementary colours
% in the Munsell system. This is also mainly for control and
% double-checking.
% -------------------------------- VERSION --------------------------------
% 2025.07.08 polished for publication [cw]

if nargin < 5
    adapt_chroma = 'halfinducer';
    if nargin < 4
        HeringRYGB = [10.6014492753623;72.1533333333333;125.133333333333;227.94];
    end
end

% CALIBRATION -------------------------------------------------------------
src_colspace = STI.src_space;
WP = STI.wp;
wp = WP.xyY;
BG = STI.bg;
bgY = BG.xyY(3);
mon_xyY = STI.mon.xyY;
mon_ldt = STI.mon.ldt;
mon_oog = STI.mon.oog;
rgb_max = STI.mon.rgb_max;
sti_n   = size(STI.inducer,1);
switch lower(src_colspace)
    case 'luv'
        colspace = 'Luv_pol';
        sens = 'ss'; cmf = '1931';
    case 'dkl'
        colspace = 'dkl_pol';
        sens = 'sp'; cmf = 'judd';
end
mon_lms = XYZ2lms(xyY2XYZ(STI.mon.xyY), sens, cmf);

% DETERMINE AFTERIMAGE STRENGTH -------------------------------------------
% Strength is determined by adjusting the chroma of the inducer.
inducer0 = [STI.inducer.(src_colspace)(:,1), STI.inducer.(colspace)]; % Inducer in polar coordinates, including L*
switch adapt_chroma
    case 'inducer'
        chroma = inducer0(:,3);
    case 'halfinducer'
        chroma = inducer0(:,3)/2;
    otherwise
        if numel(adapt_chroma) == 1
            chroma = ones(sti_n, 1)*adapt_chroma;
        elseif numel(adapt_chroma) == sti_n
            chroma = adapt_chroma;
        else
            error('which_chroma4consim does not make sense.');
        end
end
% ADAPTER DIFFERS FROM INDUCER BY CHROMA (as adaptation strength):
adapt = [inducer0(:,1:2),chroma];
% ADAPTER = colourconverter(adapt, colspace,2,...
%     wp,mon_xyY,mon_ldt,mon_oog, [], rgb_max, 1, 0);
ADAPTER = colourconverter(adapt, colspace, 2,...
    wp, sens, cmf,...
    mon_xyY,mon_ldt,mon_oog, [], rgb_max,...
    BG.xyY(3), 0);

%% LMS
%ADAPTER.lms     = XYZ2lms(ADAPTER.XYZ, sens, cmf);
BG.lms          = XYZ2lms(STI.bg.XYZ,  sens, cmf);
STI.inducer.lms = XYZ2lms(STI.inducer.XYZ,  sens, cmf);

% CONE-CONTRAST PREDICTION ------------------------------------------------
induced_lms     = afterimage_simulator(ADAPTER.lms, BG.lms);
XYZ             = lms2XYZ(induced_lms, sens, cmf);
%ConeAdapt       = colourconverter(XYZ, 'xyz',2,wp,mon_xyY,mon_ldt,mon_oog, [], rgb_max, 1, 0);
ConeAdapt       = colourconverter(XYZ, 'xyz', 2,...
    wp, sens, cmf,... 
    mon_xyY,mon_ldt,mon_oog, [], rgb_max,...
    bgY, 0);
%[ConeAdapt.dkl, ConeAdapt.dkl_pol, ConeAdapt.lms] = XYZ2dkl(XYZ, BG.lms, sens, cmf, mon_lms);

%% DKL
dkl0    = lms2dkl(STI.inducer.lms, BG.lms, mon_lms);
dkl0R   = lms2dkl(ADAPTER.lms, BG.lms, mon_lms);
[~, chr0] = cart2pol(dkl0(:,2),dkl0(:,3));
[~, chr0R] = cart2pol(dkl0R(:,2),dkl0R(:,3));
dkl_k   = chr0R./chr0; % Factor k indicating adaptation strength
dkl_k   = mean(dkl_k);
fprintf('Subtractive Adaption Strength (DKL) is: dkl_k = %.2f\n', dkl_k);

% DKL COMPLEMENTARIES -----------------------------------------------------
dkl = BG.dkl-dkl_k*dkl0;
[dkl_azi, dkl_rad] = cart2pol(dkl(:,2),dkl(:,3));
dkl_azi = rad2deg(dkl_azi);
dkl_pol = [dkl_azi(:), dkl_rad(:)];
lms = dkl2lms(dkl, BG.lms, mon_lms);
XYZ = lms2XYZ(lms, sens, cmf);
DKL = colourconverter(XYZ, 'XYZ', 2,...
    wp, sens, cmf,...
    mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
    bgY, 0);
DKL.dkl_pol = dkl_pol;

% CIECAM02 -------------------------------------------------------------
wpY = WP.xyY(3);
La = BG.xyY(3);

% ADAPTATION STRENGTH:
%kill_D = []; % Empty = Normal use of CIECAM02.
kill_D = 1; % This means D is fixed at zero for CIECAM02, no matter the surround setting.
surround = 'average'; % This does not make a difference if D is fixed at zero (kill_D)

% NORMALISE TO 1:
Y_b         = BG.xyY(3)/wpY;
XYZwp       = WP.XYZ/wpY;
bg_norm     = BG.XYZ/wpY;
indu_norm   = ADAPTER.XYZ./wpY;

% CIELUV, CIELAB, CIECAM02, and MUNSELL COMPLEMENTARIES -------------------
sti_n = size(ADAPTER,1);
comp_L = ADAPTER.Lab(:,1);
wpY = WP.xyY(3);
indu_lab = [ADAPTER.Lab(:,1), ADAPTER.Lab_pol];
oppo_lab = NaN(sti_n, 3);
oppo_mun = NaN(sti_n, 3); 
indu_mun = NaN(sti_n, 3);
Munsell = table(indu_lab, indu_mun, oppo_mun, oppo_lab);
for sti = 1:sti_n
    [Munsell.oppo_lab(sti,:), Munsell.oppo_mun(sti,:), Munsell.indu_mun(sti,:)] = munsell_opposer(CHIPS, indu_lab(sti,:));
    ADAPTED(sti,:) = colourconverter(BG.XYZ, 'XYZ', 2,...
        [ADAPTER.xyY(sti,1:2), wpY], sens, cmf,...
        mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
        bgY,0);
    CIECAM02_out(sti,:) = xyz2ciecam02(bg_norm, indu_norm(sti,:), La, Y_b, surround, kill_D);
end
wp_mun = [0.31006, 0.31616, 100]; % illuminant C from https://en.wikipedia.org/wiki/Standard_illuminant

% MUNSELL
MUN = colourconverter(...
    Munsell.oppo_lab, 'lab_pol', 2,...
    WP.xyY, sens, cmf,...
    mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
    bgY, 0);

% LUV, LAB, and CIECAM02 ACCORDING TO ADAPTATION TRANSFORM ----------------

% CIELUV
LUV = colourconverter(ADAPTED.Luv, 'Luv', 2,...
    WP.xyY, sens, cmf,...
    mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
    bgY, 0);

% CIELAB
LAB = colourconverter(ADAPTED.Lab, 'Lab', 2,...
    WP.xyY, sens, cmf,...
    mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
    bgY, 0);

% CIECAM02
XYZ = ciecam022xyz(CIECAM02_out, XYZwp, La, Y_b, surround, kill_D);
%XYZ = XYZ*wpY;
XYZ = XYZ*bgY;
CAM02 = colourconverter(XYZ, 'xyz', 2,...
    WP.xyY, sens, cmf,...
    mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
    bgY, 0);

% HSV ---------------------------------------------------------------------
rgb = ADAPTER.rgbl/rgb_max;
rgb(rgb < 0) = 0; % Artifacts
rgb(rgb > 1) = 1; % Artifacts
hsv = rgb2hsv(rgb);
%hsv = rgb2hsl(rgb); % same hue as hsv!
hsv(:,1) = mod(hsv(:,1)-0.5,1);
rgb2 = hsv2rgb(hsv)*rgb_max;
HSV = colourconverter(rgb2, 'rgbl', 2,...
    wp, sens, cmf,...
    mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
    bgY, 0);

% HERING ------------------------------------------------------------------
azi_luv = mod(ADAPTER.Luv_pol(:,1),360);
Hhue = hering_afterimage_simulater(HeringRYGB, azi_luv);
Hluv = [comp_L,Hhue,ADAPTER.Luv_pol(:,2)];
Hering = colourconverter(Hluv, 'luv_pol', 2,...
    wp, sens, cmf,...
    mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
    bgY, 0);

switch lower(src_colspace)
    case 'dkl' % Constant chroma in DKL instead of Luv
        dkl_pol(:,1) = Hering.dkl_pol(:,1);
        dkl_pol(:,2) = ones(size(Hering,1),1)*adapt_chroma;
        Hering = colourconverter([Hering.dkl(:,1), dkl_pol], 'dkl_pol', 2,...
            wp, sens, cmf,...
            mon_xyY, mon_ldt, mon_oog, [], rgb_max,...
            bgY, 0);
end

% SELECT TARGET SPACE -----------------------------------------------------

switch lower(trgt_colspace)
    case 'luv'
        colspace = 'Luv_pol';
    case 'dkl'
        colspace = 'dkl_pol';
    case 'lab'
        colspace = 'Lab_pol';
end

% HUE
ConeCon = ConeAdapt.(colspace)(:,1);
Luv     = LUV.(colspace)(:,1);
dkl     = DKL.(colspace)(:,1);
Lab     = LAB.(colspace)(:,1);
mun     = MUN.(colspace)(:,1);
hsv     = HSV.(colspace)(:,1);
RYGB    = Hering.(colspace)(:,1);
cam02   = CAM02.(colspace)(:,1);

OppoHue = table(ConeCon, Luv, dkl, hsv, Lab, cam02, mun, RYGB);

% CHROMA
ConeCon = ConeAdapt.(colspace)(:,2);
Luv     = LUV.(colspace)(:,2);
dkl     = DKL.(colspace)(:,2);
Lab     = LAB.(colspace)(:,2);
mun     = MUN.(colspace)(:,2);
hsv     = HSV.(colspace)(:,2);
RYGB    = Hering.(colspace)(:,2);
cam02   = CAM02.(colspace)(:,2);

OppoChroma = table(ConeCon, Luv, dkl, hsv, Lab, cam02, mun, RYGB);

% Luminance for Control 
ConeCon = ConeAdapt.xyY(:,3);
Luv     = LUV.xyY(:,3);
dkl     = DKL.xyY(:,3);
Lab     = LAB.xyY(:,3);
mun     = MUN.xyY(:,3);
hsv     = HSV.xyY(:,3);
RYGB    = Hering.xyY(:,3);
cam02   = CAM02.xyY(:,3);

OppoLum = table(ConeCon, Luv, dkl, hsv, Lab, cam02, mun, RYGB);

%% ***************************** SUBFUNCTIONS *****************************

%% ciecam022xyz
function XYZ = ciecam022xyz(ciecam02, XYZ_w, L_A, Y_b, surround, kill_D)
%ciecam022xyz Convert CIECAM02 values to XYZ.
%
%   out = ciecam022xyz(in,XYZ_w,L_A,Y_b,surround) converts the table of
%   CIECAM02 values to XYZ (CIE 1931 2-Degree Standard Colorimetric
%   Observer) values. in should be a table with variables J, Q, C, M, s, h,
%   and H.
%
%   XYZ_w is a three-element vector containing the CIE XYZ values
%   for the adopted white point.
%
%   L_A is the adapting luminance (in cd/m^2).
%
%   Y_b is the relative background luminance (in the range
%   [0.0,1.0].
%
%   surround is either 'average' (the typical relative luminance for
%   viewing reflection prints), 'dim' (the typical relative
%   luminance for CRT displays or televisions), or 'dark' (the
%   typical relative luminance for projected transparencies).
%
%   REFERENCE
%
%   Mark D. Fairchild, Color Appearance Models, 3rd edition, John
%   Wiley & Sons, 2013, pp. 287-302.
%
%   EXAMPLE
%
%   Convert an XYZ color to CIECAM02. (From Fairchild, Table 16.4,
%   Case 1, p. 299)
%  
%  Version 2.1.2.0 (151 KB) by Steve Eddins https://uk.mathworks.com/matlabcentral/fileexchange/64161-color-tools-for-matlab
%
%       J = 41.73;
%       Q = 195.37;
%       C = 0.10;
%       M = 0.11;
%       s = 2.36;
%       h = 219.0;
%       H = 278.1;
%       in = table(J,Q,C,M,s,h,H);
%       XYZ_w = [0.9505 1.0000 1.0888];
%       L_A = 318.31;
%       Y_b = 0.20;
%       surround = 'average';
%       out = ciecam022xyz(in,XYZ_w,L_A,Y_b,surround)
%
%   See also xyz2ciecam02.
%   Copyright MathWorks 2016-2018
% 2025.04.26 added kill_D to set D directly to 1. [cw]

if nargin < 6
    kill_D = [];
end

XYZ_w = XYZ_w * 100;
Y_w = XYZ_w(2);
Y_b = Y_b * 100;
% Calculate A_w using equations 7.1 through 7.28. CIE 159:2004, p. 9.
n = Y_b / Y_w;
z = 1.48 + sqrt(n);
switch surround
   case 'average'
      c = 0.69;
      N_c = 1.0;
      F = 1.0;
      
   case 'dim'
      c = 0.59;
      N_c = 0.9;
      F = 0.9;
      
   case 'dark'
      c = 0.525;
      N_c = 0.8;
      F = 0.8;
end
D = F*(1 - (1/3.6)*exp(-(L_A + 42)/92));

% CW 2025.04.26
if ~isempty(kill_D)
    D = 1;
end

N_bb = 0.725*(1/n)^0.2;
N_cb = N_bb;
k = 1/(5*L_A + 1);
F_L = (0.2 * k^4 * 5 * L_A) + 0.1*(1-k^4)^2 * (5*L_A)^(1/3);
M_cat02 = [ ...
   0.7328  0.4296 -0.1624
   -0.7036  1.6975  0.0061
   0.0030  0.0136  0.9834];
M_HPE = [ ...
    0.38971  0.68898  -0.07868
    -0.22981 1.18340   0.04641
    0.00000  0.00000   1.00000 ];
RGB_w = XYZ_w * M_cat02';
RGB_w_c = ((Y_w * D ./ RGB_w) + (1 - D)) .* RGB_w;
RGB_w_p = (M_HPE * (M_cat02 \ RGB_w_c'))';
f = @(RGB) sign(RGB) .* (400 * (F_L * abs(RGB) / 100).^0.42) ./ (27.13 + (F_L * abs(RGB) / 100).^0.42) + 0.1;
RGB_w_ap = f(RGB_w_p);
A_w = (2*RGB_w_ap(1) + RGB_w_ap(2) + (1/20)*RGB_w_ap(3) - 0.305) * N_bb;
% Compute t from C and J.
% Equation 8.6, CIE 159:2004, p. 10.
%
C = ciecam02.C;
J = ciecam02.J;
t = (C ./ (sqrt(J/100) * (1.64 - 0.29.^n).^0.73)) .^ (1/0.9);
% Equation 8.5, CIE 159:2004, p. 10, and Table 3, p. 8.
if ismember('H',ciecam02.Properties.VariableNames) && ...
        ~ismember('h',ciecam02.Properties.VariableNames)
    H = ciecam02.H;
    i = min(floor(H/100)+1,4);
    
    h_i = [20.14 90.00 164.25 237.53 380.14]';
    e_i = [0.8 0.7 1.0 1.2 0.8]';
    H_i = [0 100 200 300 400]';
    
    num = (H - H_i(i)) .* (e_i(i+1).*h_i(i) - e_i(i).*h_i(i+1)) - 100*h_i(i).*e_i(i+1);
    den = (H - H_i(i)) .* (e_i(i+1) - e_i(i)) - 100*e_i(i+1);
    h_p = num ./ den;
    
    h = h_p;
    h(h>360) = h(h>360) - 360;
else
    h = ciecam02.h;
end
% Equation 8.7, CIE 159:2004, p. 10.
e_t = 0.25 * (cos(h*pi/180 + 2) + 3.8);
% Equation 8.8, CIE 159:2004, p. 10.
A = A_w * (J/100).^(1/(c*z));
% Equation 8.9, CIE 159:2004, p. 10.
p_1 = (50000/13) * N_c * N_cb * e_t ./ t;
% Equation 8.10, CIE 159:2004, p. 10.
p_2 = (A/N_bb) + 0.305;
% Equation 8.11, CIE 159:2004, p. 10.
p_3 = 21/20;
% Equation 8.12, CIE 159:2004, p. 10.
h_r = h*(pi/180);
% Equations 8.13 - 8.18, CIE 159:2004, pp. 10-11.
a = zeros(size(h));
b = zeros(size(h));
sin_greater = abs(sin(h_r)) >= abs(cos(h_r));
% Compute a and b values for sin_greater
p_4 = p_1 ./ sin(h_r);
b_1_num = p_2 .* (2 + p_3) * 460 / 1403;
b_1_den = p_4 + (2 + p_3) * (220/1403) * ...
    (cos(h_r)./sin(h_r)) - (27/1403) + p_3*(6300/1403);
b_1 = b_1_num ./ b_1_den;
a_1 = b_1 .* (cos(h_r)./sin(h_r));
% Compute a and b values for ~sin_greater
p_5 = p_1 ./ cos(h_r);
a_2_num = p_2 .* (2 + p_3) * (460/1403);
a_2_den = p_5 + (2 + p_3)*(220/1403) - ((27/1403) - p_3*(6300/1403)) .* (sin(h_r)./cos(h_r));
a_2 = a_2_num ./ a_2_den;
b_2 = a_2 .* (sin(h_r)./cos(h_r));
a(sin_greater) = a_1(sin_greater);
b(sin_greater) = b_1(sin_greater);
a(~sin_greater) = a_2(~sin_greater);
b(~sin_greater) = b_2(~sin_greater);
% "Note that if t is equal to zero, then set a and b equal to zero,
% calculate A using equation 8.8, compute p_2 using equation 8.10 and go
% directly to equation 8.19." CIE 159:2004, p. 10.
t_equal_0 = (t == 0);
a(t_equal_0) = 0;
b(t_equal_0) = 0;
% Equations 8.19 - 8.21, CIE 159:2004, p. 11.
M = [ ...
    460    451    288
    460   -891   -261
    460   -220   -6300
    ] / 1403;
RGB_a_p = [p_2 a b] * M';
% Equations 8.22 - 8.24, CIE 159:2004, p. 11.
RGB_p = (100/F_L) * ( (27.13 * abs(RGB_a_p - 0.1)) ./ (400 - abs(RGB_a_p - 0.1)) ).^(1/0.42);
% "If any of the values of (Ra' - 0.1), (Ga' - 0.1) or (Ba' - 0.1) are
% negative, then the corresponding value R', G' or B' must be made
% negative."
RGB_p = RGB_p .* sign(RGB_a_p - 0.1);
% Equation 8.25, CIE 159:2004, p. 11.
RGB_c = (M_cat02 * (M_HPE \ RGB_p'))';
% Equation 8.27 - 8.29, CIE 159:2004, p. 11.
RGB = RGB_c ./ ((Y_w * D ./ RGB_w) + 1 - D);
% Equation 8.30, CIE 159:2004, p. 11.
XYZ = (M_cat02 \ RGB')';
XYZ = XYZ / 100;

%% cinterpolator
function [newy, newx] = cinterpolator(X1, Y1, X2, onecycle)
% 2012sep24-25 * [cw]

if nargin < 4
    onecycle= 360;    
end

% check input:
if size(X1,1) ~= size(Y1,1) | size(X1,2) ~= size(Y1,2)
    error('cinterpolator: x1 and y1 don''t match!');
end

% Delete NaNs:
inds1 = ~isnan(X1);
inds2 = ~isnan(Y1);
inds = inds1&inds2;
X1 = X1(inds);
Y1 = Y1(inds);

% Map to cycle:
xy(:,1) = mod(X1,onecycle);
xy(:,2) = Y1;
xy = sortrows(xy,1);
x1 = xy(:,1);
y1 = xy(:,2);
x2 = mod(X2, onecycle);
%x2 = sort(x2);

% Eliminate multiple x-values by averaging the corresponding y-values:
[x1 y1] = double_killer(x1,y1);
% [x2 dbl_n] = double_killer(x2);
% if dbl_n > 0
%     fprintf('WARNING - cinterpolator: interpolation involved double target x2; deleted doubles!\n');
% end

% Extend x and y so as to obtain a full cycle:
xex = [...
    x1-onecycle;...
    x1;...
    x1+onecycle];
yex = [y1; y1; y1];
% Actual interpolation:
y2 = interp1(xex, yex, x2, 'linear');
newx = x2;
newy = y2;

%% double_killer
function [newx newy dbl_n] = double_killer(x,y)
%2012sep24 [cw]
newx = x;
qy = 0;
if nargin > 1
    qy = 1;
    newy = y;
end
tab = tabulate(newx);
inds = find(tab(:,2) > 1);
dbls = tab(inds,1);
dbl_n = numel(dbls);
for nr = 1:dbl_n
    inds = newx == dbls(nr);
    newx(inds(2:end)) = [];
    if qy
        newy(inds(1)) = mean(y1(inds));
        newy(inds(2:end)) = [];
    end
end

if ~qy
    newy = dbl_n;
end

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

%% dkl2lms
function [lms, M] = dkl2lms(dkl, lms0, monlms)
% Inverts lms2dkl by inverting the matrix.
% 2025.04.23

L0 = lms0(1);
M0 = lms0(2);
S0 = lms0(3);

% MATRIX THAT SCALES RELATIVE TO WP (L0, M0, S0)
M = [...
    1/(L0+M0),  1/(L0+M0),      0;...
    1/L0,       -1/M0,      0;...
    -1/(L0+M0), -1/(L0+M0), 1/S0];

if nargin > 2
    [maxL, maxS] = map2mongamut(M, lms0, monlms);
    M(2,:) = M(2,:)/maxL;
    M(3,:) = M(3,:)/maxS;
end

dlms = M \ dkl';
dlms = dlms';

% INCREMENTS TO LMS:
lms = dlms+ones(size(dlms,1),1)*lms0;
lms = round(lms,12);

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

%% lms2XYZ
function XYZ = lms2XYZ(lms, fndmtls, cmf)
% 2014dec14 * [cw]
% NOT TESTED YET!!!!!

if nargin < 3
    cmf = 'judd';
    if nargin < 2
        fndmtls = 'smithpokorny';
    end
end

M = coneconversionmatrix(fndmtls, cmf);
XYZ = lms/(M'); %apparently more precise than: XYZ = lms*inv(M')

%% map2mongamut
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

% The below code, going through lms signal first, is equivalent to direct
% conversion of dlms:
% pole_lms = pole_dlms+ones(size(pole_dlms,1),1)*lms0;
% pole_lms = round(pole_lms,12);
% pole_rgb= pole_lms/(monlms(1:3,1:3));
% pole_drgb = pole_rgb-rgb0

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

%% hering_afterimage_simulater
function induced_azi = hering_afterimage_simulater(RYGB_azi, inducer_azi)
% 2022.09.25 * 

inducer_azi = mod(inducer_azi,360);
R = RYGB_azi(1,:);
Y = RYGB_azi(2,:);
G = RYGB_azi(3,:);
B = RYGB_azi(4,:);

% RED-YELLOW --------------------------------------------------------------
width1 = circularsubtractor(Y,R,360,'mindist'); % RED-YELLOW
width2 = circularsubtractor(B,G,360,'mindist'); % OPPONENT: GREEN-BLUE
inds = inducer_azi > R & inducer_azi <= Y;
rel_indu = circularsubtractor(inducer_azi(inds),R,360,'signed')./width1;
induced_azi(inds,1) = G + (rel_indu*width2);

% YELLOW-GREEN ------------------------------------------------------------
width1 = circularsubtractor(G,Y,360,'mindist'); % YELLOW-GREEN
width2 = circularsubtractor(R,B,360,'mindist'); % OPPONENT: BLUE-RED
inds = inducer_azi > Y & inducer_azi <= G;
rel_indu = circularsubtractor(inducer_azi(inds),Y,360,'signed')./width1;
induced_azi(inds,1) = B + (rel_indu*width2);

% GREEN-BLUE --------------------------------------------------------------
width1 = circularsubtractor(B,G,360,'mindist'); % GREEN-BLUE
width2 = circularsubtractor(Y,R,360,'mindist'); % OPPONENT: RED-YELLOW
inds = inducer_azi > G & inducer_azi <= B;
rel_indu = circularsubtractor(inducer_azi(inds),G,360,'signed')./width1;
induced_azi(inds,1) = R + (rel_indu*width2);

% BLUE-RED ----------------------------------------------------------------
width1 = circularsubtractor(R,B,360,'mindist'); % BLUE-RED
width2 = circularsubtractor(G,Y,360,'mindist'); % OPPONENT: YELLOW-GREEN
inds1 = inducer_azi > B & inducer_azi <= 360;
inds2 = inducer_azi >= 0 & inducer_azi <= R;
inds = inds1 | inds2;
rel_indu = circularsubtractor(inducer_azi(inds),B,360,'signed')./width1;
induced_azi(inds,1) = Y + (rel_indu*width2);


%% munsell_opposer
function [oppo_lab, oppo_mun, MunInducer] = munsell_opposer(CHIPS, lab_lhc)
% 2025.04.15 * [cw]

lab_L   = lab_lhc(:,1);
lab_hue = lab_lhc(:,2);
lab_chr = lab_lhc(:,3);

% DELETE CHIP WITH ZERO CHROMA --------------------------------------------
CHIPS(CHIPS.MChr == 0,:) = [];

% ADD POLAR COORDINATES TO TABLE WITH CHIPS -------------------------------
[azi_lab_mun, chr_mun]  = cart2pol(CHIPS.Lab(:,2),CHIPS.Lab(:,3));
CHIPS.Lab_pol = [rad2deg(azi_lab_mun), chr_mun];

% MUNSELL HUES ------------------------------------------------------------
all_hues = {...
    '2.50R','5.00R','7.50R','10.00R',...
    '2.50YR','5.00YR','7.50YR','10.00YR',...
    '2.50Y','5.00Y','7.50Y','10.00Y',...
    '2.50GY', '5.00GY', '7.50GY', '10.00GY',...
    '2.50G', '5.00G', '7.50G', '10.00G',...
    '2.50BG', '5.00BG', '7.50BG', '10.00BG',...
    '2.50B', '5.00B','7.50B','10.00B',...
    '2.50PB','5.00PB','7.50PB','10.00PB',...
    '2.50P','5.00P','7.50P','10.00P',...
    '2.50RP', '5.00RP', '7.50RP', '.00R'}; 
munhue_n = numel(all_hues);
munhue_nr = 1:munhue_n;
munhue_deg = munhue_nr/munhue_n*360; 
CHIPS.Mhue_nr = cellfun(@(x) find(strcmpi(all_hues,x)),CHIPS.Mhue);
CHIPS.Mhue_deg = mod(CHIPS.Mhue_nr/munhue_n*360,360);
[CHIPS.Mx, CHIPS.My] = pol2cart(deg2rad(CHIPS.Mhue_deg),CHIPS.MChr);

% TANSLATE INDUCER LIGHTNESS TO MUNSELL LIGHTNESS -------------------------
MunVal = (1:9)';
L = unique(CHIPS.Lab(:,1));
iMunVal = interp1(L, MunVal, lab_lhc(1), 'linear', 'extrap');
lowL = floor(iMunVal);
upL = ceil(iMunVal);
indsL1 = CHIPS.MVal == lowL;
indsL2 = CHIPS.MVal == upL;

% TRANSLATE INDUCER HUE INTO MUNSELL HUE ----------------------------------
if lab_lhc(3) < max(CHIPS.Lab_pol(indsL1,2)) & lab_lhc(3) < max(CHIPS.Lab_pol(indsL2,2))
    indsX = abs(CHIPS.Lab_pol(:,2) - lab_lhc(3)) <= 10; % Reduce range of candiate hues
else % If too saturated, choose maximum for each hue
    tmp1 = grpstats(CHIPS.MChr(indsL1), CHIPS.Mhue(indsL1), 'max');
    tmp2 = grpstats(CHIPS.MChr(indsL2), CHIPS.Mhue(indsL2), 'max');
    minimax = min([tmp1;tmp2]);
    indsX = CHIPS.MChr == minimax; 
end

% LOWER LIGHTNESS:
if isempty(CHIPS(indsL1&indsX,:))
    keyboard
end

iMunNr1 = munsell_opposition_helper1(CHIPS(indsL1&indsX,:), lab_hue);

% UPPER LIGHTNESS:
if upL < 10
    iMunNr2 = munsell_opposition_helper1(CHIPS(indsL2&indsX,:), lab_hue);
else
    iMunNr2 = [];
end

% INTERPOLATE HUE ACROSS LIGHTNESS:
iMunNrs = [iMunNr1, iMunNr2];
if abs(iMunNr1-iMunNr2) > 20
    [~,ind] = max(iMunNrs);
    iMunNrs(ind) = iMunNrs(ind)-40;
end
if upL < 10
    iMunNr = interp1([lowL, upL], iMunNrs, iMunVal);
else
    iMunNr = iMunNrs;
end
iMunNr = mod(iMunNr,40);
iMunDeg = iMunNr/40*360;

% TANSLATE INDUCER LAB CHROMA TO MUNSELL CHROMA -------------------------------
lowH    = floor(iMunNr);
upH     = ceil(iMunNr);
indsH1 = mod(CHIPS.Mhue_nr,40) == lowH;
iMunChr1_1 = interp1(CHIPS.Lab_pol(indsL1&indsH1,2), CHIPS.MChr(indsL1&indsH1,1), lab_lhc(3),'linear','extrap');
if upL < 10
    iMunChr1_2  = interp1(CHIPS.Lab_pol(indsL2&indsH1,2), CHIPS.MChr(indsL2&indsH1,1), lab_lhc(3),'linear','extrap');
    iMunChr1    = interp1([lowL, upL], [iMunChr1_1, iMunChr1_2], iMunVal,'linear','extrap');
    % if isnan([iMunChr1_1, iMunChr1_2])
    %     keyboard
    % end
else
    iMunChr1 = iMunChr1_1; 
end

if mod(lowH,40) == mod(upH,40) | upL == 10
    iMunChr = iMunChr1;
else % Two Candidate Munsell Hues
    indsH2 = mod(CHIPS.Mhue_nr,40) == mod(upH,40);
    iMunChr2_1 = interp1(CHIPS.Lab_pol(indsL1&indsH2,2), CHIPS.MChr(indsL1&indsH2,1), lab_lhc(3),'linear','extrap');
    iMunChr2_2 = interp1(CHIPS.Lab_pol(indsL2&indsH2,2), CHIPS.MChr(indsL2&indsH2,1), lab_lhc(3),'linear','extrap');
    iMunChr2 = interp1([lowL, upL], [iMunChr2_1, iMunChr2_2], iMunVal, 'linear','extrap');

    iMunChr = cinterpolator([lowH, upH], [iMunChr1, iMunChr2], iMunNr, 40);
end

MunInducer = [iMunVal, iMunDeg, iMunChr];

% COMPLEMENTARY MUNSELL COLOUR --------------------------------------------
oMunDeg  = mod(iMunDeg - 180, 360); 
oMunNr = mod(oMunDeg/360*40,40);
oppo_mun = [iMunVal, oMunDeg, iMunChr];

% CONVERSION FROM MUNSELL TO CIELAB ---------------------------------------

% OPPONENT LAB HUE --------------------------------------------------------
mini = min(abs(CHIPS.MChr-iMunChr));
indsC0 = abs(CHIPS.MChr-iMunChr)==mini;

% LOW LIGHTNESS:
oMunLab_hue1  = munsell_opposition_helper2(CHIPS(indsL1&indsC0,:), oMunNr);

% HIGH LIGHTNESS:
if upL < 10
    oMunLab_hue2  = munsell_opposition_helper2(CHIPS(indsL2&indsC0,:), oMunNr);
    
    % INTERPOLATE HUE ACROSS LIGHTNESS:
    oMunLab_hues = [oMunLab_hue1, oMunLab_hue2];
    if abs(oMunLab_hue1-oMunLab_hue2) > 180
        [~,ind] = max(oMunLab_hues);
        oMunLab_hues(ind) = oMunLab_hues(ind)-360;
    end
    oMunLab_hue = interp1([lowL, upL], oMunLab_hues, iMunVal);
else
    oMunLab_hue = oMunLab_hue1;
end
oMunLab_hue  = mod(oMunLab_hue,360);

% OPPONENT CHROMA ---------------------------------------------------------
lowMH = floor(oMunNr);
upMH = ceil(oMunNr);
indsH1 = mod(CHIPS.Mhue_nr,40) == lowMH;
indsH2 = CHIPS.Mhue_nr == upMH;
Mhues = [lowMH, upMH];
if abs(Mhues-Mhues) > 20
    [~,ind] = max(Mhues);
    Mhues(ind) = Mhues(ind)-20;
end

% LOW LIGHTNESS; LOW HUE
indsC = indsL1 & indsH1;
MunCs = CHIPS.MChr(indsC,1);
LabCs = CHIPS.Lab_pol(indsC,2);
if sum(indsC) == 0
    keyboard
end
LabC1_1 = interp1(MunCs, LabCs, iMunChr, 'linear', 'extrap');

% LOW LIGHTNESS; HIGH HUE
indsC = indsL1 & indsH2;
MunCs = CHIPS.MChr(indsC,1);
LabCs = CHIPS.Lab_pol(indsC,2);
LabC1_2 = interp1(MunCs, LabCs, iMunChr, 'linear', 'extrap');

% INTEGRATE ACROSS HUE
LabC1 = interp1(Mhues, [LabC1_1 LabC1_2], oMunNr, 'linear', 'extrap');

if upL < 10
    % HIGH LIGHTNESS; LOW HUE
    indsC = indsL2 & indsH1;
    MunCs = CHIPS.MChr(indsC,1);
    LabCs = CHIPS.Lab_pol(indsC,2);
    LabC2_1 = interp1(MunCs, LabCs, iMunChr, 'linear', 'extrap');
    
    % HIGH LIGHTNESS; HIGH HUE
    indsC = indsL2 & indsH2;
    MunCs = CHIPS.MChr(indsC,1);
    LabCs = CHIPS.Lab_pol(indsC,2);
    LabC2_2 = interp1(MunCs, LabCs, iMunChr, 'linear', 'extrap');
    
    % INTEGRATE ACROSS HUE
    LabC2 = interp1(Mhues, [LabC2_1 LabC2_2], oMunNr, 'linear', 'extrap');
end

% INTEGRATE ACROSS lightness
if upL < 10
    LabC = interp1([lowL, upL], [LabC1 LabC2], iMunVal, 'linear', 'extrap');
else
    LabC = LabC1;
end

% MAIN OUTPUT
oppo_lab = [lab_lhc(1), oMunLab_hue, LabC];

%% munsell_opposition_helper1
function mun_hue = munsell_opposition_helper1(CHIPS, lab_hue)
% 2025.04.16 * [cw]

dH = abs(circularsubtractor(CHIPS.Lab_pol(:,1),lab_hue, 360, 'mindist'));
[tmp,sort_inds] = sort(dH,'ascend');
TMP = CHIPS(sort_inds,:);
TMP.rank = (1:size(TMP,1))';
AGG = grpstats(TMP,'Mhue',"mean");
AGG = sortrows(AGG,'mean_rank', 'ascend');
[~, ind] = min(AGG.mean_rank);
if ind ~= 1 % sanity check;
    keyboard
end
munhue1 = AGG.mean_Mhue_nr(ind);
munhue2 = AGG.mean_Mhue_nr(2);
labhue1 = AGG.mean_Lab_pol(ind,1);
labhue2 = AGG.mean_Lab_pol(2,1);

munhues = [munhue1,munhue2];
if abs(munhue1-munhue2) > 20
    [~,ind] = max(munhues);
    munhues(ind) = munhues(ind)-40;
end
%munhues(munhues==0) = 40;

mun_hue = cinterpolator([labhue1, labhue2], [munhues(:,1), munhues(:,2)], lab_hue, 40);

%% munsell_opposition_helper2
function lab_hue = munsell_opposition_helper2(CHIPS, mun_hue)
% 2025.04.16 * [cw]

dH = abs(circularsubtractor(CHIPS.Mhue_nr,mun_hue, 40, 'mindist'));
[tmp,sort_inds] = sort(dH,'ascend');
TMP = CHIPS(sort_inds,:);
TMP.rank = (1:size(TMP,1))';
AGG = grpstats(TMP,'Mhue',"mean");
AGG = sortrows(AGG,'mean_rank', 'ascend');
[~, ind] = min(AGG.mean_rank);
if ind ~= 1 % sanity check;
    keyboard
end
munhue1 = AGG.mean_Mhue_nr(ind);
munhue2 = AGG.mean_Mhue_nr(2);
labhue1 = AGG.mean_Lab_pol(ind,1);
labhue2 = AGG.mean_Lab_pol(2,1);


labhues = [labhue1,labhue2];
if abs(labhue1-labhue2) > 180
    [~,ind] = max(labhues);
    labhues(ind) = labhues(ind)-360;
end
lab_hue = cinterpolator([munhue1, munhue2], [labhues(:,1), labhues(:,2)], mun_hue, 360);
lab_hue = mod(lab_hue,360);

%% xyz2ciecam02
function out = xyz2ciecam02(XYZ, XYZ_w, L_A, Y_b, surround, kill_D)
%xyz2ciecam02 Convert XYZ values to CIECAM02.
%
%   out = xyz2ciecam02(XYZ,XYZ_w,L_A,Y_b,surround) converts the Mx3
%   matrix of CIE XYZ values to a table of CIECAM02 values. XYZ is
%   assumed to contain CIE 1931 Standard Colorimetric Observer (2
%   degree) values in the range [0.0,1.0].
%
%   XYZ_w is a three-element vector containing the CIE XYZ values
%   for the adopted white point.
%
%   L_A is the adapting luminance (in cd/m^2).
%
%   Y_b is the relative background luminance (in the range
%   [0.0,1.0].
%
%   surround is either 'average' (the typical relative luminance for
%   viewing reflection prints), 'dim' (the typical relative
%   luminance for CRT displays or televisions), or 'dark' (the
%   typical relative luminance for projected transparencies).
%
%   REFERENCE
%
%   Mark D. Fairchild, Color Appearance Models, 3rd edition, John
%   Wiley & Sons, 2013, pp. 287-302.
%
%   EXAMPLE
%
%   Convert an XYZ color to CIECAM02. (From Fairchild, Table 16.4,
%   Case 1, p. 299)
%
%       XYZ = [0.1901 0.2000 0.2178];
%       XYZ_w = [0.9505 1.0000 1.0888];
%       L_A = 318.31;
%       Y_b = 0.20;
%       surround = 'average';
%       out = xyz2ciecam02(XYZ,XYZ_w,L_A,Y_b,surround)
%
%   See also xyz2ciecam02.
%   Copyright MathWorks 2016-2018
% Image Processing Toolbox uses tristimulus values in the range
% [0,1]. The CIECAM02 conversion assumes tristimulus values in the
% range [0,100].
%  Version 2.1.2.0 (151 KB) by Steve Eddins https://uk.mathworks.com/matlabcentral/fileexchange/64161-color-tools-for-matlab
% 2025.04.26 added kill_D to set D directly to 1. [cw]

if nargin < 6
    kill_D = [];
end

XYZ = 100 * XYZ;
XYZ_w = 100 * XYZ_w;
Y_b = 100 * Y_b;

% Fairchild, Color Appearance Models, 3rd ed., Wiley, 2013.
% Equation (16.1), Fairchild, p. 290.
M_cat02 = [ ...
    0.7328  0.4296 -0.1624
   -0.7036  1.6975  0.0061
    0.0030  0.0136  0.9834];
% Equation (16.2), Fairchild, p. 290.
RGB = XYZ * M_cat02';
% "The transformation must also be completed for the tristimulus
% values of the adapting stimulus." Fairchild, p. 290.
RGB_w = XYZ_w * M_cat02';
% Table 16.1, Fairchild, p. 289.
switch surround
   case 'average'
      c = 0.69;
      N_c = 1.0;
      F = 1.0;
      
   case 'dim'
      c = 0.59;
      N_c = 0.9;
      F = 0.9;
      
   case 'dark'
      c = 0.525;
      N_c = 0.8;
      F = 0.8;
end
% Equation (16.3), Fairchild, p. 290.
D = F*(1 - (1/3.6)*exp(-(L_A + 42)/92));

% CW 2025.04.26
if ~isempty(kill_D)
    D = 1;
end

% Equations (16.4) - (16.6), Fairchild, p. 290-291.
% Equations (7.4) - (7.6), CIE 159:2004,
% R_w = RGB_w(1);
% G_w = RGB_w(2);
% B_w = RGB_w(3);
% 
% R = RGB(:,1);
% G = RGB(:,2);
% B = RGB(:,3);
Y_w = XYZ_w(2);
% RGB_c = ((Y_w * D ./ RGB_w) + (1 - D)) .* RGB;
RGB_c = bsxfun(@times,((Y_w * D ./ RGB_w) + (1 - D)),RGB);
% R_c = ((100*D/R_w) + (1-D))*R;
% G_c = ((100*D/G_w) + (1-D))*G;
% B_c = ((100*D/B_w) + (1-D))*B;
% RGB_c = [R_c G_c B_c];
RGB_w_c = ((Y_w * D ./ RGB_w) + (1 - D)) .* RGB_w;
% R_w_c = ((100*D/R_w) + (1-D))*R_w;
% G_w_c = ((100*D/G_w) + (1-D))*G_w;
% B_w_c = ((100*D/B_w) + (1-D))*B_w;
% RGB_w_c = [R_w_c G_w_c B_w_c];
% Equation (16.7), Fairchild, p. 292.
k = 1/(5*L_A + 1);
% Equation (16.8), Fairchild, p. 292.
F_L = (0.2 * k^4 * 5 * L_A) + 0.1*(1-k^4)^2 * (5*L_A)^(1/3);
% Equation (16.9), Fairchild, p. 292.
n = Y_b / Y_w;
% Equation (16.10), Fairchild, p. 292.
N_bb = 0.725*(1/n)^0.2;
N_cb = N_bb;
% Equation (16.11), Fairchild, p. 292.
z = 1.48 + sqrt(n);
% Equation (16.13), Fairchild, p. 293.
M_HPE = [ ...
    0.38971  0.68898  -0.07868
    -0.22981 1.18340   0.04641
    0.00000  0.00000   1.00000 ];
% Equation (16.12), Fairchild, p. 292.
RGB_p = (M_HPE * (M_cat02 \ RGB_c'))';
RGB_w_p = (M_HPE * (M_cat02 \ RGB_w_c'))';
% Equations (16.15) - (16.17), Fairchild, p. 292.
% "If any of the values of R', G', or B' are negative, then their absolute
% values must be used, and then the quotient term in Equations 7.15 through
% 7.17 must be multiplied by a negative 1 before adding the value 0.1" CIE
% 159:2004, p. 8.
f = @(RGB) sign(RGB) .* (400 * (F_L * abs(RGB) / 100).^0.42) ./ (27.13 + (F_L * abs(RGB) / 100).^0.42) + 0.1;
RGB_ap = f(RGB_p);
RGB_w_ap = f(RGB_w_p);
% RGB_ap = (400 * (F_L * RGB_p / 100).^0.42) ./ (27.13 + (F_L * RGB_p / 100).^0.42) + 0.1;
% RGB_w_ap = (400 * (F_L * RGB_w_p / 100).^0.42) ./ (27.13 + (F_L * RGB_w_p / 100).^0.42) + 0.1;
% Initial opponent-type responses. Equations (16.18) and (16.19), Fairchild,
% p. 294.
R_ap = RGB_ap(:,1);
G_ap = RGB_ap(:,2);
B_ap = RGB_ap(:,3);
a = R_ap - 12*G_ap/11 + B_ap/11;
b = (R_ap + G_ap - 2*B_ap)/9;
% Hue angle, h, "expressed in degrees ranging from 0 to 360 measured from
% the positive a axis." Equation (16.20), Fairchild, p. 294.
h = atan2(b,a);
h(h<0) = h(h<0) + 2*pi;
h = h * (180/pi);
% Eccentricity factor, equation (16.21), Fairchild, p. 294.
e_t = (1/4) * (cos(h*pi/180 + 2) + 3.8);
% Convert from hue angle to hue quadrature. Equation (16.22) and Table 16.2,
% Fairchild, p. 294.
h_p = h;
add360 = h_p < 20.14;
h_p(add360) = h_p(add360) + 360;
i = (h_p >= 20.14) + ...
    (h_p >= 90) + ...
    (h_p >= 164.25) + ...
    (h_p >= 237.53);
h_i = [20.14 90.00 164.25 237.53 380.14]';
e_i = [0.8 0.7 1.0 1.2 0.8]';
H_i = [0 100 200 300 400]';
H = H_i(i) + (100*(h_p - h_i(i))./e_i(i)) ./ ...
    ((h_p - h_i(i))./e_i(i) + (h_i(i+1) - h_p)./e_i(i+1));
% Lightness, equations (16.23) and (16.24), Fairchild,
% p. 295.
A = (2*R_ap + G_ap + (1/20)*B_ap - 0.305) * N_bb;
A_w = (2*RGB_w_ap(1) + RGB_w_ap(2) + (1/20)*RGB_w_ap(3) - 0.305) * N_bb;
J = 100 * (A / A_w).^(c*z);
if ~isreal(J)
    keyboard
end
% Brightness, equation (16.25), Fairchild, p. 295.
Q = (4/c) * sqrt(J/100) * (A_w + 4) * F_L^0.25;
% Chroma, equations (16.26) and (16.27), Fairchild, pp. 295-296.
t = (50000/13) * N_c * N_cb * e_t .* sqrt(a.^2 + b.^2);
t = t ./ (R_ap + G_ap + (21/20)*B_ap);
C = t.^0.9 .* sqrt(J/100) .* (1.64 - 0.29^n).^0.73;
% Colorfulness, equation (16.28), Fairchild, p. 296.
M = C .* F_L^0.25;
% Saturation, equation (16.29), Fairchild, p. 296.
s = 100 * sqrt(M ./ Q);
out = table(J,Q,C,M,s,h,H);

%% xyY2XYZ
function XYZ = xyY2XYZ(xyY, dim)
%input: x, y, Y values in a row vector or matrix )
% 2012.01.20 [cw]

if nargin < 2
    dim = 2;
end

if dim == 2
    x = xyY(:,1); y = xyY(:,2); Y = xyY(:,3);
elseif dim == 3
    x = xyY(:,:,1); y = xyY(:,:,2); Y = xyY(:,:,3); 
end

X = (Y./y) .* x;
Y = Y;
Z = (Y./y) .* (1-y-x);

if dim == 2
    XYZ = [X Y Z];
elseif dim == 3
    XYZ = cat(3, X, Y, Z);    
end

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