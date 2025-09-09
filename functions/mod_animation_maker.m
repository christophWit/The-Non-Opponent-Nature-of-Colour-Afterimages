function fg = mod_animation_maker(filename_stem, iAzi, iChroma)
% 2025.09.05-07 * [cw]

if nargin < 3
    iChroma = 80;
    if nargin < 2
        iAzi = 60;
        if nargin == 0
            filename_stem = 'test';
        end
    end
end

% SETTINGS ----------------------------------------------------------------
DESIGN.fontsize = 8;
L = 70;
gif_file = sprintf('%s.gif', filename_stem);
png_file = sprintf('%s.png', filename_stem);


INDU = colourconverter([L iAzi iChroma], 'Luv_pol',2);
BG = colourconverter([L 0 0], 'Luv_pol',2);
TMP = colourconverter([L iAzi iChroma/2], 'Luv_pol',2);
sim_lms = afterimage_simulator(TMP.lms, BG.lms);
sim_XYZ = lms2XYZ(sim_lms,'ss','1931');
SIM = colourconverter(sim_XYZ, 'xyz',2,'mon');
CompChroma2 = SIM.Luv_pol(2);
COMP = colourconverter([L iAzi+180 CompChroma2], 'Luv_pol',2);

disp([BG;INDU;COMP; SIM]);

% PRODUCTION --------------------------------------------------------------
azi0 = 0:.5:360;
azi = 0:5:360;
mov.rad = 5;
mrkr_sz = 50;
azi = deg2rad(azi);
[x,y] = pol2cart(azi,mov.rad);
[x0,y0] = pol2cart(azi0,mov.rad);
sti_n = numel(azi);

fg = figure('Color',BG.rgb/255, 'InvertHardCopy','off');
azi = 0:5:360;
azi2 = 0:5:350;
azi = deg2rad(azi);
azi2 = deg2rad(azi2);
[x1,y1] = pol2cart(azi,3*mrkr_sz+1);
[x2,y2] = pol2cart(azi,mrkr_sz);
[pos_x,pos_y] = pol2cart(azi2,mrkr_sz*2);
frame_size = 300;
im_n = numel(azi2);
hold on
fill([1 -1 -1 1]*frame_size, [-1 -1 1 1]*frame_size,BG.rgb/255);
for im = 1:im_n
    r = fill(x1,y1,INDU.rgb/255, 'EdgeColor','k');

    x = x2+pos_x(im);
    y = y2+pos_y(im);
    p = fill(x,y,BG.rgb/255,'EdgeColor','none');
    
    if im > im_n*3/8 & im < im_n*5/8 
        c1 = fill(x2,y2,COMP.rgb/255, 'EdgeColor','k');
        c2 = fill(x2(1:round(end/2)),y2(1:round(end/2)),SIM.rgb/255, 'EdgeColor','k');
    else
        c1 = fill(x2,y2,BG.rgb/255, 'EdgeColor','k');
    end    
    plot(0,0,'k+');
    
    set(gca, ...
        'FontSize', DESIGN.fontsize,...
        'Units', 'pixel', 'Position', [0 0 frame_size+2 frame_size+2]);
    text(0,1,sprintf('  H %d deg / C %d',round(INDU.Luv_pol(1)), round(INDU.Luv_pol(2))),...
        'Units','Normalized', 'FontSize', DESIGN.fontsize,...
        'HorizontalAlignment','Left','VerticalAlignment','Top');
    set(gcf,'Units', 'pixel',...
        'Position',[500 200 frame_size frame_size],...
        'PaperUnits', 'centimeters',...
        'PaperPosition', [0  0 frame_size frame_size]);
    drawnow
    
    % Capture the plot as an image
    frame = getframe(fg);
    IM = frame2im(frame);
    [imind,cm] = rgb2ind(IM,65536);
    
    % Write to the GIF File
    if im == 1
        imwrite(imind,cm,gif_file,'gif', 'Loopcount',inf,'DelayTime',.01);
    else
        imwrite(imind,cm,gif_file,'gif','WriteMode','append','DelayTime',.01);
    end
end

% ADD PREDICTED HUES FOR OVERVIEW -----------------------------------------
c1 = fill(x2,y2,COMP.rgb/255, 'EdgeColor','k');
c2 = fill(x2(1:round(end/2)),y2(1:round(end/2)),SIM.rgb/255, 'EdgeColor','k');
hold off
set(gca, 'Units', 'Normalized');

% Hue angle of cone adaptation (Half-circle on top):
text(0.5,0.5,sprintf('%d deg',round(SIM.Luv_pol(1))),...
    'Units','Normalized','FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','Center','VerticalAlignment','Bottom');

% Hue angle of cone-opponency (Half-circel on bottom):
text(0.5,0.5,sprintf('%d deg',round(COMP.Luv_pol(1))),...
    'Units','Normalized','FontSize', DESIGN.fontsize,...
    'HorizontalAlignment','Center','VerticalAlignment','Top');

%% ***************************** SUBFUNCTIONS *****************************

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
        % A. Stockman, D. I. A. MacLeod, and N. E. Johnson, ‘‘Spectrallms2X
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