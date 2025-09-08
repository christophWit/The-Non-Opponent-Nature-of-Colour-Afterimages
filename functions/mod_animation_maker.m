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