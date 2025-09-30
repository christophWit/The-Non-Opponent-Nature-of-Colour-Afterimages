function figfile = fig_saver(folder, flnm, sz, saveformat, punits, res, printmode)
% 2016.09.10 * [cw]
% 2016.09.30 added printmode [cw]
% 2019.08.10 added svg, tif and pdf format [cw]
% 2021.01.06 added jpg format [cw]
% 2023.06.02 added output figfile 
% 2025.09.30 use fullfile to create paths for Apple/Win compatibility [cw] 

if nargin < 7
    printmode = '-painters';
    if nargin < 6
        res = 300;
        if nargin < 5
            punits = 'centimeters';
            if nargin < 4
                saveformat = 'png';
            end
        end
    end
end

res = sprintf('-r%d', res);

set(gcf, 'PaperUnits', punits, 'PaperPosition', [0, 0, sz(1), sz(2)], 'PaperSize', sz);
figfile0 = fullfile(folder, flnm);

switch lower(saveformat)
    case {'all', 'png'}
        figfile = [figfile0, '.png'];
        print('-dpng', res, figfile, printmode);
end
switch lower(saveformat)
    case {'eps'}
        figfile = [figfile0, '.eps'];
        print('-depsc2', '-tiff', res, figfile, printmode);
    case {'svg'}
        figfile = [figfile0, '.svg'];
        print('-dsvg', figfile, printmode);
    case {'tiff', 'tif'}
        figfile = [figfile0, '.tif'];
        print('-dtiff', res, figfile, printmode);
    case {'jpg'}
        figfile = [figfile0, '.jpg'];
        print('-djpeg', figfile, printmode);
    case {'all', 'pdf'}
        figfile = [figfile0, '.pdf'];
        print('-dpdf', figfile, printmode);
end
fprintf ('SAVE: %s\n', figfile); % Displaying the output path

