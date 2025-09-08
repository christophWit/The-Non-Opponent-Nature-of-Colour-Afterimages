function panel_identifier(rws, cls, xshift, yshift, fsz, mode, axs, design, start_with)
% 2019.02.18 * [cw]
% 2019.11.15 added design "up"
% 2019.12.07 added "k" in ids
% 2020.09.09 added option "start_with"

if nargin < 9
    start_with = 1;
if nargin < 8
    design = '.)';
if nargin < 7
    axs = [];
if nargin < 6
    mode = 'rows';
    if nargin < 5
        fsz = 10;
        if nargin < 4
            yshift = [0.1, 0.1];            
            if nargin < 3
                xshift = [0.25, 0.1];
            end
        end
    end
end
end
end
end

ids = {'a', 'b', 'c', 'd', 'e', 'f','g','h','i','j', 'k', 'l','m','n', 'o','p','q','r','s','t','u','v','w','x','y','z'};

ids = ids(start_with:end);

if numel(xshift) == 1
    xshift = ones(2,1)*xshift;
end
if numel(yshift) == 1
    yshift = ones(2,1)*yshift;
end

px1 = -xshift(1);
py1 = 1+yshift(1);
px2 = -xshift(2);
py2 = 1+yshift(2);
tot_n = rws*cls;

if isempty(axs)
    axs = get(gcf, 'Children');
    axs = flipud(axs);
end

tot_n = numel(axs);
for nr = 1:tot_n
    axes(axs(nr));
    switch lower(mode)
        case 'columns'
            nr2 = mod(nr,rws)+(mod(nr,rws)==0)*rws + (ceil(nr/rws)-1)*rws;
        otherwise
            nr2 = nr;
    end
    if mod(nr,cls) == 1
        px = px1;
    else
        px = px2;
    end
    if mod(nr,rws) == 1
        py = py1;
    else
        py = py2;
    end
    switch lower(design)
        case '.)'
            text(px, py, sprintf('%s.)', ids{nr2}),...
                'units', 'normalized',...
                'FontWeight', 'Bold', 'FontSize', fsz,...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'Bottom');
        case ')'
            text(px, py, sprintf('%s)', upper(ids{nr2})),...
                'units', 'normalized',...
                'FontWeight', 'Bold', 'FontSize', fsz,...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'Bottom');
        case '()'
            text(px, py, sprintf('(%s)', ids{nr2}),...
                'units', 'normalized',...
                'FontWeight', 'Bold', 'FontSize', fsz,...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'Bottom');
        case 'up'
            text(px, py, sprintf('%s', upper(ids{nr2})),...
                'units', 'normalized',...
                'FontWeight', 'Bold', 'FontSize', fsz,...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'Bottom');
        case '(up)'
            text(px, py, sprintf('(%s)', upper(ids{nr2})),...
                'units', 'normalized',...
                'FontWeight', 'Bold', 'FontSize', fsz,...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'Bottom');
    end
end


