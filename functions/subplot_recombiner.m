function [fh, ax] = subplot_recombiner(fhs0, sbp_n1, sbp_n2, ax0, legs)
% INPUT:
% fhs0 = original figures;
% sbp_n1 = numbert of vertical subplots of new figure
% sbp_n2 = numbert of horizontal subplots of new figure
% ax0 = Axes in original figure;
% legs = figure handle for legend
% OUTPUT:
% fh = handle of new figure;
% ax = handles of axes in new figure;
% 2015.05.29 * [cw]
% 2016.03.16 change loop (gr_n) to enable filling less than all subplots.
% 2016.04.14 added legend option.
% 2017.02.18 started to add doubleax option.
% 2020.08.14 added posibility to skip a panel by entering NaN in fhs0 [cw].
% 2025.08.17 Tidied up

% DEFAULT SETTINGS --------------------------------------------------------
if nargin < 5
    legs = [];
    if nargin < 4
        ax0(1:numel(fhs0)) = NaN;
    end
end

fh = figure('Name', 'Combined');
gr_n = numel(fhs0);
sbp_n = sbp_n1*sbp_n2;
ax = zeros(sbp_n,1);
for i = 1:gr_n
    if ~isnan(fhs0(i))
        ax(i)=subplot(sbp_n1,sbp_n2,i);
    end
end

for i = 1:gr_n
    if isnan(fhs0(i))
        continue;
    end
    figure(fhs0(i));
    h = get(fhs0(i),'Children');
    if isnan(ax0(i))
        ax0(i) = 1;
    end
%     sbp_ids = deal2doublax(ch);
%     sbp_ids = ax0(i) ;
    ax_ind = numel(h)+1-ax0(i);
    newh = copyobj(h(ax_ind),fh);
    for j = 1:length(newh)
        posnewh = get(newh(j),'Position');
        possub = get(ax(i),'Position');
        set(newh(j),'Position',...
            [posnewh(1) possub(2) posnewh(3) possub(4)]);
        set(newh(j),'Position',possub);
    end
    delete(ax(i));
    
end

% LEGENDS -----------------------------------------------------------------
if ~isempty(legs)
    for lg = numel(legs)
        if sum(fhs0 == legs(lg)) > 0
            legs0 = findobj(fhs0(legs(lg)), 'Tag', 'legend');
            legpos = get(legs0,'Position');
            legh = copyobj(legs0,fh);
            set(legh, 'Position', legpos);
        else 
            error('legend spec');
        end
    end
    uistack(legh,'top');
end

figure(fh);
ax = get(gcf, 'children');
ax = ax(end:-1:1);

%% SUBPLOTS
function [sbp_ids, indsa] = deal2doublax(ch)

pos = get(ch(:), 'Position'); 
[~, indsa, sbp_ids] = unique(cell2mat(pos), 'rows');

