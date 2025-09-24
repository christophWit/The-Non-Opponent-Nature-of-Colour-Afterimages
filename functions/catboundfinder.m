function [BOUNDARIES, BOUND_IDS] = catboundfinder(TARGETDATA, varargin)
%__________________________________________________________________________
%
%                   FUNCTION: CATEGORY BOUNDARY FINDER
%__________________________________________________________________________
%OBJECTIVE
% Find boundaries of categorical data.
%--------------------------------------------------------------------------
%INPUT
% TARGETDATA = Categorical data; format: integer, where each integer is a
% category ID; rows = cases; columns = measurement repetitions.
% OPTIONAL INPUT (varargin): 
% May be specified in parameter mode ('keyword', value) or structural mode
% (SETTINGS.keyword = value); keywords:
% 'orientationdata' = Scale of reference; per default, indexes of
% TARGETDATA are used as orientation data.
% 'onecycle' = One complete period of the periodic data, e.g. 360 for
% degrees; default: 360.
% 'cat_ids' = ensemble of possible category ids; necessary, if not all
% relevant categories are present in a given data set, but should be taken
% into accountm when structuring the responses; default = existing category
% ids in the TARGETDATA.
% 'gaps' = Define whether there may be gaps or not between adjacent
% categories; in case of 'no gaps', 2 adjacent boundaries are linearly
% interpolated with respect of categorical responses in order to result
% in one common boundary; format: 0 [no gaps] vs. 1 [gaps]; default = 0,
% i.e. no gaps.
% 'single'    = In order to allow categories with just one entry; set 0 to
% prohibit such categories (default).
% 'criterion' = not implemented 
%--------------------------------------------------------------------------
%OUTPUT
% BOUNDARIES  = The boundaries; format: 
% values: correpond to orientation data; if no orientation data specified,
% BOUNDARIES correspond to BOUND_IDS (see below). 
% rows: correspond to the category ids; if no category ids specified, they
% correspond to the integers in TARGETDATA. 
% columns: 1. col = lower boundaries; 2. col = upper boundaries.
% BOUND_IDS   = Row indices of the category boundaries within TARGETDATA.
%--------------------------------------------------------------------------
%EXAMPLE
% Run without input.
%--------------------------------------------------------------------------
%GENEALOGY
% 2008.09.02 - 2010.06.17 [cw]
% 2025.09.24 polished for publication [cw]
%__________________________________________________________________________

%% ****************************** PRELMINARIES ****************************
%% EXAMPLE
if nargin == 0
    [TARGETDATA, ORIENTATIONDATA, CAT_IDS] = example;
    [BOUNDARIES, BOUND_IDS] = catboundfinder(...
        TARGETDATA,...
        'orientationdata', ORIENTATIONDATA,...
        'cat_ids', CAT_IDS,...
        'plot_it', 1);
    return
end
%% DEFAULTS
SETTINGS.onecycle           = 360;
SETTINGS.fullcycle          = 1; % Consider targetdata to fill the whole cycle.
SETTINGS.orientationdata    = [];
SETTINGS.cat_ids            = [];
SETTINGS.criterion          = 0.5; % not implemented 
SETTINGS.gaps               = 0;
SETTINGS.plot_it            = 0;
SETTINGS.plot_colours       = [];
SETTINGS.singles            = 0;

%% INPUT SETTINGS
SETTINGS = setter(SETTINGS, varargin);

%% DERIVED SETTINGS
if isempty(SETTINGS.orientationdata)
    tr_n = size(TARGETDATA,1);
    SETTINGS.orientationdata = (1:tr_n)';
end
if isempty(SETTINGS.cat_ids)
    SETTINGS.cat_ids = unique(TARGETDATA);
end
if SETTINGS.plot_it == 1 & isempty(SETTINGS.plot_colours)
    SETTINGS.plot_colours = rand(numel(SETTINGS.cat_ids),3);
end
TARGETDATA = sortrows([SETTINGS.orientationdata,TARGETDATA],1);
TARGETDATA(:,1) = []; 

%% ******************* MAIN: BOUNDARY DETERMINATION ***********************
%% PREPARATION
% FIND MODES OF CATEGORICAL ANSWERS ---------------------------------------
[tmp1, tmp2, maxima] = mode(TARGETDATA,2);
Modes = NaN(numel(maxima),1); % Preallocating for speed.
for ind = 1:numel(maxima)
    if numel(maxima{ind})>1
        Modes(ind,1) = NaN; % If there are equally frequent values...    
    elseif numel(maxima{ind})==1 % If there is one maximum...
        Modes(ind,1) = maxima{ind};
    end
end

% SORT BY ORIENTATION DATA ------------------------------------------------
providata = sortrows([SETTINGS.orientationdata, Modes],1);
Modes = providata(:,2);
orientationdata_sorted = providata(:,1);
% INITIALISE VARIABLES ----------------------------------------------------
cat_tot = numel(SETTINGS.cat_ids); % Total number of categories
tr_n = size(Modes,1); % Total number of cases.
BOUNDARIES = NaN(cat_tot,2);
boundinds1 = NaN(cat_tot,2);
boundinds2 = NaN(cat_tot,2);
BOUND_IDS = NaN(cat_tot,2);

%% Find manifest/preliminry categories:
for cat = 1:cat_tot
    % CONVERT INTO BINARIES FOR EACH CATEGORY SEPARATELY ------------------
    target = SETTINGS.cat_ids(cat);
    inds = find(Modes == target);
    bin_responses(:,cat) = zeros(tr_n,1);
    bin_responses(inds,cat) = 1;
    % FIND INDICES OF BOUNDARIES ------------------------------------------
    % = 2 steps in order to find interval despite inhomogenities.
    % 1. step = find interval where all answers correspond to the category.
    boundinds1(cat,1:2) = intervaldetector(bin_responses(:,cat),1,0); % ACCURACY = 0 --> All must be correct.
    % 2. step = find interval where at least 0.5 of the answers correspond to the category
    if size(TARGETDATA,2) == 1
        ave_responses(:,cat) = neighbouraverager(bin_responses(:,cat), orientationdata_sorted, 4); % XXXX This should actually only be done for data, where there is only 1 response per x value; otherwise, frequencies should be taken!!!
        boundinds2(cat,1:2) = intervaldetector(ave_responses(:,cat),1,0.5);
    else
        ave_responses(:,cat) = nanmean((TARGETDATA==cat),2);
        ave_responses(isnan(Modes),:) = 0;
        if sum(bin_responses(:,cat))>1 % If the interval is larger than 1
            boundinds2(cat,1:2) = intervaldetector(ave_responses(:,cat),1,0.5); % ACCURACY = 0.5 --> 50% must be correct.
        elseif sum(bin_responses(:,cat))==1 % If it is a "SINGLE"
            if ave_responses(logical(bin_responses(:,cat)),cat) == 1 % If it is 100% right.
                boundinds2(cat,1:2) = find(bin_responses(:,cat));
            else
                boundinds2(cat,1:2) = NaN;
            end
        end
    end
    if ~isnan(boundinds2(cat,1:2)) % If interval boundaries monotonic (i.e. interval homogeneous).
        boundinds(cat,1:2) = boundinds2(cat,1:2);
    else
        boundinds(cat,1:2) = boundinds1(cat,1:2);
    end
    % DETERMINE BOUNDARIES CORRESPONDING TO THE INDICES -------------------
    if ~isnan(boundinds(cat,:))
        BOUNDARIES(cat,1:2) = [orientationdata_sorted(boundinds(cat,1)),orientationdata_sorted(boundinds(cat,2))];
    else
        BOUNDARIES(cat,1:2) = [NaN, NaN];
    end
    meanbound(cat,1) = cat;
    if SETTINGS.fullcycle & BOUNDARIES(cat,1) > BOUNDARIES(cat,2) % If the category crosses the cycle
            meanbound(cat,2) = circularmean2(BOUNDARIES(cat,:),SETTINGS.onecycle);
    else
        meanbound(cat,2) = nanmean(BOUNDARIES(cat,:));        
    end
end
sortedbounds = [meanbound, BOUNDARIES];
sortedbounds = sortrows(sortedbounds,2);
inds = ~isnan(sortedbounds(:,2));
sortedbounds = sortedbounds(inds,:);
inds = find(BOUNDARIES(:,1)>BOUNDARIES(:,2),1);
if ~isempty(inds) | SETTINGS.fullcycle == 1 % XXXX 2009apr30
    sortedbounds = [sortedbounds;sortedbounds(1,:)]; % Roundtour
end

%% CHECK FOR KEY MISTAKES IN THE TRANSITION BETWEEN CATEGORIES (Press the wrong key)
% If a value occurs within a category which corresponds neither to this
% interval-category nor to any adjacent interval-category it will be
% replaced (current implementation) by the value of the intervalcategory).
cat_n = size(sortedbounds,1);
if size(TARGETDATA,2) == 1
    for cat = 2:cat_n-1
        target    = sortedbounds(cat,1);
        successor = sortedbounds(cat+1,1);
        precursor = sortedbounds(cat-1,1);
        if isnan(boundinds1(precursor,2)) | isnan(boundinds1(successor,1)) % added 2010jun17
            continue
        end
        inds1 = find(Modes(boundinds1(precursor,2):boundinds1(target,1)) == successor);
        inds2 = find(Modes(boundinds1(target,2):boundinds1(successor,1)) == precursor);
        if numel(inds1) > 0
            bin_responses(boundinds1(precursor,2)+inds1-1,target) = 1;
        end
        if numel(inds2) > 0
            bin_responses(boundinds1(target,2)+inds2-1,target) = 1;
        end
        boundinds1(target,1:2) = intervaldetector(bin_responses(:,target),1,0);
        if size(TARGETDATA,2) == 1
            ave_responses(:,cat) = neighbouraverager(bin_responses(:,target), orientationdata_sorted, 4); % XXXX Again this should only be done for single responses per x value, for several responses, frequencies should be used directly.
            boundinds2(target,1:2) = intervaldetector(ave_responses(:,cat),1,0.5);
            if ~isnan(boundinds2(target,1:2))
                boundinds(target,1:2) = boundinds2(target,1:2);
            else
                boundinds(target,1:2) = boundinds1(target,1:2);
            end
        else
            boundinds(target,1:2) = boundinds1(target,1:2);
        end
        % RE-DETERMINE BOUNDARIES CORRESPONDING TO THE INDICES ----------------
        if ~isnan(boundinds(target,:))
            BOUNDARIES(target,1:2) = [orientationdata_sorted(boundinds(target,1)),orientationdata_sorted(boundinds(target,2))];
        else
            BOUNDARIES(target,1:2) = [NaN, NaN];
        end
    end
end
BOUND_IDS = boundinds; 

%% SINGLES (Categories with just one element)
if SETTINGS.singles
    inds_s = find(isnan(BOUND_IDS(:,1))); % Those for which there is no boundary yet.
    for nr = 1:numel(inds_s)
        ind = inds_s(nr);
        inds_out = find(Modes == ind); 
        if ~isempty(inds_out) % If there is any, that obtained a mode of responses (= candidate for a category)
            lower = find(inds_out(1) == mod(BOUND_IDS(:,1)-1, SETTINGS.onecycle)); % Modulo to apply to circular scales
            upper = find(inds_out(end) == mod(BOUND_IDS(:,2)+1, SETTINGS.onecycle));
%            if ~isempty(lower) & ~isempty(upper)
                BOUND_IDS(ind,:)  = [inds_out(1), inds_out(end)]; 
                BOUNDARIES(ind,:) = [orientationdata_sorted(inds_out(1)), orientationdata_sorted(inds_out(end))];
                ave(1,1) = ind;
                ave(1,2) = circularmean2(BOUNDARIES(ind,:),SETTINGS.onecycle);
                if SETTINGS.fullcycle
                    a = find(sortedbounds(1:end-1,2) < ave(2), 1, 'last');
                    b = find(sortedbounds(1:end-1,2) > ave(2), 1, 'first');
                else
                    a = find(sortedbounds(:,2) < ave(2), 1, 'last');
                    b = find(sortedbounds(:,2) > ave(2), 1, 'first');
                end
                sortedbounds = [sortedbounds(1:a,:); [ave, BOUNDARIES(ind,:)]; sortedbounds(b:end,:)];
%            end
        end
    end
end

%% INTERPOLATE ADJACENT BOUNDARIES (if gaps == 0)
if SETTINGS.gaps == 0
    cat_n = size(sortedbounds,1)-1;
    for cat = 1:cat_n
        id = sortedbounds(cat,1);   % INDEX OF PRECEDING CATEGORY
        up = sortedbounds(cat+1,1); % INDEX OF SUCCEDING CATEGORY
        
        % INTERPOLATE -----------------------------------------------------
        BOUNDARIES(id,2) = circularmean2([BOUNDARIES(id,2),BOUNDARIES(up,1)],SETTINGS.onecycle);    % Draw boundary in the center of the transitional interval.
        BOUNDARIES(up,1) = BOUNDARIES(id,2);                                                        % Set succeding identical to preceding.
        if ~isempty(find(orientationdata_sorted < (BOUNDARIES(id,2)),1,'last'))
            BOUND_IDS(id,2)  = find(orientationdata_sorted < (BOUNDARIES(id,2)),1,'last'); 
        else
            BOUND_IDS(id,2)  = max(orientationdata_sorted);
        end
        if ~isempty(find(orientationdata_sorted > BOUNDARIES(id,2),1,'first'))
            BOUND_IDS(up,1)  = find(orientationdata_sorted > BOUNDARIES(id,2),1,'first');
        else
            BOUND_IDS(up,1)  = min(orientationdata_sorted);
        end
    end
    if ~SETTINGS.fullcycle % If not periodical
        first   = sortedbounds(1,1);
        last    = sortedbounds(end,1);
        BOUNDARIES(first, 1) = NaN;
        BOUNDARIES(last, 2) = NaN;
    end
end

%% PLOT
if SETTINGS.plot_it == 1
    cat_n = size(BOUNDARIES,1);
    colormap(SETTINGS.plot_colours);
    rw_n = size(TARGETDATA,1);
    cl_n = size(TARGETDATA,2);
    hold on
    for cat = 1:cat_n
        inds = TARGETDATA == cat;
        Y = ones(rw_n,1)*(1:cl_n).*inds;
        Y(Y == 0) = NaN;
        X = SETTINGS.orientationdata*ones(1,cl_n);
        h = plot(X, Y, '.', 'MarkerEdgeColor', SETTINGS.plot_colours(cat,:));
        HANDLE.cat(cat,1)   = h(1);
        text(BOUNDARIES(cat,1), cl_n, num2str(BOUNDARIES(cat,1)), 'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Center')
    end
    HANDLE.bound_low    = line([BOUNDARIES(:,1), BOUNDARIES(:,1)]', (ones(cat_n,1)*[0 cl_n])', 'Color', [0.5 0.5 0.5]);
    HANDLE.bound_up     = line([BOUNDARIES(:,2), BOUNDARIES(:,2)]', (ones(cat_n,1)*[0 cl_n])', 'Color', [0 0 0]);
    hold off
    legend(HANDLE.cat, cellstr(num2str((1:cat_n)')), 'location', 'EastOutside');
    xlim([min(SETTINGS.orientationdata), max(SETTINGS.orientationdata)]);
    set(gca, 'Ytick', 1:cl_n);
    ylim([0 cl_n]);
end

%% **************************** SUBFUNCTIONS ******************************
% see original functions for details.

%% circularmean2
function AVERAGE = circularmean2(VALUES, ONECYCLE, amplitudes)
% VERSION: 2008aug18 [cw].
if nargin < 3
  amplitudes = ones(size(VALUES));
    if nargin < 2
        ONECYCLE = 2*pi;
    end
end
provivalues = (VALUES/ONECYCLE)*2*pi;
cm = atan2(sum(amplitudes.*sin(provivalues)), sum(amplitudes.*cos(provivalues)));
AVERAGE = cm * ONECYCLE/(2*pi);
AVERAGE = mod(AVERAGE, ONECYCLE);

%% intervaldetector
function [BOUNDARIES, allintervals] = intervaldetector(DATA, CRITERION, ACCURACY)
% Derived from intervaldetector.m; change: minimum number of elements = 1
% instead of 3. 
tot_n   = size(DATA, 1);
allintervals = [];
% Identify indices --------------------------------------------------------
indices = find(DATA >= CRITERION - ACCURACY &  DATA <= CRITERION + ACCURACY);

% Identify interval boundaries --------------------------------------------
ind_N=length(indices);
if (ind_N < 1) || (ind_N == tot_n)
    BOUNDARIES = [NaN NaN];
else
    CI(1) = indices(1);
    member_counter = 0;
    interval_counter = 0;
    for i=1:ind_N-1
        if indices(i+1) == indices(i)+1
            CI(2) = indices(i+1); % Redefine upper limit
            member_counter = member_counter+1;
        else
            if member_counter > 0
                interval_counter = interval_counter+1;
                allintervals(interval_counter,1:3) = [CI, member_counter+1];
            end
            CI(1) = indices(i+1); % Redefine lower limit
            member_counter = 0;
        end
    end
    if indices(1)==1 && indices(ind_N) == tot_n % If an interval at the end rebegins at the beginning of the DATA
        member_counter2 = 1;
        for j=ind_N:-1:2
            if indices(j-1) == indices(j)-1
                CI2(1) = indices(j-1); % Redefine lower limit
                member_counter2 = member_counter2+1; 
            else
                CI2(2) = indices(j-1); % Redefine upper limit
                member_counter2 = 1;
            end
        end
        member_counter = member_counter + member_counter2;
        CI(2)=CI2(2);
    end
    if member_counter > 0
        interval_counter = interval_counter+1;
        allintervals(interval_counter,1:3) = [CI, member_counter+1];
    end
    if isempty(allintervals)
        BOUNDARIES = [NaN];
    else
        ind = find(allintervals(:,3)==max(allintervals(:,3)));
        BOUNDARIES = allintervals(ind,1:2);
    end
    if size(BOUNDARIES,1) > 1 % If there are two intervals that have both the max number of members.
%        BOUNDARIES = [NaN]; % ALT
        BOUNDARIES = BOUNDARIES(1,:); % XXXX NEU --> ACHTUNG: SCHNELLE LÖSUNG: Scheint zu funktionieren XXXX
    end    
end

%% neighbouraverager
function [NEIGHBOUR_AVERAGES, ORIENTATION_DATA_NEW, IDs, UNNEIGHBOURED] = neighbouraverager(AVERAGER_DATA, ORIENTATION_DATA, NEIGHBOURS_N)
% VERSION: 2007nov [cw].
% Constants ---------------------------------------------------------------
tot_n = size(AVERAGER_DATA,1);
kernel_size = NEIGHBOURS_N+1;
conv_kernel= ones(kernel_size,1)*(1/kernel_size);
id = (1:tot_n)'; 
odds_indicator = mod(kernel_size, 2);
if odds_indicator == 0
    disp('CAUTION: YOU CHOSE AN ODD NUMBER OF NEIGHBOURS. Thus the data cannot be in the center of the neighbours.');
end
% Preparing data ----------------------------------------------------------
    data            = [id AVERAGER_DATA ORIENTATION_DATA];
    data_sorted     = sortrows(data, 3);
    data_sorted_2x  = [data_sorted; data_sorted];
% Convolusion -------------------------------------------------------------
    data_conv       = conv(data_sorted_2x(:,2), conv_kernel);
% Extract relevant data ---------------------------------------------------
    data_conv_pure  = [data_conv(1:tot_n, :), data_conv(tot_n+1:2*tot_n, :)]; % Eliminate byproducts of kernel 
% NEIGHBOUR_AVERAGES  = data_conv_pure(:,2);
    data_neighboured  = [data_conv_pure((kernel_size-1)/2+1:tot_n, :); data_conv_pure(1:(kernel_size-1)/2, :)]; % Displace kernel results to the right indices
NEIGHBOUR_AVERAGES =data_neighboured(:,2);
if nargout>=2
    ORIENTATION_DATA_NEW = data_sorted(:,3);
    if nargout>=3
        IDs = data_sorted(:,1);
        if nargout>=4
            UNNEIGHBOURED=data_neighboured(:,1);
        end
    end
end

%% setter
function [SETTINGS, Inputedsettings] = setter(SETTINGS, INPUTARGUMENTS_FROM_VARARGIN)
% VERSION:2009mar18 [cw].
Inputedsettings = fieldnames(SETTINGS);
tmp(1:numel(Inputedsettings),1) = {'default'};
Inputedsettings = [Inputedsettings, tmp];
unused_counter = 0;
unused_fields = '';
checker = 0;
while ~isempty(INPUTARGUMENTS_FROM_VARARGIN)
    if isstruct(INPUTARGUMENTS_FROM_VARARGIN{1}) % For Input by structure: Replace fields by new ones
        % Adapted from struc_updater (Version: % 2007 dec 6 by uk & cw;  for details see respective main function.
        NEWSETTINGS = INPUTARGUMENTS_FROM_VARARGIN{:};
        changedfields = fieldnames(NEWSETTINGS);
        for k = 1:numel(changedfields)
            SETTINGS.(changedfields{k}) = NEWSETTINGS.(changedfields{k});
            inds = find(strcmpi(changedfields{k}, Inputedsettings));
            if ~isempty(inds)
                Inputedsettings{inds,2} = 'input';
            else
                Inputedsettings{end+1,1} = 'useless';
                Inputedsettings{end,2} = changedfields{k};
                unused_counter = unused_counter+1;
                unused_fields = [unused_fields, ' ', changedfields{k}];
            end
        end
        INPUTARGUMENTS_FROM_VARARGIN(1) = [];
    else
        break;
    end
end
if unused_counter
    fprintf(['INPUT WARNING: There are ', num2str(unused_counter) ,' unused fields in your Input:\n']);
    fprintf([unused_fields, '\n']);
end
if ~isempty(INPUTARGUMENTS_FROM_VARARGIN)
    if ischar(INPUTARGUMENTS_FROM_VARARGIN{1}) % For Input by labeling: Convert labels into structure fields.
        input_keywords = Inputedsettings(:,1);
        inds = find(strcmpi(INPUTARGUMENTS_FROM_VARARGIN, 'availables')); % Option to display keywords;
            if ~isempty(inds) % DOES IT WORK (?) 2009mar18;
                for fld = 1:numel(input_keywords)
                    fprintf([input_keywords{fld}, '\n']);
                end
                INPUTARGUMENTS_FROM_VARARGIN(inds) = [];
                checker = 1;
            end
        if rem(size(INPUTARGUMENTS_FROM_VARARGIN),2)~=0
            error('INPUT ERROR: Input arguments are not even! Did you forget to specify the name/value of the input?');
        else
            for fld = 1:size(input_keywords)
                inds = find(strcmpi(input_keywords{fld},INPUTARGUMENTS_FROM_VARARGIN));
                if ~isempty(inds)
                    SETTINGS.(INPUTARGUMENTS_FROM_VARARGIN{inds})=INPUTARGUMENTS_FROM_VARARGIN{inds+1};
                    Inputedsettings(fld,2) = {'input'};
                    checker = 1;
                end
            end
            if ~checker
                error('INPUT ERROR: There are inadequate labels in your Input arguments!');
            end
        end
    else
        error('INPUT ERROR: Use sturcture for structure mode or keywords in '''' followed by values for parameter mode or both');
    end
end

%% example
function [TARGETDATA, ORIENTATIONDATA, CAT_IDS] = example
% Example data from colournaming by constant stimuli using the 8 chromatic
% Basic Colour Terms in DKL azimuth degree space (cn3_bcts) for subject
% nak.
clc; close all;
ORIENTATIONDATA = (0:3:360)';
CAT_IDS = 1:8;
TARGETDATA =[...
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     3     1     1     2     3     3
     3     3     3     3     3     3
     3     3     5     2     1     3
     3     3     3     3     3     3
     3     3     3     3     3     3
     3     3     3     3     3     3
     3     3     3     3     3     3
     3     3     3     3     3     3
     3     3     3     3     3     3
     3     3     3     4     4     4
     3     3     4     4     4     4
     4     4     4     4     4     4
     4     4     3     3     4     4
     4     4     4     4     4     4
     4     4     4     4     4     4
     4     4     4     4     4     4
     4     4     4     4     4     4
     4     4     4     4     4     4
     4     4     4     5     4     4
     4     4     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     5     5     5     5     5     5
     6     5     5     5     5     6
     6     6     6     5     5     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     3     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     6
     6     6     6     6     6     7
     6     7     7     7     6     7
     7     6     7     7     6     7
     7     7     7     7     7     7
     7     7     7     7     7     7
     7     7     7     7     7     7
     7     7     7     7     7     7
     7     7     7     7     5     7
     7     7     7     7     7     7
     7     7     7     7     7     7
     7     7     7     7     7     7
     7     7     7     7     7     7
     7     7     7     7     7     7
     7     7     7     5     7     7
     7     7     7     7     7     7
     7     1     7     7     7     7
     1     1     7     7     7     7
     7     7     7     7     7     1
     1     1     1     7     7     1
     1     1     1     1     7     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1];