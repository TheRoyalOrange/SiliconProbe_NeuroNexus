% plotRankedShankProfile.m
%
% Description: Multi-animal version of the channel/shank label plots (plotProbeLabels.m),
%   built from the summary table. For the chosen animals, the shanks of each animal are
%   re-centred on its top-ranked shank for a chosen shank label (ShankRank_<score> == 1):
%   that shank is column 3, its neighbouring shanks are columns 1-2 (left) and 4-5 (right).
%   Rows are Chan_Depth (position among the shank's kept channels, 1 = top), so row 1 is the
%   first non-Ch_Remove channel of each shank. One figure per chosen feature (measure
%   column), a 5 x 1 tiled layout (4 x 1 without the normalized plot):
%     tile 1: line plot, one line per animal, of the feature per column (mean of that
%       animal's channel means in the column); condition A solid, condition B dashed (same
%       colour per animal). With do_stats (and condB), a star in the animal's colour above
%       a column where A and B differ (see Statistics)
%     tile 2 (normalized plot, default on): the same lines, each divided by its own value at
%       the reference (normto): a chosen column (default 3 = top shank) or the line's peak
%       (highest column value). Each condition line is normalized to its own reference, so
%       every line is 1 at the reference; an animal whose line has no value at a chosen
%       reference column gets no normalized line
%     last 3 tiles: heatmap (5 columns x rows) of the mean over animals of each animal's
%       per-position mean over trials: condition A alone (green), or B - A (each animal's
%       difference first, then the mean; blue-white-red, white = 0). NaN cells are black.
%       Below each column: n = number of animals with data in that column (animals whose
%       top shank is near the probe edge do not have all 5 columns)
%   Optional trial filter: rows (trial x channel) are kept only where a table column meets a
%   condition (e.g. "P1_MUAAUC > 100"), before any averaging, for both conditions.
%   Statistics (do_stats, condB given): per animal and column, one value per trial = mean of
%   the feature over the column's channels in that trial (trial = Condition_FileNum + Trial);
%   Wilcoxon rank-sum test (ranksum) of the condA trials vs the condB trials (unpaired, NaN
%   ignored, at least min_trials per condition), Benjamini-Hochberg FDR correction over all
%   tests in the figure (use_fdr), alpha = 0.05. alpha, use_fdr and min_trials are set at
%   the top of the function. Returns figure handles and the plotted numbers; saves nothing.
%
% Inputs:
%   summaryTable (table, or char/string path to its csv) - tableres from
%     SummaryAnalysis_MultiMouse_allChannels.m, one row per good trial x channel. Columns used:
%       Animal_Name, Condition_Name, Region (text)
%       Chan_Shank (double) - shank (ProbeMaps column) of the channel
%       Chan_Depth (double) - position among that shank's kept channels, 1 = top; NaN rows dropped
%       ShankRank_<score> (double) - rank of the channel's shank by that ShankLabels score, 1 = highest
%       Condition_FileNum, Trial (double) - identify a trial (statistics only)
%       measure columns named <P1|P2|All>_... (double) - the features offered
%     INFERRED DATA CONTRACT: neighbouring Chan_Shank numbers are physically neighbouring
%     shanks (ProbeMaps column order = physical shank order)
%   animals (string array, optional) - Animal_Name values; omitted/empty = list dialog (multiple)
%   condA (char/string, optional) - main condition (Condition_Name); omitted/empty = list dialog
%   condB (char/string, optional) - comparison condition; heatmap becomes B - A, line plot adds
%     dashed B lines; omitted/empty = no comparison
%   rankfield (char/string, optional) - ShankRank_<score> column, with or without the
%     "ShankRank_" prefix (e.g. "SponActivity"); omitted/empty = list dialog
%   features (string array, optional) - measure columns to plot, one figure each;
%     omitted/empty = list dialog (multiple) of the measure columns with data
%   region (char/string, optional) - Region value; omitted/empty = the only region of the chosen
%     animals, or a list dialog if there is more than one
%   trialfilter (string array, optional) - conditions "<column> <op> <value>", op one of
%     > >= < <= == ~=; a numeric column takes a number, a text column only == / ~= with a text
%     value (quotes optional). Several are combined with AND. Rows where the column is NaN are
%     dropped. Omitted/empty = no filtering
%   do_stats (logical, optional) - true: per-animal, per-column A vs B test (only with condB);
%     default false
%   normto (double or char/string, optional) - normalized line plot: a column 1-5 to normalize
%     to (3 = top shank, 1/2 = -2/-1, 4/5 = +1/+2), "peak" (each line's highest column value),
%     or "none" / 0 / false (no normalized plot). Omitted/empty = 3
%
% Outputs:
%   figs (figure handle array, 1 x nFeatures) - one figure per feature
%   data (struct array, 1 x nFeatures) - the plotted numbers, fields:
%     .feature, .condA, .condB, .rankfield, .region (string), .trialfilter (string array)
%     .animals (string, 1 x nAnimals) - animals plotted (those with a rank-1 shank)
%     .topshank (double, 1 x nAnimals) - Chan_Shank of each animal's top-ranked shank
%     .mapA, .mapB (double, nAnimals x rows x 5) - per-animal mean over trials per position
%       (mapB empty without condB)
%     .heat (double, rows x 5) - heatmap values (mean of mapA, or of mapB - mapA)
%     .n (double, 1 x 5) - animals with data per column
%     .colA, .colB (double, nAnimals x 5) - line plot values
%     .normto (double or string) - reference of the normalized plot (column, "peak" or "none")
%     .colA_norm, .colB_norm (double, nAnimals x 5) - normalized line plot values (empty if none)
%     .p, .sig (nAnimals x 5) - rank-sum p-values (NaN = not tested) and significance after
%       FDR (empty without statistics)
%
% Dependencies: SummaryAnalysis_MultiMouse_allChannels.m for summaryTable, including its
%   separate section that adds Chan_Shank / Chan_Depth / ShankRank_<score> (ShankLabels from
%   ChannelShankLabelling_SponData.m / ChannelShank_activitylabelling.m). Style follows
%   plotProbeLabels.m and plotLWIntensityCompare.m. MATLAB R2019b+ (tiledlayout);
%   Statistics and Machine Learning Toolbox (ranksum) for do_stats.

function [figs, data] = plotRankedShankProfile(summaryTable, animals, condA, condB, rankfield, features, region, trialfilter, do_stats, normto)

figs = gobjects(1,0);
data = struct([]);
%statistics settings (do_stats): Wilcoxon rank-sum of A vs B trials per animal and column
alpha = 0.05;       %significance level
use_fdr = true;     %Benjamini-Hochberg FDR correction across all tests in a figure
min_trials = 3;     %minimum trials (non-NaN) per condition to test
offsets = -2:2;     %columns: shanks relative to the top-ranked shank

%% summary table: from a csv path or given directly; check the columns needed
if ischar(summaryTable) || isstring(summaryTable)
    summaryTable = readtable(char(summaryTable), 'TextType','string');
end
assert(istable(summaryTable), 'summaryTable must be a table (tableres) or the path to its csv')
vn = string(summaryTable.Properties.VariableNames);
needcols = ["Animal_Name","Condition_Name","Region","Chan_Shank","Chan_Depth"];
missingcols = needcols(~ismember(needcols, vn));
assert(isempty(missingcols), 'summaryTable is missing column(s): %s', strjoin(missingcols,', '))
T = summaryTable;
T.Animal_Name = string(T.Animal_Name);
T.Condition_Name = string(T.Condition_Name);
T.Region = string(T.Region);

%% trial filter (before anything else, so all choices below see only the kept rows)
if nargin < 8 || isempty(trialfilter)
    trialfilter = strings(1,0);
end
trialfilter = string(trialfilter);
nbefore = height(T);
for k = 1:numel(trialfilter)
    T = T(filterRows(T, trialfilter(k)), :);
end
nafter = height(T);
assert(nafter > 0, 'No rows left after the trial filter (%s)', strjoin(trialfilter,' & '))

%% animals
allanimals = unique(T.Animal_Name,'stable');
if nargin < 2 || isempty(animals)
    [sel,ok] = listdlg('ListString',cellstr(allanimals),'SelectionMode','multiple', ...
        'PromptString','Animals to include:','ListSize',[250 200]);
    assert(ok && ~isempty(sel), 'No animal selected')
    animals = allanimals(sel);
end
animals = string(animals(:)');
badanimal = setdiff(animals, allanimals, 'stable');
assert(isempty(badanimal), 'summaryTable has no rows for animal(s): %s', strjoin(badanimal,', '))
T = T(ismember(T.Animal_Name, animals), :);

%% region
regions = unique(T.Region,'stable');
if nargin < 7 || isempty(region)
    if numel(regions) == 1
        region = regions;
    else
        [sel,ok] = listdlg('ListString',cellstr(regions),'SelectionMode','single', ...
            'PromptString','Which region?','ListSize',[250 120]);
        assert(ok && ~isempty(sel), 'No region selected')
        region = regions(sel);
    end
end
region = string(region);
assert(ismember(region, regions), 'The chosen animals have no rows for region %s', region)
T = T(T.Region == region, :);

%% conditions
conds = unique(T.Condition_Name,'stable');
if nargin < 3 || isempty(condA)
    [sel,ok] = listdlg('ListString',cellstr(conds),'SelectionMode','single', ...
        'PromptString','Condition A (main condition):','ListSize',[250 200]);
    assert(ok && ~isempty(sel), 'No condition A selected')
    condA = conds(sel);
end
condA = string(condA);
if nargin < 4 || isempty(condB)
    condB = string(missing);
else
    condB = string(condB);
end
usedconds = condA;
if ~ismissing(condB), usedconds = [condA condB]; end
badcond = setdiff(usedconds, conds, 'stable');
assert(isempty(badcond), 'summaryTable has no rows for condition(s) %s (%s)', strjoin(badcond,', '), region)
T = T(ismember(T.Condition_Name, usedconds), :);

%% shank-rank label
rankcols = vn(startsWith(vn,"ShankRank_"));
assert(~isempty(rankcols), 'summaryTable has no ShankRank_<score> columns')
if nargin < 5 || isempty(rankfield)
    [sel,ok] = listdlg('ListString',cellstr(rankcols),'SelectionMode','single', ...
        'PromptString','Centre on the top shank of which label?','ListSize',[250 150]);
    assert(ok && ~isempty(sel), 'No shank label selected')
    rankfield = rankcols(sel);
end
rankfield = string(rankfield);
if ~startsWith(rankfield,"ShankRank_"), rankfield = "ShankRank_" + rankfield; end
assert(ismember(rankfield, rankcols), 'summaryTable has no column %s (available: %s)', rankfield, strjoin(rankcols,', '))

%% features: measure columns with data
measurecols = vn(~cellfun(@isempty, regexp(cellstr(vn), '^(P1|P2|All)_', 'once')));
hasdata = arrayfun(@(c) isnumeric(T.(c)) && any(~isnan(T.(c))), measurecols);
measurecols = measurecols(hasdata);
assert(~isempty(measurecols), 'summaryTable has no measure columns with data for these animals/conditions')
if nargin < 6 || isempty(features)
    [sel,ok] = listdlg('ListString',cellstr(measurecols),'SelectionMode','multiple', ...
        'PromptString','Data features to plot (one figure each):','ListSize',[300 300]);
    assert(ok && ~isempty(sel), 'No data feature selected')
    features = measurecols(sel);
end
features = string(features(:)');
badfeat = setdiff(features, measurecols, 'stable');
assert(isempty(badfeat), 'Feature(s) not in summaryTable or without data: %s', strjoin(badfeat,', '))

if nargin < 9 || isempty(do_stats)
    do_stats = false;
end
do_stats = logical(do_stats) && ~ismissing(condB);

%normalized line plot: reference column 1-5, "peak", or none
if nargin < 10 || isempty(normto)
    normto = 3; %default: normalize to the top shank (column 3)
end
if islogical(normto)
    if normto, normto = 3; else, normto = "none"; end
elseif isnumeric(normto) && normto == 0
    normto = "none";
end
if isnumeric(normto)
    assert(isscalar(normto) && ismember(normto,1:5), 'normto must be a column 1-5, "peak" or "none"')
    normtxt = sprintf('normalized to column %s', offsetLabel(offsets(normto)));
else
    normto = lower(string(normto));
    assert(ismember(normto,["peak","none"]), 'normto must be a column 1-5, "peak" or "none"')
    normtxt = 'normalized to each line''s peak column';
end
shownorm = ~(isstring(normto) && normto == "none");

%% each animal's top-ranked shank -> column of each row (3 = top shank)
%(animals without rows for condA, or condB if given, are left out)
topshank = NaN(1,numel(animals));
for a = 1:numel(animals)
    missingc = usedconds(~arrayfun(@(c) any(T.Animal_Name == animals(a) & T.Condition_Name == c), usedconds));
    if ~isempty(missingc)
        warning('plotRankedShankProfile: %s has no rows for condition(s) %s (%s) - left out', animals(a), strjoin(missingc,', '), region)
        continue
    end
    s1 = unique(T.Chan_Shank(T.Animal_Name == animals(a) & T.(rankfield) == 1));
    s1 = s1(~isnan(s1));
    if isempty(s1)
        warning('plotRankedShankProfile: %s has no shank with %s == 1 (%s) - left out', animals(a), rankfield, region)
    else
        topshank(a) = min(s1);
    end
end
keepanimal = ~isnan(topshank);
animals = animals(keepanimal);
topshank = topshank(keepanimal);
assert(~isempty(animals), 'None of the chosen animals has rows for %s and a shank with %s == 1', strjoin(usedconds,' and '), rankfield)
[~,aidx] = ismember(T.Animal_Name, animals);
T = T(aidx > 0, :);
aidx = aidx(aidx > 0);
T.col = T.Chan_Shank - topshank(aidx)' + 3;
T = T(T.col >= 1 & T.col <= 5 & ~isnan(T.Chan_Depth), :);
nrow = max(T.Chan_Depth);
na = numel(animals);

%colour maps as in plotLWIntensityCompare.m; one line colour per animal
greenmap = interp1([0 0.5 1],[0.95 0.98 0.93; 0.55 0.82 0.40; 0.18 0.60 0.22],linspace(0,1,256));
divmap = interp1([0 0.5 1],[0.13 0.40 0.75; 1 1 1; 0.80 0.15 0.15],linspace(0,1,256));
acol = lines(na);
filtertxt = '';
if ~isempty(trialfilter)
    filtertxt = sprintf('trials: %s (%d of %d table rows kept)', strjoin(trialfilter,' & '), nafter, nbefore);
end

%% one figure per feature
for fi = 1:numel(features)
    feat = features(fi);

    %per-animal mean over trials at each position (rows x 5)
    mapA = NaN(na,nrow,5);  mapB = [];
    for a = 1:na
        mapA(a,:,:) = positionMap(T, feat, animals(a), condA, nrow);
    end
    if ~ismissing(condB)
        mapB = NaN(na,nrow,5);
        for a = 1:na
            mapB(a,:,:) = positionMap(T, feat, animals(a), condB, nrow);
        end
        X = mapB - mapA;
    else
        X = mapA;
    end
    heat = reshape(mean(X,1,'omitnan'), nrow, 5);
    n = reshape(sum(any(~isnan(X),2),1), 1, 5);
    colA = reshape(mean(mapA,2,'omitnan'), na, 5);
    colB = [];
    if ~isempty(mapB), colB = reshape(mean(mapB,2,'omitnan'), na, 5); end
    %normalized lines: each line divided by its own value at the reference column / its peak
    colA_norm = [];  colB_norm = [];
    if shownorm
        colA_norm = normLines(colA, normto);
        if ~isempty(colB), colB_norm = normLines(colB, normto); end
    end

    %statistics: per animal and column, rank-sum of the A vs B trial values (trial = mean over the column's channels)
    P = [];  sig = [];
    if do_stats
        P = NaN(na,5);
        for a = 1:na
            for c = 1:5
                xa = trialValues(T, feat, animals(a), condA, c);
                xb = trialValues(T, feat, animals(a), condB, c);
                if numel(xa) >= min_trials && numel(xb) >= min_trials
                    P(a,c) = ranksum(xa, xb);
                end
            end
        end
        sig = significant(P, alpha, use_fdr);
    end

    %figure: line plot on top, heatmap below (as plotProbeLabels.m)
    if ismissing(condB)
        what = condA;
    else
        what = condB + " - " + condA;
    end
    ttl = {sprintf('%s: %s', feat, what), ...
        sprintf('centred on top shank by %s (%s), %d animals', erase(rankfield,"ShankRank_"), region, na)};
    if ~isempty(filtertxt)
        ttl{end+1} = filtertxt; %#ok<AGROW>
    end
    if do_stats
        ttl{end+1} = sprintf('* = rank-sum %s vs %s per animal, %s alpha = %.2f', condA, condB, ternary(use_fdr,'BH-FDR,',''), alpha); %#ok<AGROW>
    end
    figs(end+1) = figure('Name',sprintf('%s: %s (%s)',feat,what,erase(rankfield,"ShankRank_")),'Color','w', ...
        'Units','pixels','Position',[100 60 900 850]); %#ok<AGROW>
    tl = tiledlayout(figs(end),4+shownorm,1,'TileSpacing','compact');
    title(tl,ttl,'Interpreter','none');

    %tile 1: one line per animal (A solid, B dashed)
    ax1 = nexttile(tl,1);
    hold(ax1,'on');
    hl = gobjects(1,na);
    for a = 1:na
        hl(a) = plot(ax1,1:5,colA(a,:),'-o','LineWidth',1.5,'Color',acol(a,:),'MarkerFaceColor',acol(a,:),'DisplayName',animals(a));
        if ~isempty(colB)
            plot(ax1,1:5,colB(a,:),'--o','LineWidth',1.5,'Color',acol(a,:),'HandleVisibility','off');
        end
    end
    if ~isempty(colB)
        plot(ax1,NaN,NaN,'k-','DisplayName',"solid = " + condA);
        plot(ax1,NaN,NaN,'k--','DisplayName',"dashed = " + condB);
    end
    xlim(ax1,[0.5 5.5]);
    xticks(ax1,1:5); xticklabels(ax1,{});
    ylabel(ax1,{char(feat),'(column mean)'},'Interpreter','none');
    grid(ax1,'on'); ax1.GridAlpha = 0.15; box(ax1,'off');
    lg = legend(ax1,'Interpreter','none');
    lg.Layout.Tile = 'east'; %outside both plots (as the colorbar), so line points stay above their heatmap columns
    %stars above the higher of an animal's A/B points (slightly staggered so animals do not overlap)
    if do_stats && any(sig,'all')
        yl = ylim(ax1);
        dy = 0.06*diff(yl);
        for a = 1:na
            for c = find(sig(a,:))
                y = max([colA(a,c) colB(a,c)]) + dy;
                text(ax1, c + (a-(na+1)/2)*0.08, y, '*', 'Color',acol(a,:), 'FontSize',16, 'FontWeight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
            end
        end
        ylim(ax1,[yl(1) max(yl(2), max([colA(sig) ; colB(sig)]) + 2*dy)]);
    end
    linked = ax1;

    %tile 2: the same lines normalized to the reference (each line 1 at its reference)
    if shownorm
        axN = nexttile(tl,2);
        hold(axN,'on');
        for a = 1:na
            plot(axN,1:5,colA_norm(a,:),'-o','LineWidth',1.5,'Color',acol(a,:),'MarkerFaceColor',acol(a,:));
            if ~isempty(colB_norm)
                plot(axN,1:5,colB_norm(a,:),'--o','LineWidth',1.5,'Color',acol(a,:));
            end
        end
        yline(axN,1,':','Color',[0.4 0.4 0.4]);
        xlim(axN,[0.5 5.5]);
        xticks(axN,1:5); xticklabels(axN,{});
        ylabel(axN,{'normalized',sprintf('(%s)',erase(normtxt,'normalized '))},'Interpreter','none');
        grid(axN,'on'); axN.GridAlpha = 0.15; box(axN,'off');
        linked = [linked axN]; %#ok<AGROW>
    end

    %last 3 tiles: heatmap, n per column below
    ax2 = nexttile(tl,2+shownorm,[3 1]);
    imagesc(ax2,1:5,1:nrow,heat,'AlphaData',~isnan(heat));
    ax2.Color = 'k';
    if ismissing(condB)
        colormap(ax2,greenmap);
        lims = [min(heat,[],'all','omitnan') max(heat,[],'all','omitnan')];
        if any(isnan(lims)), lims = [0 1]; end
        if lims(1) == lims(2), lims = lims(1) + [-1 1]; end
    else
        colormap(ax2,divmap);
        dmax = max(abs(heat),[],'all','omitnan');
        if isempty(dmax) || isnan(dmax) || dmax == 0, dmax = 1; end
        lims = [-dmax dmax];
    end
    clim(ax2,lims);
    cb = colorbar(ax2);
    cb.Layout.Tile = 'east';
    cb.Label.String = sprintf('%s: %s (mean of animal means)', feat, what);
    cb.Label.Interpreter = 'none';
    xticks(ax2,1:5); yticks(ax2,1:nrow);
    %two-line tick labels: offset from the top shank, then n = animals with data in that column
    ticklab = cell(1,5);
    for c = 1:5
        ticklab{c} = sprintf('%s\\newline n = %d', offsetLabel(offsets(c)), n(c));
    end
    ax2.TickLabelInterpreter = 'tex';
    xticklabels(ax2,ticklab);
    xlabel(ax2,'Shank relative to the top-ranked shank');
    ylabel(ax2,'Chan\_Depth (1 = top kept channel)');
    linkaxes([linked ax2],'x');

    data(fi).feature = feat;
    data(fi).condA = condA;
    data(fi).condB = condB;
    data(fi).rankfield = rankfield;
    data(fi).region = region;
    data(fi).trialfilter = trialfilter;
    data(fi).animals = animals;
    data(fi).topshank = topshank;
    data(fi).mapA = mapA;
    data(fi).mapB = mapB;
    data(fi).heat = heat;
    data(fi).n = n;
    data(fi).colA = colA;
    data(fi).colB = colB;
    data(fi).normto = normto;
    data(fi).colA_norm = colA_norm;
    data(fi).colB_norm = colB_norm;
    data(fi).p = P;
    data(fi).sig = sig;
end
end

function keep = filterRows(T, f)
%rows of T meeting one condition "<column> <op> <value>" (NaN in the column = not kept)
tok = regexp(char(f), '^\s*(\w+)\s*(>=|<=|==|~=|>|<)\s*(.+?)\s*$', 'tokens', 'once');
assert(~isempty(tok), 'Trial filter "%s" is not of the form "<column> <op> <value>" (op: > >= < <= == ~=)', f)
[colname, op, valtxt] = tok{:};
assert(ismember(colname, T.Properties.VariableNames), 'Trial filter "%s": summaryTable has no column %s', f, colname)
x = T.(colname);
if isnumeric(x) || islogical(x)
    val = str2double(valtxt);
    assert(~isnan(val), 'Trial filter "%s": %s is numeric, so the value must be a number', f, colname)
    x = double(x);
    switch op
        case '>',  keep = x > val;
        case '>=', keep = x >= val;
        case '<',  keep = x < val;
        case '<=', keep = x <= val;
        case '==', keep = x == val;
        case '~=', keep = x ~= val & ~isnan(x);
    end
else
    assert(ismember(op, {'==','~='}), 'Trial filter "%s": %s is text, so only == or ~= can be used', f, colname)
    val = string(regexprep(valtxt, '^["'']|["'']$', ''));
    if strcmp(op,'==')
        keep = string(x) == val;
    else
        keep = string(x) ~= val & ~ismissing(string(x));
    end
end
end

function M = positionMap(T, feat, animal, cond, nrow)
%mean of feature feat over the trials of one animal and condition at each position (rows x 5, NaN = no data)
sub = T.Animal_Name == animal & T.Condition_Name == cond;
M = accumarray([T.Chan_Depth(sub) T.col(sub)], T.(feat)(sub), [nrow 5], @(x) mean(x,'omitnan'), NaN);
end

function x = trialValues(T, feat, animal, cond, c)
%one value per trial: mean of feature feat over the channels of column c (trial = Condition_FileNum + Trial; NaN dropped)
sub = T.Animal_Name == animal & T.Condition_Name == cond & T.col == c;
if ~any(sub)
    x = [];
    return
end
g = findgroups(T.Condition_FileNum(sub), T.Trial(sub));
x = splitapply(@(v) mean(v,'omitnan'), T.(feat)(sub), g);
x = x(~isnan(x));
end

function sig = significant(P, alpha, use_fdr)
%significant tests (logical, same size as P); Benjamini-Hochberg FDR over the tested entries if use_fdr
sig = false(size(P));
tested = find(~isnan(P));
ntested = numel(tested);
if ntested == 0, return, end
p = P(tested);
if use_fdr
    [ps,ord] = sort(p);
    q = ps .* ntested ./ (1:ntested)';
    q = flipud(cummin(flipud(q)));   %monotone adjusted p-values
    padj = zeros(size(p));  padj(ord) = min(q,1);
    sig(tested) = padj < alpha;
else
    sig(tested) = p < alpha;
end
end

function N = normLines(C, normto)
%each row of C (animals x 5) divided by its value at column normto, or by its peak (max) if normto is "peak"
if isnumeric(normto)
    ref = C(:,normto);
else
    ref = max(C,[],2,'omitnan');
end
N = C ./ ref;
end

function s = offsetLabel(off)
%column label: 'top' for the top-ranked shank, else the signed offset ('-2', '+1', ...)
if off == 0
    s = 'top';
else
    s = sprintf('%+d', off);
end
end

function out = ternary(cond, a, b)
%a if cond, else b
if cond, out = a; else, out = b; end
end
