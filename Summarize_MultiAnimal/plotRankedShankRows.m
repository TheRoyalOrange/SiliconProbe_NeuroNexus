% plotRankedShankRows.m
%
% Description: Companion of plotRankedShankProfile.m (same arguments, same re-centring):
%   for the chosen animals, the shanks of each animal are re-centred on its top-ranked shank
%   for a chosen shank label (ShankRank_<score> == 1): that shank is column 3, its
%   neighbouring shanks are columns 1-2 (left) and 4-5 (right). Rows are Chan_Depth
%   (position among the shank's kept channels, 1 = top). One figure per chosen feature
%   (measure column), a 4 x 1 tiled layout:
%     tile 1: line plot, one line per animal, of the chosen label's shank scores from
%       ProbeInfo.ShankLabels (the values of the shank activity plots of
%       ChannelShankLabelling_SponData.m / plotProbeLabels.m) for the 5 aligned shanks;
%       with scorenorm, each animal's scores divided by its top shank's score (top = 1)
%     tiles 2-4: heatmap (5 columns x rows) of the mean over animals of each animal's
%       per-position mean over trials: condition A alone (green), or B - A (each animal's
%       difference first, then the mean; blue-white-red, white = 0). NaN cells are black.
%       Below each column: n = number of animals with data in that column.
%       Overlaid on each heatmap row (normto not "none"): a line across the 5 columns of the
%       row's mean over animals for condition A (solid), and with condB also for B (dashed;
%       B itself, not B - A), each normalized to its own reference (normto: a column, default
%       3 = top shank, or the line's peak). Each animal's own value per cell is drawn as a
%       point (o = condA, x = condB, slightly left/right of the column centre when both are
%       shown), divided by the same reference as its row's line, so the line is the mean of
%       its points. All rows share one vertical scale: the reference value 1 sits at the
%       row's centre (thin dotted line), higher values go up, and the largest deviation of
%       any line or point fills 45% of a row height; the title gives the scale
%   Optional trial filter: rows (trial x channel) are kept only where a table column meets a
%   condition (e.g. "P1_MUAAUC > 100"), before any averaging, for both conditions.
%   Statistics (do_stats, condB given): per heatmap cell (row x column), the feature values
%   of all chosen animals' trials at that position (one value per trial, animals pooled):
%   Wilcoxon rank-sum test (ranksum) of condA vs condB (unpaired, NaN ignored, at least
%   min_trials per condition), Benjamini-Hochberg FDR correction over all tested cells
%   (use_fdr), alpha = 0.05; a star in the corner of each significant cell. alpha, use_fdr,
%   min_trials and probeinfo_dir are set at the top of the function. Returns figure handles
%   and the plotted numbers; saves nothing.
%
% Inputs:
%   summaryTable (table, or char/string path to its csv) - tableres from
%     SummaryAnalysis_MultiMouse_allChannels.m, one row per good trial x channel. Columns used:
%       Animal_Name, Condition_Name, Region (text)
%       Chan_Shank (double) - shank (ProbeMaps column) of the channel
%       Chan_Depth (double) - position among that shank's kept channels, 1 = top; NaN rows dropped
%       ShankRank_<score> (double) - rank of the channel's shank by that ShankLabels score, 1 = highest
%       measure columns named <P1|P2|All>_... (double) - the features offered
%     INFERRED DATA CONTRACT: neighbouring Chan_Shank numbers are physically neighbouring
%     shanks (ProbeMaps column order = physical shank order)
%   animals (string array, optional) - Animal_Name values; omitted/empty = list dialog (multiple)
%   condA (char/string, optional) - main condition (Condition_Name); omitted/empty = list dialog
%   condB (char/string, optional) - comparison condition; heatmap becomes B - A, row lines add
%     dashed B lines; omitted/empty = no comparison
%   rankfield (char/string, optional) - ShankRank_<score> column, with or without the
%     "ShankRank_" prefix (e.g. "SponActivity"); omitted/empty = list dialog. <score> must also
%     be a field of the animals' ProbeInfo.ShankLabels{probe} for the top line plot
%   features (string array, optional) - measure columns to plot, one figure each;
%     omitted/empty = list dialog (multiple) of the measure columns with data
%   region (char/string, optional) - Region value; omitted/empty = the only region of the chosen
%     animals, or a list dialog if there is more than one. Also picks the probe in ProbeInfo
%     (matched in ProbeInfo.Areas, or used as the probe number)
%   trialfilter (string array, optional) - conditions "<column> <op> <value>", op one of
%     > >= < <= == ~=; a numeric column takes a number, a text column only == / ~= with a text
%     value (quotes optional). Several are combined with AND. Rows where the column is NaN are
%     dropped. Omitted/empty = no filtering
%   do_stats (logical, optional) - true: per-cell A vs B test (only with condB); default false
%   normto (double or char/string, optional) - row lines: a column 1-5 to normalize to
%     (3 = top shank, 1/2 = -2/-1, 4/5 = +1/+2), "peak" (each line's highest column value),
%     or "none" / 0 / false (no row lines). Omitted/empty = 3
%   scorenorm (logical, optional) - true: top line plot shows each animal's shank scores
%     divided by its top shank's score; false (default): raw scores
%   <probeinfo_dir>\<animal>-ProbeInfo.mat - each animal's ProbeInfo, fields .Areas and
%     .ShankLabels{probe}.<score> (double, 1 x shanks, ProbeMaps column order). An animal
%     without it gets no top line (warning); its heatmap data are still used
%
% Outputs:
%   figs (figure handle array, 1 x nFeatures) - one figure per feature
%   data (struct array, 1 x nFeatures) - the plotted numbers, fields:
%     .feature, .condA, .condB, .rankfield, .region (string), .trialfilter (string array)
%     .animals (string, 1 x nAnimals) - animals plotted (those with a rank-1 shank)
%     .topshank (double, 1 x nAnimals) - Chan_Shank of each animal's top-ranked shank
%     .scores (double, nAnimals x 5) - top line plot values (normalized if scorenorm)
%     .scores_raw (double, nAnimals x 5) - raw shank scores of the aligned shanks
%     .scorenorm (logical) - whether the top line plot is normalized
%     .mapA, .mapB (double, nAnimals x rows x 5) - per-animal mean over trials per position
%       (mapB empty without condB)
%     .heat (double, rows x 5) - heatmap values (mean of mapA, or of mapB - mapA)
%     .n (double, 1 x 5) - animals with data per column
%     .rowA, .rowB (double, rows x 5) - mean over animals of mapA / mapB (row lines, before normalizing)
%     .normto (double or string) - reference of the row lines (column, "peak" or "none")
%     .rowA_norm, .rowB_norm (double, rows x 5) - normalized row lines (empty if none)
%     .ptsA_norm, .ptsB_norm (double, nAnimals x rows x 5) - each animal's values divided by
%       its row line's reference (the points; empty if none)
%     .p, .sig (rows x 5) - rank-sum p-values per cell (NaN = not tested) and significance
%       after FDR (empty without statistics)
%
% Dependencies: SummaryAnalysis_MultiMouse_allChannels.m for summaryTable, including its
%   separate section that adds Chan_Shank / Chan_Depth / ShankRank_<score>; ProbeInfo.ShankLabels
%   from ChannelShankLabelling_SponData.m / ChannelShank_activitylabelling.m. Same arguments
%   and re-centring as plotRankedShankProfile.m. MATLAB R2019b+ (tiledlayout); Statistics and
%   Machine Learning Toolbox (ranksum) for do_stats.

function [figs, data] = plotRankedShankRows(summaryTable, animals, condA, condB, rankfield, features, region, trialfilter, do_stats, normto, scorenorm)

figs = gobjects(1,0);
data = struct([]);
%statistics settings (do_stats): Wilcoxon rank-sum of A vs B trials per heatmap cell
alpha = 0.05;       %significance level
use_fdr = true;     %Benjamini-Hochberg FDR correction across all tested cells of a figure
min_trials = 3;     %minimum trials (non-NaN) per condition to test
offsets = -2:2;     %columns: shanks relative to the top-ranked shank
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo'; %<animal>-ProbeInfo.mat (shank scores)
linefill = 0.45;    %largest row-line deviation, as a fraction of a row height

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

%top line plot: raw shank scores, or each animal's scores relative to its top shank
if nargin < 11 || isempty(scorenorm)
    scorenorm = false; %default: raw scores, as on the shank activity plots
end
scorenorm = logical(scorenorm);

%% each animal's top-ranked shank -> column of each row (3 = top shank)
%(animals without rows for condA, or condB if given, are left out)
topshank = NaN(1,numel(animals));
for a = 1:numel(animals)
    missingc = usedconds(~arrayfun(@(c) any(T.Animal_Name == animals(a) & T.Condition_Name == c), usedconds));
    if ~isempty(missingc)
        warning('plotRankedShankRows: %s has no rows for condition(s) %s (%s) - left out', animals(a), strjoin(missingc,', '), region)
        continue
    end
    s1 = unique(T.Chan_Shank(T.Animal_Name == animals(a) & T.(rankfield) == 1));
    s1 = s1(~isnan(s1));
    if isempty(s1)
        warning('plotRankedShankRows: %s has no shank with %s == 1 (%s) - left out', animals(a), rankfield, region)
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

%% shank scores of the chosen label from each animal's ProbeInfo (top line plot), aligned to the columns
scorename = erase(rankfield,"ShankRank_");
scores = NaN(na,5);
for a = 1:na
    pf = fullfile(probeinfo_dir, animals(a) + "-ProbeInfo.mat");
    if ~isfile(pf)
        warning('plotRankedShankRows: %s not found - no shank score line for %s', pf, animals(a))
        continue
    end
    pinfo = load(pf,'ProbeInfo');
    PI = pinfo.ProbeInfo;
    %the Region value can be a region name (matched in ProbeInfo.Areas) or a probe number
    prb = find(strcmp(string(PI.Areas), region));
    if isempty(prb)
        pnum = str2double(region);
        if ~isnan(pnum) && pnum == round(pnum) && pnum >= 1 && pnum <= numel(PI.ProbeMaps)
            prb = pnum;
        end
    end
    if numel(prb) ~= 1 || ~isfield(PI,'ShankLabels') || numel(PI.ShankLabels) < prb ...
            || ~isstruct(PI.ShankLabels{prb}) || ~isfield(PI.ShankLabels{prb}, scorename)
        warning('plotRankedShankRows: %s has no ShankLabels score %s for region %s - no shank score line', animals(a), scorename, region)
        continue
    end
    sc = double(PI.ShankLabels{prb}.(scorename)(:)');
    for c = 1:5
        s = topshank(a) + offsets(c);
        if s >= 1 && s <= numel(sc)
            scores(a,c) = sc(s);
        end
    end
end
scores_raw = scores;
if scorenorm
    scores = scores ./ scores(:,3); %each animal's scores relative to its top shank (= 1)
    scorelab = {char(scorename),'shank score','(/ top shank)'};
else
    scorelab = {char(scorename),'shank score'};
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

    %row lines: mean over animals per row, each normalized to its own reference; the animals'
    %own values (points) are divided by the same reference, so each line is the mean of its points
    rowA = reshape(mean(mapA,1,'omitnan'), nrow, 5);
    rowB = [];
    if ~isempty(mapB), rowB = reshape(mean(mapB,1,'omitnan'), nrow, 5); end
    rowA_norm = [];  rowB_norm = [];  ptsA_norm = [];  ptsB_norm = [];
    if shownorm
        refA = rowRef(rowA, normto);
        rowA_norm = rowA ./ refA;
        ptsA_norm = mapA ./ reshape(refA,1,nrow);
        if ~isempty(rowB)
            refB = rowRef(rowB, normto);
            rowB_norm = rowB ./ refB;
            ptsB_norm = mapB ./ reshape(refB,1,nrow);
        end
    end

    %statistics: per cell, rank-sum of the A vs B trial values (all animals pooled)
    P = [];  sig = [];
    if do_stats
        P = NaN(nrow,5);
        for r = 1:nrow
            for c = 1:5
                atpos = T.Chan_Depth == r & T.col == c;
                xa = T.(feat)(atpos & T.Condition_Name == condA);  xa = xa(~isnan(xa));
                xb = T.(feat)(atpos & T.Condition_Name == condB);  xb = xb(~isnan(xb));
                if numel(xa) >= min_trials && numel(xb) >= min_trials
                    P(r,c) = ranksum(xa, xb);
                end
            end
        end
        sig = significant(P, alpha, use_fdr);
    end

    %row-line scale shared by all rows: the largest deviation from 1 (lines or points) fills linefill of a row
    k = 0;
    if shownorm
        maxdev = max(abs([rowA_norm(:); rowB_norm(:); ptsA_norm(:); ptsB_norm(:)] - 1),[],'omitnan');
        if ~isempty(maxdev) && ~isnan(maxdev) && maxdev > 0
            k = linefill / maxdev;
        end
    end

    %figure: shank score line plot on top, heatmap with row lines below
    if ismissing(condB)
        what = condA;
    else
        what = condB + " - " + condA;
    end
    ttl = {sprintf('%s: %s', feat, what), ...
        sprintf('centred on top shank by %s (%s), %d animals', scorename, region, na)};
    if ~isempty(filtertxt)
        ttl{end+1} = filtertxt; %#ok<AGROW>
    end
    if shownorm
        if ismissing(condB)
            styletxt = sprintf('%s: line solid, animals o', condA);
        else
            styletxt = sprintf('%s: solid, animals o; %s: dashed, animals x', condA, condB);
        end
        if k > 0
            scaletxt = sprintf('; half a row = %.0f%% change', 100*0.5/k);
        else
            scaletxt = '';
        end
        ttl{end+1} = sprintf('row lines: mean over animals, %s', normtxt); %#ok<AGROW>
        ttl{end+1} = sprintf('(%s%s)', styletxt, scaletxt); %#ok<AGROW>
    end
    if do_stats
        ttl{end+1} = sprintf('* = rank-sum %s vs %s per cell (animals pooled), %s alpha = %.2f', condA, condB, ternary(use_fdr,'BH-FDR,',''), alpha); %#ok<AGROW>
    end
    figs(end+1) = figure('Name',sprintf('%s: %s (%s), rows',feat,what,scorename),'Color','w', ...
        'Units','pixels','Position',[100 60 900 850]); %#ok<AGROW>
    tl = tiledlayout(figs(end),4,1,'TileSpacing','compact');
    title(tl,ttl,'Interpreter','none');

    %tile 1: shank scores of the chosen label, one line per animal
    ax1 = nexttile(tl,1);
    hold(ax1,'on');
    for a = 1:na
        plot(ax1,1:5,scores(a,:),'-o','LineWidth',1.5,'Color',acol(a,:),'MarkerFaceColor',acol(a,:),'DisplayName',animals(a));
    end
    xlim(ax1,[0.5 5.5]);
    xticks(ax1,1:5); xticklabels(ax1,{});
    ylabel(ax1,scorelab,'Interpreter','none');
    grid(ax1,'on'); ax1.GridAlpha = 0.15; box(ax1,'off');
    lg = legend(ax1,'Interpreter','none');
    lg.Layout.Tile = 'east'; %outside both plots (as the colorbar), so line points stay above their heatmap columns

    %tiles 2-4: heatmap, row lines on top, n per column below
    ax2 = nexttile(tl,2,[3 1]);
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
    hold(ax2,'on');
    if shownorm
        %with condB, A points sit slightly left and B points slightly right of the column centre
        dx = 0;
        if ~isempty(rowB_norm), dx = 0.1; end
        for r = 1:nrow
            plot(ax2,[0.5 5.5],[r r],':','Color',[0.35 0.35 0.35],'LineWidth',0.5); %reference value 1
            %white underlay so the black lines stay visible on dark cells
            yA = r - (rowA_norm(r,:) - 1)*k;
            plot(ax2,1:5,yA,'-','Color','w','LineWidth',3);
            if ~isempty(rowB_norm)
                yB = r - (rowB_norm(r,:) - 1)*k;
                plot(ax2,1:5,yB,'-','Color','w','LineWidth',3);
            end
            %each animal's own value: o = condA, x = condB
            for a = 1:na
                plot(ax2,(1:5)-dx, r - (squeeze(ptsA_norm(a,r,:))' - 1)*k, 'o','Color','k', ...
                    'MarkerSize',4,'MarkerFaceColor','w','LineWidth',0.8);
                if ~isempty(ptsB_norm)
                    plot(ax2,(1:5)+dx, r - (squeeze(ptsB_norm(a,r,:))' - 1)*k, 'x','Color','k', ...
                        'MarkerSize',6,'LineWidth',1.2);
                end
            end
            %mean lines on top
            plot(ax2,1:5,yA,'-','Color','k','LineWidth',1.5);
            if ~isempty(rowB_norm)
                plot(ax2,1:5,yB,'--','Color','k','LineWidth',1.5);
            end
        end
    end
    if do_stats
        [rs,cs] = find(sig);
        for i = 1:numel(rs)
            text(ax2, cs(i)+0.36, rs(i)-0.32, '*', 'Color','k', 'FontSize',14, 'FontWeight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
    xlim(ax2,[0.5 5.5]); ylim(ax2,[0.5 nrow+0.5]);
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
    linkaxes([ax1 ax2],'x');

    data(fi).feature = feat;
    data(fi).condA = condA;
    data(fi).condB = condB;
    data(fi).rankfield = rankfield;
    data(fi).region = region;
    data(fi).trialfilter = trialfilter;
    data(fi).animals = animals;
    data(fi).topshank = topshank;
    data(fi).scores = scores;
    data(fi).mapA = mapA;
    data(fi).mapB = mapB;
    data(fi).heat = heat;
    data(fi).n = n;
    data(fi).rowA = rowA;
    data(fi).rowB = rowB;
    data(fi).normto = normto;
    data(fi).rowA_norm = rowA_norm;
    data(fi).rowB_norm = rowB_norm;
    data(fi).ptsA_norm = ptsA_norm;
    data(fi).ptsB_norm = ptsB_norm;
    data(fi).scores_raw = scores_raw;
    data(fi).scorenorm = scorenorm;
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

function ref = rowRef(C, normto)
%reference value of each row of C (rows x 5): its value at column normto, or its peak (max) if normto is "peak"
if isnumeric(normto)
    ref = C(:,normto);
else
    ref = max(C,[],2,'omitnan');
end
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
