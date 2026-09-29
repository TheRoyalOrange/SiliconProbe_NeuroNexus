% plotLWIntensityCompare.m
%
% Description: Question-specific plot comparing light-only (L) and light+whisker (LW)
%   responses across light intensities, per channel on the probe map, next to the
%   probe's activity labels. For one animal, asks which region (from the unique values
%   of the table's Region column for that animal; the probe number is then found from
%   ProbeInfo.Areas and used to index ProbeMaps, Ch_Remove and the Shank/Chan labels)
%   and which data features (measure columns of the summary table) to plot, then makes one figure
%   per feature, a 5 x 3 tiled layout:
%     row 1: activity labels from ProbeInfo, drawn by plotProbeLabels.m (shank line +
%       channel heatmap, same colours as there): (1,1) spontaneous label, (1,2) and (1,3) two
%       stimulus activity scores chosen by the user (recommended: light only, whisker only)
%     rows 2-5: channel heatmaps of the feature for 4, 8, 12 and 15 (light intensity):
%       column 1 = L_<n>, column 2 = LW_<n> (green, one colour scale shared by all of
%       rows 2-5, columns 1-2), column 3 = LW_<n> minus L_<n> (blue-white-red, white = 0,
%       red > 0, blue < 0; one symmetric scale shared by rows 2-5). In column 3, a star marks
%       channels where LW and L differ significantly: Wilcoxon rank-sum test (ranksum) of the
%       channel's LW trial values vs its L trial values (independent trials, NaN ignored,
%       at least min_trials per condition), alpha = 0.05, Benjamini-Hochberg FDR correction
%       across the tested channels of that tile (use_fdr); the tile title gives
%       "<significant>/<tested> *". alpha, use_fdr and min_trials are set at the top of the function
%   Each heatmap cell is a channel in the ProbeMaps layout (columns = shanks, row 1 at the
%   top); its value is the mean of the feature over all of that animal's good trials in
%   that condition (all files), ignoring NaN. Channels with no data (e.g. Ch_Remove) are
%   black; removed channels are marked 'X'. Returns the figure handles; saves nothing.
%   Specific to this comparison (condition names are fixed); may later become a general
%   plotting function.
%
% Inputs:
%   animal (char/string) - animal ID, as in summaryTable.Animal_Name and <animal>-ProbeInfo.mat
%   summaryTable (table, or char/string path to its csv) - tableres from
%     SummaryAnalysis_MultiMouse_allChannels.m (or a single-type all-channels script), one row
%     per good trial x channel. Columns used:
%       Animal_Name, Region, Condition_Name (text), Channel_ID (double)
%       measure columns named <P1|P2|All>_<LFP|MUA|TF>... (double) - the features offered
%     Must contain rows of this animal for conditions L_4, LW_4, L_8, LW_8, L_12, LW_12,
%     L_15, LW_15, else an error
%   region (char/string/number, optional) - region to use, as stored in the Region column: a region
%     name matched in ProbeInfo.Areas (e.g. "V1") or a probe number used directly as the index (e.g. 1);
%     omitted/empty = list dialog of the unique regions for this animal in the table
%   features (string array, optional) - measure columns to plot; omitted/empty = list dialog
%     (multiple selection) of the measure columns with data for this animal and region
%   stimlabels (string array, 1 x 2, optional) - labels for row 1 columns 2 and 3; omitted/empty
%     = two list dialogs (one per column) of the probe's labels other than the spontaneous one,
%     with the light-only (name contains "LightOnly") and whisker-only (contains "Whisker" but
%     not "Light") labels pre-selected when present. Cancel leaves that tile with a note
%   do_stats (logical, optional) - true (default): run the per-channel LW vs L test and star
%     significant channels in column 3; false: no test, no stars, no count in the tile titles
%   E:\Roy\Processed Silicon Probe Data\ProbeInfo\<animal>-ProbeInfo.mat - ProbeInfo with
%     .Areas (region -> probe), .ProbeMaps, .Ch_Remove, and the activity labels in
%     .ShankLabels / .ChanLabels (or .chLabels); the row-1 column-1 label is the first whose name
%     contains "Spon". A missing label leaves its tile with a note
%
% Outputs:
%   figs (figure handle array, 1 x nFeatures) - one figure per feature
%
% Dependencies: plotProbeLabels.m (same folder); SummaryAnalysis_MultiMouse_allChannels.m
%   (or the single-type all-channels scripts) for summaryTable; ChannelShankLabelling_SponData.m /
%   ChannelShankLabelling_StimData.m for the labels. MATLAB R2020b+ (nested tiledlayout);
%   Statistics and Machine Learning Toolbox (ranksum).

function figs = plotLWIntensityCompare(animal, summaryTable, region, features, stimlabels, do_stats)

figs = gobjects(1,0);
animal = char(animal);
intensities = [4 8 12 15];
%per-channel test of LW vs L (rows 2-5, column 3): Wilcoxon rank-sum on the trial values
if nargin < 6 || isempty(do_stats)
    do_stats = true; %default: run the test and mark significant channels
end
alpha = 0.05;       %significance level
use_fdr = true;     %Benjamini-Hochberg FDR correction across the channels of each tile
min_trials = 3;     %minimum trials (non-NaN) per condition to test a channel

%% summary table: from a csv path or given directly; check the columns needed
if ischar(summaryTable) || isstring(summaryTable)
    summaryTable = readtable(char(summaryTable), 'TextType','string');
end
assert(istable(summaryTable), 'summaryTable must be a table (tableres) or the path to its csv')
needcols = ["Animal_Name","Region","Condition_Name","Channel_ID"];
missingcols = needcols(~ismember(needcols, summaryTable.Properties.VariableNames));
assert(isempty(missingcols), 'summaryTable is missing column(s): %s', strjoin(missingcols,', '))
T = summaryTable(string(summaryTable.Animal_Name) == string(animal), :);
assert(height(T) > 0, 'summaryTable has no rows for animal %s', animal)

%% which region (from the Region column; the probe is found from it in ProbeInfo.Areas below)
regions = unique(string(T.Region),'stable');
if nargin < 3 || isempty(region)
    [sel,ok] = listdlg('ListString',cellstr(regions),'SelectionMode','single', ...
        'PromptString',sprintf('%s: which region?',animal),'ListSize',[250 120]);
    assert(ok && ~isempty(sel), 'No region selected')
    region = regions(sel);
end
region = string(region);
assert(ismember(region, regions), '%s has no rows for region %s in summaryTable', animal, region)
T = T(string(T.Region) == region, :);

%% conditions needed (fixed for this comparison)
condL = "L_" + intensities;
condLW = "LW_" + intensities;
condpresent = unique(string(T.Condition_Name));
missingcond = setdiff([condL condLW], condpresent, 'stable');
assert(isempty(missingcond), '%s (%s): summaryTable has no rows for condition(s) %s', animal, region, strjoin(missingcond,', '))

%% which features: measure columns with data for this animal and region
vn = string(T.Properties.VariableNames);
measurecols = vn(~cellfun(@isempty, regexp(cellstr(vn), '^(P1|P2|All)_(LFP|MUA|TF)', 'once')));
hasdata = arrayfun(@(c) isnumeric(T.(c)) && any(~isnan(T.(c))), measurecols);
measurecols = measurecols(hasdata);
assert(~isempty(measurecols), 'summaryTable has no LFP/MUA/TF measure columns with data for %s (%s)', animal, region)
if nargin < 4 || isempty(features)
    [sel,ok] = listdlg('ListString',cellstr(measurecols),'SelectionMode','multiple', ...
        'PromptString','Data features to plot (one figure each):','ListSize',[300 300]);
    assert(ok && ~isempty(sel), 'No data feature selected')
    features = measurecols(sel);
end
features = string(features);
badfeat = setdiff(features, measurecols, 'stable');
assert(isempty(badfeat), 'Feature(s) not in summaryTable or without data for %s (%s): %s', animal, region, strjoin(badfeat,', '))

%% ProbeInfo: probe number of the chosen region (index into ProbeMaps, Ch_Remove and the
%% Shank/Chan labels), layout, removed channels, row-1 labels
pinfo = load(fullfile(['E:\Roy\Processed Silicon Probe Data\ProbeInfo\' animal '-ProbeInfo.mat']),'ProbeInfo');
ProbeInfo = pinfo.ProbeInfo;
%the Region value can be a region name (matched in ProbeInfo.Areas) or a probe number (used as the index)
prb = find(strcmp(ProbeInfo.Areas, region));
if isempty(prb)
    pnum = str2double(region);
    if ~isnan(pnum) && pnum == round(pnum) && pnum >= 1 && pnum <= numel(ProbeInfo.ProbeMaps)
        prb = pnum;
    end
end
assert(numel(prb) == 1, '%s: region "%s" is neither one of ProbeInfo.Areas (%s) nor a probe number 1-%d', ...
    animal, region, strjoin(string(ProbeInfo.Areas),', '), numel(ProbeInfo.ProbeMaps))
areaname = string(ProbeInfo.Areas{prb});
pmap = ProbeInfo.ProbeMaps{prb};
removed = false(size(pmap));
if isfield(ProbeInfo,'Ch_Remove') && numel(ProbeInfo.Ch_Remove) >= prb
    removed = ismember(pmap, ProbeInfo.Ch_Remove{prb});
end
labs = strings(1,0);
for fn = ["ShankLabels","ChanLabels","chLabels"]
    if isfield(ProbeInfo,fn) && iscell(ProbeInfo.(fn)) && numel(ProbeInfo.(fn)) >= prb && isstruct(ProbeInfo.(fn){prb})
        labs = [labs string(fieldnames(ProbeInfo.(fn){prb}))']; %#ok<AGROW>
    end
end
labs = unique(labs,'stable');
%column 1: the spontaneous label (name containing "Spon")
row1lab = [firstMatch(labs, contains(labs,'Spon','IgnoreCase',true)), string(missing), string(missing)];
%columns 2-3: stimulus activity scores chosen by the user (label names differ between recordings);
%the recommended ones (light only / whisker only) are pre-selected when their names match
stimopts = labs(~ismember(labs, row1lab(1)));
if nargin >= 5 && ~isempty(stimlabels)
    stimlabels = string(stimlabels);
    assert(numel(stimlabels) == 2, 'stimlabels must name 2 labels (row 1, columns 2 and 3)')
    badlab = setdiff(stimlabels, stimopts, 'stable');
    assert(isempty(badlab), '%s (probe %d): label(s) not in ProbeInfo: %s (available: %s)', ...
        animal, prb, strjoin(badlab,', '), strjoin(stimopts,', '))
    row1lab(2:3) = stimlabels;
elseif ~isempty(stimopts)
    rec = {find(contains(stimopts,'LightOnly','IgnoreCase',true),1), ...
           find(contains(stimopts,'Whisker','IgnoreCase',true) & ~contains(stimopts,'Light','IgnoreCase',true),1)};
    prompts = {{'Row 1, column 2: stimulus activity score','(recommended: light only)'}, ...
               {'Row 1, column 3: stimulus activity score','(recommended: whisker only)'}};
    for c = 1:2
        init = 1;
        if ~isempty(rec{c}), init = rec{c}; end
        [sel,ok] = listdlg('ListString',cellstr(stimopts),'SelectionMode','single','InitialValue',init, ...
            'PromptString',prompts{c},'ListSize',[300 150]);
        if ok && ~isempty(sel)
            row1lab(c+1) = stimopts(sel);
        end
    end
end
row1name = ["spontaneous (Spon)","column 2 stimulus","column 3 stimulus"];
fprintf('%s (%s, probe %d) row-1 labels: %s\n', animal, areaname, prb, strjoin(fillmissing(row1lab,'constant',"none"),', '))

%colour maps: green as plotProbeLabels' stimulus-driven labels; diverging blue-white-red for LW - L
greenmap = interp1([0 0.5 1],[0.95 0.98 0.93; 0.55 0.82 0.40; 0.18 0.60 0.22],linspace(0,1,256));
divmap = interp1([0 0.5 1],[0.13 0.40 0.75; 1 1 1; 0.80 0.15 0.15],linspace(0,1,256));

%% one figure per feature
for fi = 1:numel(features)
    feat = features(fi);

    %per-channel mean of the feature for each condition, in the ProbeMaps layout
    mapL = cell(1,numel(intensities));  mapLW = cell(1,numel(intensities));
    for n = 1:numel(intensities)
        mapL{n} = channelMap(T, feat, condL(n), pmap);
        mapLW{n} = channelMap(T, feat, condLW(n), pmap);
    end
    mapD = cellfun(@(lw,l) lw - l, mapLW, mapL, 'UniformOutput', false);
    %per-channel significance of LW vs L (star in the difference maps), if do_stats
    sigD = cell(1,numel(intensities));  ntested = zeros(1,numel(intensities));
    if do_stats
        for n = 1:numel(intensities)
            pmapP = channelTest(T, feat, condL(n), condLW(n), pmap, min_trials);
            [sigD{n}, ntested(n)] = significant(pmapP, alpha, use_fdr);
        end
    end
    allLLW = [cell2mat(mapL) cell2mat(mapLW)];
    climG = [min(allLLW,[],'all','omitnan') max(allLLW,[],'all','omitnan')];
    if any(isnan(climG)), climG = [0 1]; end
    if climG(1) == climG(2), climG = climG(1) + [-1 1]; end
    dmax = max(abs(cell2mat(mapD)),[],'all','omitnan');
    if isempty(dmax) || isnan(dmax) || dmax == 0, dmax = 1; end
    climD = [-dmax dmax];

    ttl = sprintf('%s - probe %d (%s): %s, L vs LW by light intensity', animal, prb, areaname, feat);
    figs(end+1) = figure('Name',ttl,'Color','w','Units','normalized','OuterPosition',[0.05 0.03 0.6 0.95]); %#ok<AGROW>
    tl = tiledlayout(figs(end),5,3,'TileSpacing','compact','Padding','compact');
    title(tl,ttl,'Interpreter','none');

    %row 1: activity labels, drawn by plotProbeLabels into tiles 1-3
    for c = 1:3
        if ismissing(row1lab(c))
            ax = nexttile(tl,c);
            text(ax,0.5,0.5,sprintf('no %s label\nin ProbeInfo',row1name(c)),'Units','normalized','HorizontalAlignment','center');
            axis(ax,'off');
        else
            plotProbeLabels(ProbeInfo, prb, row1lab(c), tl, c);
        end
    end

    %rows 2-5: L, LW and LW - L per light intensity
    for n = 1:numel(intensities)
        base = 3*n; %tile before this row
        axL = drawMap(nexttile(tl,base+1), mapL{n}, removed, greenmap, climG, condL(n));
        axLW = drawMap(nexttile(tl,base+2), mapLW{n}, removed, greenmap, climG, condLW(n));
        if do_stats
            dttl = sprintf('%s - %s (%d/%d *)', condLW(n), condL(n), nnz(sigD{n}), ntested(n));
        else
            dttl = sprintf('%s - %s', condLW(n), condL(n));
        end
        axD = drawMap(nexttile(tl,base+3), mapD{n}, removed, divmap, climD, dttl, sigD{n});
        if n == numel(intensities)
            cb = colorbar(axLW,'southoutside'); cb.Label.String = feat + " (L, LW; channel mean over trials)"; cb.Label.Interpreter = 'none';
            cb = colorbar(axD,'southoutside');  cb.Label.String = feat + " (LW - L)"; cb.Label.Interpreter = 'none';
            xlabel(axL,'Shank (ProbeMaps column)');
        end
        ylabel(axL,'ProbeMaps row');
    end
end
end

function lab = firstMatch(labs, mask)
%first label where mask is true, or missing
idx = find(mask,1);
if isempty(idx)
    lab = string(missing);
else
    lab = labs(idx);
end
end

function M = channelMap(T, feat, cond, pmap)
%mean of feature feat over the rows of condition cond, per channel, in the ProbeMaps layout (NaN = no data)
M = NaN(size(pmap));
sub = string(T.Condition_Name) == cond;
[ids,~,g] = unique(T.Channel_ID(sub));
vals = accumarray(g, T.(feat)(sub), [], @(x) mean(x,'omitnan'));
[onmap,loc] = ismember(ids, pmap);
M(loc(onmap)) = vals(onmap);
end

function P = channelTest(T, feat, condA, condB, pmap, min_trials)
%Wilcoxon rank-sum p-value of feature feat between the trials of condA and condB, per channel, in
%the ProbeMaps layout (NaN = not tested: channel missing or fewer than min_trials per condition)
P = NaN(size(pmap));
inA = string(T.Condition_Name) == condA;
inB = string(T.Condition_Name) == condB;
for ch = unique(T.Channel_ID(inA | inB))'
    a = T.(feat)(inA & T.Channel_ID == ch);  a = a(~isnan(a));
    b = T.(feat)(inB & T.Channel_ID == ch);  b = b(~isnan(b));
    loc = find(pmap == ch, 1);
    if isempty(loc) || numel(a) < min_trials || numel(b) < min_trials, continue, end
    P(loc) = ranksum(a, b);
end
end

function [sig, ntested] = significant(P, alpha, use_fdr)
%significant channels (logical, same layout as P); Benjamini-Hochberg FDR over the tested channels if use_fdr
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

function ax = drawMap(ax, M, removed, cmap, lims, ttl, sig)
%heatmap of M in the ProbeMaps layout: NaN cells black, removed channels marked 'X',
%significant channels (optional logical map sig) marked with a star
[nrow,nshnk] = size(M);
imagesc(ax,1:nshnk,1:nrow,M,'AlphaData',~isnan(M));
ax.Color = 'k';
colormap(ax,cmap);
clim(ax,lims);
xticks(ax,1:nshnk); yticks(ax,1:nrow);
title(ax,ttl,'Interpreter','none','FontWeight','normal');
[r,c] = find(removed);
for i = 1:numel(r)
    text(ax,c(i),r(i),'X','HorizontalAlignment','center','Color','w','FontSize',7);
end
if nargin >= 7 && ~isempty(sig)
    [r,c] = find(sig);
    for i = 1:numel(r)
        text(ax,c(i),r(i),'*','HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Color','k','FontSize',11,'FontWeight','bold');
    end
end
end
