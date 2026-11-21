% NNx_plotPairedCond_intensitycompare.m
%
% Description: General version of NNx_plotLWIntensityCompare.m. Compares two conditions
%   (column 1 = condition A, column 2 = condition B) per channel on the probe map, next to the
%   probe's activity labels, with one figure row per condition pair. The rows are whatever pairs
%   you give (or pick): light intensities (L_4 vs LW_4, L_8 vs LW_8, ...), light positions,
%   whisker intensities, etc. For one animal, asks which region (from the unique values of the
%   table's Region column for that animal; the probe number is then found from ProbeInfo.Areas
%   and used to index ProbeMaps, Ch_Remove and the Shank/Chan labels), which condition pairs
%   and which data features (measure columns of the summary table) to plot, then makes one
%   figure per feature, a (1 + nPairs) x 3 tiled layout:
%     row 1: activity labels from ProbeInfo, drawn by NNx_plotProbeLabels.m (shank line +
%       channel heatmap, same colours as there): (1,1) spontaneous label, (1,2) and (1,3) two
%       stimulus activity scores chosen by the user (light only / whisker only are pre-selected
%       when present, just as a suggestion)
%     rows 2 to 1+nPairs: channel heatmaps of the feature, one row per condition pair:
%       column 1 = condition A, column 2 = condition B (green, one colour scale shared by all
%       pair rows, columns 1-2), column 3 = B minus A (blue-white-red, white = 0, red > 0,
%       blue < 0; one symmetric scale shared by all pair rows). In column 3, a star marks
%       channels where B and A differ significantly: Wilcoxon rank-sum test (ranksum) of the
%       channel's B trial values vs its A trial values (independent trials, NaN ignored,
%       at least min_trials per condition), alpha = 0.05, Benjamini-Hochberg FDR correction
%       across the tested channels of that tile (use_fdr); the tile title gives
%       "<significant>/<tested> *". alpha, use_fdr and min_trials are set at the top of the function
%   Each heatmap cell is a channel in the ProbeMaps layout (columns = shanks, row 1 at the
%   top); its value is the mean of the feature over all of that animal's good trials in
%   that condition (all files), ignoring NaN. Channels with no data (e.g. Ch_Remove) are
%   black; removed channels are marked 'X'. Returns the figure handles; with save_png the
%   figures are also saved as PNGs (see Outputs), otherwise nothing is saved.
%
% Inputs:
%   animal (char/string) - animal ID, as in summaryTable.Animal_Name and <animal>-ProbeInfo.mat
%   summaryTable (table, or char/string path to its csv) - tableres from
%     NNx_SummaryAnalysis_MultiMouse_allChannels.m (or a single-type all-channels script), one row
%     per good trial x channel. Columns used:
%       Animal_Name, Region, Condition_Name (text), Channel_ID (double)
%       measure columns named <P1|P2|All>_<LFP|MUA|TF>... (double) - the features offered
%   region (char/string/number, optional) - region to use, as stored in the Region column: a region
%     name matched in ProbeInfo.Areas (e.g. "V1") or a probe number used directly as the index (e.g. 1);
%     omitted/empty = list dialog of the unique regions for this animal in the table
%   condpairs (string array, nPairs x 2, optional) - one figure row per pair, top to bottom;
%     column 1 = condition A, column 2 = condition B, named exactly as in Condition_Name, e.g.
%       ["L_4" "LW_4"; "L_8" "LW_8"; "L_12" "LW_12"; "L_15" "LW_15"]   (light intensities)
%       ["W_1" "LW_1"; "W_2" "LW_2"]                                     (whisker intensities, say)
%     omitted/empty = two list dialogs (multiple selection) of this animal + region's conditions:
%     first the column-1 conditions, then the same number of column-2 conditions; they are paired
%     in list order (the list is sorted by name prefix, then by the trailing number, so L_4 comes
%     before L_12). For any other order, pass condpairs. A condition without rows for this animal
%     is left NaN: its heatmap (and the B - A map of that row) is black with "No Data". Error only
%     if none of the conditions has rows
%   features (string array, optional) - measure columns to plot; omitted/empty = list dialog
%     (multiple selection) of the measure columns with data for this animal and region
%   stimlabels (string array, 1 x 2, optional) - labels for row 1 columns 2 and 3; omitted/empty
%     = two list dialogs (one per column) of the probe's labels other than the spontaneous one,
%     with the light-only (name contains "LightOnly") and whisker-only (contains "Whisker" but
%     not "Light") labels pre-selected when present. Cancel leaves that tile with a note
%   do_stats (logical, optional) - true (default): run the per-channel B vs A test and star
%     significant channels in column 3; false: no test, no stars, no count in the tile titles
%   save_png (logical, optional) - true: also save each figure as a PNG (see Outputs);
%     default false. Once the region is known and before anything is plotted, if the output
%     folder already exists a dialog
%     offers Overwrite (use it; same-name files are replaced), Quit (stops with an error) or
%     Edit Folder Name (a text box pre-filled with the current name; the new name is checked
%     the same way, until it is new or Overwrite / Quit is chosen). Closing the dialog = Quit
%   comparename (char/string, optional) - short name of the comparison, used in the figure titles
%     and the save folder / file names (e.g. "LvsLW_intensity", "LvsLW_position"). omitted/empty =
%     "<A prefix>vs<B prefix>" when all column-1 conditions share one prefix and all column-2
%     conditions share one (prefix = name without a trailing "_<number>", so L_4/L_8 + LW_4/LW_8
%     gives "LvsLW"), otherwise "PairedCond"
%   E:\Roy\Processed Silicon Probe Data\ProbeInfo\<animal>-ProbeInfo.mat - ProbeInfo with
%     .Areas (region -> probe), .ProbeMaps, .Ch_Remove, and the activity labels in
%     .ShankLabels / .ChanLabels (or .chLabels); the row-1 column-1 label is the first whose name
%     contains "Spon". A missing label leaves its tile with a note
%
% Outputs:
%   figs (figure handle array, 1 x nFeatures) - one figure per feature
%   with save_png: E:\Roy\Processed Silicon Probe Data\BundledAnimalData\SingleAnimal\
%     <folder>\<animal>_<area>_<comparename>_<feature>.png, one per feature, 150 dpi. <folder> is
%     <animal>_<area>_<comparename>, or the name entered under Edit Folder Name (created if missing).
%     <area> is the probe's name in ProbeInfo.Areas (also when region is given as a probe number)
%
% Dependencies: NNx_plotProbeLabels.m (same folder); NNx_SummaryAnalysis_MultiMouse_allChannels.m
%   (or the single-type all-channels scripts) for summaryTable; NNx_editProbeInfo_ChLabel_SponData.m /
%   NNx_editProbeInfo_ChLabel_StimData.m for the labels. MATLAB R2020b+ (nested tiledlayout);
%   Statistics and Machine Learning Toolbox (ranksum).

function figs = NNx_plotPairedCond_intensitycompare(animal, summaryTable, region, condpairs, features, stimlabels, do_stats, save_png, comparename)

figs = gobjects(1,0);
animal = char(animal);
%per-channel test of B vs A (column 3): Wilcoxon rank-sum on the trial values
if nargin < 7 || isempty(do_stats)
    do_stats = true; %default: run the test and mark significant channels
end

%save_png: output folder <savebase>\<animal>_<area>_<comparename>, chosen once the region is known (below)
if nargin < 8 || isempty(save_png)
    save_png = false; %default: figures are only returned
end
savebase = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\SingleAnimal';
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

%% which condition pairs (one figure row each): given, or picked from this animal + region's conditions
condpresent = unique(string(T.Condition_Name));
if nargin < 4 || isempty(condpairs)
    %sort by name prefix, then by the trailing number, so L_4, L_8, L_12, L_15 (not L_12, L_15, L_4, L_8)
    condnum = str2double(regexp(condpresent,'[\d.]+$','match','once')); %NaN if no trailing number
    condpre = regexprep(condpresent,'_?[\d.]+$','');
    [~,ord] = sortrows(table(condpre(:),condnum(:)));
    condlist = condpresent(ord);
    [selA,ok] = listdlg('ListString',cellstr(condlist),'SelectionMode','multiple', ...
        'PromptString',{'Column 1 (condition A), one per figure row:','(rows go in list order)'},'ListSize',[300 300]);
    assert(ok && ~isempty(selA), 'No column-1 condition selected')
    [selB,ok] = listdlg('ListString',cellstr(condlist),'SelectionMode','multiple', ...
        'PromptString',{sprintf('Column 2 (condition B): pick %d,',numel(selA)),'paired with column 1 in list order'},'ListSize',[300 300]);
    assert(ok && numel(selB) == numel(selA), 'Pick the same number of column-2 conditions as column-1 conditions (%d)', numel(selA))
    condpairs = [condlist(selA(:)) condlist(selB(:))];
end
condpairs = string(condpairs);
assert(size(condpairs,2) == 2, 'condpairs must be nPairs x 2 (column 1 = condition A, column 2 = condition B)')
condA = condpairs(:,1)'; %[1 x nPairs] column-1 condition of each row
condB = condpairs(:,2)'; %[1 x nPairs] column-2 condition of each row
npairs = numel(condA);
%a condition without rows is plotted as an all-NaN (black) map labelled "No Data"; error only if none has data
assert(any(ismember([condA condB], condpresent)), '%s (%s): summaryTable has no rows for any of the conditions %s', ...
    animal, region, strjoin([condA condB],', '))

%comparison name for titles and saving: "<A prefix>vs<B prefix>" when each column has one prefix
if nargin < 9 || isempty(comparename)
    preA = unique(regexprep(condA,'_?[\d.]+$',''));
    preB = unique(regexprep(condB,'_?[\d.]+$',''));
    if numel(preA) == 1 && numel(preB) == 1 && preA ~= preB
        comparename = preA + "vs" + preB;
    else
        comparename = "PairedCond";
    end
end
comparename = char(comparename);

%% which features: measure columns with data for this animal and region
vn = string(T.Properties.VariableNames);
measurecols = vn(~cellfun(@isempty, regexp(cellstr(vn), '^(P1|P2|All)_(LFP|MUA|TF)', 'once')));
hasdata = arrayfun(@(c) isnumeric(T.(c)) && any(~isnan(T.(c))), measurecols);
measurecols = measurecols(hasdata);
assert(~isempty(measurecols), 'summaryTable has no LFP/MUA/TF measure columns with data for %s (%s)', animal, region)
if nargin < 5 || isempty(features)
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

%% save_png: choose the output folder before anything is plotted (needs the area name)
savedir = '';
if save_png
    foldername = [animal '_' char(areaname) '_' comparename];
    %if the folder exists: Overwrite (use it, same-name files are replaced), Quit, or Edit Folder Name
    %(the new name is checked again, until it is new or the user chooses Overwrite / Quit)
    while isfolder(fullfile(savebase, foldername))
        choice = questdlg(sprintf('The folder\n%s\nalready exists.', fullfile(savebase, foldername)), ...
            'Folder already exists', 'Overwrite', 'Quit', 'Edit Folder Name', 'Quit');
        switch choice
            case 'Overwrite'
                break
            case 'Edit Folder Name'
                answer = inputdlg('New folder name:', 'Edit folder name', [1 60], {foldername});
                if ~isempty(answer) && ~isempty(strtrim(answer{1}))
                    foldername = strtrim(answer{1}); %checked again by the while condition
                end
            otherwise %Quit, or the dialog was closed
                error('NNx_plotPairedCond_intensitycompare: stopped by the user (output folder %s already exists)', fullfile(savebase, foldername))
        end
    end
    savedir = fullfile(savebase, foldername);
    if ~isfolder(savedir)
        mkdir(savedir);
    end
end

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
%light only / whisker only are pre-selected when their names match
stimopts = labs(~ismember(labs, row1lab(1)));
if nargin >= 6 && ~isempty(stimlabels)
    stimlabels = string(stimlabels);
    assert(numel(stimlabels) == 2, 'stimlabels must name 2 labels (row 1, columns 2 and 3)')
    badlab = setdiff(stimlabels, stimopts, 'stable');
    assert(isempty(badlab), '%s (probe %d): label(s) not in ProbeInfo: %s (available: %s)', ...
        animal, prb, strjoin(badlab,', '), strjoin(stimopts,', '))
    row1lab(2:3) = stimlabels;
elseif ~isempty(stimopts)
    rec = {find(contains(stimopts,'LightOnly','IgnoreCase',true),1), ...
           find(contains(stimopts,'Whisker','IgnoreCase',true) & ~contains(stimopts,'Light','IgnoreCase',true),1)};
    prompts = {{'Row 1, column 2: stimulus activity score','(pre-selected: light only, if present)'}, ...
               {'Row 1, column 3: stimulus activity score','(pre-selected: whisker only, if present)'}};
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

%colour maps: green as NNx_plotProbeLabels' stimulus-driven labels; diverging blue-white-red for B - A
greenmap = interp1([0 0.5 1],[0.95 0.98 0.93; 0.55 0.82 0.40; 0.18 0.60 0.22],linspace(0,1,256));
divmap = interp1([0 0.5 1],[0.13 0.40 0.75; 1 1 1; 0.80 0.15 0.15],linspace(0,1,256));

%% one figure per feature
for fi = 1:numel(features)
    feat = features(fi);

    %per-channel mean of the feature for each condition, in the ProbeMaps layout
    mapA = cell(1,npairs);  mapB = cell(1,npairs); %cell {1 x nPairs} [ProbeMaps rows x shanks]
    for n = 1:npairs
        mapA{n} = channelMap(T, feat, condA(n), pmap);
        mapB{n} = channelMap(T, feat, condB(n), pmap);
    end
    mapD = cellfun(@(b,a) b - a, mapB, mapA, 'UniformOutput', false);
    %per-channel significance of B vs A (star in the difference maps), if do_stats
    sigD = cell(1,npairs);  ntested = zeros(1,npairs);
    if do_stats
        for n = 1:npairs
            pmapP = channelTest(T, feat, condA(n), condB(n), pmap, min_trials);
            [sigD{n}, ntested(n)] = significant(pmapP, alpha, use_fdr);
        end
    end
    allAB = [cell2mat(mapA) cell2mat(mapB)];
    climG = [min(allAB,[],'all','omitnan') max(allAB,[],'all','omitnan')];
    if any(isnan(climG)), climG = [0 1]; end
    if climG(1) == climG(2), climG = climG(1) + [-1 1]; end
    dmax = max(abs(cell2mat(mapD)),[],'all','omitnan');
    if isempty(dmax) || isnan(dmax) || dmax == 0, dmax = 1; end
    climD = [-dmax dmax];

    ttl = sprintf('%s - probe %d (%s): %s, %s', animal, prb, areaname, feat, comparename);
    figs(end+1) = figure('Name',ttl,'Color','w','Units','normalized','OuterPosition',[0.05 0.03 0.6 0.95]); %#ok<AGROW>
    tl = tiledlayout(figs(end),1+npairs,3,'TileSpacing','compact','Padding','compact');
    title(tl,ttl,'Interpreter','none');

    %row 1: activity labels, drawn by NNx_plotProbeLabels into tiles 1-3
    for c = 1:3
        if ismissing(row1lab(c))
            ax = nexttile(tl,c);
            text(ax,0.5,0.5,sprintf('no %s label\nin ProbeInfo',row1name(c)),'Units','normalized','HorizontalAlignment','center');
            axis(ax,'off');
        else
            NNx_plotProbeLabels(ProbeInfo, prb, row1lab(c), tl, c);
        end
    end

    %rows 2 to 1+nPairs: A, B and B - A per condition pair
    for n = 1:npairs
        base = 3*n; %tile before this row
        axA = drawMap(nexttile(tl,base+1), mapA{n}, removed, greenmap, climG, condA(n));
        axB = drawMap(nexttile(tl,base+2), mapB{n}, removed, greenmap, climG, condB(n));
        if do_stats
            dttl = sprintf('%s - %s (%d/%d *)', condB(n), condA(n), nnz(sigD{n}), ntested(n));
        else
            dttl = sprintf('%s - %s', condB(n), condA(n));
        end
        axD = drawMap(nexttile(tl,base+3), mapD{n}, removed, divmap, climD, dttl, sigD{n});
        if n == npairs
            cb = colorbar(axB,'southoutside'); cb.Label.String = feat + " (columns 1-2; channel mean over trials)"; cb.Label.Interpreter = 'none';
            cb = colorbar(axD,'southoutside'); cb.Label.String = feat + " (column 2 - column 1)"; cb.Label.Interpreter = 'none';
            xlabel(axA,'Shank (ProbeMaps column)');
        end
        ylabel(axA,'ProbeMaps row');
    end

    %save_png: <savedir>\<animal>_<area>_<comparename>_<feature>.png
    if save_png
        exportgraphics(figs(end), fullfile(savedir, [animal '_' char(areaname) '_' comparename '_' char(feat) '.png']), 'Resolution', 150);
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
%significant channels (optional logical map sig) marked with a star; an all-NaN map
%(no data for that condition) is all black with "No Data" in the middle
[nrow,nshnk] = size(M);
imagesc(ax,1:nshnk,1:nrow,M,'AlphaData',~isnan(M));
ax.Color = 'k';
colormap(ax,cmap);
clim(ax,lims);
xticks(ax,1:nshnk); yticks(ax,1:nrow);
title(ax,ttl,'Interpreter','none','FontWeight','normal');
if all(isnan(M),'all')
    text(ax,0.5,0.5,'No Data','Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','Color','w','FontSize',12,'FontWeight','bold');
    return %no removed-channel or significance marks on an empty map
end
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
