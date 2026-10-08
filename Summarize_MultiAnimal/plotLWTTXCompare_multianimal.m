% plotLWTTXCompare_multianimal.m
%
% Description: TTX experiment version of plotLWIntensityCompare_multianimal.m. Compares
%   light-only (L) and light+whisker (LW) responses at 8 uW before and with TTX, for several
%   animals, from the summary table. Each animal's shanks are re-centred on its top-ranked
%   shank for an alignment label (ShankRank_<alignlabel> == 1; chosen in a dialog if
%   omitted): that shank is column 3, its neighbouring shanks columns 1-2 (left) and 4-5
%   (right). Rows are Chan_Depth 1..maxdepth (position among the shank's kept channels,
%   1 = top; maxdepth = 3). One figure per chosen feature, a 5 x 4 tiled layout:
%     row 1: activity labels from each animal's ProbeInfo, all aligned by the alignment label:
%       (1,1) the spontaneous label (name containing "Spon"), (1,2)-(1,3) stimlabels (default
%       LightOnly, WhiskerOnly). Each tile: a line per animal of the label's shank scores for
%       the 5 aligned shanks, and below it a 5 x maxdepth heatmap of the mean over animals of
%       the channel scores (ProbeInfo.ChanLabels, or chLabels). Columns 1-2: each animal's
%       shank and channel scores are divided by its score on the centre shank (lines meet at
%       1 in column 3). Column 3 (WhiskerOnly by default; its scores can be negative or near
%       0): divided by the animal's largest |shank score| over the whole probe, so its
%       largest shank is +1 or -1. Orange for Spon, green for stimulus labels
%     row 2: L_8 | LW_8 | LW_8 / L_8 | -
%     row 3: L_8_TTX | LW_8_TTX | LW_8_TTX / L_8_TTX | -
%     row 4: L_8_TTX / L_8 | LW_8_TTX / LW_8 | - | (LW_8_TTX / LW_8) - (L_8_TTX / L_8)
%     row 5: - | - | (LW_8_TTX / L_8_TTX) - (LW_8 / L_8) | -
%   (each ratio comparison sits next to the two ratio tiles it compares)
%   Single conditions (rows 2-3, columns 1-2): mean over animals of each animal's
%   per-position mean over trials; green, one colour scale shared by the 4 tiles; with
%   normto not "none", each also shows per row a line (mean over animals) and each animal's
%   value (o), divided by the line's reference (normto: centre column by default, or a
%   column / the line's peak), on one vertical scale shared by the 4 tiles (reference value
%   1 on the row's dotted centre line; the title gives the scale).
%   Ratio tiles (B / A): per animal and cell, the animal's mean of B divided by its mean of A
%   (one value per animal; NaN where the A mean is 0); heatmap = mean of the animal ratios.
%   Blue-white-red with white = 1, two scales symmetric around 1: one for the two LW / L
%   tiles (column 3), one for the two TTX / before tiles (row 4, columns 1-2).
%   Ratio comparisons: per animal, the difference of two of its ratio maps; heatmap = mean
%   over animals; blue-white-red with white = 0, one symmetric scale for both tiles.
%   Statistics (do_stats):
%     ratio tiles: per cell, the feature values of all chosen animals' trials at that position
%       (animals pooled): Wilcoxon rank-sum test (ranksum) of condition B vs condition A
%       (unpaired, NaN ignored, at least min_trials per condition)
%     ratio comparisons: per cell, one-sample t-test (ttest) of the animals' differences
%       against 0 (one value per animal; at least min_animals animals)
%     Benjamini-Hochberg FDR over the tested cells of each tile (use_fdr). The adjusted p (q;
%     raw p if use_fdr is false) is written in the centre of each tested cell; significant
%     cells get stars above it (* q < 0.05, ** < 0.01, *** < 0.001). Tile titles give
%     "(<significant>/<tested> *)".
%   A condition without rows for any chosen animal gives an all-black "No Data" tile (also
%   in the ratio / comparison tiles that use it); animals missing a condition simply do not
%   contribute to it. Every heatmap gives n = animals with data per column under its
%   columns. Optional trial filter as in plotRankedShankRows.m. The condition names, alpha,
%   use_fdr, min_trials, min_animals, maxdepth and the paths are set at the top of the
%   function.
%   Returns figure handles and the plotted numbers; with save_png the figures are also saved
%   as PNGs (see Outputs).
%
% Inputs (same order as plotLWIntensityCompare_multianimal.m):
%   animals (string array, optional) - Animal_Name values; omitted/empty = list dialog (multiple)
%   summaryTable (table, or char/string path to its csv) - tableres from
%     SummaryAnalysis_MultiMouse_allChannels.m, one row per good trial x channel. Columns used:
%       Animal_Name, Condition_Name, Region (text)
%       Chan_Shank (double) - shank (ProbeMaps column) of the channel
%       Chan_Depth (double) - position among that shank's kept channels, 1 = top; NaN rows dropped
%       ShankRank_<alignlabel> (double) - rank of the channel's shank by that ShankLabels score
%       measure columns named <P1|P2|All>_... (double) - the features offered
%     Conditions plotted: L_8, LW_8, L_8_TTX, LW_8_TTX (error only if none of them has rows
%     for the chosen animals and region)
%     INFERRED DATA CONTRACT: neighbouring Chan_Shank numbers are physically neighbouring
%     shanks (ProbeMaps column order = physical shank order)
%   region (char/string, optional) - Region value; omitted/empty = the only region of the chosen
%     animals, or a list dialog. Also picks each animal's probe (matched in ProbeInfo.Areas,
%     or used as the probe number)
%   features (string array, optional) - measure columns, one figure each; omitted/empty =
%     list dialog (multiple) of the measure columns with data
%   stimlabels (string array, 1 x 2, optional) - ShankLabels for row 1 columns 2-3; default
%     ["LightOnly","WhiskerOnly"]. Row 1 column 1 is the first label containing "Spon"
%   do_stats (logical, optional) - true (default): test per cell in the 4 ratio tiles and the
%     2 ratio comparisons
%   save_png (logical, optional) - true: also save each figure as a PNG (see Outputs);
%     default false. Once features and region are known and before anything is plotted, a
%     text box asks for a file-name prefix (default n<number of animals>; Cancel stops the
%     function). If any of the files already exists, a dialog offers Overwrite, Quit (stops
%     with an error) or Edit Prefix (the new prefix is checked again). Closing it = Quit
%   alignlabel (char/string, optional) - label whose top-ranked shank is column 3 (a
%     ShankRank_<label> column, prefix optional); omitted/empty = list dialog of the table's
%     ShankRank_ labels (LightOnly pre-selected)
%   normto (double or char/string, optional) - row lines on the single-condition heatmaps: a
%     column 1-5 to normalize to (3 = centre), "peak", or "none" / 0 / false (no row lines or
%     points). Default 3
%   trialfilter (string array, optional) - conditions "<column> <op> <value>", op one of
%     > >= < <= == ~=; a numeric column takes a number, a text column only == / ~= with a text
%     value. Several are combined with AND. Rows where the column is NaN are dropped.
%     Omitted/empty = no filtering
%   <probeinfo_dir>\<animal>-ProbeInfo.mat - each animal's ProbeInfo: .Areas, .ProbeMaps,
%     .Ch_Remove, .ShankLabels{probe}.<label> (1 x shanks), .ChanLabels{probe}.<label>
%     (rows x shanks, ProbeMaps layout; or the older .chLabels). An animal missing a label
%     has no line / heatmap values for it (warning)
%
% Outputs:
%   figs (figure handle array, 1 x nFeatures) - one figure per feature
%   data (struct array, 1 x nFeatures) - the plotted numbers, fields:
%     .feature, .region, .alignlabel, .stimlabels, .trialfilter, .normto
%     .animals (string, 1 x nAnimals), .topshank (double, 1 x nAnimals) - Chan_Shank of each
%       animal's top-ranked shank for alignlabel
%     .labels (struct array, 1 x 3) - row 1: .name, .scores (nAnimals x 5, normalized shank
%       scores: / centre shank for labels 1-2, / max |shank score| for label 3), .chanmaps
%       (nAnimals x maxdepth x 5, channel scores normalized the same way), .heat
%       (maxdepth x 5), .n (1 x 5)
%     .cond (struct array, 1 x 4: L_8, LW_8, L_8_TTX, LW_8_TTX) - .name, .map (nAnimals x
%       maxdepth x 5), .heat (maxdepth x 5), .n (1 x 5)
%     .ratio (struct array, 1 x 4: LW_8/L_8, LW_8_TTX/L_8_TTX, L_8_TTX/L_8, LW_8_TTX/LW_8) -
%       .name ("B / A"), .condB, .condA, .scale (1 = LW / L scale, 2 = TTX / before scale),
%       .map (nAnimals x maxdepth x 5, animal ratios), .heat (maxdepth x 5, mean of them), .n,
%       .p, .q, .sig (maxdepth x 5, trial-level rank-sum; empty without statistics)
%     .compare (struct array, 1 x 2: (LW_8_TTX/LW_8) - (L_8_TTX/L_8), (LW_8_TTX/L_8_TTX) -
%       (LW_8/L_8)) - .name, .first, .second (indices into .ratio), .map (nAnimals x maxdepth
%       x 5, animal differences), .heat, .n, .p, .q, .sig (one-sample t-test vs 0; empty
%       without statistics)
%   with save_png: <savebase>\<prefix>_<area>_LvsLW_TTX_<feature>.png, one per feature,
%     150 dpi (savebase = E:\Roy\Processed Silicon Probe Data\BundledAnimalData\MultAnimal)
%
% Dependencies: SummaryAnalysis_MultiMouse_allChannels.m for summaryTable, including its
%   separate section that adds Chan_Shank / Chan_Depth / ShankRank_<score>; ProbeInfo labels
%   from ChannelShankLabelling_SponData.m / ChannelShank_activitylabelling.m. Same layout
%   ideas and helpers as plotLWIntensityCompare_multianimal.m; heatmaps and statistics as
%   plotRankedShankRows.m. MATLAB R2020b+ (nested tiledlayout); Statistics and Machine
%   Learning Toolbox (ranksum, ttest).

function [figs, data] = plotLWTTXCompare_multianimal(animals, summaryTable, region, features, ...
    stimlabels, do_stats, save_png, alignlabel, normto, trialfilter)

figs = gobjects(1,0);
data = struct([]);
cL = "L_8";  cLW = "LW_8";  cLt = "L_8_TTX";  cLWt = "LW_8_TTX"; %conditions (fixed)
maxdepth = 3;       %rows: Chan_Depth 1..maxdepth
alpha = 0.05;       %significance level
use_fdr = true;     %Benjamini-Hochberg FDR correction across the tested cells of each ratio / comparison tile
min_trials = 3;     %minimum trials (non-NaN) per condition to test a cell
min_animals = 3;    %minimum animals with a value to t-test a ratio-comparison cell
offsets = -2:2;     %columns: shanks relative to the top-ranked shank
linefill = 0.45;    %largest row-line deviation, as a fraction of a row height
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';
savebase = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\MultAnimal';

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
if nargin < 10 || isempty(trialfilter)
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
if nargin < 1 || isempty(animals)
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
if nargin < 3 || isempty(region)
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

%% conditions (fixed): a condition without rows becomes a "No Data" tile
assert(any(ismember([cL cLW cLt cLWt], T.Condition_Name)), 'The chosen animals have no rows for any of %s (%s)', ...
    strjoin([cL cLW cLt cLWt],', '), region)
T = T(ismember(T.Condition_Name, [cL cLW cLt cLWt]), :);

%% alignment label
if nargin < 8 || isempty(alignlabel)
    %list dialog of the table's ShankRank_ labels, LightOnly pre-selected
    rankcols = vn(startsWith(vn,"ShankRank_"));
    assert(~isempty(rankcols), 'summaryTable has no ShankRank_<score> columns')
    init = find(rankcols == "ShankRank_LightOnly", 1);
    if isempty(init), init = 1; end
    [sel,ok] = listdlg('ListString',cellstr(rankcols),'SelectionMode','single','InitialValue',init, ...
        'PromptString','Align on the top shank of which label?','ListSize',[250 150]);
    assert(ok && ~isempty(sel), 'No alignment label selected')
    alignlabel = rankcols(sel);
end
alignlabel = erase(string(alignlabel),"ShankRank_");
rankfield = "ShankRank_" + alignlabel;
assert(ismember(rankfield, vn), 'summaryTable has no column %s (available: %s)', rankfield, strjoin(vn(startsWith(vn,"ShankRank_")),', '))

%% features: measure columns with data
measurecols = vn(~cellfun(@isempty, regexp(cellstr(vn), '^(P1|P2|All)_', 'once')));
hasdata = arrayfun(@(c) isnumeric(T.(c)) && any(~isnan(T.(c))), measurecols);
measurecols = measurecols(hasdata);
assert(~isempty(measurecols), 'summaryTable has no measure columns with data for these animals/conditions')
if nargin < 4 || isempty(features)
    [sel,ok] = listdlg('ListString',cellstr(measurecols),'SelectionMode','multiple', ...
        'PromptString','Data features to plot (one figure each):','ListSize',[300 300]);
    assert(ok && ~isempty(sel), 'No data feature selected')
    features = measurecols(sel);
end
features = string(features(:)');
badfeat = setdiff(features, measurecols, 'stable');
assert(isempty(badfeat), 'Feature(s) not in summaryTable or without data: %s', strjoin(badfeat,', '))

%% other options
if nargin < 5 || isempty(stimlabels)
    stimlabels = ["LightOnly","WhiskerOnly"];
end
stimlabels = string(stimlabels);
assert(numel(stimlabels) == 2, 'stimlabels must name 2 labels (row 1, columns 2 and 3)')
if nargin < 6 || isempty(do_stats)
    do_stats = true;
end
do_stats = logical(do_stats);
if nargin < 9 || isempty(normto)
    normto = 3; %default: row lines normalized to the centre column
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
if nargin < 7 || isempty(save_png)
    save_png = false;
end

%% each animal's top-ranked shank for the alignment label -> column of each row (3 = top shank)
topshank = NaN(1,numel(animals));
for a = 1:numel(animals)
    s1 = unique(T.Chan_Shank(T.Animal_Name == animals(a) & T.(rankfield) == 1));
    s1 = s1(~isnan(s1));
    if isempty(s1)
        warning('plotLWTTXCompare_multianimal: %s has no shank with %s == 1 (%s) - left out', animals(a), rankfield, region)
    else
        topshank(a) = min(s1);
    end
end
keepanimal = ~isnan(topshank);
animals = animals(keepanimal);
topshank = topshank(keepanimal);
assert(~isempty(animals), 'None of the chosen animals has a shank with %s == 1', rankfield)
[~,aidx] = ismember(T.Animal_Name, animals);
T = T(aidx > 0, :);
aidx = aidx(aidx > 0);
T.col = T.Chan_Shank - topshank(aidx)' + 3;
T = T(T.col >= 1 & T.col <= 5 & T.Chan_Depth >= 1 & T.Chan_Depth <= maxdepth, :);
na = numel(animals);

%% save_png: choose the file-name prefix before anything is plotted
if save_png
    if ~isfolder(savebase)
        mkdir(savebase);
    end
    answer = inputdlg('File-name prefix:', 'Save PNGs', [1 60], {sprintf('n%d', na)});
    if isempty(answer)
        error('plotLWTTXCompare_multianimal: stopped by the user (no file-name prefix)')
    end
    prefix = strtrim(answer{1});
    %if any file exists: Overwrite (same-name files are replaced), Quit, or Edit Prefix (checked again)
    while true
        pngnames = fullfile(savebase, prefix + "_" + region + "_LvsLW_TTX_" + features + ".png");
        existing = pngnames(isfile(pngnames));
        if isempty(existing)
            break
        end
        choice = questdlg(sprintf('%d of the PNG files already exist, e.g.\n%s', numel(existing), existing(1)), ...
            'Files already exist', 'Overwrite', 'Quit', 'Edit Prefix', 'Quit');
        switch choice
            case 'Overwrite'
                break
            case 'Edit Prefix'
                answer = inputdlg('New file-name prefix:', 'Edit prefix', [1 60], {prefix});
                if ~isempty(answer) && ~isempty(strtrim(answer{1}))
                    prefix = strtrim(answer{1}); %checked again by the loop
                end
            otherwise %Quit, or the dialog was closed
                error('plotLWTTXCompare_multianimal: stopped by the user (PNG files with prefix %s already exist)', prefix)
        end
    end
end

%% row 1: activity labels from each animal's ProbeInfo, aligned and normalized to the centre shank
labres = struct('name',{},'scores',{},'chanmaps',{},'heat',{},'n',{});
for L = 1:3
    labres(L).scores = NaN(na,5);
    labres(L).chanmaps = NaN(na,maxdepth,5);
end
labnames = strings(na,3); %label name used for each animal (the Spon label can differ between animals)
for a = 1:na
    pf = fullfile(probeinfo_dir, animals(a) + "-ProbeInfo.mat");
    if ~isfile(pf)
        warning('plotLWTTXCompare_multianimal: %s not found - no activity labels for %s', pf, animals(a))
        continue
    end
    pinfo = load(pf,'ProbeInfo');
    PI = pinfo.ProbeInfo;
    prb = find(strcmp(string(PI.Areas), region));
    if isempty(prb) %the Region value can be a probe number
        pnum = str2double(region);
        if ~isnan(pnum) && pnum == round(pnum) && pnum >= 1 && pnum <= numel(PI.ProbeMaps)
            prb = pnum;
        end
    end
    if numel(prb) ~= 1
        warning('plotLWTTXCompare_multianimal: %s has no probe for region %s - no activity labels', animals(a), region)
        continue
    end
    pmap = PI.ProbeMaps{prb};
    chremove = [];
    if isfield(PI,'Ch_Remove') && numel(PI.Ch_Remove) >= prb
        chremove = PI.Ch_Remove{prb};
    end
    %kept channels numbered 1..n from the top of each shank (as Chan_Depth in the summary script)
    depthmap = NaN(size(pmap));
    for s = 1:size(pmap,2)
        kept = pmap(:,s) > 0 & ~ismember(pmap(:,s), chremove);
        depthmap(kept,s) = 1:nnz(kept);
    end
    SL = getLabels(PI,'ShankLabels',prb);
    CL = getLabels(PI,'ChanLabels',prb);
    CLold = getLabels(PI,'chLabels',prb);
    sl = string(fieldnames(SL))';
    sponlab = sl(find(contains(sl,'Spon','IgnoreCase',true),1));
    if isempty(sponlab), sponlab = "Spon (none)"; end
    labnames(a,:) = [sponlab stimlabels];
    for L = 1:3
        lab = char(labnames(a,L));
        if ~isfield(SL,lab)
            warning('plotLWTTXCompare_multianimal: %s has no ShankLabels score %s - left out of that label', animals(a), lab)
            continue
        end
        sc = double(SL.(lab)(:)');
        if L == 3
            %row 1 column 3 (WhiskerOnly by default): scores can be negative or near 0, so scale by the
            %animal's largest |shank score| (whole probe): largest value -> +1, same value negative -> -1
            c0 = max(abs(sc),[],'omitnan');
        else
            c0 = sc(topshank(a)); %centre-shank score: the normalization reference
        end
        for c = 1:5
            s = topshank(a) + offsets(c);
            if s >= 1 && s <= numel(sc)
                labres(L).scores(a,c) = sc(s) / c0;
            end
        end
        chm = [];
        if isfield(CL,lab)
            chm = double(CL.(lab));
        elseif isfield(CLold,lab)
            chm = double(CLold.(lab));
        end
        if isempty(chm) || ~isequal(size(chm),size(pmap))
            continue %no channel scores: line only
        end
        for c = 1:5
            s = topshank(a) + offsets(c);
            if s < 1 || s > size(pmap,2), continue, end
            for d = 1:maxdepth
                r = find(depthmap(:,s) == d, 1);
                if ~isempty(r)
                    labres(L).chanmaps(a,d,c) = chm(r,s) / c0;
                end
            end
        end
    end
end
for L = 1:3
    names = unique(labnames(labnames(:,L) ~= "",L),'stable');
    if isempty(names), names = "(no label)"; end
    labres(L).name = strjoin(names,' / ');
    labres(L).heat = reshape(mean(labres(L).chanmaps,1,'omitnan'), maxdepth, 5);
    hasval = ~isnan(labres(L).scores) | reshape(any(~isnan(labres(L).chanmaps),2), na, 5); %animals x 5
    labres(L).n = sum(hasval, 1);
end

%colour maps (plotProbeLabels.m / plotLWIntensityCompare.m); one line colour per animal
orangemap = interp1([0 0.5 1],[1.00 0.96 0.90; 1.00 0.62 0.15; 0.92 0.38 0.00],linspace(0,1,256));
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

    %single conditions: per-animal means (animals x depth x 5), heatmaps, n, row lines
    conds = [cL cLW cLt cLWt];
    cnd = struct('name',cellstr(conds));
    allG = [];  devs = [];
    for i = 1:numel(conds)
        M = NaN(na,maxdepth,5);
        for a = 1:na
            M(a,:,:) = positionMap(T, feat, animals(a), conds(i), maxdepth);
        end
        cnd(i).name = conds(i);
        cnd(i).map = M;
        cnd(i).heat = reshape(mean(M,1,'omitnan'), maxdepth, 5);
        cnd(i).n = reshape(sum(any(~isnan(M),2),1), 1, 5);
        allG = [allG; cnd(i).heat(:)]; %#ok<AGROW>
        %row lines: line = mean over animals, points = each animal, both / the line's reference
        if shownorm
            [cnd(i).row, cnd(i).pts] = rowLines(M, cnd(i).heat, normto);
            devs = [devs; cnd(i).row(:); cnd(i).pts(:)]; %#ok<AGROW>
        end
    end

    %ratio tiles B / A: per animal and cell, animal's mean of B / animal's mean of A (1 value per animal);
    %heatmap = mean of the animal ratios. Order: LW/L, LW_TTX/L_TTX, L_TTX/L, LW_TTX/LW
    rspec = {cLW, cL, 1;  cLWt, cLt, 1;  cLt, cL, 2;  cLWt, cLW, 2};
    rat = struct('condB',rspec(:,1)','condA',rspec(:,2)','scale',rspec(:,3)');
    for j = 1:numel(rat)
        iA = conds == rat(j).condA;  iB = conds == rat(j).condB;
        mA = cnd(iA).map;  mA(mA == 0) = NaN; %no ratio where the denominator is 0
        R = cnd(iB).map ./ mA;
        rat(j).name = rat(j).condB + " / " + rat(j).condA;
        rat(j).map = R;
        rat(j).heat = reshape(mean(R,1,'omitnan'), maxdepth, 5);
        rat(j).n = reshape(sum(any(~isnan(R),2),1), 1, 5);
        %statistics: per cell, rank-sum of the B vs A trial values (animals pooled), FDR per tile
        rat(j).p = [];  rat(j).q = [];  rat(j).sig = [];
        if do_stats
            P = NaN(maxdepth,5);
            for r = 1:maxdepth
                for c = 1:5
                    atpos = T.Chan_Depth == r & T.col == c;
                    xa = T.(feat)(atpos & T.Condition_Name == rat(j).condA);  xa = xa(~isnan(xa));
                    xb = T.(feat)(atpos & T.Condition_Name == rat(j).condB);  xb = xb(~isnan(xb));
                    if numel(xa) >= min_trials && numel(xb) >= min_trials
                        P(r,c) = ranksum(xa, xb);
                    end
                end
            end
            if use_fdr
                Q = bhAdjust(P);
            else
                Q = P; %#ok<UNRCH> (use_fdr is a setting at the top)
            end
            rat(j).p = P;  rat(j).q = Q;  rat(j).sig = Q < alpha;
        end
    end

    %ratio comparisons: per animal, difference of two of its ratio maps; heatmap = mean over animals.
    %(1) (LW_TTX/LW) - (L_TTX/L)   (2) (LW_TTX/L_TTX) - (LW/L)
    cspec = {4, 3;  2, 1}; %{ratio index of the first term, ratio index of the subtracted term}
    cmp = struct('first',cspec(:,1)','second',cspec(:,2)');
    for j = 1:numel(cmp)
        Dr = rat(cmp(j).first).map - rat(cmp(j).second).map; %animals x depth x 5
        cmp(j).name = "(" + rat(cmp(j).first).name + ") - (" + rat(cmp(j).second).name + ")";
        cmp(j).map = Dr;
        cmp(j).heat = reshape(mean(Dr,1,'omitnan'), maxdepth, 5);
        cmp(j).n = reshape(sum(any(~isnan(Dr),2),1), 1, 5);
        %statistics: per cell, one-sample t-test of the animals' differences against 0, FDR per tile
        cmp(j).p = [];  cmp(j).q = [];  cmp(j).sig = [];
        if do_stats
            P = NaN(maxdepth,5);
            for r = 1:maxdepth
                for c = 1:5
                    x = Dr(:,r,c);  x = x(~isnan(x));
                    if numel(x) >= min_animals
                        [~,P(r,c)] = ttest(x);
                    end
                end
            end
            if use_fdr
                Q = bhAdjust(P);
            else
                Q = P; %#ok<UNRCH> (use_fdr is a setting at the top)
            end
            cmp(j).p = P;  cmp(j).q = Q;  cmp(j).sig = Q < alpha;
        end
    end

    %colour scales: green for the 4 single conditions; ratios around 1 (LW/L tiles, TTX-effect tiles);
    %ratio comparisons around 0
    climG = [min(allG,[],'omitnan') max(allG,[],'omitnan')];
    if isempty(climG) || any(isnan(climG)), climG = [0 1]; end
    if climG(1) == climG(2), climG = climG(1) + [-1 1]; end
    rmax = ones(1,2); %largest |ratio - 1| per ratio scale
    for g = 1:2
        h = [rat([rat.scale] == g).heat];
        m = max(abs(h(:) - 1),[],'omitnan');
        if ~isempty(m) && ~isnan(m) && m > 0, rmax(g) = m; end
    end
    h = [cmp.heat];
    cmax = max(abs(h(:)),[],'omitnan');
    if isempty(cmax) || isnan(cmax) || cmax == 0, cmax = 1; end
    k = 0; %row-line scale shared by the 4 single-condition tiles
    if shownorm
        maxdev = max(abs(devs - 1),[],'omitnan');
        if ~isempty(maxdev) && ~isnan(maxdev) && maxdev > 0
            k = linefill / maxdev;
        end
    end

    %figure
    ttl = {sprintf('%s (%s): L vs LW at 8 uW, before and with TTX, %d animals', feat, region, na), ...
        sprintf('aligned on each animal''s top %s shank (column "top"), Chan_Depth 1-%d', alignlabel, maxdepth)};
    if ~isempty(filtertxt)
        ttl{end+1} = filtertxt; %#ok<AGROW>
    end
    if shownorm && k > 0
        ttl{end+1} = sprintf('single-condition row lines: mean over animals, %s; o = animals; half a row = %.0f%% change', normtxt, 100*0.5/k); %#ok<AGROW>
    end
    ttl{end+1} = 'ratios B / A: mean over animals of each animal''s mean B / mean A (white = 1)'; %#ok<AGROW>
    ttl{end+1} = 'ratio comparisons: mean over animals of the difference of two animal ratios (white = 0)'; %#ok<AGROW>
    if do_stats
        ttl{end+1} = sprintf('stats: ratio tiles rank-sum B vs A trials per cell (animals pooled); comparisons one-sample t-test of the %d animal values vs 0', na); %#ok<AGROW>
        ttl{end+1} = sprintf('number = %s; * < 0.05, ** < 0.01, *** < 0.001', ternary(use_fdr,'BH-FDR adjusted p per tile','p')); %#ok<AGROW>
    end
    figs(end+1) = figure('Name',sprintf('%s (%s): L vs LW with TTX, %d animals',feat,region,na),'Color','w', ...
        'Units','pixels','Position',[40 30 1150 1250]); %#ok<AGROW> (fits a ~1300 px high screen)
    tl = tiledlayout(figs(end),5,4,'TileSpacing','compact','Padding','compact');
    title(tl,ttl,'Interpreter','none');

    %row 1: activity labels (nested 4 x 1 layout per label: line plot + heatmap), columns 1-3
    for L = 1:3
        ntl = tiledlayout(tl,4,1,'TileSpacing','compact');
        ntl.Layout.Tile = L;
        title(ntl,labres(L).name,'Interpreter','none');
        if L == 1, cmap = orangemap; else, cmap = greenmap; end
        axS = nexttile(ntl,1);
        hold(axS,'on');
        for a = 1:na
            plot(axS,1:5,labres(L).scores(a,:),'-o','LineWidth',1.5,'MarkerSize',5,'Color',acol(a,:), ...
                'MarkerFaceColor',acol(a,:),'DisplayName',animals(a));
        end
        xlim(axS,[0.5 5.5]); xticks(axS,1:5); xticklabels(axS,{});
        if L == 3
            ylabel(axS,{'shank score','(/ max |score|)'});
        else
            ylabel(axS,{'shank score','(/ top shank)'});
        end
        grid(axS,'on'); axS.GridAlpha = 0.15; box(axS,'off');
        if L == 1
            lg = legend(axS,'Interpreter','none','FontSize',7,'Orientation','horizontal','NumColumns',min(na,3));
            lg.Layout.Tile = 'north'; %above the line plot, so it covers no data
        end
        axH = nexttile(ntl,2,[3 1]);
        drawHeat(axH, labres(L).heat, cmap, [], labres(L).n);
        cb = colorbar(axH);
        cb.Layout.Tile = 'east';
        if L == 3
            cb.Label.String = 'channel score (/ max |shank score|)';
        else
            cb.Label.String = 'channel score (/ top shank)';
        end
        ylabel(axH,'Chan\_Depth');
        linkaxes([axS axH],'x');
    end

    %rows 2-3, columns 1-2: single conditions (tile = 4*(row-1) + column)
    singletile = [5 6 9 10]; %L_8, LW_8, L_8_TTX, LW_8_TTX
    axC = gobjects(1,4);
    for i = 1:4
        axC(i) = drawHeat(nexttile(tl,singletile(i)), cnd(i).heat, greenmap, climG, cnd(i).n, cnd(i).name);
        if shownorm && k > 0
            drawRowLines(axC(i), cnd(i).row, cnd(i).pts, k);
        end
    end
    ylabel(axC(1),'Chan\_Depth');  ylabel(axC(3),'Chan\_Depth');
    %ratio tiles: row 2 col 3, row 3 col 3, row 4 cols 1-2
    rattile = [7 11 13 14];
    axR = gobjects(1,4);
    for j = 1:4
        if do_stats
            rttl = sprintf('%s (%d/%d *)', rat(j).name, nnz(rat(j).sig), nnz(~isnan(rat(j).q)));
        else
            rttl = char(rat(j).name);
        end
        g = rat(j).scale;
        axR(j) = drawHeat(nexttile(tl,rattile(j)), rat(j).heat, divmap, 1 + [-rmax(g) rmax(g)], rat(j).n, rttl);
        if do_stats && ~all(isnan(rat(j).heat),'all')
            drawStats(axR(j), rat(j).q, rat(j).heat, 1, rmax(g), alpha);
        end
    end
    ylabel(axR(3),'Chan\_Depth');
    %ratio comparisons: row 4 col 4 (compares row 4 cols 2 and 1), row 5 col 3 (compares col 3 rows 3 and 2)
    cmptile = [16 19];
    axK = gobjects(1,2);
    for j = 1:2
        if do_stats
            kttl = {char(cmp(j).name), sprintf('(%d/%d *)', nnz(cmp(j).sig), nnz(~isnan(cmp(j).q)))};
        else
            kttl = char(cmp(j).name);
        end
        axK(j) = drawHeat(nexttile(tl,cmptile(j)), cmp(j).heat, divmap, [-cmax cmax], cmp(j).n, kttl);
        if do_stats && ~all(isnan(cmp(j).heat),'all')
            drawStats(axK(j), cmp(j).q, cmp(j).heat, 0, cmax, alpha);
        end
    end
    xlabel(axR(3),'Shank relative to the top-ranked shank');
    %colourbars in empty tiles (a colourbar next to a heatmap would shrink every tile of the layout)
    cbarTile(tl, 8, greenmap, climG, feat + " (mean of animal means)");
    cbarTile(tl, 12, divmap, 1 + [-rmax(1) rmax(1)], "LW / L ratio (mean of animal ratios)");
    cbarTile(tl, 15, divmap, 1 + [-rmax(2) rmax(2)], "TTX / before ratio (mean of animal ratios)");
    cbarTile(tl, 20, divmap, [-cmax cmax], "difference of ratios (mean over animals)");

    %save_png: <savebase>\<prefix>_<area>_LvsLW_TTX_<feature>.png
    if save_png
        exportgraphics(figs(end), fullfile(savebase, prefix + "_" + region + "_LvsLW_TTX_" + feat + ".png"), 'Resolution', 150);
    end

    data(fi).feature = feat;
    data(fi).region = region;
    data(fi).alignlabel = alignlabel;
    data(fi).stimlabels = stimlabels;
    data(fi).trialfilter = trialfilter;
    data(fi).normto = normto;
    data(fi).animals = animals;
    data(fi).topshank = topshank;
    data(fi).labels = labres;
    data(fi).cond = cnd;
    data(fi).ratio = rat;
    data(fi).compare = cmp;
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

function [rowN, ptsN] = rowLines(map, heat, normto)
%row lines: each row of heat (mean over animals) / its reference (column normto, or its peak);
%points: each animal's values (map, animals x rows x 5) / the same reference, so a line is the mean of its points
if isnumeric(normto)
    ref = heat(:,normto);
else
    ref = max(heat,[],2,'omitnan');
end
rowN = heat ./ ref;
ptsN = map ./ reshape(ref,1,[]);
end

function Q = bhAdjust(P)
%Benjamini-Hochberg adjusted p-values (same size as P; NaN = not tested)
Q = NaN(size(P));
tested = find(~isnan(P));
m = numel(tested);
if m == 0, return, end
[ps,ord] = sort(P(tested));
q = ps .* m ./ (1:m)';
q = min(flipud(cummin(flipud(q))), 1); %monotone, capped at 1
Q(tested(ord)) = q;
end

function ax = drawHeat(ax, M, cmap, lims, ncol, ttl)
%5-column heatmap (rows = Chan_Depth), NaN cells black, n animals per column in the tick labels;
%an all-NaN map is black with "No Data"
[nrow,ncols] = size(M);
imagesc(ax,1:ncols,1:nrow,M,'AlphaData',~isnan(M));
ax.Color = 'k';
colormap(ax,cmap);
if ~isempty(lims), clim(ax,lims); end
xlim(ax,[0.5 ncols+0.5]); ylim(ax,[0.5 nrow+0.5]);
yticks(ax,1:nrow); xticks(ax,1:ncols);
offs = -2:2;
lab = cell(1,ncols);
for c = 1:ncols
    lab{c} = sprintf('%s\\newline n=%d', offsetLabel(offs(c)), ncol(c));
end
ax.TickLabelInterpreter = 'tex';
xticklabels(ax,lab);
if nargin >= 6
    title(ax,ttl,'Interpreter','none','FontWeight','normal');
end
if all(isnan(M),'all')
    text(ax,0.5,0.5,'No Data','Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','Color','w','FontSize',12,'FontWeight','bold');
end
end

function drawRowLines(ax, rowN, ptsN, k)
%per row: dotted reference line (value 1 at the row centre), each animal's point (o), mean line on top
hold(ax,'on');
[nrow,ncols] = size(rowN);
for r = 1:nrow
    if all(isnan(rowN(r,:))), continue, end
    plot(ax,[0.5 ncols+0.5],[r r],':','Color',[0.35 0.35 0.35],'LineWidth',0.5);
    y = r - (rowN(r,:) - 1)*k;
    plot(ax,1:ncols,y,'-','Color','w','LineWidth',3); %white underlay, visible on dark cells
    for a = 1:size(ptsN,1)
        plot(ax,1:ncols, r - (squeeze(ptsN(a,r,:))' - 1)*k, 'o','Color','k','MarkerSize',3.5, ...
            'MarkerFaceColor','w','LineWidth',0.7);
    end
    plot(ax,1:ncols,y,'-','Color','k','LineWidth',1.4);
end
end

function drawStats(ax, Q, heat, centre, dmax, alpha)
%adjusted p in the centre of each tested cell; stars above it if significant; white text on dark
%cells (|value - centre| > 60% of the colour range; centre = 1 for ratios, 0 for differences)
for r = 1:size(Q,1)
    for c = 1:size(Q,2)
        if isnan(Q(r,c)), continue, end
        if abs(heat(r,c) - centre) > 0.6*dmax, tc = 'w'; else, tc = 'k'; end
        if Q(r,c) < 0.001
            ptxt = '<0.001';
        else
            ptxt = sprintf('%.3f', Q(r,c));
        end
        text(ax,c,r+0.12,ptxt,'HorizontalAlignment','center','VerticalAlignment','middle','Color',tc,'FontSize',8);
        if Q(r,c) < alpha
            if Q(r,c) < 0.001, st = '***'; elseif Q(r,c) < 0.01, st = '**'; else, st = '*'; end
            text(ax,c,r-0.2,st,'HorizontalAlignment','center','VerticalAlignment','middle','Color',tc, ...
                'FontSize',12,'FontWeight','bold');
        end
    end
end
end

function cbarTile(tl, tile, cmap, lims, label)
%horizontal colourbar alone in an empty tile of the layout (hidden axes carry the colour map and limits)
ax = nexttile(tl, tile);
colormap(ax, cmap);
clim(ax, lims);
axis(ax, 'off');
cb = colorbar(ax, 'Location', 'north');
cb.Label.String = label;
cb.Label.Interpreter = 'none';
end

function S = getLabels(ProbeInfo, fieldname, prb)
%struct of labels for probe prb from ProbeInfo.(fieldname), or an empty struct if there is none
S = struct();
if isfield(ProbeInfo,fieldname) && iscell(ProbeInfo.(fieldname)) && numel(ProbeInfo.(fieldname)) >= prb ...
        && isstruct(ProbeInfo.(fieldname){prb})
    S = ProbeInfo.(fieldname){prb};
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
