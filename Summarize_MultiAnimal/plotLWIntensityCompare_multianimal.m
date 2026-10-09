% plotLWIntensityCompare_multianimal.m
%
% Description: Multi-animal version of plotLWIntensityCompare.m, built from the summary
%   table like plotRankedShankRows.m. Compares light-only (L) and light+whisker (LW)
%   responses across light intensities for several animals. Each animal's shanks are
%   re-centred on its top-ranked shank for an alignment label (ShankRank_<alignlabel> == 1,
%   default LightOnly): that shank is column 3, its neighbouring shanks columns 1-2 (left)
%   and 4-5 (right). Rows are Chan_Depth 1..maxdepth (position among the shank's kept
%   channels, 1 = top; maxdepth = 3). One figure per chosen feature, a 5 x 3 tiled layout (row 1 is 25% taller than rows 2-5):
%     row 1: activity labels from each animal's ProbeInfo, all aligned by the alignment label:
%       (1,1) the spontaneous label (name containing "Spon"), (1,2)-(1,3) stimlabels (default
%       LightOnly, WhiskerOnly). Each tile: a line per animal of the label's shank scores for
%       the 5 aligned shanks, and below it a 5 x maxdepth heatmap of the mean over animals of
%       the channel scores (ProbeInfo.ChanLabels, or chLabels). Columns 1-2: each animal's
%       shank and channel scores are divided by its score on the centre shank (lines meet at
%       1 in column 3). Column 3 (WhiskerOnly by default; its scores can be negative or near
%       0): divided by the animal's largest |shank score| over the whole probe, so its
%       largest shank is +1 or -1 and all shanks lie within [-1, 1] (single channels can go
%       slightly beyond). Orange for Spon, green for stimulus labels (as plotProbeLabels.m)
%     rows 2-5: light intensities 4, 8, 12, 15. Each tile is a 5 x maxdepth heatmap of the
%       mean over animals of each animal's per-position mean over trials:
%       column 1 = L_<n>, column 2 = LW_<n> (green, one colour scale shared by all 8 tiles);
%       with rowlines, each of these tiles also shows per row a line of the row's values
%       across the 5 columns, stretched so the row's min / max sit at 20% / 80% of the row
%       height (shape only; the colour gives the size). Column 3 = LW_<n> - L_<n> (each
%       animal's difference first, then the mean; blue-white-red, one symmetric scale for
%       the 4 tiles)
%   Statistics (do_stats, column 3): per cell, the feature values of all chosen animals'
%   trials at that position (animals pooled): Wilcoxon rank-sum test (ranksum) of LW_<n> vs
%   L_<n> (unpaired, NaN ignored, at least min_trials per condition), Benjamini-Hochberg FDR
%   over the tested cells of that tile (use_fdr). The adjusted p (q; raw p if use_fdr is
%   false) is written in the centre of each tested cell; significant cells get stars above
%   it (* q < 0.05, ** < 0.01, *** < 0.001). Tile title: "<significant>/<tested> *".
%   A condition without rows for any chosen animal gives an all-black "No Data" tile (and
%   "No Data" in that intensity's LW - L tile); animals missing a condition simply do not
%   contribute to it. Every heatmap gives n = animals with data per column under its
%   columns. Optional trial filter as in plotRankedShankRows.m. alpha, use_fdr, min_trials,
%   maxdepth and the paths are set at the top of the function. The figure title is short,
%   "<region> Response to L vs LW: <feature> (n = <animals>)"; the run information (table,
%   alignment, animals and top shanks, filter, row lines, statistics, colour scales) is in
%   data.info and, with save_png, in a .txt file next to the PNG. Returns figure handles and
%   the plotted numbers; nothing is printed.
%
% Inputs (same order as plotLWIntensityCompare.m; alignlabel, rowlines and trialfilter added at the end):
%   animals (string array, optional) - Animal_Name values; omitted/empty = list dialog (multiple)
%   summaryTable (table, or char/string path to its csv) - tableres from
%     SummaryAnalysis_MultiMouse_allChannels.m, one row per good trial x channel. Columns used:
%       Animal_Name, Condition_Name, Region (text)
%       Chan_Shank (double) - shank (ProbeMaps column) of the channel
%       Chan_Depth (double) - position among that shank's kept channels, 1 = top; NaN rows dropped
%       ShankRank_<alignlabel> (double) - rank of the channel's shank by that ShankLabels score
%       measure columns named <P1|P2|All>_... (double) - the features offered
%     Conditions plotted: L_4, LW_4, L_8, LW_8, L_12, LW_12, L_15, LW_15 (error only if none
%     of them has rows for the chosen animals and region)
%     INFERRED DATA CONTRACT: neighbouring Chan_Shank numbers are physically neighbouring
%     shanks (ProbeMaps column order = physical shank order)
%   region (char/string, optional) - Region value; omitted/empty = the only region of the chosen
%     animals, or a list dialog. Also picks each animal's probe (matched in ProbeInfo.Areas,
%     or used as the probe number)
%   features (string array, optional) - measure columns, one figure each; omitted/empty =
%     list dialog (multiple) of the measure columns with data
%   stimlabels (string array, 1 x 2, optional) - ShankLabels for row 1 columns 2-3; default
%     ["LightOnly","WhiskerOnly"]. Row 1 column 1 is the first label containing "Spon"
%   do_stats (logical, optional) - true (default): LW vs L test per cell in column 3
%   save_png (logical, optional) - true: also save each figure as a PNG (see Outputs);
%     default false. Once features and region are known and before anything is plotted, a
%     text box shows the whole default file name, n<number of animals>_<area>_LvsLW, to edit
%     (the feature is added to it; Cancel stops the function). If any of the PNG / TXT files
%     already exists, a dialog offers Overwrite, Quit (stops with an error) or Edit Name (the
%     new name is checked again). Closing it = Quit
%   alignlabel (char/string, optional) - label whose top-ranked shank is column 3 (a
%     ShankRank_<label> column, prefix optional); default "LightOnly"
%   rowlines (logical, optional) - true (default): draw the stretched row lines on the L / LW
%     heatmaps; false: plain heatmaps
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
%     .feature, .region, .alignlabel, .stimlabels, .trialfilter, .rowlines
%     .info (string array) - the run information (written to the .txt file with save_png)
%     .animals (string, 1 x nAnimals), .topshank (double, 1 x nAnimals) - Chan_Shank of each
%       animal's top-ranked shank for alignlabel
%     .labels (struct array, 1 x 3) - row 1: .name, .scores (nAnimals x 5, normalized shank
%       scores: / centre shank for labels 1-2, / max |shank score| for label 3), .chanmaps
%       (nAnimals x maxdepth x 5, channel scores normalized the same way), .heat
%       (maxdepth x 5), .n (1 x 5)
%     .intensity (struct array, 1 x 4) - rows 2-5: .n_uW (4/8/12/15), .mapL, .mapLW (nAnimals x
%       maxdepth x 5), .heatL, .heatLW, .heatD (maxdepth x 5), .nL, .nLW, .nD (1 x 5),
%       .p, .q, .sig (maxdepth x 5; empty without statistics)
%   with save_png: <savebase>\<name>_<feature>.png (150 dpi) and <name>_<feature>.txt (data.info),
%     one pair per feature (savebase = E:\Roy\Processed Silicon Probe Data\BundledAnimalData\MultAnimal)
%
% Dependencies: SummaryAnalysis_MultiMouse_allChannels.m for summaryTable, including its
%   separate section that adds Chan_Shank / Chan_Depth / ShankRank_<score>; ProbeInfo labels
%   from ChannelShankLabelling_SponData.m / ChannelShank_activitylabelling.m. Single-animal
%   version: plotLWIntensityCompare.m; heatmaps and statistics as plotRankedShankRows.m.
%   MATLAB R2020b+ (nested tiledlayout); Statistics and Machine Learning Toolbox (ranksum).

function [figs, data] = plotLWIntensityCompare_multianimal(animals, summaryTable, region, features, ...
    stimlabels, do_stats, save_png, alignlabel, rowlines, trialfilter)

figs = gobjects(1,0);
data = struct([]);
intensities = [4 8 12 15];
maxdepth = 3;       %rows: Chan_Depth 1..maxdepth
alpha = 0.05;       %significance level
use_fdr = true;     %Benjamini-Hochberg FDR correction across the tested cells of each LW - L tile
min_trials = 3;     %minimum trials (non-NaN) per condition to test a cell
offsets = -2:2;     %columns: shanks relative to the top-ranked shank
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';
savebase = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\MultAnimal';

%% summary table: from a csv path or given directly; check the columns needed
if ischar(summaryTable) || isstring(summaryTable)
    tablesource = string(summaryTable);
    summaryTable = readtable(char(summaryTable), 'TextType','string');
else
    tablesource = "table passed in (tableres)";
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
condL = "L_" + intensities;
condLW = "LW_" + intensities;
assert(any(ismember([condL condLW], T.Condition_Name)), 'The chosen animals have no rows for any of %s (%s)', ...
    strjoin([condL condLW],', '), region)
T = T(ismember(T.Condition_Name, [condL condLW]), :);

%% alignment label
if nargin < 8 || isempty(alignlabel)
    alignlabel = "LightOnly"; %default: centre on each animal's top LightOnly shank
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
if nargin < 9 || isempty(rowlines)
    rowlines = true; %default: draw the row lines on the L / LW heatmaps
end
rowlines = logical(rowlines);
if nargin < 7 || isempty(save_png)
    save_png = false;
end

%% each animal's top-ranked shank for the alignment label -> column of each row (3 = top shank)
topshank = NaN(1,numel(animals));
for a = 1:numel(animals)
    s1 = unique(T.Chan_Shank(T.Animal_Name == animals(a) & T.(rankfield) == 1));
    s1 = s1(~isnan(s1));
    if isempty(s1)
        warning('plotLWIntensityCompare_multianimal: %s has no shank with %s == 1 (%s) - left out', animals(a), rankfield, region)
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

%% save_png: choose the file name before anything is plotted (the whole default name is shown and editable)
if save_png
    if ~isfolder(savebase)
        mkdir(savebase);
    end
    savename = sprintf('n%d_%s_LvsLW', na, region); %default: n<animals>_<area>_LvsLW
    namehelp = sprintf(['File name (the feature is added): <name>_<feature>.png and .txt\n' ...
        'Folder: %s\nWith this name: %s_%s.png'], savebase, savename, features(1));
    answer = inputdlg(namehelp, 'Save PNGs', [1 80], {savename});
    if isempty(answer) || isempty(strtrim(answer{1}))
        error('plotLWIntensityCompare_multianimal: stopped by the user (no file name)')
    end
    savename = strtrim(answer{1});
    %if any file exists: Overwrite (same-name files are replaced), Quit, or Edit Name (checked again)
    while true
        outnames = fullfile(savebase, savename + "_" + features + [".png"; ".txt"]);
        existing = outnames(isfile(outnames));
        if isempty(existing)
            break
        end
        choice = questdlg(sprintf('%d of the output files already exist, e.g.\n%s', numel(existing), existing(1)), ...
            'Files already exist', 'Overwrite', 'Quit', 'Edit Name', 'Quit');
        switch choice
            case 'Overwrite'
                break
            case 'Edit Name'
                answer = inputdlg(sprintf('New file name (the feature is added):\n<name>_<feature>.png and .txt in %s', savebase), ...
                    'Edit name', [1 80], {savename});
                if ~isempty(answer) && ~isempty(strtrim(answer{1}))
                    savename = strtrim(answer{1}); %checked again by the loop
                end
            otherwise %Quit, or the dialog was closed
                error('plotLWIntensityCompare_multianimal: stopped by the user (files named %s_<feature> already exist)', savename)
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
        warning('plotLWIntensityCompare_multianimal: %s not found - no activity labels for %s', pf, animals(a))
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
        warning('plotLWIntensityCompare_multianimal: %s has no probe for region %s - no activity labels', animals(a), region)
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
            warning('plotLWIntensityCompare_multianimal: %s has no ShankLabels score %s - left out of that label', animals(a), lab)
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

    %rows 2-5: per-animal means (animals x depth x 5), heatmaps, n, row lines, statistics
    inten = struct('n_uW',num2cell(intensities));
    allLLW = [];  allD = [];
    for n = 1:numel(intensities)
        mapL = NaN(na,maxdepth,5);  mapLW = NaN(na,maxdepth,5);
        for a = 1:na
            mapL(a,:,:) = positionMap(T, feat, animals(a), condL(n), maxdepth);
            mapLW(a,:,:) = positionMap(T, feat, animals(a), condLW(n), maxdepth);
        end
        D = mapLW - mapL;
        inten(n).mapL = mapL;  inten(n).mapLW = mapLW;
        inten(n).heatL = reshape(mean(mapL,1,'omitnan'), maxdepth, 5);
        inten(n).heatLW = reshape(mean(mapLW,1,'omitnan'), maxdepth, 5);
        inten(n).heatD = reshape(mean(D,1,'omitnan'), maxdepth, 5);
        inten(n).nL = reshape(sum(any(~isnan(mapL),2),1), 1, 5);
        inten(n).nLW = reshape(sum(any(~isnan(mapLW),2),1), 1, 5);
        inten(n).nD = reshape(sum(any(~isnan(D),2),1), 1, 5);
        allLLW = [allLLW; inten(n).heatL(:); inten(n).heatLW(:)];
        allD = [allD; inten(n).heatD(:)];
        %statistics (LW - L tile): per cell, rank-sum of the LW vs L trial values (animals pooled), FDR per tile
        inten(n).p = [];  inten(n).q = [];  inten(n).sig = [];
        if do_stats
            P = NaN(maxdepth,5);
            for r = 1:maxdepth
                for c = 1:5
                    atpos = T.Chan_Depth == r & T.col == c;
                    xa = T.(feat)(atpos & T.Condition_Name == condL(n));   xa = xa(~isnan(xa));
                    xb = T.(feat)(atpos & T.Condition_Name == condLW(n));  xb = xb(~isnan(xb));
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
            inten(n).p = P;  inten(n).q = Q;  inten(n).sig = Q < alpha;
        end
    end
    climG = [min(allLLW,[],'omitnan') max(allLLW,[],'omitnan')];
    if isempty(climG) || any(isnan(climG)), climG = [0 1]; end
    if climG(1) == climG(2), climG = climG(1) + [-1 1]; end
    dmax = max(abs(allD),[],'omitnan');
    if isempty(dmax) || isnan(dmax) || dmax == 0, dmax = 1; end
    climD = [-dmax dmax];

    %short figure title; the run information goes to data.info (and to the .txt file with save_png)
    ttl = sprintf('%s Response to L vs LW: %s (n = %d)', region, feat, na);
    info = string(ttl); %#ok<*AGROW> (info grows line by line, a dozen lines)
    info(end+1) = "created " + string(datetime('now','Format','yyyy-MM-dd HH:mm')) + " by plotLWIntensityCompare_multianimal";
    info(end+1) = "summary table: " + tablesource;
    info(end+1) = sprintf('region %s; light intensities %s uW (L_n, LW_n)', region, strjoin(string(intensities),', '));
    info(end+1) = sprintf('aligned on each animal''s top %s shank (column "top"); rows = Chan_Depth 1-%d', alignlabel, maxdepth);
    info(end+1) = "animals (top shank): " + strjoin(animals + " (" + string(topshank) + ")", ', ');
    if isempty(filtertxt)
        info(end+1) = "trials: no trial filter (" + nafter + " table rows)";
    else
        info(end+1) = filtertxt;
    end
    info(end+1) = "heatmaps: mean over animals of each animal's per-position mean over trials; n = animals with data per column";
    info(end+1) = "colour scales: L and LW tiles share one green scale; LW - L tiles share one symmetric blue-white-red scale (white = 0)";
    info(end+1) = "row 1: shank-score lines per animal and mean channel scores, / centre-shank score (" + labres(1).name + ", " + ...
        labres(2).name + "), / largest |shank score| (" + labres(3).name + ")";
    if rowlines
        info(end+1) = "L / LW row lines: per row, the mean over animals across the 5 columns, stretched so the row's min / max sit at 20% / 80% of the row height (shape only; colour gives the size)";
    else
        info(end+1) = "L / LW row lines: off";
    end
    if do_stats
        info(end+1) = sprintf('LW - L statistics: per cell, rank-sum test of LW vs L trial values (animals pooled, at least %d trials each); number = %s; * < 0.05, ** < 0.01, *** < 0.001', ...
            min_trials, ternary(use_fdr,'BH-FDR adjusted p over the cells of each tile','p'));
    else
        info(end+1) = "LW - L statistics: off";
    end
    figs(end+1) = figure('Name',sprintf('%s (%s): L vs LW, %d animals',feat,region,na),'Color','w', ...
        'Units','pixels','Position',[40 30 1250 1300]);
    %21 x 3 grid: row 1 spans 5 grid rows, each intensity row 4 (row 1 = 1.25 x the others)
    tl = tiledlayout(figs(end),21,3,'TileSpacing','compact','Padding','compact');
    title(tl,ttl,'Interpreter','none');

    %row 1: activity labels (nested 5 x 1 layout per label: line plot 2/5 + heatmap 3/5)
    for L = 1:3
        ntl = tiledlayout(tl,5,1,'TileSpacing','compact');
        ntl.Layout.Tile = L;
        ntl.Layout.TileSpan = [5 1];
        title(ntl,labres(L).name,'Interpreter','none');
        if L == 1, cmap = orangemap; else, cmap = greenmap; end
        axS = nexttile(ntl,1,[2 1]);
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
        axH = nexttile(ntl,3,[3 1]);
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

    %rows 2-5: L, LW and LW - L per light intensity
    for n = 1:numel(intensities)
        base = 3*(5 + 4*(n-1)); %tile before this row (row n+1 starts at grid row 6 + 4*(n-1))
        axL = drawHeat(nexttile(tl,base+1,[4 1]), inten(n).heatL, greenmap, climG, inten(n).nL, condL(n));
        axLW = drawHeat(nexttile(tl,base+2,[4 1]), inten(n).heatLW, greenmap, climG, inten(n).nLW, condLW(n));
        if rowlines
            drawRowLines(axL, inten(n).heatL);
            drawRowLines(axLW, inten(n).heatLW);
        end
        if do_stats
            dttl = sprintf('%s - %s (%d/%d *)', condLW(n), condL(n), nnz(inten(n).sig), nnz(~isnan(inten(n).q)));
        else
            dttl = sprintf('%s - %s', condLW(n), condL(n));
        end
        axD = drawHeat(nexttile(tl,base+3,[4 1]), inten(n).heatD, divmap, climD, inten(n).nD, dttl);
        if do_stats && ~all(isnan(inten(n).heatD),'all')
            drawStats(axD, inten(n).q, inten(n).heatD, dmax, alpha);
        end
        ylabel(axL,'Chan\_Depth');
        if n == numel(intensities)
            cb = colorbar(axLW,'southoutside'); cb.Label.String = feat + " (L, LW; mean of animal means)"; cb.Label.Interpreter = 'none';
            cb = colorbar(axD,'southoutside');  cb.Label.String = feat + " (LW - L)"; cb.Label.Interpreter = 'none';
            xlabel(axL,'Shank relative to the top-ranked shank');
        end
    end

    %save_png: <savebase>\<name>_<feature>.png and .txt (the run information)
    if save_png
        exportgraphics(figs(end), fullfile(savebase, savename + "_" + feat + ".png"), 'Resolution', 150);
        fid = fopen(fullfile(savebase, savename + "_" + feat + ".txt"), 'w');
        fprintf(fid, '%s\n', info);
        fclose(fid);
    end

    data(fi).feature = feat;
    data(fi).region = region;
    data(fi).alignlabel = alignlabel;
    data(fi).stimlabels = stimlabels;
    data(fi).trialfilter = trialfilter;
    data(fi).rowlines = rowlines;
    data(fi).info = info;
    data(fi).animals = animals;
    data(fi).topshank = topshank;
    data(fi).labels = labres;
    data(fi).intensity = inten;
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

function drawRowLines(ax, heat)
%per heatmap row: the row's values across the columns as a line, stretched so its min / max sit at
%20% / 80% of the row height (higher values up; a flat row at the row centre). Shape only.
hold(ax,'on');
[nrow,ncols] = size(heat);
for r = 1:nrow
    v = heat(r,:);
    if all(isnan(v)), continue, end
    lo = min(v,[],'omitnan');  hi = max(v,[],'omitnan');
    if hi > lo
        y = r + 0.3 - 0.6*(v - lo)/(hi - lo);
    else
        y = r + zeros(size(v));
    end
    plot(ax,1:ncols,y,'-','Color','w','LineWidth',3); %white underlay, visible on dark cells
    plot(ax,1:ncols,y,'-o','Color','k','LineWidth',1.4,'MarkerSize',3,'MarkerFaceColor','k');
end
end

function drawStats(ax, Q, heat, dmax, alpha)
%adjusted p in the centre of each tested cell; stars above it if significant; white text on dark cells
for r = 1:size(Q,1)
    for c = 1:size(Q,2)
        if isnan(Q(r,c)), continue, end
        if abs(heat(r,c)) > 0.6*dmax, tc = 'w'; else, tc = 'k'; end
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
