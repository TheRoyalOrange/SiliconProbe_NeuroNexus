% plotShankDistance_ConditionCompare.m
%
% Description: Asks whether the distance between an animal's spontaneous-activity peak shank and
%   its peak shank for an alignment label (e.g. LightOnly) changes the difference between two
%   conditions (e.g. L_8 vs LW_8). For each animal (from its ProbeInfo.ShankLabels):
%     sponpeak  = shank with the highest spontaneous score (label containing "Spon", or sponlabel)
%     alignpeak = shank with the highest alignment-label score
%     shankdist = sponpeak - alignpeak  (positive = alignment peak is LEFT of, i.e. a lower shank
%                 number than, the spon peak; e.g. spon 6, LightOnly 4 -> +2). Negative allowed.
%   Peaks: NaN scores ignored, ties -> lower shank number (same rule as the ShankRank_ columns).
%   Data points: every trial of the two conditions on the channels at Chan_Depth 1..maxdepth of the
%   alignment peak shank (the centre "top" column of the plotLWIntensityCompare_multianimal.m
%   heatmaps). One figure, a maxdepth x nFeatures tiled layout (rows = Chan_Depth 1..maxdepth,
%   columns = features). x = shankdist (each integer in its range; tick label gives n animals with
%   points in that tile). At each x: box and whisker of condition 1 (black, left) and condition 2
%   (red, right), each pooling the trials of all animals at that distance (box = quartiles, line =
%   median, whiskers = most extreme points within 1.5 IQR), with all points overlaid (small
%   horizontal jitter, same colours, one marker shape per animal). y axes are linked down each
%   feature column.
%   Statistics (do_stats): per tile (Chan_Depth x feature) and per shank distance, condition 2 vs
%   condition 1 on the trials at that distance: linear mixed model y ~ Cond + (1|Animal) (fitlme,
%   random intercept per animal) when >= 2 animals share the distance, else y ~ Cond (fitlm, no
%   animal term possible, equivalent to a t-test). p = the Cond coefficient's p (Estimate =
%   cond2 - cond1); tested only with >= min_trials non-NaN points per condition. Benjamini-Hochberg
%   FDR over the tested distances of each tile (use_fdr). The adjusted p (q; raw p if use_fdr is
%   false) is written above each tested distance, with stars if significant (* q < 0.05,
%   ** < 0.01, *** < 0.001); tile title "<significant>/<tested> *". CAVEAT: the models assume
%   roughly normal residuals, and AUC values are often skewed. INFERRED: trials of one animal at
%   one channel are treated as independent apart from the animal intercept.
%   Problems (animal without ProbeInfo / probe / label, ProbeInfo peak differing from the table's
%   ShankRank_<align> == 1 shank, a failed model fit) go into runnotes, shown at the end with
%   animalinfo and stattable; nothing is printed while it runs. The figure title is short; the run
%   information (table, animals and their peaks / distances, conditions, features, box and
%   statistics definitions) is in info and, with save_png, in a .txt file next to the PNG.
%
% Inputs:
%   summaryTable (table, workspace variable, set by the user before running - this script does not
%     make it) - tableres from SummaryAnalysis_MultiMouse_allChannels.m (or readtable of its csv,
%     or several of them stacked), one row per good trial x channel. Columns used:
%       Animal_Name, Condition_Name, Region (text)
%       Chan_Shank (double) - shank (ProbeMaps column) of the channel
%       Chan_Depth (double) - position among that shank's kept channels, 1 = top; NaN rows dropped
%       measure columns named <P1|P2|All>_... (double, e.g. All_MUAAUC) - the features offered
%       optional: Channel_ID, Trial, Condition_FileNum (copied into ptdata if present);
%         ShankRank_<align> (double) - used only for the cross-check
%   Settings section: maxdepth (3), sponlabel ("" = first ShankLabels field containing "Spon",
%     as plotLWIntensityCompare_multianimal.m), probeinfo_dir, savebase, and
%     save_png (logical) - true: also save the figure as a PNG (see Outputs); default false. Once
%     the data points are gathered and before anything is plotted, a text box shows the whole
%     default file name, n<animals>_<area>_<cond1>vs<cond2>_<align>Dist_<features>, to edit
%     (Cancel stops the script). If the PNG or TXT already exists, a dialog offers Overwrite,
%     Quit (stops with an error) or Edit Name (the new name is checked again). Closing it = Quit
%     (as plotLWIntensityCompare_multianimal.m)
%     do_stats (logical, default true), alpha (0.05), use_fdr (true), min_trials (3) - statistics
%   Dialogs: animals (multiple), region (only if more than one), the alignment label (ShankLabels
%     fields shared by all chosen animals' ProbeInfo, LightOnly preselected), exactly 2 conditions
%     (1st = black, 2nd = red), features (multiple, one column of tiles each)
%   <probeinfo_dir>\<animal>-ProbeInfo.mat - each animal's ProbeInfo: .Areas (probe found by
%     region, or region used as the probe number), .ProbeMaps, .ShankLabels{probe}.<label>
%     (double, 1 x nShanks, ProbeMaps column order)
%   INFERRED DATA CONTRACTS (not stated elsewhere):
%     - lower Chan_Shank / ProbeMaps column number = physically more LEFT shank
%     - the Spon label is the first ShankLabels field containing "Spon" (unless sponlabel is set)
%     - Chan_Shank in the table uses the same shank numbering as ProbeInfo.ShankLabels
%
% Outputs (workspace):
%   animalinfo (table, nAnimals x ...) - Animal_Name, SponLabel, SponPeak, AlignPeak, ShankDist,
%     TableTopShank (Chan_Shank with ShankRank_<align> == 1, NaN if no such column), and
%     n_<condition>_d<depth> (double) - number of points per condition and depth
%   ptdata (table, one row per plotted point) - Animal_Name, ShankDist, Condition_Name, Chan_Shank,
%     Chan_Depth, Channel_ID / Trial / Condition_FileNum (if in summaryTable), and the features
%   stattable (table, one row per tile x distance with points; empty without do_stats) - Feature,
%     Chan_Depth, ShankDist, nAnimals, n_<cond1>, n_<cond2> (points), Model ("LME", "LM", or ""
%     = not tested / fit failed), Estimate (cond2 - cond1), p, q (BH-FDR adjusted, or p), Sig
%   runnotes (string array) - problems met during the run (empty = none)
%   info (string array) - the run information (written to the .txt file with save_png)
%   fig (figure handle)
%   with save_png: <savebase>\<name>.png (150 dpi) and <name>.txt (info)
%     (savebase = E:\Roy\Processed Silicon Probe Data\BundledAnimalData\MultAnimal)
%
% Dependencies: SummaryAnalysis_MultiMouse_allChannels.m for summaryTable (including its section
%   adding Chan_Shank / Chan_Depth); ProbeInfo.ShankLabels from ChannelShankLabelling_SponData.m /
%   ChannelShank_activitylabelling.m. Same alignment as plotLWIntensityCompare_multianimal.m.
%   MATLAB R2020b+ (tiledlayout legend tile); Statistics and Machine Learning Toolbox (quantile,
%   fitlme, fitlm).

%%
summaryTable1 = readtable("E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\LvsLW_variedIntensity_n4.csv");
summaryTable2 = readtable("E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\LvsLW_TTXtreatmentS1_n3.csv");

%summaryTable1 = summaryTable1(summaryTable1.Animal_Num ~= 2,:);

summaryTable = cat(1,summaryTable1,summaryTable2);
summaryTable.All_LFPRMSsum = summaryTable.P1_LFPRMSsum + summaryTable.P2_LFPRMSsum;

%% settings
maxdepth = 3;       %rows: Chan_Depth 1..maxdepth on the alignment peak shank
sponlabel = "";     %ShankLabels field for spontaneous activity; "" = first field containing "Spon"
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';
condcolors = [0 0 0; 0.85 0.1 0.1];          %condition 1 black, condition 2 red
animalmarkers = {'o','s','^','d','v','>','<','p','h'}; %one per animal (cycled)
boxoffset = 0.18;   %condition box centres at x -/+ boxoffset
boxwidth = 0.3;     %box width (shank units)
jitterwidth = 0.14; %points spread over +/- jitterwidth/2 around the box centre
save_png = false;   %true: save the figure as <savebase>\<name>.png + <name>.txt (name chosen in a dialog)
savebase = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\MultAnimal';
do_stats = true;    %cond2 vs cond1 per tile and shank distance (mixed model, see header)
alpha = 0.05;       %significance level
use_fdr = true;     %Benjamini-Hochberg FDR across the tested distances of each tile
min_trials = 3;     %minimum points (non-NaN) per condition to test a distance

%% summary table (made by the user, not here): check the columns needed
assert(exist('summaryTable','var') == 1 && istable(summaryTable), ...
    'Put the summary table (tableres, or readtable of its csv) in the workspace as summaryTable first')
vn = string(summaryTable.Properties.VariableNames);
needcols = ["Animal_Name","Condition_Name","Region","Chan_Shank","Chan_Depth"];
missingcols = needcols(~ismember(needcols, vn));
assert(isempty(missingcols), 'summaryTable is missing column(s): %s', strjoin(missingcols,', '))
T = summaryTable;
T.Animal_Name = string(T.Animal_Name);
T.Condition_Name = string(T.Condition_Name);
T.Region = string(T.Region);
T = T(~isnan(T.Chan_Depth), :);
runnotes = strings(0,1);

%% animals
allanimals = unique(T.Animal_Name,'stable');
[sel,ok] = listdlg('ListString',cellstr(allanimals),'SelectionMode','multiple', ...
    'InitialValue',1:numel(allanimals),'PromptString','Animals to include:','ListSize',[250 200]);
assert(ok && ~isempty(sel), 'No animal selected')
animals = allanimals(sel)';
T = T(ismember(T.Animal_Name, animals), :);

%% region
regions = unique(T.Region,'stable');
if numel(regions) == 1
    region = regions;
else
    [sel,ok] = listdlg('ListString',cellstr(regions),'SelectionMode','single', ...
        'PromptString','Which region?','ListSize',[250 120]);
    assert(ok && ~isempty(sel), 'No region selected')
    region = regions(sel);
end
T = T(T.Region == region, :);

%% each animal's ShankLabels for that region's probe (from ProbeInfo)
na = numel(animals);
SLs = cell(1,na);       %ShankLabels struct of the region's probe, [] if none
labelsets = cell(1,na); %names of the 1 x nShanks numeric fields
for a = 1:na
    pf = fullfile(probeinfo_dir, animals(a) + "-ProbeInfo.mat");
    if ~isfile(pf)
        runnotes(end+1) = sprintf('%s: %s not found - animal left out', animals(a), pf); %#ok<SAGROW>
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
        runnotes(end+1) = sprintf('%s: no single probe for region %s in ProbeInfo.Areas - animal left out', animals(a), region); %#ok<SAGROW>
        continue
    end
    if ~(isfield(PI,'ShankLabels') && iscell(PI.ShankLabels) && numel(PI.ShankLabels) >= prb && isstruct(PI.ShankLabels{prb}))
        runnotes(end+1) = sprintf('%s: no ShankLabels for probe %d - animal left out', animals(a), prb); %#ok<SAGROW>
        continue
    end
    SLs{a} = PI.ShankLabels{prb};
    nshank = size(PI.ProbeMaps{prb},2);
    fn = string(fieldnames(SLs{a}))';
    isscore = arrayfun(@(f) isnumeric(SLs{a}.(char(f))) && numel(SLs{a}.(char(f))) == nshank, fn);
    labelsets{a} = fn(isscore);
end
clear pinfo PI
haslabels = ~cellfun(@isempty, labelsets);
assert(any(haslabels), 'None of the chosen animals has ShankLabels for region %s', region)

%% alignment label (shared by all animals with labels; the spon label is not offered)
commonlabels = labelsets{find(haslabels,1)};
for a = find(haslabels)
    commonlabels = intersect(commonlabels, labelsets{a}, 'stable');
end
alignchoices = commonlabels(~contains(commonlabels,'Spon','IgnoreCase',true));
if sponlabel ~= ""
    alignchoices = setdiff(alignchoices, sponlabel, 'stable');
end
assert(~isempty(alignchoices), 'No ShankLabels score shared by all chosen animals to align on')
initsel = find(alignchoices == "LightOnly", 1);
if isempty(initsel), initsel = 1; end
[sel,ok] = listdlg('ListString',cellstr(alignchoices),'SelectionMode','single','InitialValue',initsel, ...
    'PromptString','Alignment label (peak shank vs spon peak):','ListSize',[300 200]);
assert(ok && ~isempty(sel), 'No alignment label selected')
alignlabel = alignchoices(sel);

%% two conditions to compare
condchoices = unique(T.Condition_Name,'stable');
sel = [];
while numel(sel) ~= 2
    [sel,ok] = listdlg('ListString',cellstr(condchoices),'SelectionMode','multiple', ...
        'PromptString',{'Choose exactly 2 conditions','(1st = black, 2nd = red):'},'ListSize',[250 250]);
    assert(ok, 'No conditions selected')
end
conds = condchoices(sel)';
T = T(ismember(T.Condition_Name, conds), :);

%% features: measure columns with data
measurecols = vn(~cellfun(@isempty, regexp(cellstr(vn), '^(P1|P2|All)_', 'once')));
hasdata = arrayfun(@(c) isnumeric(T.(c)) && any(~isnan(T.(c))), measurecols);
measurecols = measurecols(hasdata);
assert(~isempty(measurecols), 'summaryTable has no measure columns with data for these animals/conditions')
[sel,ok] = listdlg('ListString',cellstr(measurecols),'SelectionMode','multiple', ...
    'PromptString','Data features to plot (one column each):','ListSize',[300 300]);
assert(ok && ~isempty(sel), 'No data feature selected')
features = measurecols(sel);
nf = numel(features);

%% per animal: spon peak, alignment peak and their distance in shanks
SponLabel = strings(na,1);
SponPeak = NaN(na,1);
AlignPeak = NaN(na,1);
TableTopShank = NaN(na,1);
rankfield = "ShankRank_" + alignlabel;
for a = find(haslabels)
    sl = labelsets{a};
    if sponlabel == ""
        sp = sl(find(contains(sl,'Spon','IgnoreCase',true),1));
    else
        sp = sl(sl == sponlabel);
    end
    if isempty(sp)
        if sponlabel == "", sptxt = "containing Spon"; else, sptxt = sponlabel; end
        runnotes(end+1) = sprintf('%s: no spontaneous ShankLabels score (%s) - animal left out', animals(a), sptxt); %#ok<SAGROW>
        continue
    end
    SponLabel(a) = sp;
    %highest score, NaN ignored; max returns the first (lowest shank) of tied maxima
    sc = double(SLs{a}.(char(sp))(:)');
    if ~all(isnan(sc)), [~,SponPeak(a)] = max(sc); end
    sc = double(SLs{a}.(char(alignlabel))(:)');
    if ~all(isnan(sc)), [~,AlignPeak(a)] = max(sc); end
    if isnan(SponPeak(a)) || isnan(AlignPeak(a))
        runnotes(end+1) = sprintf('%s: all-NaN %s or %s shank scores - animal left out', animals(a), sp, alignlabel); %#ok<SAGROW>
    end
    %cross-check with the table's rank (the ProbeInfo labels can have changed since the table was made)
    if ismember(rankfield, vn)
        s1 = unique(T.Chan_Shank(T.Animal_Name == animals(a) & T.(rankfield) == 1));
        s1 = s1(~isnan(s1));
        if ~isempty(s1)
            TableTopShank(a) = min(s1);
            if TableTopShank(a) ~= AlignPeak(a)
                runnotes(end+1) = sprintf('%s: ProbeInfo %s peak is shank %d but the table''s %s == 1 shank is %d (ProbeInfo used)', ...
                    animals(a), alignlabel, AlignPeak(a), rankfield, TableTopShank(a)); %#ok<SAGROW>
            end
        end
    end
end
ShankDist = SponPeak - AlignPeak; %positive = alignment peak left of (lower shank number than) the spon peak
animalinfo = table(animals', SponLabel, SponPeak, AlignPeak, ShankDist, TableTopShank, ...
    'VariableNames', ["Animal_Name","SponLabel","SponPeak","AlignPeak","ShankDist","TableTopShank"]);

%% data points: the two conditions on Chan_Depth 1..maxdepth of each animal's alignment peak shank
keepcols = ["Animal_Name","Condition_Name","Chan_Shank","Chan_Depth"];
keepcols = [keepcols, intersect(["Channel_ID","Trial","Condition_FileNum"], vn, 'stable'), features];
ptdata = table();
for a = 1:na
    if isnan(ShankDist(a)), continue, end
    rows = T.Animal_Name == animals(a) & T.Chan_Shank == AlignPeak(a) & T.Chan_Depth >= 1 & T.Chan_Depth <= maxdepth;
    if ~any(rows)
        runnotes(end+1) = sprintf('%s: no rows of %s on shank %d, Chan_Depth 1-%d', animals(a), strjoin(conds,'/'), AlignPeak(a), maxdepth); %#ok<SAGROW>
        continue
    end
    sub = T(rows, keepcols);
    sub.ShankDist = repmat(ShankDist(a), height(sub), 1);
    ptdata = [ptdata; sub]; %#ok<AGROW> (one block per animal)
end
assert(~isempty(ptdata), 'No data points: check runnotes')
ptdata = movevars(ptdata, 'ShankDist', 'After', 'Animal_Name');

%points per condition and depth, per animal
for c = 1:2
    for d = 1:maxdepth
        cname = matlab.lang.makeValidName("n_" + conds(c) + "_d" + d);
        animalinfo.(cname) = arrayfun(@(a) nnz(ptdata.Animal_Name == a & ptdata.Condition_Name == conds(c) & ptdata.Chan_Depth == d), animalinfo.Animal_Name);
    end
end

%% statistics: per tile (Chan_Depth x feature) and shank distance, condition 2 vs condition 1
%y ~ Cond + (1|Animal) (fitlme) when >= 2 animals share the distance, else y ~ Cond (fitlm);
%p of the Cond coefficient, BH-FDR over the tested distances of each tile
stattable = table();
if do_stats
    xvals = min(ptdata.ShankDist):max(ptdata.ShankDist);
    ncol1 = matlab.lang.makeValidName("n_" + conds(1));
    ncol2 = matlab.lang.makeValidName("n_" + conds(2));
    for f = 1:nf
        feat = features(f);
        for d = 1:maxdepth
            tilerows = table();
            for x = xvals
                sub = ptdata(ptdata.Chan_Depth == d & ptdata.ShankDist == x & ~isnan(ptdata.(feat)), :);
                if isempty(sub), continue, end
                n1 = nnz(sub.Condition_Name == conds(1));
                n2 = nnz(sub.Condition_Name == conds(2));
                nanim = numel(unique(sub.Animal_Name));
                model = "";  est = NaN;  p = NaN;
                if n1 >= min_trials && n2 >= min_trials
                    %condition 1 is the reference level, so the Cond coefficient = cond2 - cond1
                    tbl = table(double(sub.(feat)), categorical(sub.Condition_Name, conds), categorical(sub.Animal_Name), ...
                        'VariableNames', ["y","Cond","Animal"]);
                    try
                        if nanim >= 2
                            mdl = fitlme(tbl, 'y ~ Cond + (1|Animal)');
                            mtype = "LME";
                        else
                            mdl = fitlm(tbl, 'y ~ Cond');
                            mtype = "LM";
                        end
                        est = mdl.Coefficients.Estimate(2);
                        p = mdl.Coefficients.pValue(2);
                        model = mtype;
                    catch err
                        runnotes(end+1) = sprintf('%s, Chan_Depth %d, distance %+d: model fit failed (%s) - not tested', ...
                            feat, d, x, err.message); %#ok<SAGROW>
                    end
                end
                tilerows = [tilerows; table(feat, d, x, nanim, n1, n2, model, est, p, ...
                    'VariableNames', ["Feature","Chan_Depth","ShankDist","nAnimals",ncol1,ncol2,"Model","Estimate","p"])]; %#ok<AGROW>
            end
            if isempty(tilerows), continue, end
            %Benjamini-Hochberg adjusted p over the tested distances of this tile (NaN = not tested)
            tilerows.q = NaN(height(tilerows),1);
            tested = find(~isnan(tilerows.p));
            m = numel(tested);
            if m > 0
                if use_fdr
                    [ps,ord] = sort(tilerows.p(tested));
                    qv = ps .* m ./ (1:m)';
                    qv = min(flipud(cummin(flipud(qv))), 1); %monotone, capped at 1
                    tilerows.q(tested(ord)) = qv;
                else
                    tilerows.q = tilerows.p;
                end
            end
            tilerows.Sig = tilerows.q < alpha;
            stattable = [stattable; tilerows]; %#ok<AGROW>
        end
    end
end

%% save_png: choose the file name before anything is plotted (the whole default name is shown and editable)
plotanimals = unique(ptdata.Animal_Name,'stable')';
if save_png
    if ~isfolder(savebase)
        mkdir(savebase);
    end
    savename = sprintf('n%d_%s_%svs%s_%sDist_%s', numel(plotanimals), region, conds(1), conds(2), alignlabel, strjoin(features,'_'));
    namehelp = sprintf('File name: <name>.png and <name>.txt\nFolder: %s', savebase);
    answer = inputdlg(namehelp, 'Save PNG', [1 80], {savename});
    if isempty(answer) || isempty(strtrim(answer{1}))
        error('plotShankDistance_ConditionCompare: stopped by the user (no file name)')
    end
    savename = strtrim(answer{1});
    %if a file exists: Overwrite (same-name files are replaced), Quit, or Edit Name (checked again)
    while true
        outnames = fullfile(savebase, savename + [".png"; ".txt"]);
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
                answer = inputdlg(sprintf('New file name:\n<name>.png and .txt in %s', savebase), ...
                    'Edit name', [1 80], {savename});
                if ~isempty(answer) && ~isempty(strtrim(answer{1}))
                    savename = strtrim(answer{1}); %checked again by the loop
                end
            otherwise %Quit, or the dialog was closed
                error('plotShankDistance_ConditionCompare: stopped by the user (%s.png / .txt already exist)', savename)
        end
    end
end

%% run information (figure title stays short; this goes to info and, with save_png, the .txt file)
ttl = sprintf('%s: %s vs %s on the %s peak shank, by distance from the spon peak (n = %d)', ...
    region, conds(1), conds(2), alignlabel, numel(plotanimals));
info = string(ttl);
info(end+1) = "created " + string(datetime('now','Format','yyyy-MM-dd HH:mm')) + " by plotShankDistance_ConditionCompare";
info(end+1) = "summary table: workspace variable summaryTable (" + height(summaryTable) + " rows)";
info(end+1) = sprintf('region %s; conditions %s (black) vs %s (red); features %s', region, conds(1), conds(2), strjoin(features,', '));
info(end+1) = sprintf('shank distance = spon peak shank - %s peak shank (+ = %s peak left of / lower shank number than the spon peak); peaks = highest ProbeInfo.ShankLabels score, ties -> lower shank', ...
    alignlabel, alignlabel);
pa = animalinfo(ismember(animalinfo.Animal_Name, plotanimals), :);
info(end+1) = "animals (spon label: spon peak / " + alignlabel + " peak -> distance): " + ...
    strjoin(pa.Animal_Name + " (" + pa.SponLabel + ": " + pa.SponPeak + " / " + pa.AlignPeak + " -> " + pa.ShankDist + ")", ', ');
info(end+1) = sprintf('points: every trial on Chan_Depth 1-%d of each animal''s %s peak shank (rows = Chan_Depth); x tick n = animals with points in that tile', ...
    maxdepth, alignlabel);
info(end+1) = "boxes: each condition pools all animals at that distance; box = quartiles, line = median, whiskers = most extreme points within 1.5 IQR; points jittered horizontally, one marker shape per animal";
if do_stats && isempty(stattable)
    info(end+1) = "statistics: no distance cells with data";
elseif do_stats
    if use_fdr, qtxt = 'BH-FDR adjusted p over the tested distances of each tile'; else, qtxt = 'p'; end
    info(end+1) = sprintf(['statistics: per tile and distance, %s vs %s on its trials: y ~ Cond + (1|Animal) (fitlme) when >= 2 animals share ' ...
        'the distance, else y ~ Cond (fitlm); p of the Cond coefficient (Estimate = %s - %s); at least %d points per condition; ' ...
        'number = %s; * < 0.05, ** < 0.01, *** < 0.001 (significant: q < %g); %d of %d distance cells tested (%d LME, %d LM), %d significant; ' ...
        'assumes roughly normal residuals (AUC values are often skewed)'], ...
        conds(2), conds(1), conds(2), conds(1), min_trials, qtxt, alpha, nnz(~isnan(stattable.q)), height(stattable), ...
        nnz(stattable.Model == "LME"), nnz(stattable.Model == "LM"), nnz(stattable.Sig));
else
    info(end+1) = "statistics: off";
end
if isempty(runnotes)
    info(end+1) = "runnotes: none";
else
    info = [info, "runnotes: " + runnotes'];
end

%% plot: rows = Chan_Depth, columns = features; x = shank distance, boxes + points per condition
xvals = min(ptdata.ShankDist):max(ptdata.ShankDist);
rs = RandStream('mt19937ar','Seed',1); %fixed jitter, global rng untouched

fig = figure('Name',sprintf('%s vs %s by %s - Spon peak distance (%s)', conds(1), conds(2), alignlabel, region), ...
    'Color','w','Units','pixels','Position',[60 60 max(450*nf,700) 900]);
tl = tiledlayout(fig, maxdepth, nf, 'TileSpacing','compact','Padding','compact');
title(tl, ttl, 'Interpreter','none');
xlabel(tl, sprintf('Spon peak shank - %s peak shank  (+ = %s peak left of spon peak)', alignlabel, alignlabel), 'Interpreter','none');

axs = gobjects(maxdepth, nf);
for f = 1:nf
    feat = features(f);
    for d = 1:maxdepth
        ax = nexttile(tl, (d-1)*nf + f);
        axs(d,f) = ax;
        hold(ax,'on');
        atdepth = ptdata.Chan_Depth == d & ~isnan(ptdata.(feat));
        for x = xvals
            for c = 1:2
                xc = x + (2*c-3)*boxoffset; %c = 1 -> left, c = 2 -> right
                y = ptdata.(feat)(atdepth & ptdata.ShankDist == x & ptdata.Condition_Name == conds(c));
                if isempty(y), continue, end
                %box and whisker: quartiles, median, whiskers to the most extreme points within 1.5 IQR
                q = quantile(y, [0.25 0.5 0.75]);
                iqrange = q(3) - q(1);
                wlo = min(y(y >= q(1) - 1.5*iqrange));
                whi = max(y(y <= q(3) + 1.5*iqrange));
                patch(ax, xc + boxwidth/2*[-1 1 1 -1], [q(1) q(1) q(3) q(3)], condcolors(c,:), ...
                    'FaceAlpha',0.12, 'EdgeColor',condcolors(c,:), 'LineWidth',1.2);
                plot(ax, xc + boxwidth/2*[-1 1], [q(2) q(2)], '-', 'Color',condcolors(c,:), 'LineWidth',2);
                plot(ax, [xc xc], [q(3) whi], '-', 'Color',condcolors(c,:), 'LineWidth',1);
                plot(ax, [xc xc], [q(1) wlo], '-', 'Color',condcolors(c,:), 'LineWidth',1);
                plot(ax, xc + boxwidth/4*[-1 1], [whi whi], '-', 'Color',condcolors(c,:), 'LineWidth',1);
                plot(ax, xc + boxwidth/4*[-1 1], [wlo wlo], '-', 'Color',condcolors(c,:), 'LineWidth',1);
                %points, one marker shape per animal
                for a = 1:numel(plotanimals)
                    ya = ptdata.(feat)(atdepth & ptdata.ShankDist == x & ptdata.Condition_Name == conds(c) & ptdata.Animal_Name == plotanimals(a));
                    if isempty(ya), continue, end
                    xj = xc + (rand(rs, numel(ya), 1) - 0.5) * jitterwidth;
                    scatter(ax, xj, ya, 22, condcolors(c,:), animalmarkers{mod(a-1,numel(animalmarkers))+1}, 'LineWidth',0.8);
                end
            end
        end
        %x ticks: every distance in the range, with n animals that have points in this tile
        nx = arrayfun(@(x) numel(unique(ptdata.Animal_Name(atdepth & ptdata.ShankDist == x))), xvals);
        xlim(ax, [xvals(1)-0.5 xvals(end)+0.5]);
        xticks(ax, xvals);
        ax.TickLabelInterpreter = 'tex';
        xticklabels(ax, arrayfun(@(x,n) sprintf('%+d\\newline n=%d', x, n), xvals, nx, 'UniformOutput', false));
        ylabel(ax, sprintf('Chan\\_Depth %d\n%s', d, strrep(feat,'_','\_')));
        if d == 1
            title(ax, feat, 'Interpreter','none');
        end
        grid(ax,'on'); ax.GridAlpha = 0.15; box(ax,'off');
        hold(ax,'off');
    end
    linkaxes(axs(:,f), 'y'); %same y scale down each feature column
end

%statistics: q (and stars if significant) above each tested distance, in room added at the top of each column
if do_stats && ~isempty(stattable)
    for f = 1:nf
        yl = ylim(axs(1,f));
        ylim(axs(1,f), [yl(1) yl(2) + 0.2*diff(yl)]); %y is linked: applies to the whole column
        yl = ylim(axs(1,f));
        for d = 1:maxdepth
            st = stattable(stattable.Feature == features(f) & stattable.Chan_Depth == d, :);
            ntested = nnz(~isnan(st.q));
            for k = find(~isnan(st.q))'
                if st.q(k) < 0.001, ptxt = '<0.001'; else, ptxt = sprintf('%.3f', st.q(k)); end
                if st.Sig(k)
                    if st.q(k) < 0.001, stars = '***'; elseif st.q(k) < 0.01, stars = '**'; else, stars = '*'; end
                    lbl = {['\bf' stars], ['\rm' ptxt]}; %stars (bold) above the q
                else
                    lbl = ptxt;
                end
                text(axs(d,f), st.ShankDist(k), yl(2), lbl, 'HorizontalAlignment','center', ...
                    'VerticalAlignment','top', 'FontSize',8, 'Interpreter','tex');
            end
            if ntested > 0
                if d == 1
                    title(axs(d,f), sprintf('%s  (%d/%d *)', features(f), nnz(st.Sig), ntested), 'Interpreter','none');
                else
                    title(axs(d,f), sprintf('%d/%d *', nnz(st.Sig), ntested), 'FontWeight','normal');
                end
            end
        end
    end
end

%legend: conditions (colour) and animals (marker shape, with their shank distance)
hold(axs(1,1),'on');
lh = gobjects(1, 2 + numel(plotanimals));
for c = 1:2
    lh(c) = plot(axs(1,1), NaN, NaN, 's', 'MarkerSize',9, 'MarkerFaceColor',condcolors(c,:), ...
        'MarkerEdgeColor',condcolors(c,:), 'DisplayName',conds(c));
end
for a = 1:numel(plotanimals)
    lh(2+a) = plot(axs(1,1), NaN, NaN, animalmarkers{mod(a-1,numel(animalmarkers))+1}, 'Color',[0.4 0.4 0.4], ...
        'LineWidth',0.8, 'DisplayName', sprintf('%s (%+d)', plotanimals(a), animalinfo.ShankDist(animalinfo.Animal_Name == plotanimals(a))));
end
hold(axs(1,1),'off');
lg = legend(lh, 'Interpreter','none', 'FontSize',8);
lg.Layout.Tile = 'east';

%save_png: <savebase>\<name>.png and .txt (the run information)
if save_png
    exportgraphics(fig, fullfile(savebase, savename + ".png"), 'Resolution', 150);
    fid = fopen(fullfile(savebase, savename + ".txt"), 'w');
    fprintf(fid, '%s\n', info);
    fclose(fid);
end

%% run report
disp(animalinfo)
if do_stats
    disp(stattable)
end
if isempty(runnotes)
    disp('runnotes: none')
else
    disp('runnotes:')
    disp(runnotes)
end
