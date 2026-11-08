% SummaryPlots_run.m
%
% Description: One place to make every plot from the summary table (the csv written by
%   SummaryAnalysis_MultiMouse_allChannels.m). Each %% section sets up and calls one plotting
%   function in Functions\Plotting. The sections are NOT meant to be run top to bottom:
%   run the Setup section first (loads the table), then run only the section(s) you want
%   (Ctrl+Enter in the section). Each section has its own little input block at its top.
%   In those blocks, [] means "ask me in a dialog" (not "" - an empty string is not empty to Matlab).
%   Nothing here saves anything, except the functions that have their own save_png option.
%
% Inputs:
%   User-set config (Setup section):
%     use_tableres (0/1) - 1 = use tableres already in the workspace (straight from the summary
%       script) instead of loading csvs
%     csv_directory (char) - folder with the summary csvs
%   Dialog (Setup section, when use_tableres = 0):
%     csvfiles (cell array of char, 1 x nFiles) - full paths of the csv(s) picked from a list of
%       every .csv in csv_directory; loaded and stacked
%     probeinfo_directory (char) - folder with <animal>-ProbeInfo.mat (probe-label section)
%   <csv_directory>\<csvfiles> - tableres from SummaryAnalysis_MultiMouse_allChannels.m, one row per
%     good trial x channel. Columns used by the plotting functions:
%       Animal_Name, Condition_Name, Region (text), Channel_ID (double)
%       Chan_Shank, Chan_Depth, ShankRank_<score> (double) - from the summary script's label section
%       measure columns named <P1|P2|All>_<LFP|MUA|TF>... (double) - the plotted features
%     Several csvs are stacked into one table, so they must have the same columns
%   <probeinfo_directory>\<animal>-ProbeInfo.mat - ProbeInfo (.Areas, .ProbeMaps, .Ch_Remove,
%     .ShankLabels, .ChanLabels); the plotting functions read it themselves (see their headers)
%
% Outputs (workspace, one set per section so sections don't overwrite each other):
%   summaryTable (table) - the loaded / stacked summary table
%   figs_labels (struct) - plotProbeLabels figures, one field per mouse
%   figs_LWint (figure handles) - plotLWIntensityCompare
%   figs_LWint_multi, data_LWint_multi - plotLWIntensityCompare_multianimal
%   figs_LWTTX, data_LWTTX - plotLWTTXCompare_multianimal
%   figs_shankprof, data_shankprof - plotRankedShankProfile
%   figs_shankrowprof, data_shankrowprof - plotRankedShankRows
%   (data_* struct fields: see each function's header)
%
% Dependencies: SummaryAnalysis_MultiMouse_allChannels.m (makes the csv, including its section that
%   adds Chan_Shank / Chan_Depth / ShankRank_<score>); the plotting functions in Functions\Plotting
%   (that folder must be on the Matlab path); NNx_editProbeInfo_ChLabel_SponData.m /
%   NNx_editProbeInfo_ChLabel_StimData.m (labels in ProbeInfo).



%% Setup: Load the Summary Table (ALWAYS RUN THIS FIRST)

%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY

%Is the table still in the workspace from the summary script?
use_tableres = 0;
    % 0 = load the csv(s) below
    % 1 = use tableres from the workspace

%Where are the summary csvs?
csv_directory = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR';

%Which csv(s)? you pick them from a list of the csvs in csv_directory when this section runs
%(more than one gets stacked into one table, so they need the same columns)
%typed file names replaced by the list dialog below, so nobody has to type them:
% csvfiles = {'LvsLW_variedIntensity_n4.csv'};
%     % {} = pick them in a file dialog
%     % e.g. {'LvsLW_variedIntensity_n4.csv','LvsLW_TTXtreatmentS1_n3.csv'}

%Where are the ProbeInfo files? (only the probe-label section uses this)
probeinfo_directory = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%%%nothing else to fill in%%%

if use_tableres == 1
    summaryTable = tableres; %the table straight from SummaryAnalysis_MultiMouse_allChannels.m
else
    %replaced by the list dialog below (no typed file names anymore):
    % if isempty(csvfiles) == true
    %     csvfiles = uipickfiles('FilterSpec',csv_directory); %full paths of the picked csvs
    % else
    %     for i = 1:numel(csvfiles)
    %         csvfiles{i} = fullfile(csv_directory, csvfiles{i}); %file name -> full path
    %     end
    % end

    %list every csv in the folder and let the user pick one or more
    csvlist = dir(fullfile(csv_directory,'*.csv')); %all csvs in csv_directory
    csvnames = {csvlist.name}; %cell {1 x number of csvs} of file names
    [sel, ok] = listdlg('ListString',csvnames,'SelectionMode','multiple','ListSize',[400 300], ...
        'Name','Pick summary csv(s)','PromptString','Which csv(s)? (ctrl+click for more than one)');
    if ok == 0
        error('No csv picked - rerun the Setup section and pick at least one')
    end
    csvfiles = fullfile(csv_directory, csvnames(sel)); %cell {1 x picked csvs} of full paths

    summaryTable = table();
    for i = 1:numel(csvfiles)
        summaryTable = cat(1, summaryTable, readtable(csvfiles{i},'TextType','string')); %stack csv i under the others
    end
end
disp(['Setup: summary table loaded, ' num2str(height(summaryTable)) ' rows from ' num2str(numel(unique(summaryTable.Animal_Name))) ' animals'])




%Everything below is one plot type per section. Run only the ones you want.



%% Data Visualization: Probe activity labels, single animals (plotProbeLabels)
%one figure per label and mouse: shank scores on top, channel-score heatmap (ProbeMaps layout) below
%doesn't need the summary table, just ProbeInfo

%Who is the data?
mice = {"20260423-p12"};

%Which brain area per mouse? (picks the probe: must match one entry of that mouse's ProbeInfo.Areas)
region = {"V1"};

%Which labels?
labelnames = strings(1,0);
    % strings(1,0) = every label found for that probe
    % e.g. ["SponActivity","LightOnly_mixedIntensity"]

figs_labels = struct();
for m = 1:numel(mice)
    prbinfo = load(fullfile(probeinfo_directory, [char(mice{m}) '-ProbeInfo.mat']),'ProbeInfo');
    prb = find(strcmp(prbinfo.ProbeInfo.Areas, region{m})); %probe whose area matches the region
    figs_labels.(matlab.lang.makeValidName(char(mice{m}))) = plotProbeLabels(prbinfo.ProbeInfo, prb, labelnames);
end
clear prbinfo



%% Data Visualization: L vs LW across light intensities, one animal (plotLWIntensityCompare)
%one figure per feature: activity labels on top, then L, LW and LW - L channel maps for 4, 8, 12, 15 uW

%Who is the data? (as in the Animal_Name column)
animal = "20260423-p12";

%Which region? (as in the Region column, or a probe number)
region_LWint = []; %int stands for intensities
    % [] = pick in a dialog

%Which features? (measure columns, one figure each)
features_LWint = [];
    % [] = pick in a dialog
    % e.g. ["All_MUAAUC","P1_LFPmin"]

%Which two stimulus labels in the top row (columns 2 and 3)?
stimlabels_LWint = [];
    % [] = pick in two dialogs (light only and whisker only pre-selected)

%Star channels where LW and L differ? (Wilcoxon rank-sum, FDR)
do_stats_LWint = 1;
    % 0 = no stats
    % 1 = test and star

%Save the figures as PNGs? (asks about the folder if it already exists)
save_png_LWint = 0;
    % 0 = just show them
    % 1 = also save

figs_LWint = plotLWIntensityCompare(animal, summaryTable, region_LWint, features_LWint, ...
    stimlabels_LWint, do_stats_LWint, save_png_LWint);



%% Data Visualization: L vs LW across light intensities, multi-animal (plotLWIntensityCompare_multianimal)
%one figure per feature: shanks of each animal re-centred on its top shank for alignlabel,
%labels on top, then L, LW and LW - L heatmaps (shank x depth) for 4, 8, 12, 15 uW

%Who is the data? (Animal_Name values)
animals_LWint_multi = [];
    % [] = pick in a dialog
    % e.g. ["20260423-p12","20260425-p14"]

%Which region?
region_LWint_multi = [];
    % [] = the only region of those animals, or a dialog

%Which features?
features_LWint_multi = [];
    % [] = pick in a dialog

%Which two stimulus labels in the top row (columns 2 and 3)?
stimlabels_LWint_multi = ["LightOnly","WhiskerOnly"];

%Which label's top shank is the centre column?
alignlabel_LWint_multi = "LightOnly";

%Star cells where LW and L differ? (Wilcoxon rank-sum, FDR, animals pooled)
do_stats_LWint_multi = 1;
    % 0 = no stats
    % 1 = test and star

%Draw the stretched row lines on the L / LW heatmaps?
rowlines_LWint_multi = 1;
    % 0 = plain heatmaps
    % 1 = row lines

%Only keep some trials? ("<column> <op> <value>", all combined with AND)
trialfilter_LWint_multi = [];
    % [] = keep everything
    % e.g. ["P1_MUAAUC > 100"]

%Save the figures as PNGs (+ a txt with the run info)? (asks for the file name)
save_png_LWint_multi = 0;
    % 0 = just show them
    % 1 = also save

[figs_LWint_multi, data_LWint_multi] = plotLWIntensityCompare_multianimal(animals_LWint_multi, summaryTable, ...
    region_LWint_multi, features_LWint_multi, stimlabels_LWint_multi, do_stats_LWint_multi, ...
    save_png_LWint_multi, alignlabel_LWint_multi, rowlines_LWint_multi, trialfilter_LWint_multi);



%% Data Visualization: L vs LW before and with TTX, multi-animal (plotLWTTXCompare_multianimal)
%one figure per feature: same alignment as above, then L_8 / LW_8 before and with TTX, their ratios
%and ratio comparisons. The table needs the L_8, LW_8, L_8_TTX, LW_8_TTX conditions

%Who is the data? (Animal_Name values)
animals_LWTTX = [];
    % [] = pick in a dialog

%Which region?
region_LWTTX = [];
    % [] = the only region of those animals, or a dialog

%Which features?
features_LWTTX = [];
    % [] = pick in a dialog

%Which two stimulus labels in the top row (columns 2 and 3)?
stimlabels_LWTTX = ["LightOnly","WhiskerOnly"];

%Which label's top shank is the centre column?
alignlabel_LWTTX = [];
    % [] = pick in a dialog (LightOnly pre-selected)

%Test the ratio tiles and ratio comparisons?
do_stats_LWTTX = 1;
    % 0 = no stats
    % 1 = test and star

%Draw the stretched row lines on the single-condition heatmaps?
rowlines_LWTTX = 1;
    % 0 = plain heatmaps
    % 1 = row lines

%Only keep some trials?
trialfilter_LWTTX = [];
    % [] = keep everything
    % e.g. ["P1_MUAAUC > 100"]

%Save the figures as PNGs (+ a txt with the run info)?
save_png_LWTTX = 0;
    % 0 = just show them
    % 1 = also save

[figs_LWTTX, data_LWTTX] = plotLWTTXCompare_multianimal(animals_LWTTX, summaryTable, region_LWTTX, features_LWTTX, ...
    stimlabels_LWTTX, do_stats_LWTTX, save_png_LWTTX, alignlabel_LWTTX, rowlines_LWTTX, trialfilter_LWTTX);



%% Data Visualization: Ranked shank profile, multi-animal (plotRankedShankProfile)
%one figure per feature: per-animal line plots across the 5 shanks around each animal's top shank
%(+ normalized version), and heatmaps of condA (or condB - condA)

%Who is the data? (Animal_Name values)
animals_shankprof = [];
    % [] = pick in a dialog

%Main condition (Condition_Name)
condA_shankprof = [];
    % [] = pick in a dialog
    % e.g. "L_8"

%Comparison condition (heatmap becomes B - A, dashed B lines)
condB_shankprof = [];
    % [] = no comparison
    % e.g. "LW_8"

%Which shank ranking re-centres the shanks? (ShankRank_<score>, prefix optional)
rankfield_shankprof = [];
    % [] = pick in a dialog
    % e.g. "SponActivity"

%Which features?
features_shankprof = [];
    % [] = pick in a dialog

%Which region?
region_shankprof = [];
    % [] = the only region of those animals, or a dialog

%Only keep some trials?
trialfilter_shankprof = [];
    % [] = keep everything
    % e.g. ["P1_MUAAUC > 100"]

%Star cells where A and B differ? (only with condB)
do_stats_shankprof = 0;
    % 0 = no stats
    % 1 = test and star

%What should the normalized line plot be normalized to?
normto_shankprof = 3;
    % 1-5 = that column (3 = top shank, 1/2 = -2/-1, 4/5 = +1/+2)
    % "peak" = each line's highest value
    % "none" = no normalized plot

[figs_shankprof, data_shankprof] = plotRankedShankProfile(summaryTable, animals_shankprof, condA_shankprof, condB_shankprof, ...
    rankfield_shankprof, features_shankprof, region_shankprof, trialfilter_shankprof, do_stats_shankprof, normto_shankprof);



%% Data Visualization: Ranked shank rows, multi-animal (plotRankedShankRows)
%one figure per feature: shank scores of the ranking label on top, then condA (or condB - condA)
%heatmaps with a normalized line across each row

%Who is the data? (Animal_Name values)
animals_shankrowprof = [];
    % [] = pick in a dialog

%Main condition (Condition_Name)
condA_shankrowprof = [];
    % [] = pick in a dialog
    % e.g. "L_8"

%Comparison condition (heatmap becomes B - A, dashed B row lines)
condB_shankrowprof = [];
    % [] = no comparison
    % e.g. "LW_8"

%Which shank ranking re-centres the shanks? (must also be a ProbeInfo.ShankLabels field)
rankfield_shankrowprof = [];
    % [] = pick in a dialog
    % e.g. "SponActivity"

%Which features?
features_shankrowprof = [];
    % [] = pick in a dialog

%Which region?
region_shankrowprof = [];
    % [] = the only region of those animals, or a dialog

%Only keep some trials?
trialfilter_shankrowprof = [];
    % [] = keep everything
    % e.g. ["P1_MUAAUC > 100"]

%Star cells where A and B differ? (only with condB)
do_stats_shankrowprof = 0;
    % 0 = no stats
    % 1 = test and star

%What should the row lines be normalized to?
normto_shankrowprof = 3;
    % 1-5 = that column (3 = top shank, 1/2 = -2/-1, 4/5 = +1/+2)
    % "peak" = each line's highest value
    % "none" = no row lines

%Divide each animal's shank scores (top plot) by its top shank's score?
scorenorm_shankrowprof = 0;
    % 0 = raw scores
    % 1 = top shank = 1

[figs_shankrowprof, data_shankrowprof] = plotRankedShankRows(summaryTable, animals_shankrowprof, condA_shankrowprof, condB_shankrowprof, rankfield_shankrowprof, ...
    features_shankrowprof, region_shankrowprof, trialfilter_shankrowprof, do_stats_shankrowprof, normto_shankrowprof, scorenorm_shankrowprof);



%% Data Visualization: Shank distance vs condition difference, multi-animal (plotShankDistance_ConditionCompare)
%will be fixed later - plotShankDistance_ConditionCompare.m is a script, not a function, and its first
%section still loads its own hardcoded csvs into summaryTable (it would overwrite the one from Setup)
%and its settings section resets everything set here. Until then, open it and run it there:
%   edit plotShankDistance_ConditionCompare
%
% %How many channels down the alignment peak shank? (rows of the figure)
% maxdepth = 3;
%
% %Which ShankLabels field is spontaneous activity?
% sponlabel = "";
%     % "" = the first field containing "Spon"
%
% %Test condition 2 vs condition 1 at each shank distance? (mixed model, FDR)
% do_stats = true;
%
% %Save the figure as a PNG (+ a txt with the run info)?
% save_png = false;
%
% run('plotShankDistance_ConditionCompare.m') %uses summaryTable from Setup; dialogs ask for animals, conditions, features
