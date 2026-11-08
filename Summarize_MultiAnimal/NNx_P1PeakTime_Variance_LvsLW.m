% NNx_P1PeakTime_Variance_LvsLW.m
%
% Description: quick check of trial-to-trial jitter in the P1 peak time
%   (proxy for how well trials line up) between Light and Light+Whisker.
%   Hypothesis: adding whisker makes the response onset more consistent, so
%   P1 peak time spread should be SMALLER in LW than in L.
%   You pick which animals to use, and which Condition_Name entries are L and which are LW, from lists.
%   Same number of each -> paired in list order (e.g. L_4 <-> LW_4) and compared
%   pair by pair, so intensities and TTX never get mixed; otherwise all picked L
%   vs all picked LW as one pair. Only trials with a response (P1_MUAAUC above
%   mua_thresh) are used. Spread is computed within each animal x channel x pair
%   ("unit") so between-animal/channel latency differences don't inflate it,
%   then units are compared paired (L vs LW).
%   Optionally (heatplot) the spread is also drawn as the ranked-shank summary heatmaps
%   (NNx_plotLWIntensityCompare/NNx_plotLWTTXCompare_multianimal), shanks centred on the top shank
%   of each shank label you pick (one figure per label), each cell = mean absolute deviation
%   of the P1 peak time.
%
% Inputs:
%   one or more csvs from NNx_SummaryAnalysis_MultiMouse_allChannels.m (one row per trial x channel),
%   stacked into one table; they must have the same column names (any order), otherwise the
%   script stops with an error. Columns used:
%     P1_LFPPeakTime (double) - time of the negative LFP peak in the P1 window (75-350 ms),
%       ms post stimulus; NaN when the window has no local peak
%     Animal_Name (string) - mouse id
%     Condition_Name (string) - condition label (e.g. 'L_8', 'LW_8', 'L_8_TTX', 'W')
%     Channel_ID (double) - raw channel ID (unique across probes)
%     P1_MUAAUC (double) - integral of the firing rate over the P1 window, spikes/s*ms;
%       NaN if the spiking file was missing (those trials are dropped by the response filter)
%     P1_LFPRMSsum (double) - area under the LFP rms envelope over the P1 window, uV*ms; trials
%       above rms_max or NaN are dropped (rmsfilter)
%
% Outputs (workspace only, nothing saved):
%   tableres.Source_CSV (string) - name of the csv each row came from
%   frac_kept_respfilter (double) - fraction of picked-condition trials that passed the response filter
%   frac_kept_rmsfilter (double) - fraction of those that also passed the P1_LFPRMSsum <= rms_max filter
%   unitres (table, units x columns) - per-unit SD, IQR, n trials, no-peak fraction, Brown-Forsythe p
%   stats_paired (table, pairs+1 x columns) - per pair (and all pairs together): median SDs and
%     one-sided signrank on SD and IQR across units
%   stats_pooled (struct) - Brown-Forsythe on unit-median-centered peak times, all units pooled
%   (heatplot > 0) tableres.P1_LFPPeakTimeAbsDev (double, ms) - per trial |peak time - median| within
%     its animal x channel x condition; NaN with fewer than 2 peak trials
%   (heatplot > 0) alignlabels (string, 1 x picks) - ShankRank_ labels chosen for alignment
%   (heatplot > 0) heatfigs, heatdata (cell, 1 x alignlabels) - figure handles and plotted numbers
%     returned by NNx_plotLWIntensityCompare_multianimal.m (heatplot 1) or NNx_plotLWTTXCompare_multianimal.m
%     (heatplot 2), one per alignment label; each heatmap cell = mean absolute deviation of P1 peak time, ms
%   stats_paired also has the mixed model columns: lme_n_animals, lme_SDratio (SD_LW / SD_L estimated by
%     fitlme (REML) on log(SD) ~ Cond (+ Pair for "all") + (1 + Cond|Animal) + (1|ChanKey); <1 = LW tighter),
%     lme_SDratio_CIlo/hi (95% CI), lme_p (one-sided, LW tighter); CI and p use Satterthwaite DF.
%     lme_models (cell) - the fitted models
%   (chanmaps == 1) chanfigs (figure handles, pairs x 1) - one figure per pair, one probe-layout tile per
%     animal: per channel SD_LW - SD_L (ms), dark grey = "Not Enough Trials" (in the csv, but < 2 trials
%     left in L or LW after the filters), black = removed / not in the csv, * = Brown-Forsythe q < 0.05
%     (FDR within animal)
%   (chanmaps == 1) stats_animal (table, pairs x animals rows) - per animal: n channels, median SDs,
%     channels LW tighter, one-sided signrank p across its channels, n significant channels,
%     n Not Enough Trials channels
%   (chanmaps == 1) chanmap_notes (string) - problems hit while drawing (missing ProbeInfo, unmapped channels)
%
% Dependencies: NNx_SummaryAnalysis_MultiMouse_allChannels.m csv output (heatmaps also need its
%   Chan_Shank, Chan_Depth, Region and ShankRank_<label> columns, and each animal's ProbeInfo).
%   NNx_plotLWIntensityCompare_multianimal.m / NNx_plotLWTTXCompare_multianimal.m for the heatmaps.
%   <probeinfo_dir>\<animal>-ProbeInfo.mat (.Areas, .ProbeMaps) for the per-animal channel maps.
%   Statistics Toolbox (vartestn, signrank, iqr, ranksum, fitlme). MATLAB R2022a+ (clim).
%
% Inferred data contract: Channel_ID values are the numbers in ProbeInfo.ProbeMaps{probe}
%   (true for the 64-channel V1 probes checked; channels not found go into chanmap_notes)
%

%% Setup: pick the csv and conditions

%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY

%where are the csv files?
csv_directory = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR';

%commented out: you now pick the L and LW conditions from a list after the csv loads
%what do the light and light+whisker condition names START with in Condition_Name?
%condi_L = "L";
%condi_LW = "LW";
condi_L = "L"; %just the labels used in tables and plots now
condi_LW = "LW";

%%%Analysis details%%%
%only use trials with a response? (MUA in the P1 window above a threshold)
respfilter = 1;
    % 0 = use all trials
    % 1 = only trials with P1_MUAAUC > mua_thresh
mua_thresh = 100; %spikes/s*ms, 0 = at least one spike in the P1 window

%drop trials with a huge LFP in the P1 window? (artifacts)
rmsfilter = 0;
    % 0 = keep them
    % 1 = drop trials with P1_LFPRMSsum > rms_max or NaN
rms_max = 1500; %area under the LFP rms envelope over the P1 window (uV*ms)

%Want the ranked-shank summary heatmaps of the P1 peak time spread too?
heatplot = 1;
    % 0 = no heatmaps
    % 1 = intensity layout (NNx_plotLWIntensityCompare_multianimal: L/LW 4, 8, 12, 15)
    % 2 = TTX layout (NNx_plotLWTTXCompare_multianimal: L/LW 8 with and without TTX)
%save the heatmaps as PNG (+ txt of run info) in the MultAnimal folder?
save_heatpng = 0;
    % 0 = no
    % 1 = yes (you'll get a file name box)

%Want the per-animal channel maps? (one figure per L/LW pair, one probe-layout tile per animal)
chanmaps = 1;
    % 0 = no
    % 1 = yes
%where are the ProbeInfo files? (for each animal's probe layout)
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';


%commented out: the allChannels csv finds P1 as a true local peak and gives NaN when there
%isn't one, so "max() landed on the window edge" can't happen anymore. NaN trials are dropped below instead
%P1 window edges in ms after stim (to flag peaks that just hit the edge)
%P1_edges = [75 350];

%drop trials whose peak hit the window edge? (those are usually "no real peak" trials)
%dropedges = 1;
    % 0 = keep them
    % 1 = drop them

%%%nothing else to fill in%%%

%csv_file = uipickfiles('FilterSpec',csv_directory,'Prompt','Pick the summary csv','REFilter','\.csv$','NumFiles',1); %one csv only, now you can stack several
%tableres = readtable(csv_file{1},'TextType','string'); %one row per trial x channel
csv_file = uipickfiles('FilterSpec',csv_directory,'Prompt','Pick the summary csv(s)','REFilter','\.csv$');

%load each csv and stack them into one table (same columns needed; their order can differ)
tableres = [];
for i = 1:numel(csv_file)
    filei = readtable(csv_file{i},'TextType','string'); %one row per trial x channel
    [~,csvname] = fileparts(csv_file{i});
    filei.Source_CSV = repmat(string(csvname),height(filei),1); %which csv the row came from
    if i == 1
        tableres = filei;
    else
        if isempty(setxor(tableres.Properties.VariableNames,filei.Properties.VariableNames)) == false
            error(['These csvs do not have the same columns: ' csvname ' vs the first one picked'])
        end
        tableres = [tableres; filei(:,tableres.Properties.VariableNames)]; %put columns in the first csv's order, then stack
    end
end

%which animals to use? (all pre-selected)
animal_list = unique(tableres.Animal_Name,'stable');
idxA = listdlg('ListString',animal_list,'SelectionMode','multiple','Name','Animals', ...
    'PromptString','Which animals should I use?','ListSize',[300 200],'InitialValue',1:numel(animal_list));
pick_animals = animal_list(idxA);
tableres = tableres(ismember(tableres.Animal_Name,pick_animals),:); %everything below only sees these animals


%% Data PreProcessing: keep only L and LW rows, and pair them up

%tableres = tableres(tableres.Condition_Name == condi_L | tableres.Condition_Name == condi_LW,:); %exact names don't exist (L_8, LW_8_TTX...), match by start instead
%commented out: matching by name start replaced by picking the conditions from a list
%isLWrow = startsWith(tableres.Condition_Name,condi_LW);
%isLrow = startsWith(tableres.Condition_Name,condi_L) & isLWrow == false; %starts with L but not LW
%tableres = tableres(isLrow | isLWrow,:);
%tableres.LvsLW = repmat(condi_L,height(tableres),1);
%tableres.LvsLW(startsWith(tableres.Condition_Name,condi_LW)) = condi_LW;
%tableres.pair = extractAfter(tableres.Condition_Name,strlength(tableres.LvsLW)); %strip the L/LW off the front

%which conditions are in this csv?
condi_list = unique(tableres.Condition_Name,'stable');

%pick L and LW conditions. If you pick the same number of each, they get paired in list order
%(1st L with 1st LW, ...) and compared pair by pair. Otherwise all L vs all LW as one pair
ishsokay = 0;
while ishsokay == 0
    idxL = listdlg('ListString',condi_list,'SelectionMode','multiple','Name','Light', ...
        'PromptString','Which conditions are LIGHT (L)?','ListSize',[300 200]);
    idxLW = listdlg('ListString',condi_list,'SelectionMode','multiple','Name','Light + Whisker', ...
        'PromptString','Which conditions are LIGHT + WHISKER (LW)?','ListSize',[300 200]);
    pick_L = condi_list(idxL);
    pick_LW = condi_list(idxLW);

    if numel(pick_L) == numel(pick_LW)
        pairnames = pick_L + " vs " + pick_LW; %one pair per L/LW match
    else
        pairnames = strjoin(pick_L',"+") + " vs " + strjoin(pick_LW',"+"); %everything lumped together
    end

    quest = ['Comparing these:' newline newline char(strjoin(pairnames',newline)) newline newline 'Look right?'];
    btn1 = 'Yes, happy as a clam';
    btn2 = 'No, lemme pick again';
    answer = questdlg(quest,'Check your pairs',btn1,btn2,btn1);
    switch answer
        case btn1
            ishsokay = 1;
        case btn2
            ishsokay = 0;
    end
end

%keep only the picked conditions, label each row L or LW, and which pair it belongs to
tableres = tableres(ismember(tableres.Condition_Name,[pick_L; pick_LW]),:);
tableres.LvsLW = repmat(condi_L,height(tableres),1);
tableres.LvsLW(ismember(tableres.Condition_Name,pick_LW)) = condi_LW;
tableres.pair = repmat(pairnames(1),height(tableres),1); %right if everything is lumped into one pair
if numel(pick_L) == numel(pick_LW)
    for j = 1:numel(pick_L)
        tableres.pair(tableres.Condition_Name == pick_L(j) | tableres.Condition_Name == pick_LW(j)) = pairnames(j);
    end
end
pairs = pairnames(:);

%which animal x channel x pair combos are in the csv before the trial filters (for the
%"Not Enough Trials" squares in the per-animal channel maps)
prefilter_keys = unique(tableres.Animal_Name + "_ch" + string(tableres.Channel_ID) + "_" + tableres.pair);

%only trials with a response (some MUA in the P1 window)
ntrials_before = height(tableres);
if respfilter == 1
    tableres = tableres(tableres.P1_MUAAUC > mua_thresh,:); %NaN MUA (no spiking file) fails this too
end
frac_kept_respfilter = height(tableres) / ntrials_before; %for the report

%drop artifact trials (huge LFP in the P1 window) and trials with no P1_LFPRMSsum (NaN)
ntrials_before = height(tableres);
if rmsfilter == 1
    tableres = tableres(tableres.P1_LFPRMSsum <= rms_max,:); %NaN fails this too, so it's dropped
end
frac_kept_rmsfilter = height(tableres) / ntrials_before; %for the report, out of the trials left after the response filter

%commented out: see setup, edge peaks don't happen in the allChannels csv
%tableres.isedge = tableres.P1_PeakTime <= P1_edges(1) | tableres.P1_LFPPeakTime >= P1_edges(2); %peak sat on the window edge
%if dropedges == 1
%    tableres = tableres(tableres.isedge == 0,:);
%end
tableres.nopeak = isnan(tableres.P1_LFPPeakTime); %no P1 peak found on this trial, dropped per unit below

%a "unit" is one animal x channel x pair
tableres.unit = tableres.Animal_Name + "_ch" + string(tableres.Channel_ID) + "_" + tableres.pair;
units = unique(tableres.unit,'stable');
nunits = numel(units);


%% Data Analysis: per-unit spread of P1 peak time

sd_L = nan(nunits,1); sd_LW = nan(nunits,1); %standard deviation, ms
iqr_L = nan(nunits,1); iqr_LW = nan(nunits,1); %interquartile range, ms (less bothered by outliers)
n_L = zeros(nunits,1); n_LW = zeros(nunits,1); %number of trials with a peak
nopeak_L = nan(nunits,1); nopeak_LW = nan(nunits,1); %fraction of trials with no P1 peak (NaN)
p_bf = nan(nunits,1); %Brown-Forsythe p (equal spread test) within the unit
unit_animal = strings(nunits,1); unit_chan = zeros(nunits,1); unit_pair = strings(nunits,1);

tableres.centered = nan(height(tableres),1); %peak time minus that unit+condition's median, for pooling later

for i = 1:nunits
    isunit = tableres.unit == units(i);
    isL = isunit & tableres.LvsLW == condi_L;
    isLW = isunit & tableres.LvsLW == condi_LW;

    first = find(isunit,1);
    unit_animal(i) = tableres.Animal_Name(first);
    unit_chan(i) = tableres.Channel_ID(first);
    unit_pair(i) = tableres.pair(first);
    nopeak_L(i) = mean(tableres.nopeak(isL));
    nopeak_LW(i) = mean(tableres.nopeak(isLW));

    isL = isL & tableres.nopeak == false; %only trials with a real peak from here on
    isLW = isLW & tableres.nopeak == false;
    tL = tableres.P1_LFPPeakTime(isL); %ms
    tLW = tableres.P1_LFPPeakTime(isLW); %ms
    n_L(i) = numel(tL);
    n_LW(i) = numel(tLW);

    if n_L(i) > 1 && n_LW(i) > 1 %skip units missing a condition
        sd_L(i) = std(tL);
        sd_LW(i) = std(tLW);
        iqr_L(i) = iqr(tL);
        iqr_LW(i) = iqr(tLW);
        p_bf(i) = vartestn([tL; tLW],[repmat(condi_L,n_L(i),1); repmat(condi_LW,n_LW(i),1)], ...
            'TestType','BrownForsythe','Display','off');

        tableres.centered(isL) = tL - median(tL); %take out the unit's latency so units can be pooled
        tableres.centered(isLW) = tLW - median(tLW);
    end
end

unitres = table(units,unit_animal,unit_chan,unit_pair,n_L,n_LW,nopeak_L,nopeak_LW,sd_L,sd_LW,sd_LW-sd_L,iqr_L,iqr_LW,p_bf, ...
    'VariableNames',["Unit","Animal_Name","Channel_ID","Pair","n_L","n_LW","NoPeakFrac_L","NoPeakFrac_LW", ...
    "SD_L","SD_LW","SD_diff_LWminusL","IQR_L","IQR_LW","p_BrownForsythe"]);


%% Data Analysis: across units (paired, per pair and all together) and pooled

good = isnan(sd_L) == false;

%paired across units: is LW spread smaller than L? (one-sided, matches the hypothesis)
%one row per pair, last row = all pairs together
pairlist = [pairs; "all"];
npl = numel(pairlist);
n_units = zeros(npl,1); median_SD_L = nan(npl,1); median_SD_LW = nan(npl,1);
n_units_LWtighter = zeros(npl,1); p_SD = nan(npl,1); p_IQR = nan(npl,1);
for j = 1:npl
    if pairlist(j) == "all"
        use = good;
    else
        use = good & unit_pair == pairlist(j);
    end
    n_units(j) = sum(use);
    if n_units(j) > 0
        median_SD_L(j) = median(sd_L(use));
        median_SD_LW(j) = median(sd_LW(use));
        n_units_LWtighter(j) = sum(sd_LW(use) < sd_L(use));
        p_SD(j) = signrank(sd_LW(use),sd_L(use),'tail','left');
        p_IQR(j) = signrank(iqr_LW(use),iqr_L(use),'tail','left');
    end
end
%mixed model alongside the signrank: same per-unit SDs, but animal and channel are random effects,
%so channels from one mouse (and one channel showing up in several pairs) aren't counted as independent.
%log(SD) because a whisker effect should scale the spread (LW = some fraction of L), not shift it by fixed ms
%long format: one row per unit x condition (units with SD > 0 in both, log(0) doesn't work)
uselme = good & sd_L > 0 & sd_LW > 0;
nlme = sum(uselme);
lmetab = table([log(sd_L(uselme)); log(sd_LW(uselme))], ...
    categorical([repmat(condi_L,nlme,1); repmat(condi_LW,nlme,1)],[condi_L condi_LW]), ... %L is the reference
    categorical([unit_pair(uselme); unit_pair(uselme)]), ...
    [unit_animal(uselme); unit_animal(uselme)], ...
    [unit_animal(uselme) + "_ch" + string(unit_chan(uselme)); unit_animal(uselme) + "_ch" + string(unit_chan(uselme))], ...
    'VariableNames',["logSD","Cond","Pair","Animal","ChanKey"]);

lme_n_animals = zeros(npl,1); %animals in that model
lme_SDratio = nan(npl,1); %SD_LW / SD_L the model estimates (exp of the Cond coefficient), <1 = LW tighter
lme_SDratio_CIlo = nan(npl,1); lme_SDratio_CIhi = nan(npl,1); %95% confidence interval of that ratio
lme_p = nan(npl,1); %one-sided p (LW tighter), from the two-sided p of the Cond coefficient
lme_models = cell(npl,1); %the fitted models, if you want to look closer
for j = 1:npl
    %(1 + Cond|Animal): each animal gets its own LW effect (random slope), so the LW effect is tested
    %against how much it varies between mice, not just between channels. Without it the DF is
    %~ number of channels and the p mostly says "true for these mice's channels", not "true across mice"
    if pairlist(j) == "all"
        sub = lmetab;
        formula = 'logSD ~ Cond + Pair + (1 + Cond|Animal) + (1|ChanKey)'; %pair as a fixed effect, LW effect averaged over pairs
    else
        sub = lmetab(lmetab.Pair == pairlist(j),:);
        formula = 'logSD ~ Cond + (1 + Cond|Animal) + (1|ChanKey)';
    end
    if height(sub) > 0
        lme_n_animals(j) = numel(unique(sub.Animal));
        if pairlist(j) == "all" && numel(unique(sub.Pair)) == 1
            formula = 'logSD ~ Cond + (1 + Cond|Animal) + (1|ChanKey)'; %only one pair, nothing to put in for Pair
        end
        mdl = fitlme(sub,formula,'FitMethod','REML'); %REML: better variance estimates with few animals
        lme_models{j} = mdl;
        %Satterthwaite degrees of freedom: the default (residual DF ~ number of rows) acts like
        %every channel is independent, which is what this model is trying to avoid
        [~,~,fe] = fixedEffects(mdl,'DFMethod','satterthwaite');
        row = strcmp(fe.Name,['Cond_' char(condi_LW)]);
        est = fe.Estimate(row);
        lme_SDratio(j) = exp(est);
        lme_SDratio_CIlo(j) = exp(fe.Lower(row));
        lme_SDratio_CIhi(j) = exp(fe.Upper(row));
        if est < 0 %two-sided p -> one-sided (LW < L)
            lme_p(j) = fe.pValue(row) / 2;
        else
            lme_p(j) = 1 - fe.pValue(row) / 2;
        end
    end
end

stats_paired = table(pairlist,n_units,median_SD_L,median_SD_LW,n_units_LWtighter,p_SD,p_IQR, ...
    lme_n_animals,lme_SDratio,lme_SDratio_CIlo,lme_SDratio_CIhi,lme_p, ...
    'VariableNames',["Pair","n_units","median_SD_L","median_SD_LW","n_units_LWtighter","p_SD","p_IQR", ...
    "lme_n_animals","lme_SDratio","lme_SDratio_CIlo","lme_SDratio_CIhi","lme_p"]);

%pooled: all trials, each unit centered on its own median first
isc = isnan(tableres.centered) == false;
stats_pooled.SD_L = std(tableres.centered(isc & tableres.LvsLW == condi_L));
stats_pooled.SD_LW = std(tableres.centered(isc & tableres.LvsLW == condi_LW));
stats_pooled.p_BrownForsythe = vartestn(tableres.centered(isc),tableres.LvsLW(isc), ...
    'TestType','BrownForsythe','Display','off');

%NOTE: units aren't independent (many channels per animal, same channels across pairs),
%so treat the signrank and pooled p values as a rough check; lme_p (mixed model) is the one that
%accounts for animal and channel


%% Data Visualization: paired SD per unit + pooled distributions

colors = lines(numel(pairs));

figure()
subplot(1,2,1)
hold on
for i = find(good)'
    plot([1 2],[sd_L(i) sd_LW(i)],'-','Color',[colors(pairs == unit_pair(i),:) .15]) %one faint line per unit, colored by pair
end
h = gobjects(numel(pairs),1);
for j = 1:numel(pairs)
    h(j) = plot([1 2],[median_SD_L(j) median_SD_LW(j)],'-o','Color',colors(j,:),'LineWidth',2.5,'MarkerFaceColor',colors(j,:)); %median per pair
end
hold off
xlim([.5 2.5])
xticks([1 2]); xticklabels([condi_L condi_LW])
ylabel('P1 peak time SD (ms)')
legend(h,pairs,'Interpreter','none')
title(['per unit, all pairs signrank p = ' num2str(p_SD(end),3)])

subplot(1,2,2)
hold on
histogram(tableres.centered(isc & tableres.LvsLW == condi_L),'BinWidth',5,'Normalization','probability','DisplayStyle','stairs','LineWidth',1.5)
histogram(tableres.centered(isc & tableres.LvsLW == condi_LW),'BinWidth',5,'Normalization','probability','DisplayStyle','stairs','LineWidth',1.5)
xline(0)
hold off
legend([condi_L condi_LW])
xlabel('P1 peak time - unit median (ms)')
ylabel('fraction of trials')
title(['pooled, Brown-Forsythe p = ' num2str(stats_pooled.p_BrownForsythe,3)])
sgtitle('P1 peak time jitter: Light vs Light+Whisker')


%% Data Visualization: ranked-shank heatmaps of P1 peak time spread

%Same heatmaps as NNx_plotLWIntensityCompare/NNx_plotLWTTXCompare_multianimal, fed a per-trial spread column:
%P1_LFPPeakTimeAbsDev = |peak time - median peak time| of that trial's animal x channel x condition.
%The heatmaps average it per cell, so each cell = mean absolute deviation (ms, smaller = trials line up
%better), and their per-cell LW vs L rank-sum test on it is a test of spread (same idea as Brown-Forsythe)
if heatplot > 0
    tableres.P1_LFPPeakTimeAbsDev = nan(height(tableres),1);
    devkey = tableres.Animal_Name + "_ch" + string(tableres.Channel_ID) + "_" + tableres.Condition_Name;
    devkeys = unique(devkey(tableres.nopeak == false));
    for i = 1:numel(devkeys)
        rows = devkey == devkeys(i) & tableres.nopeak == false;
        if sum(rows) > 1 %one trial has no spread
            t = tableres.P1_LFPPeakTime(rows); %ms
            tableres.P1_LFPPeakTimeAbsDev(rows) = abs(t - median(t));
        end
    end

    %which shank label ranking(s) should the shanks be centred on? one figure per pick
    rankcols = string(tableres.Properties.VariableNames(startsWith(tableres.Properties.VariableNames,"ShankRank_")));
    firstpick = find(rankcols == "ShankRank_SponActivity",1); %spon activity ranking pre-selected
    if isempty(firstpick)
        firstpick = 1;
    end
    [sel,~] = listdlg('ListString',cellstr(rankcols),'SelectionMode','multiple','Name','Shank alignment', ...
        'PromptString','Centre the shanks on the top shank of which label(s)?','ListSize',[300 150], ...
        'InitialValue',firstpick);
    alignlabels = rankcols(sel);

    heatanimals = pick_animals; %the animals you picked at the start
    heatfigs = cell(1,numel(alignlabels)); %{1 x alignment labels}, one plot function output each
    heatdata = cell(1,numel(alignlabels));
    %with save_heatpng, every label gets the same default file name (it doesn't include the
    %alignment), so for the 2nd label on pick Edit Name in the file-exists box and add the label
    for k = 1:numel(alignlabels)
        disp(['Heatmaps: centring on ' char(alignlabels(k))])
        if heatplot == 1
            [heatfigs{k}, heatdata{k}] = NNx_plotLWIntensityCompare_multianimal(heatanimals, tableres, [], "P1_LFPPeakTimeAbsDev", ...
                [], true, save_heatpng == 1, alignlabels(k), true, []);
        elseif heatplot == 2
            [heatfigs{k}, heatdata{k}] = NNx_plotLWTTXCompare_multianimal(heatanimals, tableres, [], "P1_LFPPeakTimeAbsDev", ...
                [], true, save_heatpng == 1, alignlabels(k), true, []);
        end
        set(heatfigs{k},'Name',"P1 spread, aligned on " + alignlabels(k)) %figure titles look the same, so label the windows
    end
    %no trial filter passed: tableres is already MUA-filtered and only has your picked conditions
end


%% Data Visualization: per-animal channel maps of the LW - L spread difference

%One figure per L/LW pair, one tile per animal laid out like its probe (ProbeInfo.ProbeMaps,
%rows = channel position top to bottom, columns = shanks). Each channel with data shows
%SD_LW - SD_L of its P1 peak time (ms; blue = LW tighter), the rest are black.
%Star = that channel's Brown-Forsythe test (L vs LW trials) significant after FDR within the animal.
%Tile title = per-animal signrank across its channels (one-sided, LW tighter), like stats_paired but one mouse at a time
if chanmaps == 1
    chanfigs = gobjects(numel(pairs),1);
    chanmap_notes = strings(0,1); %problems for the report (missing ProbeInfo, channels not on the map)
    animal_res = cell(numel(pairs),1); %per pair, a table with one row per animal

    %blue-white-red, white = 0
    bwr = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
    greyc = [0.4 0.4 0.4]; %dark grey = channel in the csv but too few trials left after the filters

    for j = 1:numel(pairs)
        inpair = good & unit_pair == pairs(j);
        pair_animals = unique(unit_animal(inpair),'stable');
        nanimals = numel(pair_animals);
        dmax = max(abs(sd_LW(inpair) - sd_L(inpair))); %one symmetric colour scale for all tiles of this figure

        a_n = zeros(nanimals,1); a_medL = nan(nanimals,1); a_medLW = nan(nanimals,1); a_tighter = zeros(nanimals,1);
        a_p = nan(nanimals,1); a_nsig = zeros(nanimals,1); a_nfew = zeros(nanimals,1);

        chanfigs(j) = figure('Name',"P1 spread per animal: " + pairs(j));
        set(gcf, 'Position', get(0, 'Screensize'));
        tl = tiledlayout('flow','TileSpacing','compact');
        title(tl,"P1 peak time SD, LW - L (ms), per channel: " + pairs(j),'Interpreter','none')

        for a = 1:nanimals
            animal = pair_animals(a);
            inA = inpair & unit_animal == animal; %this animal's channels in this pair
            ax = nexttile;

            %per-animal signrank across its channels
            a_n(a) = sum(inA);
            a_medL(a) = median(sd_L(inA));
            a_medLW(a) = median(sd_LW(inA));
            a_tighter(a) = sum(sd_LW(inA) < sd_L(inA));
            a_p(a) = signrank(sd_LW(inA),sd_L(inA),'tail','left');

            %FDR (Benjamini-Hochberg) over this animal's per-channel Brown-Forsythe p values
            pch = p_bf(inA);
            qch = nan(numel(pch),1); %q per channel, NaN = not tested
            hasp = find(isnan(pch) == false);
            [ps,ord] = sort(pch(hasp));
            nt = numel(ps);
            qs = ps .* nt ./ (1:nt)';
            qs = flipud(cummin(flipud(qs))); %make them non-decreasing from the top
            qch(hasp(ord)) = min(qs,1);
            a_nsig(a) = sum(qch < 0.05);

            %this animal's probe layout
            pf = fullfile(probeinfo_dir, animal + "-ProbeInfo.mat");
            if isfile(pf) == false
                chanmap_notes(end+1) = animal + ": no ProbeInfo file, tile left empty";
                title(ax,animal)
                continue
            end
            pinfo = load(pf,'ProbeInfo');
            region = tableres.Region(find(tableres.Animal_Name == animal,1)); %which probe = the region in the table
            prb = find(strcmp(string(pinfo.ProbeInfo.Areas),region));
            if isempty(prb) %the Region value can also be a probe number
                prb = str2double(region);
            end
            pmap = pinfo.ProbeInfo.ProbeMaps{prb}; %[channel rows x shanks] of Channel_IDs

            M = nan(size(pmap)); %[channel rows x shanks], NaN = no value (black or grey)
            S = false(size(pmap)); %significant channels
            chans = unit_chan(inA);
            dvals = sd_LW(inA) - sd_L(inA);
            for ch = 1:numel(chans)
                loc = find(pmap == chans(ch));
                if isempty(loc)
                    chanmap_notes(end+1) = animal + ": channel " + chans(ch) + " not in ProbeMaps{" + prb + "}";
                else
                    M(loc) = dvals(ch);
                    S(loc) = qch(ch) < 0.05;
                end
            end

            %grey: channels that are in the csv for this pair but have no value (<2 trials left in L or LW)
            G = false(size(pmap));
            for loc = 1:numel(pmap)
                if isnan(M(loc)) && ismember(animal + "_ch" + string(pmap(loc)) + "_" + pairs(j), prefilter_keys)
                    G(loc) = true;
                end
            end
            a_nfew(a) = sum(G(:));

            h = imagesc(ax,M);
            set(h,'AlphaData',isnan(M) == false) %NaN cells see-through...
            set(ax,'Color','k') %...onto a black background (removed / not in the csv)
            colormap(ax,bwr)
            clim(ax,[-dmax dmax])
            hold(ax,'on')
            image(ax,repmat(reshape(greyc,1,1,3),size(G)),'AlphaData',double(G)) %grey squares on top, only where G
            [sr,sc] = find(S);
            text(ax,sc,sr,'*','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14,'FontWeight','bold')
            hold(ax,'off')
            axis(ax,'image')
            xticks(ax,1:size(pmap,2)); yticks(ax,1:size(pmap,1));
            xlabel(ax,'shank'); ylabel(ax,'channel (top = 1)')
            title(ax,{char(animal), sprintf('LW tighter %d/%d, p = %.3g', a_tighter(a), a_n(a), a_p(a))},'Interpreter','none')
        end
        cb = colorbar(ax);
        cb.Layout.Tile = 'east';
        cb.Label.String = 'SD_{LW} - SD_{L} (ms)';

        %legend for the grey and black squares (invisible patches just to make the legend entries)
        hold(ax,'on')
        lg1 = patch(ax,NaN,NaN,greyc,'EdgeColor','none');
        lg2 = patch(ax,NaN,NaN,'k','EdgeColor','none');
        hold(ax,'off')
        lgd = legend([lg1 lg2],{'Not Enough Trials','Removed / not in csv'},'Orientation','horizontal');
        lgd.Layout.Tile = 'south';
        lgd.Color = 'w'; %otherwise it copies the black axes background and the text disappears
        lgd.TextColor = 'k';
        lgd.FontSize = 11;

        animal_res{j} = table(pair_animals,a_n,a_medL,a_medLW,a_tighter,a_p,a_nsig,a_nfew, ...
            'VariableNames',["Animal_Name","n_channels","median_SD_L","median_SD_LW","n_channels_LWtighter", ...
            "p_SD_signrank","n_channels_sigBF_FDR","n_channels_NotEnoughTrials"]);
    end

    %all pairs in one table for the report
    stats_animal = table();
    for j = 1:numel(pairs)
        animal_res{j}.Pair = repmat(pairs(j),height(animal_res{j}),1);
        stats_animal = [stats_animal; animal_res{j}(:,["Pair" "Animal_Name" "n_channels" "median_SD_L" "median_SD_LW" ...
            "n_channels_LWtighter" "p_SD_signrank" "n_channels_sigBF_FDR" "n_channels_NotEnoughTrials"])];
    end
end


%% report

disp(['fraction of trials kept by the response filter: ' num2str(frac_kept_respfilter,3)])
disp(['fraction of those kept by the P1_LFPRMSsum filter: ' num2str(frac_kept_rmsfilter,3)])
disp(unitres)
disp(stats_paired)
disp(stats_pooled)
if chanmaps == 1
    disp(stats_animal)
    disp(chanmap_notes)
end
%NoPeakFrac = fraction of trials with no P1 peak at all. Also note the csv peak has no amplitude
%threshold, so channels with no real response still give a (random) peak time and a big SD.
