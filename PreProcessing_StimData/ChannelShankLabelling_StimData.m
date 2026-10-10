%% ChannelShankLabelling_StimData.m
%
% Description: Identifies which channels/shanks on a probe are most active
%   during stimulus-driven activity (according to some set of stimuli chosen
%   by the user, e.g. light/visual stim), so each shank or channel can be
%   labelled by these values for later location referencing (spontaneous
%   activity is scored separately in ChannelShankLabelling_SponData.m).
%   Loads stim_lfp_stimchunks from user-picked processed LFP files (trials
%   from all files are pooled), computes each trial/channel's moving-RMS
%   envelope (200-sample window), subtracts its mean pre-stim RMS (samples
%   3000-4999), and sums over samples 5000-8000. Channel score = mean of that
%   over trials; shank score = mean of the channel scores within each
%   ProbeMaps column. The channel and shank scores are added to ProbeInfo
%   (ProbeInfo.ChanLabels / ProbeInfo.ShankLabels) under a label name the
%   user types in when prompted (e.g. 'whisker'), and ProbeInfo is saved back
%   to its file. Finally plots shank and channel scores per probe.
%
% Inputs:
%   animal (string) - animal ID, e.g. '20260423-p12'
%   probes (int vector) - indices into ProbeInfo.ProbeIds/ProbeMaps to analyze
%   save_directory (string) - folder holding <animal>-ProbeInfo.mat (i.e.
%     the ProbeInfo folder itself)
%   ProbeInfo (struct) - loaded from <save_directory>\<animal>-ProbeInfo.mat.
%     Fields used:
%       .ProbeMaps (cell, 1 x probenum) - [rows x shanks] raw channel IDs;
%         each column is one shank
%       .ProbeIds (cell, 1 x probenum) - raw channel IDs per probe
%         (reshape(ProbeMaps{p},[],1))
%       .Ch_Remove (cell, 1 x probenum) - raw channel IDs to exclude per probe
%   stimdata_directory (cell of strings, picked via uipickfiles) - processed
%     LFP .mat files (from OpenEphys_BaseAnalysis), each containing:
%       stim_lfp_stimchunks (double, trials x samples x channels) - peri-
%         stimulus LFP, 1 kHz, microvolts; 3rd dim indexed by raw channel ID
%   labelname (string, typed in via inputdlg prompt) - struct field name
%     for these scores in ProbeInfo.ChanLabels / ShankLabels, e.g. 'whisker'
%     (default 'StimActivity'). Made a valid field name if needed; asks
%     before overwriting an existing label of the same name. Cancel/empty
%     = labels not added and ProbeInfo not saved.
%
% Outputs:
%   ProbeInfo (struct) - SAVED, overwriting <save_directory>\<animal>-ProbeInfo.mat.
%     Adds/updates (other existing label fields, e.g. SponActivity, are kept):
%       .ChanLabels{p}.(labelname) (double, rows x shanks, ProbeMaps layout,
%         e.g. 8x8) - channel's chActivityScore_Stim; NaN for Ch_Remove channels
%       .ShankLabels{p}.(labelname) (double, 1 x shanks, e.g. 1x8) -
%         shnkActivityScore_Stim; NaN for shanks with no scored channels
%       (ChanLabels / ShankLabels: cell, 1 x probenum, indexed like
%       ProbeMaps, each a struct of labels - see ChannelShankLabelling_SponData.m)
%   Figures (one per probe, not saved) - top: shnkActivityScore_Stim line
%     plot across shanks; bottom: heatmap of chActivityScore_Stim in the
%     ProbeMaps layout (rows x shanks), Ch_Remove channels in black ('X').
%   Workspace only (not saved):
%   chanals (int vector, Nchan x 1) - raw channel IDs analyzed: selected
%     probes' ProbeIds minus Ch_Remove, same order as oshkosh_detect_funct's
%     .chanals
%   stimdat (double, trials x samples x Nchan) - pooled LFP of those channels
%   stimdat_RMS (double, trials x Nchan) - baselined, summed post-stim
%     moving-RMS per trial and channel (a.u., uV*samples)
%   chActivityScore_Stim (cell, 1 x length(probes)) - each 1 x (number of
%     probes(p)'s channels in chanals) double: trial-mean of stimdat_RMS,
%     in chanals order
%   shnkActivityScore_Stim (cell, 1 x length(probes)) - each 1 x Nshanks
%     double: mean channel score per shank (ProbeMaps column) of probes(p)
%
% Dependencies: OpenEphys_BaseAnalysis(_Bundled).m (builds ProbeInfo and the
%   LFP files with stim_lfp_stimchunks), OpenEphys_editProbeInfo_ChRemove.m
%   (Ch_Remove), uipickfiles.
%
% NOTE (inferred data contracts, not verified):
%   - Hardcodes 1 kHz and stimulus onset at ~sample 5000, i.e. pre_second = 5
%     in OpenEphys_BaseAnalysis (onset is actually sample pre+1 = 5001), and
%     epochs at least 8000 samples long.
%   - Every picked LFP file must have the same epoch length and channel
%     layout (trials are concatenated along dim 1).
%   - ProbeInfo.Ch_Remove must have an entry for every probe in probes.
%   - Stim scores are baselined, so they can be negative (post-stim RMS below
%     pre-stim); the heatmap color scale just spans min to max.


%%

animal = '20260821-p10';
probes = [1]; %which probes in the recording to look at? 


save_directory = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%fq_range = [2 80]; %frequency band to check for oscillations in 

load(fullfile([save_directory '\' animal '-ProbeInfo.mat']))


%% Measure strongest stimulated shank (by light)

%get lfp files to use
stimdata_directory = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\LFP','Prompt', 'Choose processed LFP .mat files (containing stim_lfp_stimchunks) you want to use');

%channels to use: whole selected probes, minus removed channels (same as oshkosh_detect_funct)
chanals_raw = vertcat(ProbeInfo.ProbeIds{probes});
baddies = vertcat(ProbeInfo.Ch_Remove{probes});
chanals = chanals_raw(~ismember(chanals_raw, baddies));

%load data (from relevant channels)
stimdat = [];
for file = 1:length(stimdata_directory)
    d = load(stimdata_directory{file},'stim_lfp_stimchunks');
    stimdat = cat(1,stimdat, d.stim_lfp_stimchunks(:,:,chanals));
end
clear d

%get RMSsums, baselined by prestim activity
stimdat_RMS = zeros(size(stimdat,1),size(stimdat,3));
for tr = 1:size(stimdat,1)
    rms = sqrt(movmean(squeeze(stimdat(tr,:,:)).^2,200));
    baseline = mean(rms(3000:4999,:),1);
    baselined_rms = rms - repmat(baseline,size(rms,1),1);
    stimdat_RMS(tr,:) = sum(baselined_rms(5000:8000,:));
end

%get channel and shank specificy values
chActivityScore_Stim = cell(1,length(probes));
shnkActivityScore_Stim = cell(1,length(probes));
for prbi = 1:length(probes)
    prb = probes(prbi);
    probemask = ismember(chanals,ProbeInfo.ProbeMaps{prb}); %this probe's channels within chanals
    prb_chanals = chanals(probemask); %raw channel IDs, same order as chActivityScore_Stim{prbi}
    chActivityScore_Stim{prbi} = mean(stimdat_RMS(:,probemask),1);
    for shnk = 1:size(ProbeInfo.ProbeMaps{prb},2)
      shnkActivityScore_Stim{prbi}(shnk) = mean(chActivityScore_Stim{prbi}(ismember(prb_chanals,ProbeInfo.ProbeMaps{prb}(:,shnk))));
   end
end

%% Add activity labels to ProbeInfo and save
% Same structure as ChannelShankLabelling_SponData.m: ProbeInfo.ChanLabels{p} /
% ProbeInfo.ShankLabels{p} are one struct per probe (indexed like ProbeMaps),
% one field per label type. The user names the field for these stimulus
% scores (e.g. 'whisker'); other label fields already in ProbeInfo are kept.
labelinput = inputdlg('Name for this stimulus label field (e.g. whisker):','Label name',[1 50],{'StimActivity'});
savelabels = ~isempty(labelinput) && ~isempty(strtrim(labelinput{1}));
if ~savelabels
    warning('No label name given - labels NOT added to ProbeInfo or saved.');
else
    labelname = matlab.lang.makeValidName(strtrim(labelinput{1})); %must be a valid struct field name
    if ~strcmp(labelname,strtrim(labelinput{1}))
        warning(['Label name changed to "' labelname '" to make it a valid field name.']);
    end

    if ~isfield(ProbeInfo,'ChanLabels')
        ProbeInfo.ChanLabels = cell(1,numel(ProbeInfo.ProbeMaps));
    end
    if ~isfield(ProbeInfo,'ShankLabels')
        ProbeInfo.ShankLabels = cell(1,numel(ProbeInfo.ProbeMaps));
    end

    %ask before overwriting a label of the same name on any selected probe
    labelexists = false;
    for prbi = 1:length(probes)
        prb = probes(prbi);
        labelexists = labelexists || (numel(ProbeInfo.ChanLabels) >= prb && isfield(ProbeInfo.ChanLabels{prb},labelname)) ...
            || (numel(ProbeInfo.ShankLabels) >= prb && isfield(ProbeInfo.ShankLabels{prb},labelname));
    end
    if labelexists
        answer = questdlg(['Label "' labelname '" already exists in ProbeInfo. Overwrite it?'],'Label exists','Overwrite','Cancel','Cancel');
        savelabels = strcmp(answer,'Overwrite');
        if ~savelabels
            warning(['Label "' labelname '" NOT overwritten - ProbeInfo not saved.']);
        end
    end
end

if savelabels
    for prbi = 1:length(probes)
        prb = probes(prbi);
        pmap = ProbeInfo.ProbeMaps{prb}; %[rows x shanks] raw channel IDs
        prb_chanals = chanals(ismember(chanals,pmap)); %raw channel IDs, same order as chActivityScore_Stim{prbi}

        %channel scores rearranged into the ProbeMaps layout (NaN = no score, e.g. Ch_Remove)
        chLabelMap = NaN(size(pmap));
        [isScored,chIdx] = ismember(pmap,prb_chanals);
        chLabelMap(isScored) = chActivityScore_Stim{prbi}(chIdx(isScored));
        chLabelMap(ismember(pmap,ProbeInfo.Ch_Remove{prb})) = NaN;

        if numel(ProbeInfo.ChanLabels) < prb || isempty(ProbeInfo.ChanLabels{prb})
            ProbeInfo.ChanLabels{prb} = struct();
        end
        if numel(ProbeInfo.ShankLabels) < prb || isempty(ProbeInfo.ShankLabels{prb})
            ProbeInfo.ShankLabels{prb} = struct();
        end
        ProbeInfo.ChanLabels{prb}.(labelname) = chLabelMap;
        ProbeInfo.ShankLabels{prb}.(labelname) = shnkActivityScore_Stim{prbi}; %shanks with no scored channels are already NaN (mean of empty)
    end

    save(fullfile([save_directory '\' animal '-ProbeInfo.mat']),'ProbeInfo')
    disp(['Saved label "' labelname '" to ProbeInfo.ChanLabels / ProbeInfo.ShankLabels']);
end

%% Plot channel and shank activity scores
% One figure per probe: top = shank score per shank (mean of that shank's
% channel scores, excluding Ch_Remove channels), bottom = heatmap of
% channel scores in the same layout as ProbeMaps (columns = shanks, so each
% shank point sits directly above its heatmap column). Ch_Remove channels
% are black and marked 'X'; other cells are labelled with their channel ID.
seqmap = interp1([0 0.5 1],[0.95 0.98 0.93; 0.55 0.82 0.40; 0.18 0.60 0.22],linspace(0,1,256)); %single-hue green colormap, pale mint = low, clear green = high (never reaches black)
for prbi = 1:length(probes)
    prb = probes(prbi);
    pmap = ProbeInfo.ProbeMaps{prb}; %[rows x shanks] raw channel IDs
    [nrow,nshnk] = size(pmap);
    removed = ismember(pmap,ProbeInfo.Ch_Remove{prb});
    prb_chanals = chanals(ismember(chanals,pmap)); %raw channel IDs, same order as chActivityScore_Stim{prbi}

    %channel scores rearranged into the ProbeMaps layout (NaN = no score)
    chScoreMap = NaN(nrow,nshnk);
    [isScored,chIdx] = ismember(pmap,prb_chanals);
    chScoreMap(isScored) = chActivityScore_Stim{prbi}(chIdx(isScored));
    chScoreMap(removed) = NaN;

    figure('Name',[animal ' probe ' num2str(prb) ' stimulus-driven activity'],'Color','w');
    tl = tiledlayout(4,1,'TileSpacing','compact');
    title(tl,[animal ' - probe ' num2str(prb) ': stimulus-driven activity'],'Interpreter','none');

    %top: shank scores (mean of the channel scores in each shank / heatmap column)
    ax1 = nexttile(tl,1);
    plot(ax1,1:nshnk,shnkActivityScore_Stim{prbi},'-o','LineWidth',2,'MarkerSize',8,'Color',seqmap(end,:),'MarkerFaceColor',seqmap(end,:));
    xlim(ax1,[0.5 nshnk+0.5]);
    xticks(ax1,1:nshnk); xticklabels(ax1,{});
    ylabel(ax1,{'Shank score','(mean of channels)'});
    grid(ax1,'on'); ax1.GridAlpha = 0.15; box(ax1,'off');

    %bottom: channel heatmap (NaN cells are transparent over a black axes background)
    ax2 = nexttile(tl,2,[3 1]);
    imagesc(ax2,1:nshnk,1:nrow,chScoreMap,'AlphaData',~isnan(chScoreMap));
    ax2.Color = 'k';
    colormap(ax2,seqmap);
    cb = colorbar(ax2);
    cb.Layout.Tile = 'east'; %keep colorbar outside both tiles so the shank columns stay aligned
    cb.Label.String = 'Channel score (baselined post-stim RMS sum, a.u.)';
    xticks(ax2,1:nshnk); yticks(ax2,1:nrow);
    xlabel(ax2,'Shank (ProbeMaps column)'); ylabel(ax2,'ProbeMaps row');
    linkaxes([ax1 ax2],'x');

    %cell labels: channel ID, or 'X' for removed channels; black text on the green cells, white on black (removed) cells
    for r = 1:nrow
        for c = 1:nshnk
            if removed(r,c)
                txt = 'X';
            else
                txt = num2str(pmap(r,c));
            end
            if isnan(chScoreMap(r,c))
                txtcol = 'w';
            else
                txtcol = 'k';
            end
            text(ax2,c,r,txt,'HorizontalAlignment','center','Color',txtcol,'FontSize',8);
        end
    end
end




