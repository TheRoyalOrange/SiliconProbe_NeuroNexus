%% ChannelShankLabelling_SponData.m
%
% Description: Identifies which channels/shanks on a probe are most active
%   during spontaneous activity, so each shank or channel can be labelled by
%   these values for later location referencing (stimulus-driven activity is
%   scored separately in ChannelShankLabelling_StimData.m). Detects
%   oscillatory probe events (prbEvent) with oshkosh_detect_funct, then for
%   every prbEvent NOT overlapping a stimulus sums each channel's moving-RMS
%   envelope (200-sample window) over the event window. Shank score = mean
%   of the channel scores within each ProbeMaps column. The channel and shank
%   scores are added to ProbeInfo as the SponActivity label (in
%   ProbeInfo.ChanLabels / ProbeInfo.ShankLabels) and ProbeInfo is saved back
%   to its file.
%
% Inputs:
%   animal (string) - animal ID, e.g. '20260423-p12'
%   save_directory (string) - folder holding <animal>-ProbeInfo.mat (i.e.
%     the ProbeInfo folder itself, e.g. '...\Processed Silicon Probe Data\ProbeInfo')
%   fq_range (double, 1x2, Hz) - bandpass used for oscillation detection
%   probes (int vector) - indices into ProbeInfo.ProbeIds/ProbeMaps to analyze
%   reorder (0/1) - pre-06/2024 channel remap flag, passed to
%     oshkosh_detect_funct (if 1, a *_channel_re_index file must be loaded
%     so newchannelidx is in scope)
%   ProbeInfo (struct) - loaded from <save_directory>\<animal>-ProbeInfo.mat.
%     Fields used here:
%       .ProbeMaps (cell, 1 x probenum) - [rows x shanks] raw channel IDs;
%         each column is one shank
%       .Ch_Remove (cell, 1 x probenum) - raw channel IDs removed per probe;
%         drawn black in the heatmap (also used inside oshkosh_detect_funct)
%     (.ProbeIds, .TTLch are used inside oshkosh_detect_funct)
%   spondata_directory (cell of strings, picked via uipickfiles) - raw
%     OpenEphys recording folders used for oscillation detection
%
% Outputs:
%   ProbeInfo (struct) - SAVED, overwriting <save_directory>\<animal>-ProbeInfo.mat.
%     Adds/updates (other existing label fields, e.g. StimActivity, are kept):
%       .ChanLabels (cell, 1 x probenum, indexed like ProbeMaps) - each a struct
%         of per-channel labels; each label field has the ProbeMaps layout
%         (rows x shanks, e.g. 8x8), numeric or cell array of strings:
%         .ChanLabels{p}.SponActivity (double, rows x shanks) - channel's
%           chActivityScore_Spon; NaN for Ch_Remove channels
%       .ShankLabels (cell, 1 x probenum) - each a struct of per-shank labels,
%         each label field 1 x shanks (e.g. 1x8), numeric or cell of strings:
%         .ShankLabels{p}.SponActivity (double, 1 x shanks) -
%           shnkActivityScore_Spon; NaN for shanks with no scored channels
%       Probes not in `probes` are left untouched (empty struct if new).
%   Workspace only (not saved):
%   oscillations, snippets, prbsnippets (structs) - from oshkosh_detect_funct
%     (see its header). Used here:
%       prbsnippets.chanals (int vector, Nchan) - raw channel IDs analyzed
%         (selected probes' channels minus Ch_Remove)
%       prbsnippets.event_data{p}{event} (double, Nchan x window samples) -
%         filtered LFP (uV, 1 kHz) of ALL analyzed channels for each prbEvent
%       oscillations.prbEvent_during_stim{p} (logical/0-1 per prbEvent)
%   chActivityScore_Spon (cell, 1 x length(probes)) - each Nchan x 1 double:
%     summed moving-RMS over probes(p)'s non-stim prbEvents, for ALL analyzed
%     channels, in prbsnippets.chanals order (a.u., uV*samples)
%   shnkActivityScore_Spon (cell, 1 x length(probes)) - each 1 x Nshanks
%     double: mean channel score per shank (ProbeMaps column) of probes(p)
%   Figures (one per probe, not saved) - top: shnkActivityScore_Spon line
%     plot across shanks; bottom: heatmap of chActivityScore_Spon in the
%     ProbeMaps layout (rows x shanks), Ch_Remove channels in black ('X').
%
% Dependencies: OpenEphys_BaseAnalysis(_Bundled).m (builds ProbeInfo),
%   OpenEphys_editProbeInfo_ChRemove.m (Ch_Remove), oshkosh_detect_funct.m,
%   oscillation_plotEvents.m, uipickfiles, Signal Processing Toolbox.
%
% NOTE (inferred data contracts, not verified):
%   - ProbeInfo.ProbeMaps values are raw channel IDs, the same IDs as
%     prbsnippets.chanals (both come from ProbeMaps via
%     ProbeIds = reshape(ProbeMaps,[],1) in OpenEphys_BaseAnalysis).
%   - The probe's Ch_Remove channels are the only ProbeMaps channels missing
%     from prbsnippets.chanals; any other missing channel would also plot
%     black, but labelled with its ID instead of 'X'.


%%
animal = '20260821-p10';

save_directory = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

fq_range = [2 80]; %frequency band to check for oscillations in 
probes = [1]; %which probes in the recording to look at? 


reorder = 0;

%IN CASE REORDER IS 0: 
%load 'pre202406_64_32_channel_re_index'
%load 'pre202406_64_32_channel_re_index_v2'
%load '202408_64_64_channel_re_index_v5.mat'
%load 20240708_64_32_channel_re_index.mat



%% NEEDS INPUT FIRST: Clear the board and set the directory (change the directory path)
load(fullfile([save_directory '\'  animal '-ProbeInfo.mat']))
load_directory = 'E:\Roy\Silicon Probe Raw Data'; %don't change unless you're not Roy and have your data elsewhere

spondata_directory = uipickfiles('FilterSpec',load_directory,'Prompt', 'Choose Folders Containing Recordings you want to use [IMMEDIATE SUBFOLDERS SHOULD BE -Record Node ###-]'); %loads file names'E:\Roy\Ungrouped Whiskers\20251011-p11\20251011-p11-light4l22_40'; %for example

 % data_directory = {  'E:\Roy\Silicon Probe Raw Data\20260423-p12\20260423-p12-spon';
 %                     'E:\Roy\Silicon Probe Raw Data\20260423-p12\20260423-p12-whisker';
 %                     %'E:\Roy\Silicon Probe Raw Data\20260423-p12\20260423-p12-light22_40';
%                    'E:\Roy\Silicon Probe Raw Data\20260226-p12\20260226-p12-spon';
%                    'E:\Roy\Silicon Probe Raw Data\20260226-p12\20260226-p12-whisker';
%                    'E:\Roy\Silicon Probe Raw Data\20260226-p12\20260226-p12-light22_40';
                    %'E:\Roy\Ungrouped Whiskers\20251011-p11\20251011-p11-spon1';
                    %'E:\Roy\Ungrouped Whiskers\20251011-p11\20251011-p11-spon2';
                    %'E:\Roy\Ungrouped Whiskers\20240325-p6-raw\20240325-p6-spon';
                    %'E:\Roy\Ungrouped Whiskers\20240325-p6-raw\20240325-p6-whisker';
                    %'E:\Roy\Ungrouped Whiskers\20250308-p9\20250308-p9-light22_72';
                    %'E:\Roy\Ungrouped Whiskers\20250308-p9\20250308-p9-light22_72-n-whisker'
                       % 'E:\Roy\Silicon Probe Raw Data\20260326-p5\20260326-p5-spon';
                       % 'E:\Roy\Silicon Probe Raw Data\20260326-p5\20260326-p5-spon2';
                       % 'E:\Roy\Silicon Probe Raw Data\20260326-p5\20260326-p5-whisker';
                       % 'E:\Roy\Silicon Probe Raw Data\20260326-p5\20260326-p5-whisker10Hz';
              %     }; %for example

[oscillations, snippets, prbsnippets] = oshkosh_detect_funct(animal, spondata_directory, probes, fq_range, save_directory, reorder, ProbeInfo)
oscillation_plotEvents(ProbeInfo,oscillations, prbsnippets)

%% Measure total spontaneous activity of each channel

%compute RMSsum of each channel for each probe event
chActivityScore_Spon = cell(1,length(probes));
shnkActivityScore_Spon = cell(1,length(probes));
for prb = 1:length(probes)
    chRMSsum_total = zeros(length(prbsnippets.chanals),1);
    for snp =  1:length(prbsnippets.event_data{prb})
        if oscillations.prbEvent_during_stim{prb}(snp) == 0
        dat_snp = prbsnippets.event_data{prb}{snp};
        snp_sums = [];
        for ch = 1:size(dat_snp,1)
        snp_sums(ch) = sum(sqrt(movmean(dat_snp(ch,:).^2,200)),2);
        end
        chRMSsum_total = chRMSsum_total + snp_sums';
        else
        end
    end
    chActivityScore_Spon{prb} = chRMSsum_total;
end

for prb = 1:length(probes)
   %prb is a position in probes, so index ProbeMaps by the actual probe number probes(prb)
   for shnk = 1:size(ProbeInfo.ProbeMaps{probes(prb)},2)
      shnkActivityScore_Spon{prb}(shnk) = mean(chActivityScore_Spon{prb}(ismember(prbsnippets.chanals,ProbeInfo.ProbeMaps{probes(prb)}(:,shnk))))
   end
end

%% Add activity labels to ProbeInfo and save
% ProbeInfo.ChanLabels{p} / ProbeInfo.ShankLabels{p}: one struct per probe
% (indexed like ProbeMaps), one field per label type. Only the SponActivity
% field is (re)written here; other label fields already in ProbeInfo are kept.
if ~isfield(ProbeInfo,'ChanLabels')
    ProbeInfo.ChanLabels = cell(1,numel(ProbeInfo.ProbeMaps));
end
if ~isfield(ProbeInfo,'ShankLabels')
    ProbeInfo.ShankLabels = cell(1,numel(ProbeInfo.ProbeMaps));
end
for prb = 1:length(probes)
    pmap = ProbeInfo.ProbeMaps{probes(prb)}; %[rows x shanks] raw channel IDs

    %channel scores rearranged into the ProbeMaps layout (NaN = no score, e.g. Ch_Remove)
    chLabelMap = NaN(size(pmap));
    [isScored,chIdx] = ismember(pmap,prbsnippets.chanals);
    chLabelMap(isScored) = chActivityScore_Spon{prb}(chIdx(isScored));
    chLabelMap(ismember(pmap,ProbeInfo.Ch_Remove{probes(prb)})) = NaN;

    if numel(ProbeInfo.ChanLabels) < probes(prb) || isempty(ProbeInfo.ChanLabels{probes(prb)})
        ProbeInfo.ChanLabels{probes(prb)} = struct();
    end
    if numel(ProbeInfo.ShankLabels) < probes(prb) || isempty(ProbeInfo.ShankLabels{probes(prb)})
        ProbeInfo.ShankLabels{probes(prb)} = struct();
    end
    ProbeInfo.ChanLabels{probes(prb)}.SponActivity = chLabelMap;
    ProbeInfo.ShankLabels{probes(prb)}.SponActivity = shnkActivityScore_Spon{prb}; %shanks with no scored channels are already NaN (mean of empty)
end

save(fullfile([save_directory '\' animal '-ProbeInfo.mat']),'ProbeInfo')

%% Plot channel and shank activity scores
% One figure per probe: top = shank score per shank (mean of that shank's
% channel scores, excluding Ch_Remove channels), bottom = heatmap of
% channel scores in the same layout as ProbeMaps (columns = shanks, so each
% shank point sits directly above its heatmap column). Ch_Remove channels
% are black and marked 'X'; other cells are labelled with their channel ID.
seqmap = interp1([0 0.5 1],[1.00 0.96 0.90; 1.00 0.62 0.15; 0.92 0.38 0.00],linspace(0,1,256)); %single-hue orange colormap, pale cream = low, vivid orange = high (red kept high so it doesn't go brown; never reaches black)
for prb = 1:length(probes)
    pmap = ProbeInfo.ProbeMaps{probes(prb)}; %[rows x shanks] raw channel IDs
    [nrow,nshnk] = size(pmap);
    removed = ismember(pmap,ProbeInfo.Ch_Remove{probes(prb)});

    %channel scores rearranged into the ProbeMaps layout (NaN = no score)
    chScoreMap = NaN(nrow,nshnk);
    [isScored,chIdx] = ismember(pmap,prbsnippets.chanals);
    chScoreMap(isScored) = chActivityScore_Spon{prb}(chIdx(isScored));
    chScoreMap(removed) = NaN;

    figure('Name',[animal ' probe ' num2str(probes(prb)) ' spontaneous activity'],'Color','w');
    tl = tiledlayout(4,1,'TileSpacing','compact');
    title(tl,[animal ' - probe ' num2str(probes(prb)) ': spontaneous activity'],'Interpreter','none');

    %top: shank scores (mean of the channel scores in each shank / heatmap column)
    ax1 = nexttile(tl,1);
    plot(ax1,1:nshnk,shnkActivityScore_Spon{prb},'-o','LineWidth',2,'MarkerSize',8,'Color',seqmap(end,:),'MarkerFaceColor',seqmap(end,:));
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
    cb.Label.String = 'Channel score (summed RMS, a.u.)';
    xticks(ax2,1:nshnk); yticks(ax2,1:nrow);
    xlabel(ax2,'Shank (ProbeMaps column)'); ylabel(ax2,'ProbeMaps row');
    linkaxes([ax1 ax2],'x');

    %cell labels: channel ID, or 'X' for removed channels; black text on the orange cells, white on black (removed) cells
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

