function [oscillations, snippets, probeSnippets] = oshkosh_detect_funct(animal, data_directory, probes, fq_range, save_directory, reorder, ProbeInfo)
% oshkosh_detect_funct.m
%
% Description: Detects oscillatory bouts of spontaneous LFP activity, per
%   channel, from one or more OpenEphys recordings for a single animal,
%   restricted to one or more whole probes. For each recording folder
%   listed, loads raw data, downsamples to 1kHz and converts to
%   microvolts, bandpasses each channel into a single band (bounds set by
%   fq_range), computes a moving-RMS envelope, and detects
%   threshold-crossing events (RMS > median + 2*MAD, peak > median +
%   3*MAD) per channel
%   using median-absolute-deviation statistics. Adjacent events closer
%   than 500ms are merged; events shorter than 300ms or too close to the
%   recording edges are discarded. Flags whether each event overlaps a
%   stimulus TTL pulse (if present). Adapted from Opatz lab / Buszaki lab
%   oscillation-detection code (see getOscillations.m,
%   https://github.com/OpatzLab/HanganuOpatzToolbox/blob/master/Detect%20Oscillations/getOscillations.m).
%   Not currently part of the PreProcessing_StimData -> Functional_Connectivity
%   / Summarize_MultiAnimal pipeline order; standalone spontaneous-activity
%   analysis. Callable version of (renamed from)
%   OpenEphys_Detection_of_Oscillatory_Events.m - channel selection is now
%   probe-based via ProbeInfo instead of a hardcoded total-channel-count
%   guess; a separate function for more flexible/manual channel selection
%   is planned later. Detects events at two levels: a channel event
%   (chEvent) is what's described above, detected independently per
%   channel; a probe event (prbEvent) is the union of a probe's channels'
%   chEvent windows - a time window where at least one channel on that
%   probe has an ongoing chEvent, merging overlapping/close-together
%   chEvent windows (gap < 300ms) across all of that probe's channels.
%   NOTE: 300ms is the user-specified prbEvent merge gap; chEvent's own
%   merge gap (for adjacent threshold-crossings on the SAME channel) is
%   actually 500ms elsewhere in this function (300ms there is instead the
%   minimum-duration cutoff for keeping a chEvent) - flagging this
%   difference in case 500ms was actually intended for prbEvent too.
%   Returns three outputs, oscillations, snippets, and probeSnippets (see
%   Outputs) - does not save any of them to disk itself; that's left to
%   the caller.
%
% Inputs:
%   animal (string) - REQUIRED. Animal identifier. Used (unless ProbeInfo
%     is passed in directly) to locate this animal's ProbeInfo file.
%   data_directory (cell array of strings, Nx1) - REQUIRED. List of
%     OpenEphys recording folder paths to process for this animal, one per
%     file/condition (e.g. spon, whisker).
%   probes (int vector) - REQUIRED. Which probe(s) to analyze, as indices
%     into ProbeInfo.ProbeIds / ProbeInfo.Ch_Remove (e.g. [1] for one
%     probe, [1 2] for two). All channels belonging to these probes are
%     used; there is no per-channel override in this function.
%   fq_range (numeric, 1x2, [low high] in Hz) - REQUIRED. Passband edges
%     used to bandpass every channel's data (data_filt) before computing
%     the RMS envelope and detecting events. Replaces the old hardcoded
%     dual low/high band split - there is now only one band.
%   save_directory (string) - REQUIRED. Base folder used to locate
%     ProbeInfo, from <save_directory>\ProbeInfo\<animal>-ProbeInfo.mat
%     (same convention as ConditionalTrialRemove_callable.m), unless
%     ProbeInfo is passed in directly. This function does NOT save
%     anything to disk itself (see Outputs) - callers who want to persist
%     oscillations/snippets/probeSnippets should save() them wherever they
%     like, e.g. save(fullfile(save_directory,[animal '_Oscillations']),
%     'oscillations') or ..., 'oscillations', 'snippets', 'probeSnippets'.
%   reorder (0 or 1, optional) - Whether to remap channel order, for
%     recordings from before the 06/2024 mapping change. Defaults to 0.
%     NOTE (pre-existing, not fixed here): when set to 1, this branch
%     references a variable newchannelidx that this function does not
%     define or load - the caller must have it in scope (e.g. via a prior
%     `load` of one of the '..._channel_re_index...' files) or this will
%     error.
%   ProbeInfo (struct, optional) - This animal's ProbeInfo struct. If
%     omitted or left empty ([]), loaded automatically from
%     <save_directory>\ProbeInfo\<animal>-ProbeInfo.mat. Fields used:
%       .ProbeIds (cell array, 1 x probenum) - channel IDs belonging to
%         each probe. NOTE (inferred data contract): this field's own
%         header comment in OpenEphys_BaseAnalysis.m claims local per-probe
%         numbering ("starts at 1"), but the code that builds it
%         (probeids{i} = reshape(probemaps{i},[],1)) just reshapes
%         ProbeMaps, so for probe 2+ it actually holds the same raw/offset
%         channel IDs as ProbeMaps (e.g. 65-128), not 1-N. This function
%         relies on the real (raw-ID) values, concatenating
%         ProbeInfo.ProbeIds(probes) and using them directly as row
%         indices into the raw data matrix.
%       .Ch_Remove (cell array, 1 x probenum) - per-probe raw channel IDs
%         to exclude (bad/noisy channels), same indexing as ProbeIds.
%       .TTLch (int) - raw-data row index of the stimulus TTL channel.
%   data (OpenEphys recording, loaded via Session/continuous .samples) -
%     raw multi-channel voltage traces at native sampling rate (int16),
%     downsampled to 1kHz and converted to microvolts (x0.1950).
%
% Outputs:
%   oscillations (struct) - event detection results at two levels,
%     channel (chEvent) and probe (prbEvent) - lightweight, no waveform
%     data, just event bookkeeping. Also carries two non-event fields:
%       .chanals (int vector) - raw channel IDs used (union of the
%         selected probes' ProbeIds, minus baddies)
%       .exclude_channels (int vector) - copy of baddies (raw channel IDs
%         excluded, from ProbeInfo.Ch_Remove for the selected probes)
%       .files (cell array of strings) - copy of data_directory
%       .probes (int vector) - copy of the probes input. NOTE: added
%         beyond what was strictly asked for the chEvent->prbEvent split -
%         without it there's no way to map prbEvent_*{p} back to which
%         real probe ID p refers to.
%     chEvent fields (cell array, 1 x Nchan, indexed by position in
%     oscillations.chanals - one detected event per channel, independent
%     of every other channel):
%       .chEvent_timestamps (cell) - [start stop] sample indices (1kHz)
%         of each detected event
%       .chEvent_ext_timestamps (cell) - event timestamps extended by 1s
%         on each side
%       .chEvent_neg_timestamps (cell) - inter-event (quiescent) window
%         timestamps
%       .chEvent_peaks (cell) - sample index of the negative peak within
%         each event
%       .chEvent_peakNormedPower (cell) - RMS amplitude at the event's peak
%       .chEvent_durations (cell) - event duration in samples
%       .chEvent_ext_durations (cell) - extended event duration in samples
%       .chEvent_file_index (cell) - index into data_directory identifying
%         which recording each event came from
%       .chEvent_index (cell) - per-file event ordinal number
%       .chEvent_during_stim (cell) - logical, whether the event overlaps
%         a detected stimulus TTL pulse
%     prbEvent fields (cell array, 1 x length(probes), indexed by position
%     in oscillations.probes - one detected event per probe, being the
%     union of that probe's own channels' chEvent windows, merged whenever
%     the gap between them is under 300ms - see Description note on this
%     constant): same shape/meaning as the corresponding chEvent_ field
%     above, one level up (a "channel" is a "probe" instead):
%       .prbEvent_timestamps, .prbEvent_ext_timestamps,
%       .prbEvent_neg_timestamps, .prbEvent_durations,
%       .prbEvent_ext_durations, .prbEvent_file_index, .prbEvent_index,
%       .prbEvent_during_stim - as above, one per probe instead of per
%         channel.
%       .prbEvent_peaks / .prbEvent_peakNormedPower - computed by
%         rescanning ALL of that probe's channels' data_filt/data_filt_rms
%         over the merged window (peak = position of the single
%         most-negative point across those channels; peakNormedPower =
%         the max RMS envelope value across those channels) - i.e. not
%         copied from any one contributing chEvent. NOTE: if this
%         rescan is too slow in practice, a documented (not yet
%         implemented) fallback is to instead copy these two values from
%         whichever contributing chEvent had the highest peakNormedPower.
%   snippets (struct, optional - only computed cost is cheap, but this can
%     still be a large variable in memory/on disk since it holds waveform
%     data) - per-chEvent, all-channel filtered-signal snippets, meant for
%     downstream analysis/visualization (e.g. paging through events and
%     plotting every channel at once). Deliberately holds only bounded
%     per-event windows, not the full continuous recording. Fields:
%       .chanals (int vector) - copy of oscillations.chanals
%       .files (cell array of strings) - copy of data_directory
%       .event_data (cell, 1 x Nchan, indexed by position in
%         oscillations.chanals) - event_data{ch} is itself a cell array,
%         one cell per event detected on that channel, in the same
%         order as oscillations.chEvent_timestamps{ch}'s rows. Each cell
%         holds a [numel(chanals) x window_length] double matrix: EVERY
%         analyzed channel's filtered signal (data_filt), sliced over
%         that event's extended window (matching
%         oscillations.chEvent_ext_timestamps{ch} for the same row).
%         window_length varies per event.
%   probeSnippets (struct, optional - same caveats as snippets) - same
%     idea as snippets, one level up: per-prbEvent, all-channel
%     filtered-signal snippets. Fields:
%       .chanals (int vector) - copy of oscillations.chanals
%       .probes (int vector) - copy of oscillations.probes
%       .files (cell array of strings) - copy of data_directory
%       .event_data (cell, 1 x length(probes), indexed by position in
%         oscillations.probes) - event_data{p} is itself a cell array,
%         one cell per prbEvent detected on that probe, in the same order
%         as oscillations.prbEvent_timestamps{p}'s rows. Each cell holds
%         a [numel(chanals) x window_length] double matrix covering EVERY
%         analyzed channel (not just that probe's own channels), sliced
%         over that prbEvent's extended window. window_length varies per
%         event.
%
% Dependencies: Requires the OpenEphys MATLAB tools (Session class) and
%   Signal Processing Toolbox (bandpass, movmean). Requires
%   <animal>-ProbeInfo.mat to already exist under save_directory (built by
%   OpenEphys_BaseAnalysis.m or a variant, with Ch_Remove populated via
%   OpenEphys_editProbeInfo_ChRemove.m), unless ProbeInfo is passed in
%   directly.

if nargin < 7 || isempty(ProbeInfo)
    probeinfo_dat = load(fullfile([save_directory '\ProbeInfo\' animal '-ProbeInfo.mat']), 'ProbeInfo');
    ProbeInfo = probeinfo_dat.ProbeInfo;
end

if nargin < 6 || isempty(reorder)
    reorder = 0;
end

%% setup things (nothing to change here)
%make the list of channels to use: whole probes, minus excluded channels
chan_probe_raw = [];
for p = 1:length(probes)
    chan_probe_raw = [chan_probe_raw; repmat(p, numel(ProbeInfo.ProbeIds{probes(p)}), 1)];
end
chanals_raw = vertcat(ProbeInfo.ProbeIds{probes});
baddies = vertcat(ProbeInfo.Ch_Remove{probes});
keep_mask = ~ismember(chanals_raw, baddies);
chanals = chanals_raw(keep_mask);
chan_probe = chan_probe_raw(keep_mask); %same length/order as chanals - which entry of probes each channel belongs to

stimCh = ProbeInfo.TTLch;

oscillations.chanals = chanals;
oscillations.exclude_channels = baddies;
oscillations.probes = probes;
oscillations.files = data_directory;

snippets.chanals = chanals;
snippets.files = data_directory;
snippets.event_data = cell(length(chanals),1);
for ch = 1:length(chanals)
    snippets.event_data{ch} = {};
end

probeSnippets.chanals = chanals;
probeSnippets.probes = probes;
probeSnippets.files = data_directory;
probeSnippets.event_data = cell(length(probes),1);
for p = 1:length(probes)
    probeSnippets.event_data{p} = {};
end

%also need to make variables for relative and absolute thresholds
%(min/max), among other things. Just check the link above.

%dont touch this
 oscillations.chEvent_timestamps = cell(length(chanals),1);
 oscillations.chEvent_ext_timestamps =  cell(length(chanals),1); %clear ext_osc
 oscillations.chEvent_neg_timestamps =  cell(length(chanals),1); %clear neg_osc
 oscillations.chEvent_peaks = cell(length(chanals),1); % peaktimes
 oscillations.chEvent_peakNormedPower = cell(length(chanals),1); % amplitudes
 oscillations.chEvent_durations = cell(length(chanals),1);
 oscillations.chEvent_ext_durations = cell(length(chanals),1);
 oscillations.chEvent_file_index = cell(length(chanals),1);
 oscillations.chEvent_index = cell(length(chanals),1);
 oscillations.chEvent_during_stim = cell(length(chanals),1);

 oscillations.prbEvent_timestamps = cell(length(probes),1);
 oscillations.prbEvent_ext_timestamps = cell(length(probes),1);
 oscillations.prbEvent_neg_timestamps = cell(length(probes),1);
 oscillations.prbEvent_peaks = cell(length(probes),1);
 oscillations.prbEvent_peakNormedPower = cell(length(probes),1);
 oscillations.prbEvent_durations = cell(length(probes),1);
 oscillations.prbEvent_ext_durations = cell(length(probes),1);
 oscillations.prbEvent_file_index = cell(length(probes),1);
 oscillations.prbEvent_index = cell(length(probes),1);
 oscillations.prbEvent_during_stim = cell(length(probes),1);




%%
for file = 1:size(data_directory,1)
    disp(['Running file ' num2str(file)])
data = [];
        session = Session(data_directory{file});



% Read in data

for i = 1:size(session.recordNodes{1,1}. recordings,2)
    datums = session.recordNodes{1,1}. recordings{1,i}.continuous('Acquisition_Board-100.acquisition_board').samples;
    %datums = session.recordNodes{1,1}. recordings{1,i}.continuous('Acquisition_Board-100.Rhythm Data').samples;
    datums = downsample(datums',30)'; fs = 1000; %to 1kH, each timepoint is 1 ms
    datums = double(datums).*0.1950; %converts int16 to uV
    datums_all{i} =datums;


end
data = horzcat(datums_all{:});



%% because you messed up the channel order for all pre 06-2024 animals
if reorder == 1
    data = data(newchannelidx,:);%; chanals(7)],:);
end




%% Find stim timepoints if they exist
%plot(data(stimCh,:)'); %check how the signal looks
normbineTTL = round(data(stimCh,:)'./max(data(stimCh,:)'));
if mean(normbineTTL) < .8
    stim_times = find(normbineTTL == 1);
else
    stim_times = [];
end
%% bandpass

data = data(chanals,:);

clear data_filt
for i = 1:length(chanals)
  data_filt(i,:) = bandpass(data(i,:),fq_range,fs, 'ImpulseResponse','iir');
  disp(i);
end

%for i = 1:length(chanals)
%  data_filt_2t8(i,:) =  bandpass(data(i,:),[2 8],fs,'ImpulseResponse','iir');
%  data_filt_8t15(i,:) =  bandpass(data(i,:),[8 15],fs,'ImpulseResponse','iir');
%  data_filt_15t30(i,:) = bandpass(data(i,:),[15,30],fs,'ImpulseResponse','iir');
%  data_filt_30t80(i,:) = bandpass(data(i,:),[30,80],fs,'ImpulseResponse','iir');
%  disp(i);
%end


%% plot the main V1 and S1 L4 channels

%figure()
%hold on
%plot(data_filt(chanals(2),:),Color= 'r')
%plot(data_filt(chanals(4),:), Color= 'k')
%legend({'V1 L4' 'S1 L4'});
%hold off


%% calculate rms values and plot histogram
clear data_filt_rms
for i = 1:length(chanals)
    %data_rms(i,:) = sqrt(movmean(data_filt(i,:).^2,200));
    %data_rms_2t8(i,:) = sqrt(movmean(data_filt_2t8(i,:).^2,200));
    %data_rms_8t15(i,:) = sqrt(movmean(data_filt_8t15(i,:).^2,200));
    %data_rms_15t30(i,:) = sqrt(movmean(data_filt_15t30(i,:).^2,200));
    %data_rms_30t80(i,:) = sqrt(movmean(data_filt_30t80(i,:).^2,200));
    data_filt_rms(i,:) = sqrt(movmean(data_filt(i,:).^2,200));
end
% %% for visualizing the activity of each channel
% figure()
% title(['file ' num2str(file)])
% hold on
% plot(data_filt_rms(8,:),Color= 'r')
% plot(data_filt_rms(47,:), Color= 'k')
% legend({'V1' 'S1'});
% hold off
%
% figure()
% title(['file ' num2str(file)])
% hold on
% plot(data_filt(8,:),Color= 'r')
% plot(data_filt(47,:), Color= 'k')
% legend({'V1' 'S1'});
% hold off
%figure()
%hold on
%plot(data_rms_high(2,:),Color= 'r')
%plot(data_rms_high(5,:), Color= 'k')
%legend({'V1 L4' 'S1 L4'});
%hold off


%% fitdist to get mean and SD

clear MAD med
for ch = 1:length(chanals)
data_MAD(ch) = median(abs(data_filt_rms(ch,5000:end-5000)-median(data_filt_rms(ch,5000:end-5000),2))/0.6745,2);

data_med(ch) = median(data_filt_rms(ch,5000:end-5000),2);

end
%% In case youd like to visualize the thresholds
%chcheck = 4;
%figure()
%title('4-100Hz RMS and MADs')
%hold on
%plot(smooth(data_filt_rms(chcheck,:),10),Color= 'k')
%yline(med_low(chcheck)+MAD_low(chcheck), Color='r');
%yline(med_low(chcheck)+(MAD_low(chcheck)*3), Color='b');
%yline(med_low(chcheck)+(MAD_low(chcheck)*5), Color='g');
%legend({'filtered_rms' 'mode+MAD' 'mode+3MAD' 'mode+5MAD'});
%hold off

%% get points where rms >= mu + SD (borrow code from previous attempts)for ch = 1:length(chanals)
for ch = 1:length(chanals)

    input_data = data_filt_rms(ch,:);
    thresholded = input_data > (data_med(ch)+(data_MAD(ch)*2)); %this is effectively the lower threshold, later is a step that requires a peak above 3 MAD

%%

    start = find(diff(thresholded) > 0);
    stop = find(diff(thresholded) < 0);

% Set a stop to last oscillation if it is incomplete (end with end of
% recording)
    if length(stop) == length(start) - 1
            stop(end + 1) = length(thresholded);
    end

% Set a start to first oscillation if it is incomplete (starts with
% beginning of recording)
    if length(stop) - 1 == length(start)
            start(2 : end + 1) = start;
            start(1) = 1;
    end

% Correct special case when both first and last oscillations are incomplete
    if start(1) > stop(1)
            start(2 : end + 1) = start;
            start(1) = 1;
            stop(end + 1) = length(thresholded);
    end

    firstPass = [start',stop'];
% Merge oscillations if inter-oscillation period is too short
    minInterOscSamples = 500 / 1000 * fs; %the first number is your chosen minimum number of ms between events
            secondPass = [];
            oscillation = firstPass(1,:);
            for i = 2 : size(firstPass, 1)
                if firstPass(i,1) - oscillation(2) < minInterOscSamples
                % Merge
                    oscillation = [oscillation(1) firstPass(i,2)];
                else
                    secondPass = [secondPass ; oscillation];
                    oscillation = firstPass(i,:);
                end
            end
            secondPass = [secondPass ; oscillation];


    % Discard oscillations with a relative peak power < rel_thresholds(2)
            % and absolute peak power < abs_thresholds(2)
            thirdPass = [];
            peakNormalizedPower = [];
            %peakAbsPower = [];
            for i = 1 : size(secondPass, 1)
                maxValue_rel = max(input_data([secondPass(i, 1) : secondPass(i, 2)]));
                %maxValue_abs = max(convolvedSignal([secondPass(i, 1) : secondPass(i, 2)]));
                if maxValue_rel > data_med+(data_MAD*3) %rel_thresholds(2) || maxValue_abs > abs_thresholds(2)
                    thirdPass = [thirdPass ; secondPass(i, :)];
                    peakNormalizedPower = [peakNormalizedPower ; maxValue_rel];
                    %peakAbsPower = [peakAbsPower ; maxValue_abs];

                end
            end


            % Detect negative peak position for each oscillation
            peakPosition = zeros(size(thirdPass,1),1);
            for i = 1 : size(thirdPass, 1)
                [~, minIndex] = min(data_filt(ch,thirdPass(i, 1) : thirdPass(i, 2)));
                peakPosition(i) = minIndex + thirdPass(i,1) - 1;
            end

            if ~ isempty(thirdPass)
                % Discard oscillations that are too short
                oscillations_ch = [(thirdPass(:,1)) (peakPosition) ...
                    (thirdPass(:,2)) peakNormalizedPower];
                duration = oscillations_ch(:,3) - oscillations_ch(:,1);
                oscillations_ch(duration < 300) = NaN;%durations(2), :) = NaN;
                oscillations_ch = oscillations_ch((all((~ isnan(oscillations_ch)), 2)), :);

                %remove bouts too close to the edge
                oscillations_ch((oscillations_ch(:,1)<(10*fs)),:) = [];
                oscillations_ch((oscillations_ch(:,3)>length(input_data)-(10*fs)),:) = [];

                %extend bouts for future needs

                ext_osc(:,1) = oscillations_ch(:,1)-(1*fs);
                ext_osc(:,2) = oscillations_ch(:,3)+(1*fs);

                %create timestamps for oscillation negative periods

                neg_osc = [ext_osc(1:end-1,2)' ; ext_osc(2:end,1)']';

                neg_osc((neg_osc(:,1)<(10*fs)),:) = []; %remove ones at start
                neg_osc((neg_osc(:,2)>length(input_data)-(10*fs)),:) = []; %remove ones at end

                %all-channel filtered-signal snippet for each event, over its extended window
                for snip_i = 1:size(ext_osc,1)
                    snippets.event_data{ch}{end+1,1} = data_filt(:, ext_osc(snip_i,1):ext_osc(snip_i,2));
                end

                %% put into a structure
                osc = oscillations_ch; clear oscillations_ch

                oscillations.chEvent_timestamps{ch} = [oscillations.chEvent_timestamps{ch}; osc(:, [1 3])];
                oscillations.chEvent_ext_timestamps{ch} = [oscillations.chEvent_ext_timestamps{ch}; ext_osc];
                oscillations.chEvent_neg_timestamps{ch} = [oscillations.chEvent_neg_timestamps{ch};  neg_osc]; %clear neg_osc
                oscillations.chEvent_peaks{ch} = [oscillations.chEvent_peaks{ch}; osc(:, 2)]; % peaktimes
                oscillations.chEvent_peakNormedPower{ch} = [oscillations.chEvent_peakNormedPower{ch}; osc(:, 4)]; % amplitudes
                %oscillations.mu = mean_low;
                %oscillations.stdev = SD;
                oscillations.chEvent_durations{ch} = [oscillations.chEvent_durations{ch}; osc(:, 3)-osc(:,1)];
                oscillations.chEvent_ext_durations{ch} = [oscillations.chEvent_ext_durations{ch}; ext_osc(:,2)-ext_osc(:,1)];
                oscillations.chEvent_file_index{ch} = [oscillations.chEvent_file_index{ch}; repmat(file,size(osc,1),1)];
                %oscillations.negevent_file_index{ch} = [oscillations.negevent_file_index{ch}; repmat(file,size(neg_osc,1),1)];
                oscillations.chEvent_index{ch} = [oscillations.chEvent_index{ch}; [1:size(osc,1)]'];

                for evnt = 1:size(osc,1)
                oscillations.chEvent_during_stim{ch} = [oscillations.chEvent_during_stim{ch}; max(ismember(stim_times, osc(evnt,1):osc(evnt,3)))];
                end

                clear ext_osc neg_osc
                %oscillations.durations = duration(duration >= durations(2)); % take only those that are long enough
                %oscillations.peakAbsPower = peakAbsPower(duration >= durations(2));
                %oscillations.nnz_norm = nnz(normSignal > rel_thresholds(1)) / length(signal);
                %oscillations.nnz_abs = nnz(convolvedSignal > abs_thresholds(1)) / length(signal);
                %oscillations.fs = fs;
                %oscillations.len_rec = length(signal);
            %else
            %    oscillations.timestamps{ch} = oscillations.timestamps{ch};
            %    oscillations.peaks{ch} = [oscillations.peaks{ch}]; % peaktimes
            %    oscillations.peakNormedPower{ch} = oscillations.peakNormedPower{ch}; % amplitudes
                %oscillations.stdev{ch} = NaN;
                %oscillations.durations = NaN;
                %oscillations.peakAbsPower = NaN;
                %oscillations.nnz_norm = NaN;
                %oscillations.nnz_abs = NaN;
            end

%% here ends the channel loop
end

%% find probe events (prbEvent): union of this probe's channels' chEvent windows, merging gaps under 300ms
prb_gap_samples = 300/1000*fs; %300ms silence gap between probe events (see header NOTE on this constant)

for p = 1:length(probes)
    probe_chans_idx = find(chan_probe == p);

    all_windows = [];
    for pc = 1:length(probe_chans_idx)
        ch2 = probe_chans_idx(pc);
        all_windows = [all_windows; oscillations.chEvent_timestamps{ch2}(oscillations.chEvent_file_index{ch2}==file,:)];
    end

    if ~isempty(all_windows)
        all_windows = sortrows(all_windows,1);

        merged = all_windows(1,:);
        prbFirstPass = [];
        for k = 2:size(all_windows,1)
            if all_windows(k,1) - merged(2) < prb_gap_samples
                merged(2) = max(merged(2), all_windows(k,2)); %union of overlapping/close windows
            else
                prbFirstPass = [prbFirstPass; merged];
                merged = all_windows(k,:);
            end
        end
        prbFirstPass = [prbFirstPass; merged];

        prb_ext = [prbFirstPass(:,1)-(1*fs), prbFirstPass(:,2)+(1*fs)];

        for m = 1:size(prbFirstPass,1)
            window_filt = data_filt(probe_chans_idx, prbFirstPass(m,1):prbFirstPass(m,2));
            [~, lin_idx] = min(window_filt(:));
            [~, col] = ind2sub(size(window_filt), lin_idx);
            prb_peak = prbFirstPass(m,1) + col - 1;
            prb_peak_power = max(max(data_filt_rms(probe_chans_idx, prbFirstPass(m,1):prbFirstPass(m,2))));

            oscillations.prbEvent_timestamps{p} = [oscillations.prbEvent_timestamps{p}; prbFirstPass(m,:)];
            oscillations.prbEvent_ext_timestamps{p} = [oscillations.prbEvent_ext_timestamps{p}; prb_ext(m,:)];
            oscillations.prbEvent_peaks{p} = [oscillations.prbEvent_peaks{p}; prb_peak];
            oscillations.prbEvent_peakNormedPower{p} = [oscillations.prbEvent_peakNormedPower{p}; prb_peak_power];
            oscillations.prbEvent_durations{p} = [oscillations.prbEvent_durations{p}; prbFirstPass(m,2)-prbFirstPass(m,1)];
            oscillations.prbEvent_ext_durations{p} = [oscillations.prbEvent_ext_durations{p}; prb_ext(m,2)-prb_ext(m,1)];
            oscillations.prbEvent_file_index{p} = [oscillations.prbEvent_file_index{p}; file];
            oscillations.prbEvent_index{p} = [oscillations.prbEvent_index{p}; m];
            oscillations.prbEvent_during_stim{p} = [oscillations.prbEvent_during_stim{p}; max(ismember(stim_times, prbFirstPass(m,1):prbFirstPass(m,2)))];

            probeSnippets.event_data{p}{end+1,1} = data_filt(:, prb_ext(m,1):prb_ext(m,2)); %ALL analyzed channels
        end

        %gaps between this file's own consecutive extended prbEvent windows (mirrors chEvent_neg_timestamps)
        prb_neg = [prb_ext(1:end-1,2), prb_ext(2:end,1)];
        oscillations.prbEvent_neg_timestamps{p} = [oscillations.prbEvent_neg_timestamps{p}; prb_neg];
    end
end

%% plot all events across channels
plotdata = NaN(length(chanals),size(data_filt,2));

for ch = 1:length(chanals)

    chdata = data_filt(ch,:);
    if ~isempty(oscillations.chEvent_timestamps{ch}(oscillations.chEvent_file_index{ch}==file,:))

        events = oscillations.chEvent_timestamps{ch}(oscillations.chEvent_file_index{ch}==file,:);

        for i = 1:size(events,1)
        plotdata(ch,events(i,1):events(i,2)) = chdata(events(i,1):events(i,2));
        end
    end

end

figure();
for ch = 1:size(plotdata,1)
hold on
plot(plotdata(ch,:)-(500*(ch-1)))

end
hold off
%% here eneds the folder loop
end

end
