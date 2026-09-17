% ConditionalTrialRemove_callable.m
%
% Description: Callable, loop-friendly trial-curation utility. For one
%   animal/condition/channel, checks whether a NAMED trial-exclusion
%   condition (tr_conditional) already has a recorded trial set in this
%   file's own tr_remove_conditional table (table, columns Name/
%   TrialIdx, stored inside the LFP/spiking/TF/CSD result files
%   themselves - see Outputs). If a match already exists, does nothing
%   and returns immediately (no figures, no dialogs, no save) - this is
%   what makes it fast to call in a loop over many files, most of which
%   already have this condition recorded. If no match exists, tells the
%   user via a dialog, then lets them visually inspect LFP and
%   spike-raster traces for one channel of an already-preprocessed
%   animal/condition (output of OpenEphys_BaseAnalysis.m or a variant),
%   pick trials to exclude, and records that choice under the
%   tr_conditional name. Does NOT touch tr_keep, and does NOT modify the
%   LFP/spiking/TF/CSD result files beyond the tr_remove_conditional
%   write described in Outputs and one other exception: if tr_remove is
%   missing or empty in spiking_results.mat (old-format data that
%   predates trial-removal tracking), it is migrated in place -
%   defaulted to all-zeros (no trials removed) and written back, via
%   save(...,'-append'), to the LFP/spiking/TF/CSD result files for this
%   animal/condition, so all four stay in sync. The recorded exclusion
%   set itself is applied later, at analysis time, by whichever script
%   chooses to combine it with tr_remove into a tr_remove_local mask.
%   Separately, every call (whether it skips via a match or records a
%   new selection) also checks/updates a single master tracking file,
%   <save_directory>\tr_conditional_master.mat (table, columns Name/
%   Animal/Filename - no TrialIdx, that stays in the per-file tables) -
%   an index of which tr_conditional names exist and which animal/
%   filename pairs each has been applied to, so that can be scanned on
%   its own without opening every per-recording result file. See
%   updateTrConditionalMaster (local function, bottom of this file).
%
% Inputs:
%   animal (string) - animal name matching the folder under
%     save_directory that holds this recording's result files.
%   condition (string) - condition name only (e.g. 'whisker'), NOT
%     animal-prefixed. Internally combined with animal as
%     full_filename = [animal '-' condition] (e.g.
%     '20260423-p12-whisker'), matching the <animal>-<condition> prefix
%     that OpenEphys_BaseAnalysis.m (or a variant) actually names the
%     LFP/spiking/figure files with.
%   tr_conditional (string) - REQUIRED. Name of the trial-exclusion
%     condition to look up/record for this file (e.g.
%     'whisker_twitch_artifact'). Not optional and has no global/
%     persistent fallback - the caller is expected to already know what
%     they want to select before starting (e.g. set once in their own
%     loop/script and pass it in on every call), since MATLAB doesn't
%     allow a global or persistent variable to share a name with a
%     function parameter, which ruled out a "prompt once, remember
%     across calls" fallback under this same name.
%   fs (int, optional) - sampling rate of the original recording, in Hz
%     (used to bin spikes into 1 ms bins). Defaults to 30000.
%   ProbeInfo (struct, optional) - This animal's ProbeInfo struct (as
%     loaded from <animal>-ProbeInfo.mat). If omitted or left empty
%     ([]), loaded automatically from
%     <save_directory>\ProbeInfo\<animal>-ProbeInfo.mat. Used only for
%     .Areas (to locate the LFP overview figures below) - no longer
%     carries tr_remove_conditional (moved to per-file storage; see
%     Outputs).
%   save_directory (string, optional) - base path containing the LFP/
%     Spiking/ProbeInfo subfolders. Defaults to
%     'E:\Roy\Processed Silicon Probe Data'.
%   tr_keep (int vector, trials), tr_remove (0/1 vector, same length as
%     tr_keep), tr_remove_conditional (table, columns Name/TrialIdx;
%     defaults to an empty table if missing from the file. TrialIdx may
%     be NaN for a given Name - that means this condition was evaluated
%     for this file and the user confirmed zero trials removed (keep
%     all), not a literal trial index; any code that later reads this
%     table to build a tr_remove_local mask must treat a NaN TrialIdx as
%     a no-op rather than an index) - loaded automatically, up front, in
%     one lightweight call, from
%     <save_directory>\Spiking\<animal>\<full_filename>-spiking_results.mat.
%     If tr_remove is missing or empty (old-format file), it is
%     migrated to zeros(1,length(tr_keep)) and written back - see
%     Description. This initial load deliberately excludes
%     stim_spike_stimchunks (the heavy 30kHz array) - see chan below for
%     when/how that's read.
%   stim_lfp_stimchunks - loaded automatically from
%     <save_directory>\LFP\<animal>\<full_filename>-LFP.mat. 3D array of
%     lfp data split into trials, downsampled to 1kHz. size[trials x
%     trial length x channels]. Loaded in full (unlike
%     stim_spike_stimchunks below) because this file is not saved with
%     '-v7.3', so a matfile partial-channel read would still have to
%     pull the whole array off disk anyway - no benefit to deferring it.
%   chan (int) - NOT a function input. Prompted for via dialog box after
%     the LFP data is loaded (bounded 1 to the number of channels in
%     stim_lfp_stimchunks), so the caller doesn't need to know the
%     channel count in advance. Before prompting, opens the saved LFP
%     overview figure(s) for this animal/condition (one per probe/
%     region, from <save_directory>\AnimalFigures\<animal>\
%     <full_filename>-<ProbeInfo.Areas{prb}>_Probe-lfp_results.fig) so
%     the user can reference channel layout when choosing. Once chosen,
%     stim_spike_stimchunks is read for this channel only, via a
%     matfile object indexed as stim_spike_stimchunks(:,:,chan), so only
%     one channel's data is ever held in memory afterward (checked via
%     checkcode + a synthetic-data correctness test - bit-identical to
%     the old all-channel-then-slice approach). NOT a disk I/O win,
%     though: checked with h5info against a real saved file, MATLAB's
%     default '-v7.3' auto-chunking on this array picks a chunk shape
%     spanning ALL trials and ALL channels (e.g. [40 6 32] for a
%     [40 trials x 105000 samples x 32 channels] array - just a thin
%     slice of samples per chunk), so every chunk has to be touched
%     regardless of which channel is requested; a timed comparison
%     against a full load showed no measurable speed difference. Fixing
%     that would require re-chunking OpenEphys_BaseAnalysis.m's save
%     calls and re-saving existing files - out of scope here. Binning
%     below operates on just this channel's slice. INFERRED DATA
%     CONTRACT (pre-existing, not introduced by this change - the old
%     all-channel binning loop and this per-channel read both assume
%     it): chan
%     indexes the same channel in both stim_lfp_stimchunks and
%     stim_spike_stimchunks, i.e. the two arrays have matching channel
%     count and ordering.
%
% Outputs:
%   ProbeInfo (struct) - the input struct, returned largely unchanged
%     (still useful to the caller for .Areas etc. without reloading it
%     themselves in a loop). No longer carries tr_remove_conditional.
%   (side effect, not a return value) tr_remove_conditional - if no
%     match was found for tr_conditional, the newly chosen trial set is
%     appended (as Name/TrialIdx rows) and written, identically, into
%     all four of this animal/condition's result files:
%     <full_filename>-spiking_results.mat, <full_filename>-LFP.mat,
%     <full_filename>-CSD_results.mat, and
%     <full_filename>_TF_results.mat (note: no dash before
%     "_TF_results", unlike the other three - matches how
%     OpenEphys_BaseAnalysis.m actually names that file).
%     TESTING: these four save() calls are currently commented out, and
%     the rows that would have been written are printed to the console
%     instead - re-enable them once verified.
%   (side effect, not a return value) tr_conditional_master - the
%     Name/Animal/Filename row for this call is added if not already
%     present, in <save_directory>\tr_conditional_master.mat. Updated on
%     every call, including the "match found, skip" path. TESTING: this
%     save() call is also currently commented out, printed instead -
%     re-enable together with the four above once verified.
%
% Dependencies: Expects OpenEphys_BaseAnalysis.m (or the _Bundled /
%   _MixedTrials variant) to have already been run for this
%   animal/condition, producing the LFP.mat and spiking_results.mat
%   files this function reads, and <animal>-ProbeInfo.mat to already
%   exist. Also expects ProbeInfo.Areas (cell array of probe/region
%   names) to be populated, and (non-fatal if missing - a warning is
%   issued instead) the per-probe lfp_results.fig files under
%   AnimalFigures to have been saved for this animal/condition.

function ProbeInfo = ConditionalTrialRemove_callable(animal, condition, tr_conditional, fs, ProbeInfo, save_directory)

if nargin < 6 || isempty(save_directory)
    save_directory = 'E:\Roy\Processed Silicon Probe Data';
end

if nargin < 5 || isempty(ProbeInfo)
    probeinfo_dat = load(fullfile([save_directory '\ProbeInfo\' animal '-ProbeInfo.mat']), 'ProbeInfo');
    ProbeInfo = probeinfo_dat.ProbeInfo;
end

if ~isstruct(ProbeInfo) || ~isscalar(ProbeInfo)
    error('ConditionalTrialRemove_callable:InvalidProbeInfo', ...
        ['ProbeInfo must be a scalar struct, got a %s. This usually means a stale call ' ...
         'still passes chan as the 3rd argument - the current signature is ' ...
         'ConditionalTrialRemove_callable(animal, condition, tr_conditional, fs, ProbeInfo, save_directory), ' ...
         'with chan now chosen via dialog instead of as an input.'], class(ProbeInfo));
end

if nargin < 4 || isempty(fs)
    fs = 30000;
end

if ~isnumeric(fs) || ~isscalar(fs)
    error('ConditionalTrialRemove_callable:InvalidFs', ...
        ['fs must be a numeric scalar, got a %s. This usually means fs was omitted from a call ' ...
         'that still supplied later arguments, letting one of them (often ProbeInfo) land in the ' ...
         'fs slot instead - the current signature is ConditionalTrialRemove_callable(animal, ' ...
         'condition, tr_conditional, fs, ProbeInfo, save_directory). Pass [] for fs to use its ' ...
         'default (30000) while still supplying later arguments.'], class(fs));
end

full_filename = [animal '-' condition]; % matches the <animal>-<condition> prefix OpenEphys_BaseAnalysis.m actually writes files under
spikingFile = fullfile([save_directory '\Spiking\' animal '\' full_filename '-spiking_results.mat']);

%% light load: just enough to decide whether to skip and to migrate tr_remove if needed - deliberately excludes
%% stim_spike_stimchunks (the heavy 30kHz array), so a "match found, skip" call never touches it. See below for
%% where the chosen channel's spike data is read on its own, only once we know we're not skipping.
spiking_dat = load(spikingFile, 'tr_keep','tr_remove','tr_remove_conditional');
tr_keep = spiking_dat.tr_keep;

if ~isfield(spiking_dat,'tr_remove_conditional')
    tr_remove_conditional = table('Size',[0,2], 'VariableTypes',{'string','double'}, 'VariableNames',{'Name','TrialIdx'});
else
    tr_remove_conditional = spiking_dat.tr_remove_conditional;
end

%% if this named condition already has a recorded trial set for this file, skip everything below (loop-friendly: no dialogs, no figures, no save)
if any(strcmp(tr_remove_conditional.Name, tr_conditional))
    fprintf('Match found for ''%s'' in %s - skipping.\n', tr_conditional, full_filename);
    updateTrConditionalMaster(save_directory, tr_conditional, animal, full_filename);
    return
end
% uiwait blocks script execution until the user clicks OK (or closes the box);
% msgbox itself stays non-modal (WindowStyle 'normal'), so other figures are still clickable while it waits.
uiwait(msgbox(sprintf('No match found for ''%s'' in %s. Please select trials to remove.', tr_conditional, full_filename)));

%% migrate tr_remove from the old format (missing/empty field) to the new format (0/1 vector matching tr_keep), writing the fix back to every related result file
if ~isfield(spiking_dat,'tr_remove') || isempty(spiking_dat.tr_remove)
    tr_remove = zeros(1,length(tr_keep));
    save(spikingFile,'tr_remove','-append');
    save(fullfile([save_directory '\CSD\' animal '\' full_filename '-CSD_results.mat']),'tr_remove','-append');
    save(fullfile([save_directory '\TF\' animal '\' full_filename '_TF_results.mat']),'tr_remove','-append');
    save(fullfile([save_directory '\LFP\' animal '\' full_filename '-LFP.mat']),'tr_remove','-append');
else
    tr_remove = spiking_dat.tr_remove;
end
lfp_dat = load(fullfile([save_directory '\LFP\' animal '\' full_filename '-LFP.mat']), 'stim_lfp_stimchunks');
stim_lfp_stimchunks = lfp_dat.stim_lfp_stimchunks;

tr_keep_local = tr_keep(~logical(tr_remove)); %only offer trials not already permanently excluded

%% open the saved LFP overview figure(s) for this animal/condition (one per probe/region), for channel reference
for prb = 1:numel(ProbeInfo.Areas)
    lfpFigFile = fullfile([save_directory '\AnimalFigures\' animal '\' full_filename '-' ProbeInfo.Areas{prb} '_Probe-lfp_results.fig']);
    if isfile(lfpFigFile)
        openfig(lfpFigFile);
    else
        warning('ConditionalTrialRemove_callable:LFPFigNotFound', 'LFP overview figure not found for %s, probe ''%s'': %s', full_filename, ProbeInfo.Areas{prb}, lfpFigFile);
    end
end

%% ask which channel to inspect
% Every dialog in this function is kept non-modal (WindowStyle 'normal'
% on inputdlg calls; listdlg is non-modal by default; questdlg has no
% non-modal option at all, so its Yes/No use below is replaced with the
% nonModalConfirm() local function), so the user can always click into
% the LFP overview figure(s) or trial plots while deciding.
nchan = size(stim_lfp_stimchunks,3);
chanAnswer = inputdlg({sprintf('Which channel would you like to inspect (1-%d)?', nchan)}, 'Select channel', [1 50], {''}, struct('WindowStyle','normal'));
chan = str2double(chanAnswer{1});

%% read only the chosen channel's spike data (keeps just 1 channel in memory afterward instead of all
%% of them; NOT a disk I/O win here - see header doc re: this array's HDF5 chunk layout)
spikingMat = matfile(spikingFile);
chspikes = squeeze(spikingMat.stim_spike_stimchunks(:,:,chan));

%% bin that channel's spikes into 1ms bins for plotting
spikes_ms = [];
for batch = 1:size(chspikes,2)/(fs/1000)
    spikebatchi = sum(chspikes(:,1+((fs/1000)*(batch-1)):(fs/1000)+((fs/1000)*(batch-1))),2)>0;
    spikes_ms(:,batch) = spikebatchi;
end

%% pick trials, with a confirm/redo loop
time = [-500:3000];
confirmed = false;
while ~confirmed
    figure();
    for tr = 1:length(tr_keep_local)
        subplot(1,2,1)
        hold on
        colormap(sky(8))
        plot(time,squeeze(stim_lfp_stimchunks(tr_keep_local(tr),4500:8000,chan))'-(400*(tr-1)))
        xline(0)
        subplot(1,2,2)
        hold on
        plot(time,spikes_ms(tr_keep_local(tr),4500:8000)'-(2*(tr-1)))
        xline(0)
    end

    promt = {sprintf('Which of these trials (plot row number, 1-%d) would you like to remove for condition ''%s''?', length(tr_keep_local), tr_conditional)};
    removed_rownum = inputdlg(promt, 'Conditional trial removal', [1 50], {''}, struct('WindowStyle','normal'));
    removed_rownum = str2num(removed_rownum{1});
    chosen_trialidx = tr_keep_local(removed_rownum);

    quest = sprintf('Remove %d trial(s) under condition ''%s''?', numel(chosen_trialidx), tr_conditional);
    answer = nonModalConfirm(quest, 'Confirm conditional removal');
    close all
    if strcmp(answer, 'Yes')
        confirmed = true;
    end
end

%% record the new selection (existingRows was already confirmed empty by the match check above)
% zero trials chosen (keep all) still gets one row, TrialIdx = NaN, so this
% condition/file pair is recorded as evaluated and isn't re-prompted next time
if isempty(chosen_trialidx)
    newRows = table(string(tr_conditional), NaN, 'VariableNames', {'Name','TrialIdx'});
else
    newRows = table(repmat(string(tr_conditional),numel(chosen_trialidx),1), chosen_trialidx(:), ...
        'VariableNames', {'Name','TrialIdx'});
end
tr_remove_conditional = [tr_remove_conditional; newRows];

%% TESTING: print what would be saved instead of saving, so nothing on disk is overwritten
fprintf('\n--- ConditionalTrialRemove_callable: rows that would be added for %s, condition ''%s'' ---\n', full_filename, tr_conditional);
disp(newRows);

% save(fullfile([save_directory '\Spiking\' animal '\' full_filename '-spiking_results.mat']), 'tr_remove_conditional', '-append'); % TESTING: disabled, see note above
% save(fullfile([save_directory '\LFP\' animal '\' full_filename '-LFP.mat']), 'tr_remove_conditional', '-append'); % TESTING: disabled, see note above
% save(fullfile([save_directory '\CSD\' animal '\' full_filename '-CSD_results.mat']), 'tr_remove_conditional', '-append'); % TESTING: disabled, see note above
% save(fullfile([save_directory '\TF\' animal '\' full_filename '_TF_results.mat']), 'tr_remove_conditional', '-append'); % TESTING: disabled, see note above

updateTrConditionalMaster(save_directory, tr_conditional, animal, full_filename);

end

function answer = nonModalConfirm(question, dlgtitle)
% Minimal Yes / No, redo selection confirm dialog. questdlg hard-codes
% WindowStyle to 'modal' with no public override, so it can't be made
% non-modal directly; this reimplements just enough of it (same two
% button labels/behavior) as an ordinary, non-modal figure so the user
% can still click into the trial-plot figure while deciding.
answer = 'No, redo selection';
fig = figure('Name', dlgtitle, 'NumberTitle', 'off', 'MenuBar', 'none', ...
    'ToolBar', 'none', 'WindowStyle', 'normal', 'Resize', 'off', ...
    'Position', [400 400 340 110]);
uicontrol(fig, 'Style', 'text', 'String', question, 'Units', 'normalized', ...
    'Position', [0.05 0.5 0.9 0.4], 'HorizontalAlignment', 'center');
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Yes', 'Units', 'normalized', ...
    'Position', [0.1 0.15 0.3 0.3], 'Callback', @(~,~) respond('Yes'));
uicontrol(fig, 'Style', 'pushbutton', 'String', 'No, redo selection', 'Units', 'normalized', ...
    'Position', [0.55 0.15 0.35 0.3], 'Callback', @(~,~) respond('No, redo selection'));
uiwait(fig);
    function respond(val)
        answer = val;
        uiresume(fig);
        close(fig);
    end
end

function updateTrConditionalMaster(save_directory, tr_conditional, animal, full_filename)
% Keeps <save_directory>\tr_conditional_master.mat in sync: a single,
% small index (table, columns Name/Animal/Filename - no TrialIdx, that
% stays in the per-file tables) of which tr_conditional names exist and
% which animal/filename pairs each has been applied to, so it can be
% scanned on its own without opening every per-recording result file.
masterFile = fullfile([save_directory '\tr_conditional_master.mat']);
if isfile(masterFile)
    master_dat = load(masterFile, 'tr_conditional_master');
    tr_conditional_master = master_dat.tr_conditional_master;
else
    tr_conditional_master = table('Size',[0,3], 'VariableTypes',{'string','string','string'}, 'VariableNames',{'Name','Animal','Filename'});
end

alreadyTracked = strcmp(tr_conditional_master.Name, tr_conditional) & ...
    strcmp(tr_conditional_master.Animal, animal) & ...
    strcmp(tr_conditional_master.Filename, full_filename);
if any(alreadyTracked)
    return
end

newRow = table(string(tr_conditional), string(animal), string(full_filename), ...
    'VariableNames', {'Name','Animal','Filename'});
tr_conditional_master = [tr_conditional_master; newRow];

%% TESTING: print what would be saved instead of saving, so nothing on disk is overwritten
fprintf('--- tr_conditional_master: would add %s / %s / %s ---\n', tr_conditional, animal, full_filename);
% save(masterFile, 'tr_conditional_master'); % TESTING: disabled, see note above

end
