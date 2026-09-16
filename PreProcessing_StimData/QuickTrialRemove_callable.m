% QuickTrialRemove_callable.m
%
% Description: Callable trial-curation utility. Lets the user visually
%   inspect LFP and spike-raster traces for one channel of an
%   already-preprocessed animal/condition (output of
%   OpenEphys_BaseAnalysis.m or a variant), pick trials to exclude, and
%   record that choice as a NAMED, TEMPORARY exclusion set in
%   ProbeInfo.tr_remove_conditional (long-format table: Name/Filename/
%   TrialIdx). Does NOT touch tr_keep or tr_remove, and does NOT modify
%   the LFP/spiking/TF/CSD result files - the recorded trial set is
%   applied later, at analysis time, by whichever script chooses to
%   combine it with tr_remove into a tr_remove_local mask. On entry,
%   presents any exclusion sets already named for this filename so the
%   user can redo (overwrite) one instead of re-picking trials from
%   scratch under a new name.
%
% Inputs:
%   animal (string) - animal name matching the folder under
%     save_directory that holds this recording's result files.
%   filename (string) - animal name + condition name, matching the
%     prefix of the LFP/spiking result files to read (e.g.
%     '20260423-p12-whisker').
%   chan (int) - channel to inspect/plot.
%   fs (int) - sampling rate of the original recording, in Hz (used to
%     bin spikes into 1 ms bins).
%   ProbeInfo (struct) - PROVIDED BY CALLER. This animal's ProbeInfo
%     struct (as loaded from <animal>-ProbeInfo.mat). Must eventually
%     contain field .tr_remove_conditional (table, columns Name/
%     Filename/TrialIdx); added automatically with the correct schema
%     if missing (backward compatibility with ProbeInfo.mat files saved
%     before this field existed).
%   save_directory (string, optional) - base path containing the LFP/
%     Spiking/ProbeInfo subfolders. Defaults to
%     'E:\Roy\Processed Silicon Probe Data'.
%   stim_spike_stimchunks (double, trials x 30kHz trial length x
%     channels), tr_keep (int vector, trials), tr_remove (0/1 vector,
%     same length as tr_keep) - loaded automatically from
%     <save_directory>\Spiking\<animal>\<filename>-spiking_results.mat.
%   stim_lfp_stimchunks - loaded automatically from
%     <save_directory>\LFP\<animal>\<filename>-LFP.mat. 3D array of lfp
%     data split into trials, downsampled to 1kHz. size[trials x trial
%     length x channels].
%
% Outputs:
%   ProbeInfo (struct) - the input struct with .tr_remove_conditional
%     updated: any existing rows for this (Name, Filename) pair are
%     replaced with the newly chosen trial set. Also saved (full
%     overwrite) to <save_directory>\ProbeInfo\<animal>-ProbeInfo.mat.
%
% Dependencies: Expects OpenEphys_BaseAnalysis.m (or the _Bundled /
%   _MixedTrials variant) to have already been run for this
%   animal/condition, producing the LFP.mat and spiking_results.mat
%   files this function reads, and <animal>-ProbeInfo.mat to already
%   exist.

function ProbeInfo = QuickTrialRemove_callable(animal, filename, chan, fs, ProbeInfo, save_directory)

if nargin < 6 || isempty(save_directory)
    save_directory = 'E:\Roy\Processed Silicon Probe Data';
end

if ~isfield(ProbeInfo,'tr_remove_conditional')
    ProbeInfo.tr_remove_conditional = table('Size',[0,3], 'VariableTypes',{'string','string','double'}, 'VariableNames',{'Name','Filename','TrialIdx'});
end

%% load spiking + LFP data automatically
spiking_dat = load(fullfile([save_directory '\Spiking\' animal '\' filename '-spiking_results.mat']), 'stim_spike_stimchunks','tr_keep','tr_remove');
stim_spike_stimchunks = spiking_dat.stim_spike_stimchunks;
tr_keep = spiking_dat.tr_keep;
tr_remove = spiking_dat.tr_remove;

lfp_dat = load(fullfile([save_directory '\LFP\' animal '\' filename '-LFP.mat']), 'stim_lfp_stimchunks');
stim_lfp_stimchunks = lfp_dat.stim_lfp_stimchunks;

tr_keep_local = tr_keep(~logical(tr_remove)); %only offer trials not already permanently excluded

%% bin spikes into 1ms bins for plotting
spikes_ms = [];
for ch = 1:size(stim_spike_stimchunks,3)
    chspikes = squeeze(stim_spike_stimchunks(:,:,ch));
    spikeper = [];
    for batch = 1:size(stim_spike_stimchunks,2)/(fs/1000)
        spikebatchi = sum(chspikes(:,1+((fs/1000)*(batch-1)):(fs/1000)+((fs/1000)*(batch-1))),2)>0;
        spikeper(:,batch) = spikebatchi;
    end
    spikes_ms(:,:,ch) = spikeper;
end

%% pick an existing named condition to redo, or create a new one
existingNames = unique(ProbeInfo.tr_remove_conditional.Name(strcmp(ProbeInfo.tr_remove_conditional.Filename, filename)));

if isempty(existingNames)
    condName = inputdlg({'Name this trial-exclusion condition:'}, 'New condition', [1 50]);
    condName = condName{1};
else
    listOptions = [cellstr(existingNames); {'+ New condition'}];
    [selIdx, ok] = listdlg('ListString', listOptions, 'SelectionMode','single', ...
        'PromptString', {'Select a condition to redo (overwrite),' 'or create a new one:'}, ...
        'Name', filename);
    if ok ~= 1
        error('QuickTrialRemove_callable:Cancelled', 'No condition selected - aborting.');
    end
    if selIdx <= numel(existingNames)
        condName = char(existingNames(selIdx));
    else
        condName = inputdlg({'Name this new trial-exclusion condition:'}, 'New condition', [1 50]);
        condName = condName{1};
    end
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
        plot(time,squeeze(spikes_ms(tr_keep_local(tr),4500:8000,chan))'-(2*(tr-1)))
        xline(0)
    end

    promt = {sprintf('Which of these trials (plot row number, 1-%d) would you like to remove for condition ''%s''?', length(tr_keep_local), condName)};
    removed_rownum = inputdlg(promt, 'Conditional trial removal', [1 50]);
    removed_rownum = str2num(removed_rownum{1});
    chosen_trialidx = tr_keep_local(removed_rownum);

    quest = sprintf('Remove %d trial(s) under condition ''%s''?', numel(chosen_trialidx), condName);
    answer = questdlg(quest, 'Confirm conditional removal', 'Yes', 'No, redo selection', 'No, redo selection');
    close all
    if strcmp(answer, 'Yes')
        confirmed = true;
    end
end

%% overwrite any existing rows for this (Name, Filename) pair, then append the new selection
existingRows = strcmp(ProbeInfo.tr_remove_conditional.Name, condName) & strcmp(ProbeInfo.tr_remove_conditional.Filename, filename);
ProbeInfo.tr_remove_conditional(existingRows,:) = [];
newRows = table(repmat(string(condName),numel(chosen_trialidx),1), repmat(string(filename),numel(chosen_trialidx),1), chosen_trialidx(:), ...
    'VariableNames', {'Name','Filename','TrialIdx'});
ProbeInfo.tr_remove_conditional = [ProbeInfo.tr_remove_conditional; newRows];

save(fullfile([save_directory '\ProbeInfo\' animal '-ProbeInfo.mat']), 'ProbeInfo');

end
