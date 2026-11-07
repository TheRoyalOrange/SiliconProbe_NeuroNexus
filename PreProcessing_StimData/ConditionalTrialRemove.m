% ConditionalTrialRemove.m
%
% Description: Standalone driver script that batch-applies
%   ConditionalTrialRemove_callable.m across a manually chosen set of
%   animals and conditions, under a single tr_conditional name. Meant to
%   be run start-to-finish (F5) rather than section-by-section: it
%   blocks on its own dialogs (confirm tr_conditional, confirm the
%   animal/condition selections) before looping. As close to the very
%   start of the script as possible (right after the user-set variables,
%   before anything else runs), opens a non-modal uitable of the
%   existing master index (master_filename, normally
%   tr_conditional_master.mat) so the user can scan what conditional
%   names already exist before committing to tr_conditional.
%   After the user confirms tr_conditional, checks whether that name
%   already exists in the master index and, if so, warns and asks
%   whether to continue (since re-running the same name on files it has
%   already been applied to is a no-op - ConditionalTrialRemove_callable.m
%   skips any animal/filename pair that already has a recorded trial set
%   for that name, and there is currently no way to edit an
%   already-recorded assignment). Then, before any files are picked,
%   loads every animal's ProbeInfo and checks it has at least one probe
%   activity label (ProbeInfo.ChanLabels, made by
%   NNx_editProbeInfo_ChLabel_StimData.m / _SponData.m). If any animal has no
%   labels (or no ProbeInfo file), a questdlg lists them and offers
%   "Fuck it, we'll do it live" (carry on - those animals fall back to the
%   old LFP-figs + typed channel number pick, passed to the callable as
%   do_it_live = 1 so it doesn't re-ask per file) or "Cancel. I'll be
%   back" (stop before anything runs). For each animal in animal_inc, opens a
%   uipickfiles dialog rooted at that animal's LFP results folder so the
%   user can pick which condition files to run tr_conditional against;
%   the chosen files are parsed into condition names (stripping the
%   <animal>- prefix and -LFP suffix) and handed to
%   ConditionalTrialRemove_callable.m one at a time. For each file that
%   still needs trials picked, the callable has the user pick a channel off
%   the animal's probe activity labels (ProbeInfo.ChanLabels/ShankLabels)
%   if it has any - see that file's header. If the user hits "Cancel. I'll
%   be back" at that step, the callable returns runCancelled = true and
%   this script stops the whole run right there (files finished earlier in
%   the run stay saved). Before any of that, per file, the callable checks
%   that all four result files it saves into (spiking/LFP/CSD/TF) exist.
%   Spiking or LFP missing: "Skip this file, keep going" (added to
%   skipped_files) or "Abort. I'll be back" (runCancelled = true, run
%   stops like above). Only CSD/TF missing: "Skip saving to those, keep
%   going" (trial pick still happens, those files just aren't saved into;
%   added to nosave_files) or abort. Both lists print at the end of the run.
%
% Inputs:
%   tr_conditional (string) - PROVIDED BY USER, manually typed rather
%     than selected from a list, to reduce the chance of silently
%     reusing/mistyping an existing name. Passed straight through to
%     ConditionalTrialRemove_callable.m - see its header for the naming
%     convention (e.g. 'whisker_twitch_artifact').
%   animal_inc (cell array of strings) - PROVIDED BY USER. Each entry is
%     an animal name matching a folder under
%     <save_directory>\LFP\<animal>, i.e. the same animal argument
%     ConditionalTrialRemove_callable.m expects.
%   fs (int) - PROVIDED BY USER. Sampling rate of the original
%     recording, in Hz. Forwarded to ConditionalTrialRemove_callable.m.
%     INFERRED DATA CONTRACT: not requested explicitly, but added since
%     ConditionalTrialRemove_callable.m takes it as an input (defaults to
%     30000, matching that function's own default, if left as-is).
%   save_directory (string) - PROVIDED BY USER. Base path containing the
%     LFP/Spiking/ProbeInfo/tr_conditional_master.mat this script and
%     ConditionalTrialRemove_callable.m read from. INFERRED DATA
%     CONTRACT: not requested explicitly, but added since
%     ConditionalTrialRemove_callable.m takes it as an input; defaults to
%     'E:\Roy\Processed Silicon Probe Data', matching that function's own
%     default.
%   master_filename (string) - PROVIDED BY USER, set just below the main
%     USER-SET VARIABLES block (not inside it) right before its first
%     use. Filename (not a full path) of the master tracking file within
%     save_directory, forwarded to ConditionalTrialRemove_callable.m and
%     used by this script's own master-list display/existence-check
%     below. Normally 'tr_conditional_master.mat' (the real master
%     index) - for a test run, temporarily change it to a different
%     filename (e.g. 'themastertester.mat') so a test run cannot alter
%     the real index; see the comment at its assignment below.
%   (side effect) user picks, per animal in animal_inc, which condition
%     files under <save_directory>\LFP\<animal> to include, via
%     uipickfiles. INFERRED DATA CONTRACT: each selected file is named
%     <animal>-<condition>-LFP.mat (matching how
%     NNx_BaseAnalysis.m, or a variant, actually names LFP.mat
%     output) - the condition name fed to
%     ConditionalTrialRemove_callable.m is recovered by stripping the
%     <animal>- prefix and -LFP suffix back off the picked filename.
%
% Outputs:
%   animal_condition_files (cell array, 1 x numel(animal_inc), each
%     entry a cell array of condition-name strings) - the parsed
%     animal/condition selections gathered from uipickfiles, kept around
%     after the run for reference/debugging. Matches
%     {animals}{filenames for animal} shape.
%   skipped_files (cell array of strings, 1 x files skipped) -
%     <animal>-<condition> of every file the user chose to skip because
%     its spiking or LFP file was missing (callable's fileSkipped
%     output). Listed at the end of a completed run; stays in the
%     workspace if the run is aborted partway.
%   nosave_files (cell array of strings, 1 x files) -
%     '<animal>-<condition>: <missing file names>' for every file whose
%     trials were picked but whose CSD and/or TF file was missing, so
%     tr_remove_conditional wasn't saved there (callable's missingsaves
%     output). Listed/kept the same way as skipped_files.
%   (side effect, via ConditionalTrialRemove_callable.m) tr_remove_conditional
%     rows written into each selected animal/condition's LFP/spiking/
%     CSD/TF result files, and <save_directory>\<master_filename> updated
%     - see ConditionalTrialRemove_callable.m's own header for details.
%
% Dependencies: ConditionalTrialRemove_callable.m (same folder);
%   uipickfiles (must be on the MATLAB path - already used elsewhere in
%   this repo, e.g. Summarize_MultiAnimal\SummaryAnalysis_MultiMouse_allChannels.m);
%   expects NNx_BaseAnalysis.m (or a variant) to have already been
%   run for every animal/condition selected. For testing this script's
%   own loop/wiring without touching real data, see
%   ConditionalTrialRemove_callable_test.m (same folder) and the
%   commented-out swap in the loop below.

%% ============ USER-SET VARIABLES - fill these in by hand ============
tr_conditional = 'ActiveV1';   % name of the trial-exclusion condition to apply, e.g. 'whisker_twitch_artifact'
animal_inc = {'20260616-p7'         % animals to include this run
   % '20260423-p12'
    };
fs = 30000;             % sampling rate of the original recording (Hz)
save_directory = 'E:\Roy\Processed Silicon Probe Data';
%% =====================================================================

%% show the existing master index so the user can scan before proceeding
master_filename = 'tr_conditional_master.mat'; % the real master index - for a test run, temporarily change this to a different filename (e.g. 'themastertester.mat') so you don't alter the real one
displayTrConditionalMaster(save_directory, master_filename);

%% now validate that the user actually filled in the required variables above
if isempty(tr_conditional) || ~(ischar(tr_conditional) || isstring(tr_conditional))
    error('ConditionalTrialRemove:EmptyConditional', 'Fill in tr_conditional before running.');
end
if isempty(animal_inc) || ~iscell(animal_inc) || any(cellfun(@isempty, animal_inc))
    error('ConditionalTrialRemove:EmptyAnimalInc', 'Fill in animal_inc (cell array of animal names) before running.');
end

%% confirm tr_conditional with the user
confirmMsg = sprintf('tr_conditional is set to:\n\n''%s''\n\nHappy with this?', tr_conditional);
answer = nonModalConfirm(confirmMsg, 'Confirm tr_conditional');
if ~strcmp(answer, 'Yes')
    disp('tr_conditional not confirmed - edit the variable at the top of the script and re-run.');
    return
end

%% warn if this tr_conditional name already exists in the master index
masterFile = fullfile(save_directory, master_filename);
if isfile(masterFile)
    master_dat = load(masterFile, 'tr_conditional_master');
    existingRows = strcmp(master_dat.tr_conditional_master.Name, tr_conditional);
    if any(existingRows)
        n = sum(existingRows);
        existMsg = sprintf(['tr_conditional ''%s'' already exists in %s ' ...
            '(%d animal/filename pair(s) already recorded).\n\nNote: already-recorded pairs cannot ' ...
            'be re-edited - ConditionalTrialRemove_callable.m will just skip them.\n\nContinue anyway?'], ...
            tr_conditional, master_filename, n);
        answer = nonModalConfirm(existMsg, 'tr_conditional already exists');
        if ~strcmp(answer, 'Yes')
            disp('Aborted by user - tr_conditional already exists in the master index.');
            return
        end
    end
end

%% check every animal's ProbeInfo has activity labels (used to pick channels) before starting the run
nolabel_animals = {}; %animals whose ProbeInfo has no ChanLabels labels
noprobeinfo_animals = {}; %animals with no ProbeInfo file at all
animal_haslabels = false(1,numel(animal_inc)); %[1 x animals] 1 = has labels; 0 = no labels, so its files skip straight to the old channel pick (do_it_live)
for a = 1:numel(animal_inc)
    animal = animal_inc{a};
    probeinfoFile = fullfile([save_directory '\ProbeInfo\' animal '-ProbeInfo.mat']);
    if isfile(probeinfoFile) == false
        noprobeinfo_animals{end+1} = animal; %#ok<SAGROW> only a few animals
        continue
    end
    probeinfo_dat = load(probeinfoFile, 'ProbeInfo');
    haslabels = false;
    if isfield(probeinfo_dat.ProbeInfo, 'ChanLabels')
        for prb = 1:numel(probeinfo_dat.ProbeInfo.ChanLabels)
            if isstruct(probeinfo_dat.ProbeInfo.ChanLabels{prb}) && ~isempty(fieldnames(probeinfo_dat.ProbeInfo.ChanLabels{prb}))
                haslabels = true; %at least one label on at least one probe is enough for the channel map
            end
        end
    end
    animal_haslabels(a) = haslabels;
    if haslabels == false
        nolabel_animals{end+1} = animal; %#ok<SAGROW>
    end
end
clear probeinfo_dat

%anyone missing? say who, then do it live or bail
if ~isempty(nolabel_animals) || ~isempty(noprobeinfo_animals)
    labelMsg = 'Heads up before the run starts:';
    if ~isempty(nolabel_animals)
        labelMsg = sprintf(['%s\n\nNo probe activity labels exist for:\n%s\n\nRecommend running NNx_editProbeInfo_ChLabel_StimData.m ' ...
            'and/or NNx_editProbeInfo_ChLabel_SponData.m on these first. If you go ahead, you''ll pick their channels the old way ' ...
            '(LFP overview figs + typing the channel number).'], labelMsg, strjoin(nolabel_animals, '\n'));
    end
    if ~isempty(noprobeinfo_animals)
        labelMsg = sprintf('%s\n\nNo ProbeInfo file found at all for:\n%s\n\n(These will error out when their turn comes.)', ...
            labelMsg, strjoin(noprobeinfo_animals, '\n'));
    end
    labelAnswer = questdlg(labelMsg, 'Missing activity labels', 'Fuck it, we''ll do it live', 'Cancel. I''ll be back', 'Cancel. I''ll be back');
    if ~strcmp(labelAnswer, 'Fuck it, we''ll do it live') %cancel button, or closed the box
        disp('ConditionalTrialRemove: run cancelled before starting - go make those labels. See ya later.');
        return
    end
end

%% for each animal, pick which condition files (under its LFP folder) to run this tr_conditional against
animal_condition_files = cell(1, numel(animal_inc));
for a = 1:numel(animal_inc)
    animal = animal_inc{a};
    animalLFPdir = fullfile(save_directory, 'LFP', animal);
    if ~isfolder(animalLFPdir)
        warning('ConditionalTrialRemove:AnimalFolderNotFound', ...
            'LFP folder not found for animal ''%s'': %s - skipping.', animal, animalLFPdir);
        animal_condition_files{a} = {};
        continue
    end

    picked = uipickfiles('FilterSpec', animalLFPdir, 'Prompt', ...
        sprintf('Choose ''%s'' condition file(s) for %s', tr_conditional, animal));

    if ~iscell(picked) % uipickfiles returns 0 if the user cancels
        warning('ConditionalTrialRemove:NoFilesPicked', ...
            'No files picked for animal ''%s'' - skipping.', animal);
        animal_condition_files{a} = {};
        continue
    end

    conditions = cell(1, numel(picked));
    prefix = [animal '-'];
    suffix = '-LFP';
    for f = 1:numel(picked)
        [~, base, ~] = fileparts(picked{f});
        if startsWith(base, prefix) && endsWith(base, suffix)
            conditions{f} = base(numel(prefix)+1 : end-numel(suffix));
        else
            warning('ConditionalTrialRemove:UnexpectedFilename', ...
                ['Picked file ''%s'' does not match the expected <animal>-<condition>-LFP.mat ' ...
                 'pattern for animal ''%s'' - skipping this file.'], base, animal);
            conditions{f} = '';
        end
    end
    animal_condition_files{a} = conditions(~cellfun(@isempty, conditions));
end

%% run ConditionalTrialRemove_callable.m for every animal/condition selected above
skipped_files = {}; %<animal>-<condition> of every file skipped for a missing spiking/LFP file, for the post-run report
nosave_files = {}; %<animal>-<condition>: <missing file names> - trials picked, but not saved into those missing CSD/TF files
for a = 1:numel(animal_inc)
    animal = animal_inc{a};
    conditions = animal_condition_files{a};
    ProbeInfo = []; % reloaded once per animal (via ConditionalTrialRemove_callable.m's own optional ProbeInfo input/output), then reused across this animal's conditions instead of re-loading per condition
    for c = 1:numel(conditions)
        % TESTING: comment out the real call below and uncomment the _test call to dry-run this
        % script's loop/wiring without touching any real LFP/Spiking/CSD/TF result files - see
        % ConditionalTrialRemove_callable_test.m's header for what it does instead.
        % do_it_live = 1 for animals with no activity labels: you already said "do it live" in the pre-run check,
        % so the callable skips its own no-labels dialog and goes straight to the LFP figs + typed channel number
        do_it_live = ~animal_haslabels(a);
        [ProbeInfo, runCancelled, fileSkipped, missingsaves] = ConditionalTrialRemove_callable(animal, conditions{c}, tr_conditional, fs, ProbeInfo, save_directory, master_filename, do_it_live);
        % ProbeInfo = ConditionalTrialRemove_callable_test(animal, conditions{c}, tr_conditional, fs, ProbeInfo, save_directory, master_filename);
        % runCancelled = false; fileSkipped = false; missingsaves = {}; %uncomment along with the _test call above (the _test function has none of these outputs)

        %spiking/LFP file missing and user said skip this file - note it for the report, move on
        if fileSkipped == true
            skipped_files{end+1} = [animal '-' conditions{c}]; %#ok<SAGROW> only a few files
        %trials picked, but CSD/TF file missing so it didn't get saved there - note which for the report
        elseif runCancelled == false && ~isempty(missingsaves)
            nosave_files{end+1} = [animal '-' conditions{c} ': ' strjoin(missingsaves, ', ')]; %#ok<SAGROW>
        end

        %user hit "Cancel. I'll be back" at the channel pick, or "Abort" at a missing result file - stop the whole run here
        if runCancelled == true
            disp('ConditionalTrialRemove: run cancelled - see ya later. Files already done this run stay saved.');
            return
        end
    end
end

%post-run report: anything skipped, or not saved somewhere, for missing result files?
if ~isempty(skipped_files)
    disp('ConditionalTrialRemove: skipped these (missing spiking/LFP file) - nothing saved for them:');
    disp(skipped_files');
end
if ~isempty(nosave_files)
    disp('ConditionalTrialRemove: trials picked for these, but NOT saved into these missing files:');
    disp(nosave_files');
end
disp('ConditionalTrialRemove: done.');

%% ------------------------- local functions -------------------------

function displayTrConditionalMaster(save_directory, master_filename)
% Opens a non-modal uitable of <save_directory>\<master_filename>
% (sorted by Name, then Animal, then Filename) so the user can scan
% which tr_conditional names already exist, and on which animal/filename
% pairs, before committing to a tr_conditional above. Left open (not
% uiwait'd) so it stays visible/referenceable through the rest of the
% script's dialogs.
masterFile = fullfile(save_directory, master_filename);
if ~isfile(masterFile)
    fprintf('No %s found yet - this will be the first tr_conditional recorded.\n', master_filename);
    return
end
master_dat = load(masterFile, 'tr_conditional_master');
sortedMaster = sortrows(master_dat.tr_conditional_master, {'Name','Animal','Filename'});

fig = figure('Name', sprintf('%s (existing conditionals)', master_filename), 'NumberTitle', 'off', ...
    'MenuBar', 'none', 'ToolBar', 'none', 'Position', [200 200 640 400]);
% uitable(fig, 'Data', table2cell(sortedMaster), 'ColumnName', sortedMaster.Properties.VariableNames, ...
%     'Units', 'normalized', 'Position', [0 0 1 1]); %errors once the master file exists: its columns are string, and a figure() uitable only takes numeric/logical/char
uitable(fig, 'Data', table2cell(convertvars(sortedMaster, @isstring, 'cellstr')), 'ColumnName', sortedMaster.Properties.VariableNames, ...
    'Units', 'normalized', 'Position', [0 0 1 1]); %string columns -> char so uitable will show them
end

function answer = nonModalConfirm(question, dlgtitle)
% Minimal Yes / No confirm dialog, non-modal so the user can still click
% into the tr_conditional_master table (or anything else) while
% deciding. Duplicated (rather than shared) from the near-identical
% local function in ConditionalTrialRemove_callable.m, since MATLAB
% local functions aren't callable across files - see that file if the
% two ever need to be reconciled.
answer = 'No';
fig = figure('Name', dlgtitle, 'NumberTitle', 'off', 'MenuBar', 'none', ...
    'ToolBar', 'none', 'WindowStyle', 'normal', 'Resize', 'off', ...
    'Position', [400 400 380 140]);
uicontrol(fig, 'Style', 'text', 'String', question, 'Units', 'normalized', ...
    'Position', [0.05 0.45 0.9 0.5], 'HorizontalAlignment', 'center');
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Yes', 'Units', 'normalized', ...
    'Position', [0.1 0.1 0.3 0.3], 'Callback', @(~,~) respond('Yes'));
uicontrol(fig, 'Style', 'pushbutton', 'String', 'No', 'Units', 'normalized', ...
    'Position', [0.55 0.1 0.35 0.3], 'Callback', @(~,~) respond('No'));
uiwait(fig);
    function respond(val)
        answer = val;
        uiresume(fig);
        close(fig);
    end
end
