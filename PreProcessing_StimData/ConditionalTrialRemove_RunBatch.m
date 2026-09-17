% ConditionalTrialRemove_RunBatch.m
%
% Description: Standalone driver script that batch-applies
%   ConditionalTrialRemove_callable.m across a manually chosen set of
%   animals and conditions, under a single tr_conditional name. Meant to
%   be run start-to-finish (F5) rather than section-by-section: it
%   blocks on its own dialogs (confirm tr_conditional, confirm the
%   animal/condition selections) before looping. As close to the very
%   start of the script as possible (right after the user-set variables,
%   before anything else runs), opens a non-modal uitable of the
%   existing tr_conditional_master.mat index so the user can scan what
%   conditional names already exist before committing to tr_conditional.
%   After the user confirms tr_conditional, checks whether that name
%   already exists in the master index and, if so, warns and asks
%   whether to continue (since re-running the same name on files it has
%   already been applied to is a no-op - ConditionalTrialRemove_callable.m
%   skips any animal/filename pair that already has a recorded trial set
%   for that name, and there is currently no way to edit an
%   already-recorded assignment). For each animal in animal_inc, opens a
%   uipickfiles dialog rooted at that animal's LFP results folder so the
%   user can pick which condition files to run tr_conditional against;
%   the chosen files are parsed into condition names (stripping the
%   <animal>- prefix and -LFP suffix) and handed to
%   ConditionalTrialRemove_callable.m one at a time.
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
%   (side effect) user picks, per animal in animal_inc, which condition
%     files under <save_directory>\LFP\<animal> to include, via
%     uipickfiles. INFERRED DATA CONTRACT: each selected file is named
%     <animal>-<condition>-LFP.mat (matching how
%     OpenEphys_BaseAnalysis.m, or a variant, actually names LFP.mat
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
%   (side effect, via ConditionalTrialRemove_callable.m) tr_remove_conditional
%     rows written into each selected animal/condition's LFP/spiking/
%     CSD/TF result files, and tr_conditional_master.mat updated - see
%     ConditionalTrialRemove_callable.m's own header for details.
%
% Dependencies: ConditionalTrialRemove_callable.m (same folder);
%   uipickfiles (must be on the MATLAB path - already used elsewhere in
%   this repo, e.g. Summarize_MultiAnimal\SummaryAnalysis_MultiMouseTF.m);
%   expects OpenEphys_BaseAnalysis.m (or a variant) to have already been
%   run for every animal/condition selected.

%% ============ USER-SET VARIABLES - fill these in by hand ============
tr_conditional = '';   % name of the trial-exclusion condition to apply, e.g. 'whisker_twitch_artifact'
animal_inc = {         % animals to include this run
    ''
    };
fs = 30000;             % sampling rate of the original recording (Hz)
save_directory = 'E:\Roy\Processed Silicon Probe Data';
%% =====================================================================

%% show the existing tr_conditional_master index so the user can scan before proceeding
displayTrConditionalMaster(save_directory);

%% now validate that the user actually filled in the required variables above
if isempty(tr_conditional) || ~(ischar(tr_conditional) || isstring(tr_conditional))
    error('ConditionalTrialRemove_RunBatch:EmptyConditional', 'Fill in tr_conditional before running.');
end
if isempty(animal_inc) || ~iscell(animal_inc) || any(cellfun(@isempty, animal_inc))
    error('ConditionalTrialRemove_RunBatch:EmptyAnimalInc', 'Fill in animal_inc (cell array of animal names) before running.');
end

%% confirm tr_conditional with the user
confirmMsg = sprintf('tr_conditional is set to:\n\n''%s''\n\nHappy with this?', tr_conditional);
answer = nonModalConfirm(confirmMsg, 'Confirm tr_conditional');
if ~strcmp(answer, 'Yes')
    disp('tr_conditional not confirmed - edit the variable at the top of the script and re-run.');
    return
end

%% warn if this tr_conditional name already exists in the master index
masterFile = fullfile(save_directory, 'tr_conditional_master.mat');
if isfile(masterFile)
    master_dat = load(masterFile, 'tr_conditional_master');
    existingRows = strcmp(master_dat.tr_conditional_master.Name, tr_conditional);
    if any(existingRows)
        n = sum(existingRows);
        existMsg = sprintf(['tr_conditional ''%s'' already exists in tr_conditional_master.mat ' ...
            '(%d animal/filename pair(s) already recorded).\n\nNote: already-recorded pairs cannot ' ...
            'be re-edited - ConditionalTrialRemove_callable.m will just skip them.\n\nContinue anyway?'], ...
            tr_conditional, n);
        answer = nonModalConfirm(existMsg, 'tr_conditional already exists');
        if ~strcmp(answer, 'Yes')
            disp('Aborted by user - tr_conditional already exists in the master index.');
            return
        end
    end
end

%% for each animal, pick which condition files (under its LFP folder) to run this tr_conditional against
animal_condition_files = cell(1, numel(animal_inc));
for a = 1:numel(animal_inc)
    animal = animal_inc{a};
    animalLFPdir = fullfile(save_directory, 'LFP', animal);
    if ~isfolder(animalLFPdir)
        warning('ConditionalTrialRemove_RunBatch:AnimalFolderNotFound', ...
            'LFP folder not found for animal ''%s'': %s - skipping.', animal, animalLFPdir);
        animal_condition_files{a} = {};
        continue
    end

    picked = uipickfiles('FilterSpec', animalLFPdir, 'Prompt', ...
        sprintf('Choose ''%s'' condition file(s) for %s', tr_conditional, animal));

    if ~iscell(picked) % uipickfiles returns 0 if the user cancels
        warning('ConditionalTrialRemove_RunBatch:NoFilesPicked', ...
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
            warning('ConditionalTrialRemove_RunBatch:UnexpectedFilename', ...
                ['Picked file ''%s'' does not match the expected <animal>-<condition>-LFP.mat ' ...
                 'pattern for animal ''%s'' - skipping this file.'], base, animal);
            conditions{f} = '';
        end
    end
    animal_condition_files{a} = conditions(~cellfun(@isempty, conditions));
end

%% run ConditionalTrialRemove_callable.m for every animal/condition selected above
for a = 1:numel(animal_inc)
    animal = animal_inc{a};
    conditions = animal_condition_files{a};
    ProbeInfo = []; % reloaded once per animal (via ConditionalTrialRemove_callable.m's own optional ProbeInfo input/output), then reused across this animal's conditions instead of re-loading per condition
    for c = 1:numel(conditions)
        ProbeInfo = ConditionalTrialRemove_callable(animal, conditions{c}, tr_conditional, fs, ProbeInfo, save_directory);
    end
end

disp('ConditionalTrialRemove_RunBatch: done.');

%% ------------------------- local functions -------------------------

function displayTrConditionalMaster(save_directory)
% Opens a non-modal uitable of tr_conditional_master.mat (sorted by
% Name, then Animal, then Filename) so the user can scan which
% tr_conditional names already exist, and on which animal/filename
% pairs, before committing to a tr_conditional above. Left open (not
% uiwait'd) so it stays visible/referenceable through the rest of the
% script's dialogs.
masterFile = fullfile(save_directory, 'tr_conditional_master.mat');
if ~isfile(masterFile)
    disp('No tr_conditional_master.mat found yet - this will be the first tr_conditional recorded.');
    return
end
master_dat = load(masterFile, 'tr_conditional_master');
sortedMaster = sortrows(master_dat.tr_conditional_master, {'Name','Animal','Filename'});

fig = figure('Name', 'tr_conditional_master (existing conditionals)', 'NumberTitle', 'off', ...
    'MenuBar', 'none', 'ToolBar', 'none', 'Position', [200 200 640 400]);
uitable(fig, 'Data', table2cell(sortedMaster), 'ColumnName', sortedMaster.Properties.VariableNames, ...
    'Units', 'normalized', 'Position', [0 0 1 1]);
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
