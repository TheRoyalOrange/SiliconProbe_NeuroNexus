% ConditionalTrialRemove_callable_test.m
%
% Description: Drop-in test double for ConditionalTrialRemove_callable.m,
%   for dry-running ConditionalTrialRemove.m's looping/wiring
%   logic without touching any real LFP/Spiking/CSD/TF result files.
%   Shares the real function's call signature exactly, so
%   ConditionalTrialRemove.m can swap one call for the other with
%   no other changes (see the commented-out line in that script). Does
%   NOT load or require real preprocessed data: the first time a given
%   animal/condition pair is seen, it seeds mock tr_keep/tr_remove/
%   tr_remove_conditional (same shapes/schema
%   ConditionalTrialRemove_callable.m itself falls back to for an
%   old-format file - see that file's header) instead of reading them
%   from a real spiking_results.mat. There is no real LFP/spike data to
%   visually inspect in test mode, so unlike the real function this
%   never opens figures or a trial-picking UI - on a "no match" it just
%   shows a notice and auto-records a "keep all trials" (TrialIdx = NaN)
%   row, mirroring the real function's own NaN convention for "evaluated,
%   nothing removed". This keeps the exact same skip-if-already-recorded
%   behavior testable across repeated calls/runs, since state is
%   persisted (see Outputs), while removing the parts of the real
%   pipeline that only make sense against real data.
%
% Inputs:
%   animal, condition, tr_conditional, save_directory, master_filename -
%     identical meaning to ConditionalTrialRemove_callable.m's own inputs
%     of the same names - see that file's header.
%   fs, ProbeInfo - accepted ONLY for call-signature parity with
%     ConditionalTrialRemove_callable.m (so ConditionalTrialRemove.m's
%     loop can call either interchangeably). Neither is used: there is no
%     real spike data to bin (fs) or real LFP figures to open
%     (ProbeInfo.Areas) in test mode. ProbeInfo is returned unchanged.
%
% Outputs:
%   ProbeInfo (struct or []) - passed through unchanged from the input,
%     for call-signature parity only.
%   <save_directory>\<animal>-<condition>-results_test.mat (side effect,
%     not a return value) - holds this animal/condition pair's mock
%     tr_keep/tr_remove/tr_remove_conditional (same three variables a
%     real spiking_results.mat carries), created on first use and
%     updated on every subsequent call, in place of the four real
%     LFP/spiking/CSD/TF result files. Saved FLAT at save_directory's
%     root (not nested under Spiking\<animal>\), same placement as
%     master_filename, for easy side-by-side inspection.
%   <save_directory>\<master_filename> (side effect, not a return value)
%     - updated exactly as ConditionalTrialRemove_callable.m's own
%     updateTrConditionalMaster does (duplicated here rather than shared,
%     so this test file has zero dependency on - and can never
%     accidentally call into - the real function).
%
% Dependencies: none (standalone; does not call
%   ConditionalTrialRemove_callable.m).

function ProbeInfo = ConditionalTrialRemove_callable_test(animal, condition, tr_conditional, fs, ProbeInfo, save_directory, master_filename) %#ok<INUSD> fs unused, see header

if nargin < 7 || isempty(master_filename)
    master_filename = 'tr_conditional_master.mat';
end

if nargin < 6 || isempty(save_directory)
    save_directory = 'E:\Roy\Processed Silicon Probe Data';
end

full_filename = [animal '-' condition];
testFile = fullfile(save_directory, [full_filename '-results_test.mat']);

if isfile(testFile)
    test_dat = load(testFile, 'tr_keep', 'tr_remove', 'tr_remove_conditional');
    tr_keep = test_dat.tr_keep;
    tr_remove = test_dat.tr_remove;
    tr_remove_conditional = test_dat.tr_remove_conditional;
else
    % first time this animal/condition pair is seen in test mode - seed mock
    % trial-tracking data (10 fake trials, none removed yet)
    tr_keep = 1:10;
    tr_remove = zeros(1, 10);
    tr_remove_conditional = table('Size', [0,2], 'VariableTypes', {'string','double'}, 'VariableNames', {'Name','TrialIdx'});
end

%% if this named condition already has a recorded trial set for this (fake) file, skip - same behavior as the real function
if any(strcmp(tr_remove_conditional.Name, tr_conditional))
    fprintf('[TEST] Match found for ''%s'' in %s - skipping.\n', tr_conditional, full_filename);
    updateTrConditionalMasterTest(save_directory, master_filename, tr_conditional, animal, full_filename);
    return
end

uiwait(msgbox(sprintf(['[TEST MODE] No match found for ''%s'' in %s.\n\nNo real trial data to inspect in test mode - ' ...
    'automatically recording ''keep all trials'' for this condition.'], tr_conditional, full_filename)));

%% record a "keep all trials" (TrialIdx = NaN) row, same convention the real function uses for zero trials removed
newRow = table(string(tr_conditional), NaN, 'VariableNames', {'Name','TrialIdx'});
tr_remove_conditional = [tr_remove_conditional; newRow];

fprintf('\n--- [TEST] ConditionalTrialRemove_callable_test: row added for %s, condition ''%s'' ---\n', full_filename, tr_conditional);
disp(newRow);

save(testFile, 'tr_keep', 'tr_remove', 'tr_remove_conditional');
fprintf('[TEST] wrote %s\n', testFile);

updateTrConditionalMasterTest(save_directory, master_filename, tr_conditional, animal, full_filename);

end

function updateTrConditionalMasterTest(save_directory, master_filename, tr_conditional, animal, full_filename)
% Identical logic to updateTrConditionalMaster (local function in
% ConditionalTrialRemove_callable.m) - duplicated here, rather than
% shared, so this test file has no dependency on the real callable and
% can't accidentally touch real result files. See that file if the two
% ever need to be reconciled.
masterFile = fullfile(save_directory, master_filename);
if isfile(masterFile)
    master_dat = load(masterFile, 'tr_conditional_master');
    tr_conditional_master = master_dat.tr_conditional_master;
else
    tr_conditional_master = table('Size', [0,3], 'VariableTypes', {'string','string','string'}, 'VariableNames', {'Name','Animal','Filename'});
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

fprintf('[TEST] %s: added %s / %s / %s\n', master_filename, tr_conditional, animal, full_filename);
save(masterFile, 'tr_conditional_master');

end
