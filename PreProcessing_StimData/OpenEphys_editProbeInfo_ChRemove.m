% OpenEphys_editProbeInfo_ChRemove.m
%
% Description: Downstream, standalone utility to record which channels
%   should be excluded from an animal's probe(s) (e.g. outside brain,
%   broken). Loads that animal's existing ProbeInfo.mat, overwrites its
%   .Ch_Remove field with a user-specified per-probe channel list
%   (badchans), and saves the struct back to the same file. Does not
%   itself alter any LFP/CSD/TF/spiking result files — it only records
%   which channels downstream analyses should treat as excluded.
%
% Inputs:
%   badchans (cell array, 1 x number of probes; each cell a column
%     vector of channel IDs) - PROVIDED BY USER (edited directly in the
%     script). Channels to exclude per probe, using ChanIds from the
%     existing ProbeInfo.
%   animal (string) - PROVIDED BY USER via input dialog at runtime.
%     Animal name matching the ProbeInfo file to edit.
%   [Note: ProbeInfo has prefix animal-]
%   ProbeInfo.mat (struct, loaded from disk) - struct containing several
%     fields with information relevant to the experiment, recording
%     details, and plotting of data. Expected to already exist (created
%     by OpenEphys_BaseAnalysis.m or a variant). Has the following
%     fields:
%
%       .Animal (string) - animal name used for file labels and figures
%       .Areas (cell array of strings, 1 x probenum) - which brain areas are associated with each probe
%       .ProbeNum (int) - number of probes in the recording
%       .ProbeMaps (cell array, 1 x probenum) - 2d matrices reflecting shape of each probe and associated channel IDs in data (IDs correspond to rows of raw data)
%       .ChanIds (int vector) - vector of channel IDs in order
%       .Chans (int) - number of recording channels
%       .TTLch (int) - channel ID of TTL trigger channel. should be chans+1 unless something wierd in the GUI during recording
%       .poi (int vector) - which probes in the recording are to be analyzed?
%       .ProbeIds (cell array, 1 x probenum) - channel IDs within each probe (starts at 1 : number of channels on probe)
%       .Ch_Remove (cell array) - channels to exclude per probe; edited by this script
%
% Outputs:
%   ProbeInfo.mat (overwritten in place) - same struct as above, with
%     .Ch_Remove (cell array, 1 x number of probes; each cell a vector
%     of channel IDs) set to badchans.
%
% Dependencies: Expects OpenEphys_BaseAnalysis.m (or the _Bundled /
%   _MixedTrials variant) to have already been run for this animal,
%   producing the ProbeInfo.mat file this script loads and overwrites.
%
%
%
%
%%
%first, open whatever figures you would like to use to make this decision.


%now write the channels you don't like. Use Chan_Ids from ProbeInfo.Chan_Ids
badchans = {
            [1, 2, 3, 9, 10, 11, 12 17, 18, 19, 20, 25, 26, 27, 28, 33, 34, 35, 36, 41, 42, 43, 44, 49, 50, 51, 52, 57, 58, 59, 60]' %channels from first probe
            [65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 81, 89, 97]'    %channels from second probe
                }';


%% load probeinfo and make changes
%get the animal name 
% (dialog box prevents editting wrong animal after forgetting to change)
animal = inputdlg(['Please specify the animal whos info you would like to edit'],'Animal check')


load(fullfile(['E:\Roy\Processed Silicon Probe Data\ProbeInfo\' char(animal) '-ProbeInfo.mat']))
%list channels to remove inside array for that probe. 


ProbeInfo.Ch_Remove = badchans;
disp(ProbeInfo.Ch_Remove)
%% save it

save(fullfile(['E:\Roy\Processed Silicon Probe Data\ProbeInfo\' char(animal) '-ProbeInfo.mat']),"ProbeInfo")

