%% Notes for Claude (and you, of course, dear reader)
% INPUTS: 
%   [Note: ProbeInfo has prefix animal-]
%       
%       animal - PROVIDED BY USER. string providing animal name matching
%                that of the probeinfo file you want to edit 
%   
%       ProbeInfo.mat - struct containing several fields with information
%                       relevant to the experiment, recording details, and plotting of data.
%                       It is initialized in this script and can be further updated later. Has
%                       the following fields:
%                      
%                       .Animal = animal; %animal name used for file labels and figures
%                       .Areas = areas; %which brain areas are associated with each probe
%                       .ProbeNum = probenum; %number of probes in the recording
%                       .ProbeMaps = probemaps; %2d matrices reflecting shape of each probe and associated channel IDs in data (IDs correspond to rows of raw data)
%                       .ChanIds = chan_ids; %vector of channel IDs in order
%                       .Chans = chans; %number of recording channels
%                       .TTLch = TTLch; %channel ID of TTL trigger channel. should be chans+1 unless something wierd in the GUI during recording
%                       .poi = poi; %which probes in the recording are to be analyzed?
%                       .ProbeIds = probeids; %channel IDs within each probe (starts at 1 : number of channels on probe)
%                       .Ch_Remove = {}; %empty field to be modified later using OpenEphys_editProbeInfo_ChRemove in case some channels are not useful (e.g. outside brain/broken) 
%  
% 
%  OUTPUTS: 
%
%       ProbeInfo.Ch_Remove - an edited .Ch_Remove, with the values specified in badchans. Size {1 x number of probes}[number of channels listed for probe] 
% 
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

