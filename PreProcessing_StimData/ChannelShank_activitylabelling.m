%% 
% This is a script to identify which channels/shanks on a probe are most
% active during spontaneous activity and during stimulus-driven activity
% (according to some set of stimuli chosen by the user. Ex. Visual stim).
% This allows then labelling of each shank or channel by these values for
% later location referencing. Info to be added to probeInfo variable. 


%%
animal = '20260423-p12';
probes = [1]; %which probes in the recording to look at? 


save_directory = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%fq_range = [2 80]; %frequency band to check for oscillations in 

load(fullfile([save_directory '\' animal '-ProbeInfo.mat']))


%% Measure strongest stimulated shank (by light)

%get lfp files to use
stimdata_directory = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\LFP','Prompt', 'Choose Folders Containing Recordings you want to use [IMMEDIATE SUBFOLDERS SHOULD BE -Record Node ###-]'); %loads file names'E:\Roy\Ungrouped Whiskers\20251011-p11\20251011-p11-light4l22_40'; %for example
chanals = ~ismember(ProbeInfo.Ch_Remove{fisprb},ProbeInfo.ProbeMaps{prb})
%load data (from relevant channels)
stimdat = [];
for file = 1:length(stimdata_directory)
    d = load(stimdata_directory{file},'stim_lfp_stimchunks');
    stimdat = cat(1,stimdat, d.stim_lfp_stimchunks(:,:,prbsnippets.chanals));
end
clear d

%get RMSsums, baselined by prestim activity
for tr = 1:size(stimdat,1)
    rms = sqrt(movmean(squeeze(stimdat(tr,:,:)).^2,200));
    baseline = mean(rms(3000:4999,:),1);
    baselined_rms = rms - repmat(baseline,size(rms,1),1);
    stimdat_RMS(tr,:) = sum(baselined_rms(5000:8000,:));
end

%get channel and shank specificy values
for prbi = 1:length(probes)
    prb = probes(prbi)
    chActivityScore_Stim{prbi} = mean(stimdat_RMS(:,~ismember(ProbeInfo.Ch_Remove{prb},ProbeInfo.ProbeMaps{prb}),1));
    for shnk = 1:size(ProbeInfo.ProbeMaps{prb},2)
      shnkActivityScore_Stim{prbi}(shnk) = mean(chActivityScore_Stim{prb}(~ismember(ProbeInfo.Ch_Remove{prb},ProbeInfo.ProbeMaps{prb}(:,shnk))));
   end
end




