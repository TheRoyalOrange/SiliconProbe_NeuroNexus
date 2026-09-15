%% Notes for my dearest Claude
% 
%  INPUTS:
%        animal - PROVIDED BY USER. string providing animal name matching
%                 that of the animal who's data you'd like to examine 
%                 (at some point should probably just automatically get this from filename) 
%        filename - an extension of animal, now with condition name
%                   appended. Should match the first part of whatever animal and
%                   condition data you'd like to examine
%
%
%      [Note: the following input has prefix filname-]
%       
%       spiking_results.mat - mat file containing:
%                           stim_spike_stimchunks - 3D array of MUA data,
%                                                   where each timepoint is labelled as 0 or 1, for
%                                                   no spike or spike. Data is in original 30kHz
%                                                   sampling rate.
%                                                   size[trials x 30kHz trial length x channels]
%                                                 
%                           tr_keep  - vector of trials considered "good" according to 
%                                      user. By default, this contains all trials and
%                                      is modified later by this script
%                                      size[trials]
%                           tr_remove - the complement to tr_keep. By default this is
%                                       left empty initially and modified later by this
%                                       script size[empty]
% 
%       stim_lfp_stimchunks - loaded from LFP.mat file with prefix
%                             filename. 3D array of lfp data split into
%                             trials. Data is downsampled to 1kHz from original
%                             30kHz in the recording. size[trials x trial length x channels] 
% 
% 
%
% OUTPUTS:  
% 
%    tr_remove - modified from previous state and saved in new form for all
%                corresponding LFP.mat, TF_results.mat, and spiking_results.mat
%                files with prefix = filename. (Does not change in CSD_results.mat 
%                because the data there is saved as mean across trials and
%                would have to be recalculated accordingly. Probably should
%                fix that sometime). 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

%%
%first, double click on a spike results mat file to open it

animal = '20260423-p12'; %type this in independently to filename, just to ensure you don't change the wrong files
filename = '20260423-p12-whisker';%copy and paste the animal and condition part of the file you opened
fs = 30000; %sampling rate of original recording
chan = [13]; %the channel you'd like to observe the data for


spikes_ms = [];
for ch = 1:size(stim_spike_stimchunks,3)
    
    
    %chspikes = squeeze(stim_spike_stimchunks(tr_keep,:,ch));
    chspikes = squeeze(stim_spike_stimchunks(:,:,ch));

    
        for batch = 1:size(stim_spike_stimchunks,2)/(fs/1000)
            
            spikebatchi = sum(chspikes(:,1+((fs/1000)*(batch-1)):(fs/1000)+((fs/1000)*(batch-1))),2)>0;
            spikeper(:,batch) = spikebatchi;
        end
        spikes_ms(:,:,ch) = spikeper;
        stim_spike_avgrate(ch,:) = mean(squeeze(movmean(spikeper,50,2)*1000),1);
        
   
end 

load(fullfile(['E:\Roy\Processed Silicon Probe Data\LFP\' animal '\' filename '-LFP.mat']),'stim_lfp_stimchunks')
%%
time = [-500:3000]
figure();
for tr = 1:size(stim_spike_stimchunks,1)
    subplot(1,2,1)
    hold on 
    colormap(sky(8))
    plot(time,squeeze(stim_lfp_stimchunks(tr,4500:8000,chan))'-(400*(tr-1)))
%    plot(time,squeeze(stim_lfp_stimchunks(tr,4500:8000,chan))'-(400*(tr-1)),Color=[.5 0 .3])
    xline(0)
    subplot(1,2,2)
    hold on 
    plot(time,squeeze(spikes_ms(tr,4500:8000,chan))'-(2*(tr-1)))
    xline(0)
end
%%
tr_remove = []; %type it in 
tr_keep(tr_remove) = []; %change this one
%%


%save that in all the right spots
save(fullfile(['E:\Roy\Processed Silicon Probe Data\LFP\' animal '\' filename '-LFP.mat']), 'tr_keep','tr_remove','-append')
save(fullfile(['E:\Roy\Processed Silicon Probe Data\Spiking\' animal '\' filename '-spiking_results.mat']), 'tr_keep','tr_remove','-append')
%save(fullfile(['E:\Roy\Processed Silicon Probe Data\CSD\' animal '\' filename '-CSD_results.mat']), 'tr_keep','tr_remove','-append')
save(fullfile(['E:\Roy\Processed Silicon Probe Data\TF\' animal '\' filename '_TF_results.mat']), 'tr_keep','tr_remove','-append')

