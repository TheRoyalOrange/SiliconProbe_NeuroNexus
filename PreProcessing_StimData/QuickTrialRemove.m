% QuickTrialRemove.m
%
% Description: Downstream, standalone trial-curation utility. Lets the
%   user visually inspect LFP and spike-raster traces for one channel
%   of an already-preprocessed animal/condition (output of
%   OpenEphys_BaseAnalysis.m or a variant), manually choose trials to
%   exclude, and persist the updated tr_keep/tr_remove back into the
%   corresponding LFP.mat, spiking_results.mat, and TF_results.mat
%   files (via save(...,'-append')). Deliberately does NOT update
%   CSD_results.mat: CSD is stored as a trial-mean rather than
%   per-trial data, so removing trials there would require recomputing
%   the CSD (not done by this script).
%
% Inputs:
%   animal (string) - PROVIDED BY USER. Animal name matching the data
%     to examine (typed independently of filename as a safety check
%     against editing the wrong files).
%   filename (string) - PROVIDED BY USER. animal name + condition name,
%     matching the prefix of the LFP/spiking/TF result files to edit.
%   chan (int) - PROVIDED BY USER. Channel to inspect/plot.
%   fs (int) - PROVIDED BY USER. Sampling rate of the original
%     recording, in Hz (used to bin spikes into 1 ms bins).
%   stim_spike_stimchunks (double, trials x 30kHz trial length x
%     channels) and tr_keep (int vector, trials) - PROVIDED BY USER.
%     INFERRED DATA CONTRACT: expected to already be loaded into the
%     workspace from filename-spiking_results.mat (per the script's own
%     instruction to "double click on a spike results mat file to open
%     it" before running) — this script does not load them itself,
%     unlike stim_lfp_stimchunks below.
%   stim_lfp_stimchunks - loaded automatically by this script from
%     filename-LFP.mat. 3D array of lfp data split into trials. Data is
%     downsampled to 1kHz from original 30kHz in the recording.
%     size[trials x trial length x channels]
%
% Outputs:
%   tr_keep, tr_remove - updated trial-inclusion vectors, saved
%     (via -append) into filename-LFP.mat, filename-spiking_results.mat,
%     and filename_TF_results.mat. Not written to CSD_results.mat (see
%     Description).
%
% Dependencies: Expects OpenEphys_BaseAnalysis.m (or the _Bundled /
%   _MixedTrials variant) to have already been run for this
%   animal/condition, producing the LFP.mat, spiking_results.mat, and
%   TF_results.mat files this script loads/appends to.
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

