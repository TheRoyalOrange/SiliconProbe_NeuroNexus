%% 
% This is a script to identify which channels/shanks on a probe are most
% active during spontaneous activity and during stimulus-driven activity
% (according to some set of stimuli chosen by the user. Ex. Visual stim).
% This allows then labelling of each shank or channel by these values for
% later location referencing. Info to be added to probeInfo variable. 


%%
animal = '20260423-p12';

save_directory = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

fq_range = [2 80]; %frequency band to check for oscillations in 
probes = [1]; %which probes in the recording to look at? 


reorder = 0;

%IN CASE REORDER IS 0: 
%load 'pre202406_64_32_channel_re_index'
%load 'pre202406_64_32_channel_re_index_v2'
%load '202408_64_64_channel_re_index_v5.mat'
%load 20240708_64_32_channel_re_index.mat



%% NEEDS INPUT FIRST: Clear the board and set the directory (change the directory path)
load(fullfile([save_directory '\'  animal '-ProbeInfo.mat']))
load_directory = 'E:\Roy\Silicon Probe Raw Data'; %don't change unless you're not Roy and have your data elsewhere

spondata_directory = uipickfiles('FilterSpec',load_directory,'Prompt', 'Choose Folders Containing Recordings you want to use [IMMEDIATE SUBFOLDERS SHOULD BE -Record Node ###-]'); %loads file names'E:\Roy\Ungrouped Whiskers\20251011-p11\20251011-p11-light4l22_40'; %for example

 % data_directory = {  'E:\Roy\Silicon Probe Raw Data\20260423-p12\20260423-p12-spon';
 %                     'E:\Roy\Silicon Probe Raw Data\20260423-p12\20260423-p12-whisker';
 %                     %'E:\Roy\Silicon Probe Raw Data\20260423-p12\20260423-p12-light22_40';
%                    'E:\Roy\Silicon Probe Raw Data\20260226-p12\20260226-p12-spon';
%                    'E:\Roy\Silicon Probe Raw Data\20260226-p12\20260226-p12-whisker';
%                    'E:\Roy\Silicon Probe Raw Data\20260226-p12\20260226-p12-light22_40';
                    %'E:\Roy\Ungrouped Whiskers\20251011-p11\20251011-p11-spon1';
                    %'E:\Roy\Ungrouped Whiskers\20251011-p11\20251011-p11-spon2';
                    %'E:\Roy\Ungrouped Whiskers\20240325-p6-raw\20240325-p6-spon';
                    %'E:\Roy\Ungrouped Whiskers\20240325-p6-raw\20240325-p6-whisker';
                    %'E:\Roy\Ungrouped Whiskers\20250308-p9\20250308-p9-light22_72';
                    %'E:\Roy\Ungrouped Whiskers\20250308-p9\20250308-p9-light22_72-n-whisker'
                       % 'E:\Roy\Silicon Probe Raw Data\20260326-p5\20260326-p5-spon';
                       % 'E:\Roy\Silicon Probe Raw Data\20260326-p5\20260326-p5-spon2';
                       % 'E:\Roy\Silicon Probe Raw Data\20260326-p5\20260326-p5-whisker';
                       % 'E:\Roy\Silicon Probe Raw Data\20260326-p5\20260326-p5-whisker10Hz';
              %     }; %for example

[oscillations, snippets, prbsnippets] = oshkosh_detect_funct(animal, spondata_directory, probes, fq_range, save_directory, reorder, ProbeInfo)
oscillation_plotEvents(ProbeInfo,oscillations, prbsnippets)

%% Measure total spontaneous activity of each channel

%compute RMSsum of each channel for each probe event
for prb = 1:length(probes)
    chRMSsum_total = zeros(length(prbsnippets.chanals),1);
    for snp =  1:length(prbsnippets.event_data{prb})
        if oscillations.prbEvent_during_stim{prb}(snp) == 0
        dat_snp = prbsnippets.event_data{prb}{snp};
        snp_sums = [];
        for ch = 1:size(dat_snp,1)
        snp_sums(ch) = sum(sqrt(movmean(dat_snp(ch,:).^2,200)),2);
        end
        chRMSsum_total = chRMSsum_total + snp_sums';
        else
        end
    end
    chActivityScore_Spon{prb} = chRMSsum_total;
end

for prb = 1:length(probes)
   for shnk = 1:size(ProbeInfo.ProbeMaps{prb},2)
      shnkActivityScore_Spon{prb}(shnk) = mean(chActivityScore_Spon{prb}(ismember(prbsnippets.chanals,ProbeInfo.ProbeMaps{prb}(:,shnk))))
   end
end

