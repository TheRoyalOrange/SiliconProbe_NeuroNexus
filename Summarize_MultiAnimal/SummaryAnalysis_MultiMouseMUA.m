% SummaryAnalysis_MultiMouseMUA.m
%
% Description: Pools peri-stimulus MUA across mice and conditions for a
%   user-chosen set of channels (e.g. one channel per cortical layer). Spikes
%   are binned to 1 ms and smoothed into a firing rate (50 ms moving mean),
%   then per-trial AUC and peak rate are computed for an early (P1) and late
%   (P2) window after the stimulus. Optional sections plot condition-averaged
%   rate timeseries. Results go into one long-format table (one row per
%   trial x channel) and are written to a csv for stats/plotting in R.
%
% Inputs:
%   User-set config (edited at top of script):
%     mice (cell array of strings, nMice x 1) - animal IDs; order becomes the
%       factor order in R. List a mouse twice to analyze two channel groups
%     chans (double, nChanPerMouse x nMice) - channel IDs (rows of raw data)
%       to analyze for each mouse; row i in every column is the same
%       layer/group, so it is compared across mice
%     chan_groups (string, 1 x nChanPerMouse) - label for each row of chans (e.g. "L2/3")
%     shank_groups (string, 1 x nChanPerMouse) - shank label for each row of chans (e.g. "C")
%     region (cell array of strings, 1 x nMice) - brain area recorded per mouse (e.g. "V1")
%     condis (cell array of char, 1 x nCond) - condition names; order becomes the factor order in R
%     condiFileMice (cell array, nCond x 1) - each cell is a vector (1 x nFiles in
%       that condition) giving the index into mice for each file picked in
%       superCondis_dir{con}, in pick order
%     congrouped, grouptitles (plotting sections only) - which condition
%       indices to overlay in each figure, and the figure titles
%   Files picked per condition via uipickfiles (superCondis_dir, cell 1 x nCond),
%   from E:\Roy\Processed Silicon Probe Data\Spiking\...\*-spiking_results.mat, each containing:
%     stim_spike_stimchunks (0/1, trials x samples x channels) - peri-stimulus
%       spike indicator at 30 kHz. Stim onset is assumed at 5000 ms into the
%       trial [inferred: P1/P2 windows and plot time axis assume this]
%     tr_keep (double, 1 x nKeep) - indices of good trials
%     tr_remove (double, 1 x nRemove) - indices of rejected trials
%       NOTE: only their lengths are used; ALL trials
%       (1:numel(tr_keep)+numel(tr_remove)) are analyzed, rejected ones included
%
% Outputs:
%   superCondis (cell, 1 x nCond) - each (trials x ms x channels) 1 ms binned spike indicator
%   superCondis_rate (cell, 1 x nCond) - each (trials x ms x channels) firing rate, in spikes/s
%   superCondis_avg, superCondis_std (double, nCond x channels x ms) - trial mean / SD of rate, in spikes/s
%   tableres (table, (nTrials*nChan) x 15) - one row per trial x channel. Columns:
%     P1_AUC, P2_AUC, All_AUC (double) - integral of the rate over P1 (ms 5075:5350),
%       P2 (5350:8000) and both (5075:8000), in spikes/s*ms
%     P1_Peak, P2_Peak (double) - max rate in the P1 / P2 window, in spikes/s
%     Animal_Name (string), Animal_Num (double) - mouse ID and its index in mice
%     Condition_Name (string), Condition_Num (double) - condition name and its index in condis
%     Condition_FileNum (double) - index of the source file within its condition
%     Region (string) - from region
%     Channel_ID (double) - raw-data channel ID (from chans)
%     Channel_GroupName (string) - from chan_groups
%     Shank_GroupName (string) - from shank_groups
%     Channel_GroupNum (double) - row index into chans (1..nChanPerMouse)
%   <filename>.csv - tableres written to
%     E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\
%     (file name chosen by the user in a dialog; existing names are rejected)
%
% Dependencies: OpenEphys_BaseAnalysis*.m (writes *-spiking_results.mat);
%   optionally QuickTrialRemove.m / ConditionalTrialRemove_callable.m
%   (update tr_keep/tr_remove). Requires uipickfiles (File Exchange).

%list mice to be analyzed. Note the order as it will be treated as a factor (R style)
%if mouse has more than one group of channels to analyze, list
%it twice here
mice = {"20260508-p9"
        "20260509-p10"
        "20260821-p10"
        "20260510-p11"
        
}; 
%mice = {"20260423-p12"}; 

%list channels to use for each mouse listed (ex. groups of three here for cortical
%layers. One group per animal, in same order as 'mice'). 
%Each channel is compared with its corresponding channel in other mice, so
%each entry should have the same number of channels per mouse
chans = [11	12 13;
         28	29	30;
         27 28 29;
         21 22 23
         ]'; 
%[12 13 14 20 21 22 28 29 30; %mouse 1
%         4 5 6 12 13 14 20 21 22;
%         21 22 23 29 30 31 37 38 39%mouse 2, etc
%         ]'; 
%[12 13 14 20 21 22 28 29 30]';

%give a name to each channel per mouse (should match size of dim1 of chans)
chan_groups = ["L2/3", "L4", "L5"]%,"L2/3", "L4", "L5","L2/3", "L4", "L5"];
shank_groups = ["C","C", "C"];% ["L","L","L", "C","C","C","R","R","R"];
shank_groups_types = unique(shank_groups,'stable');
%check num of chan_groups and chans per mouse match
size(chan_groups,2) == size(chans,1) && size(chans,1) == size(shank_groups,2)

%which brain area(s) are being recorded by the channels? (order should
%match dim2 of chans and length of mice
region = {"V1","V1","V1","V1"};

condis = {'W','L','LW', 'W-TTX', 'L_TTX','LW-TTX'}; %list conditions to be included, named as you'd prefer. Note the order

%{'W','L_4','LW_4', 'L_8','LW_8','L_12','LW_12','L_15','LW_15'}; %list conditions to be included, named as you'd prefer. Note the order
% as they will be treated as factors later  
%%
for con = 1:length(condis)
    superCondis_dir{con} = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\Spiking','Prompt', ['Choose ' condis{con} ' files']); 
end
%%

%now, manually fill in an array with numbers 1 to numberofmice matching
%   order of files in superCondis_dir{con}, and repeat for all conditions
% Ex. if the first file in a condition corresponds to the first listed mouse in 'mice',
%   then the first entry in condiFileMice is 1
condiFileMice = {[1 1 2 2 3 3 3 3 3 4 4];%condition 1 files
                 [1 1 2 2 3 3 3 3 3 4 4];%condition 2 files
                 [1 1 2 2 3 3 3 3 3 4 4];
                 [1 2 3 3 3 3 3 4 4]; %condition 3
                 [1 1 2 2 3 3 3 3 3 4 4]; 
                 [1 1 2 2 3 3 3 3 3 4 4];   
                 


                 };

%check num of files and filemice match
for con = 1:length(condis)
    size(superCondis_dir{con},2) == length(condiFileMice{con})
end
%%
%load the data you will be using (and only that data)
 superCondis = [];
 superCondis_rate = [];
 superCondis_avg = [];
 superCondis_std = [];

triallabel_condiname = []; %metadata for later
triallabel_condinum = []; 
triallabel_condifile = []; 
triallabel_animalnum = [];
triallabel_animalname = [];
triallabel_region = [];
triallabel_chanID = [];
triallabel_changroup = [];
triallabel_shankgroup = [];
triallabel_channum = [];


for con = 1:numel(condis)

    superCondis_contrials = [];
  
    triallabel_condiname_con = []; 
    triallabel_condinum_con = []; 
    triallabel_condifile_con = []; 
    triallabel_animalnum_con = [];
    triallabel_animalname_con = [];
    triallabel_region_con = [];
    triallabel_chanID_con = [];
    triallabel_changroup_con = [];
    triallabel_shankgroup_con = [];
    triallabel_channum_con = [];
    
    for file = 1:numel(superCondis_dir{con})
        disp(['Running Condition ', num2str(con), ' (', condis{con}, '), File ', num2str(file)])
        dat = matfile(strjoin(superCondis_dir{con}(file)));
        %trs = 1:25;
        trs =  1:length(dat.tr_keep)+length(dat.tr_remove);
        %trs = dat.tr_keep;
               
        triallabel_condiname_con = cat(1,triallabel_condiname_con,repmat(repmat(condis{con},length(trs),1),length(chan_groups),1));
        triallabel_condinum_con = cat(1,triallabel_condinum_con,repmat(repmat(con,length(trs),1),length(chan_groups),1));
        triallabel_condifile_con = cat(1,triallabel_condifile_con,repmat(repmat(file,length(trs),1),length(chan_groups),1));
        triallabel_animalnum_con = cat(1,triallabel_animalnum_con,repmat(repmat(condiFileMice{con}(file),length(trs),1),length(chan_groups),1));
        triallabel_animalname_con = cat(1,triallabel_animalname_con,repmat(repmat(mice{condiFileMice{con}(file)},length(trs),1),length(chan_groups),1));
        triallabel_region_con = cat(1,triallabel_region_con,repmat(repmat(region{condiFileMice{con}(file)},length(trs),1),length(chan_groups),1));
        
        
        triallabel_chanID_con = cat(1,triallabel_chanID_con,reshape(repmat(chans(:,condiFileMice{con}(file)),1,length(trs))',1,[])');
        triallabel_changroup_con = cat(1,triallabel_changroup_con,reshape(repmat(chan_groups',1,length(trs))',1,[])');
        triallabel_shankgroup_con = cat(1,triallabel_shankgroup_con,reshape(repmat(shank_groups',1,length(trs))',1,[])');
        triallabel_channum_con = cat(1,triallabel_channum_con,reshape(repmat((1:size(chans,1))',1,length(trs))',1,[])');
        
        superCondis_chantrials = [];
        for shnk = 1:length(shank_groups_types)
              superCondis_chantrials = cat(3,superCondis_chantrials, ...
                  dat.stim_spike_stimchunks(trs,:,chans(shank_groups == shank_groups_types(shnk),condiFileMice{con}(file))));
        end
        superCondis_contrials = cat(1,superCondis_contrials, superCondis_chantrials);
    end
    triallabel_condiname = cat(1,triallabel_condiname,string(triallabel_condiname_con)); 
    triallabel_condinum = cat(1,triallabel_condinum,triallabel_condinum_con);  
    triallabel_condifile = cat(1,triallabel_condifile,triallabel_condifile_con); 
    triallabel_animalnum = cat(1,triallabel_animalnum,triallabel_animalnum_con); 
    triallabel_animalname = cat(1,triallabel_animalname,triallabel_animalname_con); 
    triallabel_region = cat(1,triallabel_region,string(triallabel_region_con)); 
    triallabel_chanID = cat(1,triallabel_chanID,triallabel_chanID_con); 
    triallabel_changroup = cat(1,triallabel_changroup,triallabel_changroup_con); 
    triallabel_shankgroup = cat(1,triallabel_shankgroup,triallabel_shankgroup_con); 
    triallabel_channum = cat(1,triallabel_channum,triallabel_channum_con); 

    spikes_ms = [];
    spikes_rate = [];
    spikes_rateavg = [];
    spikes_ratestd = [];
    spikeper=[];
    for ch = 1:size(superCondis_contrials,3)
        
        chspikes = squeeze(superCondis_contrials(:,:,ch));
        fs = 30000;
        
        for batch = 1:size(superCondis_contrials,2)/(fs/1000)
            
            spikebatchi = sum(chspikes(:,1+((fs/1000)*(batch-1)):(fs/1000)+((fs/1000)*(batch-1))),2)>0;
            spikeper(:,batch) = spikebatchi;
        
        end

        spikes_ms(:,:,ch) = spikeper;
        spikes_rate(:,:,ch) = squeeze(movmean(spikeper,50,2)*1000);
        spikes_rateavg(ch,:) = mean(squeeze(spikes_rate(:,:,ch)),1);
        spikes_ratestd(ch,:) = std(squeeze(spikes_rate(:,:,ch)),1);
    
   
    end 
    superCondis{con} = spikes_ms;
    superCondis_rate{con} = spikes_rate;
    superCondis_avg(con,:,:) = spikes_rateavg;
    superCondis_std(con,:,:) = spikes_ratestd;
end



%calculate auc and peak for p1 and p2 of response

for con = 1:numel(condis)
    p1auc_con = [];
    p2auc_con = [];
    allauc_con = [];
    p1peak_con = [];
    p2peak_con = [];
    p1peaktimes_con = [];
    

    for ch = 1:size(chans,1)

        %AUC
        auc = cumtrapz(superCondis_rate{con}(:,5075:5350,ch),2);
        p1auc_con = cat(1,p1auc_con,auc(:,end));
        auc = cumtrapz(superCondis_rate{con}(:,5350:8000,ch),2);
        p2auc_con = cat(1,p2auc_con,auc(:,end));
        auc = cumtrapz(superCondis_rate{con}(:,5075:8000,ch),2);
        allauc_con = cat(1,allauc_con,auc(:,end));

        %Peaks
        [peaks, tps] = max(superCondis_rate{con}(:,5075:5350,ch),[],2);
        p1peak_con = cat(1,p1peak_con,peaks(:,end));
        p1peaktimes_con = cat(1,p1peaktimes_con,tps(:,end));

        peaks = max(superCondis_rate{con}(:,5350:8000,ch),[],2);
        p2peak_con = cat(1,p2peak_con,peaks(:,end));
        %peaks = max(superCondis_rate{con}(:,5020:8000),[],2);
        %allpeak{con} = peaks;

    end
    
    p1auc{con}  = p1auc_con;
    p2auc{con}  = p2auc_con;
    allauc{con} = allauc_con;
    p1peak{con} = p1peak_con;
    p1peaktimes{con} = p1peaktimes_con;
    p2peak{con} = p2peak_con;
end



%% plot spiking timeseries of channels
wind = [4900:8000]';
yends = [0 500];
time = [-100:3000]';
smoothfactor = 100;

%colors for each group within a plot
colors = {'b' 'k' 'r'};

%which groups to plot together
congrouped = {[2,3];
              [2,5];
              [3,6];
              [5,6]
              %[1,8,9];
              %[3,7,8]
                };
%what to call the plots
grouptitles = {'V1 Response L_4 vs. LW_4 ',...
               'V1 Response L_8 vs. LW_8)',...
               'V1 Response L_12 vs. LW_12)',...
               'V1 Response L_15 vs. LW_15)',...
               %'V1 Response LW vs. LW-Wdel vs. LW-Ldel'
               };

%and now plot
for grp = 1:numel(congrouped)
    figure()
sgtitle(grouptitles{grp})


for shnk = 1: size(shank_groups_types,2)
    subplot(3,size(shank_groups_types,2),1+((shnk-1)*size(shank_groups_types,2)))
    hold on
    for con = 1:numel(congrouped{grp})
        

        plot(time,smooth(squeeze(superCondis_avg(congrouped{grp}(con),1+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor),'Color',colors{con},'LineWidth',3)
        patch([time' flipud(time)'], ...
            [smooth(squeeze(superCondis_avg(congrouped{grp}(con),1+((shnk-1)*size(shank_groups_types,2)),wind))+squeeze(superCondis_std(congrouped{grp}(con),1+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor); ...
            fliplr(smooth(squeeze(superCondis_avg(congrouped{grp}(con),1+((shnk-1)*size(shank_groups_types,2)),wind))-squeeze(superCondis_std(congrouped{grp}(con),1+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor))],...
            colors{con}, 'FaceAlpha',0.2, 'EdgeColor','none')
    end
    ylim([0 300])
    hold off
    subplot(3,size(shank_groups_types,2),2+((shnk-1)*size(shank_groups_types,2)))

    hold on
    for con = 1:numel(congrouped{grp})
        plot(time,smooth(squeeze(superCondis_avg(congrouped{grp}(con),2+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor),'Color',colors{con},'LineWidth',3)
        patch([time' flipud(time)'], ...
            [smooth(squeeze(superCondis_avg(congrouped{grp}(con),2+((shnk-1)*size(shank_groups_types,2)),wind))+squeeze(superCondis_std(congrouped{grp}(con),2+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor); ...
            fliplr(smooth(squeeze(superCondis_avg(congrouped{grp}(con),2+((shnk-1)*size(shank_groups_types,2)),wind))-squeeze(superCondis_std(congrouped{grp}(con),2+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor))],...
            colors{con}, 'FaceAlpha',0.2, 'EdgeColor','none')
    end
    ylim([0 300])
    hold off
    subplot(3,size(shank_groups_types,2),3+((shnk-1)*size(shank_groups_types,2)))

    hold on
    for con = 1:numel(congrouped{grp})
      plot(time,smooth(squeeze(superCondis_avg(congrouped{grp}(con),3+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor),'Color',colors{con},'LineWidth',3)
        patch([time' flipud(time)'], ...
            [smooth(squeeze(superCondis_avg(congrouped{grp}(con),3+((shnk-1)*size(shank_groups_types,2)),wind))+squeeze(superCondis_std(congrouped{grp}(con),3+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor); ...
            fliplr(smooth(squeeze(superCondis_avg(congrouped{grp}(con),3+((shnk-1)*size(shank_groups_types,2)),wind))-squeeze(superCondis_std(congrouped{grp}(con),3+((shnk-1)*size(shank_groups_types,2)),wind)),smoothfactor))],...
            colors{con}, 'FaceAlpha',0.2, 'EdgeColor','none')
    end
    ylim([0 300])
    hold off
end
end
%%
figure(); 
hist(p1peaktimes{4},500)
figure();
hist(p1peaktimes{9},500)
%% plot spiking timeseries of channels with all the lines
wind = [4900:8000]';
yends = [0 500];
time = [-100:3000]';
smoothfactor = 100;

%colors for each group within a plot
colors = {'b' 'k' 'r'};

%which groups to plot together
congrouped = {[1,2,3];
              [1,4,5];
              [1,6,7];
              [1,8,9]
                };
%what to call the plots
grouptitles = {'V1 Response L vs. LW (4uW flash)',...
               'V1 Response L vs. LW (8uW flash)',...
               'V1 Response L vs. LW (12uW flash)',...
               'V1 Response L vs. LW (15uW flash)',...
               };

%and now plot
for grp = 1:numel(congrouped)

figure();


sgtitle(grouptitles{grp})
subplot(3,1,1)
title('L2/3')
   hold on
    
    for con = 1:numel(congrouped{grp})
        %plot(time,superCondis_rate{con}(:,wind,1)','Color', colors{con}, 'LineWidth',.1,'LineStyle','--')
        for tr = 1:size(superCondis_rate{congrouped{grp}(con)},1)
        patch([time' flipud(time)'], ...
            [smooth(superCondis_rate{congrouped{grp}(con)}(tr,wind,1),smoothfactor)' smooth(fliplr(superCondis_rate{con}(tr,wind,1)),smoothfactor)'],...
            colors{con}, 'EdgeColor',colors{con}, 'EdgeAlpha',.15,'FaceColor','none')
        end
    end
   ylim(yends)
   
subplot(3,1,2)
title('L4')
    hold on
    
    for con = 1:numel(congrouped{grp})
        for tr = 1:size(superCondis_rate{congrouped{grp}(con)},1)
        patch([time' flipud(time)'], ...
            [smooth(superCondis_rate{congrouped{grp}(con)}(tr,wind,2),smoothfactor)' smooth(fliplr(superCondis_rate{con}(tr,wind,2)),smoothfactor)'],...
            colors{con}, 'EdgeColor',colors{con}, 'EdgeAlpha',.15,'FaceColor','none')
        end
    end
    ylim(yends)
    ylabel('MUA (spikes/sec)')

subplot(3,1,3)
title('5')
    hold on
    
    for con = 1:numel(congrouped{grp})
        for tr = 1:size(superCondis_rate{[congrouped{grp}(con)]},1)
            patch([time' flipud(time)'], ...
            [smooth(superCondis_rate{congrouped{grp}(con)}(tr,wind,3),smoothfactor)' smooth(fliplr(superCondis_rate{con}(tr,wind,3)),smoothfactor)'],...
            colors{con}, 'EdgeColor',colors{con}, 'EdgeAlpha',.15,'FaceColor','none')
        end
    end
    ylim(yends)

subplot(3,1,1)
   
    for con = 1:numel(congrouped{grp})
    plot(time,smooth(squeeze(superCondis_avg(congrouped{grp}(con),1,wind)),smoothfactor/10),'Color',colors{con},'LineWidth',2)
    end
   
subplot(3,1,2)
    for con = 1:numel(congrouped{grp})
    plot(time,smooth(squeeze(superCondis_avg(congrouped{grp}(con),2,wind)),smoothfactor/10),'Color',colors{con},'LineWidth',2)
    end

subplot(3,1,3)
    for con = 1:numel(congrouped{grp})
    plot(time,smooth(squeeze(superCondis_avg(congrouped{grp}(con),3,wind)),smoothfactor/10),'Color',colors{con},'LineWidth',2)
    end
xlabel('Time Post Stimulus (ms)')
end

%%
% %% calculate summary data (peak and AUC)
% 
% 
% %calculate auc and peak for p1 and p2 of response
% 
% for con = 1:numel(condis)
%     p1auc_con = [];
%     p2auc_con = [];
%     allauc_con = [];
%     p1peak_con = [];
%     p2peak_con = [];
% 
%     for ch = 1:size(chans,1)
% 
%         %AUC
%         auc = cumtrapz(superCondis_rate{con}(:,5150:5350,ch),2);
%         p1auc_con = cat(1,p1auc_con,auc(:,end));
%         auc = cumtrapz(superCondis_rate{con}(:,5350:8000,ch),2);
%         p2auc_con = cat(1,p2auc_con,auc(:,end));
%         auc = cumtrapz(superCondis_rate{con}(:,5020:8000,ch),2);
%         allauc_con = cat(1,allauc_con,auc(:,end));
% 
%         %Peaks
%         peaks = max(superCondis_rate{con}(:,5150:5350,ch),[],2);
%         p1peak_con = cat(1,p1peak_con,peaks(:,end));
%         peaks = max(superCondis_rate{con}(:,5350:8000,ch),[],2);
%         p2peak_con = cat(1,p2peak_con,peaks(:,end));
%         %peaks = max(superCondis_rate{con}(:,5020:8000),[],2);
%         %allpeak{con} = peaks;
% 
%     end
% 
%     p1auc{con}  = p1auc_con;
%     p2auc{con}  = p2auc_con;
%     allauc{con} = allauc_con;
%     p1peak{con} = p1peak_con;
%     p2peak{con} = p2peak_con;
% end
%%

%they were placed in cell format so that a quick plot of various things
%could be made here if desired, without much fiddling with indexing. make
%that plot here, if you so choose. NOTE THAT ALL CHANNEL DATA IS IN SAME
%VECTOR, SO YOULL HAVE TO SPLIT THAT


%% make a table that can be easily read out in R for stats and plotting
trialdata_p1auc = [];
trialdata_p2auc = [];
trialdata_allauc = [];
trialdata_p1peak = [];
trialdata_p2peak = [];

for con = 1:numel(condis)
   trialdata_p1auc = vertcat(trialdata_p1auc,p1auc{con});
   trialdata_p2auc = vertcat(trialdata_p2auc,p2auc{con});
   trialdata_allauc = vertcat(trialdata_allauc,allauc{con});
   trialdata_p1peak = vertcat(trialdata_p1peak,p1peak{con});
   trialdata_p2peak = vertcat(trialdata_p2peak,p2peak{con});
end


tableres = table(trialdata_p1auc,trialdata_p2auc,trialdata_allauc,trialdata_p1peak,trialdata_p2peak,...
      triallabel_animalname, triallabel_animalnum, triallabel_condiname,triallabel_condinum,triallabel_condifile,...
    triallabel_region,triallabel_chanID,triallabel_changroup,triallabel_shankgroup,triallabel_channum,...
    'VariableNames', ["P1_AUC","P2_AUC","All_AUC","P1_Peak","P2_Peak",...
    "Animal_Name","Animal_Num","Condition_Name","Condition_Num","Condition_FileNum","Region","Channel_ID","Channel_GroupName","Shank_GroupName","Channel_GroupNum"]);

%% save the table as a csv for R

%give the file a name

dlgtitle = 'What file name to use for the data?';
promt = {'Give the file a name (no spaces please)'};
fieldsize = [1 150];
definput = {''};
opts.Resize = 'on';
opts.WindowStyle = 'normal';
filename = inputdlg(promt,dlgtitle,fieldsize,definput,opts);

%does it exist?
exists = isfile(fullfile(['E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\' cell2mat(filename)]));

while exists == 1
    dlgtitle = 'You stupid, dementia-riddled dumbass. You already made a file called that';
    promt = {'Give the file a name (that doesnt already exist this time)'};
    fieldsize = [1 150];
    definput = {''};
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    filename = inputdlg(promt,dlgtitle,fieldsize,definput,opts);

    exists = isfile(fullfile(['E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\' cell2mat(filename)]));

end

writetable(tableres,fullfile(['E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\' cell2mat(filename) '.csv']))
disp('Saved. Go to R, traitor.')






  