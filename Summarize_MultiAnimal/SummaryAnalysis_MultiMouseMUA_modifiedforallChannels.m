%list mice to be analyzed. Note the order as it will be treated as a factor (R style)
%if mouse has more than one group of channels to analyze, list
%it twice here
mice = {"20260226-p12"};%, "20260402-p12"}; 

%list channels to use for each mouse listed (ex. groups of three here for cortical
%layers. One group per animal, in same order as 'mice'). 
%Each channel is compared with its corresponding channel in other mice, so
%each entry should have the same number of channels per mouse
chans = [1:64; %mouse 1
         %12,13,14 %mouse 2, etc
         ]'; 

%give a name to each channel per mouse (should match size of dim1 of chans)
chan_groups = [string(1:64)]; %["L2/3", "L4", "L5"];

%check num of chan_groups and chans per mouse match
size(chan_groups,2) == size(chans,1)

%which brain area(s) are being recorded by the channels? (order should
%match dim2 of chans and length of mice
region = {"V1"};

condis = {'L_4','LW_4', 'L_8','LW_8','L_12','LW_12','L_15','LW_15'}; %list conditions to be included, named as you'd prefer. Note the order
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
condiFileMice = {[1,1]; %condition 1 files
                 [1,1]; %condition 2
                 [1,1]; %condition 3
                 [1,1];   %condition 4
                 [1];   %condition 5
                 [1];   %condition 6
                 [1];   %condition 7
                 [1];   %condition 8
                 %[1,2]   %condition 9
                                 
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
    triallabel_channum_con = [];

    p1auc_con = [];
    p2auc_con = [];
    allauc_con = [];
    p1peak_con = [];
    p2peak_con = [];
    
    for file = 1:numel(superCondis_dir{con})
        disp(['Running Condition ', num2str(con), ' (', condis{con}, '), File ', num2str(file)])
        dat = matfile(strjoin(superCondis_dir{con}(file)));
        % TODO: tr_keep is now always the full trial list (never shrunk by
        % permanent removal), so using it directly here no longer applies
        % tr_remove exclusion. Update to tr_keep_local = dat.tr_keep(~logical(dat.tr_remove)).
        trs = dat.tr_keep;
               
        triallabel_condiname_con = cat(1,triallabel_condiname_con,repmat(repmat(condis{con},length(trs),1),length(chan_groups),1));
        triallabel_condinum_con = cat(1,triallabel_condinum_con,repmat(repmat(con,length(trs),1),length(chan_groups),1));
        triallabel_condifile_con = cat(1,triallabel_condifile_con,repmat(repmat(file,length(trs),1),length(chan_groups),1));
        triallabel_animalnum_con = cat(1,triallabel_animalnum_con,repmat(repmat(condiFileMice{con}(file),length(trs),1),length(chan_groups),1));
        triallabel_animalname_con = cat(1,triallabel_animalname_con,repmat(repmat(mice{condiFileMice{con}(file)},length(trs),1),length(chan_groups),1));
        triallabel_region_con = cat(1,triallabel_region_con,repmat(repmat(region{condiFileMice{con}(file)},length(trs),1),length(chan_groups),1));
        
        
        triallabel_chanID_con = cat(1,triallabel_chanID_con,reshape(repmat(chans(:,condiFileMice{con}(file)),1,length(trs))',1,[])');
        triallabel_changroup_con = cat(1,triallabel_changroup_con,reshape(repmat(chan_groups',1,length(trs))',1,[])');
        triallabel_channum_con = cat(1,triallabel_channum_con,reshape(repmat((1:size(chans,1))',1,length(trs))',1,[])');

        superCondis_contrials = cat(1,superCondis_contrials, dat.stim_spike_stimchunks(trs,:,chans(:,condiFileMice{con}(file))));
    end

    triallabel_condiname = cat(1,triallabel_condiname,string(triallabel_condiname_con)); 
    triallabel_condinum = cat(1,triallabel_condinum,triallabel_condinum_con);  
    triallabel_condifile = cat(1,triallabel_condifile,triallabel_condifile_con); 
    triallabel_animalnum = cat(1,triallabel_animalnum,triallabel_animalnum_con); 
    triallabel_animalname = cat(1,triallabel_animalname,triallabel_animalname_con); 
    triallabel_region = cat(1,triallabel_region,string(triallabel_region_con)); 
    triallabel_chanID = cat(1,triallabel_chanID,triallabel_chanID_con); 
    triallabel_changroup = cat(1,triallabel_changroup,triallabel_changroup_con); 
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

        %spikes_ms(:,:,ch) = spikeper;
        spikes_rate(:,:,ch) = squeeze(movmean(spikeper,50,2)*1000);
        spikes_rateavg(ch,:) = mean(squeeze(spikes_rate(:,:,ch)),1);
        spikes_ratestd(ch,:) = std(squeeze(spikes_rate(:,:,ch)),1);
    

        auc = cumtrapz(spikes_rate(:,5150:5350,ch),2);
        p1auc_con = cat(1,p1auc_con,auc(:,end));
        auc = cumtrapz(spikes_rate(:,5350:8000,ch),2);
        p2auc_con = cat(1,p2auc_con,auc(:,end));
        auc = cumtrapz(spikes_rate(:,5020:8000,ch),2);
        allauc_con = cat(1,allauc_con,auc(:,end));

        %Peaks
        peaks = max(spikes_rate(:,5150:5350,ch),[],2);
        p1peak_con = cat(1,p1peak_con,peaks(:,end));
        peaks = max(spikes_rate(:,5350:8000,ch),[],2);
        p2peak_con = cat(1,p2peak_con,peaks(:,end));
        %peaks = max(superCondis_rate{con}(:,5020:8000),[],2);
        %allpeak{con} = peaks;

   
    end 
    %superCondis{con} = spikes_ms;
    %superCondis_rate{con} = spikes_rate;
    superCondis_avg(con,:,:) = spikes_rateavg;
    superCondis_std(con,:,:) = spikes_ratestd;
    
    p1auc{con}  = p1auc_con;
    p2auc{con}  = p2auc_con;
    allauc{con} = allauc_con;
    p1peak{con} = p1peak_con;
    p2peak{con} = p2peak_con;
end







%% plot spiking timeseries of channels
wind = [4900:8000]';
time = [-100:3000]';
smoothfactor = 50;
colors = {'k' 'r' 'b'};
figure();
sgtitle('Heres a fucking title for ya')
subplot(3,1,1)
    hold on
    for con = 1:numel(condis)
        

        plot(time,smooth(squeeze(superCondis_avg(con,1,wind)),smoothfactor),'Color',colors{con},'LineWidth',2)
        patch([time' flipud(time)'], ...
            [smooth(squeeze(superCondis_avg(con,1,wind))+squeeze(superCondis_std(con,1,wind)),smoothfactor); ...
            fliplr(smooth(squeeze(superCondis_avg(con,1,wind))-squeeze(superCondis_std(con,1,wind)),smoothfactor))],...
            colors{con}, 'FaceAlpha',0.2, 'EdgeColor','none')
    end
    hold off
subplot(3,1,2)
    hold on
    for con = 1:numel(condis)
        plot(time,smooth(squeeze(superCondis_avg(con,2,wind)),smoothfactor),'Color',colors{con})
        patch([time' flipud(time)'], ...
            [smooth(squeeze(superCondis_avg(con,2,wind))+squeeze(superCondis_std(con,2,wind)),smoothfactor); ...
            fliplr(smooth(squeeze(superCondis_avg(con,2,wind))-squeeze(superCondis_std(con,2,wind)),smoothfactor))],...
            colors{con}, 'FaceAlpha',0.2, 'EdgeColor','none')
    end
    hold off
subplot(3,1,3)
    hold on
    for con = 1:numel(condis)
      plot(time,smooth(squeeze(superCondis_avg(con,3,wind)),smoothfactor),'Color',colors{con})
        patch([time' flipud(time)'], ...
            [smooth(squeeze(superCondis_avg(con,3,wind))+squeeze(superCondis_std(con,3,wind)),smoothfactor); ...
            fliplr(smooth(squeeze(superCondis_avg(con,3,wind))-squeeze(superCondis_std(con,3,wind)),smoothfactor))],...
            colors{con}, 'FaceAlpha',0.2, 'EdgeColor','none')
    end
    hold off
 
%%

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


%% calculate summary data (peak and AUC)


%calculate auc and peak for p1 and p2 of response

%for con = 1:numel(condis)
%    p1auc_con = [];
%    p2auc_con = [];
%    allauc_con = [];
%    p1peak_con = [];
%    p2peak_con = [];
    
%    for ch = 1:size(chans,1)

        %AUC
%        auc = cumtrapz(superCondis_rate{con}(:,5150:5350,ch),2);
%        p1auc_con = cat(1,p1auc_con,auc(:,end));
%        auc = cumtrapz(superCondis_rate{con}(:,5350:8000,ch),2);
%        p2auc_con = cat(1,p2auc_con,auc(:,end));
%        auc = cumtrapz(superCondis_rate{con}(:,5020:8000,ch),2);
%        allauc_con = cat(1,allauc_con,auc(:,end));

        %Peaks
%        peaks = max(superCondis_rate{con}(:,5150:5350,ch),[],2);
%        p1peak_con = cat(1,p1peak_con,peaks(:,end));
%        peaks = max(superCondis_rate{con}(:,5350:8000,ch),[],2);
%        p2peak_con = cat(1,p2peak_con,peaks(:,end));
        %peaks = max(superCondis_rate{con}(:,5020:8000),[],2);
        %allpeak{con} = peaks;

%    end
    
%    p1auc{con}  = p1auc_con;
%    p2auc{con}  = p2auc_con;
%    allauc{con} = allauc_con;
%    p1peak{con} = p1peak_con;
%    p2peak{con} = p2peak_con;
%end
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
    triallabel_region,triallabel_chanID,triallabel_changroup,triallabel_channum,...
    'VariableNames', ["P1_AUC","P2_AUC","All_AUC","P1_Peak","P2_Peak",...
    "Animal_Name","Animal_Num","Condition_Name","Condition_Num","Condition_FileNum","Region","Channel_ID","Channel_GroupName","Channel_GroupNum"]);

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


%save(fullfile(['E:\Roy\Processed Silicon Probe Data\BundledAnimalData\20260226-p12_wholeprobe_means']), 'superCondis_avg', 'superCondis_std','superCondis_dir')



%%
wind = [4900:8000]';
yends = [0 300];
time = [-100:3000]';
smoothfactor = 100;
%colors for each group within a plot
colors = {'k' 'r'};
congrouped = {[1,2];
              [3,4];
              [5,6];
              [7,8]
                };

for i = 1:length(ProbeInfo.poi)
    poilayout{i} = ProbeInfo.ProbeMaps{ProbeInfo.poi(i)};
    poilayoutT{i} = reshape(poilayout{i}',1,[]);
    poimapT{i} = cell2mat(ProbeInfo.ProbeMaps(ProbeInfo.poi(i)))';
end



for congroup = 1:size(congrouped,1)
   
 
 figure()
    %sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' Spike Avg'])
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
            loc = poilayoutT{prb}(chan);
            if loc > 0
            subplot(size(poimapT{prb},1),size(poimapT{prb},2),chan)
            hold on
            for con = 1:length(congrouped{congroup})
                plot(time,smooth(squeeze(superCondis_avg(congrouped{congroup}(con),poimapT{prb}(chan),wind))',smoothfactor),"Color",colors{con});
                patch([time' flipud(time)'], ...
                [smooth(squeeze(superCondis_avg(congrouped{congroup}(con),poimapT{prb}(chan),wind))+squeeze(superCondis_std(congrouped{congroup}(con),poimapT{prb}(chan),wind)),smoothfactor); ...
                fliplr(smooth(squeeze(superCondis_avg(congrouped{congroup}(con),poimapT{prb}(chan),wind))-squeeze(superCondis_std(congrouped{congroup}(con),poimapT{prb}(chan),wind)),smoothfactor))],...
                colors{con}, 'FaceAlpha',0.2, 'EdgeColor','none')
                
                xlim([time(1) time(end)]);
                title(['CH',num2str(poimapT{prb}(chan))])
                ylim(yends);
                %xline(0,'Color','c')
                title(['CH',num2str(loc)])
            end
            hold off
        end
        
   
    end
end

%% same but with box and whisker of conditions



congrouped = {[1,2];
              [3,4];
              [5,6];
              [7,8]
                };

for i = 1:length(ProbeInfo.poi)
    poilayout{i} = ProbeInfo.ProbeMaps{ProbeInfo.poi(i)};
    poilayoutT{i} = reshape(poilayout{i}',1,[]);
    poimapT{i} = cell2mat(ProbeInfo.ProbeMaps(ProbeInfo.poi(i)))';
end



for congroup = 1:size(congrouped,1)
   
 
 figure()
    %sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' Spike Avg'])
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
            loc = poilayoutT{prb}(chan);
            if loc > 0
            subplot(size(poimapT{prb},1),size(poimapT{prb},2),chan)
            hold on
            
            for con = 1:length(congrouped{congroup})
                dattoplot = tableres(tableres.Condition_Num==congrouped{congroup}(1)|...
                                        tableres.Condition_Num==congrouped{congroup}(2),:);
                boxplot(dattoplot.P2_AUC,dattoplot.Condition_Name);
                               
                title(['CH',num2str(poimapT{prb}(chan))])
               %title(['CH',num2str(loc)])
            end
            hold off
        end
        
   
    end
end

