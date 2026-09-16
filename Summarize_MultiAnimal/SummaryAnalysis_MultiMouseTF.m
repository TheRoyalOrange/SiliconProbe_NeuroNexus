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
    superCondis_dir{con} = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\TF','Prompt', ['Choose ' condis{con} ' files']); 
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

%highest frequency to check (range always starts at 1, max allowed is 151)
%freqstocheck = 80; %list the endpoint

%which time windows to get spectra of?
%P1wind = [3150:3350];
%P2wind = [3350:6000]

%what period to use for baselining power? (stim is at 3000)
%baseline_window = [800 2800];
%load the data you will be using (and only that data)
%superCondis_power = cell(numel(condis));


%triallabel_condiname = []; %metadata for later
%triallabel_condinum = []; 
%triallabel_condifile = []; 
%triallabel_animalnum = [];
%triallabel_animalname = [];
%triallabel_region = [];
%triallabel_chanID = [];
%triallabel_changroup = [];
%triallabel_channum = [];


%for con = 1:numel(condis)

    %superCondis_contrials_power = [];
  
%    triallabel_condiname_con = []; 
%    triallabel_condinum_con = []; 
%    triallabel_condifile_con = []; 
%    triallabel_animalnum_con = [];
%    triallabel_animalname_con = [];
%    triallabel_region_con = [];
%    triallabel_chanID_con = [];
%    triallabel_changroup_con = [];
%    triallabel_channum_con = [];

  
    
%    for file = 1:numel(superCondis_dir{con})
%        disp(['Running Condition ', num2str(con), ' (', condis{con}, '), File ', num2str(file)])
%        dat = matfile(strjoin(superCondis_dir{con}(file)));
%        trs = dat.tr_keep;
               
%        triallabel_condiname_con = cat(1,triallabel_condiname_con,repmat(repmat(condis{con},length(trs),1),length(chan_groups),1));
%        triallabel_condinum_con = cat(1,triallabel_condinum_con,repmat(repmat(con,length(trs),1),length(chan_groups),1));
%        triallabel_condifile_con = cat(1,triallabel_condifile_con,repmat(repmat(file,length(trs),1),length(chan_groups),1));
%        triallabel_animalnum_con = cat(1,triallabel_animalnum_con,repmat(repmat(condiFileMice{con}(file),length(trs),1),length(chan_groups),1));
%        triallabel_animalname_con = cat(1,triallabel_animalname_con,repmat(repmat(mice{condiFileMice{con}(file)},length(trs),1),length(chan_groups),1));
%        triallabel_region_con = cat(1,triallabel_region_con,repmat(repmat(region{condiFileMice{con}(file)},length(trs),1),length(chan_groups),1));
        
        
%        triallabel_chanID_con = cat(1,triallabel_chanID_con,reshape(repmat(chans(:,condiFileMice{con}(file)),1,length(trs))',1,[])');
%        triallabel_changroup_con = cat(1,triallabel_changroup_con,reshape(repmat(chan_groups',1,length(trs))',1,[])');
%        triallabel_channum_con = cat(1,triallabel_channum_con,reshape(repmat((1:size(chans,1))',1,length(trs))',1,[])');%

%        datpower = zeros(length(chans(:,condiFileMice{con}(file))),freqstocheck,3001,length(trs));
        
%        for ch = 1:length(chans(:,condiFileMice{con}(file)))
%            chpower = zeros(freqstocheck,6001,length(trs));
%            chpower_baselined = zeros(freqstocheck,3001,length(trs));

%            for frx = 1:freqstocheck
                %disp(['Condi' num2str(con) ' file' num2str(file) ' Ch' num2str(ch) ' freq' num2str(frx)])
%                chpower(frx,:,:) = abs(squeeze(dat.stim_tf(chans(ch,condiFileMice{con}(file)),frx,:,:))).^2; %compute frequency power
                %P1_chpower(frx,:) = squeeze(mean(abs(squeeze(dat.stim_tf(chans(ch,condiFileMice{con}(file)),frx,P1wind,:))).^2,1)); %compute frequency power
                %P2_chpower(frx,:) = squeeze(mean(abs(squeeze(dat.stim_tf(chans(ch,condiFileMice{con}(file)),frx,P2wind,:))).^2,1)); %compute frequency power
                
                %superCondis_contrials_itpc(ch,frx,:)      = abs( mean( exp(1i*angle(squeeze(data(chans(ch,condiFileMice{con}(file)),frx,:,:)))) ,2));
            
%            end
            
%            for tr = 1:size(chpower,3)
        
%                baseline = mean(squeeze(chpower(:,baseline_window(1):baseline_window(2),tr)),2);
                %baseline = mean(squeeze(dat.stim_tf(chans(ch,condiFileMice{con}(file)),:,baseline_window(1):baseline_window(2),tr)),2);
                %datpower_dB(:,:,tr) = 10*log10(bsxfun(@rdivide, squeeze(stim_tfpower(ch,:,:,tr)), squeeze(baseline));
%                chpower_baselined(:,:,tr) = bsxfun(@rdivide, squeeze(chpower(:,3001:end,tr)), squeeze(baseline));
    
%            end

%            datpower(ch,:,:,:) = chpower_baselined;
%        end


       %superCondis_contrials_power = cat(4,superCondis_contrials_power,datpower); %next need to do baseline and powerband calc
%     superCondis_power{con} = cat(4,superCondis_power{con},datpower);
%    end

%    triallabel_condiname = cat(1,triallabel_condiname,string(triallabel_condiname_con)); 
%    triallabel_condinum = cat(1,triallabel_condinum,triallabel_condinum_con);  
%    triallabel_condifile = cat(1,triallabel_condifile,triallabel_condifile_con); 
%    triallabel_animalnum = cat(1,triallabel_animalnum,triallabel_animalnum_con); 
%    triallabel_animalname = cat(1,triallabel_animalname,triallabel_animalname_con); 
%    triallabel_region = cat(1,triallabel_region,string(triallabel_region_con)); 
%    triallabel_chanID = cat(1,triallabel_chanID,triallabel_chanID_con); 
%    triallabel_changroup = cat(1,triallabel_changroup,triallabel_changroup_con); 
%    triallabel_channum = cat(1,triallabel_channum,triallabel_channum_con); 

    %superCondis_power{con} = superCondis_contrials_power;

%end
        



%%

%freqstocheck = 80; %list the endpoint

%which time windows to get spectra of?
P1wind = [3075:3350];
P2wind = [3350:6000];

%what period to use for baselining power? (stim is at 3000)
baseline_window = [900 2900];

%Frequency band definitions (add new ones if desired, but don't edit the standard ones) 
%UNCOMMENT THE ONES YOU WANT, AND DO THE SAME IN THE CORRESPONDING SECTION
%OF THE LOOP BELOW
    %delta = [2:4];
    %theta = [4:8];
    %alpha = [8:15];
    %beta = [15:30];
    %logamma = [30:80]; %I do not recommend doing the whole band
    %higamma = [80:151]; %I do not recommend doing the whole band
    alphabeta = [10:18]; %custom band
    lologamma = [30:50]; %custom band

    
%load the data you will be using (and only that data)
superCondis_P1pow_alphabeta = [];
superCondis_P1pow_lologamma = [];
superCondis_P2pow_alphabeta = [];
superCondis_P2pow_lologamma = [];


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

    %superCondis_contrials_power = [];
  
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
        %trs =  1:size(dat.stim_tf,4);
        % TODO: tr_keep is now always the full trial list (never shrunk by
        % permanent removal), so using it directly here no longer applies
        % tr_remove exclusion. Update to tr_keep_local = dat.tr_keep(~logical(dat.tr_remove)).
        trs = dat.tr_keep;
        %trs = 1:24;
        %datpower = zeros(length(chans(:,condiFileMice{con}(file))),length(delta),6001,length(trs));
        %deltapower = abs(dat.stim_tf(chans(:,condiFileMice{con}(file)),delta,:,trs)).^2;
        %thetapower = abs(dat.stim_tf(chans(:,condiFileMice{con}(file)),theta,:,trs)).^2;
        %alphapower = abs(dat.stim_tf(chans(:,condiFileMice{con}(file)),alpha,:,trs)).^2;
        %betapower = abs(dat.stim_tf(chans(:,condiFileMice{con}(file)),beta,:,trs)).^2;
        %logammapower = abs(dat.stim_tf(chans(:,condiFileMice{con}(file)),logamma,:,trs)).^2;
        %higammapower = abs(dat.stim_tf(chans(:,condiFileMice{con}(file)),higamma,:,trs)).^2;

        alphabetapower = zeros(length(chans(:,condiFileMice{con}(file))),length(alphabeta),6001,length(trs));
        %alphabetapower = [];
        for shnk = 1:length(shank_groups_types)
        alphabetapower_shan = abs(dat.stim_tf(chans(shank_groups == shank_groups_types(shnk),condiFileMice{con}(file)),alphabeta,:,trs)).^2;
        alphabetapower(shank_groups == shank_groups_types(shnk),:,:,:) = alphabetapower_shan;
        %trs = size(alphabetapower_shan,4);
        end
        clear alphabetapower_shan

        lologammapower = zeros(length(chans(:,condiFileMice{con}(file))),length(lologamma),6001,length(trs));
        %lologammapower = [];
        for shnk = 1:length(shank_groups_types)
        lologammapower_shan = abs(dat.stim_tf(chans(shank_groups == shank_groups_types(shnk),condiFileMice{con}(file)),lologamma,:,trs)).^2;
        lologammapower(shank_groups == shank_groups_types(shnk),:,:,:) = lologammapower_shan;
        %trs = size(lologammapower_shan,4);

        end
        clear lologammapower_shan

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

        P1pow_alphabeta_file = [];
        P2pow_alphabeta_file = [];
        P1pow_lologamma_file = [];
        P2pow_lologamma_file = [];

        %P1pow_alphabeta_trial = [];
        %P2pow_alphabeta_trial = [];
        %P1pow_lologamma_trial = [];
        %P2pow_lologamma_trial = [];
         
        for ch = 1:length(chans(:,condiFileMice{con}(file)))

            P1pow_alphabeta_trial = [];
            P2pow_alphabeta_trial = [];
            P1pow_lologamma_trial = [];
            P2pow_lologamma_trial = [];
            for tr = 1:length(trs)
            
             P1pow_alphabeta_trial(tr) = mean(alphabetapower(ch,:,P1wind,tr),'all')./mean(alphabetapower(ch,:,baseline_window(1):baseline_window(2),tr),'all');
             P2pow_alphabeta_trial(tr) = mean(alphabetapower(ch,:,P2wind,tr),'all')./mean(alphabetapower(ch,:,baseline_window(1):baseline_window(2),tr),'all');
                
             P1pow_lologamma_trial(tr) = mean(lologammapower(ch,:,P1wind,tr),'all')./mean(lologammapower(ch,:,baseline_window(1):baseline_window(2),tr),'all');
             P2pow_lologamma_trial(tr) = mean(lologammapower(ch,:,P2wind,tr),'all')./mean(lologammapower(ch,:,baseline_window(1):baseline_window(2),tr),'all');
         
            end

          %P1pow_alphabeta_file = cat(1,P1pow_alphabeta_file,P1pow_alphabeta_trial);
          %P1pow_lologamma_file = cat(1,P1pow_lologamma_file,P1pow_lologamma_trial);
          %P2pow_alphabeta_file = cat(1,P2pow_alphabeta_file,P2pow_alphabeta_trial);
          %P2pow_lologamma_file = cat(1,P2pow_lologamma_file,P2pow_lologamma_trial);

          superCondis_P1pow_alphabeta = cat(1,superCondis_P1pow_alphabeta,P1pow_alphabeta_trial'); 
          superCondis_P1pow_lologamma = cat(1,superCondis_P1pow_lologamma,P1pow_lologamma_trial');
          superCondis_P2pow_alphabeta = cat(1,superCondis_P2pow_alphabeta,P2pow_alphabeta_trial');
          superCondis_P2pow_lologamma = cat(1,superCondis_P2pow_lologamma,P2pow_lologamma_trial');

        end
 
    
    end

    triallabel_condiname = cat(1,triallabel_condiname,string(triallabel_condiname_con)); 
    triallabel_condinum = cat(1,triallabel_condinum,triallabel_condinum_con);  
    triallabel_condifile = cat(1,triallabel_condifile,triallabel_condifile_con); 
    triallabel_animalnum = cat(1,triallabel_animalnum,triallabel_animalnum_con); 
    triallabel_animalname = cat(1,triallabel_animalname,triallabel_animalname_con); 
    triallabel_region = cat(1,triallabel_region,string(triallabel_region_con)); 
    triallabel_chanID = cat(1,triallabel_chanID,triallabel_chanID_con); 
    triallabel_changroup = cat(1,triallabel_changroup,triallabel_changroup_con); 
    triallabel_shankgroup = cat(1,triallabel_shankgroup_con,triallabel_shankgroup_con); 
    triallabel_channum = cat(1,triallabel_channum,triallabel_channum_con); 

    %superCondis_power{con} = superCondis_contrials_power;

end
        
%%



%% make a table that can be easily read out in R for stats and plotting



tableres = table(superCondis_P1pow_alphabeta,superCondis_P2pow_alphabeta,superCondis_P1pow_lologamma,superCondis_P2pow_lologamma,...
      triallabel_animalname, triallabel_animalnum, triallabel_condiname,triallabel_condinum,triallabel_condifile,...
    triallabel_region,triallabel_chanID,triallabel_changroup,triallabel_channum,...
    'VariableNames', ["P1Power_AlphaBeta","P2Power_AlphaBeta","P1Power_LoloGamma","P2Power_LoloGamma",...
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

