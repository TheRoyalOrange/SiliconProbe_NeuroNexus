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
    superCondis_dir{con} = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\LFP','Prompt', ['Choose ' condis{con} ' files']); 
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


% data = [];
% for i = 1:25
%     data = cat(2,data,squeeze(stim_lfp_stimchunks(i,:,30)));
% end
% fs = 1000;
% data_banded = bandpass(data,[1 155],fs, 'ImpulseResponse','iir');
% data_env = envelope(data_banded,50,'rms');
% 
% 
% data_env_trials = reshape(data_env',15000,25)';
% 
% cumtrapz(data_env_trials,2);
%%
%load the data you will be using (and only that data)

P1wind = [5075:5350];
P2wind = [5350:8000];

  P1_ampall = []; %total rms of signal in window
  P2_ampall = [];
  P1_peak = [];
  P1_peaktime = [];
  P2_peak = [];
  P2_peaktime = [];

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
triallabel_trialnum = [];

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
    triallabel_trialnum_con = [];


    
    for file = 1:numel(superCondis_dir{con})
        disp(['Running Condition ', num2str(con), ' (', condis{con}, '), File ', num2str(file)])
        dat = matfile(strjoin(superCondis_dir{con}(file)));
        %trs = 1:24;
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
                  dat.stim_lfp_stimchunks(trs,:,chans(shank_groups == shank_groups_types(shnk),condiFileMice{con}(file))));
              
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

   
   
    for ch = 1:size(superCondis_contrials,3)
        
        ch_env = reshape(envelope(bandpass(reshape(squeeze(superCondis_contrials(:,:,ch))',1,[]),[2 150],1000, 'ImpulseResponse','iir',Steepness=[.95 .95]),25,'rms')',size(superCondis_contrials,2),size(superCondis_contrials,1))';
            
        auc = cumtrapz(ch_env(:,P1wind),2);
        P1_ampall = cat(1,P1_ampall, auc(:,end));
        auc = cumtrapz(ch_env(:,P2wind),2);
        P2_ampall = cat(1,P2_ampall, auc(:,end));

        [peak,time] = max(-squeeze(superCondis_contrials(:,P1wind,ch))');
        P1_peak = cat(1,P1_peak, peak');
        P1_peaktime = cat(1,P1_peaktime,time'+(P1wind(1)-5000));
          
        [peak,time] = max(-squeeze(superCondis_contrials(:,P2wind,ch))');
        P2_peak = cat(1,P2_peak, peak');
        P2_peaktime = cat(1,P2_peaktime,time'+(P2wind(1)-5000));
        
        end
    
    
   
    end 
%%
colors = bone(12)
figure()
for i = 1:9

    %subplot(3,3,i)
    hold on
    %plot(tableres.P1_RMSsum(tableres.Condition_Num==i & tableres.Channel_GroupNum==2),'Color',colors(i,:))
    plot(smooth(tableres.P2_RMSsum(tableres.Condition_Num==i & tableres.Channel_GroupNum==2),5),'Color',colors(i,:))

end  
hold off

figure()
for i = 1:9
    %subplot(3,3,i)
    hold on
    %plot(tableres.P1_Peak(tableres.Condition_Num==i & tableres.Channel_GroupNum==2),'Color',colors(i,:))
    plot(smooth(tableres.P2_Peak(tableres.Condition_Num==i & tableres.Channel_GroupNum==2),5),'Color',colors(i,:))

end
hold off


figure()
for i = 1:9
    %subplot(3,3,i)
    hold on
    %plot(tableres.P1_PeakTime(tableres.Condition_Num==i & tableres.Channel_GroupNum==2),'Color',colors(i,:))
    plot(smooth(tableres.P2_PeakTime(tableres.Condition_Num==i & tableres.Channel_GroupNum==2),5),'Color',colors(i,:))

end
hold off
%% make table of results



tableres = table(P1_ampall,P2_ampall,P1_peak,P2_peak,P1_peaktime,P2_peaktime,...
      triallabel_animalname, triallabel_animalnum, triallabel_condiname,triallabel_condinum,triallabel_condifile,...
    triallabel_region,triallabel_chanID,triallabel_changroup,triallabel_shankgroup,triallabel_channum,...
    'VariableNames', ["P1_RMSsum","P2_RMSsum","P1_Peak","P2_Peak","P1_PeakTime","P2_PeakTime",...
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






  