%list mice to be analyzed. Note the order as it will be treated as a factor (R style)
%if mouse has more than one group of channels to analyze, list
%it twice here
mice = {"20260618-p9"}; 

%list channels to use for each mouse listed (ex. groups of three here for cortical
%layers. One group per animal, in same order as 'mice'). 
%Each channel is compared with its corresponding channel in other mice, so
%each entry should have the same number of channels per mouse
chans = [
         36 37 38 %29 30 31 37 38 39%mouse 2, etc
         ]'; 
%[12 13 14 20 21 22 28 29 30]';

%give a name to each channel per mouse (should match size of dim1 of chans)
chan_groups = ["L2/3", "L4", "L5"];%,"L2/3", "L4", "L5","L2/3", "L4", "L5"];
shank_groups = ["C","C","C"];%["L","L","L", "C","C","C","R","R","R"];
shank_groups_types = unique(shank_groups,'stable');
%check num of chan_groups and chans per mouse match
size(chan_groups,2) == size(chans,1) && size(chans,1) == size(shank_groups,2)

%which brain area(s) are being recorded by the channels? (order should
%match dim2 of chans and length of mice
region = {"V1"};

condis = {'W','W_TTX'}%,'L48','L56'}; %list conditions to be included, named as you'd prefer. Note the order
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
condiFileMice = {[1 1];
                   [1 1]
                   %[1,1];
                   %[1,1]
                   };%condition 1 files
                  
                 
                                 
                


%check num of files and filemice match
for con = 1:length(condis)
    size(superCondis_dir{con},2) == length(condiFileMice{con})
end
%%
phases = [];
phases_oi = [];
phases_noi = [];
itpc_all = {};
phases_all = {};
%tr2check = {[1:25]
            %[1:25]
            %[1:25]
            %[1:25]
            %[1:25]
            %};
for con = 1:length(condis)     
    
    trs = [];
    phases_con = [];

    for file = 1:numel(superCondis_dir{con})
        dat = matfile(strjoin(superCondis_dir{con}(file)));
        %trs =  tr2check{file};
         disp(num2str(file))       
        for shnk = 1:length(shank_groups_types)
        phase_shank = angle(dat.stim_tf(chans(shank_groups == shank_groups_types(shnk),condiFileMice{con}(file)),:,1000:6000,:));
        phases(shank_groups == shank_groups_types(shnk),:,:,:) = phase_shank;
        end
        %phases_oi = cat(4,phases_oi,phases(:,:,:,trs));
        %phases(:,:,:,trs) = [];
        phases_con = cat(4,phases_con,phases(:,:,:,:));
        
        phases = [];
        
        %phases_noi = cat(4,phases_noi,phases);
        %phases = [];
    end
    phases_all{con} = phases_con;

    for ch = 1:size(phases_con,1)
        itpc_all{con}(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_con(ch,:,:,:)))),3));
    end
end

%%
%ITPC
        %itpc = [];
        %itpc_noi = [];

        %for ch = 1:size(phases_oi,1)
        %    itpc(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_oi(ch,:,:,1:79)))),3));
        %    itpc_noi(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_oi(ch,:,:,80:end)))),3));
        %end

     
        figure();
            subplot(3,3,1)
            contour(squeeze(itpc(1,:,900:end)),150)
            subplot(3,3,4)
            contour(squeeze(itpc(2,:,900:end)),150)
            subplot(3,3,7)
            contour(squeeze(itpc(3,:,900:end)),150)
            subplot(3,3,2)
            contour(squeeze(itpc(4,:,900:end)),150)
            subplot(3,3,5)
            contour(squeeze(itpc(5,:,900:end)),150)
            subplot(3,3,8)
            contour(squeeze(itpc(6,:,900:end)),150)
            subplot(3,3,3)
            contour(squeeze(itpc(7,:,900:end)),150)
            subplot(3,3,6)
            contour(squeeze(itpc(8,:,900:end)),150)
            subplot(3,3,9)
            contour(squeeze(itpc(9,:,900:end)),150)

        figure()
        contour([-1000:1000],[1:151],squeeze(mean(itpc(3,:,:),1)),150,'Fill', 'on')
        set(gca,'CLim',[0 1])

        figure();
            subplot(3,3,1)
            contour(squeeze(itpc_noi(1,:,1000:end)),150)
            subplot(3,3,4)
            contour(squeeze(itpc_noi(2,:,1000:end)),150)
            subplot(3,3,7)
            contour(squeeze(itpc_noi(3,:,1000:end)),150)
            subplot(3,3,2)
            contour(squeeze(itpc_noi(4,:,1000:end)),150)
            subplot(3,3,5)
            contour(squeeze(itpc_noi(5,:,1000:end)),150)
            subplot(3,3,8)
            contour(squeeze(itpc_noi(6,:,1000:end)),150)
            subplot(3,3,3)
            contour(squeeze(itpc_noi(7,:,1000:end)),150)
            subplot(3,3,6)
            contour(squeeze(itpc_noi(8,:,1000:end)),150)
            subplot(3,3,9)
            contour(squeeze(itpc_noi(9,:,1000:end)),150)

        figure()
        contour(squeeze(mean(itpc_noi(3,:,:),1)),150,'Fill', 'on')
        set(gca,'CLim',[0 1])
%% 
Xlim = [950 1250];
stimline = [1000];
clim = [0 .5];
Ylim = [1 80];

figure()
sgtitle('ITPC of Active V1 Post-Whisker Stim')
subplot(3,1,1)
hold on
        title('L2/3')
        contour(squeeze(itpc(1,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)
        xline(stimline)
       hold off
 subplot(3,1,2)
 hold on
        title('L4')
        contour(squeeze(itpc(2,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)
        xline(stimline)
        hold off
 subplot(3,1,3)
 hold on
 title('L5')
        contour(squeeze(itpc(3,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)
        xline(stimline)
        
hold off


figure()
sgtitle('ITPC of Active V1 Post-Whisker Stim w/ TTX-S1')
subplot(3,1,1)
hold on
        title('L2/3')
        contour(squeeze(itpc_noi(1,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
 subplot(3,1,2)
 hold on
        title('L4')
        contour(squeeze(itpc_noi(2,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)
        xline(stimline)
 hold off
 subplot(3,1,3)
 hold on
 title('L5')
        contour(squeeze(itpc_noi(3,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)
        xline(stimline)

hold off


figure()
sgtitle('ITPC Difference of Active V1 Post-Whisker Stim (Norm-TTX)')
subplot(3,1,1)
hold on
        title('L2/3')
        contour(squeeze(itpc(1,:,:))-squeeze(itpc_noi(1,:,:)),150,'Fill', 'on')
        set(gca,'XLim',Xlim,'YLim',Ylim)
        colorbar()
        xline(stimline)
        hold off
 subplot(3,1,2)
 hold on
        title('L4')
        contour(squeeze(itpc(2,:,:))-squeeze(itpc_noi(2,:,:)),150,'Fill', 'on')
        set(gca,'XLim',Xlim,'YLim',Ylim)
        colorbar()  
        xline(stimline)

        hold off
 subplot(3,1,3)
 hold on
 title('L5')
        contour(squeeze(itpc(3,:,:))-squeeze(itpc_noi(3,:,:)),150,'Fill', 'on')
        set(gca,'XLim',Xlim,'YLim',Ylim)
        xline(stimline)

        colorbar()
hold off

%% plot itpc of TTX trials as TTX wears off
g1 = [1:15]+72;
g2 = [16:30]+72; %choosing trial groupings
g3 = [31:45]+72;
g4 = [45:60]+72;
g5 = [60:75]+72;
        
        for ch = 1:size(phases_oi,1)
            itpc_g1(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_oi(ch,:,:,g1)))),3));
            itpc_g2(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_oi(ch,:,:,g2)))),3));
            itpc_g3(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_oi(ch,:,:,g3)))),3));
            itpc_g4(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_oi(ch,:,:,g4)))),3));
            itpc_g5(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_oi(ch,:,:,g5)))),3));

            %itpc_noi(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_oi(ch,:,:,72:148)))),3));
        end

%figure for l2/3
figure()
sgtitle('ITPC of Active V1 L2/3 Post-Whisker Stim Across TTX Recovery')
subplot(5,1,1)
hold on
        title('Trials 1:15')
        contour(squeeze(itpc_g1(1,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
 subplot(5,1,2)
hold on
        title('Trials 16:30')
        contour(squeeze(itpc_g2(1,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off

subplot(5,1,3)
hold on
        title('Trials 31:45')
        contour(squeeze(itpc_g3(1,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
        subplot(5,1,4)
hold on
        title('Trials 46:60')
        contour(squeeze(itpc_g4(1,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
        subplot(5,1,5)
hold on
        title('Trials 61:75')
        contour(squeeze(itpc_g5(1,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off

%figure for L4
figure()
sgtitle('ITPC of Active V1 L4 Post-Whisker Stim Across TTX Recovery')
subplot(5,1,1)
hold on
        title('Trials 1:15')
        contour(squeeze(itpc_g1(2,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
 subplot(5,1,2)
hold on
        title('Trials 16:30')
        contour(squeeze(itpc_g2(2,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off

subplot(5,1,3)
hold on
        title('Trials 31:45')
        contour(squeeze(itpc_g3(2,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
        subplot(5,1,4)
hold on
        title('Trials 46:60')
        contour(squeeze(itpc_g4(2,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
        subplot(5,1,5)
hold on
        title('Trials 61:75')
        contour(squeeze(itpc_g5(2,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off

        %figure for l2/3
figure()
sgtitle('ITPC of Active V1 L5 Post-Whisker Stim Across TTX Recovery')
subplot(5,1,1)
hold on
        title('Trials 1:15')
        contour(squeeze(itpc_g1(3,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
 subplot(5,1,2)
hold on
        title('Trials 16:30')
        contour(squeeze(itpc_g2(3,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off

subplot(5,1,3)
hold on
        title('Trials 31:45')
        contour(squeeze(itpc_g3(3,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
        subplot(5,1,4)
hold on
        title('Trials 46:60')
        contour(squeeze(itpc_g4(3,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off
        subplot(5,1,5)
hold on
        title('Trials 61:75')
        contour(squeeze(itpc_g5(3,:,:)),150,'Fill', 'on')
        set(gca,'CLim',clim,'XLim',Xlim,'YLim',Ylim)  
        xline(stimline)
        hold off


%% load lfp results and do stim-triggerd lfp and stim-triggered phase 

for con = 1:length(condis)
    superCondis_dir{con} = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\LFP','Prompt', ['Choose ' condis{con} ' files']); 
end
%load lfp
lfp = [];
lfp_oi = [];
for con = 1:length(condis)
    for file = 1:numel(superCondis_dir{con})
        dat = matfile(strjoin(superCondis_dir{con}(file)));
        %trs =  tr2check{file};
         disp(num2str(file))       
        for shnk = 1:length(shank_groups_types)
        lfp_shank = dat.stim_lfp_stimchunks(:,4000:8000,chans(shank_groups == shank_groups_types(shnk),condiFileMice{con}(file)));
        lfp(:,:,shank_groups == shank_groups_types(shnk)) = lfp_shank;
        end
        %phases_oi = cat(4,phases_oi,phases(:,:,:,trs));
        %phases(:,:,:,trs) = [];
        lfp_oi = cat(1,lfp_oi,lfp(:,:,:));
        lfp = [];
        
        %phases_noi = cat(4,phases_noi,phases);
        %phases = [];
end
end

figure()
subplot(2,3,1)
plot(squeeze(lfp_oi(1:72,:,1))','Color','k','LineWidth',.1)
xlim([1000 1150])

subplot(2,3,2)
plot(squeeze(lfp_oi(1:72,:,2))','Color','k','LineWidth',.1)
xlim([1000 1150])

subplot(2,3,3)
plot(squeeze(lfp_oi(1:72,:,3))','Color','k','LineWidth',.1)
xlim([1000 1150])

subplot(2,3,4)
plot(squeeze(lfp_oi(73:end,:,1))','Color','k','LineWidth',.1)
xlim([1000 1150])

subplot(2,3,5)
plot(squeeze(lfp_oi(73:end,:,2))','Color','k','LineWidth',.1)
xlim([1000 1150])

subplot(2,3,6)
plot(squeeze(lfp_oi(73:end,:,3))','Color','k','LineWidth',.1)
xlim([1000 1150])


%%
%stim triggered phase
figure()
hold on
for freq = 20:50
    plot(squeeze(phases_oi(1,freq,950:1200,1:72)) + 6*(freq-1),'Color','k','LineWidth',.01)
end
hold off
%%
%ISPC
   ispc = [];
        ispc_noi = [];
        for tr = 1:72
            ispc(:,:,tr) = abs(mean(exp(1i*(squeeze(phases_oi(:,:,:,tr)))),1));
        end
        for tr = 1:length(72:148)
            ispc_noi(:,:,tr) = abs(mean(exp(1i*(squeeze(phases_oi(:,:,:,i+72)))),1));
        end

figure()
for i = 1:75
    subplot(5,5,i)
    hold on
    contour(squeeze(ispc(:,:,i)),150)
    %set(gca,'CLim',[0 1])
    %colorbar("Limits",[0 1])
    xline(1000)
    hold off
end


map = [0.6,0.1,0.2
    0.5,0.4,0.6
    0.3,0.7,0.9]; 

[val, ind]=max(cat(1,itpc(1,1:80,1000:1500),itpc(2,1:80,1000:1500),itpc(3,1:80,1000:1500)));
figure()
contour(squeeze(ind),'Fill','on')
colormap(map)

[val, ind]=max(cat(1,itpc_noi(1,1:80,1000:1500),itpc_noi(2,1:80,1000:1500),itpc_noi(3,1:80,1000:1500)));
figure()
contour(squeeze(ind),'Fill','on')
colormap(map)

map = [0.6,0.1,0.2
    0.5,0.4,0.6
    0.3,0.7,0.9]; 

figure();
subplot(3,1,1)
hold on 
title('L2/3')
[val, ind]=max(cat(1,itpc(1,:,:),itpc_noi(1,:,:)));
contour(squeeze(ind),'Fill','on')
colormap(map)

subplot(3,1,2)
hold on 
title('L4')
[val, ind]=max(cat(1,itpc(2,:,:),itpc_noi(2,:,:)));
contour(squeeze(ind),'Fill','on')
colormap(map)

subplot(3,1,3)
hold on 
title('L5')
[val, ind]=max(cat(1,itpc(3,:,:),itpc_noi(3,:,:)));
contour(squeeze(ind),'Fill','on')
colormap(map)

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

