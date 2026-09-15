
%% THIS SCRIPT IS NOT READY TO RUN AUTONOMOUSLY. FILE SAVING IS HARD CODED AND COULD LEAD TO MISLABELLED FILES OR OVERWRITING
%% DO NOT HIT 'RUN'! THIS SCRIPT IS NOT WRITTEN FOR SUCH THINGS YET. 
%list mouse to be analyzed
animal = '20260820-p8'; 



%which brain area to plot?
region = 'V1';

%which conditions to analyze? (in order)
condis = {'W', 'W_TTX'}; %list conditions to be included, named as you'd prefer. Note the order

%{'W','L32','LW32','L40','LW40','L48','LW48','L56','LW56'}; %; %list conditions to be included, named as you'd prefer. Note the order
% as they will be treated as factors later

%group the conditions as you please
groupassign ={[1],[2]};% {[1, 3, 5, 7, 9],[2, 4, 6, 8]};

%name to groups
groupnames = {'Whisker', 'Whisker_TTX'};

%load probeinfo for animal(IT IS A GOOD IDEA TO HAVE DONE ANY LABELLING OF
%BAD CHANNELS ALREADY, TO SAVE TIME IN THIS ANALYSIS)
load(fullfile(['E:\Roy\Processed Silicon Probe Data\ProbeInfo\' animal '-ProbeInfo.mat']));
chansall = ProbeInfo.ProbeIds{find(ismember(ProbeInfo.Areas, region))};
chansuse = chansall;
chansuse(ProbeInfo.Ch_Remove{1}) = [];
%% Get directory for files to be analyzed
for con = 1:length(condis)
    superCondis_dir{con} = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\TF','Prompt', ['Choose ' condis{con} ' files']); 
end
%%


% %now, manually fill in an array with numbers 1 to numberofmice matching
% %   order of files in superCondis_dir{con}, and repeat for all conditions
% % Ex. if the first file in a condition corresponds to the first listed mouse in 'mice',
% %   then the first entry in condiFileMice is 1
% condiFileMice = {[1];
%                   [1];
%                   [1];
%                   [1];
%                   [1];
%                   [1];
%                   [1];
%                   [1];
%                    [1]
%                    %[1,1];
%                    %[1,1]
%                    };%condition 1 files
% 
% 
% 
% 
% 
% 
% %check num of files and filemice match
% for con = 1:length(condis)
%     size(superCondis_dir{con},2) == length(condiFileMice{con})
% end

%%

%how many trials from each condition to include?
%trinc =  1:25;

%get phases from each trial for good channels
phases = [];
phases_condi = cell(1,size(condis,2));



%tr2check = {[1:25]
            %[1:25]
            %[1:25]
            %[1:25]
            %[1:25]
            %};
for con = 1:length(condis)
    trs{con}= [];
    for file = 1:numel(superCondis_dir{con})
        dat = matfile(strjoin(superCondis_dir{con}(file)));
        
        disp(num2str(file))       
        %for shnk = 1:length(shank_groups_types)
        phases = angle(dat.stim_tf(chansall,:,2000:4000,:));
        
        trs{con} = cat(2,trs{con},1:size(phases,4));
        
        %phases_oi = cat(4,phases_oi,phases(:,:,:,trs));
        %phases(:,:,:,trs) = [];
        phases_condi{con} = cat(4,phases_condi{con},phases(chansuse,:,:,:));
        phases = [];
        
        %phases_noi = cat(4,phases_noi,phases);
        %phases = [];
   end
end

%% ITPC
itpc = cell(1,size(condis,2));

for con = 1:size(condis,2)

    for ch = 1:size(phases_condi{con},1)
        
        itpc{con}(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_condi{con}(ch,:,:,:)))),3));
        
    end

end

prefAngle = {};
for con = 1:size(condis,2)

    for ch = 1:size(phases_condi{con},1)
        
        for freq = 1:size(phases_condi{con},2)
        
            prefAngle{con}(ch,freq,:) = angle(mean(exp(1i*(squeeze(phases_condi{con}(ch,freq,:,:)))),2))';
        
         end

    end
end

%% ITPC grouped
itpc_grouped = cell(1,size(groupassign,2));

groups = cell(1,size(groupassign,2));
phases_grouped = cell(1,size(groupassign,2));
for grp = 1:size(groupassign,2)
    for con = 1:length(groupassign{grp}) 
    groups{grp} = cat(2,groups{grp},{condis{groupassign{grp}(con)}});
    phases_grouped{grp} = cat(4,phases_grouped{grp},phases_condi{groupassign{grp}(con)});
    end
end

for con = 1:size(itpc_grouped,2)

    for ch = 1:size(phases_condi{con},1)
        
        itpc_grouped{con}(ch,:,:) = abs(mean(exp(1i*(squeeze(phases_grouped{con}(ch,:,:,:)))),3));
        %itpc_grouped{2}(ch,:,:) = abs(mean(exp(1i*(squeeze(LWgroup(ch,:,:,:)))),3));

    end

end

prefAngle_grouped = {};
for con = 1:size(groupassign,2)

    for ch = 1:size(phases_condi{con},1)
        
        for freq = 1:size(phases_condi{con},2)
        
            prefAngle_grouped{con}(ch,freq,:) = angle(mean(exp(1i*(squeeze(phases_grouped{con}(ch,freq,:,:)))),2))';

        
         end

    end
end
%% make plots in probe layout for all conditions


frex = [1:151];
fran = [frex(1) 50];
xlim = [-25 75];
time = [-1000:1000];
clims = [0 .5];

prb = find(ismember(ProbeInfo.Areas, region));
poilayout = ProbeInfo.ProbeMaps{ProbeInfo.poi(prb)};
poilayoutT = reshape(poilayout',1,[]);
poilayoutT(ismember(poilayoutT, ProbeInfo.Ch_Remove{prb})) = 0;
poimapT = cell2mat(ProbeInfo.ProbeMaps(ProbeInfo.poi(prb)))';

for con = 1:size(condis,2)
figure()
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle([animal ' ' region '-Probe ' condis{con} ' ITPC'], 'Interpreter', 'none')


for chan = 1:numel(poilayoutT)
    loc = poilayoutT(chan);
    if loc > 0
    subplot(size(poimapT,1),size(poimapT,2),chan)
    contour(time, frex, squeeze(itpc{con}(chansuse == poimapT(chan),:,:)),25, 'Fill', 'on')
        set(gca,'ydir','normal', 'ylim', fran, 'clim',clims, 'XLim',xlim)%, 'clim',[-10 10])
        %colormap(slanCM('iceburn'))
        %colorbar('southoutside');
        title(['CH',num2str(poimapT(chan))])
        xlabel('Time (ms)')
        %yticks([1:10:50])
        %yticklabels([1 10 20 30 40 50])
        %xticks([0 500 1000])
        %xticklabels([0 1000 2000])
        xline(0,'LineWidth',1,'Color','r')
    %end
end

end
end
%% make plots in probe layout for GROUPED conditions

%choose your plotting variables
frex = [1:151]; %which frequencies are present in the data?
fran = [35 50]%[frex(1) 60]; %which frequencies to show?
xlim = [20 60]; %which timepoints (relative to stim) to show?
time = [-1000:1000]; %what are the timepoints (relative to stim) present in the data
clims = [0 .4]; %what color limits for ITPC to use?

prb = find(ismember(ProbeInfo.Areas, region));
poilayout = ProbeInfo.ProbeMaps{ProbeInfo.poi(prb)};
poilayoutT = reshape(poilayout',1,[]);
poilayoutT(ismember(poilayoutT, ProbeInfo.Ch_Remove{prb})) = 0;
poimapT = cell2mat(ProbeInfo.ProbeMaps(ProbeInfo.poi(prb)))';

for con = 1:size(groupassign,2)
figure()
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle([animal ' ' region '-Probe ' groupnames{con} ' ITPC'], 'Interpreter', 'none')


for chan = 1:numel(poilayoutT)
    loc = poilayoutT(chan);
    if loc > 0
    subplot(size(poimapT,1),size(poimapT,2),chan)
    contour(time, frex, squeeze(itpc_grouped{con}(chansuse == poimapT(chan),:,:)),50, 'Fill', 'on')
        set(gca,'ydir','normal', 'ylim', fran, 'clim',clims, 'XLim',xlim)%, 'clim',[-10 10])
        %colormap(slanCM('iceburn'))
        %colorbar('southoutside');
        title(['CH',num2str(poimapT(chan))])
        xlabel('Time (ms)')
        %yticks([1:10:50])
        %yticklabels([1 10 20 30 40 50])
        %xticks([0 500 1000])
        %xticklabels([0 1000 2000])
        xline(0,'LineWidth',1,'Color','r')
    %end
end

end
end
%% Look more closely at a specific channel for all individual conditions
chcheck = 14;
frex = [1:151];
fran = [frex(1) 85];
xlim = [-10 100];
time = [-1000:1000];
clims = [0 .5];

for con = 1:size(itpc,2)
    figure();
    contour(time, frex, squeeze(itpc{con}(chansuse == chcheck,:,:)),20, 'Fill', 'on')
        set(gca,'ydir','normal', 'ylim', fran, 'clim',clims, 'XLim',xlim)%, 'clim',[-10 10])
        %colormap(slanCM('iceburn'))
        %colorbar('southoutside');
        title(condis{con})
        xlabel('Time (ms)')
        %yticks([1:10:50])
        %yticklabels([1 10 20 30 40 50])
        %xticks([0 500 1000])
        %xticklabels([0 1000 2000])
        xline(0,'LineWidth',1,'Color','r')
    %end
end



%% Look more closely at a specific channel for GROUPED conditions

chcheck = 14;
frex = [1:151];
fran = [frex(1) 85];
xlim = [-10 100];
time = [-1000:1000];
clims = [0 .5];

for con = 1:size(itpc_grouped,2)
    figure();
    contour(time, frex, squeeze(itpc_grouped{con}(chansuse == chcheck,:,:)),50, 'Fill', 'on')
        set(gca,'ydir','normal', 'ylim', fran, 'clim',clims, 'XLim',xlim)%, 'clim',[-10 10])
        %colormap(slanCM('iceburn'))
        %colorbar('southoutside');
        title(groupnames{con})
        xlabel('Time (ms)')
        %yticks([1:10:50])
        %yticklabels([1 10 20 30 40 50])
        %xticks([0 500 1000])
        %xticklabels([0 1000 2000])
        xline(0,'LineWidth',1,'Color','r')
    %end
end


%%
% colors = bone(20);
% figure(); 
% for i = 1:20
%     hold on
% plot(squeeze(itpc_grouped{1}(14,30+i,:))','Color',colors(i,:))
% end
% hold off




% %% 
% chans = [20 21 22];
% fband = [20:45];
% twind = [1010:1025];
% corphases = [];
% %corphases = zeros(length(chans),size(phases_condi{1},4),2);
% for ch = 1:length(chans)
% 
% 
%    chdat = squeeze(cat(4,Lgroup(chans(ch),fband,twind,:),LWgroup(chans(ch),fband,twind,:)));
% 
%    for tr = 1:size(chdat,3)
% 
%         corphases(ch,tr,:) = reshape(squeeze(chdat(:,:,tr))',1,[]);
% 
%    end
% 
% end
% 
% corr_ch1 = corr(squeeze(corphases(1,:,:))');
% figure();
% heatmap(corr_ch1,"Colormap",parula(20));
% 
% t = clusterdata(squeeze(corphases(1,:,:)),maxclust=5)
% 
% 
% t = clusterdata(squeeze(corphases(1,:,:)),cutoff=5)
% 
% corr_ch2 = corr(squeeze(corphases(2,:,:))');
% figure();
% heatmap(corr_ch2,"Colormap",parula(20));
% 
% corr_ch3 = corr(squeeze(corphases(3,:,:))');
% figure();
% heatmap(corr_ch3,"Colormap",parula(20));
% 
% 
% % %% 
% figure();
% plot(squeeze(prefAngle_grouped{2}(chcheck,30,:)),'Color','k','LineWidth',3)
% hold on
% plot(squeeze(LWgroup(chcheck,30,:,:)),'Color','r')
% %plot(squeeze(phases_condi{1}(59,40,:,48)),'Color','c')
% 
% 
% figure();
% hold on
% for i=1:64
% chcheck = i
% twind = [1025:1075];
% %twind = [945:995];
% 
% for freq = 1:size(phases_condi{1},2)
% for tr = 1:size(phases_condi{1},4)
%     maxlag = round(1000/freq);
%     if maxlag > length(twind)/2
%         maxlag = round(length(twind)/2);
%     end
%     [tr_xcorr, lags] = xcorr(squeeze(prefAngle{1}(chcheck,freq,twind)),squeeze(phases_condi{1}(chcheck,freq,twind,tr)),maxlag);
%     tr_lagtoPref(freq,tr) = (lags(find(tr_xcorr==max(tr_xcorr),1)));
% end
% end
% 
%     subplot(8,8,i)
%     mean(tr_lagtoPref(30:40,:),1)
%     hist(ans,20)
%     xlim([-25 25])
% end
% heatmap(tr_lagtoPref)
% figure(); hist(ans)
% [a,b] = sort(ans)
% figure();
% heatmap(corr(tr_lagtoPref(:,:)),"Colormap",parula,'ColorLimits',[-1 1]);



%% xcorr sliding window for each ch/trial/freq

%make sure the index here is the group where ITPC is expected
refpref = prefAngle_grouped{1};


twind = 1025:1075; %time window for comparison to prefangle
fwind = 30:45; %frequency window for compatison to prefangle
tspan = 50 %span of time window in ms for which to check xcorr
tr_lagtoPref = {};
for con = 1:size(phases_grouped,2)
for ch = 1:length(chansuse)
     disp(['Ch' num2str(chansuse(ch))])
    for tr = 1:size(phases_grouped{con},4)
        lagtoPref = [];

        for freq = 1:length(fwind)
            freqi = fwind(freq);
            maxlag = round(1000/freqi);
            if maxlag > tspan/2
            maxlag = round(tspan/2);
            end
            pref = squeeze(refpref(ch,freqi,twind(1)-tspan*2:twind(end)+tspan*2))';
            tr_i = squeeze(phases_grouped{con}(ch,freqi,twind(1)-tspan*2:twind(end)+tspan*2,tr))';
            for t = 1:length(twind)
               
            tm = t+tspan*2;
            [r, lags] = xcorr(pref(tm-(tspan/2):tm+(tspan/2)),tr_i(tm-(tspan/2):tm+(tspan/2)),maxlag);
            
            lagtoPref(freqi,t) = lags(find(r==max(r),1));
            end
        end
        tr_lagtoPref{con}(ch,tr) = mean(abs(lagtoPref),'all');

    end
end
end

%ddcorr = {};
%dcorr_mean = [];
%for con = 1:size(phases_grouped,2)

%for ch = 1:size(lagtoPref{con},1)
%    for tr = 1:size(lagtoPref{con},2)
%        for freq = 1:size(lagtoPref,3)
%            dcorr_mean(ch,tr,freq) = mean(abs(lagtoPref{con}(ch,tr,freq,:)),4);
%            dd = squeeze(abs(lagtoPref{con}(ch,tr,freq,:))./dcorr_mean(ch,tr,freq));
%            ddcorr{con}(ch,tr,freq,:) = dd/max(dd);
%        end
%    end
%end

%end


%% check for significant differences between phase resetting

PhaseReset_signif = [];

for i = 1:length(chansuse)
    [p,h,stats] = ranksum(tr_lagtoPref{1}(i,:), tr_lagtoPref{2}(i,:));
    PhaseReset_signif(i,1) = p; %pval
    PhaseReset_signif(i,2) = h; %rej null 1=y, 0=n
end

%% save useful info in a struct


itpc_info.ITPC = itpc_grouped;
itpc_info.prefAngles = prefAngle_grouped;
itpc_info.dimfeatures = {'chans','freqs','time'};
itpc_info.stimtime = 1000;
itpc_info.roi_time = twind;
itpc_info.roi_freqs = fwind;
itpc_info.trialLagstoPref = tr_lagtoPref;
itpc_info.trials = trs;
itpc_info.condition_directory = superCondis_dir;
itpc_info.all_conditions = condis;
itpc_info.condition_grouping = groups;
itpc_info.condition_grouping_index = groupassign;
itpc_info.groupnames = groupnames;
itpc_info.reset_signif = PhaseReset_signif;
itpc_info.ch_index = chansuse;
itpc_info.region = region;



%%
folder = ['E:\Roy\Silcon Probe Data Followup Analyses\ITPC\' animal '_grouped']
    if isfolder(folder) == false
        mkdir(folder)
    end
save(fullfile([folder '\itpc_W_vs_W-TTX']),"itpc_info")
%%









%%

figure()
for i = 1:size(lagtoPref,2)
    subplot(9,8,i)
    contour(squeeze(ddcorr(59,i,:,:)),150,'fill','on')
end
figure(); contour(squeeze(ddcorr(59,2,:,:)),150,'fill','on')

figure()
for i = 1:size(lagtoPref,2)
    subplot(9,8,i)
    hist(reshape(squeeze(ddcorr(59,i,35:45,1025:1075)),1,[]),40)
    xline(.1)
end

figure()
hist(squeeze(mean(mean(ddcorr(59,:,35:45,1025:1075),3),4)),20)
xline(.1)

for ch = 1:size(lagtoPref,1)
ddmean_roi(ch,:) = squeeze(mean(mean(ddcorr(ch,:,35:45,1025:1075),3),4));
[B, Ind]  = sort(ddmean_roi(ch,:));
trRank_bychan(ch,:) = Ind;
end

figure();
heatmap(corr(trRank_bychan),'ColorLimits',[0 .5])

Y = pdist(trRank_bychan')
[a,b,outperm] = dendrogram(linkage(Y));

%% pca stuff

%reshape
phases_forpca = [];
for ch = 1:size(lagtoPref,1)

       chdat = squeeze(phases_condi{1}(ch,:,800:1010,:));

   for tr = 1:size(chdat,3)

        phases_forpca(ch,tr,:) = reshape(squeeze(chdat(:,:,tr))',1,[]);

   end

end


[coeff,score,latent,tsquared,explained,mu] = pca(squeeze(phases_forpca(4,:,:)));

figure();bar(explained)

pcloading_phase = [];
for pc = 1:size(coeff,2)
  pcloading_phase(pc,:,:) =  reshape(coeff(:,pc)',[],size(phases_condi{1},2))';
end

figure();
for pc = 1:size(pcloading_phase,1)
subplot(9,8,pc)
   contour(squeeze(pcloading_phase(pc,:,:)));
end


trPC_corr = [];
for pc = 1:size(pcloading_phase,1)
   trPC_corr(pc) = corr(squeeze(trRank_bychan(59,:))',squeeze(score(:,pc)));
end
figure(); bar(trPC_corr)

%reconstruct signal from positively correlated pcs
score_posPCs = score(:,trPC_corr>0);
coeff_posPCs = coeff(:,trPC_corr>0);

score_posPCs = score(:,[1,3,7,13,22]);
coeff_posPCs = coeff(:,[1,3,7,13,22]);

trials_reconstructed = (score_posPCs * coeff_posPCs') + mu;

%phases_forpca = [];
%for ch = 1:size(lagtoPref,1)

%       chdat = squeeze(phases_condi{1}(ch,:,800:1010,:));

   for tr = 1:size(trials_reconstructed,1)

        trials_recon_rect(tr,:,:) = reshape(squeeze(trials_reconstructed(tr,:)),[],size(phases_condi{1},2));

   end


figure();
for tr = 1:1:size(trials_reconstructed,1)
subplot(9,8,tr)
   contour(squeeze(trials_recon_rect(trRank_bychan(59,tr),:,1:150))',150,'fill','on');
end


%% load power and check that
power = [];
power_condi = cell(1,size(condis,2));
%tr2check = {[1:25]
            %[1:25]
            %[1:25]
            %[1:25]
            %[1:25]
            %};
for con = 1:length(condis)
    for file = 1:numel(superCondis_dir{con})
        dat = matfile(strjoin(superCondis_dir{con}(file)));
        %trs =  tr2check{file};
        disp(num2str(file))       
        %for shnk = 1:length(shank_groups_types)
        power = abs(dat.stim_tf(chansall,:,2000:4000,:)).^2;


        %phases_oi = cat(4,phases_oi,phases(:,:,:,trs));
        %phases(:,:,:,trs) = [];
        power_condi{con} = cat(4,power_condi{con},power(:,:,:,:));
        power = [];

        %phases_noi = cat(4,phases_noi,phases);
        %phases = [];
   end
end



%reshape
power_forpca = [];
for ch = 1:size(lagtoPref,1)

       chdat = squeeze(power_condi{1}(ch,1:25,800:1010,:));

   for tr = 1:size(chdat,3)

        power_forpca(ch,tr,:) = reshape(squeeze(chdat(:,:,tr))',1,[]);

   end

end


[coeff,score,latent,tsquared,explained,mu] = pca(squeeze(power_forpca(59,:,:)));

figure();bar(explained)

pcloading_power = [];
for pc = 1:size(coeff,2)
  pcloading_power(pc,:,:) =  reshape(coeff(:,pc)',[],size(chdat,1))';
end

figure();
for pc = 1:size(pcloading_power,1)
subplot(9,8,pc)
   contour(squeeze(pcloading_power(pc,:,:)),150,'fill','on');
end


trPC_corr = [];
for pc = 1:size(pcloading_phase,1)
   trPC_corr(pc) = corr(squeeze(trRank_bychan(59,:))',squeeze(score(:,pc)));
end
figure(); bar(trPC_corr)

%reconstruct signal from positively correlated pcs
score_posPCs = score(:,trPC_corr<-.1);
coeff_posPCs = coeff(:,trPC_corr<-.1);

score_posPCs = score(:,[13]);
coeff_posPCs = coeff(:,[13]);

trials_reconstructed = (score_posPCs * coeff_posPCs') + mu;

%phases_forpca = [];
%for ch = 1:size(lagtoPref,1)

%       chdat = squeeze(phases_condi{1}(ch,:,800:1010,:));
    trials_recon_rect = [];
   for tr = 1:size(trials_reconstructed,1)

        trials_recon_rect(tr,:,:) = reshape(squeeze(trials_reconstructed(tr,:)),[],size(chdat,1));

   end


figure();
for tr = 1:size(trials_reconstructed,1)
subplot(9,8,tr)
   contour(squeeze(trials_recon_rect(trRank_bychan(59,tr),:,1:25))',150,'fill','on');
end
