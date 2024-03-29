%%
clear all;
OVERWRITE_ = 0;

[fwFolder,anatomFolder,derivFolder,dataDir] = setmyenv();

subjects_set = [1:7 9:16];
rDMNDAN_AllSubj_hbo = zeros(40,40,length(subjects_set));
rDMNDAN_AllSubj_hbr = zeros(40,40,length(subjects_set));
nSubjs=length(subjects_set);

% flags and thresholds.
flags.macorrect = 'splineSG'; % 'none' or 'spline'
flags.bpfilt = 'image';% 'none' 'channel' or 'image'
flags.imagerecon = 'brain+scalp'; %'brain' or 'brain+scalp'
flags.rhoSD_ssThresh = 15;
flags.gsr = 'image';%'none','channel' or 'image'
flags.r_thresh = 0.7; % .r_thresh is the threshold for the clustering
flags.plot = 0; % .plot  flag to plot the brain correlation map
flags.p_thresh = 0; % . p_thresh is used to plot r values below that p-val (use 0 to plot all the correlations)
flags.clusteringType = 1; %0: no clustering, 1:Matlab, 2:David's algorithm
flags.task = 'WM';
flags.parcel_scheme = 'schaefer_comb';
flags.Adot_scaling_factor = 100;

pipeline_str = sprintf('WM-macor-%s_bpfilt-%s_imrec-%s_gsr-%s_clust-%i',...
    flags.macorrect,flags.bpfilt,flags.imagerecon,flags.gsr,flags.clusteringType);
fOut_map=sprintf('WM-Map_%s',pipeline_str);
pipelineDir = sprintf('%sPipeline-%s/',derivFolder,pipeline_str);
if ~exist(pipelineDir,'dir')
    mkdir(pipelineDir);
end
% to obtain the corresponding nodes in the mesh for each submask
[dmn_mask,dan_mask,mesh_brain,idx_select] = Parcellation_FC(anatomFolder,fwFolder);

%%
iSubj = 2;
%iRun = 1;
%for iSubj = 1:3
for iRun=1:8
    subject = num2str(subjects_set(iSubj));
    fprintf('Subject %s------------------------\n',subject);
    eval(sprintf("SnirfFilePathr%i = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-',flags.task,'_run-',num2str(%i),'_nirs.snirf']",iRun,iRun));
    %SnirfFilePathr2 = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-',flags.task,'_run-2_nirs.snirf'];
    if eval(sprintf("~exist(SnirfFilePathr%i,'file')",iRun))
        eval(sprintf("HbO_brain_chunkr%i = [];",iRun));
        eval(sprintf("HbR_brain_chunkr%i = [];",iRun));
        continue;
    end

    %%
    %if exist(eval(sprintf("[pipelineDir,fOut_map,'_r-',num2str(%i),'_hbr.png']",iRun)),'file')
    %    continue;
    %end
    %if ~exist(SnirfFilePathr1,'file')
    %    continue;
    %end

    %%
    eval(sprintf("[snirfObjr%i,dcObjr%i,dodObjr%i] = Preprocessing_FC(SnirfFilePathr%i,flags);",iRun,iRun,iRun,iRun));

    %%
    eval(sprintf("[HbO_brainr%i,HbR_brainr%i] = ImageReconstruction_FC(snirfObjr%i,dodObjr%i,dcObjr%i,fwFolder,flags);",iRun,iRun,iRun,iRun,iRun));
    % Do I want to visualize HbO from the Image recon?
    %checkImg_FC(fwFolder,HbO_brainr1,dodObjr1.time);
    %% Band pass filtering in image space
    if strcmp(flags.bpfilt,'image')
        %[y2] = image_BandpassFilt(y,hpf,lpf,fs)
        fs = mean(1./diff(snirfObjr1.data.time));
        eval(sprintf("HbO_brainr%i = image_BandpassFilt(HbO_brainr%i, 0.009, 0.080,fs);",iRun,iRun));
        eval(sprintf("HbR_brainr%i = image_BandpassFilt(HbR_brainr%i, 0.009, 0.080,fs);",iRun,iRun));
    end
    % Global signal regression in image space?
    if strcmp(flags.gsr,'image') && ~strcmp(flags.bpfilt,'none')
        eval(sprintf("HbO_brainr%i = GlobalRegression(HbO_brainr%i);",iRun,iRun));
        eval(sprintf("HbR_brainr%i = GlobalRegression(HbR_brainr%i);",iRun,iRun));
        % HbO_brain_chunkr2 = GlobalRegression(HbO_brain_chunkr2);
        % HbR_brain_chunkr2 = GlobalRegression(HbR_brain_chunkr2);
    end
    %%
    % use twindow.init_sec = -1 if you want to use the onset and duration
    % defined in the snirf file (assuming there is only one stimulus).
    % Otherwise, change the values (in sec.) according to your needs
    twindow.stim_name = 'baseline';
    twindow.init_sec = -1;
    twindow.offset_sec = 0;
    %twindow.init_sec = 30;
    %twindow.dur_sec = 180;
    eval(sprintf("[HbO_brain_chunkr%i] = ExctractChunk(HbO_brainr%i,snirfObjr%i,twindow,flags);",iRun,iRun,iRun));
    eval(sprintf("[HbR_brain_chunkr%i] = ExctractChunk(HbR_brainr%i,snirfObjr%i,twindow,flags);",iRun,iRun,iRun));
    %[HbO_brain_chunkr1] = ExctractChunk(HbO_brainr1,snirfObjr1,twindow,flags);
    %[HbR_brain_chunkr1] = ExctractChunk(HbR_brainr1,snirfObjr1,twindow,flags);


   
end
%%
% Do we want to concatenate runs?
[HbO_brain_r1r8] = [HbO_brain_chunkr1-mean(HbO_brain_chunkr1);
    HbO_brain_chunkr2-mean(HbO_brain_chunkr2);
    HbO_brain_chunkr3-mean(HbO_brain_chunkr3);
    HbO_brain_chunkr4-mean(HbO_brain_chunkr4);
    HbO_brain_chunkr5-mean(HbO_brain_chunkr5);
    HbO_brain_chunkr6-mean(HbO_brain_chunkr6);
    HbO_brain_chunkr7-mean(HbO_brain_chunkr7);
    HbO_brain_chunkr8-mean(HbO_brain_chunkr8);];
[HbR_brain_r1r8] = [HbR_brain_chunkr1-mean(HbR_brain_chunkr1);
    HbR_brain_chunkr2-mean(HbR_brain_chunkr2);
    HbR_brain_chunkr3-mean(HbR_brain_chunkr3);
    HbR_brain_chunkr4-mean(HbR_brain_chunkr4);
    HbR_brain_chunkr5-mean(HbR_brain_chunkr5);
    HbR_brain_chunkr6-mean(HbR_brain_chunkr6);
    HbR_brain_chunkr7-mean(HbR_brain_chunkr7);
    HbR_brain_chunkr8-mean(HbR_brain_chunkr8)];
% seedsFolder ='/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/DataRSFC_Analysis/derivatives/rsfc/Pipeline-RS-macor-spline_bpfilt-image_imrec-brain+scalp_gsr-image_clust-1/';
% fOut_seeds ='reliability_RS-macor-spline_bpfilt-image_imrec-brain+scalp_gsr-image_clust-1';
% load([seedsFolder,fOut_seeds,'.mat'],'rDMNDAN_AllSubj_hbo','rDMNDAN_AllSubj_hbr','flags', ...
%         'dmn_seeds_NP_ts','dan_seeds_NP_ts','dmn_seeds_ts','dan_seeds_ts',...
%         'dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr');
% 
% dmn_improv_hbo = dmn_improv_hbo{iSubj};
% dmn_improv_hbr = dmn_improv_hbr{iSubj};
% dan_improv_hbo = dan_improv_hbo{iSubj};
% dan_improv_hbr = dan_improv_hbr{iSubj};

%%
[dmn_mask_hbo] = Clustering_FC(HbO_brain_r1r8,dmn_mask,flags);
[dmn_mask_hbr] = Clustering_FC(HbR_brain_r1r8,dmn_mask,flags);

%%
[dan_mask_hbo] = Clustering_FC(HbO_brain_r1r8,dan_mask,flags);
[dan_mask_hbr] = Clustering_FC(HbR_brain_r1r8,dan_mask,flags);

%%
[~,dmn_improv_hbo,dan_improv_hbo,~] = ClusterSelection_FC(dmn_mask_hbo,dan_mask_hbo,HbO_brain_r1r8,flags);
[~,dmn_improv_hbr,dan_improv_hbr,~] = ClusterSelection_FC(dmn_mask_hbr,dan_mask_hbr,HbR_brain_r1r8,flags);

%%
% HbO
BrainMaps_hbo = zeros(length(idx_select),1*(length(dmn_improv_hbo)+length(dan_improv_hbo)));
dmn_hbo_ts = zeros(size(HbO_brain_r1r8,1),16);
dan_hbo_ts = zeros(size(HbO_brain_r1r8,1),4);
dmn_hbr_ts = dmn_hbo_ts;
dan_hbr_ts = dan_hbo_ts;
for iSubmask=1:length(dmn_improv_hbo)
    %obtain the correlation brain map after preprocessing and by using the
    %seed passed as the first argument.
    [dmn_hbo_ts(:,iSubmask),hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_r1r8,flags.p_thresh,flags.plot);    
    BrainMaps_hbo(:,iSubmask) = A_select1;
    fprintf('(hbo)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbo));
end
for iSubmask=1:length(dan_improv_hbo)
    [dan_hbo_ts(:,iSubmask),hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_r1r8,flags.p_thresh,flags.plot);
    BrainMaps_hbo(:,length(dmn_improv_hbo)+iSubmask) = A_select1;
    fprintf('(hbo)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbo));
end

f=figure;
imagesc(corrcoef(BrainMaps_hbo),[-1,1]);
colormap("jet");
colorbar();
title({sprintf('Subject %s DMN-DAN HbO Run 1-8',subject),twindow.stim_name,pipeline_str},'Interpreter','none');
saveas(f,[pipelineDir,fOut_map,'_Subj-',num2str(iSubj),'_WMseed_',twindow.stim_name,'_hbo.png']);
close(f);

f=figure;
plot(mean(dmn_hbo_ts,2),'-','Color',[0 0 1],'LineWidth',2.5);
hold on;
plot(mean(dan_hbo_ts,2),'-','Color',[1 0 0],'LineWidth',2.5);
plot(dmn_hbo_ts,'-','Color',[0.7 0.7 1]);
plot(dan_hbo_ts,'-','Color',[1 0.7 0.7]);
xlabel('Samples');
ylabel('\Delta HbO')
legend({'avg DMN','avg DAN'});
title({sprintf('Subject %s DMN-DAN HbO time courses',subject),twindow.stim_name,pipeline_str},'Interpreter','none');
saveas(f,[pipelineDir,fOut_map,'_Subj-',num2str(iSubj),'_WMseed_',twindow.stim_name,'_hbo_ts.png']);
close(f);

% hbr
BrainMaps_hbr = zeros(length(idx_select),1*(length(dmn_improv_hbr)+length(dan_improv_hbr)));
for iSubmask=1:length(dmn_improv_hbr)
    [dmn_hbr_ts(:,iSubmask),hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_r1r8,flags.p_thresh,flags.plot);
    BrainMaps_hbr(:,iSubmask) = A_select1;
    fprintf('(hbr)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbr));
end
for iSubmask=1:length(dan_improv_hbr)
    [dan_hbr_ts(:,iSubmask),hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_r1r8,flags.p_thresh,flags.plot);
    BrainMaps_hbr(:,length(dmn_improv_hbr)+iSubmask) = A_select1;
    fprintf('(hbr)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbr));
end

f=figure;
imagesc(corrcoef(BrainMaps_hbr),[-1,1]);
colormap("jet");
colorbar();
title({sprintf('Subject %s DMN-DAN HbR Run 1-8',subject),twindow.stim_name,pipeline_str},'Interpreter','none');
saveas(f,[pipelineDir,fOut_map,'_Subj-',num2str(iSubj),'_WMseed_',twindow.stim_name,'_hbr.png']);
close(f);

f=figure;
plot(mean(dmn_hbr_ts,2),'-','Color',[0 0 1],'LineWidth',2.5);
hold on;
plot(mean(dan_hbr_ts,2),'-','Color',[1 0 0],'LineWidth',2.5);
plot(dmn_hbr_ts,'-','Color',[0.7 0.7 1]);
plot(dan_hbr_ts,'-','Color',[1 0.7 0.7]);
xlabel('Samples');
ylabel('\Delta HbR')
legend({'avg DMN','avg DAN'});
title({sprintf('Subject %s DMN-DAN HbR time courses',subject),twindow.stim_name,pipeline_str},'Interpreter','none');
saveas(f,[pipelineDir,fOut_map,'_Subj-',num2str(iSubj),'_WMseed_',twindow.stim_name,'_hbr_ts.png']);
close(f);

%end