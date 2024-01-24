%%
clear all;
OVERWRITE_ = 1;

%[fwFolder,anatomFolder,derivFolder,dataDir] = setmyenv();
maxNumCompThreads(32);
addpath('./cbrewer/');
if isempty(which('AtlasViewerGUI'))
    mypwd=pwd;cd('../atlasviewer_repo/');setpaths;cd(mypwd);
end
if isempty(which('Homer3'))
    mypwd=pwd;cd('../homer3_repo/');setpaths;cd(mypwd);
end
fwFolder = '/usr2/postdoc/smontero/Regina_files/probe/fw/';
anatomFolder = '/usr2/postdoc/smontero/Regina_files/probe/anatomical/';
derivFolder = '/usr2/postdoc/smontero/Regina_files/Data/Jan2024/derivatives/rsfc/';
dataDir = '/usr2/postdoc/smontero/Regina_files/Data/Jan2024/';

subjects_set = [1:7 9:16];
rDMNDAN_AllSubj_hbo = zeros(40,40,length(subjects_set));
rDMNDAN_AllSubj_hbr = zeros(40,40,length(subjects_set));
nSubjs=length(subjects_set);

% flags and thresholds.
flags.macorrect = 'spline'; % 'none' or 'spline'
flags.bpfilt = 'channel';% 'none' 'channel' or 'image'
flags.imagerecon = 'brain+scalp'; %'brain' or 'brain+scalp'
flags.rhoSD_ssThresh = 15;
flags.gsr = 'image';%'none','channel' or 'image'
flags.r_thresh = 0.7; % .r_thresh is the threshold for the clustering
flags.plot=0; % .plot  flag to plot the brain correlation map
flags.p_thresh = 0; % . p_thresh is used to plot r values below that p-val (use 0 to plot all the correlations)
flags.clusteringType = 1; %1:Matlab, 2:David's algorithm
flags.task = 'noexo';
pipeline_str = sprintf('noexo-macor-%s_bpfilt-%s_imrec-%s_gsr-%s_clust-%i',...
    flags.macorrect,flags.bpfilt,flags.imagerecon,flags.gsr,flags.clusteringType);
fOut_map=sprintf('noexo-Map_%s',pipeline_str);
pipelineDir = sprintf('%sPipeline-%s/',derivFolder,pipeline_str);
if ~exist(pipelineDir,'dir')
    mkdir(pipelineDir);
end
% to obtain the corresponding nodes in the mesh for each submask
[dmn_mask,dan_mask,mesh_brain,idx_select] = Parcellation_FC(anatomFolder,fwFolder);

%%
iSubj = 1;
%iRun = 1;
%for iSubj = 1:nSubjs
for iRun=1:8
subject = num2str(subjects_set(iSubj));
fprintf('Subject %s------------------------\n',subject);
SnirfFilePathr1 = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-',flags.task,'_run-',num2str(iRun),'_nirs.snirf'];
%SnirfFilePathr2 = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-',flags.task,'_run-2_nirs.snirf'];

%%
if exist([pipelineDir,fOut_map,'_r-',num2str(iRun),'_hbr.png'],'file')
    continue;
end
if ~exist(SnirfFilePathr1,'file')
    continue;
end

%%
[snirfObjr1,dcObjr1,dodObjr1] = Preprocessing_FC(SnirfFilePathr1,flags);

%%
[HbO_brainr1,HbR_brainr1] = ImageReconstruction_FC(snirfObjr1,dodObjr1,dcObjr1,fwFolder,flags);
% Do I want to visualize HbO from the Image recon?
checkImg_FC(fwFolder,HbO_brainr1,dodObjr1.time);

%%
% use twindow.init_sec = -1 if you want to use the onset and duration
% defined in the snirf file (assuming there is only one stimulus).
% Otherwise, change the values (in sec.) according to your needs
twindow.init_sec = -1;
twindow.offset_sec = 0;
%twindow.init_sec = 30;
%twindow.dur_sec = 180;
[HbO_brain_chunkr1] = ExctractChunk(HbO_brainr1,snirfObjr1,twindow,flags);
[HbR_brain_chunkr1] = ExctractChunk(HbR_brainr1,snirfObjr1,twindow,flags);
%[HbO_brain_chunkr2] = ExctractChunk(HbO_brainr2,snirfObjr2,twindow);
%[HbR_brain_chunkr2] = ExctractChunk(HbR_brainr2,snirfObjr2,twindow);


% Band pass filtering in image space
if strcmp(flags.bpfilt,'image')
    %[y2] = image_BandpassFilt(y,hpf,lpf,fs)
    fs = mean(1./diff(snirfObjr1.data.time));
    HbO_brain_chunkr1 = image_BandpassFilt(HbO_brain_chunkr1, 0.009, 0.080,fs);
    HbR_brain_chunkr1 = image_BandpassFilt(HbR_brain_chunkr1, 0.009, 0.080,fs);
    % HbO_brain_chunkr2 = image_BandpassFilt(HbO_brain_chunkr2, 0.009, 0.080,fs);
    % HbR_brain_chunkr2 = image_BandpassFilt(HbR_brain_chunkr2, 0.009, 0.080,fs);
end
% Global signal regression in image space?
if strcmp(flags.gsr,'image') && ~strcmp(flags.bpfilt,'none')
    HbO_brain_chunkr1 = GlobalRegression(HbO_brain_chunkr1);
    HbR_brain_chunkr1 = GlobalRegression(HbR_brain_chunkr1);
    % HbO_brain_chunkr2 = GlobalRegression(HbO_brain_chunkr2);
    % HbR_brain_chunkr2 = GlobalRegression(HbR_brain_chunkr2);
end
%%
% Do we want to concatenate r1 and r2?
%[HbO_brain_r1r2] = [HbO_brain_chunkr1;HbO_brain_chunkr2];
%[HbR_brain_r1r2] = [HbR_brain_chunkr1;HbR_brain_chunkr2];
%%
[dmn_mask_hbo] = Clustering_FC(HbO_brain_chunkr1,dmn_mask,flags);
[dmn_mask_hbr] = Clustering_FC(HbR_brain_chunkr1,dmn_mask,flags);

%%
[dan_mask_hbo] = Clustering_FC(HbO_brain_chunkr1,dan_mask,flags);
[dan_mask_hbr] = Clustering_FC(HbR_brain_chunkr1,dan_mask,flags);

%%
[~,dmn_improv_hbo,dan_improv_hbo,~] = ClusterSelection_FC(dmn_mask_hbo,dan_mask_hbo,HbO_brain_chunkr1,flags);
[~,dmn_improv_hbr,dan_improv_hbr,~] = ClusterSelection_FC(dmn_mask_hbr,dan_mask_hbr,HbR_brain_chunkr1,flags);

%%
% HbO
BrainMaps_hbo = zeros(length(idx_select),1*(length(dmn_improv_hbo)+length(dan_improv_hbo)));
for iSubmask=1:length(dmn_improv_hbo)
    %obtain the correlation brain map after preprocessing and by using the
    %seed passed as the first argument.
    [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);    
    BrainMaps_hbo(:,iSubmask) = A_select1;
    fprintf('(hbo)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbo));
end
for iSubmask=1:length(dan_improv_hbo)
    [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);
    BrainMaps_hbo(:,length(dmn_improv_hbo)+iSubmask) = A_select1;
    fprintf('(hbo)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbo));
end

f=figure;
imagesc(corrcoef(BrainMaps_hbo),[-1,1]);
colormap("jet");
colorbar();
title({sprintf('Subject %s DMN-DAN HbO Run %i',subject,iRun),pipeline_str},'Interpreter','none');
saveas(f,[pipelineDir,fOut_map,'_r-',num2str(iRun),'_hbo.png']);
close(f);

% hbr
BrainMaps_hbr = zeros(length(idx_select),1*(length(dmn_improv_hbr)+length(dan_improv_hbr)));
for iSubmask=1:length(dmn_improv_hbr)
    [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
    BrainMaps_hbr(:,iSubmask) = A_select1;
    fprintf('(hbr)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbr));
end
for iSubmask=1:length(dan_improv_hbr)
    [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
    BrainMaps_hbr(:,length(dmn_improv_hbr)+iSubmask) = A_select1;
    fprintf('(hbr)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbr));
end

f=figure;
imagesc(corrcoef(BrainMaps_hbr),[-1,1]);
colormap("jet");
colorbar();
title({sprintf('Subject %s DMN-DAN HbR Run %i',subject,iRun),pipeline_str},'Interpreter','none');
saveas(f,[pipelineDir,fOut_map,'_r-',num2str(iRun),'_hbr.png']);
close(f);
end