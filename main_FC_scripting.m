%%
fwFolder = 'C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\fw\';
anatomFolder = 'C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\probe_10MPhotons\anatomical\';

subjects_set = [1:7 9:16];
rDMNDAN_AllSubj_hbo = zeros(40,40,length(subjects_set));
rDMNDAN_AllSubj_hbr = zeros(40,40,length(subjects_set));
nSubjs=length(subjects_set);

% flags and thresholds.
flags.imagerecon = 'brain'; %'brain' or 'brain+scalp' Image reconstruction
flags.rhoSD_ssThresh = 15; %threshold (mm) to identify and remove short-separation channles.
flags.gsr = 'none'; %Space to perform global signal regression: 'none','channel' or 'image'
flags.r_thresh = 0.7; % value required for Clustering 1 (Matlab clustering)
flags.plot=1;  % whether to plot or not the compuytes brain correlation map.
flags.p_thresh = 0; %is used to plot r values below .p_thresh (use 0 to plot all the correlations)
flags.clusteringType = 1; %1:Matlab, 2:David's algorithm
% to obtain the corresponding nodes in the mesh for each submask
[dmn_mask,dan_mask,mesh_brain,idx_select] = Parcellation_FC(anatomFolder,fwFolder);

iSubj = 1;
%for iSubj = 1:nSubjs
subject = num2str(subjects_set(iSubj));
fprintf('Subject %s------------------------\n',subject);
SnirfFilePathr1 = ['C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\sub-',subject,'\nirs\sub-',subject,'_task-RS_run-1_nirs.snirf'];
SnirfFilePathr2 = ['C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\sub-',subject,'\nirs\sub-',subject,'_task-RS_run-2_nirs.snirf'];

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
twindow.offset_sec = 60;
%twindow.init_sec = 30;
%twindow.dur_sec = 180;
[HbO_brain_chunkr1] = ExctractChunk(HbO_brainr1,snirfObjr1,twindow);
[HbR_brain_chunkr1] = ExctractChunk(HbR_brainr1,snirfObjr1,twindow);
%[HbO_brain_chunkr2] = ExctractChunk(HbO_brainr2,snirfObjr2,twindow);
%[HbR_brain_chunkr2] = ExctractChunk(HbR_brainr2,snirfObjr2,twindow);

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
for iSubmask=1:length(dmn_improv_hbo)
    %obtain the correlation brain map after preprocessing and by using the
    %seed passed as the first argument.
    [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);    
    fprintf('(hbo)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbo));
end
for iSubmask=1:length(dan_improv_hbo)
    [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);
    fprintf('(hbo)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbo));
end

% hbr
for iSubmask=1:length(dmn_improv_hbr)
    [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
    fprintf('(hbr)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbr));
end
for iSubmask=1:length(dan_improv_hbr)
    [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
    fprintf('(hbr)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbr));
end
