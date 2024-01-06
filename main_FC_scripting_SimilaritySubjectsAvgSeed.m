%%
fwFolder = 'C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\fw\';
anatomFolder = 'C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\probe_10MPhotons\anatomical\';

subjects_set = [1:7 9:16];
BrainMaps_hbo =[];
BrainMaps_hbr = [];
HbO_brain_r1r2_all = [];
HbR_brain_r1r2_all = [];
seed_dmn_hbo = [];
seed_dan_hbo = [];
seed_dmn_hbr = [];
seed_dan_hbr = [];
rDMNDAN_AllSubj_hbo = zeros(40,40,length(subjects_set));
rDMNDAN_AllSubj_hbr = zeros(40,40,length(subjects_set));
nSubjs=length(subjects_set);

% flags and thresholds.
% flags and thresholds
flags.imagerecon = 'brain'; %'brain' or 'brain+scalp'
flags.rhoSD_ssThresh = 15;
flags.gsr = 'image';%'none','channel' or 'image'
flags.r_thresh = 0.7; % .r_thresh is the threshold for the clustering
flags.plot=0; % .plot  flag to plot the brain correlation map
flags.p_thresh = 0; % . p_thresh is used to plot r values below that p-val (use 0 to plot all the correlations)
flags.clusteringType = 1; %1:Matlab, 2:David's algorithm
[dmn_mask,dan_mask,mesh_brain,idx_select] = Parcellation_FC(anatomFolder,fwFolder);
fnameoutput=sprintf('brainMaps_seed_dmndan_superSeed_imrecon-%s_gsr-%s_clust-%i',flags.imagerecon,flags.gsr,flags.clusteringType);

for iSubj = 1:nSubjs
    subject = num2str(subjects_set(iSubj));
    fprintf('==============================\n');
    fprintf('Subject %s\n',subject);
    fprintf('==============================\n');
    SnirfFilePathr1 = ['C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\sub-',subject,'\nirs\sub-',subject,'_task-RS_run-1_nirs.snirf'];
    SnirfFilePathr2 = ['C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\sub-',subject,'\nirs\sub-',subject,'_task-RS_run-2_nirs.snirf'];

    %%-
    [snirfObjr1,dcObjr1,dodObjr1] = Preprocessing_FC(SnirfFilePathr1,flags);
    [snirfObjr2,dcObjr2,dodObjr2] = Preprocessing_FC(SnirfFilePathr2,flags);

    %%-
    [HbO_brainr1,HbR_brainr1] = ImageReconstruction_FC(snirfObjr1,dodObjr1,dcObjr1,fwFolder,flags);
    [HbO_brainr2,HbR_brainr2] = ImageReconstruction_FC(snirfObjr2,dodObjr2,dcObjr2,fwFolder,flags);
   

    %%-
    % use twindow.init_sec = -1 if you want to use the onset and duration
    % defined in the snirf file (assuming there is only one stimulus).
    % Otherwise, change the values (in sec.) according to your needs
    twindow.init_sec = -1;
    twindow.offset_sec = 60;
    %twindow.init_sec = 30;
    %twindow.dur_sec = 180;
    [HbO_brain_chunkr1] = ExctractChunk(HbO_brainr1,snirfObjr1,twindow);
    [HbO_brain_chunkr2] = ExctractChunk(HbO_brainr2,snirfObjr2,twindow);
    [HbR_brain_chunkr1] = ExctractChunk(HbR_brainr1,snirfObjr1,twindow);
    [HbR_brain_chunkr2] = ExctractChunk(HbR_brainr2,snirfObjr2,twindow);

    %%-
    % Concatenating both runs for allsubjects
    [HbO_brain_r1r2_all] = [HbO_brain_r1r2_all;HbO_brain_chunkr1;HbO_brain_chunkr2];
    [HbR_brain_r1r2_all] = [HbR_brain_r1r2_all;HbR_brain_chunkr1;HbR_brain_chunkr2];
end
%%-
[dmn_mask_hbo] = Clustering_FC(HbO_brain_r1r2_all,dmn_mask,flags);
[dmn_mask_hbr] = Clustering_FC(HbR_brain_r1r2_all,dmn_mask,flags);

%%-
[dan_mask_hbo] = Clustering_FC(HbO_brain_r1r2_all,dan_mask,flags);
[dan_mask_hbr] = Clustering_FC(HbR_brain_r1r2_all,dan_mask,flags);

%%-
[~,dmn_improv_hbo,dan_improv_hbo,~] = ClusterSelection_FC(dmn_mask_hbo,dan_mask_hbo,HbO_brain_r1r2_all,flags);
[~,dmn_improv_hbr,dan_improv_hbr,~] = ClusterSelection_FC(dmn_mask_hbr,dan_mask_hbr,HbR_brain_r1r2_all,flags);

%%-
dmn_avg_hbo = dmn_improv_hbo(1);
dmn_avg_hbr = dmn_improv_hbr(1);

%%
for iSubmask=1:length(dmn_improv_hbo)
    submask=dmn_improv_hbo(iSubmask);
    dmn_avg_hbo.name = 'AverageSubmask';
    dmn_avg_hbo.vertices = [dmn_avg_hbo.vertices; submask.vertices];
    dmn_avg_hbo.vertices_index = [dmn_avg_hbo.vertices_index; submask.vertices_index];
    dmn_avg_hbo.groups = [dmn_avg_hbo.groups; submask.groups];
    dmn_avg_hbo.mask_subsetseed = [dmn_avg_hbo.mask_subsetseed; submask.mask_subsetseed];

    submask=dmn_improv_hbr(iSubmask);
    dmn_avg_hbr.name = 'AverageSubmask';
    dmn_avg_hbr.vertices = [dmn_avg_hbr.vertices; submask.vertices];
    dmn_avg_hbr.vertices_index = [dmn_avg_hbr.vertices_index; submask.vertices_index];
    dmn_avg_hbr.groups = [dmn_avg_hbr.groups; submask.groups];
    dmn_avg_hbr.mask_subsetseed = [dmn_avg_hbr.mask_subsetseed; submask.mask_subsetseed];
end

%dan
dan_avg_hbo = dan_improv_hbo(1);
dan_avg_hbr = dan_improv_hbr(1);

for iSubmask=1:length(dan_improv_hbo)
    submask=dan_improv_hbo(iSubmask);
    dan_avg_hbo.name = 'AverageSubmask';
    dan_avg_hbo.vertices = [dan_avg_hbo.vertices; submask.vertices];
    dan_avg_hbo.vertices_index = [dan_avg_hbo.vertices_index; submask.vertices_index];
    dan_avg_hbo.groups = [dan_avg_hbo.groups; submask.groups];
    dan_avg_hbo.mask_subsetseed = [dan_avg_hbo.mask_subsetseed; submask.mask_subsetseed];

    submask=dan_improv_hbr(iSubmask);
    dan_avg_hbr.name = 'AverageSubmask';
    dan_avg_hbr.vertices = [dan_avg_hbr.vertices; submask.vertices];
    dan_avg_hbr.vertices_index = [dan_avg_hbr.vertices_index; submask.vertices_index];
    dan_avg_hbr.groups = [dan_avg_hbr.groups; submask.groups];
    dan_avg_hbr.mask_subsetseed = [dan_avg_hbr.mask_subsetseed; submask.mask_subsetseed];
end

%%-
for iSubj = 1:nSubjs
    subject = num2str(subjects_set(iSubj));
    fprintf('==============================\n');
    fprintf('Subject %s\n',subject);
    fprintf('==============================\n');
    SnirfFilePathr1 = ['C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\sub-',subject,'\nirs\sub-',subject,'_task-RS_run-1_nirs.snirf'];
    SnirfFilePathr2 = ['C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\sub-',subject,'\nirs\sub-',subject,'_task-RS_run-2_nirs.snirf'];

    %%-
    [snirfObjr1,dc1,dodObjr1] = Preprocessing_FC(SnirfFilePathr1,flags);
    [snirfObjr2,dc2,dodObjr2] = Preprocessing_FC(SnirfFilePathr2,flags);

    %%-
    [HbO_brainr1,HbR_brainr1] = ImageReconstruction_FC(snirfObjr1,dodObjr1,dcObjr1,fwFolder,flags);
    [HbO_brainr2,HbR_brainr2] = ImageReconstruction_FC(snirfObjr2,dodObjr2,dcObjr2,fwFolder,flags);

    %%-
    % use twindow.init_sec = -1 if you want to use the onset and duration
    % defined in the snirf file (assuming there is only one stimulus).
    % Otherwise, change the values (in sec.) according to your needs
    twindow.init_sec = -1;
    twindow.offset_sec = 60;
    %twindow.init_sec = 30;
    %twindow.dur_sec = 180;
    [HbO_brain_chunkr1] = ExctractChunk(HbO_brainr1,snirfObjr1,twindow);
    [HbO_brain_chunkr2] = ExctractChunk(HbO_brainr2,snirfObjr2,twindow);
    [HbR_brain_chunkr1] = ExctractChunk(HbR_brainr1,snirfObjr1,twindow);
    [HbR_brain_chunkr2] = ExctractChunk(HbR_brainr2,snirfObjr2,twindow);

    %%-
    % Do we want to concatenate r1 and r2?
    [HbO_brain_r1r2] = [HbO_brain_chunkr1;HbO_brain_chunkr2];
    [HbR_brain_r1r2] = [HbR_brain_chunkr1;HbR_brain_chunkr2];

    %%-
    % HbO
    %BrainMaps_hbo = zeros(length(idx_select),(length(dmn_avg_hbo)+length(dan_avg_hbo)));
    [seed_dmn_hbo(:,iSubj),hmap1,A_select1] = CorrelationBrainMap_FC(dmn_avg_hbo,mesh_brain,idx_select,HbO_brain_r1r2,flags.p_thresh,flags.plot);
    BrainMaps_hbo(:,iSubj,1) = A_select1;
    [seed_dan_hbo(:,iSubj),hmap1,A_select1] = CorrelationBrainMap_FC(dan_avg_hbo,mesh_brain,idx_select,HbO_brain_r1r2,flags.p_thresh,flags.plot);
    BrainMaps_hbo(:,iSubj,2) = A_select1;
    fprintf('(hbo)DMN,DAN submask\n');

    % hbr
    %BrainMaps_hbr = zeros(length(idx_select),(length(dmn_avg_hbr)+length(dan_avg_hbr)));
    [seed_dmn_hbr(:,iSubj),hmap1,A_select1] = CorrelationBrainMap_FC(dmn_avg_hbr,mesh_brain,idx_select,HbR_brain_r1r2,flags.p_thresh,flags.plot);
    BrainMaps_hbr(:,iSubj,1) = A_select1;
    [seed_dan_hbr(:,iSubj),hmap1,A_select1] = CorrelationBrainMap_FC(dan_avg_hbr,mesh_brain,idx_select,HbR_brain_r1r2,flags.p_thresh,flags.plot);
    BrainMaps_hbr(:,iSubj,2) = A_select1;
    fprintf('(hbr)DMN,DAN submask\n');
    %
    % for iSubmask=1:length(dmn_avg_hbr)
    %     [hmap1,A_select1] = CorrelationBrainMap_FC(dmn_avg_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_r1r2,flags.p_thresh,flags.plot);
    %    % [hmap2,A_select2] = CorrelationBrainMap_FC(dmn_avg_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr2,flags.p_thresh,flags.plot);
    %     BrainMaps_hbr(:,iSubmask) = A_select1;
    %     %BrainMaps_hbr(:,length(dmn_avg_hbr)+length(dan_avg_hbo)+iSubmask) = A_select2;
    %     fprintf('(hbr)DMN submask %i of %i\n',iSubmask,length(dmn_avg_hbr));
    % end
    % for iSubmask=1:length(dan_avg_hbr)
    %     [hmap1,A_select1] = CorrelationBrainMap_FC(dan_avg_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_r1r2,flags.p_thresh,flags.plot);
    %     %[hmap2,A_select2] = CorrelationBrainMap_FC(dan_avg_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr2,flags.p_thresh,flags.plot);
    %     BrainMaps_hbr(:,length(dmn_avg_hbr)+iSubmask) = A_select1;
    %     %BrainMaps_hbr(:,2*length(dmn_avg_hbr)+length(dan_avg_hbr)+iSubmask) = A_select2;
    %     fprintf('(hbr)DAN submask %i of %i\n',iSubmask,length(dan_avg_hbr));
    % end
end
save([fnameoutput,'.mat'],'BrainMaps_hbo','BrainMaps_hbr','seed_dmn_hbo','seed_dan_hbo','seed_dmn_hbr','seed_dan_hbr','flags');
%%
% Similarity
f=figure();
subplot(2,3,1);
[r,p]=corrcoef([BrainMaps_hbo(:,:,1)]);
imagesc(r,[-1 1]);
colormap('jet');
colorbar();
title('DMN-DMN');
ylabel({'HbO','Subjects DMN'});

subplot(2,3,2);
[r,p]=corrcoef([BrainMaps_hbo(:,:,2)]);
imagesc(r,[-1 1]);
colormap('jet');
colorbar();
xlabel('Subjects DAN');
ylabel('Subjects DAN');
title('DAN-DAN');

subplot(2,3,3);
[r,p]=corrcoef([BrainMaps_hbo(:,:,1),BrainMaps_hbo(:,:,2)]);
imagesc(r(1:15,16:end),[-1 1]);
colormap('jet');
colorbar();
xlabel('Subjects DMN');
ylabel('Subjects DAN');
title('DMN-DAN');
%figure(10), imagesc(r,[-1 1]); colormap("jet")

subplot(2,3,4);
[r,p]=corrcoef([BrainMaps_hbr(:,:,1)]);
imagesc(r,[-1 1]);
colormap('jet');
colorbar();
title('DMN-DMN');
xlabel('Subjects DMN');
ylabel({'HbR','Subjects DMN'});

subplot(2,3,5);
[r,p]=corrcoef([BrainMaps_hbr(:,:,2)]);
imagesc(r,[-1 1]);
colormap('jet');
colorbar();
title('DAN-DAN');
xlabel('Subjects DAN');
ylabel('Subjects DAN');

subplot(2,3,6);
[r,p]=corrcoef([BrainMaps_hbr(:,:,1),BrainMaps_hbr(:,:,2)]);
imagesc(r(1:15,16:end),[-1 1]);
colormap('jet');
colorbar();
xlabel('Subjects DMN');
ylabel('Subjects DAN');
title('DMN-DAN');
%figure(10), imagesc(r,[-1 1]); colormap("jet")
f.Position = [10         10        1049         507];
saveas(f,[fnameoutput,'_SubjectsSimilarity.png']);
%%-
% Similarity boxplots
fooHbO_dmn = corrcoef([BrainMaps_hbo(:,:,1)]);
fooHbO_dan = corrcoef([BrainMaps_hbo(:,:,2)]);
fooHbR_dmn = corrcoef([BrainMaps_hbr(:,:,1)]);
fooHbR_dan = corrcoef([BrainMaps_hbr(:,:,2)]);
fooHbO_dmndan = corrcoef([BrainMaps_hbo(:,:,1),BrainMaps_hbo(:,:,2)]);
fooHbO_dmndan = fooHbO_dmndan(1:15,16:end);
fooHbR_dmndan = corrcoef([BrainMaps_hbr(:,:,1),BrainMaps_hbr(:,:,2)]);
fooHbR_dmndan = fooHbR_dmndan(1:15,16:end);

zHbO_dmn = 0.5 * log((1+fooHbO_dmn)./(1-fooHbO_dmn));
zHbR_dmn = 0.5 * log((1+fooHbR_dmn)./(1-fooHbR_dmn));
zHbO_dan = 0.5 * log((1+fooHbO_dan)./(1-fooHbO_dan));
zHbR_dan = 0.5 * log((1+fooHbR_dan)./(1-fooHbR_dan));
zHbO_dmndan = 0.5 * log((1+fooHbO_dmndan)./(1-fooHbO_dmndan));
zHbR_dmndan = 0.5 * log((1+fooHbR_dmndan)./(1-fooHbR_dmndan));

f=figure();
subplot(1,2,1);
lst_ = find(tril(ones(size(zHbO_dmn)),-1));
boxplot([zHbO_dmn(lst_),zHbO_dan(lst_),zHbO_dmndan(lst_)]);
[h_dmn,p_dmn,~,stats_dmn] = ttest2(zHbO_dmn(lst_),zHbO_dmndan(lst_));
[h_dan,p_dan,~,stats_dan] = ttest2(zHbO_dan(lst_),zHbO_dmndan(lst_));
fprintf('HbO - DMN tval:%.2g, dof:%i, pval:%.3g\n',stats_dmn.tstat,stats_dmn.df,p_dmn);
fprintf('HbO - DAN tval:%.2g, dof:%i, pval:%.3g\n',stats_dan.tstat,stats_dan.df,p_dan);
if p_dmn< 0.01
    line([1,3],[0.97,0.97],'Color','k','LineWidth',2.5);
end
if p_dan< 0.01
    line([2,3],[1.1,1.1],'Color','k','LineWidth',2.5);
end
ylim([-1 1.2]);
xticklabels({'DMN-DMN','DAN-DAN','DMN-DAN'});
ylabel('Fisher Z transform');
title('HbO');
subplot(1,2,2);
boxplot([zHbR_dmn(lst_),zHbR_dan(lst_),zHbR_dmndan(lst_)]);
[h_dmn,p_dmn,~,stats_dmn] = ttest2(zHbR_dmn(lst_),zHbR_dmndan(lst_));
[h_dan,p_dan,~,stats_dan] = ttest2(zHbR_dan(lst_),zHbR_dmndan(lst_));
fprintf('HbR - DMN tval:%.2g, dof:%i, pval:%.3g\n',stats_dmn.tstat,stats_dmn.df,p_dmn);
fprintf('HbR - DAN tval:%.2g, dof:%i, pval:%.3g\n',stats_dan.tstat,stats_dan.df,p_dan);
if p_dmn< 0.01
    line([1,3],[0.97,0.97],'Color','k','LineWidth',2.5);
end
if p_dan< 0.01
    line([2,3],[1.1,1.1],'Color','k','LineWidth',2.5);
end
ylim([-1 1.2]);
xticklabels({'DMN-DMN','DAN-DAN','DMN-DAN'});
title('HbR');

f.Position = [10         10        1049         507];
%saveas(f,['Group-starClusters-DMN-DAN_DavidClustering.png']);
saveas(f,[fnameoutput,'_boxplots.png']);