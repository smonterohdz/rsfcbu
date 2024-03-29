%%
% test-retest reliability using a SUPER SEED
fwFolder = 'C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\fw\';
anatomFolder = 'C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\probe_10MPhotons\anatomical\';

subjects_set = [1:7 9:16];
BrainMaps_hbo =[];
BrainMaps_hbr = [];
HbO_brain_r1r2_all = [];
HbR_brain_r1r2_all = [];
rDMNDAN_AllSubj_hbo = zeros(40,40,length(subjects_set));
rDMNDAN_AllSubj_hbr = zeros(40,40,length(subjects_set));
nSubjs=length(subjects_set);

% flags and thresholds
flags.imagerecon = 'brain+scalp'; %'brain' or 'brain+scalp'
flags.rhoSD_ssThresh = 15;
flags.gsr = 'channel';%'none','channel' or 'image'
flags.r_thresh = 0.7; % .r_thresh is the threshold for the clustering
flags.plot=0; % .plot  flag to plot the brain correlation map
flags.p_thresh = 0; % . p_thresh is used to plot r values below that p-val (use 0 to plot all the correlations)
flags.clusteringType = 1; %1:Matlab, 2:David's algorithm
fnameoutput=sprintf('rDMNDAN_AllSubj_imrecon-%s_gsr-%s_clust-%i',flags.imagerecon,flags.gsr,flags.clusteringType);
[dmn_mask,dan_mask,mesh_brain,idx_select] = Parcellation_FC(anatomFolder,fwFolder);

for iSubj = 1:nSubjs    
    subject = num2str(subjects_set(iSubj));
    fprintf('Subject %s------------------------\n',subject);
    SnirfFilePathr1 = ['C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\sub-',subject,'\nirs\sub-',subject,'_task-RS_run-1_nirs.snirf'];
    SnirfFilePathr2 = ['C:\Users\smontero\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\sub-',subject,'\nirs\sub-',subject,'_task-RS_run-2_nirs.snirf'];
 
    %%
    [snirfObjr1,dcObjr1,dodObjr1] = Preprocessing_FC(SnirfFilePathr1,flags);
    [snirfObjr2,dcObjr2,dodObjr2] = Preprocessing_FC(SnirfFilePathr2,flags);

    %%
    [HbO_brainr1,HbR_brainr1] = ImageReconstruction_FC(snirfObjr1,dodObjr1,dcObjr1,fwFolder,flags);
    [HbO_brainr2,HbR_brainr2] = ImageReconstruction_FC(snirfObjr2,dodObjr2,dcObjr2,fwFolder,flags);

    %%
    % use twindow.init_sec = -1 if you want to use the onset and duration
    % defined in the snirf file (assuming there is only one stimulus).
    % Otherwise, change the values (in sec.) according to your needs
    twindow.init_sec = -1;
    twindow.offset_sec = 60;
    [HbO_brain_chunkr1] = ExctractChunk(HbO_brainr1,snirfObjr1,twindow);
    [HbO_brain_chunkr2] = ExctractChunk(HbO_brainr2,snirfObjr2,twindow);
    [HbR_brain_chunkr1] = ExctractChunk(HbR_brainr1,snirfObjr1,twindow);
    [HbR_brain_chunkr2] = ExctractChunk(HbR_brainr2,snirfObjr2,twindow);
    %%
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

%%
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
    %%
    [HbO_brainr1,HbR_brainr1] = ImageReconstruction_FC(snirfObjr1,dodObjr1,dcObjr1,fwFolder,flags);
    [HbO_brainr2,HbR_brainr2] = ImageReconstruction_FC(snirfObjr2,dodObjr2,dcObjr2,fwFolder,flags);

    %%
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
    %%
    % Do we want to concatenate r1 and r2?
    % [HbO_brain_r1r2] = [HbO_brain_chunkr1;HbO_brain_chunkr2];
    % [HbR_brain_r1r2] = [HbR_brain_chunkr1;HbR_brain_chunkr2];
    % %%
    % [dmn_mask_hbo] = Clustering_FC(HbO_brain_r1r2,dmn_mask,flags);
    % [dmn_mask_hbr] = Clustering_FC(HbR_brain_r1r2,dmn_mask,flags);
    % 
    % %%
    % [dan_mask_hbo] = Clustering_FC(HbO_brain_r1r2,dan_mask,flags);
    % [dan_mask_hbr] = Clustering_FC(HbR_brain_r1r2,dan_mask,flags);

    % %%
    % [~,dmn_improv_hbo,dan_improv_hbo,~] = ClusterSelection_FC(dmn_mask_hbo,dan_mask_hbo,HbO_brain_r1r2,flags);
    % [~,dmn_improv_hbr,dan_improv_hbr,~] = ClusterSelection_FC(dmn_mask_hbr,dan_mask_hbr,HbR_brain_r1r2,flags);

    %%
    % HbO
    BrainMaps_hbo = zeros(length(idx_select),2*(length(dmn_improv_hbo)+length(dan_improv_hbo)));
    for iSubmask=1:length(dmn_improv_hbo)
        [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);
        [~,hmap2,A_select2] = CorrelationBrainMap_FC(dmn_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr2,flags.p_thresh,flags.plot);
        BrainMaps_hbo(:,iSubmask) = A_select1;
        BrainMaps_hbo(:,length(dmn_improv_hbo)+length(dan_improv_hbo)+iSubmask) = A_select2;
        fprintf('(hbo)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbo));
    end
    for iSubmask=1:length(dan_improv_hbo)
        [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);
        [~,hmap2,A_select2] = CorrelationBrainMap_FC(dan_improv_hbo(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr2,flags.p_thresh,flags.plot);
        BrainMaps_hbo(:,length(dmn_improv_hbo)+iSubmask) = A_select1;
        BrainMaps_hbo(:,2*length(dmn_improv_hbo)+length(dan_improv_hbo)+iSubmask) = A_select2;
        fprintf('(hbo)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbo));
    end
    
    % hbr
    BrainMaps_hbr = zeros(length(idx_select),2*(length(dmn_improv_hbr)+length(dan_improv_hbr)));
    for iSubmask=1:length(dmn_improv_hbr)
        [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
        [~,hmap2,A_select2] = CorrelationBrainMap_FC(dmn_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr2,flags.p_thresh,flags.plot);
        BrainMaps_hbr(:,iSubmask) = A_select1;
        BrainMaps_hbr(:,length(dmn_improv_hbr)+length(dan_improv_hbr)+iSubmask) = A_select2;
        fprintf('(hbr)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbr));
    end
    for iSubmask=1:length(dan_improv_hbr)
        [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
        [~,hmap2,A_select2] = CorrelationBrainMap_FC(dan_improv_hbr(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr2,flags.p_thresh,flags.plot);
        BrainMaps_hbr(:,length(dmn_improv_hbr)+iSubmask) = A_select1;
        BrainMaps_hbr(:,2*length(dmn_improv_hbr)+length(dan_improv_hbr)+iSubmask) = A_select2;
        fprintf('(hbr)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbr));
    end
    %%
    %test-retest subj
    [r,p]=corrcoef(BrainMaps_hbo);
    %figure(10), imagesc(r,[-1 1]); colormap("jet")
    rDMNDAN_AllSubj_hbo(:,:,iSubj) = r;

    [r,p]=corrcoef(BrainMaps_hbr);
    %figure(10), imagesc(r,[-1 1]); colormap("jet")
    rDMNDAN_AllSubj_hbr(:,:,iSubj) = r;
end

save([fnameoutput,'.mat'],'rDMNDAN_AllSubj_hbo','rDMNDAN_AllSubj_hbr','flags');
%%-
%test-retest group
fooHbO = rDMNDAN_AllSubj_hbo(1:20,21:end,:);
fooHbR = rDMNDAN_AllSubj_hbr(1:20,21:end,:);

zHbO = 0.5 * log((1+fooHbO)./(1-fooHbO));
zHbR = 0.5 * log((1+fooHbR)./(1-fooHbR));

f=figure();
colormap(jet(64))

subplot(2,3,1)
imagesc( mean(fooHbO,3), [-1 1] );
ylabel({'HbO';'Submasks Run1'});
title('mean R')
colorbar

subplot(2,3,4)
imagesc( mean(fooHbR,3), [-1 1] );
ylabel({'HbR';'Submasks Run1'});
colorbar
xlabel({'Submasks Run2'});



subplot(2,3,2)
foo = mean(zHbO,3);
%lst = find(eye(size(foo))==1); foo(lst) = 0; % remove diagonal
imagesc( foo, [-1 1]*max(abs(foo(:))) )
title('mean Z')
colorbar

subplot(2,3,5)
foo = mean(zHbR,3);
%lst = find(eye(size(foo))==1); foo(lst) = 0; % remove diagonal
imagesc( foo, [-1 1]*max(abs(foo(:))) )
xlabel({'Submasks Run2'});
colorbar



subplot(2,3,3)
fooM = mean(zHbO,3);
fooSE = std(zHbO,[],3) / sqrt(nSubjs);
foo = fooM./fooSE;
fooSign = sign(foo);
foo = -log10(2*(1-tcdf(abs(fooM./fooSE),nSubjs-1)));
foo = foo .* (foo>2) .* fooSign;
%lst = find(eye(size(foo))==1); foo(lst) = 0; % remove diagonal
imagesc( foo, [-1 1]*max(abs(foo(:))) )
title('-log10( p_{val}(Z T-Statistic) )')
colorbar

subplot(2,3,6)
fooM = mean(zHbR,3);
fooSE = std(zHbR,[],3) / sqrt(nSubjs);
foo = fooM./fooSE;
fooSign = sign(foo);
foo = -log10(2*(1-tcdf(abs(fooM./fooSE),nSubjs-1)));
foo = foo .* (foo>2) .* fooSign;
%lst = find(eye(size(foo))==1); foo(lst) = 0; % remove diagonal
imagesc( foo, [-1 1]*max(abs(foo(:))) )
colorbar
xlabel({'Submasks Run2'});

f.Position = [10         10        1049         507];
saveas(f,[fnameoutput,'.png']);