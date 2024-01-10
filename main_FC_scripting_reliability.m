%%
% test-retest reliability using individual imprpved seeds/submasks
clear all;

[fwFolder,anatomFolder,derivFolder,dataDir] = setmyenv();

subjects_set = [1:7 9:16];
nSubjs=length(subjects_set);
BrainMaps_hbo =[];
BrainMaps_hbr = [];
dmn_improv_hbo = cell(nSubjs,1);
dan_improv_hbo = cell(nSubjs,1);
dmn_improv_hbr = cell(nSubjs,1);
dan_improv_hbr = cell(nSubjs,1);
rDMNDAN_AllSubj_hbo = zeros(40,40,nSubjs);
rDMNDAN_AllSubj_hbr = zeros(40,40,nSubjs);
[dmn_mask,dan_mask,mesh_brain,idx_select] = Parcellation_FC(anatomFolder,fwFolder);
dmn_seeds_ts = cell(length(dmn_mask),2,nSubjs);
dan_seeds_ts = cell(length(dan_mask),2,nSubjs);


% flags and thresholds
flags.macorrect = 'spline'; % 'none' or 'spline'
flags.bpfilt = 'image';% 'none' 'channel' or 'image'
flags.imagerecon = 'brain'; %'brain' or 'brain+scalp'
flags.rhoSD_ssThresh = 15;
flags.gsr = 'image';%'none','channel' or 'image'
flags.r_thresh = 0.7; % .r_thresh is the threshold for the clustering
flags.plot=0; % .plot  flag to plot the brain correlation map
flags.p_thresh = 0; % . p_thresh is used to plot r values below that p-val (use 0 to plot all the correlations)
flags.clusteringType = 1; %1:Matlab, 2:David's algorithm
pipeline_str = sprintf('macor-%s_bpfilt-%s_imrec-%s_gsr-%s_clust-%i',...
    flags.macorrect,flags.bpfilt,flags.imagerecon,flags.gsr,flags.clusteringType);
fOut_reliability=sprintf('reliability_%s',pipeline_str);
fOut_pmap=sprintf('probMap_%s',pipeline_str);
pipelineDir = sprintf('%sPipeline-%s/',derivFolder,pipeline_str);
if ~exist(pipelineDir,'dir')
    mkdir(pipelineDir);
end

%%
if ~exist([pipelineDir,fOut_reliability,'.mat'],'file');
    for iSubj = 1:nSubjs
        subject = num2str(subjects_set(iSubj));
        fprintf('==============================\n');
        fprintf('Subject %s\n',subject);
        fprintf('==============================\n');
        SnirfFilePathr1 = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-RS_run-1_nirs.snirf'];
        SnirfFilePathr2 = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-RS_run-2_nirs.snirf'];

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
        %conctenation of non-bp-filtered data
        [HbO_brain_r1r2NoFilt] = [HbO_brain_chunkr1;HbO_brain_chunkr2];
        [HbR_brain_r1r2NoFilt] = [HbR_brain_chunkr1;HbR_brain_chunkr2];

        %%
        % Band pass filtering in image space
        if strcmp(flags.bpfilt,'image')
            %[y2] = image_BandpassFilt(y,hpf,lpf,fs)
            fs = mean(1./diff(snirfObjr1.data.time));
            HbO_brain_chunkr1 = image_BandpassFilt(HbO_brain_chunkr1, 0.009, 0.080,fs);
            HbR_brain_chunkr1 = image_BandpassFilt(HbR_brain_chunkr1, 0.009, 0.080,fs);
            HbO_brain_chunkr2 = image_BandpassFilt(HbO_brain_chunkr2, 0.009, 0.080,fs);
            HbR_brain_chunkr2 = image_BandpassFilt(HbR_brain_chunkr2, 0.009, 0.080,fs);
        end
        % Global signal regression in image space?
        if strcmp(flags.gsr,'image') && ~strcmp(flags.bpfilt,'none')
            HbO_brain_chunkr1 = GlobalRegression(HbO_brain_chunkr1);
            HbR_brain_chunkr1 = GlobalRegression(HbR_brain_chunkr1);
            HbO_brain_chunkr2 = GlobalRegression(HbO_brain_chunkr2);
            HbR_brain_chunkr2 = GlobalRegression(HbR_brain_chunkr2);
        end
        %%
        % Do we want to concatenate r1 and r2?
        [HbO_brain_r1r2] = [HbO_brain_chunkr1;HbO_brain_chunkr2];
        [HbR_brain_r1r2] = [HbR_brain_chunkr1;HbR_brain_chunkr2];
        %%
        [dmn_mask_hbo] = Clustering_FC(HbO_brain_r1r2,dmn_mask,flags);
        [dmn_mask_hbr] = Clustering_FC(HbR_brain_r1r2,dmn_mask,flags);

        %%
        [dan_mask_hbo] = Clustering_FC(HbO_brain_r1r2,dan_mask,flags);
        [dan_mask_hbr] = Clustering_FC(HbR_brain_r1r2,dan_mask,flags);

        %%
        [~,dmn_improv_hbo{iSubj},dan_improv_hbo{iSubj},~] = ClusterSelection_FC(dmn_mask_hbo,dan_mask_hbo,HbO_brain_r1r2,flags);
        [~,dmn_improv_hbr{iSubj},dan_improv_hbr{iSubj},~] = ClusterSelection_FC(dmn_mask_hbr,dan_mask_hbr,HbR_brain_r1r2,flags);

        %%
        % HbO
        fprintf('(hbo)DMN submask:');
        BrainMaps_hbo = zeros(length(idx_select),2*(length(dmn_improv_hbo{iSubj})+length(dan_improv_hbo{iSubj})));
        for iSubmask=1:length(dmn_improv_hbo{iSubj})
            [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbo{iSubj}(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);
            [~,hmap2,A_select2] = CorrelationBrainMap_FC(dmn_improv_hbo{iSubj}(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr2,flags.p_thresh,flags.plot);
            BrainMaps_hbo(:,iSubmask) = A_select1;
            BrainMaps_hbo(:,length(dmn_improv_hbo{iSubj})+length(dan_improv_hbo{iSubj})+iSubmask) = A_select2;
            dmn_seeds_ts{iSubmask,1,iSubj} = mean(HbO_brain_r1r2NoFilt(:,dmn_improv_hbo{iSubj}(iSubmask).vertices_index(dmn_improv_hbo{iSubj}(iSubmask).mask_subsetseed)),2);
            fprintf('%i\t',iSubmask);
        end
        fprintf('\n(hbo)DAN submask:');
        for iSubmask=1:length(dan_improv_hbo{iSubj})
            [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbo{iSubj}(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);
            [~,hmap2,A_select2] = CorrelationBrainMap_FC(dan_improv_hbo{iSubj}(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr2,flags.p_thresh,flags.plot);
            BrainMaps_hbo(:,length(dmn_improv_hbo{iSubj})+iSubmask) = A_select1;
            BrainMaps_hbo(:,2*length(dmn_improv_hbo{iSubj})+length(dan_improv_hbo{iSubj})+iSubmask) = A_select2;
            dan_seeds_ts{iSubmask,1,iSubj} = mean(HbO_brain_r1r2NoFilt(:,dan_improv_hbo{iSubj}(iSubmask).vertices_index(dan_improv_hbo{iSubj}(iSubmask).mask_subsetseed)),2);
            fprintf('%i\t',iSubmask);
        end

        % hbr
        BrainMaps_hbr = zeros(length(idx_select),2*(length(dmn_improv_hbr{iSubj})+length(dan_improv_hbr{iSubj})));
        fprintf('\n(hbr)DMN submask:');
        for iSubmask=1:length(dmn_improv_hbr{iSubj})
            [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbr{iSubj}(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
            [~,hmap2,A_select2] = CorrelationBrainMap_FC(dmn_improv_hbr{iSubj}(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr2,flags.p_thresh,flags.plot);
            BrainMaps_hbr(:,iSubmask) = A_select1;
            BrainMaps_hbr(:,length(dmn_improv_hbr{iSubj})+length(dan_improv_hbr{iSubj})+iSubmask) = A_select2;
            dmn_seeds_ts{iSubmask,2,iSubj} = mean(HbR_brain_r1r2NoFilt(:,dmn_improv_hbr{iSubj}(iSubmask).vertices_index(dmn_improv_hbr{iSubj}(iSubmask).mask_subsetseed)),2);
            fprintf('%i\t',iSubmask);
        end

        fprintf('\n(hbr)DAN submask:');
        for iSubmask=1:length(dan_improv_hbr{iSubj})
            [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbr{iSubj}(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
            [~,hmap2,A_select2] = CorrelationBrainMap_FC(dan_improv_hbr{iSubj}(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr2,flags.p_thresh,flags.plot);
            BrainMaps_hbr(:,length(dmn_improv_hbr{iSubj})+iSubmask) = A_select1;
            BrainMaps_hbr(:,2*length(dmn_improv_hbr{iSubj})+length(dan_improv_hbr{iSubj})+iSubmask) = A_select2;
            dan_seeds_ts{iSubmask,2,iSubj} = mean(HbR_brain_r1r2NoFilt(:,dan_improv_hbr{iSubj}(iSubmask).vertices_index(dan_improv_hbr{iSubj}(iSubmask).mask_subsetseed)),2);
            fprintf('%i\t',iSubmask);
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
    save([pipelineDir,fOut_reliability,'.mat'],'rDMNDAN_AllSubj_hbo','rDMNDAN_AllSubj_hbr','flags',...
        'dmn_seeds_ts','dan_seeds_ts','dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr');
else
    load([pipelineDir,fOut_reliability,'.mat'],'rDMNDAN_AllSubj_hbo','rDMNDAN_AllSubj_hbr','flags', ...
        'dmn_seeds_ts','dan_seeds_ts','dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr');
end
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
saveas(f,[pipelineDir,fOut_reliability,'.png']);
close(f);

%%
[hmap] = probabilityMap_FC(dmn_improv_hbo,mesh_brain,idx_select,'Probability Group DMN HbO Map');
saveas(hmap,[pipelineDir,fOut_pmap,'_dmn_hbo.png']);
close(hmap);
[hmap] = probabilityMap_FC(dan_improv_hbo,mesh_brain,idx_select,'Probability Group DAN HbO Map');
saveas(hmap,[pipelineDir,fOut_pmap,'_dan_hbo.png']);
close(hmap);
[hmap] = probabilityMap_FC(dmn_improv_hbr,mesh_brain,idx_select,'Probability Group DMN HbR Map');
saveas(hmap,[pipelineDir,fOut_pmap,'_dmn_hbr.png']);
close(hmap);
[hmap] = probabilityMap_FC(dan_improv_hbr,mesh_brain,idx_select,'Probability Group DAN HbR Map');
saveas(hmap,[pipelineDir,fOut_pmap,'_dan_hbr.png']);
close(hmap);