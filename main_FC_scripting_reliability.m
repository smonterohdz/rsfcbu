%%
% test-retest reliability using individual imprpved seeds/submasks
clear all;
OVERWRITE_ = 1;

subjects_set = [1:7 9:16];
nSubjs=length(subjects_set);
BrainMaps_hbo =[];
BrainMaps_hbr = [];
dmn_improv_hbo = cell(nSubjs,1);
dan_improv_hbo = cell(nSubjs,1);
dmn_improv_hbr = cell(nSubjs,1);
dan_improv_hbr = cell(nSubjs,1);

% flags and thresholds
flags.macorrect = 'splineSG'; % 'none' or 'spline'
flags.bpfilt = 'image';% 'none' 'channel' or 'image'
flags.imagerecon = 'brain+scalp'; %'brain' or 'brain+scalp'
flags.rhoSD_ssThresh = 15;
flags.gsr = 'image';%'none','channel' or 'image'
flags.r_thresh = 0.7; % .r_thresh is the threshold for the clustering
flags.plot = 0; % .plot  flag to plot the brain correlation map
flags.p_thresh = 0; % . p_thresh is used to plot r values below that p-val (use 0 to plot all the correlations)
flags.clusteringType = 1; %0: no clustering, 1:Matlab, 2:David's algorithm
flags.task = 'RS';
flags.parcel_scheme = 'schaefer_comb';
flags.Adot_scaling_factor = 100;
[fwFolder,anatomFolder,derivFolder,dataDir] = setmyenv(flags);

[dmn_mask,dan_mask,mesh_brain,idx_select,dmn_z,dan_z] = Parcellation_test_FC(anatomFolder,fwFolder,flags);
dmn_mask(dmn_z) = [];
dan_mask(dan_z) = [];
nparcelsdmn = length(dmn_mask);
nparcelsdan = length(dan_mask);

%dmn_mask(1).type = "seed";
%f1=plot_net_mask(mesh_brain,idx_select,dmn_mask);
%dan_mask(1).type="seed";
%f2=plot_net_mask(mesh_brain,idx_select,dan_mask);
%dmn_mask(1).type = dmn_mask(2).type;
%dan_mask(1).type = dan_mask(2).type;

dmn_seeds_ts = cell(length(dmn_mask),2,nSubjs);
dan_seeds_ts = cell(length(dan_mask),2,nSubjs);
dmn_seeds_NP_ts = cell(length(dmn_mask),2,nSubjs);
dan_seeds_NP_ts = cell(length(dan_mask),2,nSubjs);

pipeline_str = sprintf('%s-macor-%s_bpfilt-%s_imrec-%s_gsr-%s_clust-%i_parcel-%s_Asf-%i',...
    flags.task,flags.macorrect,flags.bpfilt,flags.imagerecon,flags.gsr, ...
    flags.clusteringType,flags.parcel_scheme,flags.Adot_scaling_factor);
fOut_reliability=sprintf('reliability_%s',pipeline_str);
fOut_pmap=sprintf('probMap_%s',pipeline_str);
pipelineDir = sprintf('%sPipeline-%s/',derivFolder,pipeline_str);
if ~exist(pipelineDir,'dir')
    mkdir(pipelineDir);
end
% saveas(f1,[pipelineDir,'DMN_regions_3D.png']);
% close(f1);
% saveas(f2,[pipelineDir,'DAN_regions_3D.png']);
% close(f2);

rDMNDAN_AllSubj_hbo = zeros((nparcelsdan+nparcelsdmn)*2,(nparcelsdmn+nparcelsdan)*2,nSubjs);
rDMNDAN_AllSubj_hbr = zeros((nparcelsdan+nparcelsdmn)*2,(nparcelsdmn+nparcelsdan)*2,nSubjs);
%%
if OVERWRITE_ || ~exist([pipelineDir,fOut_reliability,'.mat'],'file')
    BrainMaps_hbo = zeros(length(idx_select),2*(length(dmn_mask)+length(dan_mask)),nSubjs);
    for iSubj = 1:nSubjs
        subject = num2str(subjects_set(iSubj));
        fprintf('==============================\n');
        fprintf('Subject %s\n',subject);
        fprintf('==============================\n');
        SnirfFilePathr1 = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-',flags.task,'_run-1_nirs.snirf'];
        SnirfFilePathr2 = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-',flags.task,'_run-2_nirs.snirf'];

        %%-
        [snirfObjr1,dcObjr1,dodObjr1,dcr1NP,dodr1NP] = Preprocessing_FC(SnirfFilePathr1,flags);
        [snirfObjr2,dcObjr2,dodObjr2,dcr2NP,dodr2NP] = Preprocessing_FC(SnirfFilePathr2,flags);
        %%
        [HbO_brainr1,HbR_brainr1] = ImageReconstruction_FC(snirfObjr1,dodObjr1,dcObjr1,fwFolder,flags);
        [HbO_brainr2,HbR_brainr2] = ImageReconstruction_FC(snirfObjr2,dodObjr2,dcObjr2,fwFolder,flags);
        %fs = mean(1./diff(snirfObjr1.data.time));
        %t.start = 180; t.end = 360;
        %image_recon_plot_movie([pipelineDir,'HbO_vid.mp4'],mesh_brain,idx_select,HbO_brainr1,[],t,fs);       

        %no processing
        %[HbO_r1NP,HbR_r1NP] = ImageReconstruction_FC(snirfObjr1,dodr1NP,dcr1NP,fwFolder,flags);
        %[HbO_r2NP,HbR_r2NP] = ImageReconstruction_FC(snirfObjr2,dodr2NP,dcr2NP,fwFolder,flags);
        % Do I want to visualize HbO from the Image recon?
        %checkImg_FC(fwFolder,HbO_brainr1,dodObjr1.time);
        %%
        % use twindow.init_sec = -1 if you want to use the onset and duration
        % defined in the snirf file (assuming there is only one stimulus).
        % Otherwise, change the values (in sec.) according to your needs
        twindow.init_sec = -1;
        twindow.offset_sec = 60;
        %twindow.init_sec = 30;
        %twindow.dur_sec = 180;
        [HbO_brain_chunkr1] = ExctractChunk(HbO_brainr1,snirfObjr1,twindow,flags);
        [HbO_brain_chunkr2] = ExctractChunk(HbO_brainr2,snirfObjr2,twindow,flags);
        [HbR_brain_chunkr1] = ExctractChunk(HbR_brainr1,snirfObjr1,twindow,flags);
        [HbR_brain_chunkr2] = ExctractChunk(HbR_brainr2,snirfObjr2,twindow,flags);
        % no processing
        % [HbO_chunkr1NP] = ExctractChunk(HbO_r1NP,snirfObjr1,twindow,flags);
        % [HbO_chunkr2NP] = ExctractChunk(HbO_r2NP,snirfObjr2,twindow,flags);
        % [HbR_chunkr1NP] = ExctractChunk(HbR_r1NP,snirfObjr1,twindow,flags);
        % [HbR_chunkr2NP] = ExctractChunk(HbR_r2NP,snirfObjr2,twindow,flags);
        %% 
        %conctenation of non-bp-filtered data
        % [HbO_r1r2NP] = [HbO_chunkr1NP;HbO_chunkr2NP];
        % [HbR_r1r2NP] = [HbR_chunkr1NP;HbR_chunkr2NP];
        % no processing
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
        if flags.clusteringType > 0
            [dmn_mask_hbo] = Clustering_FC(HbO_brain_r1r2,dmn_mask,flags);
            [dmn_mask_hbr] = Clustering_FC(HbR_brain_r1r2,dmn_mask,flags);
        else
            [dmn_improv_hbo{iSubj}] = dmn_mask;
            [dmn_improv_hbr{iSubj}] = dmn_mask;
        end

        %%
        if flags.clusteringType > 0
            [dan_mask_hbo] = Clustering_FC(HbO_brain_r1r2,dan_mask,flags);
            [dan_mask_hbr] = Clustering_FC(HbR_brain_r1r2,dan_mask,flags);
        else
            [dan_improv_hbo{iSubj}] = dan_mask;
            [dan_improv_hbr{iSubj}] = dan_mask;
        end
        %%
        if flags.clusteringType > 0
            [~,dmn_improv_hbo{iSubj},dan_improv_hbo{iSubj},~] = ClusterSelection_FC(dmn_mask_hbo,dan_mask_hbo,HbO_brain_r1r2,flags);
            [~,dmn_improv_hbr{iSubj},dan_improv_hbr{iSubj},~] = ClusterSelection_FC(dmn_mask_hbr,dan_mask_hbr,HbR_brain_r1r2,flags);
        end        

        %%
        % HbO
        fprintf('\n(hbo)DMN submask:');
        %BrainMaps_hbo = zeros(length(idx_select),2*(length(dmn_improv_hbo{iSubj})+length(dan_improv_hbo{iSubj})));
        for iSubmask=1:length(dmn_improv_hbo{iSubj})
            [~,hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbo{iSubj}(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);
            [~,hmap2,A_select2] = CorrelationBrainMap_FC(dmn_improv_hbo{iSubj}(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr2,flags.p_thresh,flags.plot);
            BrainMaps_hbo(:,iSubmask,iSubj) = A_select1;
            BrainMaps_hbo(:,length(dmn_improv_hbo{iSubj})+length(dan_improv_hbo{iSubj})+iSubmask,iSubj) = A_select2;
            %dmn_seeds_NP_ts{iSubmask,1,iSubj} = mean(HbO_r1r2NP(:,dmn_improv_hbo{iSubj}(iSubmask).vertices_index(dmn_improv_hbo{iSubj}(iSubmask).mask_subsetseed)),2);
            dmn_seeds_ts{iSubmask,1,iSubj} = mean(HbO_brain_r1r2(:,dmn_improv_hbo{iSubj}(iSubmask).vertices_index(dmn_improv_hbo{iSubj}(iSubmask).mask_subsetseed)),2);
            fprintf('%i\t',iSubmask);
        end
        fprintf('\n(hbo)DAN submask:');
        for iSubmask=1:length(dan_improv_hbo{iSubj})
            [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbo{iSubj}(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr1,flags.p_thresh,flags.plot);
            [~,hmap2,A_select2] = CorrelationBrainMap_FC(dan_improv_hbo{iSubj}(iSubmask),mesh_brain,idx_select,HbO_brain_chunkr2,flags.p_thresh,flags.plot);
            BrainMaps_hbo(:,length(dmn_improv_hbo{iSubj})+iSubmask,iSubj) = A_select1;
            BrainMaps_hbo(:,2*length(dmn_improv_hbo{iSubj})+length(dan_improv_hbo{iSubj})+iSubmask,iSubj) = A_select2;
            %dan_seeds_NP_ts{iSubmask,1,iSubj} = mean(HbO_r1r2NP(:,dan_improv_hbo{iSubj}(iSubmask).vertices_index(dan_improv_hbo{iSubj}(iSubmask).mask_subsetseed)),2);
            dan_seeds_ts{iSubmask,1,iSubj} = mean(HbO_brain_r1r2(:,dan_improv_hbo{iSubj}(iSubmask).vertices_index(dan_improv_hbo{iSubj}(iSubmask).mask_subsetseed)),2);
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
            %dmn_seeds_NP_ts{iSubmask,2,iSubj} = mean(HbR_r1r2NP(:,dmn_improv_hbr{iSubj}(iSubmask).vertices_index(dmn_improv_hbr{iSubj}(iSubmask).mask_subsetseed)),2);
            dmn_seeds_ts{iSubmask,2,iSubj} = mean(HbR_brain_r1r2(:,dmn_improv_hbr{iSubj}(iSubmask).vertices_index(dmn_improv_hbr{iSubj}(iSubmask).mask_subsetseed)),2);
            fprintf('%i\t',iSubmask);
        end

        fprintf('\n(hbr)DAN submask:');
        for iSubmask=1:length(dan_improv_hbr{iSubj})
            [~,hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbr{iSubj}(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr1,flags.p_thresh,flags.plot);
            [~,hmap2,A_select2] = CorrelationBrainMap_FC(dan_improv_hbr{iSubj}(iSubmask),mesh_brain,idx_select,HbR_brain_chunkr2,flags.p_thresh,flags.plot);
            BrainMaps_hbr(:,length(dmn_improv_hbr{iSubj})+iSubmask) = A_select1;
            BrainMaps_hbr(:,2*length(dmn_improv_hbr{iSubj})+length(dan_improv_hbr{iSubj})+iSubmask) = A_select2;
            %dan_seeds_NP_ts{iSubmask,2,iSubj} = mean(HbR_r1r2NP(:,dan_improv_hbr{iSubj}(iSubmask).vertices_index(dan_improv_hbr{iSubj}(iSubmask).mask_subsetseed)),2);
            dan_seeds_ts{iSubmask,2,iSubj} = mean(HbR_brain_r1r2(:,dan_improv_hbr{iSubj}(iSubmask).vertices_index(dan_improv_hbr{iSubj}(iSubmask).mask_subsetseed)),2);
            fprintf('%i\t',iSubmask);
        end
        %%
        %test-retest subj
        [r,p]=corrcoef(BrainMaps_hbo(:,:,iSubj));
        %figure(10), imagesc(r,[-1 1]); colormap("jet")
        rDMNDAN_AllSubj_hbo(:,:,iSubj) = r;

        [r,p]=corrcoef(BrainMaps_hbr);
        %figure(10), imagesc(r,[-1 1]); colormap("jet")
        rDMNDAN_AllSubj_hbr(:,:,iSubj) = r;
    end
    A_select_mean = mean(BrainMaps_hbo(:,15,:),3);
    A = zeros(1,size(mesh_brain.vertices,1));
    A(idx_select) = A_select_mean;
    hmap = plot_connectivity_seedregion(mesh_brain,A,[],idx_select(dmn_mask(5).vertices_index), ...
        [],[],[],[]);
    hold on;
    save([pipelineDir,fOut_reliability,'.mat'],'rDMNDAN_AllSubj_hbo','rDMNDAN_AllSubj_hbr','flags',...
        'dmn_seeds_NP_ts','dan_seeds_NP_ts','dmn_seeds_ts','dan_seeds_ts',...
        'dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr','BrainMaps_hbo');
else
    load([pipelineDir,fOut_reliability,'.mat'],'rDMNDAN_AllSubj_hbo','rDMNDAN_AllSubj_hbr','flags', ...
        'dmn_seeds_NP_ts','dan_seeds_NP_ts','dmn_seeds_ts','dan_seeds_ts',...
        'dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr','BrainMaps_hbo');
end
%%
%test-retest group
fooHbO = rDMNDAN_AllSubj_hbo(1:(nparcelsdmn+nparcelsdan),(nparcelsdmn+nparcelsdan+1):end,:);
fooHbR = rDMNDAN_AllSubj_hbr(1:(nparcelsdmn+nparcelsdan),(nparcelsdmn+nparcelsdan+1):end,:);

[f1,f2]= plot_corrMat_FC(rDMNDAN_AllSubj_hbo,rDMNDAN_AllSubj_hbr,subjects_set);
saveas(f1,[pipelineDir,'HbO_Subjects_BigCorrMat.png']);
saveas(f2,[pipelineDir,'HbR_Subjects_BigCorrMat.png']);
close(f1);
close(f2);


[f1,f2]= plot_corrMat_FC(fooHbO,fooHbR,subjects_set);
saveas(f1,[pipelineDir,'HbO_Subjects_CorrMat.png']);
saveas(f2,[pipelineDir,'HbR_Subjects_CorrMat.png']);
close(f1);
close(f2);


zHbO = 0.5 * log((1+fooHbO)./(1-fooHbO));
zHbR = 0.5 * log((1+fooHbR)./(1-fooHbR));

f=figure(); 
fsize = 7;
colormap(jet(64))

subplot(2,3,1)
imagesc( mean(fooHbO,3), [-1 1] );
ax= gca();
ax.TickLabelInterpreter='none';
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
ylabel({'HbO';'Regions Run1'});
yticks(1:nparcelsdmn+nparcelsdan);
yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xticklabels({});
title('mean R')
colorbar

subplot(2,3,4)
imagesc( mean(fooHbR,3), [-1 1] );
ax= gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
ax.TickLabelInterpreter='none';
ylabel({'HbO';'Regions Run1'});
yticks(1:nparcelsdmn+nparcelsdan);
yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xticks(1:nparcelsdmn+nparcelsdan);
xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xtickangle(45);
ylabel({'HbR';'Regions Run1'});
colorbar
xlabel({'Regions Run2'});



subplot(2,3,2)
foo = mean(zHbO,3);
%lst = find(eye(size(foo))==1); foo(lst) = 0; % remove diagonal
imagesc( foo, [-1 1]*max(abs(foo(:))) )
ax = gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
title('mean Z')
xticklabels({});
yticklabels({});
colorbar

subplot(2,3,5)
foo = mean(zHbR,3);
%lst = find(eye(size(foo))==1); foo(lst) = 0; % remove diagonal
imagesc( foo, [-1 1]*max(abs(foo(:))) )
ax= gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
ax.TickLabelInterpreter='none';
xlabel({'Regions Run2'});
xticks(1:nparcelsdmn+nparcelsdan);
xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xtickangle(45);
yticklabels({});
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
ax = gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
title('-log10( p_{val}(Z T-Statistic) )')
xticklabels({});
yticklabels({});
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
ax= gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
ax.TickLabelInterpreter='none';
colorbar
xlabel({'Regions Run2'});
xticks(1:nparcelsdmn+nparcelsdan);
xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xtickangle(45);
yticklabels({});

f.Position = [10         10        980         507];
saveas(f,[pipelineDir,fOut_reliability,'.png']);
close(f);

%%
% Probability maps
[hmap] = probabilityMap_FC(dmn_improv_hbo,mesh_brain,idx_select,['Probability Group DMN HbO Map (',pipeline_str,')']);
saveas(hmap,[pipelineDir,fOut_pmap,'_dmn_hbo.png']);
close(hmap);
[hmap] = probabilityMap_FC(dan_improv_hbo,mesh_brain,idx_select,['Probability Group DAN HbO Map (',pipeline_str,')']);
saveas(hmap,[pipelineDir,fOut_pmap,'_dan_hbo.png']);
close(hmap);
[hmap] = probabilityMap_FC(dmn_improv_hbr,mesh_brain,idx_select,['Probability Group DMN HbR Map (',pipeline_str,')']);
saveas(hmap,[pipelineDir,fOut_pmap,'_dmn_hbr.png']);
close(hmap);
[hmap] = probabilityMap_FC(dan_improv_hbr,mesh_brain,idx_select,['Probability Group DAN HbR Map (',pipeline_str,')']);
saveas(hmap,[pipelineDir,fOut_pmap,'_dan_hbr.png']);
close(hmap);

%%
% Seeds/submasks plots 
fs = 10.1725;
f1=figure;
hold on;
t = tiledlayout('flow','TileSpacing','compact');
for iSubj=1:nSubjs
    %subplot(3,5,iSubj);
    nexttile;
    subject = num2str(subjects_set(iSubj));
    % for iSeed=1:length(dmn_seeds_ts(:,1,iSubj))
    %     plot(dmn_seeds_ts{iSeed,1,iSubj},'-','Color',[0.7 0.7 1]);
    % end
    plot(mean([dmn_seeds_ts{:,1,1}],2),'-','Color',[0.2 0.2 1]);    
    hold on;
    % for iSeed=1:length(dan_seeds_ts(:,1,iSubj))
    %     plot(dan_seeds_ts{iSeed,1,iSubj},'-','Color',[1.0 0.8 0.0]);
    % end
    plot(mean([dan_seeds_ts{:,1,1}],2),'-','Color',[1.0 0.5 0.0]);
    xlim([0 15000]);%round(5*60*fs)]);
    title(['S',subject]);
end
sgtitle('HbO DMN and DAN average seed time courses');
lgd=legend({'DMN','DAN'});
lgd.Layout.Tile = 'east';
f1.Position = [10         10        1049         507];
saveas(f1,[pipelineDir,'HbODMN-DAN_seedTimecourses (0-5min).png']);
close(f1);
% Seeds plots 
iSubj =1;
f=figure;
hold on;
for iSeed=1:length(dmn_seeds_ts(:,1,iSubj))
    plot(dmn_seeds_ts{iSeed,1,iSubj},'-','Color',[0.7 0.7 1]);
end
plot(mean([dmn_seeds_ts{:,1,1}],2),'-','Color',[0.2 0.2 1],'LineWidth',2.5);

iSubj =1;
hold on;
for iSeed=1:length(dan_seeds_ts(:,1,iSubj))
    plot(dan_seeds_ts{iSeed,1,iSubj},'-','Color',[1.0 0.8 0.0]);
end
plot(mean([dan_seeds_ts{:,1,1}],2),'-','Color',[1.0 0.5 0.0],'LineWidth',2.5);

%% Plot anticorrelation
f3=figure(3); 
wind_min = 5;
for iSubj=1:nSubjs
    subject = num2str(subjects_set(iSubj));
    dmn_mu = mean([dmn_seeds_ts{:,1,iSubj}],2);
    dan_mu = mean([dan_seeds_ts{:,1,iSubj}],2);
    r_seeds = [];
    window_ = round(fs*60*wind_min);
    for i=1:(length(dmn_mu)-window_)
        r=corrcoef(dmn_mu(i:(i+window_)),dan_mu(i:(i+window_)));
        r_seeds = [r_seeds;r(2)];
    end
    subplot(3,5,iSubj); 
    plot(r_seeds)
    ylim([-1 1]);
    title(['S',subject]);
end
sgtitle(['HbO DMN-DAN Anticorrelation (',num2str(wind_min),' min)']);
f3.Position = [10         10        1049         507];
saveas(f3,[pipelineDir,'HbODMN-DAN_anticorr-',num2str(wind_min),'m.png']);
close(f3);