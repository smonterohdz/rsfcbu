%%
clear all;
OVERWRITE_ = 0;

%[fwFolder,anatomFolder,derivFolder,dataDir] = setmyenv();

subjects_set = [1:7 9:16];
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

iSubj = 3;
subject = num2str(subjects_set(iSubj));
% use twindow.init_sec = -1 if you want to use the onset and duration
% defined in the snirf file (assuming there is only one stimulus).
% Otherwise, change the values (in sec.) according to your needs
twindow.stim_name = 'baseline';
twindow.init_sec = -1;
twindow.offset_sec = 0;
%twindow.init_sec = 30;
%twindow.dur_sec = 180;

pipeline_str = sprintf('%s_subj-%i_macor-%s_bpfilt-%s_imrec-%s_gsr-%s_clust-%i_parcel-%s_Asf-%i',...
    flags.task,iSubj,flags.macorrect,flags.bpfilt,flags.imagerecon,flags.gsr, ...
    flags.clusteringType,flags.parcel_scheme,flags.Adot_scaling_factor);
fOut_dmndan_mat=sprintf('DMN-DAN_mat_%s',pipeline_str);
%fOut_pmap=sprintf('probMap_%s',pipeline_str);
pipelineDir = sprintf('%sPipeline-%s/',derivFolder,pipeline_str);
if ~exist(pipelineDir,'dir')
    mkdir(pipelineDir);
end
% saveas(f1,[pipelineDir,'DMN_regions_3D.png']);
% close(f1);
% saveas(f2,[pipelineDir,'DAN_regions_3D.png']);
% close(f2);

%%
fprintf('Subject %s------------------------\n',subject);
if OVERWRITE_ || ~exist([pipelineDir,fOut_dmndan_mat,'.mat'],'file')
%%
%iRun = 1;
%for iSubj = 1:3
for iRun=1:8
    eval(sprintf("SnirfFilePathr%i = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-',flags.task,'_run-',num2str(%i),'_nirs.snirf']",iRun,iRun));
    %SnirfFilePathr2 = [dataDir,'Subj-',subject,'/nirs/sub-',subject,'_task-',flags.task,'_run-2_nirs.snirf'];
    if eval(sprintf("~exist(SnirfFilePathr%i,'file')",iRun))
        eval(sprintf("HbO_brainr%i = [];",iRun));
        eval(sprintf("HbR_brainr%i = [];",iRun));
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
end
nStim = length(snirfObjr1.stim);
BrainMaps_hbo = zeros(length(idx_select),1*(length(dmn_mask)+length(dan_mask)),nStim);
BrainMaps_hbr = zeros(length(idx_select),1*(length(dmn_mask)+length(dan_mask)),nStim);
dmn_hbo_ts = cell(nStim,1);
dan_hbo_ts = cell(nStim,1);
dmn_hbr_ts = cell(nStim,1);
dan_hbr_ts = cell(nStim,1);
dmn_improv_hbo = cell(nStim,1);
dan_improv_hbo = cell(nStim,1);
dmn_improv_hbr = cell(nStim,1);
dan_improv_hbr = cell(nStim,1);
stim_labels = {snirfObjr1.stim.name};
for iStim=1:nStim
    twindow.stim_name = snirfObjr1.stim(iStim).name;
    disp(twindow.stim_name);
    for iRun=1:8
        %%
        if eval(sprintf("~isempty(HbO_brainr%i)==1",iRun))
            eval(sprintf("[HbO_brain_chunkr%i] = ExctractChunk(HbO_brainr%i,snirfObjr%i,twindow,flags);",iRun,iRun,iRun));
            eval(sprintf("[HbR_brain_chunkr%i] = ExctractChunk(HbR_brainr%i,snirfObjr%i,twindow,flags);",iRun,iRun,iRun));
        else
            eval(sprintf("[HbO_brain_chunkr%i] = []",iRun));
            eval(sprintf("[HbR_brain_chunkr%i] = []",iRun));
        end

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

    dmn_hbo_ts{iStim} = zeros(size(HbO_brain_r1r8,1),nparcelsdmn);
    dan_hbo_ts{iStim} = zeros(size(HbO_brain_r1r8,1),nparcelsdan);
    dmn_hbr_ts{iStim} = zeros(size(HbR_brain_r1r8,1),nparcelsdmn);
    dan_hbr_ts{iStim} = zeros(size(HbR_brain_r1r8,1),nparcelsdan);

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
    [~,dmn_improv_hbo{iStim},dan_improv_hbo{iStim},~] = ClusterSelection_FC(dmn_mask_hbo,dan_mask_hbo,HbO_brain_r1r8,flags);
    [~,dmn_improv_hbr{iStim},dan_improv_hbr{iStim},~] = ClusterSelection_FC(dmn_mask_hbr,dan_mask_hbr,HbR_brain_r1r8,flags);

    %%
    % HbO
    for iSubmask=1:length(dmn_improv_hbo{iStim})
        %obtain the correlation brain map after preprocessing and by using the
        %seed passed as the first argument.
        [dmn_hbo_ts{iStim}(:,iSubmask),hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbo{iStim}(iSubmask),mesh_brain,idx_select,HbO_brain_r1r8,flags.p_thresh,flags.plot);
        BrainMaps_hbo(:,iSubmask,iStim) = A_select1;
        fprintf('(hbo)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbo{iStim}));
    end
    for iSubmask=1:length(dan_improv_hbo{iStim})
        [dan_hbo_ts{iStim}(:,iSubmask),hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbo{iStim}(iSubmask),mesh_brain,idx_select,HbO_brain_r1r8,flags.p_thresh,flags.plot);
        BrainMaps_hbo(:,length(dmn_improv_hbo{iStim})+iSubmask,iStim) = A_select1;
        fprintf('(hbo)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbo{iStim}));
    end


    % hbr
    for iSubmask=1:length(dmn_improv_hbr{iStim})
        [dmn_hbr_ts{iStim}(:,iSubmask),hmap1,A_select1] = CorrelationBrainMap_FC(dmn_improv_hbr{iStim}(iSubmask),mesh_brain,idx_select,HbR_brain_r1r8,flags.p_thresh,flags.plot);
        BrainMaps_hbr(:,iSubmask,iStim) = A_select1;
        fprintf('(hbr)DMN submask %i of %i\n',iSubmask,length(dmn_improv_hbr{iStim}));
    end
    for iSubmask=1:length(dan_improv_hbr{iStim})
        [dan_hbr_ts{iStim}(:,iSubmask),hmap1,A_select1] = CorrelationBrainMap_FC(dan_improv_hbr{iStim}(iSubmask),mesh_brain,idx_select,HbR_brain_r1r8,flags.p_thresh,flags.plot);
        BrainMaps_hbr(:,length(dmn_improv_hbr{iStim})+iSubmask,iStim) = A_select1;
        fprintf('(hbr)DAN submask %i of %i\n',iSubmask,length(dan_improv_hbr{iStim}));
    end
end
save([pipelineDir,fOut_dmndan_mat,'.mat'],'BrainMaps_hbo','BrainMaps_hbr','flags','nStim',...
    'dmn_hbo_ts','dmn_hbr_ts','dan_hbo_ts','dan_hbr_ts','iSubj','twindow','stim_labels',...
    'dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr');
else
    load([pipelineDir,fOut_dmndan_mat,'.mat'],'BrainMaps_hbo','BrainMaps_hbr','flags','nStim',...
        'dmn_hbo_ts','dmn_hbr_ts','dan_hbo_ts','dan_hbr_ts','iSubj','twindow','stim_labels',...
        'dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr');
end
r_hbo = zeros(nparcelsdmn+nparcelsdan,nparcelsdmn+nparcelsdan,nStim);
r_hbr = zeros(nparcelsdmn+nparcelsdan,nparcelsdmn+nparcelsdan,nStim);
for iStim=1:nStim
    % Correlaition matrix R
    fsize = 12;
    fpos =   [458.6000   66.6000  827.2000  686.4000];
    f=figure;
    f.Position = fpos;
%    subplot(1,2,1);
    %imagesc(corrcoef(BrainMaps_hbo(:,:,iStim)),[-1,1]);
    [r_hbo(:,:,iStim),p0_hbo] = corrcoef(BrainMaps_hbo(:,:,iStim));
    z_hbo = 0.5 * log((1+r_hbo(:,:,iStim))./(1-r_hbo(:,:,iStim)));
    z_hbo(z_hbo>6) = 6; z_hbo(z_hbo<-6) = -6;
    n = length(idx_select);
    se_hbo = 1/(sqrt(n-3));
    zscore_hbo = z_hbo ./ se_hbo;
    pval_hbo= 2*(1-normcdf(abs(zscore_hbo)));
    pval_hbolog10 = -log10(pval_hbo);
    pval_hbolog10 = pval_hbolog10.*(pval_hbolog10>2).*sign(zscore_hbo);
    imagesc(z_hbo,[-1 1].*max(abs(z_hbo).*(~eye(nparcelsdmn+nparcelsdan)),[],'all','omitnan'));
    colormap("jet");
    colorbar();
    %title({sprintf('Subject %s DMN-DAN HbO Run 1-8',subject),stim_labels{iStim},pipeline_str},'Interpreter','none');
    title({sprintf('Subject %s DMN-DAN HbO Run 1-8',subject),stim_labels{iStim},'Fisher z and p-values'});
    yticks(1:nparcelsdmn+nparcelsdan);
    yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
        {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
    xticks(1:nparcelsdmn+nparcelsdan);
    xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
        {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
    xtickangle(45);
    ax= gca();
    ax.XAxis.FontSize = fsize;
    ax.YAxis.FontSize = fsize;
    ax.TickLabelInterpreter='none';
    for iX=1:(nparcelsdan+nparcelsdmn)
        for iY=1:(nparcelsdan+nparcelsdmn)
            %text(iX,iY,num2str(pval_hbo(iY,iX)));
            if pval_hbo(iY,iX)<0.01
                text(iX,iY,'*');
            else
                text(iX,iY,sprintf('%.2f',pval_hbo(iY,iX)),'HorizontalAlignment','center');
            end
        end
    end
    saveas(f,[pipelineDir,fOut_dmndan_mat,'_Subj-',num2str(iSubj),'_WMseed_',stim_labels{iStim},'_hbo.png']);
    close(f);

    % Pvalues of previous correl matrix
    % subplot(1,2,2);
    % n = length(idx_select);
    % se_hbo = 1/(sqrt(n-3));
    % zscore_hbo = z_hbo ./ se_hbo;
    % pval_hbo= 2*(1-normcdf(abs(zscore_hbo)));
    % pval_hbolog10 = -log10(pval_hbo);
    % pval_hbolog10 = pval_hbolog10.*(pval_hbolog10>2).*sign(zscore_hbo);
    % imagesc(pval_hbo,[-1 1]*max(abs(pval_hbo(:))));
    % colormap("jet");
    % cb=colorbar();
    % title({'-log10(pval(Fisher-Z))'},'Interpreter','none');
    % % yticks(1:nparcelsdmn+nparcelsdan);
    % % yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    % %     {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
    % xticks(1:nparcelsdmn+nparcelsdan);
    % xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    %     {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
    % xtickangle(45);
    % ax= gca();
    % ax.XAxis.FontSize = fsize;
    % ax.YAxis.FontSize = fsize;
    % ax.TickLabelInterpreter='none'; 
    % saveas(f,[pipelineDir,fOut_dmndan_mat,'_Subj-',num2str(iSubj),'_WMseed_',stim_labels{iStim},'_hbo.png']);
    % close(f);

    f=figure;
    plot(mean(dmn_hbo_ts{iStim},2),'-','Color',[0 0 1],'LineWidth',2.5);
    hold on;
    plot(mean(dan_hbo_ts{iStim},2),'-','Color',[1 0 0],'LineWidth',2.5);
    plot(dmn_hbo_ts{iStim},'-','Color',[0.7 0.7 1]);
    plot(dan_hbo_ts{iStim},'-','Color',[1 0.7 0.7]);
    xlabel('Samples');
    ylabel('\Delta HbO')
    legend({'avg DMN','avg DAN'});
    title({sprintf('Subject %s DMN-DAN HbO time courses',subject),stim_labels{iStim},pipeline_str},'Interpreter','none');
    saveas(f,[pipelineDir,fOut_dmndan_mat,'_Subj-',num2str(iSubj),'_WMseed_',stim_labels{iStim},'_hbo_ts.png']);
    close(f);




    f=figure;
    f.Position = fpos;
    [r_hbr(:,:,iStim),p0_hbr] = corrcoef(BrainMaps_hbr(:,:,iStim));
    z_hbr = 0.5 * log((1+r_hbr(:,:,iStim))./(1-r_hbr(:,:,iStim)));
    z_hbr(z_hbr>6) = 6; z_hbr(z_hbr<-6) = -6;
    %n = length(idx_select);
    se_hbr = 1/(sqrt(n-3));
    zscore_hbr = z_hbr ./ se_hbr;
    pval_hbr= 2*(1-normcdf(abs(zscore_hbr)));
    pval_hbrlog10 = -log10(pval_hbr);
    pval_hbrlog10 = pval_hbrlog10.*(pval_hbrlog10>2).*sign(zscore_hbr);
    imagesc(z_hbr,[-1 1].*max(abs(z_hbr).*(~eye(nparcelsdmn+nparcelsdan)),[],'all','omitnan'));
    colormap("jet");
    colorbar();
    %title({sprintf('Subject %s DMN-DAN HbR Run 1-8',subject),stim_labels{iStim},pipeline_str},'Interpreter','none');
    title({sprintf('Subject %s DMN-DAN HbR Run 1-8',subject),stim_labels{iStim},'Fisher z and p-values'});
    yticks(1:nparcelsdmn+nparcelsdan);
    yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
        {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
    xticks(1:nparcelsdmn+nparcelsdan);
    xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
        {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
    xtickangle(45);
    ax= gca();
    ax.XAxis.FontSize = fsize;
    ax.YAxis.FontSize = fsize;
    ax.TickLabelInterpreter='none';
    for iX=1:(nparcelsdan+nparcelsdmn)
        for iY=1:(nparcelsdan+nparcelsdmn)
            %text(iX,iY,num2str(pval_hbo(iY,iX)));
            if pval_hbo(iY,iX)<0.01
                text(iX,iY,'*');
            else
                text(iX,iY,sprintf('%.2f',pval_hbo(iY,iX)),'HorizontalAlignment','center');
            end
        end
    end
    saveas(f,[pipelineDir,fOut_dmndan_mat,'_Subj-',num2str(iSubj),'_WMseed_',stim_labels{iStim},'_hbr.png']);
    close(f);

    f=figure;
    plot(mean(dmn_hbr_ts{iStim},2),'-','Color',[0 0 1],'LineWidth',2.5);
    hold on;
    plot(mean(dan_hbr_ts{iStim},2),'-','Color',[1 0 0],'LineWidth',2.5);
    plot(dmn_hbr_ts{iStim},'-','Color',[0.7 0.7 1]);
    plot(dan_hbr_ts{iStim},'-','Color',[1 0.7 0.7]);
    xlabel('Samples');
    ylabel('\Delta HbR')
    legend({'avg DMN','avg DAN'});
    title({sprintf('Subject %s DMN-DAN HbR time courses',subject),stim_labels{iStim},pipeline_str},'Interpreter','none');
    saveas(f,[pipelineDir,fOut_dmndan_mat,'_Subj-',num2str(iSubj),'_WMseed_',stim_labels{iStim},'_hbr_ts.png']);
    close(f);
end

[r_hbo_ActPas,~] = corrcoef([BrainMaps_hbo(:,:,1),BrainMaps_hbo(:,:,2)]);
[r_hbr_ActPas,~] = corrcoef([BrainMaps_hbr(:,:,1),BrainMaps_hbr(:,:,2)]);

save([pipelineDir,fOut_dmndan_mat,'.mat'],'BrainMaps_hbo','BrainMaps_hbr','flags','nStim',...
    'dmn_hbo_ts','dmn_hbr_ts','dan_hbo_ts','dan_hbr_ts','iSubj','twindow','stim_labels',...
    'dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr',...
    'r_hbo','r_hbr','r_hbo_ActPas','r_hbr_ActPas');
save([pipelineDir,'output.mat'],'BrainMaps_hbo','BrainMaps_hbr','flags','nStim',...
    'dmn_hbo_ts','dmn_hbr_ts','dan_hbo_ts','dan_hbr_ts','iSubj','twindow','stim_labels',...
    'dmn_improv_hbo','dan_improv_hbo','dmn_improv_hbr','dan_improv_hbr',...
    'r_hbo','r_hbr','r_hbo_ActPas','r_hbr_ActPas');
%% after runing Rscript
fsize=11;
load([pipelineDir,'R_output.mat']);
f=figure();
subplot(1,4,1);
%imagesc(r_hbo(:,:,iStim));
z_hbo = 0.5 * log((1+r_hbo(:,:,1))./(1-r_hbo(:,:,1)));
z_hbo(z_hbo>6) = 6; z_hbo(z_hbo<-6) = -6;
imagesc(tril(z_hbo,-1),[-1 1].*max(abs(z_hbo).*(~eye(nparcelsdmn+nparcelsdan)),[],'all','omitnan'));
colormap("jet");
cb=colorbar();
cb.Label.String = 'Fisher Z val';
title('Passive');
yticks(1:nparcelsdmn+nparcelsdan);
yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xticks(1:nparcelsdmn+nparcelsdan);
xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xtickangle(45);
ylabel(['Subject ',subject]);
ax= gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
ax.TickLabelInterpreter='none';
axis image;



subplot(1,4,2);
%imagesc(r_hbo(:,:,iStim));
z_hbo = 0.5 * log((1+r_hbo(:,:,2))./(1-r_hbo(:,:,2)));
z_hbo(z_hbo>6) = 6; z_hbo(z_hbo<-6) = -6;
imagesc(tril(z_hbo,-1),[-1 1].*max(abs(z_hbo).*(~eye(nparcelsdmn+nparcelsdan)),[],'all','omitnan'));
colormap("jet");
cb=colorbar();
cb.Label.String = 'Fisher Z val';
title('Active')
% yticks(1:nparcelsdmn+nparcelsdan);
% yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
%     {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xticks(1:nparcelsdmn+nparcelsdan);
xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xtickangle(45);
ax= gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
ax.TickLabelInterpreter='none';
axis image;

subplot(1,4,3);
imagesc(z_mat(1:10,1:10));
colormap('jet');
clim([-max(z_mat(1:10,1:10),[],'all'),max(z_mat(1:10,1:10),[],'all')]);
cb=colorbar();
cb.Label.String = 'Fisher Z val';
title('Passive-Active');
% yticks(1:nparcelsdmn+nparcelsdan);
% yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
%     {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xticks(1:nparcelsdmn+nparcelsdan);
xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xtickangle(45);
ax= gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
ax.TickLabelInterpreter='none';
axis image;


sp=subplot(1,4,4);
imagesc(-log10(p_mat(1:10,1:10)));
colormap(sp,'hot');
clim([0,2]);
cb=colorbar();
cb.Label.String = '-log_{10}(p-val)';
title('Passive-Active');
f.Position = [5    290    15392    373];
% yticks(1:nparcelsdmn+nparcelsdan);
% yticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
%     {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xticks(1:nparcelsdmn+nparcelsdan);
xticklabels(replace([{dmn_mask.name},{dan_mask.name}], ...
    {'7Networks_','Default','DorsAttn'},{'','DMN','DAN'}));
xtickangle(45);
ax= gca();
ax.XAxis.FontSize = fsize;
ax.YAxis.FontSize = fsize;
ax.TickLabelInterpreter='none';
axis image;
saveas(f,[pipelineDir,'Subj-',subject,'Passive-Active_hbo.png']);
close(f);