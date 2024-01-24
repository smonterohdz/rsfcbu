function [g, rDMNDAN,mask_improv] = findGroupOfMostCorrelatedClusters_func( cluster_idx, r_run, ts_run, flags, iSortChoice, dmn_regions, dan_regions, g )




if strcmp(flags.Mask,'DMN') % DMN
    iCluster0 = 0;
    iSubMask0 = 0;
    nSubMask = length(dmn_regions);
else % DAN
    iCluster0 = sum(cluster_idx(1:length(dmn_regions)));
    iSubMask0 = length(dmn_regions);
    nSubMask = length(dan_regions);
end


%
% NO NEED TO CHANGE ANYTHING BELOW
%
nCluster = sum(cluster_idx(iSubMask0+[1:nSubMask]));


% Get the number of vertices within each cluster
nVerticesInCluster = zeros(1,nCluster);
for iCluster = 1:nCluster
    iSubMask = find(iCluster<=cumsum(cluster_idx(iSubMask0+[1:nSubMask]) ),1);
    if strcmp(flags.Mask,'DMN')
        %eval( sprintf('foo = dmn_regions(iSubMask).clustersidx%s%d;',flags.HbX,flags.Run) );
        foo = dmn_regions(iSubMask).groups;
    else
        %eval( sprintf('foo = dan_regions(iSubMask).clustersidx%s%d;',flags.HbX,flags.Run) );
        foo = dan_regions(iSubMask).groups;
    end
    nVerticesInCluster(iCluster) = length(find( foo==(iCluster-sum(cluster_idx(iSubMask0+[1:(iSubMask-1)]))) ));
end

% Loop over all the clusters
% For each cluster, find the cluster from each of the other sub-masks that
% best correlates with it. This gives a set of clusters for the mask. We
% then calculate the mean distance (as 1-R) between all of those clusters
% within the set.
% We will nominally identify the set of clusters with the minimal mean
% distance to represent the target network within the given mask.


clDistance = zeros(nCluster,1);
clSubMask = zeros(nCluster,1);
clTimeSeries = zeros(size(ts_run,1),nCluster);
r12 = zeros(nCluster,1);

% loop over all clusters
for iCluster1 = 1:nCluster
    iSubMask1 = find(iCluster1<=cumsum(cluster_idx(iSubMask0+[1:nSubMask]) ),1);
    
    lstClusters = zeros(nSubMask,1);
    lstClusters(iSubMask1) = iCluster1;

    % loop over the sub-masks to identify the cluster within each sub-mask
    % that best correlates with the cluster iCluster1
    for iSubMask2 = 1:nSubMask
        nCluster2 = cluster_idx(iSubMask0+iSubMask2);
        i2 = sum(cluster_idx(iSubMask0+[1:(iSubMask2-1)]));
        
        r_tmp = r_run(iCluster0+iCluster1,iCluster0+i2+[1:nCluster2]);
        ii = find(r_tmp == max(r_tmp));
        
        lstClusters(iSubMask2) = i2+ii;
    end        
    
    % plot intermediate results (used for testing initially)
    if 0
        r_tmp = zeros(1,nCluster);
        r_tmp(lstClusters) = r_run(iCluster0+iCluster1,iCluster0+lstClusters);
    
        figure(6)
        subplot(2,1,1)
        imagesc( [r_run(iCluster0+iCluster1,iCluster0+[1:nCluster]); r_tmp], [-1 1])
        subplot(2,1,2)
        imagesc( r_run(iCluster0+lstClusters,iCluster0+lstClusters), [-1 1] )
        pause
    end

    % record the distance for the set of clusters derived using the
    % iCluster1 seed. Note that we can use the mean of R or of R.^2. I like
    % the mean of R.^2 as this penalizes outliers with a higher distance
    % more.
%    clDistance(iCluster1) = mean(mean( 1-r_run(iCluster0+lstClusters,iCluster0+lstClusters) ));
    clDistance(iCluster1) = mean(mean( (1-r_run(iCluster0+lstClusters,iCluster0+lstClusters)).^2 ))^0.5;
    clSubMask(iCluster1) = iSubMask1;

    % record the mean time series for this set of clusters
    % But weight each cluster time series by the number of vertices within
    % that cluster. We want the average time series across clusters to be
    % weighted by the number of vertices within each cluster.
    clTimeSeries(:,iCluster1) = sum( ts_run(:,iCluster0+lstClusters).*(ones(size(ts_run,1),1)*nVerticesInCluster(lstClusters)), 2) / sum(nVerticesInCluster(lstClusters));
    % this is just the mean
    % clTimeSeries(:,iCluster1) = sum( ts_run(:,iCluster0+lstClusters), 2) / nCluster;
    
    
    % for DAN, we want to know the correlation of the the cluster time
    % series with the chosen DMN cluster time series. As there are sets of
    % DAN clusters that positively correlate with the DMN set of clusters
    % and there are sets of DAN clusters that negatively correlate with the
    % DMN set of clusters, and they have similar distances within their
    % cluster sets, we want to chose the set of DAN clusters that
    % negatively correlates with the DMN set of clusters.
    if ~strcmp(flags.Mask,'DMN')
        foo = corrcoef(g.clTimeSeriesChosen1, clTimeSeries(:,iCluster1));
        r12(iCluster1) = foo(1,2);
    end
end

% if iSortChoice = 0 then chose the most negatively correlated set of
% clusters. 
% we add distance to the metric since we have seen the case where one set
% of clusters was more negatively correlated but it had a very poor (i.e.
% large) distance within the set of clusters. It is better to chose a set
% of clusters that is slightly less negative but has a closer distance.
if ~strcmp(flags.Mask,'DMN') & iSortChoice==0
    [r12_min,iSortChoice] = min(r12+clDistance); 
    iSortChoice = -1 * iSortChoice;
end


if 0
    foo=corrcoef(clTimeSeries);
    figure(1)
    clf
    subplot(2,1,1)
    for iSubMask1 = 1:nSubMask
        nCluster1 = cluster_idx(iSubMask0+iSubMask1);
        i1 = sum(cluster_idx(iSubMask0+[1:(iSubMask1-1)]));
        
        plot( i1+[1:nCluster1], clDistance(i1+[1:nCluster1]), '.-' )
        hold on
    end
    hold off
    xlim([1 nCluster])
    
    subplot(2,1,2)
    imagesc( foo, [-1 1]);
    colormap(jet(64))
    title('CorrCoef of mean time series of set of clusters')
end


% if strcmp(flags.Mask,'DMN')
%     figure(2)
% else
%     figure(12)
% end
% if 0
%     subplot(1,1,1)
%     imagesc( r_run( iCluster0+[1:nCluster], iCluster0+[1:nCluster] ) , [-1 1])
%     colormap(jet(64))
%     title('CorrCoef of clusters')
% else
%     lstMin = find(clDistance==min(clDistance));
%     [sortDist,lstMin] = sort(clDistance,'ascend');
% 
%     ha2=subplot(2,2,2);
%     for iSubMask1 = 1:nSubMask
%         nCluster1 = cluster_idx(iSubMask0+iSubMask1);
%         i1 = sum(cluster_idx(iSubMask0+[1:(iSubMask1-1)]));
% 
%         plot( i1+[1:nCluster1], clDistance(i1+[1:nCluster1]), '.-' )
%         hold on
%     end
%     hold off
%     xlim([1 nCluster])
%     ylabel( 'Distance' )
%     if iSortChoice>0
%         title( sprintf('Best Cluster=%d; Choosen Cluster=%d',lstMin(1),lstMin(iSortChoice)) )
%     else
%         title( sprintf('Best Cluster=%d; Choosen Cluster=%d',lstMin(1),-iSortChoice) )
%     end
% 
%     foo=corrcoef(clTimeSeries);
%     subplot(2,2,4)
%     plot( 1:nCluster, foo(:,lstMin(abs(1))), '.-' )
%     xlim([1 nCluster])
%     title('CorrCoef of mean time series of set of clusters w Best Cluster')
%     ylabel( 'CorrCoef' )
%     if ~strcmp(flags.Mask,'DMN')
%         hold on
%         hl=plot(r12);
%         hold off
%         legend('w Best Cluster','w DMN Cluster')
%     end
% 
% 
%     subplot(2,2,3)
%     hl=plot( clDistance, foo(:,lstMin(abs(1))), '.' );
%     set(hl,'markersize',16)
%     xlim(get(ha2,'ylim') )
%     ylim([-1 1])
%     ylabel( 'CorrCoef' )
%     xlabel( 'Distance' )
% 
%     subplot(2,2,1)
%     hl=plot( foo(:,lstMin(abs(1))),clDistance, '.' );
%     set(hl,'markersize',16)
%     ylim(get(ha2,'ylim') )
%     xlim([-1 1])
%     xlabel( 'CorrCoef' )
%     ylabel( 'Distance' )
%     title( flags.Mask )
% end



% plot result for best choice of clusters
lstMin = find(clDistance==min(clDistance));
[sortDist,lstMin] = sort(clDistance,'ascend');


if iSortChoice>0
    iCluster1 = lstMin(iSortChoice);
else
    iCluster1 = -iSortChoice;
end

iSubMask1 = find(iCluster1<=cumsum(cluster_idx(iSubMask0+[1:nSubMask]) ),1);

lstClusters = zeros(nSubMask,1);
lstClusters(iSubMask1) = iCluster1;

for iSubMask2 = 1:nSubMask
    nCluster2 = cluster_idx(iSubMask0+iSubMask2);
    i2 = sum(cluster_idx(iSubMask0+[1:(iSubMask2-1)]));
    
    r_tmp = r_run(iCluster0+iCluster1,iCluster0+i2+[1:nCluster2]);
    ii = find(r_tmp == max(r_tmp));
    
    lstClusters(iSubMask2) = i2+ii;
    bestCluster_id = ii; % I could use ii but it is not a descriptive var name
    if strcmp(flags.Mask,'DMN')
        %eval( sprintf('dmn_regions(iSubMask2).mask_subsetseed_%s%d = dmn_regions(iSubMask2).clustersidx%s%d== bestCluster_id;',flags.HbX,flags.Run,flags.HbX,flags.Run) );
        dmn_regions(iSubMask2).mask_subsetseed=dmn_regions(iSubMask2).groups== bestCluster_id;%added by SMH
    else
        %eval( sprintf('dan_regions(iSubMask2).mask_subsetseed_%s%d = dan_regions(iSubMask2).clustersidx%s%d== bestCluster_id;',flags.HbX,flags.Run,flags.HbX,flags.Run) );
        dan_regions(iSubMask2).mask_subsetseed=dan_regions(iSubMask2).groups== bestCluster_id;%added by SMH
    end
end

if 0
    figure(3)
    imagesc( r_run( iCluster0+lstClusters, iCluster0+lstClusters ), [-1 1]);
    colormap(jet(64))
end

[lstClusters'; nVerticesInCluster(lstClusters)]

% Store stuff for final figures
if strcmp(flags.Mask,'DMN')
    g.lstC1 = iCluster0 + lstClusters;
    g.clTimeSeries1 = clTimeSeries;
    g.clTimeSeriesChosen1 = mean( ts_run(:,iCluster0+lstClusters), 2);    
    mask_improv =dmn_regions;
else
    g.lstC2 = iCluster0 + lstClusters;
    g.clTimeSeries2 = clTimeSeries;
    g.clTimeSeriesChosen2 = mean( ts_run(:,iCluster0+lstClusters), 2);
    mask_improv=dan_regions;
end

if strcmp(flags.Mask,'DMN')
    rDMNDAN = [];
else
    rDMNDAN = r_run( [g.lstC1; g.lstC2], [g.lstC1; g.lstC2] );
end

% if ~strcmp(flags.Mask,'DMN')
%     figure(4)
%     rDMNDAN = r_run( [g.lstC1; g.lstC2], [g.lstC1; g.lstC2] );
%     imagesc( rDMNDAN, [-1 1]);
%     colormap(jet(64))
% 
%     if 0
%         figure(5)
%         imagesc( corrcoef( [g.clTimeSeries1 g.clTimeSeries2] ), [-1 1]);
%         colormap(jet(64))    
%         title('CorrCoef of mean time series of set of clusters')
%     end
% end

