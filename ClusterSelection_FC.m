function [g,dmn_improv,dan_improv,rDMNDAN] = ClusterSelection_FC(dmn_mask,dan_mask,HbXBrain,flag)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
newmask =[dmn_mask;dan_mask];
cluster_idx = zeros(1,length(newmask));
all_clusters_ts = zeros(size(HbXBrain,1),length(newmask));
iCluster = 1;
for iSubmask =1:length(newmask)
    cluster_idx(iSubmask) = length(unique(newmask(iSubmask).groups));
    for iClusterSubmask =1:cluster_idx(iSubmask)
        lst = newmask(iSubmask).groups == iClusterSubmask;
        all_clusters_ts(:,iCluster) = mean(HbXBrain(:,newmask(iSubmask).vertices_index(lst)),2);
        iCluster = iCluster+1;
    end
end
r_run = corrcoef(all_clusters_ts);

flag.Mask= 'DMN';
iSortChoice = 1;
[g,~,dmn_improv] = findGroupOfMostCorrelatedClusters_func( cluster_idx, r_run, all_clusters_ts, flag, iSortChoice, dmn_mask, dan_mask );

flag.Mask= 'DAN';
iSortChoice = 0;
[g, rDMNDAN,dan_improv] = findGroupOfMostCorrelatedClusters_func( cluster_idx, r_run, all_clusters_ts, flag, iSortChoice, dmn_mask, dan_mask, g );
end