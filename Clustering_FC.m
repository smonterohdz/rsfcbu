function [newmask] = Clustering_FC(HbXBrain,mask,flags)
%Clustering_FC Genereic function to apply clustering to each submask.

%r_thresh = 0.8;
newmask = [];
for iSubmask = 1:length(mask)
    submask = mask(iSubmask);
    seedHbX_ts = zscore(HbXBrain(:,submask.vertices_index));
    if flags.clusteringType == 1
        [newsubmask] = ImproveSeed_FC_1(seedHbX_ts,submask,flags.r_thresh);
    else
        [newsubmask] = ImproveSeed_FC_2(seedHbX_ts,submask);
    end
    newmask = [newmask;newsubmask];
end

end


