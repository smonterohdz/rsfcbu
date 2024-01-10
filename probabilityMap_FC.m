function [hmap] = probabilityMap_FC(mask,mesh_brain,idx_select,maskLabel)
%probabilityMap creates the spatial probability map for a mask from the
%improved seeds from each subect.
A = zeros(size(mesh_brain.vertices,1),1);
A_select = zeros(length(idx_select),1);

nSubj = length(mask);
nSubMasks = length(mask{1}); %assuming all the subjects contain the same number of submask.
for iSubj =1 :nSubj
    for iSM = 1:nSubMasks
        nodes_ = mask{iSubj}(iSM).vertices_index(mask{iSubj}(iSM).mask_subsetseed);
        A_select(nodes_) = A_select(nodes_) + 1;
    end
end

A(idx_select) = A_select./nSubj.*100;
hmap = plot_connectivity_seedregion(mesh_brain,A,[0 100],[],[],[],[],[]);
cb=colorbar;
cb.Label.String='Spatial Correspondence (%)';
sgtitle(maskLabel);
hold on;
end