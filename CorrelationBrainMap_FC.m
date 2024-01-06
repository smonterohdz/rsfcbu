function [seed_ts,hmap,A_select] = CorrelationBrainMap_FC(submask,mesh_brain,idx_select,HbXBrain,p_thresh,plot_flag)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

A_select = zeros(length(idx_select),1);

%submask.mask_subsetseed = submask.groups == clusterStar;
seed_ts = mean(HbXBrain(:,submask.vertices_index(submask.mask_subsetseed)),2);


non_seed = HbXBrain; non_seed(:,submask.vertices_index(submask.mask_subsetseed)) = [];
R_HbO =[];
P_HbO =[];
for kk = 1:size(non_seed,2)
    [R,P] = corrcoef(seed_ts, non_seed(:,kk));
    R_HbO(kk) = R(2);P_HbO(kk) = P(2);
end
R_HbO_f = 0.5*log((1+R_HbO)./(1-R_HbO));
hmap =[];
% HbO

A_select(submask.vertices_index(submask.mask_subsetseed)) = 0;
non_seed_index = 1:size(HbXBrain,2);
non_seed_index(submask.vertices_index(submask.mask_subsetseed)) = [];
if p_thresh ~= 0
    A_select(non_seed_index(P_HbO<p_thresh)) = R_HbO_f(P_HbO<p_thresh);
else
    A_select(non_seed_index) = R_HbO_f;
end
if plot_flag ~= 0
    A = zeros(1,size(mesh_brain.vertices,1));
    A(idx_select) = A_select;
    hmap = plot_connectivity_seedregion(mesh_brain,A,[],idx_select(submask.vertices_index(submask.mask_subsetseed)), ...
        seed_ts,[],[],[]);
    hold on;
    %title(hmap.Children.Children(5),seed_region.name,'Interpreter','none');
    %title(hmap.Children.Children(1),['Seed Hb timeseries']);
    %save_mat_file = fullfile(save_mat_path,  [filename(1:end-6), ...
    %    '_',seed_region.name,'_HbO_corrMap-starCl.png']);
    %saveas(h,save_mat_file);
    %close(h);
end
end