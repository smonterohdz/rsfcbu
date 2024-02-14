function [dmn_regions,dan_regions,mesh_brain,idx_select] = Parcellation_FC(anatomFolder,fwFolder,flags)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% Load Brain Adot
load([fwFolder 'Adot.mat'],'Adot');%
index_select = log10(sum(Adot(:,:,1),1))>=-2;
idx_select = find(index_select);

% Load brain mesh
load([fwFolder, 'mesh_brain.mat'],'mesh');
mesh_brain = mesh;


Adot_select = Adot(:,index_select,1);
Adot_new = zeros(size(Adot(:,:,1)));
Adot_new(:,index_select,:) = Adot_select(:,:,1);
A = log10(sum(Adot_new(:,:,1),1));
figure;
plot_sensitivity(mesh_brain,A, [-2 0])


%load transformation matrix and labels
T_2vol = load([anatomFolder 'labelssurf2vol.txt'],'ascii');


load([anatomFolder 'labelssurf.mat'],'aal_fv','aal_ll');

fv = struct('vertices',[],'faces',[],'idxlab',[]);
idxL = [];
mylabels = {};
for ii=1:length(aal_fv)
    if ~isempty(aal_fv{ii}.vertices)
        fv.faces    = [fv.faces; aal_fv{ii}.faces+size(fv.vertices,1)];
        aal_fv{ii}.vertices = xform_apply(aal_fv{ii}.vertices, T_2vol);
        fv.vertices = [fv.vertices; aal_fv{ii}.vertices];
        idxL = [idxL; ii*ones(size(aal_fv{ii}.vertices,1),1)];
        mylabels = [mylabels;repmat(string(aal_ll{ii}),size(aal_fv{ii}.vertices,1),1)];
        %fv.labels = {fv.labels;aal_ll{ii}}
    end
end
 fv.idxlab = idxL;
 fv.labels = mylabels;


% Finding the closest vertices between meshes
d=pdist2(mesh_brain.vertices(index_select,[2,1,3]),fv.vertices(:,[2,1,3]),'fasteuclidean');
[~,idxmin] = min(d,[],2);
brainvL = idxL(idxmin);

if 0
    save('brainvL.mat','brainvL','fv');
end


% Creating submasks

id_rois = unique(brainvL);
n_parcels = length(id_rois);
seed_regions = struct('name',[],'vertices',[],'vertices_index',[]);
dmn_regions = struct('name',[],'vertices',[],'vertices_index',[]);
brain_vertices_select = mesh_brain.vertices(index_select,:);
dan_regions = struct('name',[],'vertices',[],'vertices_index',[]);

%dmn seeds:
% my_regions_lab = {'Frontal_Med_Orb_L','Frontal_Med_Orb_R',...
%     'Frontal_Sup_Orb_L','Frontal_Sup_Orb_R',...%    'Frontal_Sup_L','Frontal_Sup_R',...
%     'Frontal_Sup_Medial_L','Frontal_Sup_Medial_R',...
%     'Angular_L','Angular_R'};
%dan seeds:
% my_regions_lab = {'Parietal_Sup_L','Parietal_Sup_R',...
%      'Precentral_L','Precentral_R'};

dmn_regions_lab = {'Frontal_Med_Orb_L','Frontal_Med_Orb_R',...
    'Frontal_Sup_Orb_L','Frontal_Sup_Orb_R',...
    'Frontal_Mid_Orb_L','Frontal_Mid_Orb_R'...
    'Frontal_Sup_L','Frontal_Sup_R',...
    'Frontal_Sup_Medial_L','Frontal_Sup_Medial_R',...
    'Temporal_Sup_L','Temporal_Sup_R',...
    'Angular_L','Angular_R',...
    'Parietal_Inf_L','Parietal_Inf_R'};

dan_regions_lab = {'Parietal_Sup_L','Parietal_Sup_R',...
    'Precentral_L','Precentral_R'};


k=1;
l=1;

for j=1:n_parcels
    if any(strcmp(aal_ll{id_rois(j)},dmn_regions_lab))
        dmn_regions(k).name = aal_ll{id_rois(j)};
        dmn_regions(k).vertices_index = find(brainvL==id_rois(j));
        dmn_regions(k).vertices = brain_vertices_select(dmn_regions(k).vertices_index,:);
        dmn_regions(k).type = "mask";
        dmn_regions(k).net_lab = "dmn";
        k = k+1;
    end
    if any(strcmp(aal_ll{id_rois(j)},dan_regions_lab))
        dan_regions(l).name = aal_ll{id_rois(j)};
        dan_regions(l).vertices_index = find(brainvL==id_rois(j));
        dan_regions(l).vertices = brain_vertices_select(dan_regions(l).vertices_index,:);
        dan_regions(l).type = "mask";
        dan_regions(l).net_lab = "dan";
        l = l+1;
    end
end


end