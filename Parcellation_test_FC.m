function [dmn_regions,dan_regions,mesh_brain,idx_select] = Parcellation_test_FC(anatomFolder,fwFolder,flags)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% Load Brain Adot
load([fwFolder 'Adot.mat'],'Adot');%
index_select = log10(sum(Adot(:,:,1),1))>=-2;
idx_select = find(index_select);

Adot_select = Adot(:,index_select,1);
Adot_new = zeros(size(Adot(:,:,1)));
Adot_new(:,index_select,:) = Adot_select(:,:,1);
A = log10(sum(Adot_new(:,:,1),1));

% Load brain mesh
load([fwFolder, 'mesh_brain.mat'],'mesh_brain');
figure;
plot_sensitivity(mesh_brain,A, [-2 2])



%% reading schaefer parcels from mat file
load([fwFolder,'schaeferParcels.mat'],'rhParcelsDMN','lhParcelsDMN','rhParcelsDAN','lhParcelsDAN');

%read left and right hemispheres and concatenate the vertices and faces.
% load schaefer parcels
[vl,fl]=read_surf([fwFolder,'..\segment\lh.pial']);
[vr,fr]=read_surf([fwFolder,'..\segment\rh.pial']);
v = [vr;vl];
f = [fr;(fl+size(vr,1))];


% %plot the brain
% figure
% trisurf(f,v(:,1),v(:,2),v(:,3),'edgecolor','flat');


nParcelsDMNr = size(rhParcelsDMN,1);
nParcelsDMNl = size(lhParcelsDMN,1);
nParcelsDANr = size(rhParcelsDAN,1);
nParcelsDANl = size(lhParcelsDAN,1);

dmn_regions = struct('name',[],'vertices',[],'vertices_index',[]);
brain_vertices_select = mesh_brain.vertices(index_select,:);
dan_regions = struct('name',[],'vertices',[],'vertices_index',[]);

zeros_r = zeros(length(vr),1);
zeros_l = zeros(length(vl),1);
k = 1;
for iParcel =1:nParcelsDMNr
        dmn_regions(k).name = rhParcelsDMN{iParcel,1};
        current_mask = [rhParcelsDMN{iParcel,2};zeros_l];
        dmn_regions(k).vertices_index = find(current_mask(index_select));
        dmn_regions(k).vertices = brain_vertices_select(dmn_regions(k).vertices_index,:);
        dmn_regions(k).mask_subsetseed = true(length(dmn_regions(k).vertices_index),1);
        dmn_regions(k).type = "mask";
        dmn_regions(k).net_lab = "dmn";
        k = k+1;
end

for iParcel =1:nParcelsDMNl
        dmn_regions(k).name = lhParcelsDMN{iParcel,1};
        current_mask = [zeros_r;lhParcelsDMN{iParcel,2}];
        dmn_regions(k).vertices_index = find(current_mask(index_select));
        dmn_regions(k).vertices = brain_vertices_select(dmn_regions(k).vertices_index,:);
        dmn_regions(k).mask_subsetseed = true(length(dmn_regions(k).vertices_index),1);
        dmn_regions(k).type = "mask";
        dmn_regions(k).net_lab = "dmn";
        k = k+1;
end

l = 1;
for iParcel =1:nParcelsDANr
        dan_regions(l).name = rhParcelsDAN{iParcel,1};
        current_mask = [rhParcelsDAN{iParcel,2};zeros_l];
        dan_regions(l).vertices_index = find(current_mask(index_select));
        dan_regions(l).vertices = brain_vertices_select(dan_regions(l).vertices_index,:);
        dan_regions(l).mask_subsetseed = true(length(dan_regions(l).vertices_index),1);
        dan_regions(l).type = "mask";
        dan_regions(l).net_lab = "dan";
        l = l+1;
end

for iParcel =1:nParcelsDANl
        dan_regions(l).name = rhParcelsDAN{iParcel,1};
        current_mask = [zeros_r;rhParcelsDAN{iParcel,2}];
        dan_regions(l).vertices_index = find(current_mask(index_select));
        dan_regions(l).vertices = brain_vertices_select(dan_regions(l).vertices_index,:);
        dan_regions(l).mask_subsetseed = true(length(dan_regions(l).vertices_index),1);      
        dan_regions(l).type = "mask";
        dan_regions(l).net_lab = "dan";
        l = l+1;
end
end