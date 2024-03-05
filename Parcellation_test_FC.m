function [dmn_regions,dan_regions,mesh_brain,idx_select,dmn_zero_v,dan_zero_v] = Parcellation_test_FC(anatomFolder,fwFolder,flags)
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
% figure;
% plot_sensitivity(mesh_brain,A, [-2 2])



%% reading schaefer parcels from mat file
load([fwFolder,'schaeferParcels.mat'],'rhParcelsDMN','lhParcelsDMN','rhParcelsDAN','lhParcelsDAN');
nParcelsDMNr = size(rhParcelsDMN,1);
nParcelsDMNl = size(lhParcelsDMN,1);
nParcelsDANr = size(rhParcelsDAN,1);
nParcelsDANl = size(lhParcelsDAN,1);
%read left and right hemispheres and concatenate the vertices and faces.
% load schaefer parcels
[vl,fl]=read_surf([fwFolder,'segment',filesep,'lh.pial']);
[vr,fr]=read_surf([fwFolder,'segment',filesep,'rh.pial']);
v = [vr;vl];
f = [fr;(fl+size(vr,1))];

dmn_regions = struct('name',[],'vertices',[],'vertices_index',[]);
brain_vertices_select = mesh_brain.vertices(index_select,:);
dan_regions = struct('name',[],'vertices',[],'vertices_index',[]);

zeros_r = zeros(length(vr),1);
zeros_l = zeros(length(vl),1);
k = 1;
l = 1;
dmn_zero_v = [];
dan_zero_v = [];

% %plot the brain
% figure
% trisurf(f,v(:,1),v(:,2),v(:,3),'edgecolor','flat');

if strcmp(flags.parcel_scheme,'schaefer')

    for iParcel =1:nParcelsDMNr
        dmn_regions(k).name = rhParcelsDMN{iParcel,1};
        current_mask = [rhParcelsDMN{iParcel,2};zeros_l];
        dmn_regions(k).vertices_index = find(current_mask(index_select));
        if isempty(dmn_regions(k).vertices_index), dmn_zero_v = [dmn_zero_v;k]; end
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
        if isempty(dmn_regions(k).vertices_index), dmn_zero_v = [dmn_zero_v;k]; end
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
        if isempty(dan_regions(l).vertices_index), dan_zero_v = [dan_zero_v;l]; end
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
        if isempty(dan_regions(l).vertices_index), dan_zero_v = [dan_zero_v;l]; end
        dan_regions(l).vertices = brain_vertices_select(dan_regions(l).vertices_index,:);
        dan_regions(l).mask_subsetseed = true(length(dan_regions(l).vertices_index),1);
        dan_regions(l).type = "mask";
        dan_regions(l).net_lab = "dan";
        l = l+1;
    end
elseif strcmp(flags.parcel_scheme,'schaefer_comb')
    % getting unique labels RH dmn 
    rhParcelsDMN_rep = cell(nParcelsDMNr,1);
    for iParcel =1:nParcelsDMNr
        pieces = strsplit(rhParcelsDMN{iParcel,1},'_');
        rhParcelsDMN_rep(iParcel) = join(pieces(1:end-1),'_');
    end
    rhParcelsDMN_unique = unique(rhParcelsDMN_rep);
    nParcelsDMNr_u = length(rhParcelsDMN_unique);
    rhParcelsDMN_id = zeros(nParcelsDMNr,1);
    for iParcel=1:nParcelsDMNr_u
        idx = strcmp(rhParcelsDMN_rep,rhParcelsDMN_unique{iParcel});
        rhParcelsDMN_id(idx) = iParcel;
    end
    % getting unique labels LH dmn
    lhParcelsDMN_rep = cell(nParcelsDMNl,1);
    for iParcel =1:nParcelsDMNl
        pieces = strsplit(lhParcelsDMN{iParcel,1},'_');
        lhParcelsDMN_rep(iParcel) = join(pieces(1:end-1),'_');
    end
    lhParcelsDMN_unique = unique(lhParcelsDMN_rep);
    nParcelsDMNl_u = length(lhParcelsDMN_unique);
    lhParcelsDMN_id = zeros(nParcelsDMNl,1);
    for iParcel=1:nParcelsDMNl_u
        idx = strcmp(lhParcelsDMN_rep,lhParcelsDMN_unique{iParcel});
        lhParcelsDMN_id(idx) = iParcel;
    end

    % creating my region structure rh dmn
    for iParcel =1:nParcelsDMNr_u
        iParcel_idx = rhParcelsDMN_id==iParcel;
        dmn_regions(k).name = rhParcelsDMN_unique{iParcel};
        current_mask = [any(logical([rhParcelsDMN{iParcel_idx,2}]),2);zeros_l];
        dmn_regions(k).vertices_index = find(current_mask(index_select));
        if isempty(dmn_regions(k).vertices_index)
            dmn_zero_v = [dmn_zero_v;k];
        end
        dmn_regions(k).vertices = brain_vertices_select(dmn_regions(k).vertices_index,:);
        dmn_regions(k).mask_subsetseed = true(length(dmn_regions(k).vertices_index),1);
        dmn_regions(k).type = "mask";
        dmn_regions(k).net_lab = "dmn";
        k = k+1;
    end

    % creating my region structure lh dmn
    for iParcel =1:nParcelsDMNl_u
        iParcel_idx = lhParcelsDMN_id==iParcel;
        dmn_regions(k).name = lhParcelsDMN_unique{iParcel};
        current_mask = [zeros_r;any(logical([lhParcelsDMN{iParcel_idx,2}]),2)];
        dmn_regions(k).vertices_index = find(current_mask(index_select));
        if isempty(dmn_regions(k).vertices_index)
            dmn_zero_v = [dmn_zero_v;k];
        end
        dmn_regions(k).vertices = brain_vertices_select(dmn_regions(k).vertices_index,:);
        dmn_regions(k).mask_subsetseed = true(length(dmn_regions(k).vertices_index),1);
        dmn_regions(k).type = "mask";
        dmn_regions(k).net_lab = "dmn";
        k = k+1;
    end
    %-----------
    % getting unique labels RH dan
    rhParcelsDAN_rep = cell(nParcelsDANr,1);
    for iParcel =1:nParcelsDANr
        pieces = strsplit(rhParcelsDAN{iParcel,1},'_');
        rhParcelsDAN_rep(iParcel) = join(pieces(1:end-1),'_');
    end
    rhParcelsDAN_unique = unique(rhParcelsDAN_rep);
    nParcelsDANr_u = length(rhParcelsDAN_unique);
    rhParcelsDAN_id = zeros(nParcelsDANr,1);
    for iParcel=1:nParcelsDANr_u
        idx = strcmp(rhParcelsDAN_rep,rhParcelsDAN_unique{iParcel});
        rhParcelsDAN_id(idx) = iParcel;
    end
    %  getting unique labels LH dan
    lhParcelsDAN_rep = cell(nParcelsDANl,1);
    for iParcel =1:nParcelsDANl
        pieces = strsplit(lhParcelsDAN{iParcel,1},'_');
        lhParcelsDAN_rep(iParcel) = join(pieces(1:end-1),'_');
    end
    lhParcelsDAN_unique = unique(lhParcelsDAN_rep);
    nParcelsDANl_u = length(lhParcelsDAN_unique);
    lhParcelsDAN_id = zeros(nParcelsDANl,1);
    for iParcel=1:nParcelsDANl_u
        idx = strcmp(lhParcelsDAN_rep,lhParcelsDAN_unique{iParcel});
        lhParcelsDAN_id(idx) = iParcel;
    end

    % creating my region structure rh dan
    for iParcel =1:nParcelsDANr_u
        iParcel_idx = rhParcelsDAN_id==iParcel;
        dan_regions(l).name = rhParcelsDAN_unique{iParcel};
        current_mask = [any(logical([rhParcelsDAN{iParcel_idx,2}]),2);zeros_l];
        dan_regions(l).vertices_index = find(current_mask(index_select));
        if isempty(dan_regions(l).vertices_index)
            dan_zero_v = [dan_zero_v;l];
        end
        dan_regions(l).vertices = brain_vertices_select(dan_regions(l).vertices_index,:);
        dan_regions(l).mask_subsetseed = true(length(dan_regions(l).vertices_index),1);
        dan_regions(l).type = "mask";
        dan_regions(l).net_lab = "dan";
        l = l+1;
    end

    % creating my region structure lh dan
    for iParcel =1:nParcelsDANl_u
        iParcel_idx = lhParcelsDAN_id==iParcel;
        dan_regions(l).name = lhParcelsDAN_unique{iParcel};
        current_mask = [zeros_r;any(logical([lhParcelsDAN{iParcel_idx,2}]),2)];
        dan_regions(l).vertices_index = find(current_mask(index_select));
        if isempty(dan_regions(l).vertices_index)
            dan_zero_v = [dan_zero_v;l];
        end
        dan_regions(l).vertices = brain_vertices_select(dan_regions(l).vertices_index,:);
        dan_regions(l).mask_subsetseed = true(length(dan_regions(l).vertices_index),1);
        dan_regions(l).type = "mask";
        dan_regions(l).net_lab = "dan";
        l = l+1;
    end
end
end