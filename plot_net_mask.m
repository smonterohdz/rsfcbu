function [h] = plot_net_mask(mesh_brain,idx_select,regions)
%plot_net_mask
%   Detailed explanation goes here
A = zeros(1,size(mesh_brain.vertices,1));
A_select = zeros(1,length(idx_select));
col_range = [1 100];

if regions(1).type == "mask"
    for i=1:length(regions)
        if isfield(regions(i),'mask_subsetseed') && ~isempty(regions(i).mask_subsetseed)
            A_select(regions(i).vertices_index(regions(i).mask_subsetseed)) = 80;            
        else
            A_select(regions(i).vertices_index) = 80;
        end
    end
end


if regions(1).type == "seed"
    col_range = [1 100];

    if length(regions)==1
        if isfield(regions(1),'mask_subsetseed') && ~isempty(regions(1).mask_subsetseed)
            col_range = [1 max(regions.groups)+1];
            A_select(regions(1).vertices_index(regions(1).mask_subsetseed)) = regions.groups(regions(1).mask_subsetseed);
        else
            A_select(regions(1).vertices_index) = 1;
        end
    else
        a = 10;
        b = 90;
        for i=1:length(regions)
            idx_val = round(a + ((i-1)*(b-a))/(length(regions)-1));
            if isfield(regions(i),'mask_subsetseed') && ~isempty(regions(i).mask_subsetseed)
                A_select(regions(i).vertices_index(regions(i).mask_subsetseed)) = idx_val;
            else
                A_select(regions(i).vertices_index) = idx_val;

            end
        end
    end
end


if regions(1).type == "combmask"
    a = 40;
    b = 90;  

    idx_dmn_submask = find([regions(:).net_lab]=="dmn");
    idx_dan_submask = find([regions(:).net_lab]=="dan");
    dmn_nreg = length(idx_dmn_submask);   
    dan_nreg = length(idx_dan_submask);

    col_range = [-100 100];
    for i=1:dmn_nreg
        idx_val = round(a + ((i-1)*(b-a))/(dmn_nreg-1));
        if isfield(regions(idx_dmn_submask(i)),'mask_subsetseed') && ~isempty(regions(idx_dmn_submask(i)).mask_subsetseed)
            A_select(regions(idx_dmn_submask(i)).vertices_index(regions(idx_dmn_submask(i)).mask_subsetseed)) = idx_val;
        else
            A_select(regions(idx_dmn_submask(i)).vertices_index) = idx_val;
        end
    end

    for i=1:dan_nreg
        idx_val = -round(a + ((i-1)*(b-a))/(dan_nreg-1));
        if isfield(regions(idx_dan_submask(i)),'mask_subsetseed') && ~isempty(regions(idx_dan_submask(i)).mask_subsetseed)
            A_select(regions(idx_dan_submask(i)).vertices_index(regions(idx_dan_submask(i)).mask_subsetseed)) = idx_val;
        else
            A_select(regions(idx_dan_submask(i)).vertices_index) = idx_val;
        end
    end
end

if regions(1).type == "groupmask"
    a = 40;
    b = 90;

    idx_dmn_submask = find([regions(:).net_lab]=="dmn");
    idx_dan_submask = find([regions(:).net_lab]=="dan");
    dmn_nreg = length(idx_dmn_submask);
    dan_nreg = length(idx_dan_submask);

    col_range = [0 16];
    for i=1:dmn_nreg
        idx_val = round(a + ((i-1)*(b-a))/(dmn_nreg-1));
        if isfield(regions(idx_dmn_submask(i)),'mask_subsetseed') && ~isempty(regions(idx_dmn_submask(i)).mask_subsetseed)
            A_select(regions(idx_dmn_submask(i)).vertices_index) = regions(idx_dmn_submask(i)).mask_subsetseed;
        else
            A_select(regions(idx_dmn_submask(i)).vertices_index) = idx_val;
        end
    end

    for i=1:dan_nreg
        idx_val = -round(a + ((i-1)*(b-a))/(dan_nreg-1));
        if isfield(regions(idx_dan_submask(i)),'mask_subsetseed') && ~isempty(regions(idx_dan_submask(i)).mask_subsetseed)
            A_select(regions(idx_dan_submask(i)).vertices_index) = regions(idx_dan_submask(i)).mask_subsetseed;
        else
            A_select(regions(idx_dan_submask(i)).vertices_index) = idx_val;
        end
    end
end

A(idx_select) = A_select;
h = plot_connectivity_seedregion(mesh_brain,A,col_range,[],[],[],[]);
hold on;

end