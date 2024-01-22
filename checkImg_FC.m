function checkImg_FC(fwFolder,HbX,xtime)
%checkImg_FC this functions helps us toto visualize the reconstructed image
%   Detailed explanation goes here

axes_order = [1,2,3];

% Load Brain Adot
load([fwFolder 'Adot.mat'],'Adot');%
index_select = log10(sum(Adot(:,:,1),1))>=-2;
idx_select = find(index_select);

% Load brain mesh
load([fwFolder, 'mesh_brain.mat'],'mesh');
f = mesh.faces;
v = mesh.vertices;
xImg = zeros(size(HbX,1),size(v,1));
xImg(:,idx_select) = HbX;

% Play a movie of the image time series
vmax = min(xImg,[],'all');
figure();
 colormap('jet');
% clim([vmax -vmax]);
% hold on;
for ii=find(xtime>60,1):floor(length(xtime))
    %ii=find(xtime>86.6,1);
    hf=trisurf( f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), xImg(ii,:), ...
        'facecolor','interp','edgecolor','flat','edgealpha',0, 'visible','on');
    campos([-2238.8, 132.0, 130.0])
    title( sprintf('Time = %0.1f',xtime(ii)) )
    axis image;  
    vmax = min(xImg(ii,:),[],'all');
    clim([vmax -vmax]);
    pause(0.005);
end  
end