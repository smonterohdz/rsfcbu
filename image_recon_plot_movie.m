function fhbx = image_recon_plot_movie(filename,mesh,index_select,Hb_brain,Hb_brain_np,t,fs)
% %main image recon
% clear all;
% close all;
%
% %
% usrname=getenv('USERNAME');
% baseDir = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching'];
% inputDir = [baseDir '\Rest_Movie_WorkingMemory\'];
% outputDir = [baseDir '\Rest_Movie_WorkingMemory\results'];
%
% %%
% data_folder = inputDir;
% load([inputDir, '\fw\Adot.mat']);%Brain
% load([inputDir, '\fw\mesh_brain.mat'])%MeshBrain
% index_select = log10(sum(Adot(:,:,1),1))>=-3;
jet256 = jet(256);
% %fs = 50;
% fs = 10.1725;
init_tp = ceil(t.start*fs)+1;
end_tp = floor(t.end*fs);
step_tp = 11;
% %filename = 'test.gif';
%
% Subject_folder = dir(fullfile(data_folder, 'sub-*'));
% Subject_folder = setdiff({Subject_folder([Subject_folder.isdir]).name},{'.','..'});
%
% %%
% for ii = 1:numel(Subject_folder)
%     % if ~strcmp(Subject_folder{ii},'sub-11')
%     %     continue
%     % end
%     fprintf('%s\n', Subject_folder{ii});
%     if ~(strcmp(Subject_folder{ii},'sub-11') || strcmp(Subject_folder{ii},'sub-14') ||...
%             strcmp(Subject_folder{ii},'sub-6') || strcmp(Subject_folder{ii},'sub-9') ||...
%             strcmp(Subject_folder{ii},'sub-10'))
%         continue;
%     end
%
%     Runs_RS = dir(fullfile(data_folder,Subject_folder{ii},'nirs', '*task-RS*.snirf'));
%     %Runs_MW = dir(fullfile(data_folder,Subject_folder{ii},'nirs', '*task-MW*.snirf'));
%     Runs_MW = dir(fullfile(data_folder,Subject_folder{ii},'*Finger_Tapping*.snirf'));
%
%     Runs = [Runs_RS; Runs_MW];
%     %Runs = Runs_RS;
%     Runs = {Runs(~[Runs.isdir]).name};
%     %savedir = fullfile('..', 'results', Subject_folder{ii},'conn');
%     savedir = fullfile(data_folder,'derivatives','homer', Subject_folder{ii},'nirs');
%
%     for chbtype = {'hbr','hbo'}
%         hbtype = chbtype{1};
%         for jj = 1:numel(Runs)
%             fprintf('%s\n', Runs{jj})
%
%             save_mat_path = savedir;
%             filename = Runs{jj};
%
%             if strcmp(hbtype,"hbo") == 1
%                 load(fullfile(save_mat_path, [filename(1:end-6),'_Hb_brain.mat']),"HbO_brain")
%                 Hb_brain = HbO_brain;
%                 clear HbO_brain;
%             else
%                 load(fullfile(save_mat_path, [filename(1:end-6),'_Hb_brain.mat']),"HbR_brain")
%                 Hb_brain = HbR_brain;
%                 clear HbR_brain;
%             end
Hbmax = max(abs(Hb_brain(:)));
col_range = [-Hbmax,Hbmax].*0.02;
%
f=mesh.faces;
v=mesh.vertices;
axes_order = [2,1,3];

n_vox_mesh = size(mesh.vertices,1);

Adot_new = zeros(1,n_vox_mesh);
fhbx = figure;
fhbx.Position = [512,90,524,400];
axis tight manual % this ensures that getframe() returns a consistent size

filenamevid = filename;%[filename '.vid'];

vid = VideoWriter([filenamevid],'MPEG-4');
vid.FrameRate =4;
open(vid);
for time_p = init_tp:step_tp:end_tp
    if time_p>size(Hb_brain,2)
        break;
    end
    Adot_new(1,index_select) = Hb_brain(time_p,:);
    %figure(fhbx);
    if ~isempty(Hb_brain_np)
        subplot(1,3,1);
    end
    h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
        Adot_new,'facecolor','interp','edgecolor','interp','edgealpha',0, 'visible','on');
    set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
    % caxis([-2 0])
    ax = gca;
    ax.CLim = col_range;
    colormap(jet256);
    colorbar
    view(-90,0)    
    camtarget([128.0, 132.0, 130.0])    
    axis image
    axis off
    camlight(0,0);
    lighting phong;
    light
    title(sprintf('HbO %.2fs',time_p/fs));

    % Non-processed
    if ~isempty(Hb_brain_np)
    Adot_new(1,index_select) = Hb_brain_np(time_p,:);
    subplot(1,3,2);
    h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
        Adot_new,'facecolor','interp','edgealpha',0, 'visible','on');
    set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
    % caxis([-2 0])
    ax = gca;
    ax.CLim = col_range;
    colormap(jet256);
    colorbar
    view(-90,0)    
    camtarget([128.0, 132.0, 130.0])    
    axis image
    axis off
    camlight(0,0);
    lighting phong;
    light
    subtitle('Non GSR');
    end

    % % view 2
    % subplot(1,3,2);
    % h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
    %     Adot_new,'facecolor','interp','edgealpha',0, 'visible','on');
    % set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
    % % caxis([-2 0])
    % ax = gca;
    % ax.CLim = col_range;
    % 
    % colormap(jet256);
    % %colorbar
    % view(90,90);
    % % camlight
    % camtarget([128.0, 132.0, 130.0])
    % % camlight
    % %campos([-2238.8, 132.0, 130.0])
    % % camlight
    % % camup([-1.0, 0.0, 0.0])
    % % camlight
    % axis image
    % axis off
    % 
    % %l = camlight;
    % %set(l,'Position',[-2000 100 100]);
    % 
    % %l2 = camlight;
    % %set(l2,'Position',[50 -100 -100]);
    % camlight(0,0);
    % 
    % lighting phong;
    % %light
    % subtitle('Anterior');
    % title(sprintf('%s (%s) %.2fs',filename(1:end-6),'HbO',time_p/fs),'Interpreter','none');
    % 
    % 
    % % view 3
    % subplot(1,3,3);
    % h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
    %     Adot_new,'facecolor','interp','edgealpha',0, 'visible','on');
    % set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
    % % caxis([-2 0])
    % ax = gca;
    % ax.CLim = col_range;
    % 
    % colormap(jet256);
    % view(-90,-90);
    % % camlight
    % camtarget([128.0, 132.0, 130.0])
    % % camlight
    % %campos([-2238.8, 132.0, 130.0])
    % % camlight
    % % camup([-1.0, 0.0, 0.0])
    % % camlight
    % axis image
    % axis off
    % 
    % %l = camlight;
    % %set(l,'Position',[-2000 100 100]);
    % 
    % %l2 = camlight;
    % %set(l2,'Position',[50 -100 -100]);
    % camlight(0,0);
    % 
    % lighting phong;
    % light
    % subtitle('Posterior');
    % cb = colorbar;
    % cb.Position = [0.92 0.25 0.02 0.5];

    drawnow;
    frame = getframe(fhbx);
    writeVideo(vid,frame);

    % im = frame2im(frame);
    % [imind,cm] = rgb2ind(im,256);
    % % Write to the GIF File
    % if time_p == init_tp
    %     imwrite(imind,cm,filenamegif,'gif', 'Loopcount',inf);
    % else
    %     imwrite(imind,cm,filenamegif,'gif','WriteMode','append');
    % end
end
close(vid);
close(fhbx);
%         end
%     end
% end
end