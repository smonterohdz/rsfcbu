function h_fig = plot_connectivity_seedregion(mesh,A,color_bar_range,seed,...
    seed_HbO,seed_HbR,HbX_time,mask_idx)
f=mesh.faces;
v=mesh.vertices;
axes_order = [2,1,3];
%axes_order = [1,2,3];

h_fig = figure;
%h_fig = figure(gcf);
if isempty(seed_HbO) && isempty(seed_HbR)
    tiledlayout(3,10)
else
    tiledlayout(5,10)
end

if isempty(HbX_time)
    if ~isempty(seed_HbO)
        HbX_time = 1:length(seed_HbO);
    else
        HbX_time = 1:length(seed_HbR);
    end
end
A(isinf(A))= 0;
AZ = [90,90,-90,0,-90];
EL = [0,90,0,90,0];
mycampos=[128.0, 2238.8, 130.0;...
        287.7, 127.4, 1144.5;...%1469.2;...
        -2238.8, 132.0, 130.0;...
        -0287.7    0134.8   -1144.5;...
        128.0, -2291.8, 130.0];
mycamtarget=[128.0, 132.0, 130.0; ...
             128.0, 132.0, 130.0;...
             128.0, 132.0, 130.0;...
             128.0, 132.0, 130.0;...
             128.0, 132.0, 130.0];
mycamup=[-1.0, 0.0, 0.0;...
         -1.0, 0.0, 0.0;...
         0.0, 0.0, 1.0;...
         -1.0, 0.0, 0.0;...
         -1.0, 0.0, 0.0];
mytitles={'Left','Anterior','Superior','Posterior','Right'};

%campos([-0287.7    0134.8   -1144.5])

%Left
for i=1:5
    nexttile([3,2])
    hold on
    % h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
    %     A,'facecolor','interp','edgecolor','flat','edgealpha',0, 'visible','on');
    if ~isempty(seed)
        if size(seed,2)==3 %if columns are 3, then it contains coordinates. plot with plot3
            plot3(seed(:,axes_order(1)),seed(:,axes_order(2)),seed(:,axes_order(3)),'g.','MarkerSize',10)
        elseif size(seed',2)==1
            A(seed) = -21;
            if ~isempty(mask_idx)
                nonmask_idx = 1:length(A);
                nonmask_idx(mask_idx) =[];
                A(nonmask_idx) = 0;
            end
        end
    end
    h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
        A,'facecolor','interp','edgecolor','flat','edgealpha',0, 'visible','on');

    set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
    if ~isempty(color_bar_range)
        clim(color_bar_range);
        if isequal(color_bar_range,[0.9,1])
            myColorMap=cbrewer('seq','YlOrRd',9);
        elseif isequal(color_bar_range,[1 100])
            myColorMap=flipud(cbrewer('qual','Paired',20));
            myColorMap = [0.8,0.8,0.8;myColorMap];
        elseif isequal(color_bar_range,[0 100])
            myColorMap=cbrewer('seq','YlGn',25);
            %myColorMap = hsv(30);
            myColorMap = [0.7,0.7,0.7;myColorMap];
        elseif isequal(color_bar_range,[-100 100])            
            myColorMap=[flipud(cbrewer('seq','GnBu',length(unique(A))));...
                cbrewer('seq','YlOrRd',length(unique(A)))];
            myColorMap(length(unique(A)),:) = [0.95 0.95 0.95];
            myColorMap(length(unique(A))+1,:) = [0.95 0.95 0.95];
        else
            myColorMap=flipud(cbrewer('div','RdBu',256));
        end
        %myColorMap(113:144,:) = 0.8;
    else
        %vmax = min(abs([min(A),max(A)]));
        vmax = max(A);
        clim([-vmax, vmax]);
        if ~isempty(seed) && size(seed',2)==1
            myColorMap=[0,1,0;flipud(cbrewer('div','RdBu',256))];
        else
            myColorMap=flipud(cbrewer('div','RdBu',256));
        end
    end
    colormap(myColorMap);
   
    
    axis image
    axis off
    view(AZ(i),EL(i));
    camtarget(mycamtarget(i,:));
    campos(mycampos(i,:));
    camup(mycamup(i,:));
   
    %axis image
    %axis off
    % hold on
    % l = camlight;
    % set(l,'Position',[-2000 100 100]);
    % l2 = camlight;
    % set(l2,'Position',[50 -100 -100]);
    % camlight(0,0);
    % lighting phong;
    % light
    
    camlight('headlight');
    lighting gouraud;
    light
    subtitle(mytitles{i});
end
%% Old version
% %Anterior
% nexttile([3,2])
% hold on
% h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
%       A,'facecolor','interp','edgealpha',0, 'visible','on'); 
% if ~isempty(seed)
%     plot3(seed(:,axes_order(1)),seed(:,axes_order(2)),seed(:,axes_order(3)),'g.','LineWidth',10)
% end
% set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
% if ~isempty(color_bar_range)
%     clim(color_bar_range)
%     if isequal(color_bar_range,[0.9,1])
%         myColorMap=cbrewer('seq','YlOrRd',9);
%     elseif isequal(color_bar_range,[1 100])
%         myColorMap=flipud(cbrewer('qual','Paired',20));
%         myColorMap = [0.7,0.7,0.7;myColorMap];
%     else
%         myColorMap=flipud(cbrewer('div','RdBu',256));
%     end
%     %myColorMap(113:144,:) = 0.8;
% else
%     vmax = max(A);
%     clim([-vmax, vmax]);
%     %myColorMap = jet(256);
%     myColorMap=flipud(cbrewer('div','RdBu',256));
%     %myColorMap(123:134,:) = 0.8;
% end
% colormap(myColorMap);
% 
% % colorbar
% view(90,90)
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
% subtitle('Anterior');
% 
% 
% %Superior
% nexttile([3,2])
% hold on
% h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
%       A,'facecolor','interp','edgealpha',0, 'visible','on'); 
% if ~isempty(seed)
%     plot3(seed(:,axes_order(1)),seed(:,axes_order(2)),seed(:,axes_order(3)),'g.','LineWidth',10);
% end
% set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
% if ~isempty(color_bar_range)
%     clim(color_bar_range)
%     if isequal(color_bar_range,[0.9,1])
%         myColorMap=cbrewer('seq','YlOrRd',9);
%     elseif isequal(color_bar_range,[1 100])
%         myColorMap=flipud(cbrewer('qual','Paired',20));
%         myColorMap = [0.7,0.7,0.7;myColorMap];
%     else
%         myColorMap=flipud(cbrewer('div','RdBu',256));
%     end
%     %myColorMap(113:144,:) = 0.8;
% else
%     vmax = max(A);
%     clim([-vmax, vmax]);
%     %myColorMap = jet(256);
%     myColorMap=flipud(cbrewer('div','RdBu',256));
%     %myColorMap(123:134,:) = 0.8;
% end
% colormap(myColorMap);
% 
% % colorbar
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([-2238.8, 132.0, 130.0])
% % axis image
% % axis off
% % hold on
% % l = camlight;
% % set(l,'Position',[-2000 100 100]);
% % l2 = camlight;
% % set(l2,'Position',[50 -100 -100]);
% % lighting phong;
% % light
% axis image
% axis off
% camlight(0,0);
% lighting phong;
% light
% subtitle('Superior');
% 
% 
% % Posterior
% nexttile([3,2])
% hold on
% h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
%       A,'facecolor','interp','edgealpha',0, 'visible','on'); 
% if ~isempty(seed)
%     plot3(seed(:,axes_order(1)),seed(:,axes_order(2)),seed(:,axes_order(3)),'g.','LineWidth',10);
% end
% set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
% 
% if ~isempty(color_bar_range)
%     clim(color_bar_range)
%    if isequal(color_bar_range,[0.9,1])
%         myColorMap=cbrewer('seq','YlOrRd',9);
%     elseif isequal(color_bar_range,[1 100])
%         myColorMap=flipud(cbrewer('qual','Paired',20));
%         myColorMap = [0.7,0.7,0.7;myColorMap];
%     else
%         myColorMap=flipud(cbrewer('div','RdBu',256));
%     end
%     %myColorMap(113:144,:) = 0.8;
% else
%     vmax = max(A);
%     clim([-vmax, vmax]);
%     %myColorMap = jet(256);
%     myColorMap=flipud(cbrewer('div','RdBu',256));
%     %myColorMap(123:134,:) = 0.8;
% end
% colormap(myColorMap);
% 
% % colorbar
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([-0287.7    0134.8   -1144.5])
% camup([-1.0, 0.0, 0.0])
% % axis image
% % axis off
% % hold on
% % l = camlight;
% % set(l,'Position',[-2000 100 100]);
% % l2 = camlight;
% % set(l2,'Position',[50 -100 -100]);
% % camlight(0,0);
% % lighting phong;
% % light
% axis image
% axis off
% camlight(0,0);
% lighting phong;
% light
% subtitle('Posterior');
% 
% 
% %Right
% nexttile([3,2])
% hold on
% h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
%       A,'facecolor','interp','edgealpha',0, 'visible','on'); 
% if ~isempty(seed)
%     plot3(seed(:,axes_order(1)),seed(:,axes_order(2)),seed(:,axes_order(3)),'g.','LineWidth',10);
% end
% set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
% 
% if ~isempty(color_bar_range)
%     clim(color_bar_range)
%     if isequal(color_bar_range,[0.9,1])
%         myColorMap=cbrewer('seq','YlOrRd',9);
%     elseif isequal(color_bar_range,[1 100])
%         myColorMap=flipud(cbrewer('qual','Paired',20));
%         myColorMap = [0.7,0.7,0.7;myColorMap];
%     else
%         myColorMap=flipud(cbrewer('div','RdBu',256));
%     end
%     %myColorMap(113:144,:) = 0.8;
% else
%     vmax = max(A);
%     clim([-vmax, vmax]);
%     %myColorMap = jet(256);
%     myColorMap=flipud(cbrewer('div','RdBu',256));
%     %myColorMap(123:134,:) = 0.8;
% end
% colormap(myColorMap);
% 
% % colorbar
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([128.0, -2291.8, 130.0])
% camup([-1.0, 0.0, 0.0])
% % axis image
% % axis off
% % hold on
% % l = camlight;
% % set(l,'Position',[-2000 100 100]);
% % l2 = camlight;
% % set(l2,'Position',[50 -100 -100]);
% % camlight(0,0);
% % lighting phong;
% % light
% axis image
% axis off
% camlight(0,0);
% lighting phong;
% light
% subtitle('Right');



%set(gcf,'position', [100 100 1000 250])
set(gcf,'position', [100 100 1280 600])
colorbar();
if ~isempty(seed_HbO) && ~isempty(seed_HbR)
    nexttile([2,10])
    plot(HbX_time,seed_HbO,'r-'); hold on;
    plot(HbX_time,seed_HbR,'b-');
    axis tight;
    xlabel('Time (s)'); ylabel('Delta HbO-HbR');
    % [pxx,f] = periodogram(seed_HbO,hamming(length(seed_HbO)),length(seed_HbO),fs,'power');
    % [pwrest,idx] = max(pxx(f<fcut_max)); % FIX Make it age-dependent
    % nexttile([2,5]);
    % plot(f,10*log10(pxx));
    % xlim([0,0.5]);
    % hold on;
    % plot(f(idx),10*log10(pxx(idx)),'*','Marker','p','MarkerSize',12,'MarkerFaceColor','r');
    % xlabel('Frequency(Hz)'); ylabel('Power (dB)');
elseif ~isempty(seed_HbO)
    nexttile([2,10])
    plot(HbX_time,seed_HbO,'r-'); %hold on;
    %plot(HbX_time,seed_HbR,'b-');
    axis tight;
    xlabel('Time (s)'); ylabel('Delta HbO');
    % [pxx,f] = periodogram(seed_HbO,hamming(length(seed_HbO)),length(seed_HbO),fs,'power');
    % [pwrest,idx] = max(pxx(f<fcut_max)); % FIX Make it age-dependent
    % nexttile([2,5]);
    % plot(f,10*log10(pxx));
    % xlim([0,0.5]);
    % hold on;
    % plot(f(idx),10*log10(pxx(idx)),'*','Marker','p','MarkerSize',12,'MarkerFaceColor','r');
    % xlabel('Frequency(Hz)'); ylabel('Power (dB)');
elseif ~isempty(seed_HbR)
    nexttile([2,10])
    %plot(HbX_time,seed_HbO,'r-'); %hold on;
    plot(HbX_time,seed_HbR,'b-');
    axis tight;
    xlabel('Time (s)');ylabel('Delta HbR');
    % [pxx,f] = periodogram(seed_HbR,hamming(length(seed_HbR)),length(seed_HbR),fs,'power');
    % [pwrest,idx] = max(pxx(f<fcut_max)); % FIX Make it age-dependent
    % nexttile([2,5]);
    % plot(f,10*log10(pxx));
    % xlim([0,0.5]);
    % hold on;
    % plot(f(idx),10*log10(pxx(idx)),'*','Marker','p','MarkerSize',12,'MarkerFaceColor','r');
    % xlabel('Frequency(Hz)'); ylabel('Power (dB)');
else
    set(gcf,'position', [100 100 1280 310]);
end
end