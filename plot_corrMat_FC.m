function [f1,f2] = plot_corrMat_FC(rDMNDAN_hbo,rDMNDAN_hbr,subjects_set)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f1 = figure();
colormap(jet(64));
t = tiledlayout('flow','TileSpacing','compact');

for iSubj=1:size(rDMNDAN_hbo,3)
    subject = num2str(subjects_set(iSubj));
    nexttile;
    imagesc(rDMNDAN_hbo(:,:,iSubj), [-1 1] );
    %ylabel({'HbO';'Submasks Run1'});
    title(['Subj ',subject]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
end
sgtitle('HbO Subject Correlation matrices');
f1.Position = [261   52  850  727];
% t.Children(1).XLabel.String='DMN,DAN R1';
% t.Children(2).XLabel.String='DMN,DAN R1';
% t.Children(3).XLabel.String='DMN,DAN R1';
% 
% t.Children(3).YLabel.String='DMN,DAN R2';
% t.Children(end).YLabel.String='DMN,DAN R2';
% t.Children(end-4).YLabel.String='DMN,DAN R2';
% t.Children(end-8).YLabel.String='DMN,DAN R2';
cb=colorbar();
cb.Layout.Tile='east';
cb.Label.String='R-value';


f2 = figure();
colormap(jet(64));
t = tiledlayout('flow','TileSpacing','compact');
for iSubj=1:size(rDMNDAN_hbr,3)
    subject = num2str(subjects_set(iSubj));
    nexttile;
    imagesc(rDMNDAN_hbr(:,:,iSubj), [-1 1] );
    %ylabel({'HbO';'Submasks Run1'});
    title(['Subj ',subject]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
end
sgtitle('HbR Subject Correlation matrices');

f2.Position = [300   52  850  727];
% t.Children(1).XLabel.String='DMN,DAN R1';
% t.Children(2).XLabel.String='DMN,DAN R1';
% t.Children(3).XLabel.String='DMN,DAN R1';
% 
% t.Children(3).YLabel.String='DMN,DAN R2';
% t.Children(end).YLabel.String='DMN,DAN R2';
% t.Children(end-4).YLabel.String='DMN,DAN R2';
% t.Children(end-8).YLabel.String='DMN,DAN R2';
cb=colorbar();
cb.Layout.Tile='east';
cb.Label.String='R-value';
end