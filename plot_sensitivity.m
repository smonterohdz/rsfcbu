function plot_sensitivity(mesh,A, color_bar_range)%,fhbx)
f=mesh.faces;
v=mesh.vertices;
%figure(fhbx);

axes_order = [2,1,3];

h = trisurf(f, v(:,axes_order(1)), v(:,axes_order(2)), v(:,axes_order(3)), ...
      A,'facecolor','interp','edgealpha',0, 'visible','on'); 
set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
% caxis([-2 0])
ax = gca;
ax.CLim = color_bar_range;
myColorMap = jet(256);
%myColorMap(1:2,:) = 0.8;
colormap(myColorMap);
colorbar
view(0,0);
% camlight
camtarget([128.0, 132.0, 130.0])
% camlight
%campos([-2238.8, 132.0, 130.0])
% camlight
% camup([-1.0, 0.0, 0.0])
% camlight
axis image
axis on
box on
hold on
xlabel('X');
ylabel('Y');
zlabel('Z');
l = camlight;
set(l,'Position',[-2000 100 100]);

l2 = camlight;
set(l2,'Position',[50 -100 -100]);

camlight(0,0);

lighting phong;
light
end