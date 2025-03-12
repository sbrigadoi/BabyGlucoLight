function plotniceimagesleep(mesh,mesh_recon,mm)
figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5])
subplot(1,2,1)
trisurf(mesh.elements,mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),full(mesh.data));
%caxis([0 mm]);
%ind = find(mesh.region>=4);
%plot3(mesh.nodes(ind,1),mesh.nodes(ind,2),mesh.nodes(ind,3),'c.');
axis equal;
view(-90,0);
colormap hot;
axis off
set(gca,'Color',[1 1 1]);
set(gcf,'Color',[1 1 1]);
alpha(0.45);
hold on
tmp = sort(mesh_recon.source.num);
if length(tmp)>1
    s1 = mesh_recon.source.coord(mesh_recon.source.num == tmp(1),:);
    s2 = mesh_recon.source.coord(mesh_recon.source.num == tmp(2),:);
    plot3(s1(:,1),s1(:,2),s1(:,3),'go',s2(:,1),s2(:,2),s2(:,3),'yo',...
        mesh_recon.source.coord(3:end,1),...
        mesh_recon.source.coord(3:end,2),...
        mesh_recon.source.coord(3:end,3),'ro','LineWidth',2,'MarkerSize',8);
else
    s1 = mesh_recon.source.coord;
    plot3(s1(:,1),s1(:,2),s1(:,3),'go','LineWidth',2,'MarkerSize',8);
end

if isfield(mesh_recon,'meas') == 1
    plot3(mesh_recon.meas.coord(:,1),...
        mesh_recon.meas.coord(:,2),...
        mesh_recon.meas.coord(:,3),'bx',...
        'LineWidth',2,'MarkerSize',8);
end


ind = reshape(mesh.data(mesh_recon.elements),[],1);
ind = reshape(ind>0,[],4);
ind = sum(ind,2);
ind = find(ind>0);
[elements,nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
hold on
%trisurf(elements,nodes(:,1),nodes(:,2),nodes(:,3),full(mesh.data),'FaceColor','r');
subplot(1,2,2)
trisurf(mesh.elements,mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),full(mesh.data));
%caxis([0 mm]);
%ind = find(mesh.region>=4);
%plot3(mesh.nodes(ind,1),mesh.nodes(ind,2),mesh.nodes(ind,3),'c.');
axis equal;
view(-180,0);
colormap hot;
axis off
set(gca,'Color',[1 1 1]);
set(gcf,'Color',[1 1 1]);
alpha(0.45)

ind = reshape(mesh.data(mesh_recon.elements),[],1);
ind = reshape(ind>0,[],4);
ind = sum(ind,2);
ind = find(ind>0);
[elements,nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
hold on
%trisurf(elements,nodes(:,1),nodes(:,2),nodes(:,3),-full(mesh.data),'FaceColor','r');

tmp = sort(mesh_recon.source.num);
if length(tmp)>1
    s1 = mesh_recon.source.coord(mesh_recon.source.num == tmp(1),:);
    s2 = mesh_recon.source.coord(mesh_recon.source.num == tmp(2),:);
    plot3(s1(:,1),s1(:,2),s1(:,3),'go',s2(:,1),s2(:,2),s2(:,3),'yo',...
        mesh_recon.source.coord(3:end,1),...
        mesh_recon.source.coord(3:end,2),...
        mesh_recon.source.coord(3:end,3),'ro','LineWidth',2,'MarkerSize',8);
else
    s1 = mesh_recon.source.coord;
    plot3(s1(:,1),s1(:,2),s1(:,3),'go','LineWidth',2,'MarkerSize',8);
end

if isfield(mesh_recon,'meas') == 1
    plot3(mesh_recon.meas.coord(:,1),...
        mesh_recon.meas.coord(:,2),...
        mesh_recon.meas.coord(:,3),'bx',...
        'LineWidth',2,'MarkerSize',8);
end