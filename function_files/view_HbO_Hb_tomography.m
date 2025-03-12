function cortex_nodes = view_HbO_Hb_tomography(HbO_c,Hb_c,HbT_c,goodch_idx,weekN)

if weekN == 30;

    load('mesh_infant30week_850nm.mat')

end



%% View
%mesh_recon = mesh;
ind = reshape(mesh.region(mesh.elements),[],1);
ind = reshape(ind>=3,[],4);
ind = sum(ind,2);
ind = find(ind==4);
[mesh3.elements,mesh3.nodes]=boundfaces(mesh.nodes,mesh.elements(ind,:),0);
%[mesh_all.elements,mesh_all.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%[meshZ.elements,meshZ.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%cortex_nodes = find(mesh.region>=3);
cortex_nodes = find(mesh.region==3);

%max_dC = max(abs([min(min([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)])) max(max([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)]))]));
%min_dC = max(abs([min(min([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)])) max(max([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)]))]));

dC_min_max = [min(min([Hb_c(cortex_nodes,:) HbO_c(cortex_nodes,:) HbT_c(cortex_nodes,:)])) max(max([Hb_c(cortex_nodes,:) HbO_c(cortex_nodes,:) HbT_c(cortex_nodes,:)]))];
% HbO_CL_lhs = HbO_lhs; %comment if using DS. 
% Hb_CL_lhs = Hb_lhs; %comment if using DS.
% HbT_CL_lhs = HbT_lhs; %comment if using DS.

i=1;
i=22;
i=45;
i=75;
i=100;
i=105;

time_points = 0:15:size(HbO_c,2);
time_points(1)=1;

for i=1:size(time_points,2)
    mesh3.data = HbO_c(:,time_points(i));
    plotniceimages_1(mesh3,mesh);
    colorbar('horiz');
    max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
    min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
    %caxis([min_HbO_cortex max_HbO_cortex]);
    caxis([dC_min_max(1) dC_min_max(2)]);
end


HbO_c_cortex = sum(HbO_c(cortex_nodes,:));



end