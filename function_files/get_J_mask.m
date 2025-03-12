function [masks] = get_J_mask(mesh,J,PD,PD_data,vol2gm,gmSurfaceMesh,coverageThresh)

J_amp = J.complete; %amp - We just need this part of the Jacobian as it describes the intensity
%EACH ROW IS A MEASUREMENT CORROSPONDING TO THE MESH.LINK
%EACH COLUMN IS A NODE CORROSPONDING TO THE MESH.NODE

% % mua Amplitude only
%J = [J_amp(goodch_idx(1:end/2,1),:)]; %amp, mua only %goodch idx is symmetrical.

mesh_recon = mesh;
ind = reshape(mesh_recon.region(mesh_recon.elements),[],1);
%ind = reshape(ind>=3,[],4);
ind = reshape(ind==3,[],4); %just GM nodes 06 05 24
ind = sum(ind,2);
ind = find(ind==4);
[mesh3.elements,mesh3.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%[mesh_all.elements,mesh_all.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%[meshZ.elements,meshZ.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%cortex_nodes = find(mesh.region>=3);
cortex_nodes = find(mesh.region==3);




all_J = sum(J_amp(:,:)); %this sums the jacobian for all measurements and nodes


%Spatially normalize
mesh3.data = all_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
max_j = max(mesh3.data(1,cortex_nodes));
min_j = min(mesh3.data(1,cortex_nodes));

% figure()
% plotniceimages_1(mesh3,mesh_recon);
% title("Spatially Normalized J - ALL channels")
% caxis([min_j max_j]);
% colorbar('horiz');

% figure()
% plotniceimages_1_Jacobian(mesh3,mesh_recon);
% title("Spatially Normalized J - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
% caxis([min_j max_j]);
% colorbar('horiz');


%creating mask.
all_j_sp_norm = all_J./mesh.support(:,1)';
mask_95_cortexnodes_allch = find(all_j_sp_norm(cortex_nodes) < min(all_j_sp_norm(cortex_nodes))*0.05);
mask_nodes_bottom_5_cortexnodes_allch = find(all_j_sp_norm(cortex_nodes) > min(all_j_sp_norm(cortex_nodes))*0.05);
mask_99_cortexnodes_allch = find(all_j_sp_norm(cortex_nodes) < min(all_j_sp_norm(cortex_nodes))*0.01);
mask_nodes_bottom_1_cortexnodes_allch = find(all_j_sp_norm(cortex_nodes) > min(all_j_sp_norm(cortex_nodes))*0.01);
mesh3.data=ones(size(all_J,2),1)*-1;
mesh3.data(cortex_nodes(mask_95_cortexnodes_allch)) = 1;

masks.mask_95_cortexnodes_allch = mask_95_cortexnodes_allch;
masks.mask_99_cortexnodes_allch = mask_99_cortexnodes_allch;

masks.mask_nodes_bottom_5_cortexnodes_allch = mask_nodes_bottom_5_cortexnodes_allch;
masks.mask_nodes_bottom_1_cortexnodes_allch = mask_nodes_bottom_1_cortexnodes_allch;


% figure()
% subplot(1,2,1)
% plotniceimages_1_greyJET(mesh3,mesh_recon);
% title("Spatially Normalized J (95% Mask) - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
% caxis([-1 1]);
% colorbar('horiz');
% subplot(1,2,2)
% mesh3.data(cortex_nodes(mask_99_cortexnodes_allch)) = 1;
% plotniceimages_1_greyJET(mesh3,mesh_recon);
% title("Spatially Normalized J (99% Mask) - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
% caxis([-1 1]);
% colorbar('horiz');

%% use gmsurface mesh
%map J.complete (Intensity J) to gmsurface
J_gmsurface = vol2gm*J_amp';
%sum this across all 64 channels
J_all_gmsurface = sum(J_gmsurface,2);
%map meshsupport to gmsurface
mesh_support_gmsurface = vol2gm*mesh.support(:,1);
% spatially normalise J all gmsurface 
J_all_gmsurface_norm=J_all_gmsurface./mesh_support_gmsurface;
%plot J all gmsurface
figure()
plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface],gmSurfaceMesh.face);view(0,90);
cb = colorbar('horiz');
title("J - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm");

%plot J all gmsurface spatially normalized
figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface_norm],gmSurfaceMesh.face);view(0,90);
cb = colorbar('horiz');
title("Spatially Normalized J - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm");

%creating mask.
%all_j_sp_norm = all_J./mesh.support(:,1)';

%J_all_gmsurface_norm
mask_95_cortexnodes_allch = find(J_all_gmsurface_norm < min(J_all_gmsurface_norm)*0.05);
mask_nodes_bottom_5_cortexnodes_allch = find(J_all_gmsurface_norm > min(J_all_gmsurface_norm)*0.05);
mask_99_cortexnodes_allch = find(J_all_gmsurface_norm < min(J_all_gmsurface_norm)*0.01);
mask_nodes_bottom_1_cortexnodes_allch = find(J_all_gmsurface_norm > min(J_all_gmsurface_norm)*0.01);

mask_j_thresh_cortexnodes_allch = find(J_all_gmsurface < -coverageThresh);
mask_nodes_bottom_j_thresh_cortexnodes_allch = find(J_all_gmsurface > -coverageThresh);

masks.mask_95_cortexnodes_gmsurface_allch = mask_95_cortexnodes_allch;
masks.mask_99_cortexnodes_gmsurface_allch = mask_99_cortexnodes_allch;
masks.mask_nodes_bottom_5_cortexnodes_gmsurface_allch = mask_nodes_bottom_5_cortexnodes_allch;
masks.mask_nodes_bottom_1_cortexnodes_gmsurface_allch = mask_nodes_bottom_1_cortexnodes_allch;

masks.mask_j_thresh_cortexnodes_allch = mask_j_thresh_cortexnodes_allch;
masks.mask_nodes_bottom_j_thresh_cortexnodes_allch = mask_nodes_bottom_j_thresh_cortexnodes_allch;

masks.mask_j_thresh_gmsurface = ones(size(J_all_gmsurface_norm,1),1)*-1;
masks.mask_j_thresh_gmsurface(masks.mask_j_thresh_cortexnodes_allch)=1;


masks.mask_95_gmsurface = ones(size(J_all_gmsurface_norm,1),1)*-1;
masks.mask_95_gmsurface(masks.mask_95_cortexnodes_gmsurface_allch)=1;


figure; plotmesh_JAC([gmSurfaceMesh.node masks.mask_95_gmsurface],gmSurfaceMesh.face);view(0,90);
cb = colorbar('horiz');
title("MASK 95% Sp Normalized J - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm");

figure; plotmesh_JAC([gmSurfaceMesh.node masks.mask_j_thresh_gmsurface],gmSurfaceMesh.face);view(0,90);
cb = colorbar('horiz');
title("MASK Threshold J - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm");



end