function [mesh] =  plotjac_infant_mesh(mesh,PD_data,J,PD)

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
%so we are left with the total sensitivity for EACH node
mesh3.data = all_J; %set mesh3.data as the total sensitivty

%getting minimum and maximum J for just the cortex nodes. 
max_j = max(mesh3.data(1,cortex_nodes));
min_j = min(mesh3.data(1,cortex_nodes));

% figure()
% plotniceimages_1(mesh3,mesh_recon); %plotting the jacobian (requires function plotniceimages_1)
% caxis([min_j max_j]);
% colorbar('horiz');
% title("Non spatially Normalized J - ALL channels")

figure()
plotniceimages_1_Jacobian(mesh3,mesh_recon); %plotting the jacobian (requires function plotniceimages_1)
caxis([min_j max_j]);
colorbar('horiz');
title("Non spatially Normalized J - ALL channels")

%get threshold of top X % of sensitivty on cortex
%max_j
thresholdJ = min_j*0.95;
thresholdJ_idx= find(mesh3.data(1,cortex_nodes)>thresholdJ)';


% figure()
% %mesh3.data = sum(J_amp);
% plotniceimages_1(mesh3,mesh_recon);
% colorbar('horiz');

%Spatially normalize
mesh3.data = all_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
max_j = max(mesh3.data(1,cortex_nodes));
min_j = min(mesh3.data(1,cortex_nodes));

% figure()
% plotniceimages_1(mesh3,mesh_recon);
% title("Spatially Normalized J - ALL channels")
% caxis([min_j max_j]);
% colorbar('horiz');

figure()
plotniceimages_1_Jacobian(mesh3,mesh_recon);
title("Spatially Normalized J - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
caxis([min_j max_j]);
colorbar('horiz');

%creating mask.
all_j_sp_norm = all_J./mesh.support(:,1)';
mask_95_cortexnodes = find(all_j_sp_norm(cortex_nodes) < min(all_j_sp_norm(cortex_nodes))*0.05);
mask_99_cortexnodes = find(all_j_sp_norm(cortex_nodes) < min(all_j_sp_norm(cortex_nodes))*0.01);
mesh3.data=ones(size(all_J,2),1)*-1;
mesh3.data(cortex_nodes(mask_95_cortexnodes)) = 1;
figure()
subplot(1,2,1)
plotniceimages_1_greyJET(mesh3,mesh_recon);
title("Spatially Normalized J (95% Mask) - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
caxis([-1 1]);
colorbar('horiz');
subplot(1,2,2)
mesh3.data(cortex_nodes(mask_99_cortexnodes)) = 1;
plotniceimages_1_greyJET(mesh3,mesh_recon);
title("Spatially Normalized J (99% Mask) - All ch's - TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
caxis([-1 1]);
colorbar('horiz');

%JUST the GOOD IDX
 goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
 mesh3.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
% 
 max_j = max(mesh3.data(1,cortex_nodes));
 min_j = min(mesh3.data(1,cortex_nodes));

figure()
plotniceimages_1_Jacobian(mesh3,mesh_recon);
title("Spatially Normalized J - goodch idx only")
caxis([min_j max_j]);
colorbar('horiz');

%create mask
%mask = sum(abs(J))./max(sum(abs(J))) > 0.05
%find(mesh3.data(1,cortex_nodes) > min_j*0.95)  

%all_J %all J for all measurements for all nodes
%mask_thresh = all_J(1,cortex_nodes)
%mask = sum(abs(J))./max(sum(abs(J))) > 0.05 -index of cortex nodes
%mask = all_J(1,cortex_nodes) ./ min(all_J(1,cortex_nodes)) > 0.05;
mask =  all_J(1,:) < 0.05*min_j; %mask is all nodes less than 5% of the max J(cortex nodes)
%for all cortex nodes MASK = 1 if > 5% of max, MASK = 0 if < 5% of max 

%get mask, but for the spatially normalized Jacobian. - index of cortex
%nodes
%mask_spatial_norm =  (all_J(1,cortex_nodes)./mesh.support(cortex_nodes,1)') ./ min( (all_J(1,cortex_nodes)./mesh.support(cortex_nodes,1)')   ) > 0.05;
mask_spatial_norm =  (all_J(1,:)./mesh.support(:,1)') < 0.05*min_j; %mask is all nodes less than 5% of the max J(cortex nodes)

mesh.mask = mask;
mesh.mask_spatial_norm = mask;

mask_spatial_norm_zero = find(mask_spatial_norm == 0); 
mesh3.data(mask_spatial_norm_zero) = max_j;
%
figure()
plotniceimages_1_Jacobian(mesh3,mesh_recon);
title("Spatially Normalized J - goodch idx only - mask "+num2str(0.05*min_j)+" TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
caxis([min_j max_j]);
colorbar('horiz');


end