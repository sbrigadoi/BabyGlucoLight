function [PD_data] = tomography_PD_data_10mins(mesh,gmSurfaceMesh,J,PD_data,PD_time,PD,savefig,vol2gm,coverageThresh)

%last arg 0 = not save figure 1 = save fig

[masks] = get_J_mask(mesh,J,PD,PD_data,vol2gm,gmSurfaceMesh,coverageThresh);
%%
J_amp = J.complete; %amp - We just need this part of the Jacobian as it describes the intensity


% use gmsurface mesh
%map J.complete (Intensity J) to gmsurface
J_gmsurface = vol2gm*J_amp';
%sum this across all 64 channels
J_all_gmsurface = sum(J_gmsurface,2);
%map meshsupport to gmsurface
mesh_support_gmsurface = vol2gm*mesh.support(:,1);
% spatially normalise J all gmsurface 
J_all_gmsurface_norm=J_all_gmsurface./mesh_support_gmsurface;

%JUST the GOOD IDX
goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
goodch_J_gmsurface = vol2gm*goodch_J';

max_j = max(goodch_J_gmsurface);
min_j = min(goodch_J_gmsurface);

%meshJ = mesh3;
%meshJ.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty

%map J.complete (Intensity J) to gmsurface
%J_gmsurface = gmSurfaceMesh.vol2gm*J.complete';
%J_gmsurface = J_gmsurface';


%EACH ROW IS A MEASUREMENT CORROSPONDING TO THE MESH.LINK
%EACH COLUMN IS A NODE CORROSPONDING TO THE MESH.NODE
%need to change remCh, to index the row number. e.g 1 2 3 4 7 11 ... not
%just 1 1 1 0 0 0 1 1 etc. 
%Ylhs_wv1 = [PD_data.dod1(:,goodch_idx(1:end/2,1))]'; %780 22 03 24
%Ylhs_wv2 = [PD_data.dod1(:,goodch_idx((end/2)+1:end,1))]'; %850 22 03 24

% %map J.complete (Intensity J) to gmsurface
% J_gmsurface = vol2gm*J.complete';
% %sum this across all 64 channels
% J_all_gmsurface = sum(J_gmsurface,2);
% %map meshsupport to gmsurface
% mesh_support_gmsurface = vol2gm*mesh.support(:,1);
% % spatially normalise J all gmsurface 
% J_all_gmsurface_norm=J_all_gmsurface./mesh_support_gmsurface;
% %plot J all gmsurface
% figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface],gmSurfaceMesh.face);view(0,90);
% cb = colorbar('horiz');
% %plot J all gmsurface spatially normalized
% figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface_norm],gmSurfaceMesh.face);view(0,90);
% cb = colorbar('horiz');

%time_points = [find(PD_data.t==779.9) find(PD_data.t==779.9+(60*2)) find(PD_data.t==779.9+(60*4)) find(PD_data.t==779.9+(60*6)) find(PD_data.t==779.9+(60*8)) find(PD_data.t==779.9+(60*10)) find(PD_data.t==779.9+(60*12))];
%time_points = 1:120:601;
%time_points = 1:12:601;
time_points = 1:1:601;

switch PD_data.time_window
    case "All"
        time_points = 1:10:size(PD_data.t,1);
end

Ylhs_wv1 = [PD_data.dod1(time_points,PD_data.goodch_idx(1:end/2,1))]'; %780
Ylhs_wv2 = [PD_data.dod1(time_points,PD_data.goodch_idx((end/2)+1:end,1))]'; %850

% % mua Amplitude only
%full mesh
J = [J_amp(PD_data.goodch_idx(1:end/2,1),:)]; %amp, mua only %goodch idx is symmetrical.
%GM surface mesh (smaller)
%J = [J_gmsurface(PD_data.goodch_idx(1:end/2,1),:)]; %amp, mua only %goodch idx is symmetrical.

wavelengths = PD_data.SD.Lambda; % wavelengths of the system
nWavs = length(PD_data.SD.Lambda); % n of wavelengths
%%
Eall = [];
for i = 1:nWavs
    Etmp = GetExtinctions(wavelengths(i));
    Etmp = Etmp(1:2); %HbO and HbR only
    Eall = [Eall; Etmp./1e7]; %This will be nWavs x 2;
end
E=Eall;

%For re running this section, set J. mua Amplitude only
J = [J_amp(PD_data.goodch_idx(1:end/2,1),:)]; %amp, mua only %goodch idx is symmetrical.
%J = [J_gmsurface(PD_data.goodch_idx(1:end/2,1),:)]; %amp, mua only %goodch idx is symmetrical.

% normalise for voxels
JTJ = sum(J(:,1:end/2).^2,1);
L1 = sqrt(JTJ + (1e-2*max(JTJ)));
JTJ = sum(J(:,end/2+1:end).^2,1);
L2 = sqrt(JTJ + (1e-2*max(JTJ)));
%Can set to 1 or normalize
L = [L1 L2]; 
%L(:) = 1; %20 03 24

J = bsxfun(@rdivide,J,L);
% normalise for data magnitude
JJT = sum(J.^2,2);
%Can set to 1 or normalize
M = sqrt(JJT + (1e-2*max(JJT)));
%M(:) = 1; %20 03 24
J = bsxfun(@rdivide,J,M);

lambda = 0.5%0.5%E-1;
% get update
Hess = J*J';
reg = eye(size(Hess)).*lambda.*sqrt(norm(Hess)); %can set to 1 or use regularisation
%higher reg = more blur , lower reg = high contrast but may overfit to
%noise
%
tmp = J'*((Hess+reg)\(diag(1./M)*Ylhs_wv1));
tmp = tmp./L';

%mus1_lhs = tmp(1:end/2,:); %for using mua and mus' in Jac
%mua1_lhs = tmp(end/2+1:end,:); %for using mua and mus' in Jac

mua1_lhs = tmp(:,:); %for just using mu a in Jac

tmp = J'*((Hess+reg)\(diag(1./M)*Ylhs_wv2));
tmp = tmp./L'; %Not included before, just added 08:26 30 06 22

%mus2_lhs = tmp(1:end/2,:);%for using mua and mus' in Jac
%mua2_lhs = tmp(end/2+1:end,:);%for using mua and mus' in Jac

mua2_lhs = tmp(:,:); %for just using mu a in Jac
%
% unmix chromophores

%%
clear all_J dodConv dodConv_og goodch_J HbO_lhs Hb_lhs HbT_lhs J JJT JTJ L L1 L2 mesh_recon mesh3
%full mesh - all time points
%HbO_lhs = zeros(size(mesh.nodes,1),size(mua1_lhs,2));
%Hb_lhs = zeros(size(mesh.nodes,1),size(mua1_lhs,2));
%HbT_lhs = zeros(size(mesh.nodes,1),size(mua1_lhs,2));

%creat full mesh size then convert to gmsurface mesh smooth - single time point method
HbO_lhs = zeros(size(mesh.nodes,1),1);
Hb_lhs = zeros(size(mesh.nodes,1),1);
HbT_lhs = zeros(size(mesh.nodes,1),1);

HbO_lhs_gmsurface = zeros(size(vol2gm,1),size(mua1_lhs,2));
Hb_lhs_gmsurface = zeros(size(vol2gm,1),size(mua1_lhs,2));
HbT_lhs_gmsurface = zeros(size(vol2gm,1),size(mua1_lhs,2));

%HbT gmsurface using vol2gm transformation
%HbT_all_m_hypo_allNodes = zeros(size(mesh.nodes,1),1);
%HbT_all_m_hypo_allNodes(cortex_nodes) = HbT_all_m_hypo(:,1,1);
%HbT_all_m_hypo_gmsurface= vol2gm*HbT_all_m_hypo_allNodes;

%map J.complete (Intensity J) to gmsurface
%J_gmsurface = vol2gm*J.complete';
%sum this across all 64 channels
%J_all_gmsurface = sum(J_gmsurface,2);
%map meshsupport to gmsurface
%mesh_support_gmsurface = vol2gm*mesh.support(:,1);
% spatially normalise J all gmsurface 
%J_all_gmsurface_norm=J_all_gmsurface./mesh_support_gmsurface;
%plot J all gmsurface
%figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface],gmSurfaceMesh.face);view(0,90);
%cb = colorbar('horiz');
%plot J all gmsurface spatially normalized
%figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface_norm],gmSurfaceMesh.face);view(0,90);
%cb = colorbar('horiz');
%title("Spatially Normalised Jacobian - GM Smooth Week30")
%xlabel("x / mm");
%ylabel("y / mm");
%zlabel("z / mm");

%HbO_lhs = zeros(size(J_gmsurface,2),size(mua1_lhs,2));
%Hb_lhs = zeros(size(J_gmsurface,2),size(mua1_lhs,2));
%HbT_lhs = zeros(size(J_gmsurface,2),size(mua1_lhs,2));

%doing inversion and just saving HbX for gmsurface smooth
for i=1:size(mua1_lhs,2)
    %full size mesh , single time point
    tmp = inv(E)*[mua1_lhs(:,i)'; mua2_lhs(:,i)']; %just using mu a
    HbO_lhs(:,1) = tmp(1,:)';
    Hb_lhs(:,1) = tmp(2,:)';
    HbT_lhs(:,1) = sum(tmp);

    %convert to gmsurface size
    %HbT_all_m_hypo_allNodes = zeros(size(mesh.nodes,1),1);
    %HbT_all_m_hypo_allNodes(cortex_nodes) = HbT_all_m_hypo(:,1,1);
    %HbT_all_m_hypo_gmsurface= vol2gm*HbT_all_m_hypo_allNodes;

    HbO_lhs_gmsurface(:,i) = vol2gm*HbO_lhs;
    Hb_lhs_gmsurface(:,i) = vol2gm*Hb_lhs;
    HbT_lhs_gmsurface(:,i) = vol2gm*HbT_lhs;
end

%full size mesh, all time points
% for i = 1 : size(mua1_lhs,2)
%     tmp = inv(E)*[mua1_lhs(:,i)'; mua2_lhs(:,i)']; %just using mu a
% 
%     %tmp = inv(E)*[mus1_lhs(:,i)'; mus2_lhs(:,i)']; %just using kappa
% 
%     %tmp = inv(E)*[mua1_lhs(:,i)' mus1_lhs(:,i)'; mua2_lhs(:,i)' mus2_lhs(:,i)']; %Using mua and mus'
%     HbO_lhs(:,i) = tmp(1,:)';
%     Hb_lhs(:,i) = tmp(2,:)';
%     HbT_lhs(:,i) = sum(tmp);
% end
clear tmp mua1_lhs mua2_lhs
%
%mesh_recon = mesh;
%ind = reshape(mesh_recon.region(mesh_recon.elements),[],1);
%%ind = reshape(ind>=3,[],4);
%ind = reshape(ind==3,[],4); %just GM nodes 06 05 24
%ind = sum(ind,2);
%ind = find(ind==4);
%[mesh3.elements,mesh3.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%%[mesh_all.elements,mesh_all.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%%[meshZ.elements,meshZ.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%%cortex_nodes = find(mesh.region>=3);
%cortex_nodes = find(mesh.region==3);

%max_dC = max(abs([min(min([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)])) max(max([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)]))]));
%min_dC = max(abs([min(min([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)])) max(max([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)]))]));


%dC_min_max = [min(min([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)])) max(max([Hb_lhs(cortex_nodes,:) HbO_lhs(cortex_nodes,:) HbT_lhs(cortex_nodes,:)]))];
%just min and max from cortex nodes that are top 95% of All Ch Normalized J
%09 05 24
%dC_min_max = [min(min([Hb_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),:) HbO_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),:) HbT_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),:)])) max(max([Hb_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),:) HbO_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),:) HbT_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),:)]))];

%min max from gmsurface HbX 25 06 24
dC_min_max = [min(min([Hb_lhs_gmsurface(:,:) HbO_lhs_gmsurface(:,:) HbT_lhs_gmsurface(:,:)])) max(max([Hb_lhs_gmsurface(:,:) HbO_lhs_gmsurface(:,:) HbT_lhs_gmsurface(:,:)]))];



% HbO_CL_lhs = HbO_lhs; %comment if using DS. 
% Hb_CL_lhs = Hb_lhs; %comment if using DS.
% HbT_CL_lhs = HbT_lhs; %comment if using DS.

%time_points = [find(PD_data.t==0) find(PD_data.t==300) find(PD_data.t==600) find(PD_data.t==900) find(PD_data.t==1200) find(PD_data.t==1500) find(PD_data.t==1800)];

% for i=1:size(time_points,2)
%     mesh3.data = HbO_lhs(:,time_points(i));
%     plotniceimages_1(mesh3,mesh_recon);
%     colorbar('horiz');
%     max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
%     min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
%     %caxis([min_HbO_cortex max_HbO_cortex]);
%     caxis([dC_min_max(1) dC_min_max(2)]);
%     title("HbO at time= "+(PD_data.t(time_points(i)))/60+" minutes")
% end

%% Find nodes that have 95% of max jacobian
%JUST the GOOD IDX
%% plotting

%Plot changes across 6 time points
%caxis_lim = [-3 3];

%making caxis lim from all time points
%caxis_lim = [dC_min_max(1) dC_min_max(2)];
%caxis_lim = [ -max(abs(caxis_lim)) max(abs(caxis_lim))]; %this makes the cbar symtertical about 0


time_points_6 = 1:120:601; %the 6 time points to plot on tomography
%get min and max from the 6 time points 09 05 24
%dC_min_max = [min(min([Hb_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),time_points_6) HbO_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),time_points_6) HbT_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),time_points_6)])) max(max([Hb_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),time_points_6) HbO_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),time_points_6) HbT_lhs(cortex_nodes(masks.mask_95_cortexnodes_allch),time_points_6)]))];

%using gmsurface mask top95% J gmsurface
dC_min_max = [min(min([Hb_lhs_gmsurface(masks.mask_j_thresh_cortexnodes_allch,time_points_6) HbO_lhs_gmsurface(masks.mask_j_thresh_cortexnodes_allch,time_points_6) HbT_lhs_gmsurface(masks.mask_j_thresh_cortexnodes_allch,time_points_6)])) max(max([Hb_lhs_gmsurface(masks.mask_j_thresh_cortexnodes_allch  ,time_points_6) HbO_lhs_gmsurface(masks.mask_j_thresh_cortexnodes_allch,time_points_6) HbT_lhs_gmsurface(masks.mask_j_thresh_cortexnodes_allch ,time_points_6)]))];

%masks.mask_j_thresh_cortexnodes_allch

caxis_lim = [dC_min_max(1) dC_min_max(2)];
caxis_lim = [ -max(abs(caxis_lim)) max(abs(caxis_lim))]; %this makes the cbar symtertical about 0


%set HbX nodes in bottom 5% of Jacobian to zero
%HbO_lhs(cortex_nodes(masks.mask_nodes_bottom_5_cortexnodes_allch  ),:) = 0;
%Hb_lhs(cortex_nodes(masks.mask_nodes_bottom_5_cortexnodes_allch  ),:) = 0;
%HbT_lhs(cortex_nodes(masks.mask_nodes_bottom_5_cortexnodes_allch  ),:) = 0;

%set HbX nodes in bottom 5% of Jacobian to zero
HbO_lhs_gmsurface(masks.mask_nodes_bottom_j_thresh_cortexnodes_allch,:) = 0;
Hb_lhs_gmsurface(masks.mask_nodes_bottom_j_thresh_cortexnodes_allch,:) = 0;
HbT_lhs_gmsurface(masks.mask_nodes_bottom_j_thresh_cortexnodes_allch,:) = 0;

%mask_j_thresh_cortexnodes_allch
%masks.mask_nodes_bottom_j_thresh_cortexnodes_allch

switch PD_data.time_window
    case "All"
        time_points_6 = 1:size(PD_data.t,1)/6:size(PD_data.t,1); %the 6 time points to plot on tomography
end

plot_n = [1 2 3 5 6 7]; %subplot numbers

%  HbO 
figure()
for i=1:size(time_points_6,2)
    subplot(2,4,plot_n(i))
    %%%mesh3.data = HbO_lhs(:,time_points(i));
    %mesh3.data = HbO_lhs(:,find(time_points==time_points_6(i)));
    %plotniceimages_1_greyJET(mesh3,mesh_recon);
    plotmesh_iso2([gmSurfaceMesh.node HbO_lhs_gmsurface(:,find(time_points==time_points_6(i)))],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    cb = colorbar('horiz'); 
    ylabel(cb,'\Delta HbO (\muM)','FontSize',10,'Rotation',0)
    %max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
    %min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
    caxis([caxis_lim(1) caxis_lim(2)]);
    title("HbO at time= "+round((PD_data.t(time_points_6(i) )/60))+" minutes")
end
subplot(2,4,4)
switch PD_data.eventType
    case "m_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(PD.eventN,1):PD.events_start_end.m_hypo(PD.eventN,2)),'b-o')
    case "S_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(PD.eventN,1):PD.events_start_end.S_hypo(PD.eventN,2)),'b-o')
end
title("TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("Time / mins")
ylabel("BGC / mg/DL")
xline(15,'g--');xline(PD_time(end)-15,'r--');
xline(PD_data.time_window_t(1),'m--');xline(PD_data.time_window_t(2),'m--');
yline(72,'k--');
subplot(2,4,8)
%JUST the GOOD IDX
%goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
%meshJ = mesh3;
%meshJ.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
%max_j = max(meshJ.data(1,cortex_nodes));
%min_j = min(meshJ.data(1,cortex_nodes));
%plotniceimages_1_greyJET(meshJ,mesh_recon);
plotmesh_iso2([gmSurfaceMesh.node goodch_J_gmsurface],gmSurfaceMesh.face);view(0,90);
%plotmesh_iso2([gmSurfaceMesh.node HbO_lhs_gmsurface(:,find(time_points==time_points_6(i)))],gmSurfaceMesh.face)
view(0,90); %TOP VIEW % used for infant week 30
title("J - goodch idx only")
caxis([min_j max_j]);
cb = colorbar('horiz'); 
ylabel(cb,'\DeltaI  \Delta\mu_a^{-1}','FontSize',10,'Rotation',0)
% Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_HbO.png")
end

figure()
%JUST the GOOD IDX
%goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
%meshJ = mesh3;
%meshJ.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
%max_j = max(meshJ.data(1,cortex_nodes));
%min_j = min(meshJ.data(1,cortex_nodes));
%plotniceimages_1_Jacobian(meshJ,mesh_recon);
plotmesh_JAC([gmSurfaceMesh.node goodch_J_gmsurface],gmSurfaceMesh.face);view(0,90);
title("J - goodch idx only")
caxis([min_j max_j]);
cb = colorbar('horiz'); 
ylabel(cb,'\DeltaI  \Delta\mu_a^{-1}','FontSize',10,'Rotation',0)
% Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_JGoodCh.png")
end


%Hb
figure()
for i=1:size(time_points_6,2)
    subplot(2,4,plot_n(i))
    %mesh3.data = HbO_lhs(:,time_points(i));
    %mesh3.data = Hb_lhs(:,find(time_points==time_points_6(i)));
    %plotniceimages_1_greyJET(mesh3,mesh_recon);
    plotmesh_iso2([gmSurfaceMesh.node HbO_lhs_gmsurface(:,find(time_points==time_points_6(i)))],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    cb = colorbar('horiz'); 
    ylabel(cb,'\Delta Hb (\muM)','FontSize',10,'Rotation',0)
    %max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
    %min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
    caxis([caxis_lim(1) caxis_lim(2)]);
    title("Hb at time= "+round((PD_data.t(time_points_6(i) )/60))+" minutes")
end
subplot(2,4,4)
switch PD_data.eventType
    case "m_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(PD.eventN,1):PD.events_start_end.m_hypo(PD.eventN,2)),'b-o')
    case "S_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(PD.eventN,1):PD.events_start_end.S_hypo(PD.eventN,2)),'b-o')
end
title("TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("Time / mins")
ylabel("BGC / mg/DL")
xline(15,'g--');xline(PD_time(end)-15,'r--');
xline(PD_data.time_window_t(1),'m--');xline(PD_data.time_window_t(2),'m--');
yline(72,'k--');
subplot(2,4,8)
%JUST the GOOD IDX
%goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
%meshJ = mesh3;
%meshJ.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
%max_j = max(meshJ.data(1,cortex_nodes));
%min_j = min(meshJ.data(1,cortex_nodes));
%plotniceimages_1_greyJET(meshJ,mesh_recon);
plotmesh_iso2([gmSurfaceMesh.node goodch_J_gmsurface],gmSurfaceMesh.face);view(0,90);
title("J - goodch idx only")
caxis([min_j max_j]);
cb = colorbar('horiz'); 
ylabel(cb,'\DeltaI  \Delta\mu_a^{-1}','FontSize',10,'Rotation',0)
% Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_Hb.png")
end
%HbT
figure()
for i=1:size(time_points_6,2)
    subplot(2,4,plot_n(i))
    %mesh3.data = HbO_lhs(:,time_points(i));
    %mesh3.data = HbT_lhs(:,find(time_points==time_points_6(i)));
    %plotniceimages_1_greyJET(mesh3,mesh_recon);
    plotmesh_iso2([gmSurfaceMesh.node HbO_lhs_gmsurface(:,find(time_points==time_points_6(i)))],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    cb = colorbar('horiz'); 
    ylabel(cb,'\Delta HbT (\muM)','FontSize',10,'Rotation',0)
    %max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
    %min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
    caxis([caxis_lim(1) caxis_lim(2)]);
    title("HbT at time= "+round((PD_data.t(time_points_6(i) )/60))+" minutes")
end
subplot(2,4,4)
switch PD_data.eventType
    case "m_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(PD.eventN,1):PD.events_start_end.m_hypo(PD.eventN,2)),'b-o')
    case "S_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(PD.eventN,1):PD.events_start_end.S_hypo(PD.eventN,2)),'b-o')
end
title("TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("Time / mins")
ylabel("BGC / mg/DL")
xline(15,'g--');xline(PD_time(end)-15,'r--');
xline(PD_data.time_window_t(1),'m--');xline(PD_data.time_window_t(2),'m--');
yline(72,'k--');
subplot(2,4,8)
%JUST the GOOD IDX
%goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
%meshJ = mesh3;
%meshJ.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
%max_j = max(meshJ.data(1,cortex_nodes));
%min_j = min(meshJ.data(1,cortex_nodes));
plotmesh_iso2([gmSurfaceMesh.node goodch_J_gmsurface],gmSurfaceMesh.face);view(0,90);
%plotniceimages_1_greyJET(meshJ,mesh_recon);
title("J - goodch idx only")
caxis([min_j max_j]);
cb = colorbar('horiz'); 
ylabel(cb,'\DeltaI  \Delta\mu_a^{-1}','FontSize',10,'Rotation',0)
% Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_HbT.png")
end
%% Get 6 regions of cortex
landmark_EEG1020 = [
	-23.57	14.32	24.47;
	-1.09	15.9	37.04;
 	21.63	12.93	26.08;
 	-32.89	-16.63	35.48;
 	-2.01	-17.37	51.07;
 	29.67	-16.27	36.75;
 	-24.76	-50.57	27.24;
 	-1.33	-52.76	39.53;
 	22.94	-49.58	29.2 ];

% plot3(landmark_EEG1020(:,1),...
%         landmark_EEG1020(:,2),...
%         landmark_EEG1020 (:,3),'cx',...
%         'LineWidth',2,'MarkerSize',8);
%mesh3.data = abs(randn(size(HbO_lhs,1),1))+100; %RANDOMISE AND SET SO IT'S AT MAX CBAR
%GETTTING CORTEX ZONES 
%F3  FZ  F4
%C3  CZ  C4
%P3  PZ  P4
%radius_LM = 20;
%replace mesh.nodes with gmSurfaceMesh.node
%landmark_EEG1020_cortexF3 = find(pdist2(landmark_EEG1020(1,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexF3) = -5;
%landmark_EEG1020_cortexFZ = find(pdist2(landmark_EEG1020(2,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexFZ) = -4;
%landmark_EEG1020_cortexF4 = find(pdist2(landmark_EEG1020(3,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexF4) = -3;
%landmark_EEG1020_cortexC3 = find(pdist2(landmark_EEG1020(4,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexC3) = -2;
%landmark_EEG1020_cortexCZ = find(pdist2(landmark_EEG1020(5,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexCZ) = -1;
%landmark_EEG1020_cortexC4 = find(pdist2(landmark_EEG1020(6,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexC4) = -0;
%landmark_EEG1020_cortexP3 = find(pdist2(landmark_EEG1020(7,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexP3) = 1;
%landmark_EEG1020_cortexPZ = find(pdist2(landmark_EEG1020(8,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexPZ) = 2;
%landmark_EEG1020_cortexP4 = find(pdist2(landmark_EEG1020(9,:),gmSurfaceMesh.nodes)<radius_LM);
%mesh3.data(landmark_EEG1020_cortexP4) = 3;

%figure()
%plotniceimages_1_greyJET(mesh3,mesh_recon);
%cb = colorbar('horiz'); 
%ylabel(cb,'EEG LM Zone','FontSize',10,'Rotation',0)
%max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
%min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
%caxis([caxis_lim(1) caxis_lim(2)]);
%title("9 EEG zones on cortex F3z4 C3z4 P3z4")

%Intersect between the nodes in the EEG LM SPHERES and Cortex nodes
%landmark_EEG1020_cortexF3 = intersect(landmark_EEG1020_cortexF3,cortex_nodes);
%landmark_EEG1020_cortexFZ = intersect(landmark_EEG1020_cortexFZ,cortex_nodes);
%landmark_EEG1020_cortexF4 = intersect(landmark_EEG1020_cortexF4,cortex_nodes);
%landmark_EEG1020_cortexC3 = intersect(landmark_EEG1020_cortexC3,cortex_nodes);
%landmark_EEG1020_cortexCZ = intersect(landmark_EEG1020_cortexCZ,cortex_nodes);
%landmark_EEG1020_cortexC4 = intersect(landmark_EEG1020_cortexC4,cortex_nodes);
%landmark_EEG1020_cortexP3 = intersect(landmark_EEG1020_cortexP3,cortex_nodes);
%landmark_EEG1020_cortexPZ = intersect(landmark_EEG1020_cortexPZ,cortex_nodes);
%landmark_EEG1020_cortexP4 = intersect(landmark_EEG1020_cortexP4,cortex_nodes);

ROI_data = ones(size(gmSurfaceMesh.node,1),1)*4 ;%  abs(randn(size(gmSurfaceMesh.node,1),1))+100; %RANDOMISE AND SET SO IT'S AT MAX CBAR
%GETTTING CORTEX ZONES 
%F3  FZ  F4
%C3  CZ  C4
%P3  PZ  P4
radius_LM = 20;
landmark_EEG1020_cortexF3 = find(pdist2(landmark_EEG1020(1,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexF3) = -5;
landmark_EEG1020_cortexFZ = find(pdist2(landmark_EEG1020(2,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexFZ) = -4;
landmark_EEG1020_cortexF4 = find(pdist2(landmark_EEG1020(3,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexF4) = -3;
landmark_EEG1020_cortexC3 = find(pdist2(landmark_EEG1020(4,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexC3) = -2;
landmark_EEG1020_cortexCZ = find(pdist2(landmark_EEG1020(5,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexCZ) = -1;
landmark_EEG1020_cortexC4 = find(pdist2(landmark_EEG1020(6,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexC4) = -0;
landmark_EEG1020_cortexP3 = find(pdist2(landmark_EEG1020(7,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexP3) = 1;
landmark_EEG1020_cortexPZ = find(pdist2(landmark_EEG1020(8,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexPZ) = 2;
landmark_EEG1020_cortexP4 = find(pdist2(landmark_EEG1020(9,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexP4) = 3;

EEG_lm_9ROI = [];

figure()
plotmesh_iso2([gmSurfaceMesh.node ROI_data],gmSurfaceMesh.face)
clim([-5 4])
cb = colorbar('horiz');
ylabel(cb,'EEG LM Zone','FontSize',10,'Rotation',0)
title("9 EEG zones on cortex F3z4 C3z4 P3z4")

PD.data.landmarkEEGnodes.landmark_EEG1020_cortexF3 = landmark_EEG1020_cortexF3;
PD.data.landmarkEEGnodes.landmark_EEG1020_cortexFZ = landmark_EEG1020_cortexFZ;
PD.data.landmarkEEGnodes.landmark_EEG1020_cortexF4 = landmark_EEG1020_cortexF4;
PD.data.landmarkEEGnodes.landmark_EEG1020_cortexC3 = landmark_EEG1020_cortexC3;
PD.data.landmarkEEGnodes.landmark_EEG1020_cortexCZ = landmark_EEG1020_cortexCZ;
PD.data.landmarkEEGnodes.landmark_EEG1020_cortexC4 = landmark_EEG1020_cortexC4;
PD.data.landmarkEEGnodes.landmark_EEG1020_cortexP3 = landmark_EEG1020_cortexP3;
PD.data.landmarkEEGnodes.landmark_EEG1020_cortexPZ = landmark_EEG1020_cortexPZ;
PD.data.landmarkEEGnodes.landmark_EEG1020_cortexP4 = landmark_EEG1020_cortexP4;

%GET MIN MAX AXIS
EGG1020_HbX_max = max([max(max([max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexF3,:))) max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexFZ,:))) max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexF4,:)));
max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexC3,:))) max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexCZ,:))) max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexC4,:)));
max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexP3,:))) max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexPZ,:))) max(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexP4,:)));])) ;

max(max([max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexF3,:))) max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexFZ,:))) max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexF4,:)));
max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexC3,:))) max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexCZ,:))) max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexC4,:)));
max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexP3,:))) max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexPZ,:))) max(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexP4,:)));])) ;

max(max([max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexF3,:))) max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexFZ,:))) max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexF4,:)));
max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexC3,:))) max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexCZ,:))) max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexC4,:)));
max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexP3,:))) max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexPZ,:))) max(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexP4,:)));])) ]);

EGG1020_HbX_min =min([min(min([min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexF3,:))) min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexFZ,:))) min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexF4,:)));
min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexC3,:))) min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexCZ,:))) min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexC4,:)));
min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexP3,:))) min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexPZ,:))) min(mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexP4,:)));])) ;

min(min([min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexF3,:))) min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexFZ,:))) min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexF4,:)));
min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexC3,:))) min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexCZ,:))) min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexC4,:)));
min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexP3,:))) min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexPZ,:))) min(mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexP4,:)));])) ;

min(min([min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexF3,:))) min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexFZ,:))) min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexF4,:)));
min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexC3,:))) min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexCZ,:))) min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexC4,:)));
min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexP3,:))) min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexPZ,:))) min(mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexP4,:)));])) ]);

figure()
subplot(3,3,1)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexF3,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexF3,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexF3,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("F3")
subplot(3,3,2)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexFZ,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexFZ,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexFZ,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("FZ")
subplot(3,3,3)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexF4,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexF4,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexF4,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("F4")
subplot(3,3,4)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexC3,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexC3,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexC3,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("C3")
subplot(3,3,5)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexCZ,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexCZ,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexCZ,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("CZ")
subplot(3,3,6)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexC4,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexC4,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexC4,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("C4")
subplot(3,3,7)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexP3,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexP3,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexP3,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("P3")
subplot(3,3,8)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexPZ,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexPZ,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexPZ,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("PZ")
subplot(3,3,9)
plot(PD_data.t,mean(HbO_lhs_gmsurface(landmark_EEG1020_cortexP4,:)),'r-')
hold on
plot(PD_data.t,mean(Hb_lhs_gmsurface(landmark_EEG1020_cortexP4,:)),'b-')
plot(PD_data.t,mean(HbT_lhs_gmsurface(landmark_EEG1020_cortexP4,:)),'g-')
ylim([EGG1020_HbX_min EGG1020_HbX_max])
xlabel("Time / s")
ylabel("\Delta HbX / \muM")
title("P4")
sgtitle("TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
% Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_SpectEEG1020HbX.png")
end
%% Normalized HbT, HbO, Hb (div by max of Hbo Hb HbT) -1 to 1
normal_factor = max([abs(dC_min_max(1)) abs(dC_min_max(2))]);

HbT_lhs_norm = HbT_lhs_gmsurface./normal_factor;
HbO_lhs_norm = HbO_lhs_gmsurface./normal_factor;
Hb_lhs_norm = Hb_lhs_gmsurface./normal_factor;

%HbO
figure()
for i=1:size(time_points_6,2)
    subplot(2,4,plot_n(i))
    %mesh3.data = HbO_lhs(:,time_points(i));
    %mesh3.data = HbO_lhs_norm(:,find(time_points==time_points_6(i)));
    %plotniceimages_1_greyJET(mesh3,mesh_recon);
    plotmesh_iso2([gmSurfaceMesh.node HbO_lhs_norm(:,find(time_points==time_points_6(i)))],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    cb = colorbar('horiz'); 
    ylabel(cb,'\Delta HbO (Norm. A.U)','FontSize',10,'Rotation',0)
    %max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
    %min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
    caxis([-1 1]);
    title("Norm. HbO at time= "+round((PD_data.t(time_points_6(i) )/60))+" minutes")
end
subplot(2,4,4)
switch PD_data.eventType
    case "m_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(PD.eventN,1):PD.events_start_end.m_hypo(PD.eventN,2)),'b-o')
    case "S_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(PD.eventN,1):PD.events_start_end.S_hypo(PD.eventN,2)),'b-o')
end
title("TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("Time / mins")
ylabel("BGC / mg/DL")
xline(15,'g--');xline(PD_time(end)-15,'r--');
xline(PD_data.time_window_t(1),'m--');xline(PD_data.time_window_t(2),'m--');
yline(72,'k--');
subplot(2,4,8)
%JUST the GOOD IDX
%goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
%meshJ = mesh3;
%meshJ.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
%max_j = max(meshJ.data(1,cortex_nodes));
%min_j = min(meshJ.data(1,cortex_nodes));
%plotniceimages_1_greyJET(meshJ,mesh_recon);
plotmesh_iso2([gmSurfaceMesh.node goodch_J_gmsurface],gmSurfaceMesh.face);view(0,90);
title("J - goodch idx only")
caxis([min_j max_j]);
cb = colorbar('horiz'); 
ylabel(cb,'\DeltaI  \Delta\mu_a^{-1}','FontSize',10,'Rotation',0)
% Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_HbO_NORM.png")
end

%Hb
figure()
for i=1:size(time_points_6,2)
    subplot(2,4,plot_n(i))
    %mesh3.data = HbO_lhs(:,time_points(i));
    %mesh3.data = Hb_lhs_norm(:,find(time_points==time_points_6(i)));
    %plotniceimages_1_greyJET(mesh3,mesh_recon);
    plotmesh_iso2([gmSurfaceMesh.node Hb_lhs_norm(:,find(time_points==time_points_6(i)))],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    cb = colorbar('horiz'); 
    ylabel(cb,'\Delta Hb (Norm. A.U)','FontSize',10,'Rotation',0)
    %max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
    %min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
    caxis([-1 1]);
    title("Norm. Hb at time= "+round((PD_data.t(time_points_6(i) )/60))+" minutes")
end
subplot(2,4,4)
switch PD_data.eventType
    case "m_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(PD.eventN,1):PD.events_start_end.m_hypo(PD.eventN,2)),'b-o')
    case "S_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(PD.eventN,1):PD.events_start_end.S_hypo(PD.eventN,2)),'b-o')
end
title("TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("Time / mins")
ylabel("BGC / mg/DL")
xline(15,'g--');xline(PD_time(end)-15,'r--');
xline(PD_data.time_window_t(1),'m--');xline(PD_data.time_window_t(2),'m--');
yline(72,'k--');
subplot(2,4,8)
%JUST the GOOD IDX
%goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
%meshJ = mesh3;
%meshJ.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
%max_j = max(meshJ.data(1,cortex_nodes));
%min_j = min(meshJ.data(1,cortex_nodes));
%plotniceimages_1_greyJET(meshJ,mesh_recon);
plotmesh_iso2([gmSurfaceMesh.node goodch_J_gmsurface],gmSurfaceMesh.face);view(0,90);
title("J - goodch idx only")
caxis([min_j max_j]);
cb = colorbar('horiz'); 
ylabel(cb,'\DeltaI  \Delta\mu_a^{-1}','FontSize',10,'Rotation',0)
% Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_Hb_NORM.png")
end


%HbT
figure()
for i=1:size(time_points_6,2)
    subplot(2,4,plot_n(i))
    %mesh3.data = HbO_lhs(:,time_points(i));
    %mesh3.data = HbT_lhs_norm(:,find(time_points==time_points_6(i)));
    %plotniceimages_1_greyJET(mesh3,mesh_recon);
    plotmesh_iso2([gmSurfaceMesh.node HbT_lhs_norm(:,find(time_points==time_points_6(i)))],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    cb = colorbar('horiz'); 
    ylabel(cb,'\Delta HbT (Norm. A.U)','FontSize',10,'Rotation',0)
    %max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
    %min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
    caxis([-1 1]);
    title("Norm. HbT at time= "+round((PD_data.t(time_points_6(i) )/60))+" minutes")
end
subplot(2,4,4)
switch PD_data.eventType
    case "m_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(PD.eventN,1):PD.events_start_end.m_hypo(PD.eventN,2)),'b-o')
    case "S_hypo"
        plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(PD.eventN,1):PD.events_start_end.S_hypo(PD.eventN,2)),'b-o')
end
title("TimeW "+PD_data.time_window+" PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+". nirs")
xlabel("Time / mins")
ylabel("BGC / mg/DL")
xline(15,'g--');xline(PD_time(end)-15,'r--');
xline(PD_data.time_window_t(1),'m--');xline(PD_data.time_window_t(2),'m--');
yline(72,'k--');
subplot(2,4,8)
%JUST the GOOD IDX
%goodch_J = sum(J_amp(PD_data.goodch_idx(1:end/2),:)); %this sums the jacobian for all measurements and nodes
% %so we are left with the total sensitivity for EACH node
%meshJ = mesh3;
%meshJ.data = goodch_J./mesh.support(:,1)'; %set mesh3.data as the total sensitivty
%max_j = max(meshJ.data(1,cortex_nodes));
%min_j = min(meshJ.data(1,cortex_nodes));
%plotniceimages_1_greyJET(meshJ,mesh_recon);
plotmesh_iso2([gmSurfaceMesh.node goodch_J_gmsurface],gmSurfaceMesh.face);view(0,90);
title("J - goodch idx only")
caxis([min_j max_j]);
cb = colorbar('horiz'); 
ylabel(cb,'\DeltaI  \Delta\mu_a^{-1}','FontSize',10,'Rotation',0)
% Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_HbT_NORM.png")
end


%% PLOT Z score
% Z score: Run this if analysing delta C from start to end of event
%start_baseline_ds = floor(start_baseline/ds);
%end_baseline_ds = floor(end_baseline/ds);
%start_peak_event_ds = floor(start_peak_event/ds);
%end_peak_event_ds = floor(end_peak_event/ds);

%mean of activation - mean of baseline / STD baseline
% HbO_lhs_Z = (mean(  HbO_lhs(:,12:end)') - mean(HbO_lhs(:,1:11)')   ./ std( HbO_lhs(:,1:11)') ) ;
% HbO_lhs_Z = HbO_lhs_Z';

%old 09 05 24
%HbO_lhs_Z = (mean(  HbO_lhs(:,121:end)') - mean(HbO_lhs(:,1:120)')   ./ std( HbO_lhs(:,1:120)') ) ;
%HbO_lhs_Z = HbO_lhs_Z';

%([5;5]-[4;3])./[0.05 ;0.03]
%Z = (Mean stim - mean baseline) / SD baseline
%Wan Chu Su 2023 fnirs

HbO_lhs_Z = ( mean(  HbO_lhs_gmsurface(:,121:end)') - mean(HbO_lhs_gmsurface(:,1:120)') )   ./ std( HbO_lhs_gmsurface(:,1:120)')  ;
HbO_lhs_Z = HbO_lhs_Z';

Hb_lhs_Z = (mean( Hb_lhs_gmsurface(:,121:end)') - mean(Hb_lhs_gmsurface(:,1:120)') )  ./ std(Hb_lhs_gmsurface(:,1:120)') ;
Hb_lhs_Z = Hb_lhs_Z';

HbT_lhs_Z = (mean(HbT_lhs_gmsurface(:,121:end)') - mean(HbT_lhs_gmsurface(:,1:120)') )  ./ std(HbT_lhs_gmsurface(:,1:120)') ;
HbT_lhs_Z = HbT_lhs_Z';

%set nodes of the 5% mask criteria to 0 (if mask = 0, node is in bottom 5%
%of good ch spatial norm J)
%HbO_lhs_Z(find(mesh.mask_spatial_norm == 0),1) = 0;
%Hb_lhs_Z(find(mesh.mask_spatial_norm == 0),1) = 0;
%HbT_lhs_Z(find(mesh.mask_spatial_norm == 0),1) = 0;

%HbO_lhz_Z(cortex_nodes(masks.mask_nodes_bottom_5_cortexnodes_allch),:) = 0;
%Hb_lhz_Z(cortex_nodes(masks.mask_nodes_bottom_5_cortexnodes_allch),:) = 0;
%HbT_lhz_Z(cortex_nodes(masks.mask_nodes_bottom_5_cortexnodes_allch),:) = 0;

HbO_lhz_Z(masks.mask_nodes_bottom_j_thresh_cortexnodes_allch,:) = 0;
Hb_lhz_Z(masks.mask_nodes_bottom_j_thresh_cortexnodes_allch,:) = 0;
HbT_lhz_Z(masks.mask_nodes_bottom_j_thresh_cortexnodes_allch,:) = 0;

%mask_j_thresh_cortexnodes_allch
%masks.mask_nodes_bottom_j_thresh_cortexnodes_allch

%HbT_lhs_Z = (mean(HbT_lhs(:,27:end)') - mean(HbT_lhs(:,1:27)') )  ./ std(HbT_lhs(:,1:27)') ;
%HbT_lhs_Z = HbT_lhs_Z';

%Z_min_max = [min(min([Hb_lhs_Z(cortex_nodes,:) HbO_lhs_Z(cortex_nodes,:) HbT_lhs_Z(cortex_nodes,:)])) max(max([Hb_lhs_Z(cortex_nodes,:) HbO_lhs_Z(cortex_nodes,:) HbT_lhs_Z(cortex_nodes,:)]))];

%Z_min_max_HbO = [min(min([HbO_lhs_Z(cortex_nodes,:)])) max(max([ HbO_lhs_Z(cortex_nodes,:) ]))];
%Z_min_max_Hb = [min(min([Hb_lhs_Z(cortex_nodes,:)])) max(max([ Hb_lhs_Z(cortex_nodes,:) ]))];
%Z_min_max_HbT = [min(min([HbT_lhs_Z(cortex_nodes,:)])) max(max([ HbT_lhs_Z(cortex_nodes,:) ]))];

Z_min_max_HbO = [min(min([HbO_lhs_Z(:,:)])) max(max([ HbO_lhs_Z(:,:) ]))];
Z_min_max_Hb = [min(min([Hb_lhs_Z(:,:)])) max(max([ Hb_lhs_Z(:,:) ]))];
Z_min_max_HbT = [min(min([HbT_lhs_Z(:,:)])) max(max([ HbT_lhs_Z(:,:) ]))];

Z_min_max_HbO = [ -max(abs(Z_min_max_HbO)) max(abs(Z_min_max_HbO))]; %this makes the cbar symtertical about 0
Z_min_max_Hb = [ -max(abs(Z_min_max_Hb)) max(abs(Z_min_max_Hb))]; %this makes the cbar symtertical about 0
Z_min_max_HbT = [ -max(abs(Z_min_max_HbT)) max(abs(Z_min_max_HbT))]; %this makes the cbar symtertical about 0

%meshZ=mesh3;
%meshZ.data = HbO_lhs_Z;


    figure()
    hold on
    subplot(1,3,1)
    %plotniceimages_1_greyJET(meshZ,mesh_recon);
    plotmesh_iso2([gmSurfaceMesh.node HbO_lhs_Z],gmSurfaceMesh.face);view(0,90);
    cb = colorbar('horiz'); 
    ylabel(cb,'Z score / A.U','FontSize',10,'Rotation',0)
    %caxis([-5 5]);
    caxis([Z_min_max_HbO(1) Z_min_max_HbO(2)]);
    title("Z score HbO, Subject "+num2str(PD.subjectN)+" "+PD.eventType+" event N "+num2str(PD.eventN)+"")
%meshZ.data = Hb_lhs_Z;
    
    subplot(1,3,2)
    %plotniceimages_1_greyJET(meshZ,mesh_recon);
    plotmesh_iso2([gmSurfaceMesh.node Hb_lhs_Z],gmSurfaceMesh.face);view(0,90);
    cb = colorbar('horiz'); 
    ylabel(cb,'Z score / A.U','FontSize',10,'Rotation',0)
    caxis([Z_min_max_Hb(1) Z_min_max_Hb(2)]);
    title("Z score Hb, Subject "+num2str(PD.subjectN)+" "+PD.eventType+" event N "+num2str(PD.eventN)+"")
%meshZ.data = HbT_lhs_Z;
    
    subplot(1,3,3)
    plotmesh_iso2([gmSurfaceMesh.node HbT_lhs_Z],gmSurfaceMesh.face);view(0,90);
    cb = colorbar('horiz'); 
    ylabel(cb,'Z score / A.U','FontSize',10,'Rotation',0)
    caxis([Z_min_max_HbT(1) Z_min_max_HbT(2)]);
    title("Z score HbT, Subject "+num2str(PD.subjectN)+" "+PD.eventType+" event N "+num2str(PD.eventN)+"")
    % Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    %set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_Z.png")
    end
    hold off
    


    %return values
    PD_data.HbO = HbO_lhs_gmsurface;
    PD_data.Hb = Hb_lhs_gmsurface;
    PD_data.HbT = HbT_lhs_gmsurface;
    PD_data.goodch_J_gmsurface = goodch_J_gmsurface;


end