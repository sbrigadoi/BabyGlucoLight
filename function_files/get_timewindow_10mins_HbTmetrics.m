function [HbT_metrics_PD_events] = get_timewindow_10mins_HbTmetrics(good_PD_events_subN,good_PD_events,time_window,eventType,HbT_all_gmsurface,gmSurfaceMesh,masks,eventType_subN , eventType_subN_idx )


%% metric 0 %n nodes x n events - T val with respect to 0 change

%HbT_all_T_val_gmsurface = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));

t_val = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
h_val_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
ci_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
pval_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));

for i=1:size(HbT_all_gmsurface,1);
    for j=1:size(HbT_all_gmsurface,3)
        [t_test(i,j),pval_t(i,j)] = ttest(HbT_all_gmsurface(i,:,j));
        [h_val_t(i,j),pval_t(i,j),ci_t,stats] = ttest(HbT_all_gmsurface(i,:,j));
        t_val(i,j) = stats.tstat;
    end
end 

mask_j_thresh_allnodes_allch_gmsurface = zeros(size(HbT_all_gmsurface,1),1);
mask_j_thresh_allnodes_allch_gmsurface(masks.mask_j_thresh_cortexnodes_allch,1)=1;

adj_p_t = ones(size(HbT_all_gmsurface,3),size(mask_j_thresh_allnodes_allch_gmsurface,1));

for i=1:size(HbT_all_gmsurface,3)
    [h_t, crit_p_t, adj_ci_cvrg_t, adj_p_mask_t(:,i)]=fdr_bh(pval_t(logical(mask_j_thresh_allnodes_allch_gmsurface)));
end

for i=1:size(adj_p_mask_t,2)
    adj_p_t(i,logical(mask_j_thresh_allnodes_allch_gmsurface)) = adj_p_mask_t(:,i);
end

t_val(~logical(mask_j_thresh_allnodes_allch_gmsurface),:)=0;
pval_t(~logical(mask_j_thresh_allnodes_allch_gmsurface),:)=-1;
adj_p_t(:,~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;


%% metric 1 %n nodes x n events
HbT_all_max_gmsurface = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));

    for i=1:size(HbT_all_gmsurface,1)
        [HbT_all_max HbT_all_max_idx] = (max(abs(HbT_all_gmsurface(i,121:end,:)))); %get max HbT from T=2:10mins
        HbT_all_max_idx = reshape(HbT_all_max_idx(1,1,:), 1, []);
    
        HbT_all_max = HbT_all_gmsurface(i,121+HbT_all_max_idx-1,:); %get max HbT from T=2:10mins
        HbT_all_max = diag(reshape(HbT_all_max,[],size(HbT_all_max,3))) ;%rehsape into sqaure, and get diag

        HbT_all_max_gmsurface(i,:) = HbT_all_max;

    end

%% metric 2 %n nodes x n events | Variance and range of HbT at each node across time - for each event
    
VAR_HbT_all_gmsurface = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
range_HbT_all_gmsurface = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));

    for i=1:size(HbT_all_gmsurface,1)
        for j=1:size(HbT_all_gmsurface,3)
            VAR_HbT_all_gmsurface(i,j) = var(HbT_all_gmsurface(i,:,j));
            range_HbT_all_gmsurface(i,j) = range(HbT_all_gmsurface(i,:,j));

        end

    end

%% metric 3 - mean total HbT each hemisphere

hemisphere_HbT_all_gmsurface = zeros(2,size(HbT_all_gmsurface,3));

%left and right hemisphere
LHS = find(gmSurfaceMesh.node(:,1)<0);
RHS = find(gmSurfaceMesh.node(:,1)>0);

for j=1:size(HbT_all_gmsurface,3)
    hemisphere_HbT_all_gmsurface(1,j) = mean(sum(HbT_all_gmsurface(LHS,121:end,j))) / (size(HbT_all_gmsurface,2) -121);
    hemisphere_HbT_all_gmsurface(2,j) = mean(sum(HbT_all_gmsurface(RHS,121:end,j))) / (size(HbT_all_gmsurface,2) -121);

end

%% metric 3.5 mean total HbT each node

mean_change_HbT_all_gmsurface = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
for j=1:size(HbT_all_gmsurface,3)
   mean_change_HbT_all_gmsurface(:,j) = mean(HbT_all_gmsurface(:,121:end,j),2); %/ (size(HbT_all_gmsurface,2) -121);
end

%% metric 4 - total HbT in 9 ROI
% Get 6 regions of cortex
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

EEG_9_ROI_HbT_all_gmsurface = zeros(9,size(HbT_all_gmsurface,3));

for j=1:size(HbT_all_gmsurface,3)
    EEG_9_ROI_HbT_all_gmsurface(1,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexF3,121:end,j))) / (size(HbT_all_gmsurface,2) -121);
    EEG_9_ROI_HbT_all_gmsurface(2,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexFZ,121:end,j))) / (size(HbT_all_gmsurface,2) -121);
    EEG_9_ROI_HbT_all_gmsurface(3,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexF4,121:end,j))) / (size(HbT_all_gmsurface,2) -121);

    EEG_9_ROI_HbT_all_gmsurface(4,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexC3,121:end,j))) / (size(HbT_all_gmsurface,2) -121);
    EEG_9_ROI_HbT_all_gmsurface(5,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexCZ,121:end,j))) / (size(HbT_all_gmsurface,2) -121);
    EEG_9_ROI_HbT_all_gmsurface(6,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexC4,121:end,j))) / (size(HbT_all_gmsurface,2) -121);

    EEG_9_ROI_HbT_all_gmsurface(7,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexP3,121:end,j))) / (size(HbT_all_gmsurface,2) -121);
    EEG_9_ROI_HbT_all_gmsurface(8,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexPZ,121:end,j))) / (size(HbT_all_gmsurface,2) -121);
    EEG_9_ROI_HbT_all_gmsurface(9,j) = mean(sum(HbT_all_gmsurface(landmark_EEG1020_cortexP4,121:end,j))) / (size(HbT_all_gmsurface,2) -121);

end

%% put metrics together

HbT_metrics_PD_events.HbT_all_max_gmsurface = HbT_all_max_gmsurface;
HbT_metrics_PD_events.VAR_HbT_all_gmsurface = VAR_HbT_all_gmsurface;
HbT_metrics_PD_events.range_HbT_all_gmsurface = range_HbT_all_gmsurface;
HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface = hemisphere_HbT_all_gmsurface;
HbT_metrics_PD_events.EEG_9_ROI_HbT_all_gmsurface = EEG_9_ROI_HbT_all_gmsurface;
HbT_metrics_PD_events.mean_change_HbT_all_gmsurface = mean_change_HbT_all_gmsurface;

HbT_metrics_PD_events.t_val = t_val;
HbT_metrics_PD_events.pval_t = pval_t;
HbT_metrics_PD_events.adj_p_t = adj_p_t;



end