function [event_jac_threshold_idx] = event_jac_threshold_check(headVolumeMesh,gmSurfaceMesh,J_gmsurface,good_ch_idx_events,masks,coverageThresh)



event_jac_threshold_idx = zeros(size(good_ch_idx_events,2),1);


%J_gmsurface
sum_J_gmsurface = zeros(size(good_ch_idx_events,2), size(J_gmsurface,1));
mask_J_threshold = zeros(size(good_ch_idx_events,2), size(J_gmsurface,1));
n_nodes_below_J_thresh = zeros(size(good_ch_idx_events,2),1);
mean_p_dist_gmnodes_J_thresh =  zeros(size(good_ch_idx_events,2),1);

for i=1:size(good_ch_idx_events,2)
    good_ch_idx_event_i = nonzeros(good_ch_idx_events(:,i));

    sum_J_gmsurface(i,:) = sum(J_gmsurface(:,good_ch_idx_event_i(1:end/2,1))');

    for j=1:size(J_gmsurface,1)
        if abs(sum_J_gmsurface(i,j)) > abs(coverageThresh)
            mask_J_threshold(i,j) = 1;
        end
    end

    %N nodes within top 95% J mask and below the J thresh mask
    %n_nodes_below_J_thresh(i) = size(masks.mask_95_cortexnodes_gmsurface_allch,1) - nnz(mask_J_threshold(i,masks.mask_95_cortexnodes_gmsurface_allch));
    n_nodes_below_J_thresh(i) = size(masks.mask_j_thresh_cortexnodes_allch,1) - nnz(mask_J_threshold(i,masks.mask_j_thresh_cortexnodes_allch));

%masks.mask_j_thresh_cortexnodes_allch

    %pdist between nodes within top 95% J mask and below the J thresh mask
    mean_p_dist_gmnodes_J_thresh(i) = mean(pdist(gmSurfaceMesh.node(find(mask_J_threshold(i,masks.mask_j_thresh_cortexnodes_allch)==0),:)));

    figure()
    plotmesh_iso2([gmSurfaceMesh.node mask_J_threshold(i,:)'] ,gmSurfaceMesh.face)
    xlabel('x / mm')
    ylabel('y / mm')
    zlabel('z / mm')
    view(0,90);
    cb = colorbar('horiz');
    title(num2str(i))

    %show 9ROI's

end

%%pdist between nodes within top 95% J mask and below the J thresh mask
%%mean(pdist(gmSurfaceMesh.node(find(mask_J_threshold(i,masks.mask_95_cortexnodes_gmsurface_allch)==0),:)));

%%N nodes within top 95% J mask and below the J thresh mask
%%n_nodes_below_J_thresh(i) = size(masks.mask_95_cortexnodes_gmsurface_allch,1) - nnz(mask_J_threshold(i,masks.mask_95_cortexnodes_gmsurface_allch))

% figure()
% plotmesh_iso2([gmSurfaceMesh.node mask_J_threshold(i,:)'] ,gmSurfaceMesh.face)
% xlabel('x / mm')
% ylabel('y / mm')
% zlabel('z / mm')
% view(0,90);
% cb = colorbar('horiz');
                                                        %row number was 2,
                                                        %now is 1,
                                                        %mask_J_threshold(2
    dist_pdist = pdist(gmSurfaceMesh.node(find(mask_J_threshold(1,masks.mask_j_thresh_cortexnodes_allch)==0),:));
    figure()
    hist(dist_pdist)

%% 9ROIS
%% Get 6 regions of cortex

%this should be changed for each week size
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

ROI_data(masks.mask_nodes_bottom_j_thresh_cortexnodes_allch) = 5;
%GETTTING CORTEX ZONES 
%F3  FZ  F4
%C3  CZ  C4
%P3  PZ  P4

%radius LM should be changed for each week size
radius_LM = 20;
radius_LM = 18;

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

%EEG_lm_9ROI = [];

figure()
plotmesh_iso2([gmSurfaceMesh.node ROI_data],gmSurfaceMesh.face)
clim([-5 5])
cb = colorbar('horiz');
ylabel(cb,'EEG LM Zone','FontSize',10,'Rotation',0)
title("9 EEG zones on cortex F3z4 C3z4 P3z4")

%find EEG 9LM that are within J thresh
%mask_J_threshold %1 or 0
%masks.mask_j_thresh_cortexnodes_allch %node idx
%landmark_EEG1020_cortexF3 %node idx

intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexF3);


%number of nodes overlapping between the 9EEG landmark ROI and nodes above
%J threshold from all channels (i.e best case scenario)
N_nodes_overlap_jthresh_all_9EEG = [size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexF3),1);
size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexFZ),1);
size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexF4),1);
size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexC3),1);
size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexCZ),1);
size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexC4),1);
size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexP3),1);
size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexPZ),1);
size(intersect(masks.mask_j_thresh_cortexnodes_allch, landmark_EEG1020_cortexP4),1)];


%for each event
for i=1:size(good_ch_idx_events,2)
    %number of nodes overlapping between the 9EEG landmark ROI and nodes above
    %J threshold from good channels idx
    N_nodes_overlap_jthresh_goodchidx_9EEG = [size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexF3),2);
    size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexFZ),2);
    size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexF4),2);
    size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexC3),2);
    size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexCZ),2);
    size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexC4),2);
    size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexP3),2);
    size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexPZ),2);
    size(intersect(find(mask_J_threshold(i,:)==1), landmark_EEG1020_cortexP4),2)];

    frac_N_nodes_jthresh_goodch_to_jthresh_all_9EEG = N_nodes_overlap_jthresh_goodchidx_9EEG./N_nodes_overlap_jthresh_all_9EEG;

    %if non of the fractions < 0.5 for any of the ROI, then consider it a good
    %event
    %was 0.5 , now set to 0.33 23 10 24
    if size(find(frac_N_nodes_jthresh_goodch_to_jthresh_all_9EEG<0.33),1)==0
        event_jac_threshold_idx(i) = 1;
    end

end


end