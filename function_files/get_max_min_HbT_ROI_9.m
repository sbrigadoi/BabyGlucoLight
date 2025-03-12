function [max_max_HbT_ROI_9] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexAB)

%get min and max vlues and reshapre in 2D matrix
max_HbT_ROI_9 = max(max(HbT_all_gmsurface(landmark_EEG1020_cortexAB,121:end,:)));
min_HbT_ROI_9 = min(min(HbT_all_gmsurface(landmark_EEG1020_cortexAB,121:end,:)));

max_HbT_ROI_9 = reshape(max_HbT_ROI_9(1,1,:),1,[]);
min_HbT_ROI_9 = reshape(min_HbT_ROI_9(1,1,:),1,[]);

%see if the max or min value has bigger magnitude
idx_max = abs(max_HbT_ROI_9)>=abs(min_HbT_ROI_9);
idx_min = abs(max_HbT_ROI_9)<abs(min_HbT_ROI_9);

max_max_HbT_ROI_9(idx_max) = max_HbT_ROI_9(idx_max);
max_max_HbT_ROI_9(idx_min) = min_HbT_ROI_9(idx_min);

end