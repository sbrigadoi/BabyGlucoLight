function [] = R_val_2_metrics(subjectN,eventType,savefig,stacked,mild_and_severe,time_window,PD,vol2gm,masks,gmSurfaceMesh,HbT_all_gmsurface,peak_glucose_all,good_PD_events,eventType_all,PD_data,good_PD_events_subN,BGC_metrics_PD_events,HbT_metrics_PD_events,eventType_subN , eventType_subN_idx )

%%
% run this - r and t analysis
if mild_and_severe ==  1
    eventType= "S_m_hypo";
end

%for stacking subjects
if stacked == 1; %yes we are stacking
    switch time_window
        case "A"
            %subjectN = 820253349;
            %subjectN = 81015253334394144495556585960;
            switch eventType
                case "m_hypo"
                        subjectN = 860;
                case "S_hypo"
                        subjectN = 1549;
            end
        case "D"
            %subjectN = 81520253349;
            %subjectN = 81015253334394144495556585960;
           switch eventType
                case "m_hypo"
                        subjectN = 860;
                case "S_hypo"
                        subjectN = 1549;
            end

        case "F"
            %subjectN = 781520253349;
            %subjectN = 81015253334394144495556585960;
           switch eventType
                case "m_hypo"
                        subjectN = 860;
                case "S_hypo"
                        subjectN = 1549;
            end

    end
else
    good_PD_events_subN = [ones(1,size(good_PD_events,2))*subjectN];
end

PD.eventType = eventType;
PD.subjectN = subjectN;

%Mapping mask to gmSurface
mask_j_thresh_allnodes_allch_gmsurface = zeros(size(vol2gm,1),1);
mask_j_thresh_allnodes_allch_gmsurface(masks.mask_j_thresh_cortexnodes_allch,1)=1;

%top 95% JAC All Ch = 1, bottom 5% = 0;
% figure; plotmesh_iso2([gmSurfaceMesh.node mask_95_allnodes_allch_gmsurface],gmSurfaceMesh.face)
% cb = colorbar('horiz');

figure; plotmesh_iso2([gmSurfaceMesh.node mask_j_thresh_allnodes_allch_gmsurface],gmSurfaceMesh.face)
cb = colorbar('horiz');

% run this
rho = zeros(size(HbT_all_gmsurface,1),1);
pval = zeros(size(HbT_all_gmsurface,1),1);

%t_val = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
%h_val_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
%ci_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
%pval_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
%% HbT Max
metric_name= "HbT-max-vs-local-BCG_range";
BCG_metric = BGC_metrics_PD_events.local_BCG_range;
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-max-vs-global-BCG_range";
BCG_metric = BGC_metrics_PD_events.global_BCG_range;
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-max-vs-local-BCG_var";
BCG_metric = BGC_metrics_PD_events.local_BCG_var;
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-max-vs-global-BCG_var";
BCG_metric = BGC_metrics_PD_events.global_BCG_var;
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-max-vs-global-BCG-minima";
BCG_metric = BGC_metrics_PD_events.global_minima; 
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-max-vs-global-BCG-mean";
BCG_metric = BGC_metrics_PD_events.global_BCG_mean; 
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-max-vs-local-BCG-mean";
BCG_metric = BGC_metrics_PD_events.local_BCG_mean; 
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-max-vs-global-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.global_BCG_T_val; 
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

%B = rmmissing(A)
metric_name= "HbT-max-vs-local-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.local_BCG_T_val; 
HbT_metric = HbT_metrics_PD_events.HbT_all_max_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)
% HbT range
metric_name= "HbT-range-vs-local-BCG-range";
BCG_metric = BGC_metrics_PD_events.local_BCG_range;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-range-vs-global-BCG-range";
BCG_metric = BGC_metrics_PD_events.global_BCG_range;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-range-vs-local-BCG-mean";
BCG_metric = BGC_metrics_PD_events.local_BCG_mean;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-range-vs-global-BCG-mean";
BCG_metric = BGC_metrics_PD_events.global_BCG_mean;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-range-vs-local-BCG-var";
BCG_metric = BGC_metrics_PD_events.local_BCG_var;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-range-vs-global-BCG-var";
BCG_metric = BGC_metrics_PD_events.global_BCG_var;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-range-vs-local-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.local_BCG_T_val;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-range-vs-global-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.global_BCG_T_val;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-range-vs-global-BCG-minima";
BCG_metric = BGC_metrics_PD_events.global_minima;
HbT_metric = HbT_metrics_PD_events.range_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)


% HbT VAR
metric_name= "HbT-VAR-vs-global-BCG-range";
BCG_metric = BGC_metrics_PD_events.global_BCG_range;
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-VAR-vs-local-BCG-range";
BCG_metric = BGC_metrics_PD_events.local_BCG_range;
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-VAR-vs-global-BCG-var";
BCG_metric = BGC_metrics_PD_events.global_BCG_var;
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-var-vs-local-BCG-var";
BCG_metric = BGC_metrics_PD_events.local_BCG_var; 
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-VAR-vs-global-BCG-mean";
BCG_metric = BGC_metrics_PD_events.global_BCG_mean;
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-var-vs-local-BCG-mean";
BCG_metric = BGC_metrics_PD_events.local_BCG_mean; 
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-VAR-vs-global-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.global_BCG_T_val;
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-var-vs-local-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.local_BCG_T_val; 
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-VAR-vs-global-BCG-minima";
BCG_metric = BGC_metrics_PD_events.global_minima;
HbT_metric = HbT_metrics_PD_events.VAR_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

% HbT Mean
metric_name= "HbT-mean-vs-local-BCG-range";
BCG_metric = BGC_metrics_PD_events.local_BCG_range; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-mean-vs-global-BCG-range";
BCG_metric = BGC_metrics_PD_events.global_BCG_range; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-mean-vs-global-BCG-minima";
BCG_metric = BGC_metrics_PD_events.global_minima; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-mean-vs-global-BCG-mean";
BCG_metric = BGC_metrics_PD_events.global_BCG_mean; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-mean-vs-local-BCG-mean";
BCG_metric = BGC_metrics_PD_events.local_BCG_mean; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)


metric_name= "HbT-mean-vs-global-BCG-var";
BCG_metric = BGC_metrics_PD_events.global_BCG_var; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-mean-vs-local-BCG-var";
BCG_metric = BGC_metrics_PD_events.local_BCG_var; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-mean-vs-global-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.global_BCG_T_val; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-mean-vs-local-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.local_BCG_T_val; 
HbT_metric = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

% HbT Tval
metric_name= "HbT-Tval-vs-local-BCG-mean";
BCG_metric = BGC_metrics_PD_events.local_BCG_mean; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-Tval-vs-local-BCG-range";
BCG_metric = BGC_metrics_PD_events.local_BCG_range; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-Tval-vs-local-BCG-var";
BCG_metric = BGC_metrics_PD_events.local_BCG_var; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-Tval-vs-global-BCG-mean";
BCG_metric = BGC_metrics_PD_events.global_BCG_mean; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-Tval-vs-global-BCG-range";
BCG_metric = BGC_metrics_PD_events.global_BCG_range; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-Tval-vs-global-BCG-var";
BCG_metric = BGC_metrics_PD_events.global_BCG_var; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-Tval-vs-global-BCG-minima";
BCG_metric = BGC_metrics_PD_events.global_minima; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-Tval-vs-local-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.local_BCG_T_val; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)

metric_name= "HbT-Tval-vs-global-BCG-Tval";
BCG_metric = BGC_metrics_PD_events.global_BCG_T_val; 
HbT_metric = HbT_metrics_PD_events.t_val;
calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)


%%
% %% Get R val %one each for each node
% for i=1:size(HbT_all_gmsurface,1)
%     %Get R val
%     if size(HbT_all_gmsurface,3) > 2 %b fit only works with > 2 events 
% 
%         %[b, stats] = robustfit(peak_glucose_all,HbT_all_max(1,:)'   );
%         [b, stats] = robustfit(BGC_metrics_PD_events.local_BCG_range ,HbT_metrics_PD_events.HbT_all_max_gmsurface(i,:) );
% 
%         %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   b(1) + b(2)*peak_glucose_all     );
%         [rho(i),pval(i)] = corr(HbT_metrics_PD_events.HbT_all_max_gmsurface(i,:)',   b(1) + b(2)*BGC_metrics_PD_events.local_BCG_range);
% 
%     else
%         %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   peak_glucose_all     );
%         [rho(i),pval(i)] = corr(HbT_metrics_PD_events.HbT_all_max_gmsurface(i,:)',   BGC_metrics_PD_events.local_BCG_range    );
% 
%     end
% 
%     %length of glucose
%     %[b, stats] = robustfit(length_glucose_all_m_S_hypo,HbT_all_m_S_hypo_max(1,:)'   );
%     %[rho(i),pval(i)] = corr(HbT_all_m_S_hypo_max(1,:)',   b(1) + b(2)*length_glucose_all_m_S_hypo     );
%     %for j=1:size(HbT_all_gmsurface,3)
%     %    [t_test(i,j),pval_t(i,j)] = ttest(HbT_all_gmsurface(i,:,j));
%     %    [h_val_t(i,j),pval_t(i,j),ci_t,stats] = ttest(HbT_all_gmsurface(i,:,j));
%     %    t_val(i,j) = stats.tstat;
%     %end
% end 
% 
% %% FDR correct P vals and plot
% adj_p = ones(1,size(mask_j_thresh_allnodes_allch_gmsurface,1));
% %adj_p_t = ones(size(HbT_all_gmsurface,3),size(mask_j_thresh_allnodes_allch_gmsurface,1));
% 
% %adj_p_mask_t = ones(size(mask_95_allnodes_allch_gmsurface,1),size(HbT_all_m_hypo_gmsurface,3));
% % FDR correction
% [h, crit_p, adj_ci_cvrg, adj_p_mask]=fdr_bh(pval(logical(mask_j_thresh_allnodes_allch_gmsurface)));
% 
% %set nodes in top 95% mask to be corrected p vals.
% adj_p(logical(mask_j_thresh_allnodes_allch_gmsurface)) = adj_p_mask;
% 
% %max(rho)
% %set R P to be zero for nodes outside of mask.
% %rho(masks.mask_nodes_bottom_5_cortexnodes_allch)=0;
% %pval(masks.mask_nodes_bottom_5_cortexnodes_allch)=-1;
% 
% rho(~logical(mask_j_thresh_allnodes_allch_gmsurface))=0;
% pval(~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;
% adj_p(~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;
% 
% if size(HbT_all_gmsurface,3) > 1 %max rho > 1 events 
%     max_rho = max([max(rho) abs(min(rho))]);
% else
%     max_rho = 1;
% end
% 
% max_pval = max([max(pval) abs(min(pval))]);
% max_adj_p = max([max(adj_p) abs(min(adj_p))]);
% 
% % run this plot mesh (iso2mesh)
% %rho=rho';
% %pval=pval';
% figure()
% subplot(1,3,1)
% plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
% view(0,90); %TOP VIEW % used for infant week 30
% clim([-max_rho max_rho]);
% cb = colorbar('horiz');
% ylabel(cb,"R value")
% title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
% subplot(1,3,2)
% plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.1 0.1]);
% cb = colorbar('horiz');
% ylabel(cb,"P value")
% title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% subplot(1,3,3)
% plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
% view(0,90);
% clim([-max_adj_p max_adj_p]);
% cb = colorbar('horiz');
% ylabel(cb,"FDR cor. P value")
% title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% 
% figure()
% subplot(1,3,1)
% plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
% view(0,90); %TOP VIEW % used for infant week 30
% %clim([-max_rho max_rho]);
% clim([-1 1]);
% cb = colorbar('horiz');
% ylabel(cb,"R value")
% title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
% subplot(1,3,2)
% plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.1 0.1]);
% cb = colorbar('horiz');
% ylabel(cb,"P value")
% title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% subplot(1,3,3)
% plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
% view(0,90);
% clim([-max_adj_p max_adj_p]);
% cb = colorbar('horiz');
% ylabel(cb,"FDR cor. P value")
% title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% if savefig == 1
%     saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_val.png")
% end
% 
% figure()
% subplot(1,3,1)
% plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
% view(0,90); %TOP VIEW % used for infant week 30
% %clim([-max_rho max_rho]);
% clim([-1 1]);
% cb = colorbar('horiz');
% ylabel(cb,"R value")
% title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
% ax = gca;
% ax.FontSize = 10;
% subplot(1,3,2)
% plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.05 0.05]);
% cb = colorbar('horiz');
% ylabel(cb,"P value")
% title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% ax = gca;
% ax.FontSize = 10;
% subplot(1,3,3)
% plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
% view(0,90);
% clim([-max_adj_p max_adj_p]);
% cb = colorbar('horiz');
% ylabel(cb,"FDR cor. P value")
% title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% if savefig == 1
%     saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_val_p0_05.png")
% end
% ax = gca;
% ax.FontSize = 10;




end