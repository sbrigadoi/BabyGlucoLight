function [] = t_val_metric(subjectN,eventType,savefig,stacked,mild_and_severe,time_window,PD,vol2gm,masks,gmSurfaceMesh,HbT_all_gmsurface,peak_glucose_all,good_PD_events,eventType_all,PD_data,good_PD_events_subN,BGC_metrics_PD_events ,HbT_metrics_PD_events,eventType_subN , eventType_subN_idx )


% run this - r and t analysis
if mild_and_severe ==  1
    eventType= "S_m_hypo";
end

%for stacking subjects
if stacked == 1; %yes we are stacking
  
            switch eventType
                case "m_hypo"
                        subjectN = 860;
                case "S_hypo"
                        subjectN = 1549;
            end
    
end
PD.subjectN = subjectN;
%HbT_metrics_PD_events.t_val
%HbT_metrics_PD_events.adj_p_t

t_val = HbT_metrics_PD_events.t_val;

pval_t =HbT_metrics_PD_events.adj_p_t';

t_val_pos_mask = zeros(size(t_val,1),size(t_val,2));
t_val_neg_mask = zeros(size(t_val,1),size(t_val,2));
t_val_zero_mask = zeros(size(t_val,1),size(t_val,2));

t_val_pos_neg_mask_net = zeros(size(t_val,1),1);

for i=1:size(t_val,2)
    t_val_pos = find(t_val(:,i)>0);
    t_val_neg = find(t_val(:,i)<0);
    t_val_zero = find(t_val(:,i)==0);


    t_val_pos_mask(t_val_pos,i) = 1;
    t_val_neg_mask(t_val_neg,i) = 1;
    t_val_zero_mask(t_val_zero,i) = 1;
end

%for each node, sum if mask is pos or neg, 
%pos value t_val_pos_neg_mask_net means more events HbT increase
%neg value t_val_pos_neg_mask_net means more events HbT decrease

%only consider nodes with p<0.05 for the tvalue


%find(pval_t(:,1)>0.05);

%set nodes which have P val (t) > 0.05 to zero
for i=1:size(t_val,2)
    t_val_pos_mask(find(pval_t(:,i)>0.05),i) = 0;
    t_val_neg_mask(find(pval_t(:,i)>0.05),i) = 0;
end

t_val_pos_mask_net = sum(t_val_pos_mask,2);
t_val_neg_mask_net = sum(t_val_neg_mask,2);

t_val_pos_neg_mask_net = t_val_pos_mask_net - t_val_neg_mask_net;

%t_val_pos_neg_mask_net_norm = t_val_pos_neg_mask_net/max(abs(t_val_pos_neg_mask_net));

t_val_pos_neg_mask_net_norm = t_val_pos_neg_mask_net/size(t_val,2);
%% plot
figure()
%subaxis(1,3,1,'sh',0.04) 
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
axis off
view(1.5,90); %TOP VIEW % used for infant week 30
clim([-1 1]);
%cb = colorbar('horiz');
%ylabel(cb,"\pm Net T-Val")
ax = gca;
ax.FontSize = 20;
sgtitle("PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+" TimeW "+time_window+" Net T-val N-events="+num2str(size(t_val,2))+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_T-val-net.png")
end








end