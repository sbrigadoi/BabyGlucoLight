function [] = L_R_hemisphere_meta_data(HbT_metrics_PD_events,PD_meta_data,good_PD_events_subN,good_PD_events)


%find(PD_meta_data.PD_sex(:,2)==0)
%%

n_female = size(find(PD_meta_data.PD_sex(good_PD_events_subN',2)==0),1);
n_male = size(find(PD_meta_data.PD_sex(good_PD_events_subN',2)==1),1);

n_unblind = size(find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==0),1);
n_blind = size(find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==1),1);



UB_HbT_Hemisphere = [HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==0) ) HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==0) )]; 
Blinded_HbT_Hemisphere = [HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==1) ) HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==1) )];

[h,p,ci,stats] = ttest2(UB_HbT_Hemisphere,Blinded_HbT_Hemisphere);

UB_HbT_gmsurface = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface(:,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==0));
B_HbT_gmsurface = HbT_metrics_PD_events.mean_change_HbT_all_gmsurface(:,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==1));


[h,p,ci,stats] = ttest2(mean(UB_HbT_gmsurface'),mean(B_HbT_gmsurface'));

% FDR correction
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p);

%%
figure()
subplot(2,2,1)
sgtitle("PD m/s hypo only - no hyper events")
plot(zeros(1,n_female),HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,find(PD_meta_data.PD_sex(good_PD_events_subN',2)==0) ),'rx')
hold on
plot(ones(1,n_female)*0.5,HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,find(PD_meta_data.PD_sex(good_PD_events_subN',2)==0) ),'rx')
plot(ones(1,n_male),HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,find(PD_meta_data.PD_sex(good_PD_events_subN',2)==1) ),'bx')
plot(ones(1,n_male)*1.5,HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,find(PD_meta_data.PD_sex(good_PD_events_subN',2)==1) ),'bx')
xlabel('LHS F - RHS F - LHS M - RHS M')
ylabel('Mean hemisphere change in HbT')
xlim([-0.5 2])
legend('LHS F','RHS F','LHS M','RHS M')
title("Sex")

subplot(2,2,2)
plot(zeros(1,n_unblind),HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==0) ),'rx')
hold on
plot(ones(1,n_unblind)*0.5,HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==0) ),'rx')
plot(ones(1,n_blind),HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==1) ),'bx')
plot(ones(1,n_blind)*1.5,HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==1) ),'bx')
xlabel('LHS UB - RHS UB - LHS B - RHS B')
ylabel('Mean hemisphere change in HbT')
xlim([-0.5 2])
legend('LHS UB','RHS UB','LHS B','RHS B')
title("Blinded/Unblinded")

subplot(2,2,3)
plot(PD_meta_data.PD_bw(good_PD_events_subN',2),HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,:),'rx')
hold on
plot(PD_meta_data.PD_bw(good_PD_events_subN',2),HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,:),'bx')
ylabel('Mean hemisphere change in HbT / \mu M')
xlabel('Birth weight / g')
legend('LHS','RHS')
title("Birth Weight")

subplot(2,2,4)
plot(PD_meta_data.PD_ga(good_PD_events_subN',2),HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,:),'rx')
hold on
plot(PD_meta_data.PD_ga(good_PD_events_subN',2),HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,:),'bx')
ylabel('Mean hemisphere change in HbT / \mu M')
xlabel('gest. age / days')
legend('LHS','RHS')
title("Gest. Age")



figure()
plot(ones(1,n_unblind)*0.5,HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==0) ),'ro')
hold on
plot(ones(1,n_unblind)*0.5,HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==0) ),'rx')
plot(ones(1,n_blind)*1.5,HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(1,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==1) ),'bo')
plot(ones(1,n_blind)*1.5,HbT_metrics_PD_events.hemisphere_HbT_all_gmsurface(2,find(PD_meta_data.PD_blinded(good_PD_events_subN',2)==1) ),'bx')
xlabel('LHS UB - RHS UB - LHS B - RHS B')
ylabel('Mean hemisphere change in HbT')
xlim([-0.5 2])
legend('LHS UB','RHS UB','LHS B','RHS B')
title("Blinded/Unblinded")






end