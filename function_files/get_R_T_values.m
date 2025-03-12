function [] = get_R_T_values(subjectN,eventType,savefig,stacked,mild_and_severe,time_window,PD,vol2gm,masks,gmSurfaceMesh,HbT_all_gmsurface,peak_glucose_all,good_PD_events,eventType_all,PD_data,good_PD_events_subN);

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
            subjectN = 860;
        case "D"
            %subjectN = 81520253349;
            %subjectN = 81015253334394144495556585960;
           subjectN = 860;

        case "F"
            %subjectN = 781520253349;
            %subjectN = 81015253334394144495556585960;
           subjectN = 860;

    end
else
    good_PD_events_subN = [ones(1,size(good_PD_events,2))*subjectN];
end



%subjectN = 820253349;
%subjectN = 81520253349;
PD.eventType = eventType;
PD.subjectN = subjectN;

%Mapping mask to gmSurface
% mask_95_allnodes_allch = zeros(size(vol2gm,2),1);
% mask_95_allnodes_allch(cortex_nodes(masks.mask_95_cortexnodes_allch)) = 1;
% mask_95_allnodes_allch_gmsurface = vol2gm*mask_95_allnodes_allch;
% 
% mask_95_allnodes_allch_gmsurface(find(mask_95_allnodes_allch_gmsurface~=0))=1;

%mask_j_thresh_cortexnodes_allch
%masks.mask_nodes_bottom_j_thresh_cortexnodes_allch

mask_j_thresh_allnodes_allch_gmsurface = zeros(size(vol2gm,1),1);
mask_j_thresh_allnodes_allch_gmsurface(masks.mask_j_thresh_cortexnodes_allch,1)=1;

%top 95% JAC All Ch = 1, bottom 5% = 0;
% figure; plotmesh_iso2([gmSurfaceMesh.node mask_95_allnodes_allch_gmsurface],gmSurfaceMesh.face)
% cb = colorbar('horiz');

figure; plotmesh_iso2([gmSurfaceMesh.node mask_j_thresh_allnodes_allch_gmsurface],gmSurfaceMesh.face)
cb = colorbar('horiz');

% nnz(mask_95_allnodes_allch_gmsurface)
% 
% cortex_nodes(masks.mask_95_cortexnodes_allch);
% cortex_nodes;

% run this
rho = zeros(size(HbT_all_gmsurface,1),1);
pval = zeros(size(HbT_all_gmsurface,1),1);

t_val = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
h_val_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
ci_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
pval_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));

% Max-min BGC vs max HbT
%[BGC_values_PD_events] = get_timewindow_10mins_BGCvalues(good_PD_events_subN,good_PD_events,time_window);


% for i=1:size(HbT_all_gmsurface,1)
%     %Max HbT
%     [HbT_all_max HbT_all_max_idx] = (max(abs(HbT_all_gmsurface(i,121:end,:)))); %get max HbT from T=2:10mins
%     HbT_all_max_idx = reshape(HbT_all_max_idx(1,1,:), 1, []);
% 
%     HbT_all_max = HbT_all_gmsurface(i,121+HbT_all_max_idx-1,:); %get max HbT from T=2:10mins
%     HbT_all_max = diag(reshape(HbT_all_max,[],size(HbT_all_max,3))) ;%rehsape into sqaure, and get diag
% 
%     %max - min glucose
%     if size(peak_glucose_all,1) > 2 %b fit only works with > 2 events 
% 
%         %[b, stats] = robustfit(peak_glucose_all,HbT_all_max(1,:)'   );
%         [b, stats] = robustfit(peak_glucose_all,HbT_all_max   );
% 
%         %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   b(1) + b(2)*peak_glucose_all     );
%         [rho(i),pval(i)] = corr(HbT_all_max,   b(1) + b(2)*peak_glucose_all     );
% 
%     else
%         %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   peak_glucose_all     );
%         [rho(i),pval(i)] = corr(HbT_all_max,   peak_glucose_all     );
% 
%     end
% 
% 
% end

% R val. Peak BCG vs max HbT
%for i=1:size(cortex_nodes,1)
for i=1:size(HbT_all_gmsurface,1)

    %HbT_all_S_hypo_max = (max(HbT_all_S_hypo(i,121:end,:)));    
    %HbT_all_max = (max(HbT_all(i,121:end,:))); %get max HbT from T=2:10mins
    %used pre 31 05 24
    %HbT_all_max = (max(HbT_all_m_hypo_gmsurface(i,121:end,:))); %get max HbT from T=2:10mins
    %used 31 05 24
    %HbT_all_max = (max(abs(HbT_all_gmsurface(i,121:end,:)))); %get max HbT from T=2:10mins


    [HbT_all_max HbT_all_max_idx] = (max(abs(HbT_all_gmsurface(i,121:end,:)))); %get max HbT from T=2:10mins
    HbT_all_max_idx = reshape(HbT_all_max_idx(1,1,:), 1, []);
    
    HbT_all_max = HbT_all_gmsurface(i,121+HbT_all_max_idx-1,:); %get max HbT from T=2:10mins
    HbT_all_max = diag(reshape(HbT_all_max,[],size(HbT_all_max,3))) ;%rehsape into sqaure, and get diag

  
    %[a b]= max(.....) %one output will be index, one output will be the
    %value (from input i.e 121:end)
    %find index of abs max
    %but keep information of + or - value
    % updated 26 08 24

    %get b and rho values for each node
    %[b, stats] = robustfit(peak_glucose_all_S_hypo,HbT_all_S_hypo_max(1,:)'   );
    %[rho(i),pval(i)] = corr(HbT_all_S_hypo_max(1,:)',   b(1) + b(2)*peak_glucose_all_S_hypo     );
    
    %peak glucose
    if size(peak_glucose_all,1) > 2 %b fit only works with > 2 events 

        %[b, stats] = robustfit(peak_glucose_all,HbT_all_max(1,:)'   );
        [b, stats] = robustfit(peak_glucose_all,HbT_all_max   );

        %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   b(1) + b(2)*peak_glucose_all     );
        [rho(i),pval(i)] = corr(HbT_all_max,   b(1) + b(2)*peak_glucose_all     );

    else
        %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   peak_glucose_all     );
        [rho(i),pval(i)] = corr(HbT_all_max,   peak_glucose_all     );

    end

    %length of glucose
    %[b, stats] = robustfit(length_glucose_all_m_S_hypo,HbT_all_m_S_hypo_max(1,:)'   );
    %[rho(i),pval(i)] = corr(HbT_all_m_S_hypo_max(1,:)',   b(1) + b(2)*length_glucose_all_m_S_hypo     );
    for j=1:size(HbT_all_gmsurface,3)
        [t_test(i,j),pval_t(i,j)] = ttest(HbT_all_gmsurface(i,:,j));
        [h_val_t(i,j),pval_t(i,j),ci_t,stats] = ttest(HbT_all_gmsurface(i,:,j));
        t_val(i,j) = stats.tstat;
    end
end 

%[h,p,ci,stats] = ttest(HbT_all_m_hypo_gmsurface(1,:,1));
%t_val = stats.tstat
%
% run this
%pval(logical(mask_95_allnodes_allch_gmsurface))
% set up adj p matrix of ones
%mask_j_thresh_allnodes_allch_gmsurface

adj_p = ones(1,size(mask_j_thresh_allnodes_allch_gmsurface,1));
adj_p_t = ones(size(HbT_all_gmsurface,3),size(mask_j_thresh_allnodes_allch_gmsurface,1));

%adj_p_mask_t = ones(size(mask_95_allnodes_allch_gmsurface,1),size(HbT_all_m_hypo_gmsurface,3));
% FDR correction
[h, crit_p, adj_ci_cvrg, adj_p_mask]=fdr_bh(pval(logical(mask_j_thresh_allnodes_allch_gmsurface)));

%fix this
for i=1:size(HbT_all_gmsurface,3)
    [h_t, crit_p_t, adj_ci_cvrg_t, adj_p_mask_t(:,i)]=fdr_bh(pval_t(logical(mask_j_thresh_allnodes_allch_gmsurface)));
end

%set nodes in top 95% mask to be corrected p vals.
adj_p(logical(mask_j_thresh_allnodes_allch_gmsurface)) = adj_p_mask;

for i=1:size(adj_p_mask_t,2)
    adj_p_t(i,logical(mask_j_thresh_allnodes_allch_gmsurface)) = adj_p_mask_t(:,i);
end

%max(rho)
%set R P to be zero for nodes outside of mask.
%rho(masks.mask_nodes_bottom_5_cortexnodes_allch)=0;
%pval(masks.mask_nodes_bottom_5_cortexnodes_allch)=-1;

rho(~logical(mask_j_thresh_allnodes_allch_gmsurface))=0;
pval(~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;
adj_p(~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;

t_val(~logical(mask_j_thresh_allnodes_allch_gmsurface),:)=0;
pval_t(~logical(mask_j_thresh_allnodes_allch_gmsurface),:)=-1;
adj_p_t(:,~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;

if size(peak_glucose_all,1) > 1 %max rho > 1 events 
    max_rho = max([max(rho) abs(min(rho))]);
else
    max_rho = 1;
end

max_pval = max([max(pval) abs(min(pval))]);
max_adj_p = max([max(adj_p) abs(min(adj_p))]);

max_t_val = max([max(t_val) abs(min(t_val))]);
max_pval_t = max([max(pval_t) abs(min(pval_t))]);
max_adj_p_t = max([max(adj_p_t) abs(min(adj_p_t))]);

%% 
% run this plot mesh (iso2mesh)
%rho=rho';
%pval=pval';
figure()
subplot(1,3,1)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(0,90); %TOP VIEW % used for infant week 30
clim([-max_rho max_rho]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
subplot(1,3,2)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(0,90);
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(1,3,3)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(0,90);
clim([-max_adj_p max_adj_p]);
cb = colorbar('horiz');
ylabel(cb,"FDR cor. P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")

figure()
subplot(1,3,1)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(0,90); %TOP VIEW % used for infant week 30
%clim([-max_rho max_rho]);
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
subplot(1,3,2)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(0,90);
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(1,3,3)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(0,90);
clim([-max_adj_p max_adj_p]);
cb = colorbar('horiz');
ylabel(cb,"FDR cor. P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_val.png")
end

figure()
subplot(1,3,1)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(0,90); %TOP VIEW % used for infant week 30
%clim([-max_rho max_rho]);
clim([-0.05 0.05]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
ax = gca;
ax.FontSize = 10;
subplot(1,3,2)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(0,90);
clim([-0.05 0.05]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
ax = gca;
ax.FontSize = 10;
subplot(1,3,3)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(0,90);
clim([-max_adj_p max_adj_p]);
cb = colorbar('horiz');
ylabel(cb,"FDR cor. P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_val_p0_05.png")
end
ax = gca;
ax.FontSize = 10;



%average across p<0.05 nodes to get scatter plot
p_val_0_05_R = intersect(find(pval<=0.05),find(pval>=0));
max_HbT_p_val_0_05_R_all = zeros(size(p_val_0_05_R,1),size(good_PD_events,2));
%max_HbT_p_val_0_05_R;
%peak_glucose_all
for i=1:size(p_val_0_05_R,1)
    %Abs HbT
    %max_HbT_p_val_0_05_R(i,:) = max(abs(HbT_all_gmsurface(p_val_0_05_R(i),121:end,:)));

    [max_HbT_p_val_0_05_R max_HbT_p_val_0_05_R_idx] = (max(abs(HbT_all_gmsurface(p_val_0_05_R(i),121:end,:)))); %get max HbT from T=2:10mins
    max_HbT_p_val_0_05_R_idx = reshape(max_HbT_p_val_0_05_R_idx(1,1,:), 1, []);
    
    max_HbT_p_val_0_05_R = HbT_all_gmsurface(p_val_0_05_R(i),121+max_HbT_p_val_0_05_R_idx-1,:); %get max HbT from T=2:10mins
    max_HbT_p_val_0_05_R = diag(reshape(max_HbT_p_val_0_05_R,[],size(max_HbT_p_val_0_05_R,3)))';%rehsape into sqaure, and get diag
    max_HbT_p_val_0_05_R_all(i,:) = max_HbT_p_val_0_05_R;

end

c_max_HbT_p_val_0_05_R = polyfit(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:)),1);
fit_y_c_p_glucose_max_r_HbT = c_max_HbT_p_val_0_05_R(1)*peak_glucose_all + c_max_HbT_p_val_0_05_R(2);


figure()
plot(peak_glucose_all,max_HbT_p_val_0_05_R_all(:,:),'bo')
hold on
plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'r--')
xlabel("Peak Glucose Value / mg/dL")
ylabel("Max HbT / \mu M")
title("Combined plot R of P<0.05- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" N CORTEX NODES "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_p_0_05nodes_scatter.png")
end

%c_max_HbT_p_val_0_05_R = polyfit(peak_glucose_all,mean(max_HbT_p_val_0_05_R(:,:)),1);
%fit_y_c_p_glucose_max_r_HbT = c_max_HbT_p_val_0_05_R(1)*peak_glucose_all + c_max_HbT_p_val_0_05_R(2);


figure()
plot(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:)),'bo')
hold on
plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'r--')
xlabel("Peak Glucose Value / mg/dL")
ylabel("Max HbT / \mu M")
title("Mean Combined plot R of P<0.05- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" N CORTEX NODES "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_p_0_05nodes_scatter_mean.png")
end
% testing coloured plot/scatter
%figure()
%plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'r--','LineWidth',1)
%hold on
%figure()
%hold on
%for i=1:size(good_PD_events_subN,2)
%get colour
%good_PD_events_subN(i);  %'MarkerFaceColor',[0.5,0.5,0.5]
%plot(peak_glucose_all(i),mean(max_HbT_p_val_0_05_R_all(:,i)),'x','LineWidth',7,'MarkerFaceColor',[good_PD_events_subN(i)/100,good_PD_events_subN(i)/100, good_PD_events_subN(i)/100]);
%s1 = plot(peak_glucose_all(i),mean(max_HbT_p_val_0_05_R_all(:,i)),'x','LineWidth',7);
%s1.MarkerFaceColor = [good_PD_events_subN(i)/100 good_PD_events_subN(i)/100 good_PD_events_subN(i)/100];
%end
%xlabel("Peak Glucose Value / mg/dL")
%ylabel("Max HbT / \mu M")
%title("Mean Combined plot R of P<0.05- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" N CORTEX NODES "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
%testing colour scatter options
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%legend('fit',split(num2str(good_PD_events_subN)))
%legend('x','c')
%legend(num2str(good_PD_events_subN))
%newStr = split(num2str(good_PD_events_subN))
%scatter plot colours

%figure()
%plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'r--','LineWidth',1)
%hold on

%1:1:size(unique(good_PD_events_subN'),1)
%good_PD_events_subN
%rng(2,"twister")
%seed = rng;

%RGB_col = [102/255 197/255 204/255;
%246/255 207/255 113/255;
%248/255 156/255 116/255;
%220/255 176/255 242/255;
%135/255 197/255 95/255;
%179/255 179/255 179/255];

RGB_col = [241/255 71/255 37/255;
241/255 245/255 33/255;
17/255 251/255 56/255;
255/255 144/255 37/255;
236/255 137/255 132/255;
202/255 125/255 187/255;
136/255 8/255 165/255;
63/255 73/255 254/255;
16/255 81/255 165/255;
13/255 8/255 241/255;
10/255 18/255 202/255;
54/255 208/255 182/255;
67/255 34/255 10/255;
200/255 168/255 175/255;
188/255 143/255 16/255;
251/255 1/255 165/255;];

%col_scat = [1.6*good_PD_events_subN'/100, 1.5*good_PD_events_subN'/100, 0.8*good_PD_events_subN'/100];

col_scat_idx = zeros(size(good_PD_events_subN,2),1);
col_scat = zeros(size(good_PD_events_subN,2),1);
for i=1:size(good_PD_events_subN,2)
    col_scat_idx(i,1) = find(good_PD_events_subN(i) == unique(good_PD_events_subN));
end

[b, stats] = robustfit(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:))   );

%[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   b(1) + b(2)*peak_glucose_all     );
[rho_HbT_pval_0_05_R_all,pvalHbT_pval_0_05_R_all] = corr(mean(max_HbT_p_val_0_05_R_all(:,:))',   b(1) + b(2)*peak_glucose_all     );

%plot(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:)),'bo')
%
%timeWA psoter fNIRS 04 09 24


%%
if stacked == 1
    switch time_window
        case "A"
            col_scat_idx(1:3,1) = 1;
            col_scat_idx(4,1) = 2;
            col_scat_idx(5,1) = 3;
            col_scat_idx(6:8,1) = 4;
            col_scat_idx(9:19,1) = 6;
            col_scat_idx(20,1) = 7;
            col_scat_idx(21,1) = 8;
            col_scat_idx(22,1) = 9;
            col_scat_idx(23:36,1) = 10;
            col_scat_idx(37:38,1) = 11;
            col_scat_idx(39,1) = 12;
            col_scat_idx(40,1) = 14;
            col_scat_idx(41:44,1) = 15;
            %time w A all fnirs 24 poster
            % col_scat_idx(1:3,1) = 2;
            % col_scat_idx(4:8,1) = 4;
            % col_scat_idx(9:11,1) = 5;
            % col_scat_idx(12:22,1) = 6;
            % col_scat_idx(23:31,1) = 7;

        case "D"
            col_scat_idx(1:3,1) = 1;
            col_scat_idx(4,1) = 2;
            col_scat_idx(5:6,1) = 3;
            col_scat_idx(7,1) = 4;
            col_scat_idx(8:20,1) = 5;
            col_scat_idx(21,1) = 6;
            col_scat_idx(22,1) = 7;
            col_scat_idx(23,1) = 8;
            col_scat_idx(24,1) = 9;
            col_scat_idx(25:39,1) = 10;
            col_scat_idx(40,1) = 11;
            col_scat_idx(41:42,1) = 12;
            col_scat_idx(43,1) = 13;
            col_scat_idx(44:45,1) = 14;
            col_scat_idx(46:48,1) = 15;
            %time W D
            % col_scat_idx(1:3,1) = 2;
            % col_scat_idx(4:5,1) = 3;
            % col_scat_idx(6:13,1) = 4;
            % col_scat_idx(14,1) = 5;
            % col_scat_idx(15:27,1) = 6;
            % col_scat_idx(28:38,1) = 7;

        case "F"
            col_scat_idx(1:5,1) = 1;
            col_scat_idx(5,1) = 2;
            col_scat_idx(7:8,1) = 3;
            col_scat_idx(9:10,1) = 4;
            col_scat_idx(11:20,1) = 5;
            
            col_scat_idx(21,1) = 7;
            col_scat_idx(22,1) = 8;
            col_scat_idx(23,1) = 9;
            col_scat_idx(24:38,1) = 10;
            col_scat_idx(39:40,1) = 11;
            col_scat_idx(41:42,1) = 12;
            col_scat_idx(43,1) = 13;
            col_scat_idx(44:45,1) = 14;
            col_scat_idx(46:50,1) = 15;
            %time W F
            %7
            % col_scat_idx(1,1) = 1;
            % col_scat_idx(2:6,1) = 2;
            % col_scat_idx(7:8,1) = 3;
            % col_scat_idx(9:13,1) = 4;
            % col_scat_idx(14:15,1) = 5;
            % col_scat_idx(16:25,1) = 6;
            % col_scat_idx(26:37,1) = 7;
    end
end

col_scat = RGB_col(col_scat_idx(:,1),:);

%savefig=1;

figure()
plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'k-','LineWidth',8)
hold on
scatter(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:)),130,col_scat,'filled')
xlim([40 72])
ylim([-13 13])
xlabel("Peak Glucose Value / mg/dL")
ylabel("Max HbT / \mu M")
title("Mean Combined R of P<0.05- Glc vs Max HbT. R="+num2str(rho_HbT_pval_0_05_R_all)+" P="+num2str(pvalHbT_pval_0_05_R_all)+" Sub. "+num2str(subjectN)+" "+eventType+" N CORTEX NODES "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
ax = gca;
ax.FontSize = 16;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_p_0_05nodes_scatter_mean_colour.png")
end

%savefig=0;

figure()
subplot(2,3,1)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(90,0); %TOP VIEW % used for infant week 30
%clim([-max_rho max_rho]);
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
subplot(2,3,4)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(-90,0);
%clim([-max_rho max_rho]);
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
subplot(2,3,2)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(90,0); %TOP VIEW % used for infant week 30
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(2,3,5)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(-90,0);
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(2,3,3)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(90,0); %TOP VIEW % used for infant week 30
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(2,3,6)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(-90,0);
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_val_side.png")
end

% run this
%find(rho==max(rho));
%pval(find(rho==max(rho)));

%max_rho_HbT = max(HbT_all(find(rho==max(rho)),121:end,:));
%max_rho_HbT = max(HbT_all_m_hypo_gmsurface(find(rho==max(rho)),121:end,:));

%max rho HbT using max abs HbT all pre 26 08 24
%max_rho_HbT = max(abs(HbT_all_gmsurface(find(rho==max(rho)),121:end,:)));

%max rho HbT using + and - polarity
%max_rho_HbT = HbT_all_gmsurface(find(rho==max(rho)),121:end,:);

[max_rho_HbT max_rho_HbT_idx] = (max(abs(HbT_all_gmsurface(find(rho==max(rho)),121:end,:)))); %get max HbT from T=2:10mins
max_rho_HbT_idx = reshape(max_rho_HbT_idx(1,1,:), 1, []);
    
max_rho_HbT = HbT_all_gmsurface(find(rho==max(rho)),121+max_rho_HbT_idx-1,:); %get max HbT from T=2:10mins
max_rho_HbT = diag(reshape(max_rho_HbT,[],size(max_rho_HbT,3)));%rehsape into sqaure, and get diag

%max_rho_HbT = max(HbT_all(33044,121:end,:));
%max_rho_HbT = max_rho_HbT(1,:)'; %pre 26 08 24 
%peak_glucose_all;

switch subjectN
     case 7
        figure()
        plot(peak_glucose_all,max_rho_HbT,'x')
        legend('m hypo')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+"")
 
     case 8
        figure()
        plot(peak_glucose_all,max_rho_HbT,'x')
        legend('m hypo')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+"")
 
    case 20
        figure()
        plot(peak_glucose_all,max_rho_HbT,'x')
        legend('m hypo')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+"")
    case 33
        figure()
        plot(peak_glucose_all(1:end-2),max_rho_HbT(1:end-2),'x')
        hold on
        plot(peak_glucose_all(end-1:end),max_rho_HbT(end-1:end),'rx')
        legend('m hypo','S hypo')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+"")
end

c_p_glucose_max_r_HbT = polyfit(peak_glucose_all,max_rho_HbT,1);
fit_y_c_p_glucose_max_r_HbT = c_p_glucose_max_r_HbT(1)*peak_glucose_all + c_p_glucose_max_r_HbT(2);

figure()
plot(peak_glucose_all,max_rho_HbT,'x')
hold on
        plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'o--')
        legend('m hypo','p1 fit','Location','northwest')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+" p1 "+num2str(c_p_glucose_max_r_HbT(1))+" p0 "+num2str(c_p_glucose_max_r_HbT(2))+" ")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_maxR_node.png")
end

%savefig=1;
%colour
figure()
scatter(peak_glucose_all,max_rho_HbT,65,col_scat,'filled')
hold on
plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'k-','LineWidth',1)
xlim([40 72])
ylim([-13 13])
%legend('m hypo','p1 fit','Location','northwest')
xlabel("Peak Glucose Value / mg/dL")
ylabel("Max HbT / \mu M")
title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+" p1 "+num2str(c_p_glucose_max_r_HbT(1))+" p0 "+num2str(c_p_glucose_max_r_HbT(2))+" ")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_maxR_node_colour.png")
end
%savefig=0;
%
%t test, t val.
% switch subjectN
%     case 7
%         eventN_t_test = good_PD_events;
%     case 8
%         eventN_t_test = good_PD_events;
%     case 20
%         eventN_t_test = good_PD_events;
%     case 25
%         eventN_t_test = good_PD_events;
%     case 33
%         eventN_t_test = good_PD_events;
% end
%% 
% t val stuff below
eventN_t_test = good_PD_events;

for i=1:size(t_val,2)
    figure()
    subplot(1,3,1)
    plotmesh_iso2([gmSurfaceMesh.node t_val(:,i)],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    clim([-max_t_val max_t_val]);
    cb = colorbar('horiz');
    ylabel(cb,"T value")
    title("T val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" "+num2str(eventN_t_test(i))+" TW "+time_window+"")
    subplot(1,3,2)
    plotmesh_iso2([gmSurfaceMesh.node pval_t(:,i)],gmSurfaceMesh.face)
    view(0,90);
    clim([-0.1 0.1]);
    cb = colorbar('horiz');
    ylabel(cb,"P value")
    title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" "+num2str(eventN_t_test(i))+" TW "+time_window+"")
    subplot(1,3,3)
    plotmesh_iso2([gmSurfaceMesh.node adj_p_t(i,:)'],gmSurfaceMesh.face)
    view(0,90);
    %clim([-max_adj_p_t max_adj_p_t]);
    clim([-0.1 0.1]);
    cb = colorbar('horiz');
    ylabel(cb,"FDR cor. P value")
    title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" "+num2str(eventN_t_test(i))+" TW "+time_window+"")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val"+num2str(eventN_t_test(i))+".png")
    end
end

% t test DICE
%t_val
%mask_95_allnodes_allch_gmsurface

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

figure()
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_mask(:,1)],gmSurfaceMesh.face)
cb = colorbar('horiz');

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

%size(t_val,2)
figure()
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(0,90); %TOP VIEW % used for infant week 30
title("T val NET (p<0.05) (N Pos - N Neg), PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_net.png")
end

figure()
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(0,90); %TOP VIEW % used for infant week 30
title("T val NET Norm (p<0.05) (N Pos - N Neg)/N Events, PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_net_norm.png")
end

figure()
subplot(2,2,1)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(-90,0); %TOP VIEW % used for infant week 30
subplot(2,2,2)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(90,0); %TOP VIEW % used for infant week 30
subplot(2,2,3)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(0,0); %TOP VIEW % used for infant week 30
subplot(2,2,4)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(180,0); %TOP VIEW % used for infant week 30
sgtitle("T val NET (p<0.05) (N Pos - N Neg), PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_net_side.png")
end

figure()
subplot(2,2,1)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(-90,0); %TOP VIEW % used for infant week 30
subplot(2,2,2)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(90,0); %TOP VIEW % used for infant week 30
subplot(2,2,3)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(0,0); %TOP VIEW % used for infant week 30
subplot(2,2,4)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(180,0); %TOP VIEW % used for infant week 30
sgtitle("T val NET Norm (p<0.05) (N Pos - N Neg)/N Events, PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_net_norm_side.png")
end

for i=1:size(t_val,2)
    for j=1:size(t_val,2)
        t_val_dice_pos(i,j) = dice(t_val_pos_mask(:,i),t_val_pos_mask(:,j));
        t_val_dice_neg(i,j) = dice(t_val_neg_mask(:,i),t_val_neg_mask(:,j));
    end
end

figure()
subplot(1,2,1)
imagesc(t_val_dice_pos); 
xlabel("Event N")
xlabel("Event N")
title("T val DICE (Pos), PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
clim([0 1]);
colorbar;
subplot(1,2,2)
imagesc(t_val_dice_neg); 
xlabel("Event N")
xlabel("Event N")
title("T val DICE (Neg), PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
clim([0 1]);
colorbar;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_matrix.png")
end

% Average HbT across all events

    t_val_AVG = mean(t_val')';
    t_val_pc_pos = size(find(t_val_AVG>0),1)/nnz(mask_j_thresh_allnodes_allch_gmsurface)*100;
    t_val_pc_neg = size(find(t_val_AVG<0),1)/nnz(mask_j_thresh_allnodes_allch_gmsurface)*100;

    figure()
    plotmesh_iso2([gmSurfaceMesh.node t_val_AVG],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    %clim([-max_t_val max_t_val]);
    clim([-max(abs(t_val_AVG)) max(abs(t_val_AVG))]);
    cb = colorbar('horiz');
    ylabel(cb,"T value")
    title("T val. AVG HbT, PD"+num2str(subjectN)+" "+eventType_all+" "+num2str(eventN_t_test(i))+" TW "+time_window+" "+num2str(t_val_pc_pos)+"% + "+num2str(t_val_pc_neg)+" % -")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_avg.png")
    end

%HbT_all_m_hypo_gmsurface_AVG = mean(HbT_all_m_hypo_gmsurface);

%create overlap between active R value nodes and Net T value nodes

%savefig = 1;
%pick three colours
%red - R val active only 
%blue - Net T active only
%green - Both R and Net V active
%active could be defined by >0.5 or <-0.5 in each metric
%%
%r 
rho;
r_thresh_active = 0.5;
r_thresh_active_pval = 0.05;
t_thresh_active = 0.25;

%rho_active_mask = find(abs(rho)>=r_thresh_active);
rho_active_mask = find(abs(pval)<=r_thresh_active_pval);

%t net
t_val_pos_neg_mask_net_norm;

t_val_active_mask = find(abs(t_val_pos_neg_mask_net_norm)>=t_thresh_active);

r_t_val_active_mask = zeros(size(HbT_all_gmsurface,1),1);

r_t_val_active_mask(rho_active_mask,1) = r_t_val_active_mask(rho_active_mask,1)+1;
r_t_val_active_mask(t_val_active_mask,1) = r_t_val_active_mask(t_val_active_mask,1)+2;

figure()
plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
view(0,90)
clim([0 3])
cb = colorbar('horiz');
ylabel(cb,"Active R(1) T(2) RT(3)")
title("Active R T val HbT, PD"+num2str(subjectN)+" "+eventType_all+" N events = "+num2str(size(eventN_t_test,2))+" TW "+time_window+" R thresh pval< "+num2str(r_thresh_active_pval)+" T thresh "+num2str(t_thresh_active)+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
     saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_active_R_T_val.png")
end

% figure()
% plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
% view(0,90)
% clim([0 3])
% cb = colorbar('horiz');
% ylabel(cb,"Active R(1) T(2) RT(3)")
%     title("Active R T val HbT, PD"+num2str(subjectN)+" "+eventType_all+" N events = "+num2str(size(eventN_t_test,2))+" TW "+time_window+" R thresh pval "+num2str(r_thresh_active_pval)+" T thresh "+num2str(t_thresh_active)+"")
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     if savefig == 1
%         saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_active_R_T_val.png")
%     end

figure()
subplot(1,3,1)
plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
view(0,90)
clim([0 3])
cb = colorbar('horiz');
ylabel(cb,"Active R(1) T(2) RT(3)")
title("Active R T val HbT, PD"+num2str(subjectN)+" "+eventType_all+" N events = "+num2str(size(eventN_t_test,2))+" TW "+time_window+" R thresh pval "+num2str(r_thresh_active_pval)+" T thresh "+num2str(t_thresh_active)+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(1,3,2)
plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
view(-90,0)
clim([0 3])
cb = colorbar('horiz');
ylabel(cb,"Active R(1) T(2) RT(3)")
title("Left Hemishpere")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(1,3,3)
plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
view(90,0)
clim([0 3])
cb = colorbar('horiz');
ylabel(cb,"Active R(1) T(2) RT(3)")
title("Right Hemishpere")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

if savefig == 1
     saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_active_R_T_val_CLRview.png")
end


end