function calculate_R_val_2_metrics(BCG_metric,HbT_metric,metric_name,HbT_all_gmsurface,PD,mask_j_thresh_allnodes_allch_gmsurface,subjectN,eventType,time_window,gmSurfaceMesh,savefig,good_PD_events_subN,stacked)


% run this
rho = zeros(size(HbT_all_gmsurface,1),1);
pval = zeros(size(HbT_all_gmsurface,1),1);
Rsq1 = zeros(size(HbT_all_gmsurface,1),1);

        % %SB code from whatsapp 28 11 24
        % [rho,pVal(j,iCh)] = corr(mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF),blockAvgTime(mean_dc(:,1,chInt(iCh),iBHRF)==0,iB));
        % 
        % %subplot(length(chInt),3,k)
        % if pVal(j,iCh) <= .05
        %     scatter(mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF),blockAvgTime(mean_dc(:,1,chInt(iCh),iBHRF)==0,iB),'r')
        % else
        %     scatter(mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF),blockAvgTime(mean_dc(:,1,chInt(iCh),iBHRF)==0,iB),'filled')
        % end
        % 
        % [b,stats] = robustfit(mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF),blockAvgTime(mean_dc(:,1,chInt(iCh),iBHRF)==0,iB));
        % 
        % Rsq1 = corr(blockAvgTime(mean_dc(:,1,chInt(iCh),iBHRF)==0,iB),b(1)+b(2)*mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF))^2;
        % hold on;
        % plot(mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF),b(1)+b(2)*mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF),'r','LineWidth',2)
        % 


%% Get R val %one each for each node
for i=1:size(HbT_all_gmsurface,1)
    %Get R val
    if size(HbT_all_gmsurface,3) > 2 %b fit only works with > 2 events 

        %[b, stats] = robustfit(peak_glucose_all,HbT_all_max(1,:)'   );
        [b, stats] = robustfit(BCG_metric ,HbT_metric(i,:) );
    
        pval(i) = stats.p(2); %the pval is from the pval of the robustfit gradient b(2)

        %if grad is pos, keep R^2 pos, if grad is neg, make R^2 neg
        if b(2)<0
            Rsq1(i) = - (corr(  HbT_metric(i,:)',b(1)+b(2)*BCG_metric  )^2);
        else
            Rsq1(i) = corr(  HbT_metric(i,:)',b(1)+b(2)*BCG_metric  )^2;
        end

        %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   b(1) + b(2)*peak_glucose_all     );

        %used pre 27 11 24
        %[rho(i),pval(i)] = corr(HbT_metric(i,:)',   b(1) + b(2)*BCG_metric);
        
        %[rho(i),pval(i)] = corr(HbT_metric(i,:)', BCG_metric);
        %[rho(i),pval(i)] = corr(BCG_metric, HbT_metric(i,:)');


        


        % figure()
        % plot( BCG_metric,HbT_metric(i,:)','rx','MarkerSize',12)
        % hold on
        % plot(b(1) + b(2)*BCG_metric , HbT_metric(i,:)','bx','MarkerSize',12)
        % xlabel('BGC metric')
        % ylabel('HbT metric')
        % 
        % figure()
        % plot( BCG_metric,HbT_metric(i,:)','rx','MarkerSize',12)
        % hold on
        % plot(b(1) + b(2)*BCG_metric , HbT_metric(i,:)','b-','MarkerSize',12)
        % xlabel('BGC metric')
        % ylabel('HbT metric')
        % HbT_metric(i,:)',   b(1) + b(2)*BCG_metric



    else
        %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   peak_glucose_all     );
        [rho(i),pval(i)] = corr(HbT_metric(i,:)',   BCG_metric    );


        %rho(2570)
        %rho(2168)

    end

    %length of glucose
    %[b, stats] = robustfit(length_glucose_all_m_S_hypo,HbT_all_m_S_hypo_max(1,:)'   );
    %[rho(i),pval(i)] = corr(HbT_all_m_S_hypo_max(1,:)',   b(1) + b(2)*length_glucose_all_m_S_hypo     );
    %for j=1:size(HbT_all_gmsurface,3)
    %    [t_test(i,j),pval_t(i,j)] = ttest(HbT_all_gmsurface(i,:,j));
    %    [h_val_t(i,j),pval_t(i,j),ci_t,stats] = ttest(HbT_all_gmsurface(i,:,j));
    %    t_val(i,j) = stats.tstat;
    %end
end 

% %%
% 
%         Rsq1 = corr(  HbT_metric(i,:)',b(1)+b(2)*BCG_metric  )^2;
%         hold on;
%         plot(mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF),b(1)+b(2)*mean_dc(mean_dc(:,1,chInt(iCh),iBHRF)==0,1,iCh,iBHRF),'r','LineWidth',2)
% 
% 
% i = find(Rsq1 == max(Rsq1))
% min(Rsq1)
% 
%         figure()
%         plot( BCG_metric,HbT_metric(i,:)','rx','MarkerSize',12)
%         hold on
%         %plot(b(1) + b(2)*BCG_metric , HbT_metric(i,:)','bx','MarkerSize',12)
%         plot(BCG_metric,b(1) + b(2)*BCG_metric,'cx-','MarkerSize',12)
% 
%         xlabel('BGC metric')
%         ylabel('HbT metric')


%% FDR correct P vals and plot
adj_p = ones(1,size(mask_j_thresh_allnodes_allch_gmsurface,1));
%adj_p_t = ones(size(HbT_all_gmsurface,3),size(mask_j_thresh_allnodes_allch_gmsurface,1));

%adj_p_mask_t = ones(size(mask_95_allnodes_allch_gmsurface,1),size(HbT_all_m_hypo_gmsurface,3));
% FDR correction
[h, crit_p, adj_ci_cvrg, adj_p_mask]=fdr_bh(pval(logical(mask_j_thresh_allnodes_allch_gmsurface)));

%set nodes in top 95% mask to be corrected p vals.
adj_p(logical(mask_j_thresh_allnodes_allch_gmsurface)) = adj_p_mask;

%max(rho)
%set R P to be zero for nodes outside of mask.
%rho(masks.mask_nodes_bottom_5_cortexnodes_allch)=0;
%pval(masks.mask_nodes_bottom_5_cortexnodes_allch)=-1;

rho(~logical(mask_j_thresh_allnodes_allch_gmsurface))=0;
pval(~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;
adj_p(~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;
Rsq1(~logical(mask_j_thresh_allnodes_allch_gmsurface))=0;

%%%%%%%%%%%%%% Save fig 1.Rsqr,pval, fdr pval
figure()
subaxis(1,3,1,'sh',0.04) 
plotmesh_iso2([gmSurfaceMesh.node Rsq1],gmSurfaceMesh.face)
axis off
view(0,90); %TOP VIEW % used for infant week 30
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"\pm R^2")
ax = gca;
ax.FontSize = 20;
subaxis(1,3,2,'sh',0.04) 
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
axis off
view(0,90); %TOP VIEW % used for infant week 30
clim([-0.05 0.05]);
cb = colorbar('horiz');
ylabel(cb,"p")
ax = gca;
ax.FontSize = 20;
subaxis(1,3,3,'sh',0.04) 
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
axis off
view(0,90); %TOP VIEW % used for infant week 30
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"fdr p")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
ax = gca;
ax.FontSize = 20;
sgtitle("PD"+num2str(PD.subjectN)+" "+PD.eventType+" "+num2str(PD.eventN)+" TimeW "+time_window+" Rsqr val "+metric_name+"")
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_Rsqr_val"+metric_name+".png")
end


% figure()
% subplot(1,3,1)
% plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
% axis off
% view(0,90); %TOP VIEW % used for infant week 30
% clim([-1 1]);
% cb = colorbar('horiz');
% ylabel(cb,"\pm R^2")
% title("R val. "+metric_name+" , PD"+num2str(subjectN)+" "+eventType+"")
% subplot(1,3,2)
% plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.1 0.1]);
% cb = colorbar('horiz');
% ylabel(cb,"P value")
% title("P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% subplot(1,3,3)
% plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.05 0.05]);
% cb = colorbar('horiz');
% ylabel(cb,"FDR cor. P value")
% title("FDR adj P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% figure()
% plot(BCG_metric, HbT_metric(find(Rsq1==max(Rsq1)),:), 'rx')
% hold on
% plot(BCG_metric, HbT_metric(find(Rsq1==min(Rsq1)),:), 'bx')
% 
% find(Rsq1==max(Rsq1))
% 
% find(Rsq1==min(Rsq1))


if size(HbT_all_gmsurface,3) > 1 %max rho > 1 events 
    max_rho = max([max(rho) abs(min(rho))]);
else
    max_rho = 1;
end

max_pval = max([max(pval) abs(min(pval))]);
max_adj_p = max([max(adj_p) abs(min(adj_p))]);

% run this plot mesh (iso2mesh)
%rho=rho';
%pval=pval';

% figure()
% subplot(1,3,1)
% plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
% view(0,90); %TOP VIEW % used for infant week 30
% clim([-max_rho max_rho]);
% cb = colorbar('horiz');
% ylabel(cb,"R value")
% title("R val. "+metric_name+" , PD"+num2str(subjectN)+" "+eventType+"")
% subplot(1,3,2)
% plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.1 0.1]);
% cb = colorbar('horiz');
% ylabel(cb,"P value")
% title("P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% subplot(1,3,3)
% plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
% view(0,90);
% clim([-max_adj_p max_adj_p]);
% cb = colorbar('horiz');
% ylabel(cb,"FDR cor. P value")
% title("FDR adj P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% figure()
% subplot(1,3,1)
% plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
% view(0,90); %TOP VIEW % used for infant week 30
% %clim([-max_rho max_rho]);
% clim([-1 1]);
% cb = colorbar('horiz');
% ylabel(cb,"R value")
% title("R val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+"")
% subplot(1,3,2)
% plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.1 0.1]);
% cb = colorbar('horiz');
% ylabel(cb,"P value")
% title("P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% subplot(1,3,3)
% plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
% view(0,90);
% clim([-max_adj_p max_adj_p]);
% cb = colorbar('horiz');
% ylabel(cb,"FDR cor. P value")
% title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% if savefig == 1
%     saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_val"+metric_name+".png")
% end

% figure()
% subplot(1,3,1)
% plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
% view(0,90); %TOP VIEW % used for infant week 30
% %clim([-max_rho max_rho]);
% clim([-1 1]);
% cb = colorbar('horiz');
% ylabel(cb,"R value")
% title("R val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+"")
% ax = gca;
% ax.FontSize = 10;
% subplot(1,3,2)
% plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.05 0.05]);
% cb = colorbar('horiz');
% ylabel(cb,"P value")
% title("P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% ax = gca;
% ax.FontSize = 10;
% subplot(1,3,3)
% plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
% view(0,90);
% clim([-max_adj_p max_adj_p]);
% cb = colorbar('horiz');
% ylabel(cb,"FDR cor. P value")
% title("FDR adj P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% if savefig == 1
%     saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_val_p0_05"+metric_name+".png")
% end
%ax = gca;
%ax.FontSize = 10;

% figure()
% subplot(1,3,1)
% plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
% view(0,90); %TOP VIEW % used for infant week 30
% %clim([-max_rho max_rho]);
% clim([-1 1]);
% cb = colorbar('horiz');
% ylabel(cb,"R value")
% title("R val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+"")
% ax = gca;
% ax.FontSize = 10;
% subplot(1,3,2)
% plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.05 0.05]);
% cb = colorbar('horiz');
% ylabel(cb,"P value")
% title("P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% ax = gca;
% ax.FontSize = 10;
% subplot(1,3,3)
% plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
% view(0,90);
% clim([-0.1 0.1]);
% cb = colorbar('horiz');
% ylabel(cb,"FDR cor. P value")
% title("FDR adj P val. "+metric_name+", PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% if savefig == 1
%     saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_val_p0_05_adjp_0_1"+metric_name+".png")
% end
% ax = gca;
% ax.FontSize = 10;

%% plot graph
RGB_col = [241/255 71/255 37/255;
241/255 245/255 33/255;
17/255 251/255 56/255;
255/255 144/255 37/255;
236/255 137/255 132/255;
202/255 125/255 187/255;
136/255 8/255 165/255;
100/255 73/255 200/255;  %pd41
16/255 81/255 165/255;
78/255 34/255 205/255;   %45
50/255 220/255 60/255;   %48
13/255 8/255 241/255;
50/255 200/255 100/255; %pd55
54/255 208/255 182/255;
67/255 34/255 10/255;
200/255 168/255 175/255;
188/255 143/255 16/255;
251/255 1/255 165/255;];

col_scat_idx = zeros(size(good_PD_events_subN,2),1);
col_scat = zeros(size(good_PD_events_subN,2),1);
for i=1:size(good_PD_events_subN,2)
    col_scat_idx(i,1) = find(good_PD_events_subN(i) == unique(good_PD_events_subN));
end
%% get rval of pval <0.05

%average across p<0.05 nodes to get scatter plot
p_val_0_05_R = intersect(find(pval<=0.05),find(pval>=0)); %pval
p_val_fdr_0_05_R = intersect(find(adj_p'<=0.05),find(adj_p'>=0)); %FDR CORR PVAL!!

%find nodes with neg and pos R^2
Rsq_pos = find(Rsq1>0); %pos R^2
Rsq_neg = find(Rsq1<0); %neg R^2


% max_HbT_p_val_0_05_R_all = zeros(size(p_val_0_05_R,1),size(good_PD_events_subN,2));
% max_HbT_p_val_fdr_0_05_R_all = zeros(size(p_val_fdr_0_05_R,1),size(good_PD_events_subN,2));
% max_HbT_p_val_0_05_R_all = mean(HbT_metric(p_val_0_05_R,:));
% max_HbT_p_val_fdr_0_05_R_all = mean(HbT_metric(p_val_fdr_0_05_R,:));
% max_HbT_p_val_0_05_R_all = zeros(size(p_val_0_05_R,1),size(good_PD_events_subN,2));
% max_HbT_p_val_fdr_0_05_R_all = zeros(size(p_val_fdr_0_05_R,1),size(good_PD_events_subN,2));

%mean HbT from Pos R^2
%max_HbT_p_val_0_05_R_all_pos = mean(HbT_metric(intersect(p_val_0_05_R,Rsq_pos),:));
%max_HbT_p_val_fdr_0_05_R_all_pos = mean(HbT_metric(intersect(p_val_fdr_0_05_R,Rsq_pos),:));

%max Rsqr val of Pos R^2
max_HbT_p_val_0_05_R_all_pos = HbT_metric(intersect(find(Rsq1==max(Rsq1)),p_val_0_05_R),:);
max_HbT_p_val_fdr_0_05_R_all_pos = HbT_metric(intersect(find(Rsq1==max(Rsq1)),p_val_fdr_0_05_R),:);

node_max_Rsqr_p_val_0_05_pos = intersect(find(Rsq1==max(Rsq1)),p_val_0_05_R);
node_max_Rsqr_p_val_fdr_0_05_pos = intersect(find(Rsq1==max(Rsq1)),p_val_fdr_0_05_R);

max_Rsq1_pos_p_val_0_05 = Rsq1(intersect(find(Rsq1==max(Rsq1)),p_val_0_05_R));
max_Rsq1_pos_p_val_fdr_0_05 = Rsq1(intersect(find(Rsq1==max(Rsq1)),p_val_fdr_0_05_R));

pval_max_Rsq1_pos_p_val_0_05 = pval(intersect(find(Rsq1==max(Rsq1)),p_val_0_05_R));
pval_max_Rsq1_pos_p_val_fdr_0_05 = adj_p(intersect(find(Rsq1==max(Rsq1)),p_val_fdr_0_05_R));

%if size p_val_fdr_0_05_R<0, i.e there are fdr pvals<0.05
if size(p_val_fdr_0_05_R,1)>0
    %if the max R^2 val has fdr p>0.05, find highest R^2 val. with p<0.05
    if size(max_HbT_p_val_fdr_0_05_R_all_pos,1)==0
        max_HbT_p_val_fdr_0_05_R_all_pos = HbT_metric(find(Rsq1==max(Rsq1(p_val_fdr_0_05_R))),:);
        node_max_Rsqr_p_val_fdr_0_05_pos = find(Rsq1==max(Rsq1(p_val_fdr_0_05_R)));
        max_Rsq1_pos_p_val_fdr_0_05 = Rsq1(find(Rsq1==max(Rsq1(p_val_fdr_0_05_R))));
        pval_max_Rsq1_pos_p_val_fdr_0_05 = adj_p(find(Rsq1==max(Rsq1(p_val_fdr_0_05_R))));
    end
end

%mean HbT from neg R^2
%max_HbT_p_val_0_05_R_all_neg = mean(HbT_metric(intersect(p_val_0_05_R,Rsq_neg),:));
%max_HbT_p_val_fdr_0_05_R_all_neg = mean(HbT_metric(intersect(p_val_fdr_0_05_R,Rsq_neg),:));

%max Rsqr val of neg R^2
max_HbT_p_val_0_05_R_all_neg = HbT_metric(intersect(find(Rsq1==min(Rsq1)),p_val_0_05_R),:);
max_HbT_p_val_fdr_0_05_R_all_neg = HbT_metric(intersect(find(Rsq1==min(Rsq1)),p_val_fdr_0_05_R),:);

max_Rsq1_neg_p_val_0_05 = Rsq1(intersect(find(Rsq1==min(Rsq1)),p_val_0_05_R));
max_Rsq1_neg_p_val_fdr_0_05 = Rsq1(intersect(find(Rsq1==min(Rsq1)),p_val_fdr_0_05_R));

pval_max_Rsq1_neg_p_val_0_05 = pval(intersect(find(Rsq1==min(Rsq1)),p_val_0_05_R));
pval_max_Rsq1_neg_p_val_fdr_0_05 = adj_p(intersect(find(Rsq1==min(Rsq1)),p_val_fdr_0_05_R));

node_max_Rsqr_p_val_0_05_neg = intersect(find(Rsq1==min(Rsq1)),p_val_0_05_R);
node_max_Rsqr_p_val_fdr_0_05_neg = intersect(find(Rsq1==min(Rsq1)),p_val_fdr_0_05_R);

%if size p_val_fdr_0_05_R<0, i.e there are fdr pvals<0.05
if size(p_val_fdr_0_05_R,1)>0
%if the max R^2 val has fdr p>0.05, find highest R^2 val. with p<0.05
    if size(max_HbT_p_val_fdr_0_05_R_all_neg,1)==0
        max_HbT_p_val_fdr_0_05_R_all_neg = HbT_metric(find(Rsq1==min(Rsq1(p_val_fdr_0_05_R))),:);
        node_max_Rsqr_p_val_fdr_0_05_neg = find(Rsq1==min(Rsq1(p_val_fdr_0_05_R)));
        max_Rsq1_neg_p_val_fdr_0_05 = Rsq1(find(Rsq1==min(Rsq1(p_val_fdr_0_05_R))));
        pval_max_Rsq1_neg_p_val_fdr_0_05 = adj_p(find(Rsq1==min(Rsq1(p_val_fdr_0_05_R))));
    end
end
% figure()
% plot(BCG_metric,HbT_metric(find(Rsq1==max(Rsq1)),:),'rx-')
% hold on
% plot(BCG_metric,HbT_metric(find(Rsq1==min(Rsq1)),:),'bx-')
% find(Rsq1==max(Rsq1))
% find(Rsq1==min(Rsq1))


%find(Rsq_pos==max(Rsq_pos))

%bins = conncomp(G)

%figure()
%plot(BCG_metric,max_HbT_p_val_0_05_R_all_pos,'rx')
%hold on
%plot(BCG_metric,max_HbT_p_val_0_05_R_all_neg,'bx')

% figure()
% plot(BCG_metric,HbT_metric(intersect(p_val_0_05_R,Rsq_pos),:),'rx--')
% hold on
% figure()
% plot(BCG_metric,HbT_metric(intersect(p_val_0_05_R,Rsq_neg),:),'bx')
% 
% HbT_metric(intersect(p_val_0_05_R,Rsq_pos),:)
% HbT_metric(intersect(p_val_0_05_R,Rsq_neg),:)

N_nodes_Rsq_pos_pval_0_05 = size(intersect(p_val_0_05_R,Rsq_pos),1);
N_nodes_Rsq_neg_pval_0_05 = size(intersect(p_val_0_05_R,Rsq_neg),1);

N_nodes_Rsq_pos_pval_fdr_0_05 = size(intersect(p_val_fdr_0_05_R,Rsq_pos),1);
N_nodes_Rsq_neg_pval_fdr_0_05 = size(intersect(p_val_fdr_0_05_R,Rsq_neg),1);

% find(rho==max(rho))
% HbT_metric(find(rho==max(rho)),:)
% figure()
% plot(BCG_metric,HbT_metric(find(rho==max(rho)),:),'rx')
% title("rho ="+num2str(max(rho))+"")

% find(rho==min(rho))
% HbT_metric(find(rho==min(rho)),:)
% figure()
% plot(BCG_metric,HbT_metric(find(rho==min(rho)),:),'rx')
% title("rho ="+num2str(min(rho))+"")


%max_HbT_p_val_0_05_R;
%peak_glucose_all
% for i=1:size(p_val_0_05_R,1)
%     %Abs HbT
%     %max_HbT_p_val_0_05_R(i,:) = max(abs(HbT_all_gmsurface(p_val_0_05_R(i),121:end,:)));
% 
%     [max_HbT_p_val_0_05_R max_HbT_p_val_0_05_R_idx] = (max(abs(HbT_all_gmsurface(p_val_0_05_R(i),121:end,:)))); %get max HbT from T=2:10mins
%     max_HbT_p_val_0_05_R_idx = reshape(max_HbT_p_val_0_05_R_idx(1,1,:), 1, []);
% 
%     max_HbT_p_val_0_05_R = HbT_all_gmsurface(p_val_0_05_R(i),121+max_HbT_p_val_0_05_R_idx-1,:); %get max HbT from T=2:10mins
%     max_HbT_p_val_0_05_R = diag(reshape(max_HbT_p_val_0_05_R,[],size(max_HbT_p_val_0_05_R,3)))';%rehsape into sqaure, and get diag
%     max_HbT_p_val_0_05_R_all(i,:) = max_HbT_p_val_0_05_R;
% 
%     %HbT_metric
% 
% end

%c_max_HbT_p_val_0_05_R = polyfit(BCG_metric,mean(max_HbT_p_val_0_05_R_all(:,:)),1);


%HbT_metric(find(Rsq1==max(Rsq1)),:)

%%
%pos and neg line of best fit from p<0.05 of R^2
%%%%% for now 03 12 24 - use single node of highesr + or - rsqr value

%USE ROBUST FIT TO GET LINE OF BEST FIT

if size(max_HbT_p_val_0_05_R_all_pos,1) == size(BCG_metric,2)
    %c_max_HbT_p_val_0_05_R = polyfit(BCG_metric,max_HbT_p_val_0_05_R_all(:,:),1);

    %c_max_HbT_p_val_0_05_R_pos = polyfit(BCG_metric,max_HbT_p_val_0_05_R_all_pos(:,:),1);
    c_max_HbT_p_val_0_05_R_pos  = robustfit(BCG_metric ,max_HbT_p_val_0_05_R_all_pos );
    %c_max_HbT_p_val_0_05_R_pos = polyfit(BCG_metric,HbT_metric(find(Rsq1==max(Rsq1)),:),1);


    %%%c_max_HbT_p_val_fdr_0_05_R = polyfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all(:,:),1);

    %fit_y_c_p_glucose_max_r_HbT = c_max_HbT_p_val_0_05_R(1)*BCG_metric + c_max_HbT_p_val_0_05_R(2);
    fit_y_c_p_glucose_max_r_HbT_pos = c_max_HbT_p_val_0_05_R_pos(2)*BCG_metric + c_max_HbT_p_val_0_05_R_pos(1);
    %%%fit_y_c_p_fdr_glucose_max_r_HbT = c_max_HbT_p_val_fdr_0_05_R(1)*BCG_metric + c_max_HbT_p_val_fdr_0_05_R(2);
end
if size(max_HbT_p_val_0_05_R_all_neg,1) == size(BCG_metric,2)
    %c_max_HbT_p_val_0_05_R = polyfit(BCG_metric,max_HbT_p_val_0_05_R_all(:,:),1);

    %c_max_HbT_p_val_0_05_R_neg = polyfit(BCG_metric,max_HbT_p_val_0_05_R_all_neg(:,:),1);
   c_max_HbT_p_val_0_05_R_neg  = robustfit(BCG_metric ,max_HbT_p_val_0_05_R_all_neg );

    %%%c_max_HbT_p_val_fdr_0_05_R = polyfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all(:,:),1);

    %fit_y_c_p_glucose_max_r_HbT = c_max_HbT_p_val_0_05_R(1)*BCG_metric + c_max_HbT_p_val_0_05_R(2);
    fit_y_c_p_glucose_max_r_HbT_neg = c_max_HbT_p_val_0_05_R_neg(2)*BCG_metric + c_max_HbT_p_val_0_05_R_neg(1);
    %%%fit_y_c_p_fdr_glucose_max_r_HbT = c_max_HbT_p_val_fdr_0_05_R(1)*BCG_metric + c_max_HbT_p_val_fdr_0_05_R(2);
end

%pos and neg line of best fit from FDR p<0.05 of R^2
if size(max_HbT_p_val_fdr_0_05_R_all_pos,1) == size(BCG_metric,2)
    %c_max_HbT_p_val_0_05_R = polyfit(BCG_metric,max_HbT_p_val_0_05_R_all(:,:),1);
    %c_max_HbT_p_val_fdr_0_05_R_pos = polyfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all_pos(:,:),1);
        c_max_HbT_p_val_fdr_0_05_R_pos  = robustfit(BCG_metric ,max_HbT_p_val_fdr_0_05_R_all_pos );

    %%%c_max_HbT_p_val_fdr_0_05_R = polyfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all(:,:),1);

    %fit_y_c_p_glucose_max_r_HbT = c_max_HbT_p_val_0_05_R(1)*BCG_metric + c_max_HbT_p_val_0_05_R(2);
    fit_y_c_p_fdr_glucose_max_r_HbT_pos = c_max_HbT_p_val_fdr_0_05_R_pos(2)*BCG_metric + c_max_HbT_p_val_fdr_0_05_R_pos(1);
    %%%fit_y_c_p_fdr_glucose_max_r_HbT = c_max_HbT_p_val_fdr_0_05_R(1)*BCG_metric + c_max_HbT_p_val_fdr_0_05_R(2);
end
if size(max_HbT_p_val_fdr_0_05_R_all_neg,1) == size(BCG_metric,2)
    %c_max_HbT_p_val_0_05_R = polyfit(BCG_metric,max_HbT_p_val_0_05_R_all(:,:),1);
    %c_max_HbT_p_val_fdr_0_05_R_neg = polyfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all_neg(:,:),1);
    c_max_HbT_p_val_fdr_0_05_R_neg  = robustfit(BCG_metric ,max_HbT_p_val_fdr_0_05_R_all_neg );

    
    %%%c_max_HbT_p_val_fdr_0_05_R = polyfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all(:,:),1);

    %fit_y_c_p_glucose_max_r_HbT = c_max_HbT_p_val_0_05_R(1)*BCG_metric + c_max_HbT_p_val_0_05_R(2);
    fit_y_c_p_fdr_glucose_max_r_HbT_neg = c_max_HbT_p_val_fdr_0_05_R_neg(2)*BCG_metric + c_max_HbT_p_val_fdr_0_05_R_neg(1);
    %%%fit_y_c_p_fdr_glucose_max_r_HbT = c_max_HbT_p_val_fdr_0_05_R(1)*BCG_metric + c_max_HbT_p_val_fdr_0_05_R(2);
end



% if size(max_HbT_p_val_fdr_0_05_R_all,2) == size(BCG_metric,1)
%     c_max_HbT_p_val_fdr_0_05_R = polyfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all(:,:),1);
%     fit_y_c_p_fdr_glucose_max_r_HbT = c_max_HbT_p_val_fdr_0_05_R(1)*BCG_metric + c_max_HbT_p_val_fdr_0_05_R(2);
% end

%%%[b, stats] = robustfit(BCG_metric,mean(max_HbT_p_val_0_05_R_all(:,:))   );

%[b, stats] = robustfit(BCG_metric,max_HbT_p_val_0_05_R_all(:,:)   );


%if size(p_val_0_05_R,1)>2
%    [b, stats] = robustfit(BCG_metric,max_HbT_p_val_0_05_R_all(:,:)   );
%end
%%%%%%%get b val for robust fit of pos r^2 for p and pfdr p

% if size(max_HbT_p_val_0_05_R_all_pos,2)>2
%     [b_pos, stats_pos] = robustfit(BCG_metric,max_HbT_p_val_0_05_R_all_pos(:,:)   );
%     [rho_HbT_pval_0_05_R_all_pos,pvalHbT_pval_0_05_R_all_pos] = corr(max_HbT_p_val_0_05_R_all_pos(:,:)',   b_pos(1) + b_pos(2)*BCG_metric     );
% 
% end
% if size(max_HbT_p_val_fdr_0_05_R_all_pos,2)>2
%     [b_fdr_pos, stats_fdr_pos] = robustfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all_pos(:,:)   );
%     [rho_HbT_pval_fdr_0_05_R_all_pos,pvalHbT_pval_fdr_0_05_R_all_pos] = corr(max_HbT_p_val_fdr_0_05_R_all_pos(:,:)',   b_fdr_pos(1) + b_fdr_pos(2)*BCG_metric     );
% end
% %get b val for robust fit of neg r^2 for p and pfdr p
% if size(max_HbT_p_val_0_05_R_all_neg,2)>2
%     [b_neg, stats_neg] = robustfit(BCG_metric,max_HbT_p_val_0_05_R_all_neg(:,:)   );
%     [rho_HbT_pval_0_05_R_all_neg,pvalHbT_pval_0_05_R_all_neg] = corr(max_HbT_p_val_0_05_R_all_neg(:,:)',   b_neg(1) + b_neg(2)*BCG_metric     );
% end
% if size(max_HbT_p_val_fdr_0_05_R_all_neg,2)>2
%     [b_fdr_neg, stats_fdr_neg] = robustfit(BCG_metric,max_HbT_p_val_fdr_0_05_R_all_neg(:,:)   );
%     [rho_HbT_pval_fdr_0_05_R_all_neg,pvalHbT_pval_fdr_0_05_R_all_neg] = corr(max_HbT_p_val_fdr_0_05_R_all_neg(:,:)',   b_fdr_neg(1) + b_fdr_neg(2)*BCG_metric     );
% 
% end

%[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   b(1) + b(2)*peak_glucose_all     );
%[rho_HbT_pval_0_05_R_all,pvalHbT_pval_0_05_R_all] = corr(mean(max_HbT_p_val_0_05_R_all(:,:))',   b(1) + b(2)*BCG_metric     );

%[rho_HbT_pval_0_05_R_all,pvalHbT_pval_0_05_R_all] = corr(max_HbT_p_val_0_05_R_all(:,:)',   b(1) + b(2)*BCG_metric     );



% % 
% % if size(p_val_0_05_R,1)>2
% %     [rho_HbT_pval_0_05_R_all,pvalHbT_pval_0_05_R_all] = corr(max_HbT_p_val_0_05_R_all(:,:)',   b(1) + b(2)*BCG_metric     );
% % end
% % 
% % 
% % if size(p_val_fdr_0_05_R,1)>2
% %     [rho_HbT_pval_fdr_0_05_R_all,pvalHbT_pval_fdr_0_05_R_all] = corr(max_HbT_p_val_fdr_0_05_R_all(:,:)',   b(1) + b(2)*BCG_metric     );
% % end
%%
if stacked == 1
    switch time_window
        case "A"
            switch eventType
                case "m_hypo"
            col_scat_idx(1:3,1) = 1;
            col_scat_idx(4,1) = 2;
            col_scat_idx(5,1) = 3;
            col_scat_idx(6:8,1) = 4;
            col_scat_idx(9:19,1) = 6;
            col_scat_idx(20,1) = 7;
            col_scat_idx(21,1) = 8;
            col_scat_idx(22,1) = 9;
            col_scat_idx(23,1) = 10; %PD45
            col_scat_idx(24,1) = 11; %PD48
            col_scat_idx(23+2:36+2,1) = 10+2;
            col_scat_idx(37+2:38+2,1) = 11+2;
            col_scat_idx(39+2,1) = 12+2;
            col_scat_idx(40+2,1) = 14+2;
            col_scat_idx(41+2:44+2,1) = 15+2;
            %time w A all fnirs 24 poster
            % col_scat_idx(1:3,1) = 2;
            % col_scat_idx(4:8,1) = 4;
            % col_scat_idx(9:11,1) = 5;
            % col_scat_idx(12:22,1) = 6;
            % col_scat_idx(23:31,1) = 7;

             case "S_hypo"
            col_scat_idx(1,1) = 3;
            col_scat_idx(2,1) = 4;
            col_scat_idx(3:5,1) = 5;
            col_scat_idx(6:9,1) = 10+2;%+2 as added PD45 and 48
            %col_scat_idx(6:8,1) = 10;
            %'recomment col scat 6:9 when inc pd49'

            case "S_m_hypo"
            col_scat_idx(1:3,1) = 1;
            col_scat_idx(4,1) = 2;
            col_scat_idx(5,1) = 3;
            col_scat_idx(6:8,1) = 4;
            col_scat_idx(9:19,1) = 6;
            col_scat_idx(20,1) = 7;
            col_scat_idx(21,1) = 8;
            col_scat_idx(22,1) = 9;
            col_scat_idx(23,1) = 10;
            col_scat_idx(24,1) = 11;
            col_scat_idx(23+2:36+2,1) = 10+2;%+2 as added PD45 and 48;
            col_scat_idx(37+2:38+2,1) = 11+2;%+2 as added PD45 and 48;
            col_scat_idx(39+2,1) = 12+2;%+2 as added PD45 and 48;
            col_scat_idx(40+2,1) = 14+2;%+2 as added PD45 and 48;
            col_scat_idx(41+2:44+2,1) = 15+2;%+2 as added PD45 and 48;
           
            col_scat_idx(44+1+2,1) = 3;
            col_scat_idx(44+2+2,1) = 4;
            col_scat_idx(44+3+2:44+5+2,1) = 5;
            col_scat_idx(44+6+2:44+9+2,1) = 10+2;
           end

        case "D"
            switch eventType
                case "m_hypo"
            col_scat_idx(1:3,1) = 1;
            col_scat_idx(4,1) = 2;
            col_scat_idx(5:6,1) = 3;
            col_scat_idx(7,1) = 4;
            col_scat_idx(8:20,1) = 5;
            col_scat_idx(21,1) = 6;
            col_scat_idx(22,1) = 7;
            col_scat_idx(23,1) = 8;
            col_scat_idx(24,1) = 9;
            col_scat_idx(25:26,1) = 10;%pd45(2e)
            col_scat_idx(25+2:39+2,1) = 10+2;
            col_scat_idx(40+2,1) = 11+2;
            col_scat_idx(41+2:42+2,1) = 12+2;
            col_scat_idx(43+2,1) = 13+2;
            col_scat_idx(44+2:45+2,1) = 14+2;
            col_scat_idx(46+2:48+2,1) = 15+2;

            case "S_hypo"
            col_scat_idx(1,1) = 3;
            col_scat_idx(2,1) = 4;
            col_scat_idx(3:5,1) = 5;
            col_scat_idx(6,1) = 10+2;
            %'recomment col scat idx 6,1 when inc pd49'
            

            case "S_m_hypo"
            col_scat_idx(1:3,1) = 1;
            col_scat_idx(4,1) = 2;
            col_scat_idx(5:6,1) = 3;
            col_scat_idx(7,1) = 4;
            col_scat_idx(8:20,1) = 5;
            col_scat_idx(21,1) = 6;
            col_scat_idx(22,1) = 7;
            col_scat_idx(23,1) = 8;
            col_scat_idx(24,1) = 9;
            col_scat_idx(25:26,1) = 10;%pd45(2e)
            col_scat_idx(25+2:39+2,1) = 10+2;
            col_scat_idx(40+2,1) = 11+2;
            col_scat_idx(41+2:42+2,1) = 12+2;
            col_scat_idx(43+2,1) = 13+2;
            col_scat_idx(44+2:45+2,1) = 14+2;
            col_scat_idx(46+2:48+2,1) = 15+2;

            col_scat_idx(1+48+2,1) = 3;
             col_scat_idx(2+48+2,1) = 4;
             col_scat_idx(3+48+2:5+48+2,1) = 5;
             col_scat_idx(6+48+2,1) = 10+2;
            %time W D
            % col_scat_idx(1:3,1) = 2;
            % col_scat_idx(4:5,1) = 3;
            % col_scat_idx(6:13,1) = 4;
            % col_scat_idx(14,1) = 5;
            % col_scat_idx(15:27,1) = 6;
            % col_scat_idx(28:38,1) = 7;
            end
        case "F"
            switch eventType
                case "m_hypo"
            col_scat_idx(1:5,1) = 1;
            col_scat_idx(6,1) = 2;
            col_scat_idx(7:8,1) = 3;
            col_scat_idx(9:10,1) = 4;
            col_scat_idx(11:20,1) = 5;
            
            col_scat_idx(21,1) = 7;
            col_scat_idx(22,1) = 8;
            col_scat_idx(23,1) = 9;
            col_scat_idx(24:25,1) = 10;%pd45(e2)
            col_scat_idx(24+2:38+2,1) = 10+2;
            col_scat_idx(39+2:40+2,1) = 11+2;
            col_scat_idx(41+2:42+2,1) = 12+2;
            col_scat_idx(43+2,1) = 13+2;
            col_scat_idx(44+2:45+2,1) = 14+2;
            col_scat_idx(46+2:50+2,1) = 15+2;
            %time W F
            %7
            % col_scat_idx(1,1) = 1;
            % col_scat_idx(2:6,1) = 2;
            % col_scat_idx(7:8,1) = 3;
            % col_scat_idx(9:13,1) = 4;
            % col_scat_idx(14:15,1) = 5;
            % col_scat_idx(16:25,1) = 6;
            % col_scat_idx(26:37,1) = 7;
                case "S_hypo"
                col_scat_idx(1,1) = 3;
                col_scat_idx(2,1) = 4;
                col_scat_idx(3:5,1) = 5;
                col_scat_idx(6:7,1) = 10+2;
                
                case "S_m_hypo"
                col_scat_idx(1:5,1) = 1;
            col_scat_idx(6,1) = 2;
            col_scat_idx(7:8,1) = 3;
            col_scat_idx(9:10,1) = 4;
            col_scat_idx(11:20,1) = 5;
            col_scat_idx(21,1) = 7;
            col_scat_idx(22,1) = 8;
            col_scat_idx(23,1) = 9;
            col_scat_idx(24:25,1) = 10;%pd45(e2)
            col_scat_idx(24+2:38+2,1) = 10+2;
            col_scat_idx(39+2:40+2,1) = 11+2;
            col_scat_idx(41+2:42+2,1) = 12+2;
            col_scat_idx(43+2,1) = 13+2;
            col_scat_idx(44+2:45+2,1) = 14+2;
            col_scat_idx(46+2:50+2,1) = 15+2;

                col_scat_idx(50+1+2,1) = 3;
                col_scat_idx(50+2+2,1) = 4;
                col_scat_idx(50+3+2:50+5+2,1) = 5;
                col_scat_idx(50+6+2:50+7+2,1) = 10+2;

            

            end
    end
end

col_scat = RGB_col(col_scat_idx(:,1),:);

%%
%N_nodes_Rsq_neg_pval_fdr_0_05
if size(max_HbT_p_val_0_05_R_all_pos,1)>0
    figure()
    plot(BCG_metric,fit_y_c_p_glucose_max_r_HbT_pos,'k-','LineWidth',8)
    hold on
    scatter(BCG_metric,max_HbT_p_val_0_05_R_all_pos(:,:),220,col_scat,'filled')
    %xlim([40 72])
    %ylim([-13 13])
    xlabel(""+metric_name+" BGC")
    ylabel(""+metric_name+" HbT")
    title("Max+ R^2 P<0.05"+metric_name+". R^2="+num2str(max_Rsq1_pos_p_val_0_05)+" P="+num2str(pval_max_Rsq1_pos_p_val_0_05)+" S "+num2str(subjectN)+" "+eventType+" N Nds "+num2str(N_nodes_Rsq_pos_pval_0_05)+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R_pos(2))+" p0 "+num2str(c_max_HbT_p_val_0_05_R_pos(1))+" nd "+num2str(node_max_Rsqr_p_val_0_05_pos)+"")
    ax = gca;
    ax.FontSize = 26;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_p_0_05nodes_scatter_maxRsqr_colour_POS"+metric_name+".png")
    end
end

if size(max_HbT_p_val_fdr_0_05_R_all_pos,1)>0
    figure()
    plot(BCG_metric,fit_y_c_p_fdr_glucose_max_r_HbT_pos,'k-','LineWidth',8)
    hold on
    scatter(BCG_metric,max_HbT_p_val_fdr_0_05_R_all_pos(:,:),220,col_scat,'filled') %130
    %xlim([40 72])
    %ylim([-13 13])
    xlabel(""+metric_name+" BGC")
    ylabel(""+metric_name+" HbT")
    title("Max+ R^2 fdr P<0.05"+metric_name+". R^2="+num2str(max_Rsq1_pos_p_val_fdr_0_05)+" P="+num2str(pval_max_Rsq1_pos_p_val_fdr_0_05)+" S "+num2str(subjectN)+" "+eventType+" N Nds "+num2str(N_nodes_Rsq_pos_pval_fdr_0_05)+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_fdr_0_05_R_pos(2))+" p0 "+num2str(c_max_HbT_p_val_fdr_0_05_R_pos(1))+" nd "+num2str(node_max_Rsqr_p_val_fdr_0_05_pos)+"")
    ax = gca;
    ax.FontSize = 26; %20
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_p_fdr_0_05nodes_scatter_maxRsqr_colour_POS"+metric_name+".png")
    end
end

if size(max_HbT_p_val_0_05_R_all_neg,1)>0
    figure()
    plot(BCG_metric,fit_y_c_p_glucose_max_r_HbT_neg,'k-','LineWidth',8)
    hold on
    scatter(BCG_metric,max_HbT_p_val_0_05_R_all_neg(:,:),220,col_scat,'filled')
    %xlim([40 72])
    %ylim([-13 13])
    xlabel(""+metric_name+" BGC")
    ylabel(""+metric_name+" HbT")
    title("Max- R^2 P<0.05"+metric_name+". R^2="+num2str(max_Rsq1_neg_p_val_0_05)+" P="+num2str(pval_max_Rsq1_neg_p_val_0_05)+" S "+num2str(subjectN)+" "+eventType+" N Nds "+num2str(N_nodes_Rsq_neg_pval_0_05)+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R_neg(2))+" p0 "+num2str(c_max_HbT_p_val_0_05_R_neg(1))+" nd "+num2str(node_max_Rsqr_p_val_0_05_neg)+"")
    ax = gca;
    ax.FontSize = 26;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_p_0_05nodes_scatter_maxRsqr_colour_NEG"+metric_name+".png")
    end
   
end

if size(max_HbT_p_val_fdr_0_05_R_all_neg,1)>0
    figure()
    plot(BCG_metric,fit_y_c_p_fdr_glucose_max_r_HbT_neg,'k-','LineWidth',8)
    hold on
    scatter(BCG_metric,max_HbT_p_val_fdr_0_05_R_all_neg(:,:),220,col_scat,'filled')
    %xlim([40 72])
    %ylim([-13 13])
    xlabel(""+metric_name+" BGC")
    ylabel(""+metric_name+" HbT")
    title("Max- R^2 fdr P<0.05"+metric_name+". R^2="+num2str(max_Rsq1_neg_p_val_fdr_0_05)+" P="+num2str(pval_max_Rsq1_neg_p_val_fdr_0_05)+" S "+num2str(subjectN)+" "+eventType+" N Nds "+num2str(N_nodes_Rsq_neg_pval_fdr_0_05)+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_fdr_0_05_R_neg(2))+" p0 "+num2str(c_max_HbT_p_val_fdr_0_05_R_neg(1))+" nd "+num2str(node_max_Rsqr_p_val_fdr_0_05_neg)+"")
    ax = gca;
    ax.FontSize = 26;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_p_fdr_0_05nodes_scatter_maxRsqr_colour_NEG"+metric_name+".png")
    end
end




%% old pre 28 11 24
% if size(p_val_0_05_R,1)>2
%     figure()
%     plot(BCG_metric,fit_y_c_p_glucose_max_r_HbT,'k-','LineWidth',8)
%     hold on
%     scatter(BCG_metric,max_HbT_p_val_0_05_R_all(:,:),130,col_scat,'filled')
%     %xlim([40 72])
%     %ylim([-13 13])
%     xlabel(""+metric_name+" BGC")
%     ylabel(""+metric_name+" HbT")
%     title("Mean Com. R of P<0.05"+metric_name+". R="+num2str(rho_HbT_pval_0_05_R_all)+" P="+num2str(pvalHbT_pval_0_05_R_all)+" Sub. "+num2str(subjectN)+" "+eventType+" N CTX Nds "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
%     ax = gca;
%     ax.FontSize = 16;
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     if savefig == 1
%         saveas(gcf,"PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_p_0_05nodes_scatter_mean_colour"+metric_name+".png")
%     end
% end
% 
% if size(p_val_fdr_0_05_R,1)>2
%     figure()
%     plot(BCG_metric,fit_y_c_p_fdr_glucose_max_r_HbT,'k-','LineWidth',8)
%     hold on
%     scatter(BCG_metric,max_HbT_p_val_fdr_0_05_R_all(:,:),130,col_scat,'filled')
%     %xlim([40 72])
%     %ylim([-13 13])
%     xlabel(""+metric_name+" BGC")
%     ylabel(""+metric_name+" HbT")
%     title("Mean Com. R of FDR P<0.05"+metric_name+". R="+num2str(rho_HbT_pval_fdr_0_05_R_all)+" P="+num2str(pvalHbT_pval_fdr_0_05_R_all)+" Sub. "+num2str(subjectN)+" "+eventType+" N CTX Nds "+num2str(size(p_val_fdr_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_fdr_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_fdr_0_05_R(2))+" ")
%     ax = gca;
%     ax.FontSize = 16;
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     if savefig == 1
%         saveas(gcf,"PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(PD.eventN)+"_TimeW_"+time_window+"_R_FDR_p_0_05nodes_scatter_mean_colour"+metric_name+".png")
%     end
% 
% end
%%

end