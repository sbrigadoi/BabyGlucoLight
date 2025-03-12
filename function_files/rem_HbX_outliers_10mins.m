function [PD_data] = rem_HbX_outliers_10mins(PD_data,std_factor,t_pointsN)
% %% Remove outlyers from HbO  signal
% %consider just good ch idx
% 
HbO_1_mean = zeros(size(PD_data.HbO_1,1),1);
HbO_1_std = zeros(size(PD_data.HbO_1,1),1);

Hb_1_mean = zeros(size(PD_data.Hb_1,1),1);
Hb_1_std = zeros(size(PD_data.Hb_1,1),1);

HbT_1_mean = zeros(size(PD_data.HbT_1,1),1);
HbT_1_std = zeros(size(PD_data.HbT_1,1),1);

for i=1:size(PD_data.HbO_1,1)
    HbO_1_mean(i,1) = mean(PD_data.HbO_1(i,PD_data.goodch_idx(1:end/2)));
    HbO_1_std(i,1) = std(PD_data.HbO_1(i,PD_data.goodch_idx(1:end/2)));

    Hb_1_mean(i,1) = mean(PD_data.Hb_1(i,PD_data.goodch_idx(1:end/2)));
    Hb_1_std(i,1) = std(PD_data.Hb_1(i,PD_data.goodch_idx(1:end/2)));

    HbT_1_mean(i,1) = mean(PD_data.HbT_1(i,PD_data.goodch_idx(1:end/2)));
    HbT_1_std(i,1) = std(PD_data.HbT_1(i,PD_data.goodch_idx(1:end/2)));
end

figure()
plot(PD_data.t,PD_data.HbO_1(:,PD_data.goodch_idx(1:end/2)))
hold on
plot(PD_data.t,HbO_1_mean,'-','LineWidth',6,'Color',[1 0 0 0.5])
plot(PD_data.t,HbO_1_mean+(std_factor*HbO_1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
plot(PD_data.t,HbO_1_mean-(std_factor*HbO_1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
ylabel('d HbO (M)')
xlabel('Time / s')
xline(900,'g--');
title("HbO "+std_factor+" STD - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+"")

figure()
plot(PD_data.t,PD_data.Hb_1(:,PD_data.goodch_idx(1:end/2)))
hold on
plot(PD_data.t,Hb_1_mean,'-','LineWidth',6,'Color',[1 0 0 0.5])
plot(PD_data.t,Hb_1_mean+(std_factor*Hb_1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
plot(PD_data.t,Hb_1_mean-(std_factor*Hb_1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
ylabel('d Hb (M)')
xlabel('Time / s')
xline(900,'g--');
title("Hb "+std_factor+" STD - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+"")


%1 = timep is outlier, 0 = timep not outlier
outlier_more = PD_data.HbO_1(:,PD_data.goodch_idx(1:end/2))> HbO_1_mean+(std_factor*HbO_1_std); %above 2STD
outlier_less = PD_data.HbO_1(:,PD_data.goodch_idx(1:end/2))< HbO_1_mean-(std_factor*HbO_1_std); %below 2STD

%find N time pts that are outliers
N_outlier_more = zeros(1,size(outlier_more,2));
N_outlier_less = zeros(1,size(outlier_more,2));
for i=1:size(outlier_more,2)
    N_outlier_more(1,i) = nnz(outlier_more(:,i));
    N_outlier_less(1,i) = nnz(outlier_less(:,i));   
end
N_outlier_total = N_outlier_more + N_outlier_less;
%if X% of timepts are outliers, disregard ch
%outlier ch given 0idx , goodch given 1dx
%t_pointsN is the factor for how many time pts need to be outside the
%boundary to be an outlier ch
ch_outlier = N_outlier_total< (size(PD_data.t,1)*t_pointsN);

%update goodch_idx
goodch_idx = [PD_data.goodch_idx(ch_outlier==1) ;PD_data.goodch_idx(ch_outlier==1)+64];
PD_data.goodch_idx = goodch_idx;
%update goodch idx, remch and SD meas.
PD_data.SD.MeasListAct(goodch_idx) = 1;

figure()
plot(PD_data.t,PD_data.HbO_1(:,PD_data.goodch_idx(1:end/2)))
hold on
plot(PD_data.t,HbO_1_mean,'-','LineWidth',6,'Color',[1 0 0 0.5])
plot(PD_data.t,HbO_1_mean+(std_factor*HbO_1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
plot(PD_data.t,HbO_1_mean-(std_factor*HbO_1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
ylabel('d HbO (M)')
xlabel('Time / s')
xline(900,'g--');
title("HbO Outliers rm "+std_factor+" STD-fraction "+t_pointsN+" t pts outlier - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+"")


figure()
plot(PD_data.t,PD_data.Hb_1(:,PD_data.goodch_idx(1:end/2)))
hold on
plot(PD_data.t,Hb_1_mean,'-','LineWidth',6,'Color',[1 0 0 0.5])
plot(PD_data.t,Hb_1_mean+(std_factor*Hb_1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
plot(PD_data.t,Hb_1_mean-(std_factor*Hb_1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
ylabel('d Hb (M)')
xlabel('Time / s')
xline(900,'g--');
title(" Hb Outliers rm "+std_factor+" STD-fraction "+t_pointsN+" t pts outlier - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+"")


end