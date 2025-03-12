function [PD_data] = rem_OD_outliers_10mins(PD_data,std_factor,t_pointsN)

% Remove outlyers from DOD signal
%consider just good ch idx

dod1_mean = zeros(size(PD_data.dod1,1),1);
dod1_std = zeros(size(PD_data.dod1,1),1);

for i=1:size(PD_data.dod1,1)
    dod1_mean(i,1) = mean(PD_data.dod1(i,PD_data.goodch_idx));
    dod1_std(i,1) = std(PD_data.dod1(i,PD_data.goodch_idx));
end

%outlier factor
%std_factor = 1.5;

figure()
plot(PD_data.t,PD_data.dod1(:,PD_data.goodch_idx(1:end/2)))
hold on
plot(PD_data.t,dod1_mean,'-','LineWidth',6,'Color',[1 0 0 0.5])
plot(PD_data.t,dod1_mean+(std_factor*dod1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
plot(PD_data.t,dod1_mean-(std_factor*dod1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
ylabel('Dod (M)')
xlabel('Time / s')
xline(900,'g--');
title("Dod 1 - with "+num2str(std_factor)+ " STD - goodch idx - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+" ")

%1 = timep is outlier, 0 = timep not outlier
outlier_more = PD_data.dod1(:,PD_data.goodch_idx)> dod1_mean+(std_factor*dod1_std); %above XSTD
outlier_less = PD_data.dod1(:,PD_data.goodch_idx)< dod1_mean-(std_factor*dod1_std); %below XSTD

%find N time pts that are outliers
for i=1:size(outlier_more,2)
    N_outlier_more(1,i) = nnz(outlier_more(:,i));
    N_outlier_less(1,i) = nnz(outlier_less(:,i));   
end
N_outlier_total = N_outlier_more + N_outlier_less;
%if X% (t_pointsN) of timepts are outliers, disregard ch
%outlier ch given 0idx , goodch given 1dx

ch_outlier = N_outlier_total< (size(PD_data.t,1)*t_pointsN);

%update goodch_idx
goodch_idx = [PD_data.goodch_idx(ch_outlier==1) ]; %;goodch_idx(ch_outlier==1)+64];
%goodch_idx = nonzeros(goodch_idx);
%update goodch idx, remch and SD meas.
    for i=1:size(goodch_idx,1)
        if (ismember(goodch_idx(i)+64,goodch_idx)) %| ismember(goodch_idx(i)-64,goodch_idx)
            %goodch_idx(i)
            goodch_idx(i) = goodch_idx(i);
            %remCh(i)=1;
        else
            goodch_idx(i) = 0;
        end
    end
    
    goodch_idx = nonzeros(goodch_idx);
    goodch_idx = [goodch_idx; goodch_idx+64];
    
    PD_data.SD.MeasListAct(goodch_idx) = 1;
    PD_data.goodch_idx = goodch_idx;

figure()
plot(PD_data.t,PD_data.dod1(:,goodch_idx))
hold on
plot(PD_data.t,dod1_mean,'-','LineWidth',6,'Color',[1 0 0 0.5])
plot(PD_data.t,dod1_mean+(std_factor*dod1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
plot(PD_data.t,dod1_mean-(std_factor*dod1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
ylabel('d HbO (M)')
xlabel('Time / s')
xline(900,'g--');
title("Outliers rm "+std_factor+" STD-fraction "+t_pointsN+" t pts outlier - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+"")


% perform outlier again - but for first 60s %%%%%%%%%%%%%%%%%%%%%%%%%5
clear var N_outlier_more N_outlier_less N_outlier_total ch_outlier
dod1_mean = zeros(size(PD_data.dod1,1),1);
dod1_std = zeros(size(PD_data.dod1,1),1);

for i=1:size(PD_data.dod1,1)
    dod1_mean(i,1) = mean(PD_data.dod1(i,PD_data.goodch_idx));
    dod1_std(i,1) = std(PD_data.dod1(i,PD_data.goodch_idx));
end

%1 = timep is outlier, 0 = timep not outlier
outlier_more = PD_data.dod1(1:((PD_data.fs*60)+1),PD_data.goodch_idx)> dod1_mean(1:((PD_data.fs*60)+1))+(std_factor*dod1_std(1:((PD_data.fs*60)+1))); %above XSTD
outlier_less = PD_data.dod1(1:((PD_data.fs*60)+1),PD_data.goodch_idx)< dod1_mean(1:((PD_data.fs*60)+1))-(std_factor*dod1_std(1:((PD_data.fs*60)+1))); %below XSTD

%find N time pts that are outliers
for i=1:size(outlier_more,2)
    N_outlier_more(1,i) = nnz(outlier_more(:,i));
    N_outlier_less(1,i) = nnz(outlier_less(:,i));   
end
N_outlier_total = N_outlier_more + N_outlier_less;
%if X% (t_pointsN) of timepts are outliers, disregard ch
%outlier ch given 0idx , goodch given 1dx
%t pts x 0.1 since we are looking at first 1min, not 10 min.

ch_outlier = N_outlier_total< (size(PD_data.t,1)*(t_pointsN*0.1));

%update goodch_idx
goodch_idx = [PD_data.goodch_idx(ch_outlier==1) ]; %;goodch_idx(ch_outlier==1)+64];
%goodch_idx = nonzeros(goodch_idx);
%update goodch idx, remch and SD meas.
    for i=1:size(goodch_idx,1)
        if (ismember(goodch_idx(i)+64,goodch_idx)) %| ismember(goodch_idx(i)-64,goodch_idx)
            %goodch_idx(i)
            goodch_idx(i) = goodch_idx(i);
            %remCh(i)=1;
        else
            goodch_idx(i) = 0;
        end
    end
    
    goodch_idx = nonzeros(goodch_idx);
    goodch_idx = [goodch_idx; goodch_idx+64];
    
    PD_data.SD.MeasListAct(goodch_idx) = 1;
    PD_data.goodch_idx = goodch_idx;


figure()
plot(PD_data.t,PD_data.dod1(:,goodch_idx))
hold on
plot(PD_data.t,dod1_mean,'-','LineWidth',6,'Color',[1 0 0 0.5])
plot(PD_data.t,dod1_mean+(std_factor*dod1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
plot(PD_data.t,dod1_mean-(std_factor*dod1_std),'-','LineWidth',6,'Color',[0 1 0 0.5])
ylabel('d HbO (M)')
xlabel('Time / s')
xline(900,'g--');
title("OutL rm from 1st 60s"+std_factor+" STD-fraction "+t_pointsN+" t pts outlier - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+"")

end