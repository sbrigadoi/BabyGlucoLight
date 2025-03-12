function [PD_data] = Iteration10minDOT(subjectN,eventType,time_window,i,J,mesh,gmSurfaceMesh,t_points,vol2gm,coverageThresh,event_good)

switch eventType
    case "m_hypo"
        eventN = event_good(i);
    case "S_hypo"
        eventN = event_good(i);
end


% switch subjectN
%     case 7 %add other subject Ns
%         switch eventType
%             case "m_hypo"
%                 PD_m_hypo_events = [8 9 10 11 12 13 14 15 16 17 18 19];
%                 PD_m_hypo_events = [8 9 10 11 13 14 17];
%                 eventN = PD_m_hypo_events(i);
%         end
%     case 8 %add other subject Ns
%         switch eventType
%             case "m_hypo"
%                 PD_m_hypo_events = [1 2 3 4 5 6 7 8 9 10 11];
%                 eventN = PD_m_hypo_events(i);
%         end
%     case 15
%         switch eventType
%             case "m_hypo"
%                 PD_m_hypo_events = [1 2 3 4 5 6 7 8];
%                 eventN = PD_m_hypo_events(i);
%             case "S_hypo"
%                 PD_S_hypo_events = [1 2];
%                 eventN = PD_S_hypo_events(i);
%         end
%     case 20 %add other subject Ns
%         switch eventType
%             case "m_hypo"
%                 switch time_window
%                     case "A"
%                       PD_m_hypo_events = [3 4 5 6 8 9 11];
% 
% 
%                     case "D"
%                         PD_m_hypo_events = [3 4 5 6 7 8 9 10 11 12];
%                         PD_m_hypo_events = [3 4 5 6 7 9 11 12];
%                 end
%             eventN = PD_m_hypo_events(i);
%         end
%     case 33 %add other subject Ns
%         switch eventType
%             case "m_hypo"
%                 switch time_window
%                     case "A"
%                     PD_m_hypo_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
%                     case "D"
%                     PD_m_hypo_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
%                     PD_m_hypo_events = [1 3 4 6 7 9 10 11 12 13 15 19]; %based on good idx 15 05 24
%                 end
%                 eventN = PD_m_hypo_events(i);
%             case "S_hypo"
%                 PD_S_hypo_events = [1 2 3]; %S hypo
%                 PD_S_hypo_events = [2 3]; %S hypo  %based on good idx 15 05 24
%                 eventN = PD_S_hypo_events(i);
%         end
% end


PD_data = load("H:\PadovaPostDoc\BabyGluCo\NIRS_data\PD"+num2str(subjectN)+"\formatted\PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(eventN)+".nirs",'-mat');
PD_data = getfield(PD_data,"PD"+num2str(subjectN)+"_"+eventType+"_X");
%PD_data = PD_data.PD33_m_hypo_X; %change this
PD = load("H:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(subjectN)+".mat");
%PD = load("E:\PadovaPostDoc\BabyGluCo\formatted_PD_GuyPerkins\PD"+num2str(subjectN)+".mat");

PD = getfield(PD,"PD"+num2str(subjectN)+"");
%PD = PD.PD33; %change this
switch eventType
        case "m_hypo"
            PD_time = 0:5:5*PD.events_length.m_hypo(eventN);
        case "S_hypo"
            PD_time = 0:5:5*PD.events_length.S_hypo(eventN);
end

PD_data.time_window = time_window;

PD.subjectN = subjectN;
PD.eventN = eventN;
PD.eventType = eventType;

PD_data.subjectN = subjectN;
PD_data.eventN = eventN;
PD_data.eventType = eventType;

% Crop data. Just want 2mins before, 8 mins after (so total 10)
[PD_data time_window_t] = get_timewindow_10mins(PD_data,PD,PD_time,time_window);


fs=10;
PD_data.fs = fs;

% % for poster SfNIRS 2024 04 09 24
% figure()
%     %plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(eventN,1):PD.events_start_end.m_hypo(eventN,2)),'b-o','LineWidth',8)
%     plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(eventN,1):PD.events_start_end.S_hypo(eventN,2)),'b-o','LineWidth',8)
% 
%     %axes('FontSize', 24);
%     %title("Time W"+time_window+" PD "+num2str(subjectN)+" "+eventType+" "+num2str(eventN)+". nirs")
%     xlabel("Time (Minutes)","FontSize",24)
%     ylabel("BGC (mg/dL)","FontSize",24)
%     xline(15,'--','LineWidth',6,'Color',[255/255, 144/255, 37/255]);xline(PD_time(end)-15,'-.','LineWidth',6,'Color',[255/255, 144/255, 37/255]);
%     xline(time_window_t(1),'m--','LineWidth',6);xline(time_window_t(2),'m--','LineWidth',6);
%     xline(0,'m--','LineWidth',6);xline(10,'m--','LineWidth',6);    
%     %xline(53,'m--','LineWidth',6);xline(63,'m--','LineWidth',6);
%     xline(PD_time(end)-17,'m--','LineWidth',6);xline(PD_time(end)-7,'m--','LineWidth',6);
%     yline(72,'k--','LineWidth',6);
%     yline(47,'k-.','LineWidth',6);
% ax = gca;
% ax.FontSize = 26;  % Font Size of 20,   15
%legend()



if eventType == "m_hypo"; %S_hypo
    figure()
    plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(eventN,1):PD.events_start_end.m_hypo(eventN,2)),'b-o')
    title("Time W"+time_window+" PD "+num2str(subjectN)+" "+eventType+" "+num2str(eventN)+". nirs")
    xlabel("Time / mins")
    ylabel("BGC / mg/DL")
    xline(15,'g--');xline(PD_time(end)-15,'r--');
    xline(time_window_t(1),'m--');xline(time_window_t(2),'m--');
    yline(72,'k--');
end
if eventType == "S_hypo"; %S_hypo
    figure()
    plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(eventN,1):PD.events_start_end.S_hypo(eventN,2)),'b-o')
    title("Time W"+time_window+" PD "+num2str(subjectN)+" "+eventType+" "+num2str(eventN)+". nirs")
    xlabel("Time / mins")
    ylabel("BGC / mg/DL")
    xline(15,'g--');xline(PD_time(end)-15,'r--');
    xline(time_window_t(1),'m--');xline(time_window_t(2),'m--');
    yline(72,'k--');
end

% Step 2. Signal Qaulity check
dRange = [5E-4 3];
SNRrange = 2;
[goodch_idx, PD_data] = find_good_ch_15mins(PD_data,dRange,SNRrange,1);

goodch_idx_step2 = PD_data.goodch_idx;

% Step 3. Get DoD
meanValue = mean(PD_data.d);
dodConv = -log(abs(PD_data.d)./meanValue);
dodConv_og = -log(abs(PD_data.d)./meanValue);
PD_data.dod = dodConv;

figure()
plot(PD_data.t,PD_data.dod(:,PD_data.goodch_idx))
xlabel('Time / S')
ylabel('\Delta OD')
title("Dod goodch idx - PD "+num2str(subjectN)+" event "+eventType+" N "+num2str(eventN)+" ")

if size(goodch_idx,1) < 12

    %PD_data.HbO = zeros(size(mesh.nodes,1),t_points);
    %PD_data.Hb = zeros(size(mesh.nodes,1),t_points);
    %PD_data.HbT = zeros(size(mesh.nodes,1),t_points);

    PD_data.HbO = zeros(size(vol2gm,1),t_points);
    PD_data.Hb = zeros(size(vol2gm,1),t_points);
    PD_data.HbT = zeros(size(vol2gm,1),t_points);

else


    %%%% Step 4. Motion corr. OUTPUT dod1
    methodN =1;
    PD_data = PD_data_motioncorr_15min(PD_data,methodN,1);

    goodch_idx_step4 = PD_data.goodch_idx;
    goodch_idx_step4_gvtd = PD_data.goodch_idx_gvdt;


    %poster SfNIRS 04 09 24
    % figure()
    % plot(PD_data.t,PD_data.dod(:,PD_data.goodch_idx(13)),'Color',[255/255,144/255,37/255],'LineWidth',3)
    % hold on
    % plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(13)),'Color',[87/255,104/255,210/255],'LineWidth',3)
    % xlabel('Time / S')
    % ylabel('\Delta OD')
    % ax = gca;
    % ax.FontSize = 16;


    figure()
    plot(PD_data.t,PD_data.dod(:,PD_data.goodch_idx(1)))
    hold on
    plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(1)))
    xlabel('Time / S')
    ylabel('\Delta OD')
    title("Dod1 (Mo corr.) goodch idx - PD "+num2str(subjectN)+" event "+eventType+" N "+num2str(eventN)+" ")

    figure()
    plot(PD_data.t,PD_data.dod1(:,PD_data.goodch_idx))
    xlabel('Time / S')
    ylabel('\Delta OD')
    title("Dod1 (Mo corr.) goodch idx - PD "+num2str(subjectN)+" event "+eventType+" N "+num2str(eventN)+" ")

    %%% Baseline subtraction
    PD_data.dod1(:,:) = PD_data.dod1(:,:) - mean(PD_data.dod1(1:2*fs*60,:));
    PD_data.dod1_noc(:,:) = PD_data.dod1_noc(:,:) - mean(PD_data.dod1_noc(1:2*fs*60,:));

    figure()
    plot(PD_data.t,PD_data.dod1(:,PD_data.goodch_idx))
    xlabel('Time / S')
    ylabel('\Delta OD')
    title("Baseline subtraction goodch idx - Subject"+num2str(subjectN)+" event "+eventType+" N "+num2str(eventN)+" ")
    % Remove outliers from dod1

    std_factor = 2;%1.5; %boundary of N STD's from mean
    t_pointsN = 0.25;%0.4; %fraction of timepoints outside of boundary to be an outlier channel
    [PD_data] = rem_OD_outliers_10mins(PD_data,std_factor,t_pointsN);
    goodch_idx_step4_OD_outliers = PD_data.goodch_idx;


    % Step 5. Spectoscropy
    close all;
    [PD_data] = spectroscopy_DPF_10mins(PD_data);

    % Remove outliers from HbO (don't have to do this step)
    [PD_data] = rem_HbX_outliers_10mins(PD_data,std_factor,t_pointsN);
        goodch_idx_step5_HbX_outliers = PD_data.goodch_idx;

    %currently only uses HbO to find outliers. 
    %%% Downsample before tomography (by a factor of fs, fs is 10Hz, so will DS to 1Hz)
    [PD_data] = PD_data_downsample(PD_data,fs);

    % %%% Step 6. Image Recon
    % w=0;
    % %use this 26 07 23
    % load('Jacobian_infant30week_850nm.mat')
    % load('mesh_infant30week_850nm.mat')
    [mesh] = plotjac_infant_mesh(mesh,PD_data,J,PD);
    
    %%% Tomography
    [PD_data] = tomography_PD_data_10mins(mesh,gmSurfaceMesh,J,PD_data,PD_time,PD,1,vol2gm,coverageThresh);
end

switch PD_data.eventType
            case "m_hypo"
                %find min glucose pos (index will corrospond to inddex in PD_time)
                %min_glucose_pos = find(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)))  );
                
                peak_glucose = min(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)));
                length_glucose = PD.events_length.m_hypo(eventN);

                %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                %make it -2mins before point and 8mins after =10mins total
                %%time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];
                %can change the (end) to (1) . (end) means it will always look at
                %the LATEST glucoce MINIMUM, (1) means it looks at the FIRST
                %glucose minimum
            case "S_hypo"
                %min_glucose_pos = find(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)))  );
                %%time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];

                peak_glucose = min(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)));
                length_glucose = PD.events_length.S_hypo(eventN);
         end

PD_data.peak_glucose = peak_glucose;
PD_data.length_glucose = length_glucose;

end