function [PD_data] = view_all_events_subject(subjectN,eventType,time_window)
savefig = 1;

fs = 10;
PD_data.fs = fs;

PD_data.subjectN = subjectN;
%PD_data.eventN = eventN;
PD_data.eventType = eventType;
PD_data.time_window = time_window;

%PD 10 1:1
%PD 13 5:5
%PD 15 1:8
%PD 20 3:12
%PD 25 1:5
%PD33 1:16
% file_numbers = dir('E:\PadovaPostDoc\BabyGluCo\NIRS_data\PD33\formatted');
% start_of_NIRS = strfind(cellstr(NIRS_file_names_char),'.nirs');
% E:\PadovaPostDoc\BabyGluCo\NIRS_data\PD33\formatted

switch PD_data.subjectN
    case 7 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [8 9 10 11 12 13 14 15 16 17 18 19];
        end
    case 8 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2 3 4 5 6 7 8 9 10 11];
        end
    case 10 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1];    
        end
    case 13 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [5];
        end
    case 15 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2 3 4 5 6 7 8];
            case "S_hypo"
                PD_S_hypo_events = [1 2]; %S hypo 
        end
    case 20 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [3 4 5 6 7 8 9 10 11 12];
        end
    case 25 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2 3 4 5];
            case "S_hypo"
                PD_S_hypo_events = [1]; %S hypo 
        end
    case 28 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2];    
        end
    case 33 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
            case "S_hypo"
                PD_S_hypo_events = [1 2 3]; %S hypo 
        end
    case 34 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1];    
        end

    case 39 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2];
            
        end
    case 41 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [2 3 4];
            
        end
    case 44 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1];
            
        end
    case 45 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2];
            
        end

    case 48 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2];
            
        end

    case 49 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];
            case "S_hypo"
                PD_S_hypo_events = [1 2 3 4 5]; %S hypo 
        end

    case 55 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2 3 4];
            
                
        end

    case 56 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2 3 4 5];    
            case "S_hypo"
                PD_S_hypo_events = [1];    
        end

    case 58 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1];    
        end

    case 59 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2];    
        end

    case 60 %add other subject Ns
        switch PD_data.eventType
            case "m_hypo"
                PD_m_hypo_events = [1 2 3 4 5 6 7 8 9 10 11];    
        end

end

switch eventType
    case "m_hypo"
        N_events = size(PD_m_hypo_events,2);
    case "S_hypo"
        N_events = size(PD_S_hypo_events,2);
end



subplot_col = 5;
if N_events < 5
   subplot_col =  N_events;
end

subplot_row = ceil(N_events/subplot_col);

figure()
for i=1:N_events
    switch eventType
        case "m_hypo"
            eventN = PD_m_hypo_events(i);
        case "S_hypo"
            eventN = PD_S_hypo_events(i);
    end
    
    PD_data = load("H:\PadovaPostDoc\BabyGluCo\NIRS_data\PD"+num2str(subjectN)+"\formatted\PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(eventN)+".nirs",'-mat');
    PD_data = getfield(PD_data,"PD"+num2str(subjectN)+"_"+eventType+"_X");
    PD_data.eventN = eventN;
    PD_data.subjectN = subjectN;
    PD_data.time_window = time_window;
    PD_data.eventType = eventType;
    
    PD = load("H:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(subjectN)+".mat");
    switch PD_data.subjectN
        case 7
            PD=PD.PD7;
        case 8
            PD=PD.PD8;
        case 10
            PD = PD.PD10;
        case 13
            PD=PD.PD13;
        case 15
            PD=PD.PD15;
        case 20
            PD=PD.PD20;
        case 25
            PD=PD.PD25;
        case 28
            PD = PD.PD28;
        case 33         %%add other subject Ns
            PD = PD.PD33; %change this
        case 34
            PD = PD.PD34;
        case 39
            PD = PD.PD39;
        case 41
            PD = PD.PD41;
        case 44
            PD = PD.PD44;
         case 45
            PD = PD.PD45;
        case 48
            PD = PD.PD48;
        case 49         %%add other subject Ns
            PD = PD.PD49;
        case 55         %%add other subject Ns
            PD = PD.PD55;
        case 56
            PD = PD.PD56;
        case 58
            PD = PD.PD58;
        case 59
            PD = PD.PD59;
        case 60
            PD = PD.PD60;

    end


    switch eventType
        case "m_hypo"
            PD_time = 0:5:5*PD.events_length.m_hypo(eventN);
        case "S_hypo"
            PD_time = 0:5:5*PD.events_length.S_hypo(eventN);
    end

    %get time window
    [PD_data time_window_t] = get_timewindow_10mins(PD_data,PD,PD_time,time_window);

    % Step 2. Signal Qaulity check
    dRange = [5E-4 3];%3];
    SNRrange = 2; %2;
    [goodch_idx, PD_data] = find_good_ch_15mins(PD_data,dRange,SNRrange,0);
    if size(goodch_idx,1) <2 %less than 2good chs, i.e 0 or 1
        goodch_idx = [1;65];
    end
    PD_data.goodch_idx = goodch_idx;


    subplot(subplot_row,subplot_col,i)
    plot(PD_data.t,PD_data.d(:,goodch_idx))
    xlabel('Time / S')
    ylabel('D / A.U')
    title("PD "+PD_data.subjectN+" "+PD_data.eventType+" "+num2str(eventN)+" N good ch = "+num2str(size(goodch_idx,1))+" TimeW "+time_window+" ")
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD_data.subjectN)+"_"+PD_data.eventType+"_"+num2str(PD_data.eventN)+"_TimeW_"+PD_data.time_window+"_raw_data.png")
end


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

%plot all the glucose events
figure()
for i=1:N_events
    switch eventType
        case "m_hypo"
            eventN = PD_m_hypo_events(i);
        case "S_hypo"
            eventN = PD_S_hypo_events(i);
    end
    PD_data = load("H:\PadovaPostDoc\BabyGluCo\NIRS_data\PD"+num2str(subjectN)+"\formatted\PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(eventN)+".nirs",'-mat');
    PD_data = getfield(PD_data,"PD"+num2str(subjectN)+"_"+eventType+"_X");
    PD_data.eventN = eventN;
    PD_data.subjectN = subjectN;
    PD_data.time_window = time_window;
    PD_data.eventType = eventType;
    PD = load("H:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(subjectN)+".mat");
    switch PD_data.subjectN
        case 7
            PD=PD.PD7;
        case 8
            PD=PD.PD8;
        case 10
            PD=PD.PD10;
        case 13
            PD=PD.PD13;
        case 15
            PD=PD.PD15;
        case 20
            PD=PD.PD20;
        case 25
            PD=PD.PD25;
        case 28
            PD=PD.PD28;
        case 33         %%add other subject Ns
            PD = PD.PD33; %change this
        case 34
            PD=PD.PD34;
        case 39
            PD = PD.PD39;
        case 41
            PD = PD.PD41;
        case 44
            PD = PD.PD44;
        case 45
            PD = PD.PD45;
         case 48
            PD = PD.PD48;
        case 49         %%add other subject Ns
            PD = PD.PD49; %change this
        case 55         %%add other subject Ns
            PD = PD.PD55; %change this
        case 56
            PD=PD.PD56;
        case 58
            PD=PD.PD58;
        case 59
            PD=PD.PD59;
        case 60
            PD=PD.PD60;

    end
    %PD_time = 0:5:5*PD.events_length.m_hypo(eventN);
    switch eventType
        case "m_hypo"
            PD_time = 0:5:5*PD.events_length.m_hypo(eventN);
        case "S_hypo"
            PD_time = 0:5:5*PD.events_length.S_hypo(eventN);
    end
    %get time window
    [PD_data time_window_t] = get_timewindow_10mins(PD_data,PD,PD_time,time_window);
    subplot(subplot_row,subplot_col,i)
    switch eventType
        case "m_hypo"
            plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(eventN,1):PD.events_start_end.m_hypo(eventN,2)),'b-o')
        case "S_hypo"
            plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(eventN,1):PD.events_start_end.S_hypo(eventN,2)),'b-o')
    end
    title("Time W"+time_window+" PD "+num2str(subjectN)+" "+eventType+" "+num2str(eventN)+". nirs")
    xlabel("Time / mins")
    ylabel("BGC / mg/DL")
    xline(15,'g--');xline(PD_time(end)-15,'r--');
    xline(time_window_t(1),'m--');xline(time_window_t(2),'m--');
    yline(72,'k--');
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD_data.subjectN)+"_"+PD_data.eventType+"_"+num2str(PD_data.eventN)+"_TimeW_"+PD_data.time_window+"_raw_glucose.png")
end


end