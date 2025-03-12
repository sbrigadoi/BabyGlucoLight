function [PD_data time_window_t] = get_timewindow_10mins(PD_data,PD,PD_time,time_window)

switch time_window %time window - time in minutes
    case "A" 
        time_window_t=[0 10];
    case "B"
        time_window_t=[13 23];

    case "C" % 5Min before min glucose
        switch PD_data.eventType
            case "m_hypo"
                %find min glucose pos (index will corrospond to inddex in PD_time)
                min_glucose_pos = find(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)))  );
                %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                %make it -2mins before point and 8mins after =10mins total
                time_window_t=[PD_time(min_glucose_pos(1)-1 )-2 PD_time(min_glucose_pos(1)-1 )+8];
            case "S_hypo"
                min_glucose_pos = find(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)))  );
                time_window_t=[PD_time(min_glucose_pos(1)-1 )-2 PD_time(min_glucose_pos(1)-1 )+8];
        end

    case "D" % min glucose
        switch PD_data.eventType
            case "m_hypo"
                %find min glucose pos (index will corrospond to inddex in PD_time)
                min_glucose_pos = find(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)))  );
                %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                %make it -2mins before point and 8mins after =10mins total
                time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];
                %can change the (end) to (1) . (end) means it will always look at
                %the LATEST glucoce MINIMUM, (1) means it looks at the FIRST
                %glucose minimum
            case "S_hypo"
                min_glucose_pos = find(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)))  );
                time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];
         end

    case "E" % 5 min post min glucose
        switch PD_data.eventType
            case "m_hypo"
                %find min glucose pos (index will corrospond to inddex in PD_time)
                min_glucose_pos = find(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)))  );
                %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                %make it -2mins before point and 8mins after =10mins total
                time_window_t=[PD_time(min_glucose_pos(end)+1 )-2 PD_time(min_glucose_pos(end)+1 )+8];
                %could change the (1) to (end) so that it's always looking at the
                %increase from the min
            case "S_hypo"
                 min_glucose_pos = find(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)))  );
                time_window_t=[PD_time(min_glucose_pos(end)+1 )-2 PD_time(min_glucose_pos(end)+1 )+8];
         end
    case "F" %start of end baseline
        time_window_t=[PD_time(end)-17 PD_time(end)-7];

    %case "G" %last point pre hypo glucose

    case "All"
        time_window_t=[PD_time(1) PD_time(end)];
    end

    PD_data.time_window_t = time_window_t;

    % Crop data. Just want 2mins before, 8 mins after (so total 10)
    fs=10;
    switch time_window
        case "All"
            PD_data.t = PD_data.t((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60),:);
            PD_data.s = PD_data.s((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60),:);
            PD_data.d = PD_data.d((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60),:);
            PD_data.aux = PD_data.aux((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60),:);

        otherwise
            PD_data.t = PD_data.t((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60)+1,:);
            PD_data.s = PD_data.s((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60)+1,:);
            PD_data.d = PD_data.d((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60)+1,:);
            PD_data.aux = PD_data.aux((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60)+1,:);
    end

end