function [BGC_metrics_PD_events]  = get_timewindow_10mins_BGCmetrics(good_PD_events_subN,good_PD_events,time_window,eventType,eventType_subN , eventType_subN_idx )

%glucose metrics
%for global we consider the entire single Event i.e time 0 to time end of
%event

%for local, we consider just the 10 minutes (3 samples) of BCG for the
%given A-F time window!

%variance
local_BCG_var = zeros(size(good_PD_events_subN,2),1);
global_BCG_var = zeros(size(good_PD_events_subN,2),1);

%max - min
local_BCG_range = zeros(size(good_PD_events_subN,2),1);
global_BCG_range = zeros(size(good_PD_events_subN,2),1);

%global minima
global_minima = zeros(size(good_PD_events_subN,2),1);

%Mean - local and global
local_BCG_mean = zeros(size(good_PD_events_subN,2),1);
global_BCG_mean = zeros(size(good_PD_events_subN,2),1);

%T val - local and global
local_BCG_T_val = zeros(size(good_PD_events_subN,2),1);
global_BCG_T_val = zeros(size(good_PD_events_subN,2),1);

%length of event

length_BGC_event = zeros(size(good_PD_events_subN,2),1);

t_val = zeros(size(good_PD_events_subN,2),1);
h_val_t = zeros(size(good_PD_events_subN,2),1);
ci_t = zeros(size(good_PD_events_subN,2),1);
pval_t = zeros(size(good_PD_events_subN,2),1);


for i=1:size(good_PD_events_subN,2)
    %load glucose data from subject N
    PD_N = load("H:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(good_PD_events_subN(i))+"");
    PD_N = getfield(PD_N,"PD"+num2str(good_PD_events_subN(i))+"");

    %get start and end time points for glucose event
    %PD_N.events_start_end.m_hypo(good_PD_events(i),:); %eventN = good_PD_events(i)

    %BCG values during entire event (GLOBAL BCG metrics)
    % NEED TO distinguish between m or S hypo
    switch eventType_subN_idx(i)
        case 1 %mhypo
            global_BCG = PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1) : PD_N.events_start_end.m_hypo(good_PD_events(i),2));
        case 2 %Shypo
            global_BCG = PD_N.glucose(PD_N.events_start_end.S_hypo(good_PD_events(i),1) : PD_N.events_start_end.S_hypo(good_PD_events(i),2));
        case 3 %mhyper
            global_BCG = PD_N.glucose(PD_N.events_start_end.m_hyper(good_PD_events(i),1) : PD_N.events_start_end.m_hyper(good_PD_events(i),2));
        case 4 %Shyper
            global_BCG = PD_N.glucose(PD_N.events_start_end.S_hyper(good_PD_events(i),1) : PD_N.events_start_end.S_hyper(good_PD_events(i),2));
    end

    global_BCG_var(i,1) = var(global_BCG);
    if eventType == "m_hypo"
        global_BCG_range(i,1) = range(global_BCG)*-1; %make it minus for hypo 
    elseif eventType == "S_hypo"
        global_BCG_range(i,1) = range(global_BCG)*-1; %make it minus for hypo 
    else 
        global_BCG_range(i,1) = range(global_BCG) %keep it plus for hyper
    end
    global_BCG_mean(i,1) = mean(global_BCG);
    global_minima(i,1) = min(global_BCG);

    length_BGC_event(i,1) = size(global_BCG,1);

        [t_test(i,1),pval_t(i,1)] = ttest(global_BCG,global_BCG(1));
        [h_val_t(i,1),pval_t(i,1),ci_t,stats] = ttest(global_BCG,global_BCG(1));
        global_BCG_T_val(i,1) = stats.tstat;

    

    %Get infex of time window
    switch time_window %time window - time in minutes
        case "A"
            time_window_t=[0 10];
            %time_window_t_BCG=[1 3]; %sample 1 to 3 covers first 10 mins 0,5,10 min
           time_window_t_BCG=[1 4]; %sample 1 to 3 covers first 15 mins 0,5,10,15 min

            %time_window_t_BCG=[2 4]; %sample 1 to 3 covers first 10 mins 0,5,10 min

        case "B"
            time_window_t=[13 23]; %sample 1 to 3 covers first 10 mins 0,5,10 min
            time_window_t_BCG =[3 5]; %sample 3 to 5 covers 10 to 20 mins 10,15,20 min

        case "C" % 5Min before min glucose
            switch eventType_subN_idx(i)
                case 1 %mhypo
                    %find min glucose pos (index will corrospond to inddex in PD_time)
                    min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)))  );
                    %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                    %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                    %make it -2mins before point and 8mins after =10mins total
                    %time_window_t=[PD_time(min_glucose_pos(1)-1 )-2 PD_time(min_glucose_pos(1)-1 )+8];
                    
                    time_window_t_BCG =[min_glucose_pos-1 min_glucose_pos+1]; %time sample of min glucose pos -1 (5mins)
                case 2 %shypo
                    min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.S_hypo(good_PD_events(i),1):PD_N.events_start_end.S_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.S_hypo(good_PD_events(i),1):PD_N.events_start_end.S_hypo(good_PD_events(i),2)))  );
                    %time_window_t=[PD_time(min_glucose_pos(1)-1 )-2 PD_time(min_glucose_pos(1)-1 )+8];
                    time_window_t_BCG =[min_glucose_pos-1 min_glucose_pos+1]; %time sample of min glucose pos -1 (5mins)

            end

        case "D" % min glucose
            switch eventType_subN_idx(i)
                case 1 %mhypo
                    %find min glucose pos (index will corrospond to inddex in PD_time)
                    min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)))  );
                    %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                    %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                    %make it -2mins before point and 8mins after =10mins total
                    
                    %time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];
                    
                    %time_window_t_BCG=[min_glucose_pos(end) min_glucose_pos(end)+2];%min to min plus 10 mins 

                    time_window_t_BCG=[min_glucose_pos(end) min_glucose_pos(end)+3]; %min to min plus 15 mins
                    %can change the (end) to (1) . (end) means it will always look at
                    %the LATEST glucoce MINIMUM, (1) means it looks at the FIRST
                    %glucose minimum
                case 2 %shypo
                    min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.S_hypo(good_PD_events(i),1):PD_N.events_start_end.S_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.S_hypo(good_PD_events(i),1):PD_N.events_start_end.S_hypo(good_PD_events(i),2)))  );
                    %time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];
                    %time_window_t_BCG=[min_glucose_pos(end) min_glucose_pos(end)+2];
                    time_window_t_BCG=[min_glucose_pos(end) min_glucose_pos(end)+3];

            end

        case "E" % 5 min post min glucose
            switch eventType(i)
                case 1 %"m_hypo"
                    %find min glucose pos (index will corrospond to inddex in PD_time)
                    min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)))  );
                    %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                    %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                    %make it -2mins before point and 8mins after =10mins total
                    %time_window_t=[PD_time(min_glucose_pos(end)+1 )-2 PD_time(min_glucose_pos(end)+1 )+8];
                    time_window_t_BCG =[min_glucose_pos+1 min_glucose_pos+3]; %time sample of min glucose pos -1 (5mins)

                    %could change the (1) to (end) so that it's always looking at the
                    %increase from the min
                case 2% "S_hypo"
                    min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.S_hypo(good_PD_events(i),1):PD_N.events_start_end.S_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.S_hypo(good_PD_events(i),1):PD_N.events_start_end.S_hypo(good_PD_events(i),2)))  );
                    %time_window_t=[PD_time(min_glucose_pos(end)+1 )-2 PD_time(min_glucose_pos(end)+1 )+8];
                    time_window_t_BCG =[min_glucose_pos+1 min_glucose_pos+3]; %time sample of min glucose pos -1 (5mins)

            end
        case "F" %start of end baseline
            %time_window_t=[PD_time(end)-17 PD_time(end)-7];
            %incorrectF, similar to G
            %time_window_t_BCG =[size(global_BCG,1)-2 size(global_BCG,1)]; %last-2 lastlast time sample %time win F on 30 10 24 was this by mistake
            
            %correct F
            %time_window_t_BCG =[size(global_BCG,1)-3 size(global_BCG,1)-1]; %last-3 last-1 last time sample
            
            %F to end
            time_window_t_BCG =[size(global_BCG,1)-3 size(global_BCG,1)]; %last-3 last last time sample

            %F(-20 mins from end to 5mins from end)
            %time_window_t_BCG =[size(global_BCG,1)-4 size(global_BCG,1)-1]; %last-3 last last time sample


            %case "G" %last point pre hypo glucose
        case "G"
            time_window_t_BCG =[size(global_BCG,1)-2 size(global_BCG,1)]

        case "All"
            %time_window_t=[PD_time(1) PD_time(end)];
            time_window_t_BCG =[1 size(global_BCG,1)]; %1: all time sample N

    end


        %global_BCG(time_window_t_BCG(1):time_window_t_BCG(2))
        local_BCG_var(i,1) = var(global_BCG(time_window_t_BCG(1):time_window_t_BCG(2)));
        %local_BCG_range(i,1) = range(global_BCG(time_window_t_BCG(1):time_window_t_BCG(2)));

        if global_BCG(time_window_t_BCG(1)) > global_BCG(time_window_t_BCG(2))
            local_BCG_range(i,1) = range(global_BCG(time_window_t_BCG(1):time_window_t_BCG(2)))*-1; %make it minus for decreasing local BGC range 
        else 
            local_BCG_range(i,1) = range(global_BCG(time_window_t_BCG(1):time_window_t_BCG(2))); %keep it plus for for increasing local BGC range 
        end

        local_BCG_mean(i,1) = mean(global_BCG(time_window_t_BCG(1):time_window_t_BCG(2)));
        
        
        %if the four samples in the local BGC are the same, set tval as 0.
        if size(unique(global_BCG(time_window_t_BCG(1):time_window_t_BCG(2))),1)==1
            %ttest_vals = [global_BCG(time_window_t_BCG(1):time_window_t_BCG(2)-1) ; global_BCG(time_window_t_BCG(1))-1];
            local_BCG_T_val(i,1) = 0;
        else 
            ttest_vals = global_BCG(time_window_t_BCG(1):time_window_t_BCG(2));
            [t_test(i,1),pval_t(i,1)] = ttest(ttest_vals,ttest_vals(1) );
            [h_val_t(i,1),pval_t(i,1),ci_t,stats] = ttest(ttest_vals,ttest_vals(1));
            local_BCG_T_val(i,1) = stats.tstat;

        end

        %ttest_vals = global_BCG(time_window_t_BCG(1):time_window_t_BCG(2));


        




    % PD_data.time_window_t = time_window_t;
    % 
    % % Crop data. Just want 2mins before, 8 mins after (so total 10)
    % fs=10;
    % switch time_window
    %     case "All"
    %         PD_data.t = PD_data.t((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60),:);
    %         PD_data.s = PD_data.s((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60),:);
    %         PD_data.d = PD_data.d((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60),:);
    %         PD_data.aux = PD_data.aux((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60),:);
    % 
    %     otherwise
    %         PD_data.t = PD_data.t((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60)+1,:);
    %         PD_data.s = PD_data.s((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60)+1,:);
    %         PD_data.d = PD_data.d((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60)+1,:);
    %         PD_data.aux = PD_data.aux((time_window_t(1)*fs*60)+1:(time_window_t(2)*fs*60)+1,:);
    % end






end

BGC_metrics_PD_events.local_BCG_var = local_BCG_var;
BGC_metrics_PD_events.global_BCG_var = global_BCG_var;
BGC_metrics_PD_events.local_BCG_range = local_BCG_range;
BGC_metrics_PD_events.global_BCG_range = global_BCG_range;

BGC_metrics_PD_events.global_minima = global_minima;


BGC_metrics_PD_events.global_BCG_mean = global_BCG_mean;
BGC_metrics_PD_events.local_BCG_mean = local_BCG_mean;

BGC_metrics_PD_events.local_BCG_T_val = local_BCG_T_val;
BGC_metrics_PD_events.global_BCG_T_val = global_BCG_T_val;

BGC_metrics_PD_events.length_BGC_event = length_BGC_event;



end