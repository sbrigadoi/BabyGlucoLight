function [events_analysis_results,eu_analysis_results] = pax_selection_glycemic_analysis(PD,idx_PD,threshold_25_perc,result_disp,result_plot);
% Perform the events analysis and extract the euglycemia intervals
%
% Input  : PD                = struct with PD table, infogap and n_gap
%          idx_PD            = current PD index
%          threshold_25_perc = 25 percentile of all datset event duration,
%                              threshold for critical Euglycemia gap
%          result_disp       = 1 display result in Command Window
%          result_plot       = plot CGM results, 1 = <=R2022  3 = R2023
%
% Output : events_analysis_results = detailed events results + CGM signal
%                                    interpolated
%
%          eu_analysis_results     = detailed events results + critical gap
%                                    indexes

%% INITIALIZATION OF THE CGM DATA

% Extraction of the CGM parameters
PD_table = PD.PD;
cgm_value = PD_table.("mg/dL"); %[mg/dL]
cgm_date = PD_table.Time;

% Length of the final CGM signal (after gap correction)
n_cgm = length(cgm_value);

% Sampling Grid
t = [1:1:n_cgm]'; %[Sample]

% Time Grid
Ts_cgm = 5; %[min]
t_cgm = [0:Ts_cgm:Ts_cgm*(n_cgm-1)]';

%% CGM INTERPOLATION 

% In order to correct the non critical gap it is computed the 'CGM gap
% interpolation' suggested in:
% 'Minding the Gaps in Continuous Glucose Monitoring: A Method to Repair 
% Gaps to Achieve More Accurate Glucometrics' by Fonda et. al.

% Gap information
infogap = PD.infogap;

% Removal of samples within the time gaps
idx_to_delete = [];

for i = 1:1:PD.n_gap
    to_add        = [infogap.idx_start_gap(i):1:infogap.idx_end_gap(i)]';
    idx_to_delete = [idx_to_delete ; to_add];
end

% CGM data with time gaps
cgm_value_s                = cgm_value;
cgm_value_s(idx_to_delete) = [];
% Sampling grid
ts                = t;
ts(idx_to_delete) = [];

% Linear Interpolation (Fonda et. al.)
cgm_value_lin_interp = cgm_value;
for j=1:1:PD.n_gap
    gd_total = cgm_value(infogap.idx_end_gap(j)+1)-cgm_value(infogap.idx_start_gap(j)-1);
    for k = infogap.idx_start_gap(j):1:infogap.idx_end_gap(j)
        to_add = cgm_value_lin_interp(k-1)+(gd_total/infogap.n_samples(j));
        cgm_value_lin_interp(k) = to_add;
    end
end

%% GLYCEMIC EVENT DETECTION
% Initialization of support variables
th_eu_inf       = 72;  %[mg/dL]
th_eu_sup       = 144; %[mg/dL]
idx_start_event = [];  %Event start index
idx_end_event   = [];  %Event end index
idx_first       = [];  %First index out of eu range 
idx_last        = [];  %Last index out of eu range
f_eu_first      = 0;   %Flag for first absolute euglycemic sample
f_start         = 0;   %Flag event started
f_end           = 1;   %Flag event ended

a = 0; %Flag start finding event start index
b = 0; %Flag start finding event end index
p = 0; %Moving window support varibale for event start
q = 0; %Moving window support variable for event end

for i=1:1:n_cgm
    
    % Search for the event start point if an event is not already started
    % and the CGM signal went to euglycemia at least one time
    if (cgm_value(i)<th_eu_inf | cgm_value(i)>th_eu_sup) & f_start==0 & f_end==1        
        % Dimensional check, break if it is not possible to check the two
        % consecutive samples
        if i+2>n_cgm
            break
        end

        % Check if the two consecutive samples after the current are
        % outside the euglycemia range
        c = 0;
        for j=i+1:i+2
            if (cgm_value(j)<th_eu_inf | cgm_value(j)>th_eu_sup)
                c = c+1;
            end
        end
        
        % If the current sample and the next two are outside the euglycemia 
        % range, the search of the event starting point begins. 
        if c==2
            % Flag to point out that an event start and then it is
            % neccessary to find the event ending point
            f_start = 1;
            f_end   = 0;
            % Flag to start search the event starting point 
            a = 0;
            % Save of the current index as the first outside euglycemia
            % range
            idx_first = [idx_first; i];
            
            % Searching for the event starting point 
            while a==0
                
                % Dimensional check. If the moving window exceeds from the
                % CGM acquisition interval, the first event starting point
                % is outside the CGM acquisition. However, it is still 
                % considered valid since it satisfies the definition of
                % event. The starting point is placed at the first acquired 
                % sample.
                
                if i-p-4<1
                    idx_start_event = [idx_start_event; 1];                    
                    f_start = 1;
                    a = 1;
                    p = 0;
                    break
                end
                
                % Computation of the absolute CGM difference in a 15' 
                % moving window
                d_g = abs(cgm_value(i-p-4)-cgm_value(i-p-1));
                
                % If the absolute GCM difference is above 15 mg/dL and the
                % first sample of the moving window is in euglycemia, that
                % index is taken as the event starting point. Otherwise,
                % the moving window is shifted one step backwards

               if d_g<15 & cgm_value(i-p-4:1:i-p-1)>=th_eu_inf & cgm_value(i-p-4:1:i-p-1)<=th_eu_sup
                   
                   idx_start_event = [idx_start_event; i-p-4];
                   % Reset of moving window parameter
                   p = 0;
                   % Flag for the end of the searching of the event
                   % starting point
                   a = 1;
                else
                    % Backwards translation of the moving window
                    p = p+1;
                end
            end
        end
    end
    
    % Search for the event end point if an event is already started and the
    % CGM signal went to euglycemia at least one time
    if (cgm_value(i)>=th_eu_inf & cgm_value(i)<=th_eu_sup) & f_start==1 & f_end==0        
        % Dimensional check, break if it is not possible to check the two
        % consecutive samples. In this case, the last sample of the event
        % is the one before the first in euglycemia and the event end point
        % is the last one acquired.
        if i+2>n_cgm
            idx_last = [idx_last; i-1];
            idx_end_event = [idx_end_event; n_cgm];
            break
        end
        
        % Check if the two consecutive samples after the current are
        % inside the euglycemia range
        c = 0;
        for j=i+1:i+2
            if (cgm_value(j)>=th_eu_inf & cgm_value(j)<=th_eu_sup)
                c = c+1;
            end
        end
        
        % If the current sample and the next two are inside the euglycemia 
        % range, the search of the event ending point begins. 
        if c==2
            % Flag to point out that an event end and then, if necessary, 
            % it is required to find the start of a new event
            f_end = 1;
            f_start = 0;
            % Flag to start search the event ending point
            b = 0;
            
            % Save of the previous index as the last outside euglycemia 
            % range
            idx_last = [idx_last; i-1];
            
            % Searching for the event ending point 
            while b == 0
                
                % Dimensional check. If the moving window exceeds from the
                % CGM acquisition interval, the last sample is taken as
                % event starting point if it is not in euglycemia (the real 
                % ending point is outside the GCM acquisition)
                 if i+q+3>n_cgm
                    idx_end_event = [idx_end_event; n_cgm];
                    q = 0;
                    b = 1;
                    break
                 end

                % Computation of the absolute CGM difference in a 15' 
                % moving window. Pay attention that the i-th sample is
                % the first inside the euglycemia
                d_g = abs(cgm_value(i+q)-cgm_value(i+q+3));
                
                % If the absolute GCM difference is above 15 mg/dL and all
                % samples of the moving window is in euglycemia, that
                % index is taken as the event ending point. Otherwise,
                % the moving window is shifted one step forward

                if d_g<15 & (cgm_value(i+q:1:i+q+3)>=th_eu_inf & cgm_value(i+q:1:i+q+3)<=th_eu_sup)
                    
                    idx_end_event = [idx_end_event; i+q+3];
                    % Reset of the moving window parameter
                    q = 0;
                    % Flag for the end of the searching of the event
                    % ending point
                    b = 1;
                else 
                    % Farward translation of the moving window
                    q = q+1;
                end
            end
        end
    end
end

% If the dimension of the vector with the event end indexes is lower than
% the one of the vector with the event start indexes, it means that the
% event ends outside the CGM acquisition. Therefore, the index of the last
% sample is taken as the last point inside the event and as event end index
if length(idx_end_event) < length(idx_start_event)
    idx_end_event = [idx_end_event; n_cgm];
    idx_last      = [idx_last; n_cgm];
end

% The same consideration is not necessary with regard to the start of an 
% event since the constraints have already been implemented in the code
% (start of an event only in euglycemia, first dimensional check)

% Events that start and end at the same time are merged
to_remove = [];
for k=1:1:length(idx_start_event)-1
    if idx_start_event(k)==idx_start_event(k+1)
        to_remove = [to_remove k+1];
    end
end

idx_start_event(to_remove) = [];
idx_end_event(to_remove)   = [];
idx_first(to_remove)       = [];
idx_last(to_remove)        = [];

% If two consecutive events overlap, they are merged
to_remove_start = [];
to_remove_end   = [];
for i = 1:1:length(idx_start_event)-1
    if idx_start_event(i+1) <= idx_end_event(i)
        to_remove_start = [to_remove_start ; i+1];
        to_remove_end   = [to_remove_end ; i];
    end
end

idx_start_event(to_remove_start) = [];
idx_end_event(to_remove_end)     = [];
idx_first(to_remove_start)       = [];
idx_last(to_remove_end)          = [];

%% CRITICAL GAPS IDENTIFICATION AND CORRECTION

% Variables initialization
n_gap_tot        = length(infogap.idx_end_gap);
n_events         = length(idx_start_event);
date_start_event = PD_table.Time(idx_start_event);
date_end_event   = PD_table.Time(idx_end_event);
duration_event   = minutes(date_end_event-date_start_event); % [Min]
f_critic_event   = zeros(n_events,1); % 1 critical event 0 event ok
n_gap_corr_ev    = 0;

% Critical gap indentification and correction
cgm_value_lin_interp_correct   = cgm_value;
idx_critic_gap_start           = [];
idx_critic_gap_end             = [];

% Gap identification and correction inside each event
for i = 1:1:n_events
    
    % Start and end indexes of gap inside each event
    start_gap_event = [];
    end_gap_event   = [];
    c = 1;
    
    for j = 1:1:n_gap_tot
   
        % If there is no gap inside an event, check the next one.
        % Otherwise, the current gap is assumed all inside the event
        if (infogap.idx_start_gap(j) < idx_start_event(i) & infogap.idx_end_gap(j) < idx_start_event(i)) | (infogap.idx_start_gap(j) > idx_end_event(i) & infogap.idx_end_gap(j) > idx_end_event(i)) 
            continue
        else
            start_gap_event(c,1) = infogap.idx_start_gap(j);
            end_gap_event(c,1)   = infogap.idx_end_gap(j);
            c = c+1;
        end
        
    end
    
    % Erase possible 0 values
    start_gap_event(start_gap_event==0) = [];
    end_gap_event(end_gap_event==0)     = [];
    
    % Number gap inside the current event
    n_gap_tmp = length(start_gap_event);
    % The temporal duration of a gap is defined as the number of missing
    % samples multiplied by 5 minutes. 
    duration_gap         = (end_gap_event-start_gap_event+1)*5;  %[Min]
    duration_gap_perc    = (duration_gap/duration_event(i))*100; %[Perc]
    
    % Censoring events
    threshold_30  = 30; %[Perc]
    threshold_60  = 60; %[Min]
    threshold_cgm = 15; %[mg/dL]
    
    % Critical events criteria:
    % -) At least one gap lasts more than 30% of the duration of the event
    % -) At least one gap lasts more than 60 minutes and the absolute 
    %    difference between the cgm values before or after 30 minutes from 
    %    the start or end of the gap is greater than 15 mg/dL
    
    % Check for critical gap inside the current event
    for k=1:1:n_gap_tmp
        
        % Critical gap flag
        f_critic_gap = 0;
        
        % If the current gap lasts more than 30% of the duration of the 
        % event, the event and the gap are flagged as critical
        if duration_gap_perc(k)>threshold_30
            idx_critic_gap_start = [idx_critic_gap_start ; start_gap_event(k)];
            idx_critic_gap_end   = [idx_critic_gap_end ; end_gap_event(k)];
            f_critic_event(i)    = 1;
            f_critic_gap         = 1;
        end
        
        % Gap duration over 60'
        if duration_gap(k)>threshold_60
            
            % Dimensional check and absolute difference before gap
            if start_gap_event-6<0
                tmp_before = abs(cgm_value(start_gap_event)-cgm_value(1));
            else
                tmp_before = abs(cgm_value(start_gap_event)-cgm_value(start_gap_event-6));
            end
            
            % Dimensional check and absolute difference after gap
            if end_gap_event+6>n_cgm
                tmp_after = abs(cgm_value(n_cgm)-cgm_value(end_gap_event));
            else
                tmp_after = abs(cgm_value(end_gap_event+6)-cgm_value(end_gap_event));
            end
            
            % If the current gap is critical, the event is flagged.
            if tmp_after>15 | tmp_before>15
                idx_critic_gap_start = [idx_critic_gap_start ; start_gap_event(k)];
                idx_critic_gap_end   = [idx_critic_gap_end ; end_gap_event(k)];
                f_critic_event(i)    = 1;
                f_critic_gap         = 1;
            end
        end
        
        % If the current gap has an acceptable duration, it is replaced 
        % with the interpolation
        if f_critic_gap == 0
            cgm_value_lin_interp_correct(start_gap_event(k):1:end_gap_event(k)) = cgm_value_lin_interp(start_gap_event(k):1:end_gap_event(k));
            n_gap_corr_ev = n_gap_corr_ev + 1;
        end
    end
end

% Replacement of the original CGM signal with the correct one (Fonda et.
% al.)
cgm_value = cgm_value_lin_interp_correct;

% If an event starts or ends inside a critical gap, the start or the end is 
% shifted to the index of the first correct sample
for i=1:1:length(find(f_critic_event==1))
    
    % Critical events indexes
    idx_events_critic = find(f_critic_event==1);
    tmp = idx_events_critic(i);
    
    
    for k = 1:1:length(idx_critic_gap_start)
        
        % Current critical gap outside the current event
        if (idx_critic_gap_start(k) < idx_start_event(tmp) & idx_critic_gap_end(k) < idx_start_event(tmp)) | (idx_critic_gap_start(k) > idx_end_event(tmp) & idx_critic_gap_end(k) > idx_end_event(tmp)) 
            continue
        end
        
        % The current critical event ends inside a critical gap
        if (idx_start_event(tmp) <= idx_critic_gap_start(k)) & (idx_end_event(tmp) <= idx_critic_gap_end(k))
            idx_end_event(tmp) = idx_critic_gap_start(k)-1;
        end
        
        % The current critical event starts inside a critical gap
        if (idx_end_event(tmp) >= idx_critic_gap_end(k)) & (idx_start_event(tmp) >= idx_critic_gap_start(k))
            idx_start_event(tmp) = idx_critic_gap_end(k)+1;
        end
        
    end
end

% Control if CGM acquisition begins or ends inside an event
f_start_inside_ev = 0;
f_end_inside_ev   = 0;

% If the CGM signal starts inside an event, it is flagged
if idx_start_event(1,:) == 1
    
    f_start_inside_ev = 1;
    
    if f_critic_event(1) == 0
        f_critic_event(1) = 2;
    end
end

% If the CGM signal ends inside an event, it is flagged
if idx_end_event(end,:) == n_cgm 
    f_critic_event(end) = 2;
    f_end_inside_ev     = 1;
end

% Number of critical events 
n_events_critic  = length(find(f_critic_event==1));

%% EXTRACTION AND UPDATE OF USEFUL GLYCEMIC EVENT PARAMETERS 

n_events          = length(idx_start_event);
date_start_event  = PD_table.Time(idx_start_event);
date_end_event    = PD_table.Time(idx_end_event);
date_first_out    = PD_table.Time(idx_first);
date_last_out     = PD_table.Time(idx_last);
time_first_sample = timeofday(PD_table.Time(idx_first));
time_last_sample  = timeofday(PD_table.Time(idx_last));

duration_event            = minutes(date_end_event-date_start_event); % [Min]
duration_effective_out_eu = minutes(date_last_out-date_first_out);

duration_events_tot    = sum(duration_event);
median_duration_events = median(duration_event);

% These values ​​are defined between the first and last sample outisde the
% euglycemia range
max_value_event    = [];
min_value_event    = [];
range_value_event  = [];
mean_value_event   = [];
std_value_event    = [];
median_value_event = [];

for i=1:1:n_events
    
    max_value_event(i,1)    = max(cgm_value(idx_first(i):idx_last(i)));
    min_value_event(i,1)    = min(cgm_value(idx_first(i):idx_last(i)));
    range_value_event(i,1)  = max_value_event(i)-min_value_event(i);
    mean_value_event(i,1)   = round(mean(cgm_value(idx_first(i):idx_last(i))));
    median_value_event(i,1) = median(cgm_value(idx_first(i):idx_last(i)));
    std_value_event(i,1)    = std(cgm_value(idx_first(i):idx_last(i)));
    
end

%% EVENT CLASSIFICATION

% Type of events:
% -) Severe Hypo  = at least three sample under 47 mg/dL
% -) Severe Hyper = at least three sample over 180 mg/dL
% -) Mild Hypo    = [47 - 71] mg/dL and no Severe Hypo detect
% -) Mild Hyper   = [145 - 180] mg/dL and no Severe Hyper detect
% -) Hypo/Hyper   = if a hypo peak was detect righth before an hyper peak
%                   (or viceversa)

% Glycemic thresholds
th_mild_hypo_inf    = 47;
th_mild_hypo_sup    = 71;
th_severe_hypo_sup  = 46;
th_mild_hyper_inf   = 145;
th_mild_hyper_sup   = 180;
th_severe_hyper_inf = 181;

% Result vector
type_event = [];
n_severe_hypo  = 0;
n_mild_hypo    = 0;
n_severe_hyper = 0;
n_mild_hyper   = 0;
n_hypo_hyper   = 0;

% For each event
for i=1:n_events
    
    tmp_start_ev     = idx_start_event(i);
    tmp_end_ev       = idx_end_event(i);
    tmp_severe_hyper = 0;
    tmp_severe_hypo  = 0;
    tmp_mild_hypo    = 0;
    tmp_mild_hyper   = 0;
    
    % For each samples inside an event
    for k = tmp_start_ev:1:(tmp_end_ev-3)
        
        severe_hyper_count = 0;
        mild_hyper_count   = 0;
        severe_hypo_count  = 0;
        mild_hypo_count    = 0;
        
        % Test three consecutive samples
        for j = 1:1:3
            
            if cgm_value(k+j)>=th_severe_hyper_inf
                severe_hyper_count = severe_hyper_count + 1;
            end
            
            if cgm_value(k+j)<=th_severe_hypo_sup
                severe_hypo_count = severe_hypo_count + 1;
            end
            
            if cgm_value(k+j)>=th_mild_hypo_inf & cgm_value(k+j)<=th_mild_hypo_sup
                mild_hypo_count = mild_hypo_count + 1;
            end
            
            if cgm_value(k+j)>=th_mild_hyper_inf & cgm_value(k+j)<=th_mild_hyper_sup
                mild_hyper_count = mild_hyper_count + 1;
            end
            
        end
        
        % If at least three samples are in severe range, the event is
        % classified as severe 
        if severe_hyper_count == 3
            tmp_severe_hyper = tmp_severe_hyper + 1;
        end
        
        if severe_hypo_count == 3
            tmp_severe_hypo = tmp_severe_hypo + 1;
        end
        
        % Mild control in order to detect possible 'Hypo/Hyper'
        if mild_hyper_count == 1
            tmp_mild_hyper = tmp_mild_hyper + 1;
        end

        if mild_hypo_count == 1
            tmp_mild_hypo = tmp_mild_hypo + 1;
        end
    end
    
    % Hypo/Hyper classification
    if ((tmp_severe_hyper>0 | tmp_mild_hyper>0) & (tmp_severe_hypo>0 | tmp_mild_hypo>0))...
            | ((tmp_severe_hypo>0 | tmp_mild_hypo>0) & (tmp_severe_hyper>0 | tmp_mild_hyper>0))
        
        type_event   = [type_event ; string('Hypo / Hyper')];
        n_hypo_hyper = n_hypo_hyper + 1;
        continue
    end
    
    % Severe Hypo classification
    if tmp_severe_hypo > 0 
        type_event    = [type_event ; string('Severe Hypo')];
        n_severe_hypo = n_severe_hypo + 1;
        continue
    end
    
    % Severe Hyper classification
    if tmp_severe_hyper > 0 
        type_event     = [type_event ; string('Severe Hyper')];
        n_severe_hyper = n_severe_hyper + 1;
        continue
    end

    % Mild Hypo classification
    if tmp_mild_hypo > 0 
        type_event  = [type_event ; string('Mild Hypo')];
        n_mild_hypo = n_mild_hypo + 1;
    end

    % Mild Hyper classification
    if tmp_mild_hyper > 0 
        type_event   = [type_event ; string('Mild Hyper')];
        n_mild_hyper = n_mild_hyper + 1;
    end
end

%% SAVE OF THE GLYCEMIC EVENTS RESULTS

% Compact results display
event_pax_report = table(idx_PD,n_events,f_start_inside_ev,f_end_inside_ev,n_events_critic,n_mild_hypo,n_severe_hypo,n_mild_hyper,n_severe_hyper,n_hypo_hyper);

% Vector with progressing number of event
n_ev_prog = [1:1:n_events]';

% Table with glycemic events information
glycemic_events = table(n_ev_prog,f_critic_event,date_start_event,date_end_event,duration_event,type_event,max_value_event,...
    min_value_event,range_value_event,mean_value_event,median_value_event,std_value_event);

% Table with glycemic events useful index
idx_events = table(idx_start_event,idx_end_event,idx_first,idx_last);

%% EUGLYCEMIA INTERVAL

% The flag_event vector contains 1 if an index belong to an event and 0
% otherwise (1 = event)
flag_event = zeros(n_cgm,1);
for j = 1:1:n_events
    flag_event(idx_start_event(j):idx_end_event(j)) = 1;
end

% The flag_eu vector contains 1 if an index belongs to an euglycemic sample
% and 0 if it belongs to an event (1 = eu ; 0 = event)
flag_eu = double(not(flag_event));

% Extraction of start and end indexes of Euglycemia intervals
start        = 0;
idx_start_eu = [];
idx_end_eu   = [];

% Euglycemia interval identification
for k = 1:1:n_cgm
    
    % If the current flag value is 1 (EU) and a start of an euglycemia
    % interval is not already found, the current index is the start of the
    % euglycemia interval and until it is not found the end of that inerval
    % no other starting point are allowed
    if flag_eu(k) == 1 & start == 0
        idx_start_eu = [idx_start_eu ; k];
        start        = 1;
        continue
    end

    % If the current flag value is 0 (Event) and the start of an euglycemia
    % interval was already found, the previous index is the end of the
    % euglycemia interval and until it is not found the start of a new 
    % inerval no other ending point are allowed
    if flag_eu(k) == 0 & start == 1;
        idx_end_eu = [idx_end_eu ; k-1];
        start      = 0;
        continue
    end
        
end

% If the acquisition end in euglycemia, the last end index is added
% manually
if length(idx_start_eu)>length(idx_end_eu)
    idx_end_eu  = [idx_end_eu ; n_cgm];
end

%% GAP IDENTIFICATION INSIDE EUGLYCEMIC INTERVALS

% Identification of possible critical gap inside events 
n_eu_interval   = length(idx_start_eu);
date_start_eu   = PD_table.Time(idx_start_eu);
date_end_eu     = PD_table.Time(idx_end_eu);
duration_eu_int = minutes(date_end_eu-date_start_eu); % [Min]
check_eu        = ones(n_eu_interval,1);  % 0 delete 1 event ok
n_gap_corr_eu   = 0;

% Indexes of critical Euglycemia gap
idx_start_critical_gap_eu = [];
idx_end_critical_gap_eu   = [];

% The threshold is set equal to the 25-th percentile of the total duration 
% of glycemic events computes over all the dataset.
% The total duration of an event is defined as the time between the start
% and the end of an event excluding the first and last 15 minutes, which by
% definition must be in euglycemia. 

% Initialization of the new indexes vectors  
idx_start_eu_new = idx_start_eu;
idx_end_eu_new   = idx_end_eu;

% Critical gap indentification and correction
cgm_value_lin_interp_correct   = cgm_value;

% Test every gap 
for i = 1:1:n_eu_interval
    
    % Vector with start and end indexes of gap inside the current
    % euglycemic interval
    start_gap_eu = [];
    end_gap_eu   = [];
    c            = 1;
    
    % For every gap inside the CGM signal
    for j = 1:1:n_gap_tot
        
        % If there is no gap inside an euglycemia (EU) interval, check the
        % next one. Otherwise, the current gap is assumed all inside the EU
        % interval.
        if (infogap.idx_start_gap(j) < idx_start_eu(i) & infogap.idx_end_gap(j) < idx_start_eu(i)) | (infogap.idx_start_gap(j) > idx_end_eu(i) & infogap.idx_end_gap(j) > idx_end_eu(i)) 
            continue
        else
            start_gap_eu(c,1) = infogap.idx_start_gap(j);
            end_gap_eu(c,1)   = infogap.idx_end_gap(j);
            c                 = c+1;
        end
        
    end
    
    % Number gap inside the current EU interval
    n_gap_eu = length(start_gap_eu);
    % The temporal duration of a gap is defined as the number of missing
    % samples multiplied by 5 minutes
    duration_gap = (end_gap_eu-start_gap_eu+1)*5; %[Min]
    
    % Censoring EU interval
    for k=1:1:n_gap_eu
        
        % An EU interval is divided if it contains a gap with duration 
        % greater than the 25-th percentile of the EU total time.
        % Otherwise, if the gap is not critical inside the EU interval it 
        % is corrected by the relative portion of the interpolated CGM 
        % signal computed before.
        
        % Split of critical gap inside EU interval
        if duration_gap(k)>=threshold_25_perc
            
            % Save of Eu critical gap indexes
            idx_start_critical_gap_eu = [idx_start_critical_gap_eu ; start_gap_eu(k)];
            idx_end_critical_gap_eu   = [idx_end_critical_gap_eu ; end_gap_eu(k)];
            
            % Flag the current euglycemia interval as critical
            check_eu(i)  = 0;
            
            % Gap contains all the EU interval -> EU interval is completely
            % eliminated
            if start_gap_eu(k)<=idx_start_eu(i) & end_gap_eu(k)>=idx_end_eu(i)
                idx_start_eu_new(i) = NaN;
                idx_end_eu_new(i)   = NaN;
            end
            
            % If a gap is fully inside the EU interval -> EU interval is
            % divided in order to avoide the critical gap 
            if start_gap_eu(k)>idx_start_eu(i) & end_gap_eu(k)<idx_end_eu(i)
                idx_start_eu_new   = [idx_start_eu_new ; (end_gap_eu(k)+1)];
                idx_end_eu_new     = [idx_end_eu_new ; (start_gap_eu(k)-1)];
            end
            
            % If a gap starts outside the EU interval but end inside -> the
            % new divided EU interval starts immediately after the end of
            % the gap
            if start_gap_eu(k)<=idx_start_eu(i) & end_gap_eu(k)<idx_end_eu(i)
                idx_start_eu_new(i) = end_gap_eu(k)+1;
            end
            
            % If a gap starts inside the EU interval but end outside -> the
            % new divided EU interval ends immediately before the start of
            % the gap
            if start_gap_eu(k)>idx_start_eu(i) & end_gap_eu(k)>=idx_end_eu(i)
                idx_end_eu_new(i) = start_gap_eu(k)-1;
            end
        
        % Correction of the non critical gap 
        % (Linear Interpolation, Fonda et. al.)
        else
            cgm_value_lin_interp_correct(start_gap_eu(k):1:end_gap_eu(k)) = cgm_value_lin_interp(start_gap_eu(k):1:end_gap_eu(k));
            n_gap_corr_eu = n_gap_corr_eu + 1;
        end
    end
end

% Replacement of the original CGM signal with the correct one (Fonda et.
% al.)
cgm_value = cgm_value_lin_interp_correct;

% Number of critical EU gap identified
n_critical_gap_removed = length(find(check_eu==0));

% Duplication of the all EU indexes
idx_start_eu_raw = idx_start_eu;
idx_end_eu_raw   = idx_end_eu;

% Correction of the new EU indexes
idx_start_eu_new(isnan(idx_start_eu_new)) = [];
idx_end_eu_new(isnan(idx_end_eu_new))     = [];

% Substitution of the new EU indexes in the old variables
idx_start_eu = sort(idx_start_eu_new);
idx_end_eu   = sort(idx_end_eu_new);

%% EXTRACTION OF USEFULL EUGLYCEMIA INTERVAL PARAMETERS

% General Euglycemia Interval information
n_eu_interval = length(idx_start_eu);
date_start_eu = PD_table.Time(idx_start_eu);
date_end_eu   = PD_table.Time(idx_end_eu);

% Euglycemia Interval duration metrics
duration_eu_int    = minutes(date_end_eu-date_start_eu); %[Min]
duration_eu_tot    = sum(duration_eu_int); %[Min]
mean_duration_eu   = round(mean(duration_eu_int)); %[Min]
median_duration_eu = round(median(duration_eu_int)); %[Min]
MAD_duration_eu    = round(mad(duration_eu_int,1)); %0 = mean; 1 = median
max_duration_eu    = max(duration_eu_int); %[Min]
min_duration_eu    = min(duration_eu_int); %[Min]
 
% Number of 5' windows inside each Euglycemic interval and inside the first
% and last Euglycemia interval
dur_window_eu      = 5; %[Min]
n_windows_eu       = floor(duration_eu_int/dur_window_eu); 
n_windows_eu_tot   = sum(n_windows_eu);

% Glycemia value metrics
max_value_eu    = [];
min_value_eu    = [];
range_value_eu  = [];
mean_value_eu   = [];
std_value_eu    = [];
median_value_eu = [];

for i=1:1:n_eu_interval
    
    max_value_eu(i,1)    = max(cgm_value(idx_start_eu(i):idx_end_eu(i)));
    min_value_eu(i,1)    = min(cgm_value(idx_start_eu(i):idx_end_eu(i)));
    range_value_eu(i,1)  = max_value_eu(i)-min_value_eu(i);
    mean_value_eu(i,1)   = round(mean(cgm_value(idx_start_eu(i):idx_end_eu(i))));
    median_value_eu(i,1) = median(cgm_value(idx_start_eu(i):idx_end_eu(i)));
    std_value_eu(i,1)    = std(cgm_value(idx_start_eu(i):idx_end_eu(i)));
    
end

%% SAVE OF THE EUGLYCEMIC INTERVALS RESULTS

% Vector with progressing number of euglycemic intervals
n_eu_prog = [1:1:n_eu_interval]';

% Table with euglycemic interval information
euglycemic_int = table(n_eu_prog,date_start_eu,date_end_eu,duration_eu_int,max_value_eu,...
    min_value_eu,range_value_eu,mean_value_eu,median_value_eu,std_value_eu,n_windows_eu);

% Compact euglycemia interval results
eu_pax_report = table(idx_PD,n_eu_interval,duration_eu_tot,median_duration_eu,max_duration_eu,min_duration_eu,MAD_duration_eu,n_critical_gap_removed,n_gap_corr_eu);

% Table with euglycemic interval useful index
idx_eu = table(idx_start_eu,idx_end_eu);

% Table with critical gap indexes
idx_start_critic_gap = sort(unique([idx_critic_gap_start ; idx_start_critical_gap_eu]));
idx_end_critic_gap   = sort(unique([idx_critic_gap_end ; idx_end_critical_gap_eu]));
idx_critical_gap = table(idx_start_critic_gap,idx_end_critic_gap);

%% FINAL OUTPUT OF THE FUNCTION

% Events analysis output
events_analysis_results = struct('event_pax_report',event_pax_report,'glycemic_events',glycemic_events,'idx_events',idx_events,'cgm_value',cgm_value);

% Euglycemia analysis output
eu_analysis_results = struct('eu_pax_report',eu_pax_report,'euglycemic_int',euglycemic_int,'idx_eu',idx_eu,'idx_critical_gap',idx_critical_gap);

%% PLOT AND DISPLAY FINAL RESULTS

if result_disp == 1
    disp('==============================================================================')
    disp('PATIENT GLYCEMIC EVENTS REPORT')
    disp(' ')
    disp(event_pax_report)
    disp('==============================================================================')
    disp('GLYCEMIC EVENTS DETAILS ')
    disp(' ')
    disp(glycemic_events)
    disp('==============================================================================')
    disp('PATIENT EUGLYCEMIA INTERVALS REPORT')
    disp(' ')
    disp(eu_pax_report)
    disp('==============================================================================')
    disp('GLYCEMIC EVENTS DETAILS ')
    disp(' ')
    disp(euglycemic_int)
    disp('==============================================================================')
end

% MATLAB <= R2022
if result_plot == 1
    
    figure()
    
    hold on

    plot(t_cgm,cgm_value)

    yline(th_eu_inf,'--r'); yline(th_eu_sup,'--r')
    yline(th_mild_hypo_inf,'--r');yline(th_mild_hyper_sup,'--r')

    if length(idx_start_eu)>0
        plot(t_cgm(idx_start_eu),cgm_value(idx_start_eu),'g.','MarkerSize',20)
        plot(t_cgm(idx_end_eu),cgm_value(idx_end_eu),'gx','LineWidth',2)   
    end
    
    if length(find(f_critic_event==1))>0
        plot(t_cgm(idx_start_event(f_critic_event==1)),cgm_value(idx_start_event(f_critic_event==1)),'k.','MarkerSize',20)
        plot(t_cgm(idx_end_event(f_critic_event==1)),cgm_value(idx_end_event(f_critic_event==1)),'kx','LineWidth',2)
    end
    
    if length(find(f_critic_event==0))>=0
        plot(t_cgm(idx_start_event(f_critic_event==0)),cgm_value(idx_start_event(f_critic_event==0)),'r.','MarkerSize',20)
        plot(t_cgm(idx_end_event(f_critic_event==0)),cgm_value((idx_end_event(f_critic_event==0))),'rx','LineWidth',2)
    end
    
    f=get(gca,'Children');
    legend([f(2),f(1)],'Start Event','End Event')
    
    if width(PD.infogap)>0
        plot(t_cgm(PD.infogap.idx_start_gap),cgm_value(PD.infogap.idx_start_gap),'k^','LineWidth',2)
        if length(PD.infogap.idx_start_gap)>0
            f=get(gca,'Children');
            legend([f(1),f(3),f(2)],'Start Gap','Start Event','End Event')
        end
    end

    title(['CGM DATA PD ',num2str(idx_PD)])
    xlabel('Time [Min]')
    ylabel('[mg/dL]')
    
    if max(cgm_value)>240 & min(cgm_value)>30
        axis([0 t_cgm(end) 30 (max(cgm_value)+10)])

    elseif max(cgm_value)<240 & min(cgm_value)<30
        axis([0 t_cgm(end) min(cgm_value)-10 240])

    elseif max(cgm_value)>240 & min(cgm_value)<30
        axis([0 t_cgm(end) min(cgm_value)-10 (max(cgm_value)+10)])
        
    else
        axis([0 t_cgm(end) 30 240])
    end
end

% MATLAB R2023
if result_plot == 3
    
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));

    hold on

    plot(t_cgm,cgm_value,'k')

    % yline(th_eu_inf,'--r','LineWidth',2.5); yline(th_eu_sup,'--r','LineWidth',2.5)
    yline(th_eu_inf,'--r'); yline(th_eu_sup,'--r')


    yline(th_eu_inf,'--r'); yline(th_eu_sup,'--r')
    yline(th_mild_hypo_inf,'--r');yline(th_mild_hyper_sup,'--r')

    if length(find(f_critic_event==0))>0
        xregion(t_cgm(idx_start_event(f_critic_event==0)),t_cgm(idx_end_event(f_critic_event==0)),"FaceColor",'r')
    end

    if length(idx_start_eu)>0
        xregion(t_cgm(idx_start_eu),t_cgm(idx_end_eu),"FaceColor",'g')
    end

    if length(find(f_critic_event==1))>0
        xregion(t_cgm(idx_start_event(f_critic_event==1)),t_cgm(idx_end_event(f_critic_event==1)),"FaceColor",[0.9290 0.6940 0.1250])
    end

    if f_critic_event(1) == 2
        xregion(t_cgm(idx_start_event(1)),t_cgm(idx_end_event(1)))
    end

    if f_critic_event(end) == 2
        xregion(t_cgm(idx_start_event(end)),t_cgm(idx_end_event(end)))        
    end

    plot(t_cgm(idx_start_event),cgm_value(idx_start_event),'r.','MarkerSize',20)
    plot(t_cgm(idx_end_event),cgm_value(idx_end_event),'rx','LineWidth',2)
    f=get(gca,'Children');
    legend([f(2),f(1)],'Start Event','End Event')
    
    if width(PD.infogap)>0
        plot(t_cgm(PD.infogap.idx_start_gap),cgm_value(PD.infogap.idx_start_gap),'k^','LineWidth',2)
        if length(PD.infogap.idx_start_gap)>0
            f=get(gca,'Children');
            legend([f(1),f(3),f(2)],'Start Gap','Start Event','End Event')
        end
    end

    title(['CGM DATA PD ',num2str(idx_PD)])
    xlabel('Time [Min]')
    ylabel('[mg/dL]')
    % set(gca,'fontweigth','bold')
    % set(gca,'TickLabelInterpreter','none');
    % set(gca,'fontweight','bold','fontsize',14); 
    % h = gca;
    % h.YAxis.FontWeight = 'bold';
    % h.YAxis.FontSize = 18;
    
    % axesH = gca;
    % axesH.XAxis.TickLabelInterpreter = 'latex'
    % axesH.XAxis.TickLabelFormat      = '\\textbf{%g}';
    % yticks([72 144])
    ylabel('Glucose Concentration [mg/dL]','FontWeight','normal','FontSize',12)
    
    if max(cgm_value)>240 & min(cgm_value)>30
        axis([0 t_cgm(end) 30 (max(cgm_value)+10)])

    elseif max(cgm_value)<240 & min(cgm_value)<30
        axis([0 t_cgm(end) min(cgm_value)-10 240])

    elseif max(cgm_value)>240 & min(cgm_value)<30
        axis([0 t_cgm(end) min(cgm_value)-10 (max(cgm_value)+10)])
        
    else
        axis([0 t_cgm(end) 30 240])
        % axis([0 t_cgm(end) 30 160])
    end

    % Salvataggio HD dei plot

    % name = sprintf('PD%d.fig',idx_PD)
    % print(name,'-djpeg','-r600')
    % close all
    % savefig(name)

end

