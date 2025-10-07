function [eu_analysis_results] = pax_selection_euglycemia_analysis(PD,events_analysis_results,threshold_25_perc,idx_PD,result_disp,result_plot)


%% INITIALIZATION OF THE CGM DATA

% Extraction of the CGM parameters
PD_table = PD.PD;
cgm_value = PD_table.("mg/dL"); %[mg/dL]
cgm_date = PD_table.Time;

% idx_start_event = events_analysis_results.idx_raw_events.idx_start_event_raw;
% idx_end_event = events_analysis_results.idx_raw_events.idx_end_event_raw;
idx_start_event = events_analysis_results.idx_events.idx_start_event;
idx_end_event = events_analysis_results.idx_events.idx_end_event;
% idx_first_eu = events_analysis_results.idx_first_eu;

% Length of the final CGM signal (after gap correction)
n_cgm = length(cgm_value);

% Sampling Grid
t = [1:1:n_cgm]'; %[Sample]

% Time Grid
Ts_cgm = 5; %[min]
t_cgm = [0:Ts_cgm:Ts_cgm*(n_cgm-1)]';

%% CGM INTERPOLATION 

% In order to correct the non critic gap it is computed the 'CGM gap
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

%% EUGLYCEMIA INTERVAL

% The flag_event vector contains 1 if an index belong to an event and 0
% otherwise (1 = event)

flag_event   = zeros(n_cgm,1);
n_events_raw = length(idx_start_event);

for j = 1:1:n_events_raw
    flag_event(idx_start_event(j):idx_end_event(j)) = 1;
end

% The flag_eu vector contains 1 if an index belongs to an euglycemic sample
% and 0 if it belongs to an event (1 = eu ; 0 = event)
flag_eu = double(not(flag_event));

% If acquisition starts outside euglycemia
% for i = 1:(idx_first_eu-1)
%     flag_eu(i) = 0;
% end

% Extraction of start and end indexes of Euglycemia intervals
start        = 0;
idx_start_eu = [];
idx_end_eu   = [];

for k = 1:1:n_cgm
    
    % If the current flag value is 1 (eu) and a start of an euglycemia
    % interval is not already found, the current index is the start of the
    % euglycemia interval and until it is not found the end of that inerval
    % no other starting point are allowed
    if flag_eu(k) == 1 & start == 0
        idx_start_eu = [idx_start_eu ; k];
        start        = 1;
        continue
    end

    % If the current flag value is 0 (event) and the start of an euglycemia
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

%% CRITICAL GAP IDENTIFICATION INSIDE EUGLYCEMIC INTERVALS

% The threshold is set equal to the 25-th percentile of the total duration 
% of glycemic events computes over all the dataset.
% The total duration of an event is defined as the time between the start
% and the end of an event excluding the first and last 15 minutes, which by
% definition must be in euglycemia. 

% Identification of possible critical gap inside events 
n_gap_tot       = length(infogap.idx_end_gap);
n_eu_interval   = length(idx_start_eu);
date_start_eu   = PD_table.Time(idx_start_eu);
date_end_eu     = PD_table.Time(idx_end_eu);
duration_eu_int = minutes(date_end_eu-date_start_eu); % [Min]
check_eu        = ones(n_eu_interval,1);  % 0 delete 1 event ok
n_gap_corr_eu   = 0;

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
        % greater than the 25-th percentile of the EU total time or is 
        % greater than 60'. Otherwise, if the gap is not critical inside 
        % the EU interval it is corrected by the relative portion of the 
        % interpolated CGM signal computed before.
        
        % Split of critical gap inside EU interval
        if duration_gap(k)>=threshold_25_perc
            
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
            % new divided EU interval end immediately before the start of
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

% Function output
eu_analysis_results = struct('eu_pax_report',eu_pax_report,'euglycemic_int',euglycemic_int,'idx_eu',idx_eu);

%% DISPLAY OF EUGLYCEMIA ANALYSIS RESULTS

if result_disp == 1
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
%% PLOT OF THE FINAL RESULTS

th_eu_inf           = 72;  %[mg/dL]
th_eu_sup           = 144; %[mg/dL]
th_mild_hypo_inf    = 47;  %[mg/dL]
th_mild_hypo_sup    = 71;  %[mg/dL]
th_severe_hypo_sup  = 46;  %[mg/dL]
th_mild_hyper_inf   = 145; %[mg/dL]
th_mild_hyper_sup   = 180; %[mg/dL]
th_severe_hyper_inf = 181; %[mg/dL]

if result_plot == 1
    % Risultati finali con rimozione eventi e correzione eu (non interp)
    figure()

    hold on

    plot(t,cgm_value)

    yline(th_eu_inf,'--r')
    yline(th_eu_sup,'--r')
    yline(47,'--r');yline(72,'--r');yline(144,'--r');yline(180,'--r')

    plot(idx_start_event,cgm_value(idx_start_event),'r.','MarkerSize',20)
    plot(idx_end_event,cgm_value(idx_end_event),'rx','LineWidth',2)
    
    plot(idx_start_eu,cgm_value(idx_start_eu),'k.','MarkerSize',20)
    plot(idx_end_eu,cgm_value(idx_end_eu),'kx','LineWidth',2)

    if width(PD.infogap)>0
        plot(t(PD.infogap.idx_start_gap),cgm_value(PD.infogap.idx_start_gap),'k^','LineWidth',2)
    end

    title(['CGM DATA PD ',num2str(idx_PD)])
    xlabel('Time [Min]')
    ylabel('[mg/dL]')
    f=get(gca,'Children');
    legend([f(5),f(4),f(3),f(2),f(1)],'Start Event','End Event','Start Eu','End Eu','Start gap')
    axis([0 length(t) 40 240])
end