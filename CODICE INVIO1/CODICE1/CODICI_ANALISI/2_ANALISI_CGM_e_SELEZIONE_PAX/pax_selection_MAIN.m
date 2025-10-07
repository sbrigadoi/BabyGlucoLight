%% PATIENT SELECTION

close all
clear all
clc

%% SETTING OF FORLDERS PATH 

base_path     = pwd;
code_path     = fullfile(base_path,'CODICE','DEFINITIVO','CGM_SELEZIONE_PAX_NIRS');
data_path     = fullfile(base_path,'DATI','CGM_INTERP/');
function_path = fullfile(base_path,'CODICE','DEFINITIVO','CGM_SELEZIONE_PAX_NIRS');

addpath(genpath(pwd))
addpath(code_path)
addpath(data_path)
addpath(function_path)

%% ALL-DATASET GLYCEMIC EVENTS ANALYSIS

file_name      = 'PD%d.mat';
pax_idx_vector = [1:1:17 19:1:20 22:1:26 28 33:1:45 47:1:55 59]';
n_pax          = length(pax_idx_vector);

% Glycemic Analysis: Display Setting
result_disp       = 0;  % 1 = show results on command window for each PD
result_plot       = 3;  % 1 = show plot R2021  3 = show plot R2023 for each PD
threshold_25_perc = 45; % [Min] Critical Euglycemia Gap threshold % RULE 1

% Display of the preliminar information
disp('ALL-DATASET GLICEMIC EVENTS ANALYSIS')
disp(' ')
disp(['Number of patients: ',num2str(n_pax)])
disp(' ')
disp('START')
disp(' ')

% Glycemic events analysis for each patient inside the dataset

% Variables initialization
time_before_first_eu = [];
pax_event_report     = [];
idx_start_first_eu_w   = [];
idx_end_first_eu_w     = [];
idx_start_last_eu_w    = [];
idx_end_last_eu_w      = [];
TIR_first = [];
TIR_last  = [];
TIR_eu    = [];
TOR_eu    = [];
TIR_eu_perc = [];
TOR_eu_perc = [];
n_critical_event = [];
n_critical_eu    = [];
n_critic_gap     = [];
duration_cr_gap  = [];
f_first_eu_gap   = [];
f_last_eu_gap    = [];
median_t_cr_gap  = [];
max_t_cr_gap = [];
min_t_cr_gap = [];
dur_gap_perc = [];
duration_cgm = [];
dur_gap      = [];
dur_CGM      = [];
TAR = [];
TBR = [];
pax_eu_report = [];
%%
for i = 1:1:n_pax
    
    idx_pax = pax_idx_vector(i);
    disp('==============================================================================')
    disp(['PATIENT ',num2str(idx_pax)])
    
    % Load of the current patient data
    to_load     = sprintf(file_name,idx_pax);
    
    PD_struct   = load(to_load);
    table_name  = string(fieldnames(PD_struct));
    PD          = getfield(PD_struct,table_name);
    
    % GLYCEMIC ANALYSIS
    [events_analysis_results,eu_analysis_results] = pax_selection_glycemic_analysis(PD,idx_pax,threshold_25_perc,result_disp,result_plot);
    
    pax_event_report = [pax_event_report ; events_analysis_results.event_pax_report];
    pax_eu_report    = [pax_eu_report ; eu_analysis_results.eu_pax_report];

    % SAVE OF THE CGM SIGNAL INTERPOLATED %was commented out pre 07 03 25
     %PD.PD.("mg/dL") = events_analysis_results.cgm_value(:);
     %save_PD_path_name   = fullfile(base_path,'DATI','CGM_INTERP',to_load);
     %save (save_PD_path_name,'PD')
    
    % PATIENT SELECTION METRICS
    
    % Duration of total CGM acquisition
    tmp_dur = minutes(PD.PD.Time(end)-PD.PD.Time(1));
    dur_CGM = [dur_CGM ; tmp_dur];
    
    % Time elapsed between the start of the acquisition and the first 
    % Euglycemia interval
    time_before_first_eu = [time_before_first_eu ; minutes(eu_analysis_results.euglycemic_int.date_start_eu(1) - PD.PD.Time(1))];
    
    % saveas(gcf,['PD' num2str(idx_pax) '.fig']);
    
    % First and last Eu interval indexes
    idx_start_first_eu_w = [idx_start_first_eu_w ; eu_analysis_results.idx_eu.idx_start_eu(1)];
    idx_end_first_eu_w   = [idx_end_first_eu_w ; eu_analysis_results.idx_eu.idx_end_eu(1)];
    idx_start_last_eu_w  = [idx_start_last_eu_w ; eu_analysis_results.idx_eu.idx_start_eu(end)];
    idx_end_last_eu_w    = [idx_end_last_eu_w ; eu_analysis_results.idx_eu.idx_end_eu(end)];

    % Duration of the first and last Euglycemia interval
    TIR_first = [TIR_first ; eu_analysis_results.euglycemic_int.duration_eu_int(1)];
    TIR_last  = [TIR_last  ; eu_analysis_results.euglycemic_int.duration_eu_int(end)];
    
    % 'Total time in Eu range' (TIR) between the first and last Euglycemia 
    % interval (limits included)
    TIR_eu = [TIR_eu ; sum(eu_analysis_results.euglycemic_int.duration_eu_int)];
    
    % 'Total time outside Eu range' (TOR) between the first and last
    % Euglycemia interval(limits included). Is the sum of the duration of
    % all events.
    idx_end_first_eu  = eu_analysis_results.idx_eu.idx_end_eu(1);
    idx_start_last_eu = eu_analysis_results.idx_eu.idx_end_eu(end);
    
    % Index of the first event after the first Euglycemia interval
    for i = 1:1:length(events_analysis_results.idx_events.idx_start_event)
        if events_analysis_results.idx_events.idx_start_event(i)>idx_end_first_eu
            n_first_event = i;
            break
        end
    end

    % Index of the last event before the last Euglycemia interval
    for i = 1:1:length(events_analysis_results.idx_events.idx_start_event)
        if events_analysis_results.idx_events.idx_start_event(i)>idx_start_last_eu
            n_last_event = i-1;
            break
        else
            n_last_event = i;
        end
    end
    
    TOR_eu = [TOR_eu ; sum(events_analysis_results.glycemic_events.duration_event(n_first_event:n_last_event))];
    
    % TAR = Time Above Range (Hyper)
    idx_TAR = find(events_analysis_results.glycemic_events.f_critic_event == 0 & (events_analysis_results.glycemic_events.type_event == 'Mild Hyper' | events_analysis_results.glycemic_events.type_event == 'Severe Hyper'));  
    TAR     = [TAR ; sum(events_analysis_results.glycemic_events.duration_event(idx_TAR))];

    % TBR = Time Below Range (Hypo)
    idx_TBR = find(events_analysis_results.glycemic_events.f_critic_event == 0 & (events_analysis_results.glycemic_events.type_event == 'Mild Hypo' | events_analysis_results.glycemic_events.type_event == 'Severe Hypo'));  
    TBR     = [TBR ; sum(events_analysis_results.glycemic_events.duration_event(idx_TBR))];

    % TIR and TOR as a percentage
    tmp_duration_cgm = minutes(PD.PD.Time(eu_analysis_results.idx_eu.idx_end_eu(end))- PD.PD.Time(eu_analysis_results.idx_eu.idx_start_eu(1)));
    duration_cgm     = [duration_cgm ; tmp_duration_cgm];
    TIR_eu_perc      = [TIR_eu_perc ; round((TIR_eu(end)/tmp_duration_cgm)*100)];
    TOR_eu_perc      = [TOR_eu_perc ; round((TOR_eu(end)/tmp_duration_cgm)*100)];

    
    % Number of critical events between the start of the cgm acquisition 
    % and the end of the last Euglycemia interval
    tmp_critic_ev    = length(find(events_analysis_results.glycemic_events.f_critic_event(1:n_last_event)==1));
    n_critical_event = [n_critical_event ; tmp_critic_ev];
    
    % Number of critical gap inside Euglycemia interval
    tmp_critic_eu = eu_analysis_results.eu_pax_report.n_critical_gap_removed;
    n_critical_eu = [n_critical_eu ; tmp_critic_eu];
    
    % Total duration of critical gap (defined as number of missing samples 
    % times 5') between the start of the acquisition and the end of the
    % last Euglycemia interval.
    
    info_critic_gap     = eu_analysis_results.idx_critical_gap;    
    idx_end_analysis    = eu_analysis_results.idx_eu.idx_end_eu(end);
    
    idx_critic_gap      = find(eu_analysis_results.idx_critical_gap.idx_start_critic_gap<idx_end_analysis);
    tmp_n_critic_gap    = length(idx_critic_gap);
    tmp_duration_cr_gap = ((info_critic_gap.idx_end_critic_gap(idx_critic_gap)-info_critic_gap.idx_start_critic_gap(idx_critic_gap))+1)*5;    
    tmp_tot_dur_cr_gap  = sum(tmp_duration_cr_gap);
    
    % Gap duration percentage
    dur_gap      = [dur_gap ; tmp_tot_dur_cr_gap];
    dur_gap_perc = [dur_gap_perc ; round((tmp_tot_dur_cr_gap/tmp_duration_cgm)*100)];
    
    % Median, Max and Min critic gap duration
    if length(tmp_duration_cr_gap)>0
        tmp_median_dur_cr_g = median(tmp_duration_cr_gap);
        tmp_max = max(tmp_duration_cr_gap);
        tmp_min = min(tmp_duration_cr_gap);
    else
        tmp_median_dur_cr_g = 0;
        tmp_max = 0;
        tmp_min = 0;
    end
        
    duration_cr_gap = [duration_cr_gap ; tmp_tot_dur_cr_gap];
    median_t_cr_gap = [median_t_cr_gap ; tmp_median_dur_cr_g];
    max_t_cr_gap    = [max_t_cr_gap ; tmp_max];
    min_t_cr_gap    = [min_t_cr_gap ; tmp_min];
    n_critic_gap    = [n_critic_gap ; tmp_n_critic_gap];
    
    % Flag the last Euglycemia interval if it end right before the start of
    % a critic gap
    if length(find(info_critic_gap.idx_start_critic_gap == (idx_end_analysis+1)))==1
        f_last_eu_gap = [f_last_eu_gap ; 1];
    else
        f_last_eu_gap = [f_last_eu_gap ; 0];
    end
    
    % Flag the first Euglycemia interval if it begin right before the start
    % of a critic gap
    if length(find(info_critic_gap.idx_end_critic_gap == (eu_analysis_results.idx_eu.idx_start_eu(1)-1)))==1
        f_first_eu_gap = [f_first_eu_gap ; 1];
    else
        f_first_eu_gap = [f_first_eu_gap ; 0];
    end

end

close all
%% RESULTS TABLE

pax_analysis = table(pax_idx_vector,time_before_first_eu,TIR_first,TIR_last,duration_cgm,TIR_eu,TIR_eu_perc,TOR_eu,TOR_eu_perc,TAR,TBR,dur_gap,dur_gap_perc);
gap_analysis = table(pax_idx_vector,n_critic_gap,f_first_eu_gap,f_last_eu_gap,duration_cr_gap,median_t_cr_gap,min_t_cr_gap,max_t_cr_gap);

sincro_results = table(pax_idx_vector,idx_start_first_eu_w,idx_end_first_eu_w,idx_start_last_eu_w,idx_end_last_eu_w);

tot_CGM_duration_results = table(pax_idx_vector,dur_CGM);

%% MAX and MIN duration of CGM signal (Congress)

min_cgm     = min(tot_CGM_duration_results.dur_CGM);
idx_min_cgm = pax_idx_vector(find(tot_CGM_duration_results.dur_CGM == min_cgm));
min_cgm_tab = table(idx_min_cgm,min_cgm);
disp(min_cgm_tab);

max_cgm     = max(tot_CGM_duration_results.dur_CGM);
idx_max_cgm = pax_idx_vector(find(tot_CGM_duration_results.dur_CGM == max_cgm));
max_cgm_tab = table(idx_max_cgm,max_cgm);
disp(max_cgm_tab);

%% PATIENTS SELECTION

%% CRITICAL GAP SELECTION

% Patients with no critical gap between the start of the CGM acqusition and
% the end of the last Euglycemia Interval. If a critical gap starts 
% immediately after the end of the last Euglycemia Interval, the patient is
% removed as well.

idx_no_gap   = [];

for i=1:1:n_pax
    if gap_analysis.n_critic_gap(i)==0 & gap_analysis.f_last_eu_gap(i)==0
        idx_no_gap = [idx_no_gap ; i];
    end
end

disp('==============================================================================')
disp('PATIENTS WITH NO CRITICAL GAP')
disp(num2str(pax_idx_vector(idx_no_gap)'))

pax_analysis_no_gap = pax_analysis(idx_no_gap,:);

%% TOTAL TIME OUTSIDE EUGLYCEMIA SELECTION

th_25_perc_TOR = quantile(pax_analysis_no_gap.TOR_eu_perc(:),0.25)

figure('Position', [100, 100, 700, 400])
boxplot(pax_analysis_no_gap.TOR_eu_perc(:),'Widths',0.5)
title(['Time Outside Euglycemia Range [%] , 25-th percentile = ',num2str(th_25_perc_TOR),'%'])
title(['Threshold for the third criterion of patient selection (Q1 of TOR%): ',num2str(th_25_perc_TOR),'%'])
ylabel('Time [%]')
xticklabels('Overall Time Outside euglycemia Range [TOR%]')
% print('FIG_23_BP_EU_SOGLIA','-djpeg','-r600')

pax_analysis_after_TOR = pax_analysis_no_gap(pax_analysis_no_gap.TOR_eu_perc(:)>=th_25_perc_TOR,:);

pax_selected = pax_analysis_after_TOR.pax_idx_vector;

disp('==============================================================================')
disp('PATIENTS SELECTED')
disp(num2str(pax_selected'))

%% TOTAL TIME INSIDE EUGLYCEMIA SELECTION (NOT CONSIDERED FOR MY THESIS)

% th_75_perc_TIR = quantile(pax_analysis_after_TOR.TIR_eu_perc(:),0.75)
% 
% figure()
% boxplot(pax_analysis_after_TOR.TIR_eu_perc(:))
% title(['Time Inside Euglycemia Range [%] , 75-th percentile = ',num2str(th_75_perc_TIR),'%'])
% 
% pax_analysis_after_TIR = pax_analysis_after_TOR(pax_analysis_after_TOR.TIR_eu_perc(:)<=th_75_perc_TIR,:);
% 
% pax_tenuti = pax_analysis_after_TIR.pax_idx_vector

%% SAVE OF THE RESULTS

 save pax_analysis_results_no_corr.mat pax_analysis
 save pax_event_report_no_corr.mat pax_event_report
close all