% PATIENT SELECTION

close all
clear all
clc

%% SETTING OF FORLDERS PATH 

base_path     = pwd;
code_path     = fullfile(base_path,'CODICE','DEFINITIVO','CGM_SELEZIONE_PAX_NIRS');
data_path     = fullfile(base_path,'DATI','CGM_INTERP');
function_path = fullfile(base_path,'CODICE','DEFINITIVO','CGM_SELEZIONE_PAX_NIRS');
cgm_results_path = fullfile(base_path,'RISULTATI DEFINITIVI/');

addpath(genpath(pwd))
addpath(code_path)
addpath(data_path)
addpath(function_path)
addpath(cgm_results_path)


%% SYNCHRONIZATION INFO

% DATE PRESE DA PDF INFO SYNCHRO, SONO LE DATE DI INIZIO E FINE
% SINCRONIZZAZIONE NIRS-CGM

% PD4_s = datetime('12-Apr-2020 18:20:00');
% PD4_e = datetime('16-Apr-2020 12:50:00');
% 
% PD5_s = datetime('17-Apr-2020 00:00:00');
% PD5_e = datetime('20-Apr-2020 09:00:00');
% 
% PD8_s = datetime('02-May-2020 19:34:00');
% PD8_e = datetime('06-May-2020 16:29:00');
% 
% PD9_s = datetime('09-May-2020 01:24:00');
% PD9_e = datetime('12-May-2020 18:49:00');
% 
% PD11_s = datetime('28-May-2020 21:18:00');
% PD11_e = datetime('29-May-2020 17:23:00');
% 
% PD14_s = datetime('23-Jun-2020 20:43:00');
% PD14_e = datetime('27-Jun-2020 00:52:00');
% 
% PD15_s = datetime('16-Jul-2020 14:05:00');
% PD15_e = datetime('19-Jul-2020 13:10:00');
% 
% PD19_s = datetime('19-Sep-2020 17:08:00');
% PD19_e = datetime('21-Sep-2020 21:38:00');
% 
% PD47_s = datetime('10-Sep-2022 20:02:00');
% PD47_e = datetime('13-Sep-2022 13:32:00');
% 
% date_cgm_1 = table();
% idx_PD   = [4;5;8;9;11;14;15;19;47];
% 
% date_cgm_1.idx_PD = idx_PD;
% date_cgm_1.start  = [PD4_s;PD5_s;PD8_s;PD9_s;PD11_s;PD14_s;PD15_s;PD19_s;PD47_s];
% date_cgm_1.end    = [PD4_e;PD5_e;PD8_e;PD9_e;PD11_e;PD14_e;PD15_e;PD19_e;PD47_e];
% date_cgm_1.dur    = minutes(date_cgm_1.end - date_cgm_1.start);

% PD_s = datetime();
% PD_e = datetime();

%% GLYCEMIC METRICS ANALYSIS

file_name = 'PD%d.mat';
%idx_PD    = [4;5;8;9;11;14;15;19;47];

idx_PD    = [3;4;5;8;9;10;11;14;15;19;25;47];
n_pax     = length(idx_PD);


% Indexes of first and last Euglycemia window extracted from the patient
% selection phase.
load idx_eu_results.mat

% Load of patient analysis and event report already corrected
%load pax_analysis_results_for_cca.mat
%load pax_event_report_for_cca.mat

load pax_analysis_results_with_corr.mat
load pax_event_report_with_corr.mat

%get rid of PD59
pax_analysis = pax_analysis(1:47,:);
pax_event_report = pax_event_report(1:47,:);

%just use the 12PD we will want to use
pax_analysis = pax_analysis(idx_PD,:);
pax_event_report = pax_event_report(idx_PD,:); 

%%
% Initialization of usefull variables
Ts = 5;  %[Min]
date_cgm = table();
date_cgm.idx_PD = idx_PD;
cgm_metrics = table();
% cgm_metrics.PD = idx_PD;

% Display of the preliminar information
disp('EXTRACTION OF CGM METRICS')
disp(' ')
disp(['PD index = ',num2str(idx_PD')])
disp(' ')
disp(['Number of patients: ',num2str(n_pax)])
disp(' ')
disp('START')
disp(' ')

% Extracion of CGM metrics for each selected patients
dur_hour = [];
warning('off')
for i = 1:1:n_pax
    
    idx_curr_PD = idx_PD(i);
    disp('==============================================================================')
    disp(['PATIENT ',num2str(idx_curr_PD)])
    
    % Load of the current patient data
    to_load     = sprintf(file_name,idx_curr_PD);
    
    PD_struct   = load(to_load);
    table_name  = string(fieldnames(PD_struct));
    PD          = getfield(PD_struct,table_name);
    
    % EXTRACTION OF CGM DATA 

    % Total Time of acquisition 
    % Total time is defined as the time between the start of the first
    % euglycemia window and the end of the last euglycemia window.
    
    % The start and end date is extracted from the 'idx_eu_results.mat'
    % obtained from the PD selection
    
    % Pay attention that for PD5 it was used as last euglycemia window the
    % second-last since the last event was flagged as critic (only critic event
    % for that patient)
    idx_tab   = find(idx_eu_int.pax_idx_vector == idx_curr_PD);
    tmp_start = PD.PD.Time(idx_eu_int.idx_start_first_eu_w(idx_tab));
    tmp_end   = PD.PD.Time(idx_eu_int.idx_end_last_eu_w(idx_tab));

    date_cgm.start(i) = tmp_start;
    date_cgm.end(i)   = tmp_end;

    % PD 5 is analyzed separately
    if idx_curr_PD == 5
        date_cgm.end(i) = datetime('20-Apr-2020 09:00:00');
        idx_eu_int.idx_start_last_eu_w(idx_tab) = find(PD.PD.Time == datetime('19-Apr-2020 22:05:00'));
        idx_eu_int.idx_end_last_eu_w(idx_tab) = find(PD.PD.Time == date_cgm.end(i));
    end

    % Duration of the CGM acquisition in minutes
    date_cgm.dur(i)   = minutes(date_cgm.end(i)-date_cgm.start(i));

    % Start and end index 
    tmp_idx_s = find(PD.PD.Time == PD.PD.Time(idx_eu_int.idx_start_first_eu_w(idx_tab)));
    tmp_idx_e = find(PD.PD.Time == PD.PD.Time(idx_eu_int.idx_end_last_eu_w(idx_tab)));

    % MEAN, MAD, COEFFICIENT OF VARIATION AND STANDARD DEVIATION 
    cgm_metrics.mean(i,1)   = mean(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e));
    cgm_metrics.sd(i,1)     = std(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e));
    cgm_metrics.median(i,1) = median(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e));
    cgm_metrics.mad(i,1)    = mad(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e),1);
    cgm_metrics.cv(i,1)     = std(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e))/mean(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e));

    % First and third quartile of the cgm signal
    cgm_metrics.q1(i,1) = quantile(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e),0.25);
    cgm_metrics.q3(i,1) = quantile(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e),0.75);
    % Interquartile range
    cgm_metrics.iqr(i,1) = cgm_metrics.q3(i)-cgm_metrics.q1(i);
  
    % RATE OF CHANGE (ROC) AND AVERAGE ABSOLUTE RATE OF CHANGE (AARC)
    % ROC = [a(t2)-a(t1)]/(T2-t2)
    pos = 0;
    for j = tmp_idx_s:1:(tmp_idx_e-1)
        pos = pos+1;
        tmp_ROC(pos) = (PD.PD.("mg/dL")(j+1)-PD.PD.("mg/dL")(j))/Ts;
    end
    % AARC = (sum|ROC|)/x    --> x = number of osservation
    cgm_metrics.AARC(i,1) = mean(abs(tmp_ROC));

    % SD RATE OF CHANGE (SDRC)
    % Standard deviation of rate of change values ROC
    cgm_metrics.SDRC(i,1) = std(tmp_ROC);
    

    % CONGAn (Continuous Overall Net Glycemic Action over n*1h)
    % Standard deviation between the current value and the one observed n
    % hours before
    
    % Minutes and number of samples for each n-hours interval of CONGAn
    n_1  = 60;    %[Min]
    n_2  = 60*2;  %[Min]
    n_4  = 60*4;  %[Min]
    n_24 = 60*24; %[Min]
    n_sample_1  = n_1/Ts;
    n_sample_2  = n_2/Ts;
    n_sample_4  = n_4/Ts;
    n_sample_24 = n_24/Ts;
    
    % Number of samples between the start and the end of the GCM interval
    % under analysis
    n_dur = length(PD.PD.Time(tmp_idx_s:tmp_idx_e));

    % CONGA1
    % Standard deviation computed over the vector with differences of
    % glucose from the current value and the value 1 hour before (starting
    % from one hour after the start of CGM interval under analysis)
    pos = 0;
    for k = (tmp_idx_s+(n_sample_1)):1:n_dur
        pos = pos+1;
        tmp_diff_CONGA1(pos) = PD.PD.("mg/dL")(k)-PD.PD.("mg/dL")(k-n_sample_1); 
    end
    cgm_metrics.CONGA1(i,1) = std(tmp_diff_CONGA1);

    % CONGA2
    % Standard deviation computed over the vector with differences of
    % glucose from the current value and the value 2 hour before (starting
    % from two hour after the start of CGM interval under analysis)
    pos = 0;
    for k = (tmp_idx_s+(n_sample_2)):1:n_dur
        pos = pos+1;
        tmp_diff_CONGA2(pos) = PD.PD.("mg/dL")(k)-PD.PD.("mg/dL")(k-n_sample_2); 
    end
    cgm_metrics.CONGA2(i,1) = std(tmp_diff_CONGA2);

    % CONGA4
    % Standard deviation computed over the vector with differences of
    % glucose from the current value and the value 4 hour before (starting
    % from four hour after the start of CGM interval under analysis)
    pos = 0;
    for k = (tmp_idx_s+(n_sample_4)):1:n_dur
        pos = pos+1;
        tmp_diff_CONGA4(pos) = PD.PD.("mg/dL")(k)-PD.PD.("mg/dL")(k-n_sample_4); 
    end
    cgm_metrics.CONGA4(i,1) = std(tmp_diff_CONGA4);

    % MODD = Mean Of Daily Difference (interday glycemic variability)             CALCOLO SOLO SU DUE GIORNI O DI PIU'?
    % This parameter is calculated as the mean of the absolute differences
    % between glucose values at the same time on two consecutive days.
    % Attention PD11 doesn't reach the 24 hours of acquisition --> 0
    pos = 0;
    n_two_day = (tmp_idx_s+(n_sample_24*2));
    for k = (tmp_idx_s+(n_sample_24)):1:n_two_day
        pos = pos+1;
        if k > length(PD.PD.("mg/dL"))
            % tmp_diff_MODD(pos) = Inf;
            tmp_diff_MODD(pos) = 0;

            continue
        end
        tmp_diff_MODD(pos) =  abs(PD.PD.("mg/dL")(k)-PD.PD.("mg/dL")(k-n_sample_24));
    end
    cgm_metrics.MODD(i,1) = mean(tmp_diff_MODD);

    % GMI = Glucose Management Indicator (NOT USED)
    % cgm_metrics.GMI(i,1) = 3.31 + (0.02392*cgm_metrics.mean(i));

    % e_alc = Estimate HbA1c (emoglobina glicosilata, biomarker for the 
    % risk of developing chronic complications of diabetes). Now it is 
    % used the GMI parameter. (NOT USED)

    % cgm_metrics.e_alc(i) = (46.7 + cgm_metrics.mean(i))/28.7;

    % AUC = Area Under the Curve (trapz method)                            
    tmp_AUC = trapz(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e));
    % AUC_n = Normalized Area Under the Curve (over total acquisition time)
    cgm_metrics.AUC_n(i) = tmp_AUC/(date_cgm.dur(i)); 


    % Min and Max CGM value
    cgm_metrics.min(i) = min(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e));
    cgm_metrics.max(i) = max(PD.PD.("mg/dL")(tmp_idx_s:tmp_idx_e));

end
warning('on')

% Add TIR% TOR% TBR% TAR% to the cgm_metrics table
cgm_metrics.TIR_perc = pax_analysis.TIR_eu_perc;
cgm_metrics.TOR_perc = pax_analysis.TOR_eu_perc;
cgm_metrics.TBR_perc = pax_analysis.TBR_perc;
cgm_metrics.TAR_perc = pax_analysis.TAR_perc;

% Add number of normalized hypo/hyper glycemic event to the cgm_metrics 
% table (over total number of events)
for i = 1:1:length(idx_PD)
    tmp_n_hypo_n  = (pax_event_report.n_mild_hypo(i)+pax_event_report.n_severe_hypo(i))/pax_event_report.n_events(i);
    tmp_n_hyper_n = (pax_event_report.n_mild_hyper(i)+pax_event_report.n_severe_hyper(i))/pax_event_report.n_events(i);
    cgm_metrics.hypo_n(i) = tmp_n_hypo_n;
    cgm_metrics.hyper_n(i) = tmp_n_hyper_n;
end

%%
% to_save = fullfile(base_path,'RISULTATI DEFINITIVI/CGM_METRICS/prova_cgm_metrics.mat');
% save(to_save, 'cgm_metrics')
% disp('CGM Metrics saved')

disp('CGM Metrics done')