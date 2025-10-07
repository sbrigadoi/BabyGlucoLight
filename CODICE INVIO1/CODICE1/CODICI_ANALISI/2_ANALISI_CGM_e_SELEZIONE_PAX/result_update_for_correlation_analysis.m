% MANUAL CORRECTION OF CGM ANALYSIS RESULTS

% -) Remove of the critic events on PD5 with update of the relative metrics
% -) Remove from all metrics the last event if it is after the last 
%    euglycemia interval
% -) Correction of the duration of CGM defined as the sum of the time spent
%    outside and inside the euglycemia interval
% -) Correction of the TOR% e TIR% (now the sum is 100%)
% -) Computation of TAR% e TBR%

clear all
close all
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

%%
% Load of automated glycemic analysis results
%load pax_analysis_results.mat
%load pax_event_report.mat

%these have had no correction 17 03 25
load('pax_analysis_results_no_corr.mat')
load('pax_event_report_no_corr.mat')

%idx_PD    = [4;5;8;9;11;14;15;19;47];
%new 01 03 25 so have to +1 to5 and 9, +2 to 14 15 19
%no 59, as we can't use 59
idx_PD    = [3;4;5;8;9;10;11;14;15;19;25;47];


% Manual correction to remove the event outside the last euglycemia

pd_n=5;
% PD 5
pax_event_report.n_mild_hypo(pax_event_report.idx_PD == pd_n) = pax_event_report.n_mild_hypo(pax_event_report.idx_PD == pd_n)-1;
pax_event_report.n_events(pax_event_report.idx_PD == pd_n) = pax_event_report.n_events(pax_event_report.idx_PD == pd_n)-1;
pax_event_report.n_events_critic(pax_event_report.idx_PD == pd_n) = 0;

pax_analysis.TOR_eu(pax_analysis.pax_idx_vector == pd_n) = pax_analysis.TOR_eu(pax_analysis.pax_idx_vector == pd_n)-75;
pax_analysis.TIR_eu(pax_analysis.pax_idx_vector == pd_n) = pax_analysis.TIR_eu(pax_analysis.pax_idx_vector == pd_n)-135;
pax_analysis.TIR_last(pax_analysis.pax_idx_vector == pd_n) = 655;
pax_analysis.duration_cgm(pax_analysis.pax_idx_vector == pd_n) = 4860;
pax_analysis.dur_gap(pax_analysis.pax_idx_vector == pd_n) = 0;
pax_analysis.dur_gap_perc(pax_analysis.pax_idx_vector == pd_n) = 0;

pd_n = 9;
% PD 9
pax_event_report.n_severe_hyper(pax_event_report.idx_PD == pd_n) = pax_event_report.n_severe_hyper(pax_event_report.idx_PD == pd_n)-1;
pax_event_report.n_events(pax_event_report.idx_PD == pd_n) = pax_event_report.n_events(pax_event_report.idx_PD == pd_n)-1;

pd_n = 14;
% PD 14
pax_event_report.n_severe_hyper(pax_event_report.idx_PD == pd_n) = pax_event_report.n_severe_hyper(pax_event_report.idx_PD == pd_n)-1;
pax_event_report.n_events(pax_event_report.idx_PD == pd_n) = pax_event_report.n_events(pax_event_report.idx_PD == pd_n)-1;

pd_n = 15;
% PD 15
pax_event_report.n_severe_hypo(pax_event_report.idx_PD == pd_n) = pax_event_report.n_severe_hypo(pax_event_report.idx_PD == pd_n)-1;
pax_event_report.n_events(pax_event_report.idx_PD == pd_n) = pax_event_report.n_events(pax_event_report.idx_PD == pd_n)-1;

pd_n = 19;
% PD 19
pax_event_report.n_severe_hyper(pax_event_report.idx_PD == pd_n) = pax_event_report.n_severe_hyper(pax_event_report.idx_PD == pd_n)-1;
pax_event_report.n_events(pax_event_report.idx_PD == pd_n) = pax_event_report.n_events(pax_event_report.idx_PD == pd_n)-1;

for i = 1:1:length(idx_PD)
    pax_analysis.duration_cgm(i) = pax_analysis.TIR_eu(i)+pax_analysis.TOR_eu(i);
    pax_analysis.TIR_eu_perc(i)  = round((pax_analysis.TIR_eu(i)/pax_analysis.duration_cgm(i))*100,2);
    pax_analysis.TOR_eu_perc(i)  = round((pax_analysis.TOR_eu(i)/pax_analysis.duration_cgm(i))*100,2);
end

pax_analysis.time_before_first_eu = [];
pax_analysis.dur_gap = [];
pax_analysis.dur_gap_perc = [];
pax_event_report.f_start_inside_ev = [];
pax_event_report.f_end_inside_ev = [];
pax_event_report.n_events_critic = [];

pax_analysis.TAR_perc = round((pax_analysis.TAR./pax_analysis.duration_cgm)*100,2);
pax_analysis.TBR_perc = round((pax_analysis.TBR./pax_analysis.duration_cgm)*100,2);

%% SAVE OF THE RESULTS

 destination_pax_an  = fullfile(base_path,'RISULTATI DEFINITIVI/CGM_METRICS/pax_analysis_results_for_cca.mat');
 destinantion_ev_rep = fullfile(base_path,'RISULTATI DEFINITIVI/CGM_METRICS/pax_event_report_for_cca.mat');
% 
 save (destination_pax_an,'pax_analysis');
 save (destinantion_ev_rep,'pax_event_report');
 disp('Results saved')

disp('DONE')