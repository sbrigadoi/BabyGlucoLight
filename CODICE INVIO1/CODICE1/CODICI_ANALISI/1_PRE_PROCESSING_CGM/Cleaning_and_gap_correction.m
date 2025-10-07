% CLEANING AND CORRECTION OF THE RAW DATASET

close all
clear all
clc

%% SETTING OF FORLDERS PATH 
base_path     = pwd;
data_path     = fullfile(base_path,'DATI','CGM_ESTRATTO');
function_path = fullfile(base_path,'CODICE','DEFINITIVO/PRE_PROCESSING/');

addpath(genpath(pwd))
addpath(data_path)
addpath(function_path)

%% DATASET CLEANING AND GAP CORRECTION

% Initialization of support variables 
file_name        = 'PD%d.mat';
pax_idx_vector   = [1:1:17 19:1:20 22:1:26 28 33:1:45 47:1:55]';
n_pax            = length(pax_idx_vector);
PD_idx           = 0;
n_gap            = 0;
duration_tot_gap = 0;

% Display of the preliminar information
disp('CLEANING AND GAP CORRECTION')
disp(' ')
disp(['Number of patients: ',num2str(n_pax)])
disp(' ')
disp('START')
disp(' ')

% Cleaning and gap correction for each patient inside the dataset
for i = 1:1:n_pax
    
    pax_idx = pax_idx_vector(i);
    disp('==============================================================================')
    disp(['PATIENT ',num2str(pax_idx)])
    disp('==============================================================================')
    
    % Load of the current patient raw data
    to_load    = sprintf(file_name,pax_idx);
    
    PD_struct  = load(to_load);
    table_name = string(fieldnames(PD_struct));
    PD_raw     = getfield(PD_struct,table_name);
    
    % When the 'Calibration' flag column is present, it is necessary to
    % remove the calibration samples. 
    variable_names = PD_raw.Properties.VariableNames;
    
    if ismember('Calibration',variable_names)
        % Calibration correction
        PD_clean = calibration_cleaning(PD_raw);
        cgm_value = PD_clean.mgdL;
        cgm_date  = PD_clean.Time;
    else
        % No calibration correction is needed
        PD_clean = PD_raw;
        cgm_value = PD_clean.mgdL;
        cgm_date  = PD_clean.Time;
        % Manual extraction of the start and finish date of the acquisition
        date_CGM_start = cgm_date(1);
        date_CGM_end   = cgm_date(end);

        % Manual display of the results
        disp(' ')
        disp('CGM INTERVAL DATE')
        disp(['Start : ',datestr(date_CGM_start,0)])
        disp(['End   : ',datestr(date_CGM_end,0)])
        disp(' ')
        disp('CGM CALIBRATION CORRECTIONS: NOT REQUIRED')
    end
    
    % Temporal gap correction
    PD_full   = gap_correction(cgm_value,cgm_date);
    
    % Current patient identification number and total number of gap
    PD_idx(i,1) = pax_idx;
    n_gap(i,1)  = PD_full.n_gap;
    
    % Total gap duration for the current patient
    if PD_full.n_gap > 0
        duration_tot_gap(i,1) = sum((PD_full.infogap.n_samples)*5);
    else
        duration_tot_gap(i,1) = 0;
    end
    
    % Save of the final result (!!PAY ATTENTION WHEN REMOVE THE COMMENT!!)
    
%     PD = PD_full;
%     save_PD_path_name   = fullfile(base_path,'DATI','CGM_CLEAN','DEFINITIVI',to_load);
%     save (save_PD_path_name,'PD')
            
end

% Display of the last information
disp(' ')
disp(' ')
disp('END')
disp('==============================================================================')
disp('==============================================================================')
disp(' ')

% Creation of a table with the patient identification number, the number of
% gap and the total gap duration.
PD_gap_report = table(PD_idx,n_gap,duration_tot_gap);

% Display of the table
disp('PD GAP REPORT TABLE')
disp(' ')
disp(PD_gap_report)

% Save of the gap report table (!!PAY ATTENTION WHEN REMOVE THE COMMENT!!)
% save_gap_report_path   = fullfile(base_path,'DATI','CGM_CLEAN','DEFINITIVI','PD_GAP_REPORT_TABLE.mat');
% save(save_gap_report_path,'PD_gap_report')