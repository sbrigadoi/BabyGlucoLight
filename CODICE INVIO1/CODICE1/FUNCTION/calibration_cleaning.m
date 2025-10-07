function [PD_clean] = calibration_cleaning(PD_raw)
% Cleaning of the dataset from calibration samples
%
% [PD_clean] = calibration_cleaning(PD_raw)
%
% Input: PD_raw : table which contains at least one column called 
%                 'Calibration' with the flag of calibration in categorical
%                 format and one column called 'Time' with the date of
%                 acquisition in 'datetime' format.
% 
% Output: PD_clean : the same imput table whitout the row which corrisponds
%                    to the calibration samples.

% Extraction of the column with the calibration flag
calibration = PD_raw.Calibration;
% Extracting the name of the calibration flag
name_calibration_flag = categories(calibration);

% Calibration removal
to_remove = [];

for i = 1:1:length(calibration)
    if calibration(i) == name_calibration_flag
       to_remove = [to_remove i];
    end
end

PD_raw(to_remove,:) = [];

% Output of the function
PD_clean = PD_raw;

% Extraction of the start and finish date of the CGM acquisition
date_CGM = PD_clean.Time;
date_CGM_start = date_CGM(1);
date_CGM_end = date_CGM(end);

% Display of the results
% disp('==============================================================================')
disp(' ')
disp('CGM INTERVAL DATE')
disp(['Start : ',datestr(date_CGM_start,0)])
disp(['End   : ',datestr(date_CGM_end,0)])
disp(' ')
disp(['CGM CALIBRATION CORRECTIONS: ',num2str(length(to_remove))])