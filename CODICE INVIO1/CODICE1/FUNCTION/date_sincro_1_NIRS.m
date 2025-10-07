function [date_start_sincro,date_end_sincro] = date_sincro_1_NIRS(date_CGM,date_start,date_end,f_disp)
% Input:  Date of CGM acquisition, starting and ending date for NIRS 
%         dataset. f_disp = 1 display results in Command Window
% Output: Interval date for CGM-NIRS synchronization

% Extraction of the dates of CGM start and end
date_CGM_start = date_CGM(1);
date_CGM_end = date_CGM(end);

% Setting variables
tmp = 0;
tmp_s = 0;
NIRS_gap_s = 0; % Flag start
NIRS_gap_e = 0; % Flag end

% Search for the date of the start of synchronization

% If synchronization is not possible the flag f_s is set to 1
f_s = 0;

% CGM starts before the start of the first NIRS dataset
if date_CGM_start <= date_start
    tmp_s = find(date_CGM >= date_start,1,'first');
    date_start_sincro = date_CGM(tmp_s);
    NIRS_gap_s = 1;
end

% CGM starts inside the first NIRS dataset, presence of NIRS gap 
if (date_CGM_start >= date_start & date_CGM_start<= date_end)
    date_start_sincro = date_CGM_start;
    NIRS_gap_s = 1;
end

% CGM starts after the end of the second (and last) dataset
if date_CGM_start > date_end
    disp('SYNCHRONIZATION NOT POSSIBLE')
    f_s = 1;
    date_start_sincro = NaN;
end

% Search for the date of the end of synchronization

% CGM ends inside the first NIRS dataset
if (date_CGM_end >= date_start & date_CGM_end <= date_end)  
    date_end_sincro = date_CGM_end;
end

% CGM ends date is greater than the end of the NIRS acquisition
if date_CGM_end > date_end
    tmp = find(date_CGM<=date_end,1,'last');
    date_end_sincro = date_CGM(tmp);
    NIRS_gap_e = 1;
end

% CGM ends date is before the starts of the NIRS acquisition
if date_CGM_end < date_start
    disp('SYNCHRONIZATION NOT POSSIBLE')
    f_s = 1;
    date_end_sincro = NaN;
end

% Display of the synchronization date
if f_disp == 1
    disp('==============================================================================')
    disp('CGM - NIRS SYNCHRONIZATION INTERVAL DATE')
    disp(' ')
    if f_s == 0
        disp(['Start : ',datestr(date_start_sincro)])
        disp(['End   : ',datestr(date_end_sincro)])
    else
        disp('SYNCHRONIZATION NOT POSSIBLE')
        date_start_sincro = NaN;
        date_end_sincro = NaN;
    end
end
