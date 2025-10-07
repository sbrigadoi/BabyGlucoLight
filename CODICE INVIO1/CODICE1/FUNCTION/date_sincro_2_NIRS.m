function [date_start_sincro,date_end_sincro] = date_sincro_2_NIRS(date_CGM,date_start_first,date_end_first,date_start_second,date_end_second,f_disp)
% Input:  Date of CGM acquisition, starting and ending date for both NIRS
%         dataset. f_disp = 1 display results in Command Window
% Output: Interval date for CGM-NIRS synchronization

%%
% Extraction of the dates of CGM start and end
date_CGM_start = date_CGM(1);
date_CGM_end = date_CGM(end);

% Setting variables
tmp = 0;
tmp_s = 0;
NIRS_gap_s = 0; % Flag start
NIRS_gap_e = 0; % Flag end

% Search for the date of the start of synchronization

% CGM starts before the start of the first NIRS dataset
if date_CGM_start <= date_start_first & date_CGM_end>date_end_first
    tmp_s = find(date_CGM >= date_start_first,1,'first');
    date_start_sincro = date_CGM(tmp_s);
    NIRS_gap_s = 1;
end

if date_CGM_start <= date_start_first & date_CGM_end<=date_end_first
    tmp_s = find(date_CGM >= date_start_first,1,'first');
    date_start_sincro = date_CGM(tmp_s);
end

% CGM starts inside the first NIRS dataset, presence of NIRS gap 
if (date_CGM_start >= date_start_first & date_CGM_end> date_end_first)
    date_start_sincro = date_CGM_start;
    NIRS_gap_s = 1;
end

if (date_CGM_start >= date_start_first & date_CGM_end <= date_end_first)
    date_start_sincro = date_CGM_start;
end

% CGM starts inside the second NIRS dataset
if date_CGM_start >= date_start_second & date_CGM_start <= date_end_second
        date_start_sincro = date_CGM_start;
end
    
% CGM starts inside the NIRS acquisition break
if date_CGM_start > date_end_first & date_CGM_start < date_start_first
    tmp_s = find(date_CGM >= date_start_second,1,'first');
    date_start_sincro = date_CGM(tmp_s);
end

% CGM starts after the end of the second (and last) dataset
if date_CGM_start > date_end_second
    disp('SYNCHRONIZATION NOT POSSIBLE')
end
%%
% Search for the date of the end of synchronization

% CGM ends inside the first NIRS dataset
if (date_CGM_end >= date_start_first & date_CGM_end <= date_end_first)  
    date_end_sincro = date_CGM_end;
end

% CGM ends inside the second NIRS dataset
if date_CGM_end <= date_end_second & date_CGM_start >= date_start_second
    date_end_sincro = date_CGM_end;
end

% CGM ends inside the second NIRS dataset
if date_CGM_end <= date_end_second & date_CGM_start < date_start_second
    date_end_sincro = date_CGM_end;
    NIRS_gap_e = 1;
end

% CGM ends inside the NIRS acquisition break
if date_CGM_end > date_end_first & date_CGM_end < date_start_second
    tmp = find(date_CGM<=date_end_first,1,'last');
    date_end_sincro = date_CGM(tmp);
end

% CGM ends date is greater than the end of the second NIRS acquisition (PD15
% case)
if date_CGM_end > date_end_second
    tmp = find(date_CGM<=date_end_second,1,'last');
    date_end_sincro = date_CGM(tmp);
    NIRS_gap_e = 1;
end

% CGM ends date is before the starts of the NIRS acquisition
if date_CGM_end < date_start_first
    disp('SYNCHRONIZATION NOT POSSIBLE')
end

% If NIRS gap is present, the synchronization interval is divided

% DA CONTROLLARE

% if NIRS_gap_s == 1 & NIRS_gap_e == 1
%     date_start_sincro = [date_start_sincro , date_start_second];
%     date_end_sincro = [date_end_first , date_end_sincro];
% end

%%
% Display of the synchronization date
if f_disp == 1 
    disp('==============================================================================')
    disp('CGM - NIRS SYNCHRONIZATION INTERVAL DATE')
    disp(' ')
    disp('Start : ')
    disp(datestr(date_start_sincro))
    disp('End   : ')
    disp(datestr(date_end_sincro))
    if NIRS_gap_s == 1 & NIRS_gap_e == 1
        disp(['NIRS gap: Present'])
    else
        disp(['NIRS gap: Not Present'])
    end
end