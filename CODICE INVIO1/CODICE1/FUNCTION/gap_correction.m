function [PD_full]=gap_correction(cgm_value,cgm_date)
% Temporal Gap identification and correction
%
% [PD_clean]=gap_correction(cgm_value,cgm_date)
%
% Inputs: cgm_value = vector of the acquired CGM values (double)
%         cgm_date  = vector of the acquisition date (datetime)
%
% Output: PD_full   = PD_full.PD : table with the date and the glycemic
%                                  value after the gap correction
%
%                     PD_full.infogap : matrix where in the first column 
%                                       there are the indexes of the gap 
%                                       starting point (as the first with 
%                                       corrections), in the second column
%                                       the indexes of the gap ending point
%                                       (as the last with corrections) in
%                                       the third column the number of 
%                                       missing sample inside the gap.                                  

%% GAP IDENTIFICATION

% Extraction of variables under analysis and initialization of support
% variables

date_start_gap = [];
date_end_gap   = [];
gap_duration   = [];
n_CGM          = length(cgm_date);
Ts             = 5; %[minutes]
n_gap          = 0;

% Number of missing sample inside the temporal gap
n_missing = [];
% Matrix which contains inside the first column the starting index of the
% gap and in the second column the duration of the gap in terms of sampling
% points
idx_gap = [];

% Computation of the difference between contiguous time points
for i = 1:n_CGM
    
    % Dimensional check
    if i+1 > n_CGM
        break
    end
    
    % Temporal difference between two contiguous samples in double format
    diff_CGM = minutes(time(caldiff([cgm_date(i) cgm_date(i+1)])));
    
    % If the difference is over or equal to 10' there is a temporal gap
    if diff_CGM >= 10
        % Increment of the number of gap 
        n_gap = n_gap + 1;
        % The start date of the gap is 5' after the last sample acquired
        date_start_gap = [date_start_gap dateshift(cgm_date(i),'start','minutes',5)];
        % The end date is the start date translated of the difference
        % between two contiguous samples minus 10'.
        date_end_gap = [date_end_gap dateshift(date_start_gap(end),'start','minutes',diff_CGM-5)];
        % Gap duration in sampling points
        n_missing = (floor(diff_CGM/5)-1);
        % Gap duration is defined as the number of missing samples times 
        % multiplied the CGM sampling time 5'(in minutes)
        gap_duration = [gap_duration n_missing*Ts];
        % Save of the starting index and the gap duration in terms of 
        % sampling points
        idx_gap(n_gap,1) = i;
        idx_gap(n_gap,2) = n_missing;
    end
end

% Display of the preliminary inspection of the CGM dataset
% disp('==============================================================================')
disp(' ')
disp('CGM DATA GAP')
disp(['Number of gap     : ',num2str(n_gap)])
disp('Gap date interval :')
for i = 1:length(date_start_gap)
    disp([datestr(date_start_gap(i),0),' - ',datestr(date_end_gap(i),0),' : ',num2str(gap_duration(i)),' minutes']);
end

%% GAP CORRECTION
% If there are at least one gap 
if n_gap > 0
    % In order to fix the temporal gap it is used the mean value between the
    % two samples among the gap.
    
    cgm_value_full = [];
    cgm_date_full  = [];
    
    n_gap = size(idx_gap,1);
    c     = 1;
    
    % For each CGM value, if the current index is equal to one of the gap start
    % index, in the new full vector is copied the value of the current index 
    % and then the number of missing sample inside the gap times the mean value
    % between the to samples among the gap. Otherwise it is copied the same
    % value. The same approach it was used to create the vector with the date
    % corrisponding to the full vector of CGM samples

    for i = 1:length(cgm_value)

        if i == idx_gap(c,1)
            cgm_value_full = [cgm_value_full ; cgm_value(i)];
            mean_gap       = round(mean(cgm_value(i:i+1)));
            cgm_value_full = [cgm_value_full ; mean_gap*ones(idx_gap(c,2),1)];

            start_date     = cgm_date(idx_gap(c,1));
            stop_date      = start_date + minutes(idx_gap(c,2)*5);
            to_add_date    = linspace(start_date,stop_date,(idx_gap(c,2)+1));
            cgm_date_full  = [cgm_date_full ; to_add_date'];


            if c ~= n_gap
                c = c+1;
            end
        else
            cgm_value_full = [cgm_value_full ; cgm_value(i)];
            cgm_date_full  = [cgm_date_full ; cgm_date(i)];
        end
    end

    % It is created the vector with the new gap start indexes
    idx_gap_new = [idx_gap(1,1)];
    
    for j = 2:1:n_gap
        idx_gap_new(j,1) = idx_gap(j,1) + sum(idx_gap(1:j-1,2));
    end
    % The start gap index is the one corrisponding to the first sample with
    % the mean value of the gap
    idx_start_gap = idx_gap_new(:,1) + 1;
    idx_end_gap   = idx_gap_new(:,1) + (idx_gap(:,2));
    n_samples     = idx_gap(:,2);
    infogap       = table(idx_start_gap,idx_end_gap,n_samples);

    % Variable name change in order to mantain the same nomenclature
    Time          = cgm_date_full;
    PD_full_table = table(cgm_date_full,cgm_value_full);
    PD_full_table.Properties.VariableNames = {'Time','mg/dL'};
    
    % Output structure
    PD_full.PD      = PD_full_table; 
    PD_full.infogap = infogap;
    PD_full.n_gap   = n_gap;
    
else
    % If no gap were found no correction is needed
    PD_full_table = table(cgm_date,cgm_value);
    PD_full_table.Properties.VariableNames = {'Time','mg/dL'};
    
    % In order to maintain consistency between the two cases, if no gap is
    % detected the same empty variable are created
    idx_start_gap = [];
    idx_end_gap   = [];
    n_samples     = [];
    infogap       = table(idx_start_gap,idx_end_gap,n_samples);
    
    % Output structure
    PD_full.PD      = PD_full_table;
    PD_full.infogap = infogap;
    PD_full.n_gap   = n_gap;
end
    
