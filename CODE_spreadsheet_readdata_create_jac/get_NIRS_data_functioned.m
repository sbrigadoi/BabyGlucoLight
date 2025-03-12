%get NIRS data functioned 
clear all;close all;clc;

%Set params.
subjectN = 48;
hard_drive = "H";
year_of_data = '2022'; %double check each time

event_type = "m_hypo";
%event_type = "S_hypo";
%event_type = "m_hyper";
%event_type = "S_hyper";


% Load folder and pull characters from text file in dir.
NIRS_data_folder = dir(""+hard_drive+":\PadovaPostDoc\BabyGluCo\NIRS_data\PD"+num2str(subjectN)+"");
PDX = load(""+hard_drive+":\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(subjectN)+".mat")
PDX = getfield(PDX,"PD"+num2str(subjectN)+"");


NIRS_file_names_char = [NIRS_data_folder.name];

start_of_PD = strfind(cellstr(NIRS_file_names_char),'PD');
start_of_PD = start_of_PD{1, 1}  ;

start_of_NIRS = strfind(cellstr(NIRS_file_names_char),'.nirs');
start_of_NIRS  = start_of_NIRS {1, 1}  ;

start_of_year = strfind(cellstr(NIRS_file_names_char),year_of_data);
start_of_year  = start_of_year {1, 1}  ;

%% find time and date
%if there are BIS files, split into two parts, non bis and bis
%date_start = 18; %end is +10
%day_start = 14; %end is +2
%time_start = 5; %end is +7

[date_start,day_start,time_start] = get_date_day_time_start(subjectN);

for i=1:size(start_of_PD,2)
    date(i,:) = NIRS_file_names_char(start_of_PD(i)+date_start:start_of_PD(i)+date_start+10);
    day(i,:) = NIRS_file_names_char(start_of_PD(i)+day_start:start_of_PD(i)+day_start+2);
    time(i,:) = NIRS_file_names_char(start_of_PD(i)+time_start:start_of_PD(i)+time_start+7);
end
%% change . to : for data and time
for i=1:size(time,1)
    time_rows(i,:)= strrep(time(i,:),'.',':');
end
% convert time char to datetime time
time_rows = datetime(time_rows, 'Format', 'HH:mm:ss');
%
date_rows = datetime(date, 'Format', 'dd MMM yyyy');
%
DateAndTime_rows = date_rows + timeofday(time_rows);
DateAndTime_rows.Format = [date_rows.Format ' ' time_rows.Format];
%% Create Nirs_file_number indicies
[nirs_file_number_matlab,nirs_file_number,N_seperate_recordings,rec_length] = get_nirs_file_number_matlab(subjectN);

% Sort
for i=1:size(DateAndTime_rows,1)
    DateAndTime_rows_adjusted(i,:) = DateAndTime_rows(i) + (  nirs_file_number(i)*seconds(3000)) ;
end
%didn't run these two lines. can only run these if there is 1 measurenemt.
%PD20 has 3 seperate measurements.
%sorted_nirs_file_number = sort(nirs_file_number) ;
%sorted_DateAndTime_rows_adjusted = sort(DateAndTime_rows_adjusted) ;
%need to add +1 so thffe 0 index can be sorted, then the +1 is taken away
nirs_file_number(:,1)=nirs_file_number(:,1)+1;
sorted_nirs_file_number= sortrows(nirs_file_number,[2,1]); %sort by col2. then by col1.
nirs_file_number(:,1)=nirs_file_number(:,1)-1; %take away the +1, so net +0
sorted_nirs_file_number(:,1)= sorted_nirs_file_number(:,1)-1; %take away the +1, so net +0

%%%%%%%
switch N_seperate_recordings
    case 1
    sorted_DateAndTime_rows_adjusted_1 = sortrows(DateAndTime_rows_adjusted(1:sum(rec_length(1:1))+1,1));

    sorted_DateAndTime_rows_adjusted = [sorted_DateAndTime_rows_adjusted_1];

    case 2
    sorted_DateAndTime_rows_adjusted_1 = sortrows(DateAndTime_rows_adjusted(1:sum(rec_length(1:1))+1,1));
    sorted_DateAndTime_rows_adjusted_2 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:1))+2:sum(rec_length(1:2))+2,1));

    sorted_DateAndTime_rows_adjusted = [sorted_DateAndTime_rows_adjusted_1;sorted_DateAndTime_rows_adjusted_2];

    case 3
    sorted_DateAndTime_rows_adjusted_1 = sortrows(DateAndTime_rows_adjusted(1:sum(rec_length(1:1))+1,1));
    sorted_DateAndTime_rows_adjusted_2 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:1))+2:sum(rec_length(1:2))+2,1));
    sorted_DateAndTime_rows_adjusted_3 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:2))+3:sum(rec_length(1:3))+3,1));

    sorted_DateAndTime_rows_adjusted = [sorted_DateAndTime_rows_adjusted_1;sorted_DateAndTime_rows_adjusted_2;sorted_DateAndTime_rows_adjusted_3];

    case 4
    sorted_DateAndTime_rows_adjusted_1 = sortrows(DateAndTime_rows_adjusted(1:sum(rec_length(1:1))+1,1));
    sorted_DateAndTime_rows_adjusted_2 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:1))+2:sum(rec_length(1:2))+2,1));
    sorted_DateAndTime_rows_adjusted_3 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:2))+3:sum(rec_length(1:3))+3,1));
    sorted_DateAndTime_rows_adjusted_4 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:3))+4:sum(rec_length(1:4))+4,1));

    sorted_DateAndTime_rows_adjusted = [sorted_DateAndTime_rows_adjusted_1;sorted_DateAndTime_rows_adjusted_2;sorted_DateAndTime_rows_adjusted_3;sorted_DateAndTime_rows_adjusted_4];
    
    case 5
    sorted_DateAndTime_rows_adjusted_1 = sortrows(DateAndTime_rows_adjusted(1:sum(rec_length(1:1))+1,1));
    sorted_DateAndTime_rows_adjusted_2 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:1))+2:sum(rec_length(1:2))+2,1));
    sorted_DateAndTime_rows_adjusted_3 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:2))+3:sum(rec_length(1:3))+3,1));
    sorted_DateAndTime_rows_adjusted_4 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:3))+4:sum(rec_length(1:4))+4,1));
    sorted_DateAndTime_rows_adjusted_5 = sortrows(DateAndTime_rows_adjusted(sum(rec_length(1:4))+5:sum(rec_length(1:5))+5,1));

    sorted_DateAndTime_rows_adjusted = [sorted_DateAndTime_rows_adjusted_1;sorted_DateAndTime_rows_adjusted_2;sorted_DateAndTime_rows_adjusted_3;sorted_DateAndTime_rows_adjusted_4;sorted_DateAndTime_rows_adjusted_5];

end

%% Go into m-hypo

%event_type

PDX.DateAndTime_rows.m_hypo = PDX.date_of_events.m_hypo + timeofday(PDX.time_of_events.m_hypo);
PDX.DateAndTime_rows.m_hypo.Format = [PDX.date_of_events.m_hypo.Format ' ' PDX.time_of_events.m_hypo.Format];

%PDX.DateAndTime_rows.m_hyper = PDX.date_of_events.m_hyper + timeofday(PDX.time_of_events.m_hyper);
%PDX.DateAndTime_rows.m_hyper.Format = [PDX.date_of_events.m_hyper.Format ' ' PDX.time_of_events.m_hyper.Format];

for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    %find time diff between glucose eve and NIRS meas.
    delta_time = (PDX.DateAndTime_rows.m_hypo(i,1) - sorted_DateAndTime_rows_adjusted);
    %find only the positive delta. i.e NIRS event is before event
    positiveIndexes = (PDX.DateAndTime_rows.m_hypo(i,1) - sorted_DateAndTime_rows_adjusted) > 0;
    %find smallest pos delta. I.e NIRS event closest to glucose event
    min_pos_delta_time = min(delta_time(positiveIndexes));

    if size(min_pos_delta_time,1) == 0
        min_pos_delta_time = max(delta_time);
    end

    %find index of this.
    start_nirs_file_event(i,1) = find(delta_time==min_pos_delta_time);
    start_PD_nirs_file_event(i,1:2) = sorted_nirs_file_number(find(delta_time==min_pos_delta_time),:);

    %find time diff between glucose eve and NIRS meas.
    delta_time = (PDX.DateAndTime_rows.m_hypo(i,2) - sorted_DateAndTime_rows_adjusted);
    positiveIndexes = (PDX.DateAndTime_rows.m_hypo(i,2) - sorted_DateAndTime_rows_adjusted) > 0;
    min_pos_delta_time = min(delta_time(positiveIndexes));

    if size(min_pos_delta_time,1) == 0
        min_pos_delta_time = max(delta_time);
    end

    end_nirs_file_event(i,1) = find(delta_time==min_pos_delta_time);
    end_PD_nirs_file_event(i,1:2) = sorted_nirs_file_number(find(delta_time==min_pos_delta_time),:);

end

%
%get the length of each glucose event. I.e if it is over 1 or 2 or 3.. NIRS
%file.
for i=1:size(start_PD_nirs_file_event,1)
    size_start_to_end_events(i,1) = size(start_PD_nirs_file_event(i,:):end_PD_nirs_file_event(i,:),2);
end

start_to_end_PD_nirs_file_event = zeros(  size(start_PD_nirs_file_event,1),max(size_start_to_end_events)*2   );

%create a matrix which contains the numbers for the consecutive NIRS files.
%I.e if a glucose events starts at N=3 and ends at N=5, then you will need
%file numbers 3 4 5 for the event.
for i=1:size(start_PD_nirs_file_event,1)
    if start_PD_nirs_file_event(i,1) < end_PD_nirs_file_event(i,1) %to avoid problem of PD49
        start_to_end_PD_nirs_file_event(i,1:size_start_to_end_events(i)) = start_PD_nirs_file_event(i,1):end_PD_nirs_file_event(i,1);
        start_to_end_PD_nirs_file_event(i,size_start_to_end_events(i)+1:end) = start_PD_nirs_file_event(i,2):end_PD_nirs_file_event(i,2);
    end
end
%for PD49, for m hypo5. The start of the event is in file 20 - 2, the end
%of the file is in 1 - 1 , but the : indexing doesn't allow for the end
%index to be a lower number, i.e 1<2. 
%so manually setting PD49 mhypo5.

%%% need to do the same, but for the loc_nirs_file numbers. 
%the loc nirs file start/end
%are the indices that are the start/end of the nirs file order as read in
%matlab. 
%i.e the matlab order is 1 10 11 12 13.... and these indices below point to
%the start/end of the index in the matlab list. 

%intersect(find(nirs_file_number_matlab(:,1) == start_to_end_PD_nirs_file_event(i,1)),find(nirs_file_number_matlab(:,2) == start_to_end_PD_nirs_file_event(i,4)))
%max(size_start_to_end_events)

%start and end of loc file
for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    loc_nirs_file_start(i,1) = intersect(find(nirs_file_number_matlab(:,1) == start_PD_nirs_file_event(i,1)),find(nirs_file_number_matlab(:,2) == start_PD_nirs_file_event(i,2)));
    loc_nirs_file_end(i,1) = intersect(find(nirs_file_number_matlab(:,1) == end_PD_nirs_file_event(i,1)),find(nirs_file_number_matlab(:,2) == end_PD_nirs_file_event(i,2)));
end   

%the below fills in the gaps of the start-end loc files.
%creating matrix to hold the indices of the loc nirs file start end
%(normally N events rows x number of NIRS files cols (12 x 3)
loc_nirs_file_start_end = zeros(size(PDX.DateAndTime_rows.m_hypo,1),max(size_start_to_end_events));

for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    for j=1:max(size_start_to_end_events)
        if start_to_end_PD_nirs_file_event(i,1) > 0
            %loc_nirs_file_start_end
            loc_nirs_file_start_end(i,j) = intersect(find(nirs_file_number_matlab(:,1) == start_to_end_PD_nirs_file_event(i,j)),find(nirs_file_number_matlab(:,2) == start_to_end_PD_nirs_file_event(i,j+max(size_start_to_end_events))));
        end
        
    end
end
% X.ALL
%delta_start = zeros(size(PD20.DateAndTime_rows.m_hypo,1),1);
%delta_end = zeros(size(PD20.DateAndTime_rows.m_hypo,1),1);
for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    delta_start(i,1) = PDX.DateAndTime_rows.m_hypo(i,1) - sorted_DateAndTime_rows_adjusted(start_nirs_file_event(i,1));
    delta_end(i,1) = PDX.DateAndTime_rows.m_hypo(i,2) - sorted_DateAndTime_rows_adjusted(end_nirs_file_event(i,1));
    
    frame_start(i,1) = floor(delta_start(i,1)/seconds(3000)*30001);
    frame_end(i,1) =floor(delta_end(i,1)/seconds(3000)*30001);
end
%name_NIRS_start
%join(name_NIRS_start,name_NIRS_end)
%NIRS_data_folder_onlydata = NIRS_data_folder(3:end,:);
%loc_nirs_file_start_end
%NIRS_data_folder_onlydata.name(loc_nirs_file_start_end(i,j),1)
name_times_N = unique(time,'rows','stable');
name_dates_N = unique(date,'rows','stable');
name_days_N = unique(day,'rows','stable');
%
%name_start = 'E:\PadovaPostDoc\BabyGluCo\NIRS_data\PD20\PD20';
name_start = ''+hard_drive+':\PadovaPostDoc\BabyGluCo\NIRS_data\PD'+num2str(subjectN)+'\PD'+num2str(subjectN)+'';

switch subjectN
    case 5
        name_start = ''+hard_drive+':\PadovaPostDoc\BabyGluCo\NIRS_data\PD'+num2str(subjectN)+'\MS PD'+num2str(subjectN)+'';
end

name_start = convertStringsToChars(name_start);


names_glucose_event = cell(size(PDX.DateAndTime_rows.m_hypo,1),max(size_start_to_end_events));
for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    if size_start_to_end_events(i) > 0 %if not any overlapping files, skip
        if start_to_end_PD_nirs_file_event(i,1) > 0 % if not valid/missing overlap, skip
    
            for j=1:max(size_start_to_end_events)
                name_time = name_times_N(start_to_end_PD_nirs_file_event(i,j+max(size_start_to_end_events)),:); %was +3 08 05 24
                name_date = name_dates_N(start_to_end_PD_nirs_file_event(i,j+max(size_start_to_end_events)),:);
                name_day = name_days_N(start_to_end_PD_nirs_file_event(i,j+max(size_start_to_end_events)),:);
                name_pdN = num2str(start_to_end_PD_nirs_file_event(i,j));
                name_end = '.nirs';
            %names_glucose_event(i,j) = [name_start ' ' name_time ' ' name_date ' ' name_pdN name_end];
                names_glucose_event{i,j} = [name_start ' ' name_time ' ' name_day ' ' name_date ' ' name_pdN name_end];
            end
        end
    end
end

% Load events PD X
%test_load = load(names_glucose_event{1,1},'-mat');

%Need to select the events that are actually covered by the NIRS data.
%I.e. For the 1 or 2 or 3 or 4 PD numbers, they should be listed as
%consectuive and increasing. 
%i.e
%PD 3 4 5 is a legit coverage of 3 NIRS PD Numbers (3,4,5)
%PD 3 4 1 is coverage of 2 NIRS PD numbers (3 and 4) 
%PD 9 3 3 is coverage of NO NIRS PD Numbers

%best way to check this is by looking at 
%start_PD_nirs_file_event and
%end_PD_nirs_file_event

%startPD +1 or 2 should = endPD

%start_to_end_PD_nirs_file_event

%if start_to_end_PD_nirs_file_event(i,j) > start_to_end_PD_nirs_file_event(i,j+2)

%end

%if start_to_end_PD_nirs_file_event(i,j) > start_to_end_PD_nirs_file_event(i,j+2)

%end

%best way to check if the start PD or endPD number are covered by the NIRS
%data are that if the delta time for the start/end are less than 00:50:00
%(50 mins) as this is the longest possible delta time, since each NIRS
%reading is a max of 50 mins. 
%and the start/end seperate PD lists will have the first (start) and last
%(end) NIRS numbers in there.

length_NIRS_PD = duration(0,50,0);

% valid_start_end = zeros(size(PD20.DateAndTime_rows.m_hypo,1),2);
% for i=1:size(PD20.DateAndTime_rows.m_hypo,1)
%     if delta_start(i) < length_NIRS_PD
%         valid_start_end(i,1) = 1;
%     end
%     if delta_end(i)< length_NIRS_PD
%         valid_start_end(i,2) = 1;
%     end        
% end

valid_start_end = zeros(size(PDX.DateAndTime_rows.m_hypo,1),2);
for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    if (delta_start(i) < length_NIRS_PD) && (delta_start(i) > 0) 
        valid_start_end(i,1) = 1;
    end
    if (delta_end(i)< length_NIRS_PD) && (delta_start(i) > 0)
        valid_start_end(i,2) = 1;
    end        
end


%start_end_pd 
%test to see how many of the PD numbers contribute to the event
%i.e if it's 1, 2, 3 or 4 of the PD numbers. 

%start_to_end_PD_nirs_file_event(1,1) < start_to_end_PD_nirs_file_event(1,3)
                                       
%if start_to_end_PD_nirs_file_event(i,j)+1 = start_to_end_PD_nirs_file_event(i,j+1)
%    valid_start_end_pdN =

%start_to_end_PD_nirs_file_event

valid_start_end_pdN = zeros(size(PDX.DateAndTime_rows.m_hypo,1),max(size_start_to_end_events));
%
% 1= valid PD number
% 0= not valid PD number
for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    for j=1:max(size_start_to_end_events)-1 %was 2
        if start_to_end_PD_nirs_file_event(i,j)+1 == start_to_end_PD_nirs_file_event(i,j+1);
            valid_start_end_pdN(i,j+1)=1;
        end
    end
end

for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
     if valid_start_end_pdN(i,2) == 1;
            valid_start_end_pdN(i,1)=1;
     end
end
%if start_to_end_PD_nirs_file_event(i,j)+2 == start_to_end_PD_nirs_file_event(i,j+2);
%    valid_start_end_pdN(i,j+2)=1;
%end
%
for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    if sum(valid_start_end(i,:)) == 0
        valid_start_end_pdN(i,:) = 0; 
    end
end
%% SAVE load nirs files mhypo PDX CHANGE NAME AT END

%names_glucose_event
%we need to check that the NIRS files are all complete, i.e they have
%30,001 frames in them, and that the end NIRS file has at least frame_end
%number of frames
%If not, then we can set valid_start_end_pdN(i,:) to zero for that event
%(row) so that the all NIRS data isn't created and save operation not done.


for i=1:size(PDX.DateAndTime_rows.m_hypo,1)
    %for j=1:size(valid_start_end_pdN,2)
        if sum(valid_start_end_pdN(i,:)) == 1
            nirs_start = load(names_glucose_event{i,1},'-mat');
            if size(nirs_start.t,1) < 30001 %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_start.d(frame_end(i):end,:) ];
            end
        end
    %end
        if sum(valid_start_end_pdN(i,:)) == 2
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_end = load(names_glucose_event{i,2},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_end.d(1:frame_end(i),:)];
            end
        end
        
        if sum(valid_start_end_pdN(i,:)) == 3
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle = load(names_glucose_event{i,2},'-mat');
            nirs_end = load(names_glucose_event{i,3},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle.t,1) < 30001 |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle.d ; nirs_end.d(1:frame_end(i),:)];
            end
        end
        
        if sum(valid_start_end_pdN(i,:)) == 4
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_end = load(names_glucose_event{i,4},'-mat');
             if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
             else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d  ; nirs_end.d(1:frame_end(i),:)];
             end
        end
    
        if sum(valid_start_end_pdN(i,:)) == 5
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_middle_3 = load(names_glucose_event{i,4},'-mat');
            nirs_end = load(names_glucose_event{i,5},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  | size(nirs_middle_3.t,1) <30001 |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d ; nirs_middle_3.d  ; nirs_end.d(1:frame_end(i),:)];
            end
         end

        if sum(valid_start_end_pdN(i,:)) == 6
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_middle_3 = load(names_glucose_event{i,4},'-mat');
            nirs_middle_4 = load(names_glucose_event{i,5},'-mat');
            nirs_end = load(names_glucose_event{i,6},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  | size(nirs_middle_3.t,1) <30001 | size(nirs_middle_4.t,1) <30001  |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d ;nirs_middle_3.d ;nirs_middle_4.d ; nirs_end.d(1:frame_end(i),:)];
            end
        end

        if sum(valid_start_end_pdN(i,:)) == 7
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_middle_3 = load(names_glucose_event{i,4},'-mat');
            nirs_middle_4 = load(names_glucose_event{i,5},'-mat');
            nirs_middle_5 = load(names_glucose_event{i,6},'-mat');
            nirs_end = load(names_glucose_event{i,7},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  | size(nirs_middle_3.t,1) <30001 | size(nirs_middle_4.t,1) <30001 | size(nirs_middle_5.t,1) <30001 |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d ;nirs_middle_3.d ;nirs_middle_4.d ;nirs_middle_5.d ; nirs_end.d(1:frame_end(i),:)];
            end
        end
        
        if sum(valid_start_end_pdN(i,:)) == 8
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_middle_3 = load(names_glucose_event{i,4},'-mat');
            nirs_middle_4 = load(names_glucose_event{i,5},'-mat');
            nirs_middle_5 = load(names_glucose_event{i,6},'-mat');
            nirs_middle_6 = load(names_glucose_event{i,7},'-mat');
            nirs_end = load(names_glucose_event{i,8},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  | size(nirs_middle_3.t,1) <30001 |size(nirs_middle_4.t,1) <30001 |size(nirs_middle_5.t,1) <30001 |size(nirs_middle_6.t,1) <30001    |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d ;nirs_middle_3.d ;nirs_middle_4.d ;nirs_middle_5.d ;nirs_middle_6.d ; nirs_end.d(1:frame_end(i),:)];
            end
        end

        if sum(valid_start_end_pdN(i,:)) == 9
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_middle_3 = load(names_glucose_event{i,4},'-mat');
            nirs_middle_4 = load(names_glucose_event{i,5},'-mat');
            nirs_middle_5 = load(names_glucose_event{i,6},'-mat');
            nirs_middle_6 = load(names_glucose_event{i,7},'-mat');
            nirs_middle_7 = load(names_glucose_event{i,8},'-mat');
            nirs_end = load(names_glucose_event{i,9},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  | size(nirs_middle_3.t,1) <30001 |size(nirs_middle_4.t,1) <30001 |size(nirs_middle_5.t,1) <30001 |size(nirs_middle_6.t,1) <30001 |size(nirs_middle_7.t,1) <30001 |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d ;nirs_middle_3.d ;nirs_middle_4.d ;nirs_middle_5.d ;nirs_middle_6.d;nirs_middle_7.d ; nirs_end.d(1:frame_end(i),:)];
            end
        end

        if sum(valid_start_end_pdN(i,:)) == 10
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_middle_3 = load(names_glucose_event{i,4},'-mat');
            nirs_middle_4 = load(names_glucose_event{i,5},'-mat');
            nirs_middle_5 = load(names_glucose_event{i,6},'-mat');
            nirs_middle_6 = load(names_glucose_event{i,7},'-mat');
            nirs_middle_7 = load(names_glucose_event{i,8},'-mat');
            nirs_middle_8 = load(names_glucose_event{i,9},'-mat');
            nirs_end = load(names_glucose_event{i,10},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  | size(nirs_middle_3.t,1) <30001 |size(nirs_middle_4.t,1) <30001 |size(nirs_middle_5.t,1) <30001 |size(nirs_middle_6.t,1) <30001 |size(nirs_middle_7.t,1) <30001 |size(nirs_middle_8.t,1) <30001 |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d ;nirs_middle_3.d ;nirs_middle_4.d ;nirs_middle_5.d ;nirs_middle_6.d;nirs_middle_7.d;nirs_middle_8.d   ; nirs_end.d(1:frame_end(i),:)];
            end
        end

    if sum(valid_start_end_pdN(i,:)) == 11
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_middle_3 = load(names_glucose_event{i,4},'-mat');
            nirs_middle_4 = load(names_glucose_event{i,5},'-mat');
            nirs_middle_5 = load(names_glucose_event{i,6},'-mat');
            nirs_middle_6 = load(names_glucose_event{i,7},'-mat');
            nirs_middle_7 = load(names_glucose_event{i,8},'-mat');
            nirs_middle_8 = load(names_glucose_event{i,9},'-mat');
            nirs_middle_9 = load(names_glucose_event{i,10},'-mat');
            nirs_end = load(names_glucose_event{i,11},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  | size(nirs_middle_3.t,1) <30001 |size(nirs_middle_4.t,1) <30001 |size(nirs_middle_5.t,1) <30001 |size(nirs_middle_6.t,1) <30001 |size(nirs_middle_7.t,1) <30001 |size(nirs_middle_8.t,1) <30001 |size(nirs_middle_9.t,1) <30001 |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d ;nirs_middle_3.d ;nirs_middle_4.d ;nirs_middle_5.d ;nirs_middle_6.d;nirs_middle_7.d;nirs_middle_8.d ;nirs_middle_9.d ; nirs_end.d(1:frame_end(i),:)];
            end
        end

    if sum(valid_start_end_pdN(i,:)) == 14
            nirs_start = load(names_glucose_event{i,1},'-mat');
            nirs_middle_1 = load(names_glucose_event{i,2},'-mat');
            nirs_middle_2 = load(names_glucose_event{i,3},'-mat');
            nirs_middle_3 = load(names_glucose_event{i,4},'-mat');
            nirs_middle_4 = load(names_glucose_event{i,5},'-mat');
            nirs_middle_5 = load(names_glucose_event{i,6},'-mat');
            nirs_middle_6 = load(names_glucose_event{i,7},'-mat');
            nirs_middle_7 = load(names_glucose_event{i,8},'-mat');
            nirs_middle_8 = load(names_glucose_event{i,9},'-mat');
            nirs_middle_9 = load(names_glucose_event{i,10},'-mat');
            nirs_middle_10 = load(names_glucose_event{i,11},'-mat');
            nirs_middle_11 = load(names_glucose_event{i,12},'-mat');
            nirs_middle_12 = load(names_glucose_event{i,13},'-mat');
            nirs_end = load(names_glucose_event{i,14},'-mat');
            if size(nirs_start.t,1) < 30001 | size(nirs_middle_1.t,1) < 30001 |size(nirs_middle_2.t,1) <30001  | size(nirs_middle_3.t,1) <30001 |size(nirs_middle_4.t,1) <30001|size(nirs_middle_5.t,1) <30001 |size(nirs_middle_6.t,1) <30001 |size(nirs_middle_7.t,1) <30001 |size(nirs_middle_8.t,1) <30001 |size(nirs_middle_9.t,1) <30001 |size(nirs_middle_10.t,1) <30001 |size(nirs_middle_11.t,1) <30001 |size(nirs_middle_12.t,1) <30001      |size(nirs_end.t,1) < frame_end(i) %if file is not complete
                valid_start_end_pdN(i,:) = 0; %make it invalid
            else
                all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle_1.d ; nirs_middle_2.d ;nirs_middle_3.d ;nirs_middle_4.d ;nirs_middle_5.d ;nirs_middle_6.d;nirs_middle_7.d;nirs_middle_8.d  ;nirs_middle_9.d ;nirs_middle_10.d;nirs_middle_11.d;nirs_middle_12.d; nirs_end.d(1:frame_end(i),:)];
            end
        end

    %all_NIRS_data = [nirs_start.d(frame_start(i):end,:) ; nirs_middle.d ; nirs_end.d(1:frame_end(i),:)];
    
    if sum(valid_start_end_pdN(i,:)) > 0;
        PDX_m_hypo_1 = nirs_start;
        PDX_m_hypo_1.d =all_NIRS_data;
        PDX_m_hypo_1.t= [0:1:(size(all_NIRS_data,1))-1]'*0.1;
        PDX_m_hypo_1.s = zeros(size(all_NIRS_data,1),1);
        PDX_m_hypo_1.aux = zeros(size(all_NIRS_data,1),8);

        %change names

        switch subjectN
            case 1
            PD1_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 2
            PD2_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 5
            PD5_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 9
            PD9_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 10
            PD10_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 11
            PD11_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 28
            PD28_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 34
            PD34_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 35
            PD35_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 36
            PD36_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 39
            PD39_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 41
            PD41_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 43
            PD43_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 44
            PD44_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 45
            PD45_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 46
            PD46_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 47
            PD47_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 48
            PD48_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX


            case 55
            PD55_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 56
            PD56_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 57
            PD57_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 58
            PD58_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 59
            PD59_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX

            case 60
            PD60_m_hypo_X = PDX_m_hypo_1; %always keep RHS as PDX


        end
        save("PD"+num2str(subjectN)+"_"+event_type+"_"+i+".nirs","PD"+num2str(subjectN)+"_"+event_type+"_X","-mat")
    

    end
end

%% Test load

%load("PD41_m_hypo_2.nirs","-mat")


