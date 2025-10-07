% % CGM - NIRS SYNCHRONIZATION
% 
% CASO SPECIALE: Ha 2 database, il secondo ha i numeri sbagliati. Il file
% 5 è in realtà il file 0. E così via. Si è deciso di farlo a mano per
% evitare problemi
% 
% close all
% clear all
% clc
% 
% idx_PD = 4;
% 
%% SETTING OF FOLDER PATH 
% 
% base_path      = pwd;
% data_path_cgm  = fullfile(base_path,'DATI','CGM_INTERP');
% data_path_nirs = fullfile('E:\DATI\NIRS\PD4');
% function_path  = fullfile(base_path,'CODICE','FUNCTION');
% 
% addpath(genpath(pwd))
% addpath(data_path_cgm)
% addpath(function_path)
% addpath(data_path_nirs)

%% PRIMA INTERVALLO DI EUGLICEMIA

% % Load of cgm data
% PD_struct = load ('PD4.mat');
% date_CGM  = PD_struct.PD.PD.Time;
% cgm_value = PD_struct.PD.PD.("mg/dL");
% 
% % Extraction of the start and finish date of the CGM acquisition
% date_CGM_start = date_CGM(1);
% date_CGM_end = date_CGM(end);
% 
% % Display of the preliminary inspection of the CGM dataset
% disp('==============================================================================')
% disp('CGM INTERVAL DATE')
% disp(['Start : ',datestr(date_CGM_start)])
% disp(['End   : ',datestr(date_CGM_end)])
% 
% %% 
% % Load of the usefull glycemic analysis result
% load idx_eu_results.mat
% 
% PD_results = idx_eu_int(idx_eu_int.pax_idx_vector==idx_PD, : );
% 
% % Date interval of the first and last Euglycemia window
% date_start_first_eu = date_CGM(PD_results.idx_start_first_eu_w);
% date_end_first_eu   = date_CGM(PD_results.idx_end_first_eu_w);
% date_start_last_eu  = date_CGM(PD_results.idx_start_last_eu_w);
% date_end_last_eu    = date_CGM(PD_results.idx_end_last_eu_w);
% dur_first_eu        = minutes(date_end_first_eu - date_start_first_eu);
% dur_last_eu         = minutes(date_end_last_eu - date_start_last_eu);
% 
% % Display of the first and last Euglycemia window dates
% disp(' ')
% disp('FIRST EUGLYCEMIA INTERVAL')
% disp(['Start    : ',datestr(date_start_first_eu)])
% disp(['End      : ',datestr(date_end_first_eu)])
% disp(['Duration : ',num2str(dur_first_eu)])
% disp(' ')
% disp('LAST EUGLYCEMIA INTERVAL')
% disp(['Start    : ',datestr(date_start_last_eu)])
% disp(['End      : ',datestr(date_end_last_eu)])
% disp(['Duration : ',num2str(dur_last_eu)])
% 
% %% NIRS DATASET ANALYSES
% 
% % It is analysed as automatically as possible the NIRS dataset in order to
% % extract the acquisition information in order to proceed with the CGM-NIRS
% % synchronization.
% 
% % Count of NIRS files for each dataset.
% % It is generated the path name of the folder which contains the NIRS
% % files. The command 'dir' extracts the informations of all files which
% % contains inside the name the words passed as argument.
% tmp = fullfile(data_path_nirs,'*PD4 18.16*.nirs');
% info_first_dataset = dir(tmp);
% n_dataset_1 = length(info_first_dataset);
% 
% % NIRS sampling period
% Ts_NIRS = 0.1; % [s]
% % Standard length of a NIRS file
% n_standard = 30001;
% 
% %% EXTRACTION OF ACQUISITION INTERVAL DATE
% 
% % Extraction of the acquisition start date
% 
% % Extraction of the file name
% file_name_first = info_first_dataset.name;
% % Separation of the information contained inside the file name and
% % conversion in DateTime format
% sep = split(string(file_name_first));
% time_first = split(sep(2),'.');
% datetime_format = strcat(sep(4),'-',sep(5),'-',sep(6),{' '},time_first(1),':',time_first(2),':',time_first(3));
% date_start_first = datetime(datetime_format);
% 
% % In order to get the length of the last acqusition file it is necessary to
% % load it
% load('PD4 18.16.25 Sun 12 Apr 2020 113.nirs','-mat');
% n_last = length(t);
% 
% % Total number of samples acquired in the first dataset and convertion in
% % seconds
% n1_complete = (n_dataset_1-1)*n_standard;
% n1_tot = n1_complete + n_last;
% duration_first = seconds(Ts_NIRS*n1_tot); %[s]
% 
% % End date of first dataset acquisition
% date_end_first = date_start_first + duration_first;
% 
% % Display of the NIRS acqusition interval dates 
% disp('==============================================================================')
% disp('NIRS INTERVAL DATE')
% disp(' ')
% disp('NIRS 1 :')
% disp(['Start : ',datestr(date_start_first)])
% disp(['End   : ',datestr(date_end_first)])
% 
% a = 0;
% %% CGM AND NIRS SYNCHRONIZATION
% 
% % By the use of the 'date_sincro_1_NIRS' it is extracted the intreval date
% % of the CGM-NIRS synchronization in case of 1 NIRS dataset.
% 
% f_disp = 0; % Display Command Window results
% [date_start_sincro,date_end_sincro] = date_sincro_1_NIRS(date_CGM,date_start_first,date_end_first,f_disp);
% 
% % The same approach is used to determine whether the start and the end of
% % the first and last Euglycemia window is inside the NIRS interval
% 
% f_disp = 1; % NO display Command Window results
% [date_start_sincro_first,date_end_sincro_first] = date_sincro_1_NIRS(date_CGM(PD_results.idx_start_first_eu_w:PD_results.idx_end_first_eu_w),date_start_first,date_end_first,f_disp);
% % [date_start_sincro_last,date_end_sincro_last]  = date_sincro_1_NIRS(date_CGM(PD_results.idx_start_last_eu_w:PD_results.idx_end_last_eu_w),date_start_first,date_end_first,f_disp);
% a = 0;
% 
% %% ESTRAZIONE INDICI DI INIZIO E FINE SINCRONIZZAZIONE
% 
% % Si verifica dove si colloca la data di inizio sincronizzazione e si
% % calcola l'intervallo di tempo che intercorre tra tale data e l'inizio 
% % dell'acquisizione NIRS. Si salva il nome del dataset NIRS nel quale cade
% % la data di inizio sincronizzazione.
% 
% if date_start_sincro < date_end_first 
%     diff_Nirs_sincro = date_start_sincro - date_start_first;
%     date_start_Nirs = date_start_first;
%     dataset_name = 'PD4 18.16.25 Sun 12 Apr 2020 %d.nirs';
% end
% 
% % Durata in secondi di ogni file NIRS campionato completo (tranne l'ultimo)
% dur_Nirs_file_sec = seconds(n_standard*0.1); %[s]
% 
% % Si calcola la data di inizio acquisizione per ogni file NIRS del dataset
% % in analisi.
% date_acq_Nirs_file = [];
% for i=0:(n_dataset_1-1)
%     date_acq_Nirs_file = [date_acq_Nirs_file ; date_start_Nirs+(i*dur_Nirs_file_sec)];
% end
% 
% % Si estrae il nome del file NIRS contenente l'istante di inizio
% % sincronizzazione
% conv = 0;
% i = 0;
% % Si scorre il file contenente la data di inizio acquisizione per ogni file
% % NIRS del dataset contenente la data di inizio sincronizzazione
% while conv == 0
%     i = i+1;
%     % Se la data di inizio acquisizione del file i-esimo è maggiore della
%     % data di inizio sincronizzazione allora l'inidice del file che
%     % contiene la data ricercata è il precedente 
%     if date_acq_Nirs_file(i) > date_start_sincro
%         idx_file_Nirs_start_sincro = i - 1;
%         conv = 1;
%     end
%     % Se la data di inizio acquisizione del file i-esimo è uguale alla data
%     % di inizio sincronizzazione allora l'inidice del file che contiene la 
%     % data ricercata è il corrente
%     if date_acq_Nirs_file(i) == date_start_sincro
%         idx_file_Nirs_start_sincro = i;
%         conv = 1;
%     end
%     % Se la data di inizio acquisizione del file i-esimo è minore della 
%     % data di inizio sincronizzazione e si è arrivati all'ultimo file NIRS
%     % allora l'inidice del file che contiene la data ricercata è il 
%     % corrente    
%     if date_acq_Nirs_file(i) < date_start_sincro & i == length(date_acq_Nirs_file)
%         idx_file_Nirs_start_sincro = i;
%         conv = 1;
%     end
% 
%     % Una volta trovato l'indice del file contenente il punto di inizio
%     % sincronizzazione si crea la variabile che contiene il nome del file
%     if conv == 1
%         % Data file Nirs per inizio sincro
%         time_file_Nirs_start_sincro = date_acq_Nirs_file(idx_file_Nirs_start_sincro);
%         % Nome file Nirs per inizio sincro, si sottrae uno dell'indice
%         % trovato in precedenza poichè la nomenclatura dei file inizia da
%         % zero.
%         name_file_Nirs_start_sincro = sprintf(dataset_name, idx_file_Nirs_start_sincro-1);
%     end
% end
% 
% % Si effettua lo stesso procedimento per identificare il nome del file
% % contenente la data di fine sincronizzazione
% conv = 0;
% i = 0;
% while conv == 0
%     i = i+1;
%     if date_acq_Nirs_file(i) > date_end_sincro
%         idx_file_Nirs_end_sincro = i - 1;
%         conv = 1;
%     end
%     if date_acq_Nirs_file(i) == date_end_sincro
%         idx_file_Nirs_end_sincro = i;
%         conv = 1;
%     end
%     if date_acq_Nirs_file(i) < date_end_sincro & i == length(date_acq_Nirs_file)
%         idx_file_Nirs_end_sincro = i;
%         conv = 1;
%     end
% 
%     if conv == 1
%         time_file_Nirs_end_sincro = date_acq_Nirs_file(idx_file_Nirs_end_sincro);
%         name_file_Nirs_end_sincro = sprintf(dataset_name, idx_file_Nirs_end_sincro-1);
%     end
% end
% 
% % Si calcola l'indice di inizio sincronizzazione all'interno del file NIRS 
% % contentenete la data di inizio sincronizzazione
% diff_start_Nirs = date_start_sincro - time_file_Nirs_start_sincro;
% diff_start_Nirs_s = seconds(diff_start_Nirs); %[s]
% idx_start_sincro_Nirs = round(diff_start_Nirs_s/Ts_NIRS); %[time point]
% 
% % Si calcola l'indice di fine sincronizzazione all'interno del file NIRS 
% % contentenete la data di fine sincronizzazione 
% diff_end_Nirs = date_end_sincro - time_file_Nirs_end_sincro;
% diff_end_Nirs_s = seconds(diff_end_Nirs);
% idx_end_sincro_Nirs = round(diff_end_Nirs_s/Ts_NIRS);
% 
% % Display informazioni sulla sincronizzazione lato database NIRS
% disp('==============================================================================')
% disp('NIRS - CGM SYNCHRONIZATION')
% disp(' ')
% disp(['NIRS FILE INFO FOR SINCRO FROM ',datestr(date_start_sincro),' /TO/ ',datestr(date_end_sincro)])
% disp(' ')
% disp('START')
% disp(['Name   : ',name_file_Nirs_start_sincro])
% disp(['Time   : ',datestr(time_file_Nirs_start_sincro)])
% disp(['Sample : ',num2str(idx_start_sincro_Nirs)])
% disp(' ')
% disp('END')
% disp(['Name   : ',name_file_Nirs_end_sincro])
% disp(['Time   : ',datestr(time_file_Nirs_end_sincro)])
% disp(['Sample : ',num2str(idx_end_sincro_Nirs)])
% 
% %% ESTRAZIONE INDICI DI INIZIO E FINE SINCRONIZZAZIONE PRIMA EU
% 
% % Si verifica dove si colloca la data di inizio sincronizzazione e si
% % calcola l'intervallo di tempo che intercorre tra tale data e l'inizio 
% % dell'acquisizione NIRS. Si salva il nome del dataset NIRS nel quale cade
% % la data di inizio sincronizzazione.
% if date_start_sincro_first < date_end_first 
%     diff_Nirs_sincro = date_start_sincro_first - date_start_first;
%     date_start_Nirs = date_start_first;
%     dataset_name = 'PD4 18.16.25 Sun 12 Apr 2020 %d.nirs';
% end
% 
% % Durata in secondi di ogni file NIRS campionato completo (tranne l'ultimo)
% dur_Nirs_file_sec = seconds(n_standard*0.1); %[s]
% 
% % Si calcola la data di inizio acquisizione per ogni file NIRS del dataset
% % in analisi.
% date_acq_Nirs_file = [];
% for i=0:(n_dataset_1-1)
%     date_acq_Nirs_file = [date_acq_Nirs_file ; date_start_Nirs+(i*dur_Nirs_file_sec)];
% end
% 
% % Si estrae il nome del file NIRS contenente l'istante di inizio
% % sincronizzazione
% conv = 0;
% i = 0;
% % Si scorre il file contenente la data di inizio acquisizione per ogni file
% % NIRS del dataset contenente la data di inizio sincronizzazione
% while conv == 0
%     i = i+1;
%     % Se la data di inizio acquisizione del file i-esimo è maggiore della
%     % data di inizio sincronizzazione allora l'inidice del file che
%     % contiene la data ricercata è il precedente 
%     if date_acq_Nirs_file(i) > date_start_sincro_first
%         idx_file_Nirs_start_sincro = i - 1;
%         conv = 1;
%     end
%     % Se la data di inizio acquisizione del file i-esimo è uguale alla data
%     % di inizio sincronizzazione allora l'inidice del file che contiene la 
%     % data ricercata è il corrente
%     if date_acq_Nirs_file(i) == date_start_sincro_first
%         idx_file_Nirs_start_sincro = i;
%         conv = 1;
%     end
%     % Se la data di inizio acquisizione del file i-esimo è minore della 
%     % data di inizio sincronizzazione e si è arrivati all'ultimo file NIRS
%     % allora l'inidice del file che contiene la data ricercata è il 
%     % corrente    
%     if date_acq_Nirs_file(i) < date_start_sincro_first & i == length(date_acq_Nirs_file)
%         idx_file_Nirs_start_sincro = i;
%         conv = 1;
%     end
% 
%     % Una volta trovato l'indice del file contenente il punto di inizio
%     % sincronizzazione si crea la variabile che contiene il nome del file
%     if conv == 1
%         % Data file Nirs per inizio sincro
%         time_file_Nirs_start_sincro_first = date_acq_Nirs_file(idx_file_Nirs_start_sincro);
%         % Nome file Nirs per inizio sincro, si sottrae uno dell'indice
%         % trovato in precedenza poichè la nomenclatura dei file inizia da
%         % zero.
%         name_file_Nirs_start_sincro_first = sprintf(dataset_name, idx_file_Nirs_start_sincro-1);
%     end
% end
% 
% % Si effettua lo stesso procedimento per identificare il nome del file
% % contenente la data di fine sincronizzazione
% conv = 0;
% i = 0;
% while conv == 0
%     i = i+1;
%     if date_acq_Nirs_file(i) > date_end_sincro_first
%         idx_file_Nirs_end_sincro = i - 1;
%         conv = 1;
%     end
%     if date_acq_Nirs_file(i) == date_end_sincro_first
%         idx_file_Nirs_end_sincro = i;
%         conv = 1;
%     end
%     if date_acq_Nirs_file(i) < date_end_sincro_first & i == length(date_acq_Nirs_file)
%         idx_file_Nirs_end_sincro = i;
%         conv = 1;
%     end
% 
%     if conv == 1
%         time_file_Nirs_end_sincro_first = date_acq_Nirs_file(idx_file_Nirs_end_sincro);
%         name_file_Nirs_end_sincro_first = sprintf(dataset_name, idx_file_Nirs_end_sincro-1);
%     end
% end
% 
% % Si calcola l'indice di inizio sincronizzazione all'interno del file NIRS 
% % contentenete la data di inizio sincronizzazione che nel PD15 corrisponde 
% % alla data del primo campione CGM valido (16-Jul_2020 14:05:00) 
% diff_start_Nirs = date_start_sincro_first - time_file_Nirs_start_sincro_first;
% diff_start_Nirs_s = seconds(diff_start_Nirs); %[s]
% idx_start_sincro_Nirs_first = round(diff_start_Nirs_s/Ts_NIRS); %[time point]
% 
% % Si calcola l'indice di fine sincronizzazione all'interno del file NIRS 
% % contentenete la data di fine sincronizzazione che nel PD15 corrisponde 
% % alla data dell'ultimo campione CGM valido e contenuto contemporaneamente
% % in un file NIRS (20-Jul_2020 11:40:00)
% diff_end_Nirs = date_end_sincro_first - time_file_Nirs_end_sincro_first;
% diff_end_Nirs_s = seconds(diff_end_Nirs);
% idx_end_sincro_Nirs_first = round(diff_end_Nirs_s/Ts_NIRS);
% 
% % Display informazioni sulla sincronizzazione lato database NIRS
% disp('==============================================================================')
% disp('FIRST EUGLYCEMIA INTERVAL SYNCHRONIZATION')
% disp(' ')
% disp(['NIRS FILE INFO FOR SINCRO FROM ',datestr(date_start_sincro_first),' /TO/ ',datestr(date_end_sincro_first)])
% disp(' ')
% disp('START')
% disp(['Name   : ',name_file_Nirs_start_sincro_first])
% disp(['Time   : ',datestr(time_file_Nirs_start_sincro_first)])
% disp(['Sample : ',num2str(idx_start_sincro_Nirs_first)])
% disp(' ')
% disp('END')
% disp(['Name   : ',name_file_Nirs_end_sincro_first])
% disp(['Time   : ',datestr(time_file_Nirs_end_sincro_first)])
% disp(['Sample : ',num2str(idx_end_sincro_Nirs_first)])
% 
% 
% a= 0;
% %% LOAD AND SAVE OF REQUIRED NIRS DATA 
% 
% % Da cambiare per altri PAX
% %                                                                           SPECIFIC FOR PD4
% idx_file_Nirs_start_sincro_first = 5;
% idx_file_Nirs_end_sincro_first = 23;
% 
% % idx_file_Nirs_start_sincro_last = 110;
% % idx_file_Nirs_end_sincro_last = 115;
% 
% % First EU
% if date_start_first_eu < date_end_first 
%     dataset_name_first = 'PD4 18.16.25 Sun 12 Apr 2020 %d.nirs';
% end
% 
% a = 0;
% %% FIRST EUGLYCEMIA INTERVAL
% 
% % Load of the NIRS file for the first Euglycemia interval
% disp('===================================================================')
% disp('Load NIRS file first Euglycemia interval')
% disp(' ')
% 
% % Inizializzo 
% d_tot = [];
% t_tot = [];
% aux_tot = [];
% s_tot = [];
% 
% name_file_Nirs_start_sincro_first = 'PD4 18.16.25 Sun 12 Apr 2020 5.nirs';
% % % Carico il primo file
% first_start_file_name = name_file_Nirs_start_sincro_first;
% disp(['Load: ',first_start_file_name])
% load(first_start_file_name,'-mat');
% d_tot   = vertcat(d_tot,d(idx_start_sincro_Nirs_first:end,:));
% t_tot   = vertcat(t_tot,t(idx_start_sincro_Nirs_first:end,:));
% aux_tot = vertcat(aux_tot,aux(idx_start_sincro_Nirs_first:end,:));
% s_tot   = vertcat(s_tot,s(idx_start_sincro_Nirs_first:end,:));
% 
% % Carico i restanti file interi
% for k = (idx_file_Nirs_start_sincro_first+1):(idx_file_Nirs_end_sincro_first-1)
%     tmp = sprintf(dataset_name_first, k);
%     disp(['Load: ',tmp])
%     load(tmp,'-mat');
%     d_tot = vertcat(d_tot,d);
%     t_tot = vertcat(t_tot,t);
%     aux_tot = vertcat(aux_tot,aux);
%     s_tot = vertcat(s_tot,s);
% end
% 
% % Carico l'ultimo file fino all'indice trovato prima
% name_file_Nirs_end_sincro_first = 'PD4 18.16.25 Sun 12 Apr 2020 23.nirs'
% first_end_file_name = name_file_Nirs_end_sincro_first;
% disp(['Load: ',first_end_file_name])
% load(first_end_file_name,'-mat');
% d_tot   = vertcat(d_tot,d(1:idx_end_sincro_Nirs_first,:));
% t_tot   = vertcat(t_tot,t(1:idx_end_sincro_Nirs_first,:));
% aux_tot = vertcat(aux_tot,aux(1:idx_end_sincro_Nirs_first,:));
% s_tot   = vertcat(s_tot,s(1:idx_end_sincro_Nirs_first,:));
% 
% % Save of the NIRS data under analysis
% % (!! Be carefull when uncomment the following code !!)
% 
% % clear d t aux s
% % d = d_tot;
% % t = t_tot;
% % aux = aux_tot;
% % s = s_tot;
% % 
% % fold_path = fullfile(base_path,'CODICE/','DEFINITIVO/','NIRS/','PD4/');
% % mkdir(fold_path,'PD4_NIRS_DATA_EU/');
% % 
% % to_save = fullfile(base_path,'CODICE/','DEFINITIVO/','NIRS/','PD4/','PD4_NIRS_DATA_EU/','PD4_NIRS_first_eu_tot_2.mat');
% % save (to_save, 'd','t','aux','s','SD')
% 
% a = 0;
% %%
% 
% synchro_result_first_eu = struct();
% synchro_result_first_eu.file_start_first_eu = name_file_Nirs_start_sincro_first;
% synchro_result_first_eu.time_start_first_eu_file = time_file_Nirs_start_sincro_first;
% synchro_result_first_eu.idx_start_first_eu = idx_start_sincro_Nirs_first;
% synchro_result_first_eu.file_end_first_eu = name_file_Nirs_end_sincro_first;
% synchro_result_first_eu.time_end_first_eu_file = time_file_Nirs_end_sincro_first;
% synchro_result_first_eu.idx_end_first_eu = idx_end_sincro_Nirs_first;
% 
% synchro_result = struct();
% synchro_result.idx_PD   = idx_PD;
% synchro_result.first_eu = synchro_result_first_eu;
% a = 0;
% %%
% % Save of the results (!!Be carefull when uncomment the next line!!)
% % save PD4_synchro_result_2.mat synchro_result



%% ULTIMO INTERVALLO DI EUGLICEMIA 

% CASO SPECIALE RICOMINCIO DA CAPO

close all
clear all
clc

idx_PD = 25;


%% SETTING OF FOLDER PATH 

base_path      = pwd;
data_path_cgm  = fullfile(base_path,'DATI','CGM_INTERP');
%data_path_nirs = fullfile('E:\DATI\NIRS\PD59\');
data_path_nirs = fullfile(base_path,'DATI','NIRS','PD25');
function_path  = fullfile(base_path,'CODICE','FUNCTION');

addpath(genpath(pwd))
addpath(data_path_cgm)
addpath(function_path)
addpath(data_path_nirs)

%%
% Load of cgm data
PD_struct = load ('PD25.mat');
date_CGM  = PD_struct.PD.PD.Time;
cgm_value = PD_struct.PD.PD.("mg/dL");

% Extraction of the start and finish date of the CGM acquisition
date_CGM_start = date_CGM(1);
date_CGM_end = date_CGM(end);

% Display of the preliminary inspection of the CGM dataset
disp('==============================================================================')
disp('CGM INTERVAL DATE')
disp(['Start : ',datestr(date_CGM_start)])
disp(['End   : ',datestr(date_CGM_end)])

%% 
% Load of the usefull glycemic analysis result
load idx_eu_results.mat

%did this and saved new variable 14 03 25
% idx_eu_int.idx_end_first_eu_w(3) = 916-51;
%idx_eu_int.idx_start_last_eu_w(3) = 929+51;
%idx_eu_int.idx_end_first_eu_w(10) = 351-51;
%idx_eu_int.idx_start_last_eu_w(10) = 374+51;

PD_results = idx_eu_int(idx_eu_int.pax_idx_vector==idx_PD, : );

% Date interval of the first and last Euglycemia window
date_start_first_eu = date_CGM(PD_results.idx_start_first_eu_w);
date_end_first_eu   = date_CGM(PD_results.idx_end_first_eu_w);
date_start_last_eu  = date_CGM(PD_results.idx_start_last_eu_w);
date_end_last_eu    = date_CGM(PD_results.idx_end_last_eu_w);
dur_first_eu        = minutes(date_end_first_eu - date_start_first_eu);
dur_last_eu         = minutes(date_end_last_eu - date_start_last_eu);

% Display of the first and last Euglycemia window dates
disp(' ')
disp('FIRST EUGLYCEMIA INTERVAL')
disp(['Start    : ',datestr(date_start_first_eu)])
disp(['End      : ',datestr(date_end_first_eu)])
disp(['Duration : ',num2str(dur_first_eu),' Min'])
disp(' ')
disp('LAST EUGLYCEMIA INTERVAL')
disp(['Start    : ',datestr(date_start_last_eu)])
disp(['End      : ',datestr(date_end_last_eu)])
disp(['Duration : ',num2str(dur_last_eu),' Min'])

%% NIRS DATASET ANALYSES

% It is analysed as automatically as possible the NIRS dataset in order to
% extract the acquisition information in order to proceed with the CGM-NIRS
% synchronization.

% Count of NIRS files for each dataset.
% It is generated the path name of the folder which contains the NIRS
% files. The command 'dir' extracts the informations of all files which
% contains inside the name the words passed as argument.
tmp = fullfile(data_path_nirs,'*PD25 18.17.*.nirs');
info_first_dataset = dir(tmp);
n_dataset_1 = length(info_first_dataset);

% NIRS sampling period
Ts_NIRS = 0.1; % [s]
% Standard length of a NIRS file
n_standard = 30001;

%% EXTRACTION OF ACQUISITION INTERVAL DATE

% Extraction of the acquisition start date

% Extraction of the file name
file_name_first = info_first_dataset.name;
% Separation of the information contained inside the file name and
% conversion in DateTime format
sep = split(string(file_name_first));
time_first = split(sep(2),'.');
datetime_format = strcat(sep(4),'-',sep(5),'-',sep(6),{' '},time_first(1),':',time_first(2),':',time_first(3));
date_start_first = datetime(datetime_format);

% In order to get the length of the last acqusition file it is necessary to
% load it
load('PD25 18.17.02 Wed 03 Mar 2021 141.nirs','-mat');
n_last = length(t);

% Total number of samples acquired in the first dataset and convertion in
% seconds
n1_complete = (n_dataset_1-1)*n_standard;
n1_tot = n1_complete + n_last;
duration_first = seconds(Ts_NIRS*n1_tot); %[s]

% End date of first dataset acquisition
date_end_first = date_start_first + duration_first;

% Display of the NIRS acqusition interval dates 
disp('==============================================================================')
disp('NIRS INTERVAL DATE')
disp(' ')
disp('NIRS 1 :')
disp(['Start : ',datestr(date_start_first)])
disp(['End   : ',datestr(date_end_first)])

%% CGM AND NIRS SYNCHRONIZATION

% By the use of the 'date_sincro_1_NIRS' it is extracted the intreval date
% of the CGM-NIRS synchronization in case of 1 NIRS dataset.

f_disp = 0; % Display Command Window results
[date_start_sincro,date_end_sincro] = date_sincro_1_NIRS(date_CGM,date_start_first,date_end_first,f_disp);

% The same approach is used to determine whether the start and the end of
% the first and last Euglycemia window is inside the NIRS interval

f_disp = 1; % NO display Command Window results
%CHANGE FOR START AND LAST 4 03 25 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [date_start_sincro_first,date_end_sincro_first] = date_sincro_1_NIRS(date_CGM(PD_results.idx_start_first_eu_w:PD_results.idx_end_first_eu_w),date_start_first,date_end_first,f_disp);
[date_start_sincro_last,date_end_sincro_last]  = date_sincro_1_NIRS(date_CGM(PD_results.idx_start_last_eu_w:PD_results.idx_end_last_eu_w),date_start_first,date_end_first,f_disp);

a = 0;
%% ESTRAZIONE INDICI DI INIZIO E FINE SINCRONIZZAZIONE

% Si verifica dove si colloca la data di inizio sincronizzazione e si
% calcola l'intervallo di tempo che intercorre tra tale data e l'inizio 
% dell'acquisizione NIRS. Si salva il nome del dataset NIRS nel quale cade
% la data di inizio sincronizzazione.

if date_start_sincro < date_end_first 
    diff_Nirs_sincro = date_start_sincro - date_start_first;
    date_start_Nirs = date_start_first;
    dataset_name = 'PD25 18.17.02 Wed 03 Mar 2021 %d.nirs';
end

% Durata in secondi di ogni file NIRS campionato completo (tranne l'ultimo)
dur_Nirs_file_sec = seconds(n_standard*0.1); %[s]

% Si calcola la data di inizio acquisizione per ogni file NIRS del dataset
% in analisi.
date_acq_Nirs_file = [];
for i=0:(n_dataset_1-1)
    date_acq_Nirs_file = [date_acq_Nirs_file ; date_start_Nirs+(i*dur_Nirs_file_sec)];
end

% Si estrae il nome del file NIRS contenente l'istante di inizio
% sincronizzazione
conv = 0;
i = 0;
% Si scorre il file contenente la data di inizio acquisizione per ogni file
% NIRS del dataset contenente la data di inizio sincronizzazione
while conv == 0
    i = i+1;
    % Se la data di inizio acquisizione del file i-esimo è maggiore della
    % data di inizio sincronizzazione allora l'inidice del file che
    % contiene la data ricercata è il precedente 
    if date_acq_Nirs_file(i) > date_start_sincro
        idx_file_Nirs_start_sincro = i - 1;
        conv = 1;
    end
    % Se la data di inizio acquisizione del file i-esimo è uguale alla data
    % di inizio sincronizzazione allora l'inidice del file che contiene la 
    % data ricercata è il corrente
    if date_acq_Nirs_file(i) == date_start_sincro
        idx_file_Nirs_start_sincro = i;
        conv = 1;
    end
    % Se la data di inizio acquisizione del file i-esimo è minore della 
    % data di inizio sincronizzazione e si è arrivati all'ultimo file NIRS
    % allora l'inidice del file che contiene la data ricercata è il 
    % corrente    
    if date_acq_Nirs_file(i) < date_start_sincro & i == length(date_acq_Nirs_file)
        idx_file_Nirs_start_sincro = i;
        conv = 1;
    end
    
    % Una volta trovato l'indice del file contenente il punto di inizio
    % sincronizzazione si crea la variabile che contiene il nome del file
    if conv == 1
        % Data file Nirs per inizio sincro
        time_file_Nirs_start_sincro = date_acq_Nirs_file(idx_file_Nirs_start_sincro);
        % Nome file Nirs per inizio sincro, si sottrae uno dell'indice
        % trovato in precedenza poichè la nomenclatura dei file inizia da
        % zero.
        name_file_Nirs_start_sincro = sprintf(dataset_name, idx_file_Nirs_start_sincro-1);
    end
end

% Si effettua lo stesso procedimento per identificare il nome del file
% contenente la data di fine sincronizzazione
conv = 0;
i = 0;
while conv == 0
    i = i+1;
    if date_acq_Nirs_file(i) > date_end_sincro
        idx_file_Nirs_end_sincro = i - 1;
        conv = 1;
    end
    if date_acq_Nirs_file(i) == date_end_sincro
        idx_file_Nirs_end_sincro = i;
        conv = 1;
    end
    if date_acq_Nirs_file(i) < date_end_sincro & i == length(date_acq_Nirs_file)
        idx_file_Nirs_end_sincro = i;
        conv = 1;
    end
    
    if conv == 1
        time_file_Nirs_end_sincro = date_acq_Nirs_file(idx_file_Nirs_end_sincro);
        name_file_Nirs_end_sincro = sprintf(dataset_name, idx_file_Nirs_end_sincro-1);
    end
end

% Si calcola l'indice di inizio sincronizzazione all'interno del file NIRS 
% contentenete la data di inizio sincronizzazione
diff_start_Nirs = date_start_sincro - time_file_Nirs_start_sincro;
diff_start_Nirs_s = seconds(diff_start_Nirs); %[s]
idx_start_sincro_Nirs = round(diff_start_Nirs_s/Ts_NIRS); %[time point]

% Si calcola l'indice di fine sincronizzazione all'interno del file NIRS 
% contentenete la data di fine sincronizzazione 
diff_end_Nirs = date_end_sincro - time_file_Nirs_end_sincro;
diff_end_Nirs_s = seconds(diff_end_Nirs);
idx_end_sincro_Nirs = round(diff_end_Nirs_s/Ts_NIRS);

% Display informazioni sulla sincronizzazione lato database NIRS
disp('==============================================================================')
disp('NIRS - CGM SYNCHRONIZATION')
disp(' ')
disp(['NIRS FILE INFO FOR SINCRO FROM ',datestr(date_start_sincro),' /TO/ ',datestr(date_end_sincro)])
disp(' ')
disp('START')
disp(['Name   : ',name_file_Nirs_start_sincro])
disp(['Time   : ',datestr(time_file_Nirs_start_sincro)])
disp(['Sample : ',num2str(idx_start_sincro_Nirs)])
disp(' ')
disp('END')
disp(['Name   : ',name_file_Nirs_end_sincro])
disp(['Time   : ',datestr(time_file_Nirs_end_sincro)])
disp(['Sample : ',num2str(idx_end_sincro_Nirs)])

%% ESTRAZIONE INDICI DI INIZIO E FINE SINCRONIZZAZIONE PRIMA EU

% Si verifica dove si colloca la data di inizio sincronizzazione e si
% calcola l'intervallo di tempo che intercorre tra tale data e l'inizio 
% dell'acquisizione NIRS. Si salva il nome del dataset NIRS nel quale cade
% la data di inizio sincronizzazione.
if date_start_sincro_first < date_end_first 
    diff_Nirs_sincro = date_start_sincro_first - date_start_first;
    date_start_Nirs = date_start_first;
dataset_name = 'PD25 18.17.02 Wed 03 Mar 2021 %d.nirs';
end

% Durata in secondi di ogni file NIRS campionato completo (tranne l'ultimo)
dur_Nirs_file_sec = seconds(n_standard*0.1); %[s]

% Si calcola la data di inizio acquisizione per ogni file NIRS del dataset
% in analisi.
date_acq_Nirs_file = [];
for i=0:(n_dataset_1-1)
    date_acq_Nirs_file = [date_acq_Nirs_file ; date_start_Nirs+(i*dur_Nirs_file_sec)];
end

% Si estrae il nome del file NIRS contenente l'istante di inizio
% sincronizzazione
conv = 0;
i = 0;
% Si scorre il file contenente la data di inizio acquisizione per ogni file
% NIRS del dataset contenente la data di inizio sincronizzazione
while conv == 0
    i = i+1;
    % Se la data di inizio acquisizione del file i-esimo è maggiore della
    % data di inizio sincronizzazione allora l'inidice del file che
    % contiene la data ricercata è il precedente 
    if date_acq_Nirs_file(i) > date_start_sincro_first
        idx_file_Nirs_start_sincro = i - 1;
        conv = 1;
    end
    % Se la data di inizio acquisizione del file i-esimo è uguale alla data
    % di inizio sincronizzazione allora l'inidice del file che contiene la 
    % data ricercata è il corrente
    if date_acq_Nirs_file(i) == date_start_sincro_first
        idx_file_Nirs_start_sincro = i;
        conv = 1;
    end
    % Se la data di inizio acquisizione del file i-esimo è minore della 
    % data di inizio sincronizzazione e si è arrivati all'ultimo file NIRS
    % allora l'inidice del file che contiene la data ricercata è il 
    % corrente    
    if date_acq_Nirs_file(i) < date_start_sincro_first & i == length(date_acq_Nirs_file)
        idx_file_Nirs_start_sincro = i;
        conv = 1;
    end

    % Una volta trovato l'indice del file contenente il punto di inizio
    % sincronizzazione si crea la variabile che contiene il nome del file
    if conv == 1
        % Data file Nirs per inizio sincro
        time_file_Nirs_start_sincro_first = date_acq_Nirs_file(idx_file_Nirs_start_sincro);
        % Nome file Nirs per inizio sincro, si sottrae uno dell'indice
        % trovato in precedenza poichè la nomenclatura dei file inizia da
        % zero.
        name_file_Nirs_start_sincro_first = sprintf(dataset_name, idx_file_Nirs_start_sincro-1);
    end
end

% Si effettua lo stesso procedimento per identificare il nome del file
% contenente la data di fine sincronizzazione
conv = 0;
i = 0;
while conv == 0
    i = i+1;
    if date_acq_Nirs_file(i) > date_end_sincro_first
        idx_file_Nirs_end_sincro = i - 1;
        conv = 1;
    end
    if date_acq_Nirs_file(i) == date_end_sincro_first
        idx_file_Nirs_end_sincro = i;
        conv = 1;
    end
    if date_acq_Nirs_file(i) < date_end_sincro_first & i == length(date_acq_Nirs_file)
        idx_file_Nirs_end_sincro = i;
        conv = 1;
    end

    if conv == 1
        time_file_Nirs_end_sincro_first = date_acq_Nirs_file(idx_file_Nirs_end_sincro);
        name_file_Nirs_end_sincro_first = sprintf(dataset_name, idx_file_Nirs_end_sincro-1);
    end
end

% Si calcola l'indice di inizio sincronizzazione all'interno del file NIRS 
% contentenete la data di inizio sincronizzazione che nel PD15 corrisponde 
% alla data del primo campione CGM valido (16-Jul_2020 14:05:00) 
diff_start_Nirs = date_start_sincro_first - time_file_Nirs_start_sincro_first;
diff_start_Nirs_s = seconds(diff_start_Nirs); %[s]
idx_start_sincro_Nirs_first = round(diff_start_Nirs_s/Ts_NIRS); %[time point]

% Si calcola l'indice di fine sincronizzazione all'interno del file NIRS 
% contentenete la data di fine sincronizzazione che nel PD15 corrisponde 
% alla data dell'ultimo campione CGM valido e contenuto contemporaneamente
% in un file NIRS (20-Jul_2020 11:40:00)
diff_end_Nirs = date_end_sincro_first - time_file_Nirs_end_sincro_first;
diff_end_Nirs_s = seconds(diff_end_Nirs);
idx_end_sincro_Nirs_first = round(diff_end_Nirs_s/Ts_NIRS);

% Display informazioni sulla sincronizzazione lato database NIRS
disp('==============================================================================')
disp('FIRST EUGLYCEMIA INTERVAL SYNCHRONIZATION')
disp(' ')
disp(['NIRS FILE INFO FOR SINCRO FROM ',datestr(date_start_sincro_first),' /TO/ ',datestr(date_end_sincro_first)])
disp(' ')
disp('START')
disp(['Name   : ',name_file_Nirs_start_sincro_first])
disp(['Time   : ',datestr(time_file_Nirs_start_sincro_first)])
disp(['Sample : ',num2str(idx_start_sincro_Nirs_first)])
disp(' ')
disp('END')
disp(['Name   : ',name_file_Nirs_end_sincro_first])
disp(['Time   : ',datestr(time_file_Nirs_end_sincro_first)])
disp(['Sample : ',num2str(idx_end_sincro_Nirs_first)])

%% ESTRAZIONE INDICI DI INIZIO E FINE SINCRONIZZAZIONE ULTIMA EU

% Si verifica dove si colloca la data di inizio sincronizzazione e si
% calcola l'intervallo di tempo che intercorre tra tale data e l'inizio 
% dell'acquisizione NIRS. Si salva il nome del dataset NIRS nel quale cade
% la data di inizio sincronizzazione.
if date_start_sincro_last < date_end_first 
    diff_Nirs_sincro = date_start_sincro_last - date_start_first;
    date_start_Nirs = date_start_first;
dataset_name = 'PD25 18.17.02 Wed 03 Mar 2021 %d.nirs';
end

% Durata in secondi di ogni file NIRS campionato completo (tranne l'ultimo)
dur_Nirs_file_sec = seconds(n_standard*0.1); %[s]

% Si calcola la data di inizio acquisizione per ogni file NIRS del dataset
% in analisi.
date_acq_Nirs_file = [];
for i=0:(n_dataset_1-1)
    date_acq_Nirs_file = [date_acq_Nirs_file ; date_start_Nirs+(i*dur_Nirs_file_sec)];
end

% Si estrae il nome del file NIRS contenente l'istante di inizio
% sincronizzazione
conv = 0;
i = 0;
% Si scorre il file contenente la data di inizio acquisizione per ogni file
% NIRS del dataset contenente la data di inizio sincronizzazione
while conv == 0
    i = i+1;
    % Se la data di inizio acquisizione del file i-esimo è maggiore della
    % data di inizio sincronizzazione allora l'inidice del file che
    % contiene la data ricercata è il precedente 
    if date_acq_Nirs_file(i) > date_start_sincro_last
        idx_file_Nirs_start_sincro = i - 1;
        conv = 1;
    end
    % Se la data di inizio acquisizione del file i-esimo è uguale alla data
    % di inizio sincronizzazione allora l'inidice del file che contiene la 
    % data ricercata è il corrente
    if date_acq_Nirs_file(i) == date_start_sincro_last
        idx_file_Nirs_start_sincro = i;
        conv = 1;
    end
    % Se la data di inizio acquisizione del file i-esimo è minore della 
    % data di inizio sincronizzazione e si è arrivati all'ultimo file NIRS
    % allora l'inidice del file che contiene la data ricercata è il 
    % corrente    
    if date_acq_Nirs_file(i) < date_start_sincro_last & i == length(date_acq_Nirs_file)
        idx_file_Nirs_start_sincro = i;
        conv = 1;
    end
    
    % Una volta trovato l'indice del file contenente il punto di inizio
    % sincronizzazione si crea la variabile che contiene il nome del file
    if conv == 1
        % Data file Nirs per inizio sincro
        time_file_Nirs_start_sincro_last = date_acq_Nirs_file(idx_file_Nirs_start_sincro);
        % Nome file Nirs per inizio sincro, si sottrae uno dell'indice
        % trovato in precedenza poichè la nomenclatura dei file inizia da
        % zero.
        name_file_Nirs_start_sincro_last = sprintf(dataset_name, idx_file_Nirs_start_sincro-1);
    end
end

% Si effettua lo stesso procedimento per identificare il nome del file
% contenente la data di fine sincronizzazione
conv = 0;
i = 0;
while conv == 0
    i = i+1;
    if date_acq_Nirs_file(i) > date_end_sincro_last
        idx_file_Nirs_end_sincro = i - 1;
        conv = 1;
    end
    if date_acq_Nirs_file(i) == date_end_sincro_last
        idx_file_Nirs_end_sincro = i;
        conv = 1;
    end
    if date_acq_Nirs_file(i) < date_end_sincro_last & i == length(date_acq_Nirs_file)
        idx_file_Nirs_end_sincro = i;
        conv = 1;
    end
    
    if conv == 1
        time_file_Nirs_end_sincro_last = date_acq_Nirs_file(idx_file_Nirs_end_sincro);
        name_file_Nirs_end_sincro_last = sprintf(dataset_name, idx_file_Nirs_end_sincro-1);
    end
end

% Si calcola l'indice di inizio sincronizzazione all'interno del file NIRS 
% contentenete la data di inizio sincronizzazione che nel PD15 corrisponde 
% alla data del primo campione CGM valido (16-Jul_2020 14:05:00) 
diff_start_Nirs = date_start_sincro_last - time_file_Nirs_start_sincro_last;
diff_start_Nirs_s = seconds(diff_start_Nirs); %[s]
idx_start_sincro_Nirs_last = round(diff_start_Nirs_s/Ts_NIRS); %[time point]

% Si calcola l'indice di fine sincronizzazione all'interno del file NIRS 
% contentenete la data di fine sincronizzazione che nel PD15 corrisponde 
% alla data dell'ultimo campione CGM valido e contenuto contemporaneamente
% in un file NIRS (20-Jul_2020 11:40:00)
diff_end_Nirs = date_end_sincro_last - time_file_Nirs_end_sincro_last;
diff_end_Nirs_s = seconds(diff_end_Nirs);
idx_end_sincro_Nirs_last = round(diff_end_Nirs_s/Ts_NIRS);

% Display informazioni sulla sincronizzazione lato database NIRS
disp('==============================================================================')
disp('LAST EUGLYCEMIA INTERVAL SYNCHRONIZATION')
disp(' ')
disp(['NIRS FILE INFO FOR SINCRO FROM ',datestr(date_start_sincro_last),' /TO/ ',datestr(date_end_sincro_last)])
disp(' ')
disp('START')
disp(['Name   : ',name_file_Nirs_start_sincro_last])
disp(['Time   : ',datestr(time_file_Nirs_start_sincro_last)])
disp(['Sample : ',num2str(idx_start_sincro_Nirs_last)])
disp(' ')
disp('END')
disp(['Name   : ',name_file_Nirs_end_sincro_last])
disp(['Time   : ',datestr(time_file_Nirs_end_sincro_last)])
disp(['Sample : ',num2str(idx_end_sincro_Nirs_last)])

a = 0;
%%
% synchro_result_first_eu = struct();
% synchro_result_first_eu.file_start_first_eu = name_file_Nirs_start_sincro_first;
% synchro_result_first_eu.time_start_first_eu_file = time_file_Nirs_start_sincro_first;
% synchro_result_first_eu.idx_start_first_eu = idx_start_sincro_Nirs_first;
% synchro_result_first_eu.file_end_first_eu = name_file_Nirs_end_sincro_first;
% synchro_result_first_eu.time_end_first_eu_file = time_file_Nirs_end_sincro_first;
% synchro_result_first_eu.idx_end_first_eu = idx_end_sincro_Nirs_first;

synchro_result_last_eu = struct();
synchro_result_last_eu.file_start_last_eu = name_file_Nirs_start_sincro_last;
synchro_result_last_eu.time_start_last_eu_file = time_file_Nirs_start_sincro_last;
synchro_result_last_eu.idx_start_last_eu = idx_start_sincro_Nirs_last;
synchro_result_last_eu.file_end_last_eu = name_file_Nirs_end_sincro_last;
synchro_result_last_eu.time_end_last_eu_file = time_file_Nirs_end_sincro_last;
synchro_result_last_eu.idx_end_last_eu = idx_end_sincro_Nirs_last;

synchro_result = struct();
synchro_result.idx_PD   = idx_PD;


%synchro_result.first_eu = synchro_result_first_eu;
synchro_result.last_eu  = synchro_result_last_eu;

% Save of the results (!!Be carefull when uncomment the next line!!)

%save PD25_synchro_result_first_eu.mat synchro_result
save PD25_synchro_result_last_eu.mat synchro_result


%% LOAD AND SAVE OF REQUIRED NIRS DATA
% CGM - NIRS SYNCHRONIZATION INTERVAL DATE
% 
% Start : 10-Apr-2020 19:05:00
% End   : 10-Apr-2020 23:25:00
% FIRST EUGLYCEMIA INTERVAL SYNCHRONIZATION
% % 
% NIRS FILE INFO FOR SINCRO FROM 04-Mar-2021 19:12:00 /TO/ 04-Mar-2021 21:42:00
% 
% START
% Name   : PD25 18.17.02 Wed 03 Mar 2021 29.nirs
% Time   : 04-Mar-2021 18:27:04
% Sample : 26951
% 
% END
% Name   : PD25 18.17.02 Wed 03 Mar 2021 32.nirs
% Time   : 04-Mar-2021 20:57:05
% Sample : 26948

% LAST EUGLYCEMIA INTERVAL SYNCHRONIZATION
% 
% NIRS FILE INFO FOR SINCRO FROM 08-Mar-2021 08:42:00 /TO/ 08-Mar-2021 15:22:00
% 
% START
% Name   : PD25 18.17.02 Wed 03 Mar 2021 132.nirs
% Time   : 08-Mar-2021 08:17:15
% Sample : 14848
% 
% END
% Name   : PD25 18.17.02 Wed 03 Mar 2021 140.nirs
% Time   : 08-Mar-2021 14:57:16
% Sample : 14840

% Da cambiare per altri PAX
% SPECIFIC FOR PD10
idx_file_Nirs_start_sincro_first = 29;
idx_file_Nirs_end_sincro_first = 32;

idx_file_Nirs_start_sincro_last = 132;
idx_file_Nirs_end_sincro_last = 140;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% change for first and last
% % First EU
if date_start_first_eu < date_end_first 
    dataset_name = 'PD25 18.17.02 Wed 03 Mar 2021 %d.nirs';
    dataset_name_first = dataset_name;
end

% Last Eu
 if date_start_last_eu < date_end_first 
     dataset_name_last = 'PD25 18.17.02 Wed 03 Mar 2021 %d.nirs';
 end


% In teoria già sopra mi fa il controllo se per caso la fine della finestra
% cade nell'ingervallo di tempo dove non ho dati NIRS


%% FIRST EUGLYCEMIA INTERVAL

% Load of the NIRS file for the first Euglycemia interval
disp('===================================================================')
disp('Load NIRS file first Euglycemia interval')
disp(' ')

% Inizializzo 
d_tot = [];
t_tot = [];
aux_tot = [];
s_tot = [];

% % Carico il primo file
first_start_file_name = name_file_Nirs_start_sincro_first;
disp(['Load: ',first_start_file_name])
load(first_start_file_name,'-mat');
d_tot   = vertcat(d_tot,d(idx_start_sincro_Nirs_first:end,:));
t_tot   = vertcat(t_tot,t(idx_start_sincro_Nirs_first:end,:));
aux_tot = vertcat(aux_tot,aux(idx_start_sincro_Nirs_first:end,:));
s_tot   = vertcat(s_tot,s(idx_start_sincro_Nirs_first:end,:));

% Carico i restanti file interi
for k = (idx_file_Nirs_start_sincro_first+1):(idx_file_Nirs_end_sincro_first-1)
    tmp = sprintf(dataset_name_first, k);
    disp(['Load: ',tmp])
    load(tmp,'-mat');
    d_tot = vertcat(d_tot,d);
    t_tot = vertcat(t_tot,t);
    aux_tot = vertcat(aux_tot,aux);
    s_tot = vertcat(s_tot,s);
end

% Carico l'ultimo file fino all'indice trovato prima
first_end_file_name = name_file_Nirs_end_sincro_first;
disp(['Load: ',first_end_file_name])
load(first_end_file_name,'-mat');
d_tot   = vertcat(d_tot,d(1:idx_end_sincro_Nirs_first,:));
t_tot   = vertcat(t_tot,t(1:idx_end_sincro_Nirs_first,:));
aux_tot = vertcat(aux_tot,aux(1:idx_end_sincro_Nirs_first,:));
s_tot   = vertcat(s_tot,s(1:idx_end_sincro_Nirs_first,:));

%% Save of the NIRS data under analysis
% (!! Be carefull when uncomment the following code !!)

clear d t aux s
d = d_tot;
t = t_tot;
aux = aux_tot;
s = s_tot;

first_path = fullfile(base_path,'CODICE/','DEFINITIVO/','NIRS/','PD25/','PD25_NIRS_DATA_EU/','PD25_NIRS_first_eu_tot.mat');
save (first_path, 'd','t','aux','s','SD')


% 
% fold_path = fullfile(base_path,'CODICE/','DEFINITIVO/','NIRS/','PD19/');
% mkdir(fold_path,'PD19_NIRS_DATA_EU/');
% 
% to_save = fullfile(base_path,'CODICE/','DEFINITIVO/','NIRS/','PD19/','PD19_NIRS_DATA_EU/','PD19_NIRS_first_eu_tot.mat');
% save (to_save, 'd','t','aux','s','SD')

%% LAST EU INTERVAL

% Load of the NIRS file for the Last Euglycemia interval
disp('===================================================================')
disp('Load NIRS file last Euglycemia interval')
disp(' ')

% Inizializzo 
d_tot = [];
t_tot = [];
aux_tot = [];
s_tot = [];

% % Carico il primo file
last_start_file_name = name_file_Nirs_start_sincro_last;
disp(['Load: ',last_start_file_name])
load(last_start_file_name,'-mat');
d_tot   = vertcat(d_tot,d(idx_start_sincro_Nirs_last:end,:));
t_tot   = vertcat(t_tot,t(idx_start_sincro_Nirs_last:end,:));
aux_tot = vertcat(aux_tot,aux(idx_start_sincro_Nirs_last:end,:));
s_tot   = vertcat(s_tot,s(idx_start_sincro_Nirs_last:end,:));

% Carico i restanti file interi
for k = (idx_file_Nirs_start_sincro_last+1):(idx_file_Nirs_end_sincro_last-1)
    tmp = sprintf(dataset_name_last, k);
    disp(['Load: ',tmp])
    load(tmp,'-mat');
    d_tot = vertcat(d_tot,d);
    t_tot = vertcat(t_tot,t);
    aux_tot = vertcat(aux_tot,aux);
    s_tot = vertcat(s_tot,s);
end

% Carico l'ultimo file fino all'indice trovato prima
last_end_file_name = name_file_Nirs_end_sincro_last;
disp(['Load: ',last_end_file_name])
load(last_end_file_name,'-mat');
d_tot   = vertcat(d_tot,d(1:idx_end_sincro_Nirs_last,:));
t_tot   = vertcat(t_tot,t(1:idx_end_sincro_Nirs_last,:));
aux_tot = vertcat(aux_tot,aux(1:idx_end_sincro_Nirs_last,:));
s_tot   = vertcat(s_tot,s(1:idx_end_sincro_Nirs_last,:));

%% SAVE LAST EU TOT MAT
% Save of the NIRS data under analysis
% (!! Be carefull when uncomment the following code !!)

clear d t aux s
d = d_tot;
t = t_tot;
aux = aux_tot;
s = s_tot;

last_path = fullfile(base_path,'CODICE/','DEFINITIVO/','NIRS/','PD25/','PD25_NIRS_DATA_EU/','PD25_NIRS_last_eu_tot.mat');
save (last_path, 'd','t','aux','s','SD')






