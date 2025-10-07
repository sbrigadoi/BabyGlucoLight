% PATIENT SELECTION

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
pax_idx_vector = [1:1:17 19:1:20 22:1:26 28 33:1:45 47:1:55]';
% pax_idx_vector = [4;5;8;9;11;14;15;19;47];
% pax_idx_vector = 4;
n_pax          = length(pax_idx_vector);

% Glycemic Analysis: Display Setting
result_disp       = 0;  % 1 = show results on command window for each PD
result_plot       = 0;  % 1 = show plot R2021  3 = show plot R2023 for each PD
threshold_25_perc = 45; % [Min] Critical Euglycemia Gap threshold

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

n_hypo_mild = [];
n_hypo_sev  = [];
n_hyper_mild = [];
n_hyper_sev = [];
n_hypo_hyper = [];
n_tot_eventi = [];
n_event_removed = [];
f_start_in = [];
f_end_in = [];
n_critici = [];
n_inizio_fine_evento = [];
n_eventi_buoni = [];
dur_hypo_m = [];
dur_hypo_s = [];
dur_hyper_m = [];
dur_hyper_s = [];
dur_hypo_hyper = [];
TOR = [];

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

    %% TESIIIIIII

    % ESTRAGGO SOLO QUELLI NON CRITICI !
    idx_non_critic = find(events_analysis_results.glycemic_events.f_critic_event==0);
    % Classifico solo i non critici
    n_hypo_mild_tmp  = length(find(events_analysis_results.glycemic_events.type_event(idx_non_critic) == 'Mild Hypo'));
    n_hypo_sev_tmp   = length(find(events_analysis_results.glycemic_events.type_event(idx_non_critic) == 'Severe Hypo'));
    n_hyper_mild_tmp = length(find(events_analysis_results.glycemic_events.type_event(idx_non_critic) == 'Mild Hyper'));
    n_hyper_sev_tmp  = length(find(events_analysis_results.glycemic_events.type_event(idx_non_critic) == 'Severe Hyper'));
    n_hypo_hyper_tmp = length(find(events_analysis_results.glycemic_events.type_event(idx_non_critic) == 'Hypo / Hyper'));

    % SALVO RISULTATI
    tmp_n_inizio_fine_ev = events_analysis_results.event_pax_report.f_start_inside_ev + events_analysis_results.event_pax_report.f_end_inside_ev;

    n_critici = [n_critici ; length(find(events_analysis_results.glycemic_events.f_critic_event==1))];
    n_critici_inizio_fine = events_analysis_results.event_pax_report.n_events-length(idx_non_critic);
    n_tot_eventi = [n_tot_eventi ; events_analysis_results.event_pax_report.n_events];
    n_inizio_fine_evento = [n_inizio_fine_evento ; tmp_n_inizio_fine_ev];
    n_event_removed = [n_event_removed ; n_critici_inizio_fine];

    n_hypo_mild = [n_hypo_mild ; n_hypo_mild_tmp];
    n_hypo_sev  = [n_hypo_sev ; n_hypo_sev_tmp];
    n_hyper_mild = [n_hyper_mild ; n_hyper_mild_tmp];
    n_hyper_sev = [n_hyper_sev ; n_hyper_sev_tmp];
    n_hypo_hyper = [n_hypo_hyper ; n_hypo_hyper_tmp];

    f_start_in = [f_start_in ; events_analysis_results.event_pax_report.f_start_inside_ev];
    f_end_in = [f_end_in ; events_analysis_results.event_pax_report.f_end_inside_ev];
    n_eventi_buoni = [n_eventi_buoni ; length(idx_non_critic)];

    % DURATE 
    tmp_TOR = events_analysis_results.glycemic_events.duration_event(idx_non_critic);
    TOR = [TOR ; sum(tmp_TOR)];

    idx_hypo_m = find(events_analysis_results.glycemic_events.type_event == 'Mild Hypo' & events_analysis_results.glycemic_events.f_critic_event==0);
    dur_hypo_m = [dur_hypo_m ; sum(events_analysis_results.glycemic_events.duration_event(idx_hypo_m))];

    idx_hypo_s = find(events_analysis_results.glycemic_events.type_event == 'Severe Hypo'& events_analysis_results.glycemic_events.f_critic_event==0);
    dur_hypo_s = [dur_hypo_s ; sum(events_analysis_results.glycemic_events.duration_event(idx_hypo_s))];

    idx_hyper_m = find(events_analysis_results.glycemic_events.type_event == 'Mild Hyper'& events_analysis_results.glycemic_events.f_critic_event==0);
    dur_hyper_m = [dur_hyper_m ; sum(events_analysis_results.glycemic_events.duration_event(idx_hyper_m))];

    idx_hyper_s = find(events_analysis_results.glycemic_events.type_event == 'Severe Hyper'& events_analysis_results.glycemic_events.f_critic_event==0);
    dur_hyper_s = [dur_hyper_s ; sum(events_analysis_results.glycemic_events.duration_event(idx_hyper_s))];
    
    idx_hypo_hyper = find(events_analysis_results.glycemic_events.type_event== 'Hypo / Hyper'& events_analysis_results.glycemic_events.f_critic_event==0);
    dur_hypo_hyper = [dur_hypo_hyper ; sum(events_analysis_results.glycemic_events.duration_event(idx_hypo_hyper))];

    %% VECCHIE !!!!
    % SAVE OF THE CGM SIGNAL INTERPOLATED
%     PD.PD.("mg/dL") = events_analysis_results.cgm_value(:);
%     save_PD_path_name   = fullfile(base_path,'DATI','CGM_INTERP',to_load);
%     save (save_PD_path_name,'PD')
    
    % PATIENT SELECTION METRICS
    
    % Duration of total CGM acquisition
    % tmp_dur = minutes(PD.PD.Time(end)-PD.PD.Time(1));
    % dur_CGM = [dur_CGM ; tmp_dur];
    % 
    % % Time elapsed between the start of the acquisition and the first 
    % % Euglycemia interval
    % time_before_first_eu = [time_before_first_eu ; minutes(eu_analysis_results.euglycemic_int.date_start_eu(1) - PD.PD.Time(1))];
    % 
    % % saveas(gcf,['PD' num2str(idx_pax) '.fig']);
    % 
    % % First and last Eu interval indexes
    % idx_start_first_eu_w = [idx_start_first_eu_w ; eu_analysis_results.idx_eu.idx_start_eu(1)];
    % idx_end_first_eu_w   = [idx_end_first_eu_w ; eu_analysis_results.idx_eu.idx_end_eu(1)];
    % idx_start_last_eu_w  = [idx_start_last_eu_w ; eu_analysis_results.idx_eu.idx_start_eu(end)];
    % idx_end_last_eu_w    = [idx_end_last_eu_w ; eu_analysis_results.idx_eu.idx_end_eu(end)];
    % 
    % % Duration of the first and last Euglycemia interval
    % TIR_first = [TIR_first ; eu_analysis_results.euglycemic_int.duration_eu_int(1)];
    % TIR_last  = [TIR_last  ; eu_analysis_results.euglycemic_int.duration_eu_int(end)];
    % 
    % 'Total time in Eu range' (TIR) between the first and last Euglycemia 
    % interval (limits included)
    TIR_eu = [TIR_eu ; sum(eu_analysis_results.euglycemic_int.duration_eu_int)];
    % 
    % % 'Total time outside Eu range' (TOR) between the first and last
    % % Euglycemia interval(limits included). Is the sum of the duration of
    % % all events.
    % idx_end_first_eu  = eu_analysis_results.idx_eu.idx_end_eu(1);
    % idx_start_last_eu = eu_analysis_results.idx_eu.idx_end_eu(end);
    % 
    % % Index of the first event after the first Euglycemia interval
    % for i = 1:1:length(events_analysis_results.idx_events.idx_start_event)
    %     if events_analysis_results.idx_events.idx_start_event(i)>idx_end_first_eu
    %         n_first_event = i;
    %         break
    %     end
    % end
    % 
    % % Index of the last event before the last Euglycemia interval
    % for i = 1:1:length(events_analysis_results.idx_events.idx_start_event)
    %     if events_analysis_results.idx_events.idx_start_event(i)>idx_start_last_eu
    %         n_last_event = i-1;
    %         break
    %     else
    %         n_last_event = i;
    %     end
    % end
    % 
    % TOR_eu = [TOR_eu ; sum(events_analysis_results.glycemic_events.duration_event(n_first_event:n_last_event))];
    % 
    % % TAR = Time Above Range (Hyper)
    % idx_TAR = find(events_analysis_results.glycemic_events.f_critic_event == 0 & (events_analysis_results.glycemic_events.type_event == 'Mild Hyper' | events_analysis_results.glycemic_events.type_event == 'Severe Hyper'));  
    % TAR = [TAR ; sum(events_analysis_results.glycemic_events.duration_event(idx_TAR))];
    % 
    % % TBR = Time Below Range (Hypo)
    % idx_TBR = find(events_analysis_results.glycemic_events.f_critic_event == 0 & (events_analysis_results.glycemic_events.type_event == 'Mild Hypo' | events_analysis_results.glycemic_events.type_event == 'Severe Hypo'));  
    % TBR = [TBR ; sum(events_analysis_results.glycemic_events.duration_event(idx_TBR))];
    % 
    % % TIR and TOR as a percentage
    % tmp_duration_cgm = minutes(PD.PD.Time(eu_analysis_results.idx_eu.idx_end_eu(end))- PD.PD.Time(eu_analysis_results.idx_eu.idx_start_eu(1)));
    % duration_cgm = [duration_cgm ; tmp_duration_cgm];
    % TIR_eu_perc  = [TIR_eu_perc ; round((TIR_eu(end)/tmp_duration_cgm)*100)];
    % TOR_eu_perc  = [TOR_eu_perc ; round((TOR_eu(end)/tmp_duration_cgm)*100)];
    % 
    % 
    % % Number of critical events between the start of the cgm acquisition 
    % % and the end of the last Euglycemia interval
    % tmp_critic_ev    = length(find(events_analysis_results.glycemic_events.f_critic_event(1:n_last_event)==1));
    % n_critical_event = [n_critical_event ; tmp_critic_ev];
    % 
    % % Number of critical gap inside Euglycemia interval
    % tmp_critic_eu = eu_analysis_results.eu_pax_report.n_critical_gap_removed;
    % n_critical_eu = [n_critical_eu ; tmp_critic_eu];
    % 
    % % Total duration of critical gap (defined as number of missing samples 
    % % times 5') between the start of the acquisition and the end of the
    % % last Euglycemia interval.
    % 
    % info_critic_gap     = eu_analysis_results.idx_critical_gap;    
    % idx_end_analysis    = eu_analysis_results.idx_eu.idx_end_eu(end);
    % 
    % idx_critic_gap      = find(eu_analysis_results.idx_critical_gap.idx_start_critic_gap<idx_end_analysis);
    % tmp_n_critic_gap    = length(idx_critic_gap);
    % tmp_duration_cr_gap = ((info_critic_gap.idx_end_critic_gap(idx_critic_gap)-info_critic_gap.idx_start_critic_gap(idx_critic_gap))+1)*5;    
    % tmp_tot_dur_cr_gap  = sum(tmp_duration_cr_gap);
    % 
    % % Gap duration percentage
    % dur_gap      = [dur_gap ; tmp_tot_dur_cr_gap];
    % dur_gap_perc = [dur_gap_perc ; round((tmp_tot_dur_cr_gap/tmp_duration_cgm)*100)];
    % 
    % % Median, Max and Min critic gap duration
    % if length(tmp_duration_cr_gap)>0
    %     tmp_median_dur_cr_g = median(tmp_duration_cr_gap);
    %     tmp_max = max(tmp_duration_cr_gap);
    %     tmp_min = min(tmp_duration_cr_gap);
    % else
    %     tmp_median_dur_cr_g = 0;
    %     tmp_max = 0;
    %     tmp_min = 0;
    % end
    % 
    % duration_cr_gap = [duration_cr_gap ; tmp_tot_dur_cr_gap];
    % median_t_cr_gap = [median_t_cr_gap ; tmp_median_dur_cr_g];
    % max_t_cr_gap    = [max_t_cr_gap ; tmp_max];
    % min_t_cr_gap    = [min_t_cr_gap ; tmp_min];
    % n_critic_gap    = [n_critic_gap ; tmp_n_critic_gap];
    % 
    % % Flag the last Euglycemia interval if it end right before the start of
    % % a critic gap
    % if length(find(info_critic_gap.idx_start_critic_gap == (idx_end_analysis+1)))==1
    %     f_last_eu_gap = [f_last_eu_gap ; 1];
    % else
    %     f_last_eu_gap = [f_last_eu_gap ; 0];
    % end
    % 
    % % Flag the first Euglycemia interval if it begin right before the start
    % % of a critic gap
    % if length(find(info_critic_gap.idx_end_critic_gap == (eu_analysis_results.idx_eu.idx_start_eu(1)-1)))==1
    %     f_first_eu_gap = [f_first_eu_gap ; 1];
    % else
    %     f_first_eu_gap = [f_first_eu_gap ; 0];
    % end

end

%% TABELLA RISULTATI TESI
tesi = table(pax_idx_vector,n_tot_eventi,n_event_removed,n_hypo_mild,n_hypo_sev,n_hyper_mild,n_hyper_sev,n_hypo_hyper,f_start_in,f_end_in,n_critici);

% save risultati_cgm_analisi.mat tesi
% writetable(tesi,'risultati_cgm_analisi.xlsx')

a = 0;
%%
X = ["Ipoglicemia Moderata" "Ipoglicemia Severa" "Iperglicemia Moderata" "Iperglicemia Severa" "Ipo/Iper"];
tot_hypo_m = sum(n_hypo_mild);
std_hypo_m = std(n_hypo_mild);
tot_hypo_s = sum(n_hypo_sev);
std_hypo_s = std(n_hypo_sev);
tot_hyper_m = sum(n_hyper_mild);
std_hyper_m = std(n_hyper_mild);
tot_hyper_s = sum(n_hyper_sev);
std_hyper_s = std(n_hyper_sev);
tot_hypo_hyper = sum(n_hypo_hyper);
std_hypo_hyper = std(n_hypo_hyper);

Y = [tot_hypo_m,tot_hypo_s,tot_hyper_m,tot_hyper_s,tot_hypo_hyper];
std_Y = [std_hypo_m,std_hypo_s,std_hyper_m,std_hyper_s,std_hypo_hyper];

% BARPLOT NUMERO DI EVENTI

figure()

set(gcf, 'Position', get(0, 'Screensize'));

b = bar(Y,'FaceColor','flat')

% Colori
b.CData(1,:) = [0.3010 0.7450 0.9330];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.3010 0.7450 0.9330];
b.CData(4,:) = [0 0.4470 0.7410];
b.CData(5,:) = [0.8500 0.3250 0.0980];

% Numero sopra bar
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center','fontweight','bold'); 

% Grassetto asse x e asse y
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(gca,'FontWeight','bold')
% Label asse x e y
xticklabels(X)
ylabel('Numero di eventi')
% Rimuovo tick asse x
h=gca; h.XAxis.TickLength = [0 0];
% Limite y
ylim([0 200])
%Titolo
title('NUMERO TOTALE DI EVENTI RILEVATI')

% saveas(gcf,'BARPLOT_NUMERO_EVENTI.jpg')

% % print('FIG_42_BAR_NUM_EV','-djpeg','-r600')


%%
X = ["Ipoglicemia" "Ipoglicemia"];

n_tot_hypo = n_hypo_mild + n_hypo_sev;
n_tot_hyper = n_hyper_mild + n_hyper_sev;

hyper = [n_hyper_mild, n_hyper_sev];
hypo = [n_hypo_mild, n_hypo_sev];
hypo_hyper_m = [n_hypo_mild, n_hyper_mild];
hypo_hyper_s = [n_hypo_sev, n_hyper_sev];

% BOXPLOT HYPO VS HYPER

figure()

set(gcf, 'Position', get(0, 'Screensize'));

subplot(2,2,1)
boxplot(hypo,'Labels',{'IPOGLICEMIA MOD','IPOGLICEMIA SEV'})
title('DISTRIBUZIONE IPOGLICEMIA MODERATA VS SEVERA')
set(gca, 'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Numero Eventi')
ylim([-3 35])

subplot(2,2,2)
boxplot(hyper,'Labels',{'IPERGLICEMIA MOD','IPERGLICEMIA SEV'})
title('DISTRIBUZIONE IPERGLICEMIA MODERATA VS SEVERA')
set(gca, 'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Numero Eventi')
ylim([-3 20])

subplot(2,2,3)
boxplot(hypo_hyper_m,'Labels',{'IPOGLICEMIA MOD','IPERGLICEMIA MOD'})
title('DISTRIBUZIONE IPO VS IPER MODERATA')
set(gca, 'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Numero Eventi')
ylim([-3 35])

subplot(2,2,4)
boxplot(hypo_hyper_s,'Labels',{'IPOGLICEMIA SEV','IPERGLICEMIA SEV'})
title('DISTRIBUZIONE IPO VS IPER SEVERA')
set(gca, 'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Numero Eventi')
ylim([-3 7])

% saveas(gcf,'BOXPLOT_hypo_hyper.jpg')
% print(['FIG_43_BOX_PLOT_NUM_EV'],'-djpeg','-r600')


% % Grassetto asse x e asse y
% set(get(gca, 'XAxis'), 'FontWeight', 'bold');
% set(gca,'FontWeight','bold')
% % Label asse x e y
% % xticklabels(X)
% % ylabel('Numero di eventi')
% % Rimuovo tick asse x
% h=gca; h.XAxis.TickLength = [0 0];
% % Limite y
% % ylim([0 200])
% %Titolo
% title('NUMERO TOTALE DI EVENTI RILEVATI')

%% WILCOXON TEST
% clc
% [p_1,h_1] = ranksum(n_hypo_mild,n_hypo_sev)
% [p_2,h_2] = ranksum(n_hyper_mild,n_hyper_sev)
% [p_3,h_3] = ranksum(n_hypo_mild,n_hyper_mild)
% [p_4,h_4] = ranksum(n_hypo_sev,n_hyper_sev);
[p_ms_ipo,h_ms_ipo] = ranksum(n_hypo_mild,n_hypo_sev)
[p_ms_iper,h_ms_iper] = ranksum(n_hyper_mild,n_hyper_sev)
[p_mm_ipo_iper,h_mm_ipo_iper] = ranksum(n_hypo_mild,n_hyper_mild)
[p_ss_ipo_iper,h_ss_ipo_iper] = ranksum(n_hypo_sev,n_hyper_sev)

%%

%%
clear Y
X = [399,62]
Y = [sum(n_eventi_buoni),sum(n_critici) sum(n_inizio_fine_evento)];
labels = {'Eventi completi (337)','Eventi critici rimossi (31)','Eventi etichettati rimossi (31)'};
xlab =["TOTALE RILEVATI" "RIMOSSI"];
explode = [0 1 1];

figure()

set(gcf, 'Position', get(0, 'Screensize'));

subplot(1,2,1)
bb = bar(X,'FaceColor','flat');
bb.CData(2,:) = [0.8500 0.3250 0.0980];
text(1:length(X),X,num2str(X'),'vert','bottom','horiz','center','fontweight','bold');
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(gca,'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];

ylim([0 430])
ylabel('Numero di eventi')
xticklabels(xlab)
title('EVENTI TOTALI vs EVENTI RIMOSSI')

subplot(1,2,2)
pp = pie(Y,explode,labels)
xlm = xlim;
title('COMPOSIZIONE EVENTI RILEVATI','Position',[0 1.5])
set(findobj(pp,'type','text'),'FontWeight','bold') 
% print('FIG_41_PIE','-djpeg','-r600')

% saveas(gcf,'BAR_PIE_REMOVED.jpg')


%% TOR
duration_ev = table(TOR,dur_hypo_m,dur_hypo_s,dur_hyper_m,dur_hyper_s,dur_hypo_hyper);
%%
TAR = dur_hyper_m + dur_hyper_s;
TBR = dur_hypo_m + dur_hypo_s;
a = 0;
%%

X = ["Ipoglicemia Moderata" "Ipoglicemia Severa" "Iperglicemia Moderata" "Iperglicemia Severa" "Ipo/Iper"]
Y = [sum(dur_hypo_m) sum(dur_hypo_s) sum(dur_hyper_m) sum(dur_hyper_s) sum(dur_hypo_hyper)];

figure()
set(gcf, 'Position', get(0, 'Screensize'));
b = bar(Y,'FaceColor','flat')
b.CData(1,:) = [0.3010 0.7450 0.9330];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.3010 0.7450 0.9330];
b.CData(4,:) = [0 0.4470 0.7410];
b.CData(5,:) = [0.8500 0.3250 0.0980];
xticklabels(X)
ylabel('Tempo [Min]')
ylim([0 19000])
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center','fontweight','bold');
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(gca,'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];
title('TEMPO TOTALE PER EVENTO GLICEMICO')
% print('FIG_44_BAR_TEMPO','-djpeg','-r600')

% saveas(gcf,'TEMPO_TOTALE_EVENTI_GLICEMICI.jpg')
%%
figure()
set(gcf, 'Position', get(0, 'Screensize'));
boxplot([TBR,TAR],'Labels',{'IPOGLICEMIA','IPERGLICEMIA'})
title('TEMPO SPESO IN IPOGLICEMIA VS IPERGLICEMIA')
set(gca, 'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Tempo [Min]')
ylim([-100 4750])
% print('FIG_45_BOX_TEMPO','-djpeg','-r600')

% saveas(gcf,'BOXPLOT_TBR_TAR.jpg')

%%
[p,h,stats] = ranksum(TBR,TAR)

%%
figure()
subplot(1,2,2)
set(gcf, 'Position', get(0, 'Screensize'));
boxplot([TOR,TIR_eu],'Labels',{'EVENTO','EUGLICEMIA'})
title('TEMPO SPESO IN UN EVENTO VS EUGLICEMIA')
set(gca, 'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Tempo [Min]')
ylim([-100 7500])


subplot(1,2,1)
set(gcf, 'Position', get(0, 'Screensize'));
boxplot([TBR,TAR],'Labels',{'IPOGLICEMIA','IPERGLICEMIA'})
title('TEMPO SPESO IN IPOGLICEMIA VS IPERGLICEMIA')
set(gca, 'FontWeight','bold')
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Tempo [Min]')
ylim([-100 4500])
% saveas(gcf,'BOXPLOT_TBR_TAR_TIR_TOR.jpg')
print('FIG_45_BOX_TEMPO','-djpeg','-r600')

%%
[p,h,stats] = ranksum(TOR,TIR_eu)

%%
a = table(TOR,TIR_eu);
%%
% er = errorbar(x,data,errlow,errhigh);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 
% er = errorbar(Y,std_Y,'LineWidth',2);
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  

%%
% clc
% [p,h] = ranksum(n_hypo_mild,n_hyper_mild)


%% RESULTS TABLE

pax_analysis = table(pax_idx_vector,time_before_first_eu,TIR_first,TIR_last,duration_cgm,TIR_eu,TIR_eu_perc,TOR_eu,TOR_eu_perc,TAR,TBR,dur_gap,dur_gap_perc);
gap_analysis = table(pax_idx_vector,n_critic_gap,f_first_eu_gap,f_last_eu_gap,duration_cr_gap,median_t_cr_gap,min_t_cr_gap,max_t_cr_gap);

sincro_results = table(pax_idx_vector,idx_start_first_eu_w,idx_end_first_eu_w,idx_start_last_eu_w,idx_end_last_eu_w);

tot_CGM_duration_results = table(pax_idx_vector,dur_CGM);

%% MAX and MIN duration of CGM signal (Congress)

min_cgm = min(tot_CGM_duration_results.dur_CGM);
idx_min_cgm = pax_idx_vector(find(tot_CGM_duration_results.dur_CGM == min_cgm));
min_cgm_tab = table(idx_min_cgm,min_cgm);
disp(min_cgm_tab);

max_cgm = max(tot_CGM_duration_results.dur_CGM);
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
% title(['Time Outside Euglycemia Range [%] , 25-th percentile = ',num2str(th_25_perc_TOR),'%'])
title(['Threshold for the third criterion of patient selection (Q1 of TOR%): ',num2str(th_25_perc_TOR),'%'])
ylabel('Time [%]')
xticklabels('Overall Time Outside euglycemia Range [TOR%]')

pax_analysis_after_TOR = pax_analysis_no_gap(pax_analysis_no_gap.TOR_eu_perc(:)>=th_25_perc_TOR,:);

pax_selected = pax_analysis_after_TOR.pax_idx_vector;

disp('==============================================================================')
disp('PATIENTS SELECTED')
disp(num2str(pax_selected'))
%% TOTAL TIME INSIDE EUGLYCEMIA SELECTION (NOT CONSIDERED AT THE MOMENT)

% th_75_perc_TIR = quantile(pax_analysis_after_TOR.TIR_eu_perc(:),0.75)
% 
% figure()
% boxplot(pax_analysis_after_TOR.TIR_eu_perc(:))
% title(['Time Inside Euglycemia Range [%] , 75-th percentile = ',num2str(th_75_perc_TIR),'%'])
% 
% pax_analysis_after_TIR = pax_analysis_after_TOR(pax_analysis_after_TOR.TIR_eu_perc(:)<=th_75_perc_TIR,:);
% 
% pax_tenuti = pax_analysis_after_TIR.pax_idx_vector

%%
% save pax_analysis_results.mat pax_analysis
% save pax_event_report.mat pax_event_report