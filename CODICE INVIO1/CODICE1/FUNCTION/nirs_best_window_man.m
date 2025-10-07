function [d_best,t_best,C_thresh,tab_win_after_good,best_five] = nirs_best_window_man(synchro_result,d,t,rmch,spect_param,J_850,mesh_850,f_disp,disp_single,disp_result,idx_man)
% Return the best 5' time window 

% INPUT: idx_PD --> Index of patient under analysis
%        synchro_result --> info of current euglycemia interval
%        d --> All raw NIRS data
%        t --> Time vector for raw NIRS data
%        rmch --> Cell array with: dRange ; SNRrange
%        spect_param --> Cell array with: th_power_1 ; th_power_2 ; f_inf_1
%                        f_sup_1 ; f_inf_2 ; f_sup_2 (parameters for
%                        spectral analysis)
%        J_850 --> Jacobian matrix for wavelength of 850 nm
%        mesh_850 --> Home-made head mesh used for Jacobian computation
%        f_disp --> Display results
%        disp_single --> Plot for each channel at each wavelength the
%                        spectral analysis (1 = ON)
%        disp_result --> Plot spectral analysis final results (1 = ON)
%        idx_man --> Manual best window

% OUTPUT: d_best   --> NIRS signal of the best 5' time window
%         t_best   --> time vector of bet 5' time window  
%         c_thresh --> Coverage threshold
%         tab_win_after_good --> info of all possible best windows (n_good 
%                                e active node)
%         best_five --> structure with info of best 5' time window 

%% SEARCH FOR BEST 5' TIME WINDOW 

% 3 criterion:

% 1) SNR over 90 percentile 
% 2) Lowest number of bad channels (Amplitude and spectral analysis)
% 3) Highest number of active nodes which are well covered from the array
%    sensitivity (Coverage threshold (*))


% (*) = 'Array Designer: automated optimized array design for functional
% near-infrared spectroscopy' - Brigadoi et. al.

%% FIRST CRITERION : SNR THRESHOLD

% Initialization of support variables
n_NIRS = length(t);
n_ch = size(d,2)/2;
Ts = round((t(2)-t(1)),4);
fs = 1/Ts;

% Samples inside 5' window (300 sec)
w_length = ((5*60)/Ts);
% Samples for 1' shift (60 sec)
s_shift = (60/Ts);
% Row = channels ; Coloumn = time windows
SNR_first_win = [];
idx_start_window = [];
idx_end_window = [];

% 5' moving window shifted of 1' : STD and SNR

for i = 1:s_shift:(n_NIRS-w_length)

    idx_start_window = [idx_start_window ; i];
    idx_end_window   = [idx_end_window ; i+w_length];
    d_tmp = d(i:i+w_length,:);

    SNR_first_ch = mean(mean(d_tmp)./std(d_tmp,[],1));
    
    SNR_first_win = [SNR_first_win SNR_first_ch];
end

% SNR THRESHOLD = 90 PERCENTILE OF SNR (Compute over all SNR for each
% window)
SNR_thr_first = quantile(SNR_first_win,0.90);
idx_over_snr_thr_first = find(SNR_first_win>SNR_thr_first);
SNR_over_snr_th_first = SNR_first_win(idx_over_snr_thr_first);

% Index of start and end of window above 90 percentile
idx_start_over_SNR = idx_start_window(idx_over_snr_thr_first);
idx_end_over_SNR = idx_end_window(idx_over_snr_thr_first);

% Number of good channel in time window above 90 percentile
n_pos_win = length(idx_start_over_SNR);
disp(['NUMBER OF TIME WINDOWS ABOVE 90 PERCENTILE OF SNR : ',num2str(n_pos_win)])

% Vector with number of window
n_win = length(SNR_first_win);
win_vect = [1:1:n_win];


%% SECOND CRITERION : BAD CHANNALES AND SPECTRAL ANALYSIS

% Bad channels removal and spectral analysis of possibile good window 
% (above 90 percentile SNR).

% Number of good channels
n_good = [];
% Vector 1 x nCh_tot: 1 good channel ; 0 Bad channel
remove_ch_tot = [];

warning('off')
for i = 1:1:n_pos_win

    d_tmp = d(idx_start_over_SNR(i):idx_end_over_SNR(i),:);
    
    [tmp_rem] = spectral_analysis(d_tmp,Ts,rmch{1},rmch{2},spect_param{1},spect_param{2},spect_param{3},spect_param{4},spect_param{5},spect_param{6},disp_single,disp_result);

    remove_ch_tot(:,i) = tmp_rem;
    n_good = [n_good ; sum(tmp_rem)];

end
warning('on')

% Vector of indexes of possible best window
idx_pos_win = [1:1:n_pos_win]';

% Median over all good channel 
th_good_ch = median(n_good);

% Keep only the possible best window above the median of the number of good
% channel
idx_win_after_good = find(n_good>=th_good_ch);                              % SE ERRORE CON UNA SOLA FINESTRA MAGGIORE E UGUALE 


% Table with index and number of good channels
tab_win_after_good        = table();
tab_win_after_good.idx    = idx_win_after_good;
tab_win_after_good.n_good = n_good(idx_win_after_good);

%% THIRD CRITERION : SPARSIFICATION

%% JACOBIAN ANALYSIS

% We decided arbitrarily to use the Jacobian for 850 nm wavelength. 

% Seperate Intensity and phase parts of the Jacobian
J_amp = J_850.complete; 
                                                                                                            
% EACH ROW IS A MEASUREMENT CORROSPONDING TO THE MESH.LINK
% EACH COLUMN IS A NODE CORROSPONDING TO THE MESH.NODE

%% MESH GENERATION

% Generating a mesh 'mesh_recon', 'mesh3' so that we can plot the Jacobian   
% on atals recon here

mesh_recon = mesh_850;
ind = reshape(mesh_recon.region(mesh_recon.elements),[],1);
ind = reshape(ind>=3,[],4);
ind = sum(ind,2);
ind = find(ind==4);
[mesh3.elements,mesh3.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
[mesh_all.elements,mesh_all.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
cortex_nodes = find(mesh_850.region==3);

%% ALL CHANNELS JACOBIAN PLOT

% Utility variable for jacobian plot
% This sums the jacobian for all measurements and nodes so we are left with 
% the total sensitivity for each node
all_J = sum(J_amp(:,:)); 
% Set mesh3.data as the total sensitivty
mesh3.data = all_J;

% Getting minimum and maximum J for just the cortex nodes. 
max_j = max(mesh3.data(1,cortex_nodes));
min_j = min(mesh3.data(1,cortex_nodes));

if f_disp == 1

    % It requires function plotniceimages_1
    plotniceimages_1(mesh3,mesh_recon); 
    caxis([min_j max_j])
    colorbar('horiz');
    title('ALL CHANNELS FOR \lambda = 850 nm')
end

%% COVERAGE THRESHOLD

% REF: 'Array Designer: automated optimized array design for functional
% near-infrared spectroscopy' - Brigadoi et. al.

% Compute the Coverage threshold as reported in the paper.
% From each possible window remove from the Jacobian the rows (channels)
% which corresponds to bad channels, sum along coloumns in order to obtain
% one value for each node to compare to the Coverage Threshold. If the node
% value is above the threhsold, that node is count as active.

% Coverage threshold computation
p_th    = 1; % [1 %]
act_vol = 1000;  % [1 cm^3 = 1000 mm^3]
V_all   = nodevolume(mesh_850.nodes,mesh_850.elements); % [Median Voronoi volume across all nodes mm^3]
V_gm    = median(V_all(mesh_850.region==3)); % [Median Voronoi volume across nodes of GM mm^3]

d_muA   = 0.001; % [0.001 mm^-1]

C_thresh = log((100+p_th)/100)/((act_vol / V_gm)*d_muA);                       

% Possible windows analysis
count_over_th = [];

for i = 1:1:n_pos_win

    % Bad channels removal
    tmp_J = J_amp;
    tmp_J(remove_ch_tot(1:n_ch,i)==0,:) = [];
    
    % Node value = sum along coloumn of J
    int_node = abs(sum(tmp_J(:,:)));
    
    % Coverage threshold comparison
    tmp_count = 0;
    for jN = 1:size(int_node,2);
        if int_node(jN) > C_thresh
            tmp_count = tmp_count + 1;
        end
    end

    % Number of active nodes for each window
    count_over_th = [count_over_th ; tmp_count];
end

% Analysis of the windows which have a number of good channels above the
% median
count_over_th_def = count_over_th(idx_win_after_good);
tab_win_after_good.spars = count_over_th_def;

%% BEST 5' TIME WINDOW

% Best 5' window is the one with the highest number of active nodes
% [max_count idx_max_mat] = max(tab_win_after_good.spars);

% Indexes and results
idx_best_win = tab_win_after_good.idx(idx_man);
n_best_SNR = SNR_over_snr_th_first(idx_best_win);
n_best_good_ch = tab_win_after_good.n_good(idx_man);
n_best_spars = tab_win_after_good.spars(idx_man);

%% JACOBIAN OF THE BEST 5' WINDOW

remCh_best = remove_ch_tot(:,idx_best_win);
% This sums the jacobian for all measurements and nodes so we are left with 
% the total sensitivity for EACH node
all_J_best = sum(J_amp(remCh_best(1:n_ch)==1,:)); 
% Set mesh3.data as the total sensitivty
mesh3.data = all_J_best; %

% Getting minimum and maximum J for just the cortex nodes. 
max_j = max(mesh3.data(1,cortex_nodes));
min_j = min(mesh3.data(1,cortex_nodes));

% Plotting the jacobian (requires function plotniceimages_1)

if f_disp == 1

    plotniceimages_1(mesh3,mesh_recon); 
    caxis([min_j max_j])
    colorbar('horiz');
    title('BEST 5'' WINDOW FOR \lambda = 850 nm')
end

%% BEST 5' WINDOW DISPLAY OF THE RESULTS

idx_start_best = idx_start_over_SNR(idx_best_win);
idx_end_best   = idx_end_over_SNR(idx_best_win);

d_best = d(idx_start_best:idx_end_best,:);
% t_best = t(idx_start_best:idx_end_best);                                 Works bad with concatenated NIRS data
t_best = [1:1:size(d_best,1)]';

% date_start_sincro  = (synchro_result.first_eu.time_start_first_eu_file + seconds(synchro_result.first_eu.idx_start_first_eu*Ts))
% date_end_sincro    = (synchro_result.first_eu.time_end_first_eu_file + seconds(synchro_result.first_eu.idx_end_first_eu*Ts));

date_start_sincro  = (synchro_result.time_start + seconds(synchro_result.idx_start*Ts));
date_end_sincro    = (synchro_result.time_end + seconds(synchro_result.idx_end*Ts));
date_start_first_win = date_start_sincro + seconds(idx_start_best*Ts);
date_end_first_win   = date_start_sincro + seconds(idx_end_best*Ts);

% date_start_eu  = PD.PD.Time(idx_eu_int.idx_start_first_eu_w(idx_PD));
% date_end_eu    = PD.PD.Time(idx_eu_int.idx_end_first_eu_w(idx_PD));
% date_start_first_win = date_start_eu + seconds(idx_start_best*Ts);
% date_end_first_win   = date_start_eu + seconds(idx_end_best*Ts);

% Display of the results
disp('===================================================================')
disp('BEST 5'' EUGLYCEMIA WINDOW')
disp(' ')
disp(['DATE START SYNCHRO : ',datestr(date_start_sincro)])
disp(['DATE END SYNCHRO   : ',datestr(date_end_sincro)])
disp(' ')
disp(['SNR THRESHOLD (90-th PERCENTILE) = ',num2str(SNR_thr_first)])
disp(['GOOD CHANNELS THRESHOLD (MEDIAN) = ',num2str(th_good_ch)])
disp(['COVERAGE THRESHOLD               = ',num2str(C_thresh)])
disp(' ')
disp(['SNR BEST 5'' WINDOW        = ',num2str(n_best_SNR)])
disp(['GOOD CHANNELS BEST 5'' WIN = ',num2str(n_best_good_ch/2),'/',num2str(size(d,2)/2)])
disp(['ACTIVE NODE BEST 5'' WIN   = ',num2str(n_best_spars),'/',num2str(size(int_node,2))])
disp(' ')
disp(['STARTING INDEX OF THE BEST 5'' WINDOW = ',num2str(idx_start_best)])
disp(['ENDING INDEX OF THE BEST 5'' WINDOW   = ',num2str(idx_end_best)])
disp(' ')
disp(['DATE START BEST 5'' WINDOW : ',datestr(date_start_first_win)])
disp(['DATE END BEST 5'' WINDOW   : ',datestr(date_end_first_win)])

% Output structure with information of best 5' time window
best_five = struct();
best_five.date_start_eu = date_start_sincro;
best_five.date_end_eu = date_end_sincro;
best_five.SNR_thresh = SNR_thr_first;
best_five.median_good_ch = th_good_ch;
best_five.C_thresh = C_thresh;
best_five.SNR = n_best_SNR;
best_five.good_ch = n_best_good_ch;
best_five.good_ch_single_w = n_best_good_ch/2;
best_five.active_node = n_best_spars;
best_five.idx_start = idx_start_best;
best_five.idx_end = idx_end_best;
best_five.date_start = date_start_first_win;
best_five.date_end = date_end_first_win;

%% DISPLAY OF RESULTS


if f_disp == 1

    %% Display SNR for each window and for the selected window
    figure()
    plot(win_vect,SNR_first_win)
    hold on
    plot(idx_over_snr_thr_first,SNR_first_win(idx_over_snr_thr_first),'*')
    plot(idx_over_snr_thr_first(tab_win_after_good.idx(1:end)),SNR_first_win(idx_over_snr_thr_first(tab_win_after_good.idx(1:end))),'*g')
    yline(SNR_thr_first,'--r')
    legend('SNR','Over th','90 Perc th')
    title('SNR for each Window')
    xlabel('N win')
    ylabel('SNR')

    %% NUMBER OF ACTIVE NODES FOR EACH POSSIBLE 5' TIME WINDOW
    figure()
    plot(count_over_th)
    hold on
    plot(idx_win_after_good,tab_win_after_good.spars,'*')
    legend('All win C.I.','Possible best win C.I.')
    title('Coverage Index over all 5'' windows')
    xlabel('Idx Window')
    ylabel('C.I.')

end

