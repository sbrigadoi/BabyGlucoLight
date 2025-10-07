% NIRS ANALYSIS

% PD 19 LAST EUGLYCEMIA INTERVAL

close all
clear all
clc

idx_PD = 19;

% RemoveNoisyChannel FUNCTION PARAMETERS
% SNRrange was set after visual inspection of channels which survived from
% the spectral analysis. Some 'weird' channels survived but they have very
% low SNR compared to the others. 
dRange = [5e-4 3];
SNRrange = 3;                                                              

% SPECTRAL ANALYSIS PARAMETERS :                                            

th_power_1 = 0.01;
th_power_2 = 0.01;

f_inf_1 = 0.05;
f_sup_1 = 2;

f_inf_2 = 3;
f_sup_2 = 4;

%% SETTING OF FOLDER PATH 

base_path = pwd;
data_path_cgm = fullfile(base_path,'DATI','CGM_INTERP');
data_path_nirs = fullfile(base_path,'CODICE','DEFINITIVO/','NIRS/','PD19');
function_path = fullfile(base_path,'CODICE','FUNCTION');
dati_ext_disc = fullfile('E:\DATI');
head_model_ext_disc = fullfile('E:\DATI\4DNeonatalHeadModel');

addpath(genpath(pwd))
addpath(genpath(data_path_cgm))
addpath(genpath(function_path))
addpath(genpath(data_path_nirs))
addpath(genpath(dati_ext_disc))
addpath(genpath(head_model_ext_disc))
addpath('E:\DATI\4DNeonatalHeadModel\MODEL')

%% LOAD SYNCHRONIZATION RESULTS
% Load of the synchronization results in order to display usefull
% information about the start and end time of the best 5' time window
load PD9_synchro_result.mat

last_synchro_result = struct();
last_synchro_result.time_start = synchro_result.last_eu.time_start_last_eu_file;
last_synchro_result.idx_start  = synchro_result.last_eu.idx_start_last_eu;
last_synchro_result.time_end   = synchro_result.last_eu.time_end_last_eu_file;
last_synchro_result.idx_end    = synchro_result.last_eu.idx_end_last_eu;

%% LAST EUGLYCEMIA INTERVAL 
% Load of the NIRS data for last Euglycemia interval
load PD19_NIRS_last_eu_tot.mat 

n_ch = size(d,2)/2;
Ts = round((t(2)-t(1)),4);
fs = 1/Ts;

%% LOAD OF JACOBIAN MATRIX

% Load of Jacobian matrix for each wavelength

% PD19 = 35+1 weeks + IUGR --> 29 gestation weeks

% Patient 19 is IUGR with a Birth Weigth and Cranial Circunference in the 0
% percentile for it's age. Therefore, we decided to use the smallest head
% model available (29 weeks).

load J_780_29w.mat
load mesh_780_29w.mat
mesh_780_29w = mesh;

load J_850_29w.mat
load mesh_850_29w.mat
mesh_850_29w = mesh;

disp('Jacobian and mesh completely loaded')

%% IDENTIFICATION OF THE BEST 5' TIME WINDOW 

% RemoveNoisyChannel parameters
rmch = cell(2,1);
rmch{1} = dRange;
rmch{2} = SNRrange;

% Spectral analysis parameters

spect_param = cell(6,1);
spect_param{1} = th_power_1;
spect_param{2} = th_power_2;
spect_param{3} = f_inf_1;
spect_param{4} = f_sup_1;
spect_param{5} = f_inf_2;
spect_param{6} = f_sup_2;

% Usefull variables
f_disp = 1;
disp_single = 0;
disp_result = 0;

% Best 5' time window 
[d_best,t_best,C_thresh,tab_win_after_good,best_five] = nirs_best_window(last_synchro_result,d,t,rmch,spect_param,J_850,mesh_850_29w,f_disp,disp_single,disp_result);

%% RAW NIRS SIGNAL DISPLAY FOR THE BEST 5' WINDOW 

figure()
subplot(2,1,1)
plot(t_best,d_best(:,1:n_ch))
title(['RAW NIRS DATA: LAST EUGLYCEMIA INTERVAL \lambda = ',num2str(SD.Lambda(1)),' nm'])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
axis tight

subplot(2,1,2)
plot(t_best,d_best(:,n_ch+1:end))
title(['RAW NIRS DATA: LAST EUGLYCEMIA INTERVAL \lambda = ',num2str(SD.Lambda(2)),' nm'])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
axis tight

%% BEST 5' WINDOW ANALYSIS AND PRE-PROCESSING

% th_power_1 = 0.12;
% th_power_2 = 0.01;
% f_inf_1 = 0.05;
% f_sup_1 = 2;
% f_inf_2 = 3;
% f_sup_2 = 4;

% Remove bad channel and spectral analysis selected window

disp_single = 0;
disp_result = 1;

[removeChSpect] = spectral_analysis(d_best,Ts,dRange,SNRrange,th_power_1,th_power_2,f_inf_1,f_sup_1,f_inf_2,f_sup_2,disp_single,disp_result);


% 1 GOOD ; 0 REMOVED
d_best_good = d_best(:,removeChSpect==1);
SD.MeasListAct = removeChSpect;
n_ch = size(d_best_good,2)/2;

% NIRS SIGNAL DISPLAY
figure()

subplot(2,1,1)
plot(t_best,d_best_good(:,1:n_ch))
title(['RAW NIRS DATA: BEST 5'' WIN, GOOD CH \lambda = ',num2str(SD.Lambda(1)),' nm'])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
axis tight

subplot(2,1,2)
plot(t_best,d_best_good(:,n_ch+1:end))
title(['RAW NIRS DATA: BEST 5'' WIN, GOOD CH \lambda = ',num2str(SD.Lambda(2)),' nm'])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
axis tight

% Channel by channel SNR visual inspection in order to adjust the 
% RemoveNoisyChannel threshold
SNRch = [];

for i=1:1:size(d_best_good,2)
    meanData = mean(d_best_good(:,i));
    stdData  = std(d_best_good(:,i),[],1);
    SNRch(i)  = meanData./stdData;
end

figure()
plot(SNRch,'-o')
title('SNR ch-by-ch')
xlabel('Channel')
ylabel('SNR')

% Display SNR minimum value
min_snr = min(SNRch);
max_snr = max(SNRch);
disp(['MIN SNR = ',num2str(min_snr)])
disp(['MAX SNR = ',num2str(max_snr)])

%% CHANNELS ARRAY CONFIGURATION

figure()
plot(SD.SrcPos(:,1), SD.SrcPos(:,2),'.r','MarkerSize',10)
hold on
plot(SD.DetPos(:,1), SD.DetPos(:,2),'.b','MarkerSize',10)
xlabel('X')
ylabel('Y')
zlabel('Z')

% Connection of source and detectors for each channel between good and bad
% channels

nCh = size(d,2)/2;
for iCh = 1:nCh
    % MeasList: first coloumn index of source second coloumn index of
    % detector
    % SrcPos / DetPos: 3D coordinates
    srcPos = SD.SrcPos(SD.MeasList(iCh,1),:);
    detPos = SD.DetPos(SD.MeasList(iCh,2),:);
    title('Channels array configuration')

    % 3D channels connection
    % NERI TENUTI     BLU RIMOSSI
    hold on
    if removeChSpect(iCh) == 1
        plot([srcPos(1) detPos(1)],[srcPos(2) detPos(2)],'k')
    else
        plot([srcPos(1) detPos(1)],[srcPos(2) detPos(2)],'c')
    end
end

% Euclidean distance between source-detector for each channels

for iCh=1:nCh
    % 3D channels coordinates
    srcPos = SD.SrcPos(SD.MeasList(iCh,1),:);
    detPos = SD.DetPos(SD.MeasList(iCh,2),:);
    % Euclidean distance between source-detector for each channels
    distCh(iCh,1) = sqrt(sum((srcPos-detPos).^2));
end

% Histogram (10 bins)
figure()
histogram(distCh,10)
title('Euclidean distance of source-detector')
xlabel(['Distance [',SD.SpatialUnit,']'])
ylabel('Number of channel')

%% SAVE OF THE BEST 5' WINDOW NIRS SIGNAL

% % Save of the NIRS signal: ALL CHANNELS, bad channels inof are inside SD

% save PD19_NIRS_best_last.mat t_best d_best d_best_good SD
% disp('Saving complete')

% load('PD19_NIRS_best_last.mat')

%% DoD CONVERTION

% All channel
meanValue = mean(d_best);
dodConv_all = -log(abs(d_best)./meanValue);

% Good channels only
nCh = size(d_best_good,2)/2;
meanValue_last_good = mean(d_best_good);
dodConv_last_good = -log(abs(d_best_good)./meanValue_last_good);

figure()
subplot(2,1,1)
plot(t_best,dodConv_last_good(:,1:nCh))
title(['Optical Density Change Best window only good channels, \lambda = ',num2str(SD.Lambda(1)),' nm'])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')
axis tight

subplot(2,1,2)
plot(t_best,dodConv_last_good(:,nCh+1:end))
title(['Optical Density Change Best window only good channels, \lambda = ', num2str(SD.Lambda(2)),' nm'])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')
axis tight

%% WAVELET MOTION ARTIFACTS CORRECTION

% Since Yang method may introduces some bias in MA correction it is decided 
% to use the Wavelet correction. Yang's methods uses MA identification that
% leads to different results with different MA identification parameters.
% Therefore, it introduces bias in the results, not only in terms of
% different strength of correlation but also in terms of significative
% correlation pattern. Wavelet MA correction works extremely well with
% spikes but it tends to not correct shift of baseline. However our data
% seems not to be damage with a lot of this artifact. It is necessary to 
% set an higher iqr in order to preserve low frequency informative data.
% Moreover, a visual inspection of the result channel by channel is 
% required. 

% WAVELET MOTION ARTIFACT CORRECTION

% Set an high iqr in order to save information at low frequencies
iqr = 1.2;
% Run wavelet motion correction
disp('START WAVELET')
dodWavelet = hmrMotionCorrectWavelet(dodConv_all,SD,iqr);
disp('END WAVELET')

%%
% Visual inspection of the wavelet results channel by channel
% t_best = [1:1:size(dodWavelet,1)];

% for i = 1:size(dodWavelet,2)
%     if removeChSpect(i)==1
%         figure()
%         plot(t_new,dodConv_all(:,i))
%         hold on
%         plot(t_new,dodWavelet(:,i))
%         axis tight
%         legend('RAW','WAVELET')
%         if i<(size(dodWavelet,2)/2)
%             title(['RAW vs WAVELET - ch: ',num2str(i),' first'])
%         else
%             title(['RAW vs WAVELET - ch: ',num2str(i),' second'])
%         end
%         xlabel('Time [s]')
%         ylabel('\DeltaOD [A.U.]')
%         set(gcf, 'Position', get(0, 'Screensize'));
%         pause
%         close                   
%     end
% end

%%
% Visualize dod signal before and after wavelet only for good channels

dodConv_good = dodConv_all(:,SD.MeasListAct==1);
dodWavelet_good = dodWavelet(:,SD.MeasListAct==1);

% save PD_19_LAST_Wavelet.mat dodWavelet dodWavelet_good t_best
% disp('Save wavelet complete')
% load PD_8_LAST_Wavelet.mat

figure()
subplot(2,1,1)
plot(t_best,dodConv_good(:,1:n_ch))
title('RAW NIRS FIRST WAVELENGTH')
axis tight
ylabel('\DeltaOD [A.U.]')
xlabel('Time [s]')

subplot(2,1,2)
plot(t_best,dodWavelet_good(:,1:nCh))
title('WAVELET NIRS FIRST WAVELENGTH')
axis tight
ylabel('\DeltaOD [A.U.]')
xlabel('Time [s]')


figure()
subplot(2,1,1)
plot(t_best,dodConv_good(:,n_ch+1:end))
title('RAW NIRS SECOND WAVELENGTH')
axis tight
ylabel('\DeltaOD [A.U.]')
xlabel('Time [s]')

subplot(2,1,2)
plot(t_best,dodWavelet_good(:,nCh+1:end))
title('WAVELET NIRS SECOND WAVELENGTH')
axis tight
ylabel('dod')
xlabel('Time [s]')


%% YANG ET. AL. MOTION ARTIFACTS CORRECTION 

% % QUESTO PASSAGGIO E' STATO FATTO DENTRO LA FUNCTION
% % 'hmrMotionArtifactbyChannel_byYang_Giacomo'
% 
% 
% % MOTION ART = Se la Mean SD dei dati in una finestra mobile supera una
% % certa soglia.
% 
% % Soglia = media sui valori SD più piccoli calcolati per ogni finestra
% % mobile. Valori piccoli corrispondono a porzioni noise free se calcolate
% % in un intervallo sufficientemente piccolo
% 
% % PIPELINE:
% % 1) OD originale
% % 2) Passa banda = 0-0.75 Hz (rimuovo le slow drift per calcolo STD)
% % 3) Per ogni finestra mobile di 4 secondi shiftata di 2 secondi calcolo
% %    STD (vedi formula paper)
% % 4) Ordino in maniera crescente i valorio STD per ogni finestra
% % 5) Calcolo la media sui 30% dei valori più piccoli (STD piccole si
% %    riferiscono a porzioni noise free)
% 
% %% MOTION ARTIFACT RECOGNITION
% 
% % PROVA (SEMBRA BUONO MA VEDERE CANALI REGOLARI SUL FONDO)
% tIncMan = ones(size(dodConv_all,1),1);
% tMotion = 0.5;
% tMask = 0.75;
% STDEVthresh = 9;
% AMPthresh = 7;
% fs = 1/Ts;
% 
% [tInc,tIncCh] = byYang_prova_gp(dodConv_all, fs, SD, tIncMan, tMotion, tMask, STDEVthresh, AMPthresh);
% 
% 
% % Visualize detected motion artifacts in good channels at both wavelength
% dodConvGood = dodConv_all(:,removeChSpect==1);
% t = [1:Ts:(5*60)+1]';
% 
% figure()
% plot(t,dodConvGood(:,:))
% title(['DOD Good Channel MA Identification BOTH Wavelength',num2str(STDEVthresh)])
% hold on;
% for i = 1:length(tInc)
%     if tInc(i) == 0 
%         lh = plot([t(i) t(i)],[min(min(dodConvGood)) max(max(dodConvGood))],'r','LineWidth',0.5); 
%         lh.Color = [lh.Color 0.07];
%     end
% end
% axis tight
% 
% 
% % SPLINE INTERPOLATION
% 
% p_hbo = 0.99; % Paper (almost a cubic spline, if p=0 LLS)
% dodSpline = hmrMotionCorrectSpline(dodConv_all,t_best,SD,tIncCh,p_hbo);
% 
% % Plot original optical density data in all good channels of first wavelength
% 
% figure()
% subplot(2,1,1)
% plot(t,dodConvGood(:,1:end/2))
% xlabel('Time [s]')
% ylabel('Optical density [A.U.]')
% xlim([t(1) t(end)])
% title('Wavelength 1: Original Good Ch')
% 
% % Plot spline corrected optical density data in all good channels of first wavelength
% dodSplineGood = dodSpline(:,removeChSpect == 1);
% subplot(2,1,2)
% plot(t,dodSplineGood(:,1:end/2))
% xlabel('Time [s]')
% ylabel('Optical density corrected [A.U.]')
% xlim([t(1) t(end)])
% title('Wavelength 1: Spline Correction Good Ch')
% 
% % Second wavelength result
% figure()
% subplot(2,1,1)
% plot(t,dodConvGood(:,(end/2)+1:end))
% xlabel('Time [s]')
% ylabel('Optical density [A.U.]')
% xlim([t(1) t(end)])
% title('Wavelength 2: Original Good Ch')
% 
% % Plot spline corrected optical density data in all good channels of second wavelength
% subplot(2,1,2)
% plot(t,dodSplineGood(:,(end/2)+1:end))
% xlabel('Time [s]')
% ylabel('Optical density corrected [A.U.]')
% xlim([t(1) t(end)])
% title('Wavelength 2: Spline Correction Good Ch')
% 
% %% RESULTS CHANNEL BY CHANNEL
% 
% % Compare uncorrected vs. spline corrected data at each good channel with
% % superimposed vertical lines for where artifacts were detected
% 
% % tIncCh of only good channels
% % tIncCh_good = tIncCh(:,removeChSpect==1);
% % 
% % n_tot = size(d_best,2); 
% % 
% % % Both Wavelength
% % for iCh = 1:n_tot
% %     if removeChSpect(iCh) == 1 % display only good channels
% % 
% %         figure;
% %         plot(t,dodConv_all(:,iCh))
% %         hold on;
% %         plot(t,dodSpline(:,iCh))
% %         if iCh<=n_tot/2
% %             title([num2str(iCh),' - First Wavelength'])
% %         else
% %             title([num2str(iCh),' - Second Wavelength'])
% %         end
% %         % legend('dod','dod Spline corrected')
% %         xlim([t(1) t(end)])
% %         hold on;
% %         for i = 1:size(tIncCh,1) 
% %             if tIncCh(i,iCh) == 0 % here we use tIncCh since we are displaying channels one at a time and we are interested in evaluating whether spline was able to correct the artifacts detected specifically in each single channel
% %                 lh = plot([t(i) t(i)],[-0.5 0.5],'r','LineWidth',0.5);
% %                 lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
% %             end
% %         end
% % 
% %         pause
% %         close
% %     end
% % end
% 
% %%  SECOND MOTION ARTIFACT RECOGNITION
% 
% % PROVA (SEMBRA ANDARE BENE MA VEDERE SE CONVIENE TOGLIERE ALTRI CANALI
% % RUMOROSI)
% tIncMan = ones(size(dodConv_all,1),1); 
% tMotion = 0.5;
% tMask = 0.75;
% STDEVthresh = 9;
% AMPthresh = 6;
% 
% 
% [tInc_GN,tIncCh_GN] = byYang_prova_gp(dodSpline, fs, SD, tIncMan, tMotion, tMask, STDEVthresh, AMPthresh);
% 
% 
% % Visualize detected motion artifacts in good channels at both wavelength
% t = [1:Ts:(5*60)+1]';
% 
% figure()
% plot(t,dodSplineGood(:,:))
% title('Second DOD Good Channel MA identification BOTH wavelength')
% hold on;
% for i = 1:length(tInc_GN) 
%     if tInc_GN(i) == 0
%         lh = plot([t(i) t(i)],[min(min(dodSplineGood)) max(max(dodSplineGood))],'r','LineWidth',0.5); 
%         lh.Color = [lh.Color 0.07]; 
%     end
% end
% axis tight
% 
% 
% % MOTION ARTIFACT CORRECTION PT.2: GAUSSIAN WHITE NOISE SUBSTITUCTION
% 
% dodGWN = hmrMotionCorrectGWN_Giacomo(dodSpline, t, SD, tIncCh_GN);
% 
% dodGWN_good = dodGWN(:,removeChSpect==1);
% 
% 
% figure()
% subplot(2,1,1)
% plot(t,dodSplineGood(:,1:end/2))
% title('DOD after spline, Good Channel, First Wavelength')
% axis tight
% subplot(2,1,2)
% plot(t,dodGWN_good(:,1:end/2))
% title('DOD after GWN substitution, Good Channel, First Wavelength')
% axis tight
% 
% figure()
% subplot(2,1,1)
% plot(t,dodSplineGood(:,(end/2)+1:end))
% title('DOD after spline, Good Channel, Second Wavelength')
% axis tight
% subplot(2,1,2)
% plot(t,dodGWN_good(:,(end/2)+1:end))
% title('DOD after GWN substitution, Good Channel, Second Wavelength')
% axis tight
% 
% 
% %% Plot signal after GWN substitution for both wavelength
% figure()
% subplot(2,1,1)
% plot(t,dodGWN_good(:,1:end/2))
% title('DOD after GWN substitution, Good Channel, First Wavelength')
% axis tight
% 
% subplot(2,1,2)
% plot(t,dodGWN_good(:,(end/2)+1:end))
% title('DOD after GWN substitution, Good Channel, Second Wavelength')
% axis tight
% 
% %% RESULTS CHANNEL BY CHANNEL
% 
% % Compare spline vs. GWN corrected data at each good channel with
% % superimposed vertical lines for where artifacts were detected
% % 
% % n_tot = size(d_best,2);
% % 
% % for iCh = 1:n_tot
% %     if removeChSpect(iCh) == 1 % display only good channels
% %         figure;
% %         plot(t,dodSpline(:,iCh))
% %         hold on;
% %         plot(t,dodGWN(:,iCh))    
% %         if iCh<=n_tot/2
% %             title([num2str(iCh),' - First Wavelength'])
% %         else
% %             title([num2str(iCh),' - Second Wavelength'])
% %         end
% %         % legend('dod Spline','dod GWN')
% %         xlim([t(1) t(end)])
% %         hold on;
% %         for i = 1:size(tIncCh_GN,1)  
% %             if tIncCh_GN(i,iCh) == 0 % here we use tIncCh since we are displaying channels one at a time and we are interested in evaluating whether spline was able to correct the artifacts detected specifically in each single channel
% %                 % lh = plot([t(i) t(i)],[min(dodSplineGood(:,iCh)) max(dodSplineGood(:,iCh))],'r','LineWidth',0.5);
% %                 lh = plot([t(i) t(i)],[-0.7 0.7],'r','LineWidth',0.5);
% %                 lh.Color = [lh.Color 0.07]; % the  first three entries define color, the fourth one transparency
% %                 axis tight
% % 
% %             end
% %         end
% % 
% %         pause
% %         close
% %     end
% % end
% 

% END YANG METHOD
%% BAND-PASS FILTER

% Apply band-pass filter
lowerCutOff = 0.009;
higherCutOff = 0.08;
dodFilt = hmrBandpassFilt(dodWavelet,fs,lowerCutOff,higherCutOff);

dodFilt_good = dodFilt(:,removeChSpect == 1);

figure()
plot(t_best,dodFilt_good(:,1:(end/2)))
title('DOD after Band Pass filter Hz - First Wavelength')
axis tight

figure()
plot(t_best,dodFilt_good(:,(end/2)+1:end))
title('DOD after Band Pass filter (0.01 - 0.5) Hz - Second Wavelength')
axis tight

%% SAVE OF THE DOD AFTER PRE-PROCESSING AND BAND PASS FILTER

% save PD19_NIRS_preprocessed_BP_last.mat dodFilt dodFilt_good t_best SD 
% disp('Saving complete')

% load('PD19_NIRS_preprocessed_BP_last.mat')
% disp('Load of pre-processed BP filtered NIRS signal complete')

%% IMAGES RECONSTRUCTION (INVERSE PROBLEM)

% Remove bad channels from the two Jacobian 
for i = 1:length(SD.Lambda) 

    if i == 1
        tmp = J_780.complete(:,:);
    end

    if i == 2
        tmp = J_850.complete(:,:);
    end
        JCropped{i} = tmp(SD.MeasListAct(SD.MeasList(:,4)==i)==1,:);
end

% Tikhonov regularization parameter
lambda1 = 0.1;                                                              

% Load age matched meshes (IUGR)
addpath('E:\DATI\4DNeonatalHeadModel\MODEL\29-weeks')
load('AllMeshes_29weeks.mshs','-mat')

% Load 'vol2gm' matrix, it maps values from volumetric to GM surface (IUGR)
addpath 'E:\DATI\4DNeonatalHeadModel\JACOBIAN_AND_MESH\29_WEEKS_JACOBIAN_AND_MESH'
load vol2gm_29w.mat

% Since vol2gm changes only with gestational week, it can be easely loaded

% vol2gm = Vol2GM_Transform(headVolumeMesh.node(:,1:3),gmSurfaceMesh.node,3,0);
% save vol2gm_29w.mat vol2gm
% disp('Vol2gm computed and saved')

% save_r = 1 --> Save results 
save_r = 0;

% Perform images reconstruction

[img] = img_recon(lambda1,JCropped,SD,dodFilt,headVolumeMesh,gmSurfaceMesh,vol2gm,save_r);
hbo = img.hbo;
hbr = img.hbr;

% load img_recon_PD19_last_eu.mat
% disp('Load of HbO and HbR complete')

%% RESTING STATE FUNCTIONAL CONNECTIVITY

% It was created an home-made function in order to perform the rsFC 

f_res = 1; % = 1 show results
r_ROI = 15;  % [mm]

% r_ROI = 8;   % [mm] NO OVERLAPPING

J_amp = J_850.complete;
all_J = sum(J_amp(:,:)); 

[ROI_table,ROI_active_table,tenTwenty,FC_struct] = rsFC(tenFive,gmSurfaceMesh,vol2gm,all_J,C_thresh,r_ROI,hbo,hbr,f_res);

% save PD_19_LAST_rsFC.mat FC_struct ROI_table
% disp('rsFC saved')

%% END OF rsFC ANALYSIS
