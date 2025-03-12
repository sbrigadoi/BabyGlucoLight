%% Simulate HRF and eval
clear all;clc;close all;

%1. Get some Euc data dod1
%2 Conv this to HbO HbR HbO1
%3 Sim HRF  HbO 1 sim
%4Recover OD data dod2
%5 apply corr. for MA dod3
%6 Recover HRF HbO 2
%7 see difference 
%% 1 Get Euc data. 
load('E:\PadovaPostDoc\BabyGluCo\NIRS_data\PD10\PD10 11.45.26 Wed 13 May 2020 20.nirs','-mat')
%load('E:\PadovaPostDoc\BabyGluCo\formatted_PD_GuyPerkins\PD10.mat')
PD_data.d = d;
PD_data.s = s;
PD_data.t = t;
PD_data.SD = SD;
PD_data.aux = aux;
clear d s t SD aux;

load('E:\PadovaPostDoc\BabyGluCo\NIRS_data\PD10\PD10 11.45.26 Wed 13 May 2020 21.nirs','-mat')
PD_data.d = [PD_data.d ;d];
PD_data.s = [PD_data.s ;s];
PD_data.t = [PD_data.t; t+max(PD_data.t)+0.1];
PD_data.aux = [PD_data.aux ;aux];
clear d s t SD aux;
fs = 10; %sampling rate (Hz)
%PD_data.glucose = getfield(PD10,'glucose',subject);
% 2 Signal Qaulity check and find MA trains
%Step 2. Signal Qaulity check
dRange = [5E-4 3];
SNRrange = 2; %2;
[goodch_idx, PD_data] = find_good_ch(PD_data,dRange,SNRrange);
PD_data.goodch_idx = goodch_idx;

% Step 3. Find M.A Train (1=no MA Train, 0 =MA train) (1= no MA, 0=MA)
%(if ch is in MA train idx, it doesn't have a MA train)
meanValue = mean(PD_data.d);
dodConv = -log(abs(PD_data.d)./meanValue);
dodConv_og = -log(abs(PD_data.d)./meanValue);
PD_data.dod = dodConv;
%
[MA_train_idx , tIncCh, IncCh_train, IncCh_train_loc] = find_MA_train(PD_data,dodConv,goodch_idx);
PD_data.tIncCh = tIncCh;
PD_data.IncCh_train_loc = IncCh_train_loc;
% Step 3.5 Show channels
ch_idx = goodch_idx;
data_type =1;
%show_data(PD_data,dodConv_og,ch_idx,data_type)
% DO interpolation MA corr.
% PD_data = MA_inter(PD_data); 
% MA_locations_true = [find(PD_data.tIncCh(:,i)==0)];
% MA_train_locations_true= find(PD_data.IncCh_train_loc(:,i)==0);
% 
% figure()
% plot(PD_data.inter_dod_new(:,1))
% 
% figure()
% plot(PD_data.t,PD_data.dod(:,1))
% xlabel('Time / s')
% ylabel('dOD / A.U')
% 
% figure()
% plot(PD_data.t,dodConv_og(:,1))
% hold on
% plot(PD_data.t(MA_locations_true),dodConv_og(MA_locations_true,1),'m.','MarkerSize',6)
% plot(PD_data.t(MA_train_locations_true),dodConv_og(MA_train_locations_true,1),'y.','MarkerSize',6)
% xlabel('Time / s')
% ylabel('dOD / A.U')
% 
% figure()
% plot(PD_data.t,PD_data.inter_dod_new(:,13))
% hold on
% plot(PD_data.t(MA_locations_true),PD_data.inter_dod_new(MA_locations_true,1),'m.','MarkerSize',6)
% plot(PD_data.t(MA_train_locations_true),PD_data.inter_dod_new(MA_train_locations_true,1),'y.','MarkerSize',6)
% xlabel('Time / s')
% ylabel('dOD / A.U')
% 
% figure()
% plot(PD_data.t,PD_data.dod(:,13))
% hold on
% plot(PD_data.t(MA_locations_true),PD_data.dod(MA_locations_true,13),'m.','MarkerSize',6)
% xlabel('Time / s')
% ylabel('dOD / A.U')
% 
% figure()
% plot(PD_data.t,PD_data.inter_dod_new(:,1))
% %sum(tIncCh(:,1)==0)
%sum(IncCh_train_loc(:,1)==0)
%std(dodConv_og(:,1))
%
% figure()
% subplot(2,1,1)
% plot(PD_data.t,dodConv_og(:,1))
% subplot(2,1,2)
% plot(PD_data.t,dodConv_og(:,1))
% hold on
% xline([find(tIncCh(:,1)==0)]/10);
% xline([find(IncCh_train_loc(:,1)==0)]/10 , 'r-');
% xlabel('Time / s');
% ylabel('DOD / A.U');
% 
% figure()
% plot(PD_data.t,dodConv_og(:,1))
% hold on
% xline([find(tIncCh(:,1)==0)]/10);
% xline([find(IncCh_train_loc(:,1)==0)]/10 , 'r-');
% xlabel('Time / s');
% ylabel('DOD / A.U');

%sum(tIncCh(:,1)==0)

% figure()
% plot(PD_data.t,PD_data.d(:,1))
% hold on
% xline([find(IncCh_train_loc(:,1)==0)]/fs);

%Bandpass Filter
%lowerCutOff = 0;
%higherCutOff = 0.01;
fs = 10; %sampling rate (Hz)
%dodFilt = hmrBandpassFilt(dodConv,fs,lowerCutOff,higherCutOff);

% Step 4.Y block avg/baseline
%dod = dodFilt; %this should be the final DOD form of the processed data.
dod=dodConv;
PD_data.dod1 = dod;
PD_data.dod = dod;
data_type=1;
show_data(PD_data,dod,ch_idx,data_type)

% Step 5. Downsample
% Reshape and down-sample data to 0.1 Hz.
%ds = fs*60; %1 min block avg (10frames per second * 60s = N frames in 1 min)
%PD_data = downsample_data(PD_data,ds);

% Step 6. Block Avg.
%Use first 5 minutes as baseline to look at DOD across all time %USE THIS
%baseline = [1 5]; %start and end of baseline
%PD_data = block_avg_data(PD_data,baseline);
%%%%%%data_type =2;
%data_type =1;
%show_data(PD_data,dodConv_og,ch_idx,data_type)

%figure()
%plot(PD_data.t,PD_data.dod(:,1))
%hold on
%xline([find(tIncCh(:,1)==0)]/fs);
%xline([find(IncCh_train_loc(:,1)==0)]/fs , 'r-');
%xlabel('Time / s');
%ylabel('DOD / A.U');

% Step 7. Recover HbO Hb HbT
%weekN = 30; %week age of infant
%[HbO_c,Hb_c,HbT_c] = recon_HbO_Hb(PD_data,weekN,goodch_idx);
% Step 7. Recover HbO Hb HbT using spectroscopy - ch by ch

%780 850 nm
%DPF 744 5.11+0.106 A^0.723
%DPF 807 4.99+0.067 A^0.814
%DPF 832 4.67+0.062 A^0.819
%A = years from full term

%Recently, Scholkmann and Wolf (2013) modeled the DPF as a function of age and wavelength, deriving the following equation (Equation 5) for calculation of DPF at any wavelength and age A:
%DPF(λ,A)=223.3+0.05624A0.8493−5.723∗10−7λ3
%+0.001245λ2−0.9025λ.    (5)

% ppf_780 = 4.99 + 0.067*( (1-(10/40))^(0.814));
% ppf_850 = 4.67 + 0.062*( (1-(10/40))^(0.819));
% PD_data.ppf = [ppf_780 ppf_850]; % DPF factor from Duncan 1998
% 
% age = weekN
% PD_data.SD.Lambda 

A=0; %age = 0Years
%DPF(λ,A)=
ppf_780(1,1) = (223.3)+ (0.05624*A^(0.8493));
ppf_780(1,2) = -(5.723*10^-7)*(PD_data.SD.Lambda(1)^(3));
ppf_780(1,3) =0.001245*(PD_data.SD.Lambda(1)^(2));
ppf_780(1,4) =-(0.9025*PD_data.SD.Lambda(1));

ppf_850(1,1) = (223.3)+ (0.05624*A^(0.8493));
ppf_850(1,2) = -(5.723*10^-7)*(PD_data.SD.Lambda(2)^(3));
ppf_850(1,3) =0.001245*(PD_data.SD.Lambda(2)^(2));
ppf_850(1,4) =-(0.9025*PD_data.SD.Lambda(2));

ppf_780 = sum(ppf_780);
ppf_850 = sum(ppf_850);
PD_data.ppf = [ppf_780 ppf_850];

dc_1 = hmrOD2Conc( PD_data.dod1, PD_data.SD, PD_data.ppf );

[HbO_x,Hb_x,HbT_x] = dc2HbX(dc_1);
PD_data.HbO_1 = HbO_x;
PD_data.Hb_1 = Hb_x;
PD_data.HbT_1 = HbT_x;
% Step 8. View Tomography
%cortex_nodes = view_HbO_Hb_tomography(HbO_c,Hb_c,HbT_c,goodch_idx,weekN);
%figure()
%plot(PD_data.t_ds,mean(HbO_c(cortex_nodes,:)))
% view changes
% data_type=3;
% show_data(PD_data,PD_data.HbO_1,ch_idx(1:end/2),data_type)
% 
% show_data(PD_data,PD_data.HbO_1,ch_idx(1),data_type)
% 
% figure()
% plot(PD_data.t,PD_data.HbO_1(:,1))
% hold on
% xline([find(tIncCh(:,1)==0)]/fs);
% xline([find(IncCh_train_loc(:,1)==0)]/fs , 'r-');
% xlabel('Time / s');
% ylabel('DOD / A.U');
% Step 9. Simulate HRF
%INPUT:
%               t = time
%               Fc = sampling frequency
%               data_RS = real data (one channel) in concentration changes
%               duration_hrf = length of the wanted HRF
%               nHRF = how many HRFs to add in this channel
%               distance = inter-trial interval (distance between two stimuli)
%               nChrom = 1 for HbO and 2 for HbR
%               block = 1 for event-related (standard HRF) and 2 for block design (conv with rect)
%               blockDuration = duration of the block
% OUTPUT:
%               yout = semisimulated signal (HRF+real data)
%               vett_hrf = vector with the added simulated HRFs
%               u = stimulus vector (the s of Homer)
%               t_hrf = timing of HRF
%               hrf_avg = average HRF across the nHRF simulated HRFs in this channel

%PD_data = simHRF(PD_data);
%
PD_data = simHRF_010224(PD_data);
figure()
plot(PD_data.t./60,PD_data.simulation(1).HRF.HbO(:,1))
hold on
plot(PD_data.t./60,PD_data.simulation(2).HRF.HbO(:,1))
plot(PD_data.t./60,PD_data.simulation(3).HRF.HbO(:,1))
plot(PD_data.t./60,PD_data.simulation(1).HRF.HbR(:,1))
plot(PD_data.t./60,PD_data.simulation(2).HRF.HbR(:,1))
plot(PD_data.t./60,PD_data.simulation(3).HRF.HbR(:,1))
%plot(PD_data.t./60,PD_data.HbO_1(:,1),'r.')
%plot(PD_data.t./60,PD_data.Hb_1(:,1),'b.')
xline([900/60],'b--')
xline([5100/60],'b--')
xlim([0 100])
xlabel('Time / Min');
ylabel('Sim HRF / M');
legend('Low SNR - HbO','Medium SNR - HbO','High SNR - HbO','Low SNR - Hb','Medium SNR - Hb','High SNR - Hb');
%legend('Low SNR - HbO','Medium SNR - HbO','High SNR - HbO','Low SNR - Hb','Medium SNR - Hb','High SNR - Hb','HbO 1','Hb 1');
% legend('Low SNR - Duration Short','Medium SNR - Duration Short','High SNR - Duration Short','Low SNR - Duration Long','Medium SNR - Duration Long','High SNR - Duration Long')

%HRF Sim ONLY, without the added data
figure()
plot(PD_data.t./60,PD_data.simulationONLY(1).HRF.HbO(:,1))
hold on
plot(PD_data.t./60,PD_data.simulationONLY(2).HRF.HbO(:,1))
plot(PD_data.t./60,PD_data.simulationONLY(3).HRF.HbO(:,1))
plot(PD_data.t./60,PD_data.simulationONLY(1).HRF.HbR(:,1))
plot(PD_data.t./60,PD_data.simulationONLY(2).HRF.HbR(:,1))
plot(PD_data.t./60,PD_data.simulationONLY(3).HRF.HbR(:,1))
xline([900/60],'b--')
xline([5100/60],'b--')
xlim([0 100])
xlabel('Time / Min');
ylabel('Sim HRF / M');
legend('Low SNR - HbO','Medium SNR - HbO','High SNR - HbO','Low SNR - Hb','Medium SNR - Hb','High SNR - Hb');
%legend('Low SNR - HbO','Medium SNR - HbO','High SNR - HbO','Low SNR - Hb','Medium SNR - Hb','High SNR - Hb','HbO 1','Hb 1');
% legend('Low SNR - Duration Short','Medium SNR - Duration Short','High SNR - Duration Short','Low SNR - Duration Long','Medium SNR - Duration Long','High SNR - Duration Long')
%
PD_data.HbO_sim1_lowSNR= PD_data.simulation(1).HRF.HbO;
PD_data.HbO_sim1_medSNR= PD_data.simulation(2).HRF.HbO;
PD_data.HbO_sim1_highSNR= PD_data.simulation(3).HRF.HbO;

PD_data.Hb_sim1_lowSNR= PD_data.simulation(1).HRF.HbR;
PD_data.Hb_sim1_medSNR= PD_data.simulation(2).HRF.HbR;
PD_data.Hb_sim1_highSNR= PD_data.simulation(3).HRF.HbR;

PD_data.HbT_sim1_lowSNR= PD_data.simulation(1).HRF.HbO + PD_data.simulation(1).HRF.HbR;
PD_data.HbT_sim1_medSNR= PD_data.simulation(2).HRF.HbO + PD_data.simulation(2).HRF.HbR;
PD_data.HbT_sim1_highSNR= PD_data.simulation(3).HRF.HbO + PD_data.simulation(3).HRF.HbR;

%dc_sim1 = [PD_data.HbO_sim1_lowSNR, PD_data.Hb_sim1_lowSNR, (PD_data.HbO_sim1_lowSNR+PD_data.Hb_sim1_lowSNR)];
%dc_sim1_lowSNR = cat(3,PD_data.HbO_sim1_lowSNR, PD_data.Hb_sim1_lowSNR, (PD_data.HbO_sim1_lowSNR+PD_data.Hb_sim1_lowSNR));
% PD_data.dc_1_sim_lowSNR = cat(3,PD_data.HbO_sim1_lowSNR, PD_data.Hb_sim1_lowSNR, (PD_data.HbO_sim1_lowSNR+PD_data.Hb_sim1_lowSNR));
% PD_data.dc_1_sim_medSNR = cat(3,PD_data.HbO_sim1_medSNR, PD_data.Hb_sim1_medSNR, (PD_data.HbO_sim1_medSNR+PD_data.Hb_sim1_medSNR));
% PD_data.dc_1_sim_highSNR = cat(3,PD_data.HbO_sim1_highSNR, PD_data.Hb_sim1_highSNR, (PD_data.HbO_sim1_highSNR+PD_data.Hb_sim1_highSNR));

PD_data.dc_1_sim_lowSNR(:,1,:) = PD_data.HbO_sim1_lowSNR;
PD_data.dc_1_sim_lowSNR(:,2,:) = PD_data.Hb_sim1_lowSNR;
PD_data.dc_1_sim_lowSNR(:,3,:) = PD_data.HbT_sim1_lowSNR;

PD_data.dc_1_sim_medSNR(:,1,:) = PD_data.HbO_sim1_medSNR;
PD_data.dc_1_sim_medSNR(:,2,:) = PD_data.Hb_sim1_medSNR;
PD_data.dc_1_sim_medSNR(:,3,:) = PD_data.HbT_sim1_medSNR;

PD_data.dc_1_sim_highSNR(:,1,:) = PD_data.HbO_sim1_highSNR;
PD_data.dc_1_sim_highSNR(:,2,:) = PD_data.Hb_sim1_highSNR;
PD_data.dc_1_sim_highSNR(:,3,:) = PD_data.HbT_sim1_highSNR;
% Turn concentration changes back to dOD
% INPUTS:
% dc: the concentration data (#time points x 3 x #SD pairs
%     3 concentrations are returned (HbO, HbR, HbT)
% SD:  the SD structure
% ppf: partial pathlength factors for each wavelength. If there are 2
%      wavelengths of data, then this is a vector ot 2 elements.
%      Typical value is ~6 for each wavelength if the absorption change is 
%      uniform over the volume of tissue measured. To approximate the
%      partial volume effect of a small localized absorption change within
%      an adult human head, this value could be as small as 0.1.
% OUTPUTS:
% dod: the change in OD (#time points x #channels)
% 5) Convert the semi-simulated data back to attenuation changes
PD_data.dod2_medSNR = hmrConc2OD( PD_data.dc_1_sim_medSNR, PD_data.SD, PD_data.ppf );
PD_data.dod2_lowSNR = hmrConc2OD( PD_data.dc_1_sim_lowSNR, PD_data.SD, PD_data.ppf );
PD_data.dod2_highSNR = hmrConc2OD( PD_data.dc_1_sim_highSNR, PD_data.SD, PD_data.ppf );

figure()
plot(PD_data.t, PD_data.dod1(:,3))
hold on
plot(PD_data.t, PD_data.dod2_lowSNR(:,3))
plot(PD_data.t, PD_data.dod2_medSNR(:,3))
plot(PD_data.t, PD_data.dod2_highSNR(:,3))
xlabel('Time / s');
ylabel('DOD / A.U')
legend('DOD 1','DOD 2 lowSNR','DOD 2 medSNR','DOD 2 highSNR')

clear A dc_1 dc_2_lowSNR dc_2_medSNR dc_2_highSNR dod dodConv dodConv_og Hb_x HbO_x HbT_x ;
%% test SplineSG
p=0.99;
FrameSize_sec = 0.5; %1/fs;
[dod2_SG_medSNR ,tIncCh_baseline_dod2_medSNR,tInc_baseline_dod2_medSNR] = hmrMotionCorrectSplineSG(PD_data.dod2_medSNR, PD_data.dod2_medSNR, PD_data.t, PD_data.SD, p, FrameSize_sec,1);

figure()
subplot(4,1,1)
plot(PD_data.t,PD_data.dod2_medSNR(:,1))
hold on
plot(PD_data.t,dod2_SG_medSNR(:,1))
plot(PD_data.t(find(tIncCh_baseline_dod2_medSNR(:,1)==0)),PD_data.dod2_medSNR(find(tIncCh_baseline_dod2_medSNR(:,1)==0),1),'m.','MarkerSize',1)
legend('dod2 medSNR','dod2 SG medSNR','BL shifts')
subplot(4,1,2)
plot(PD_data.t,PD_data.dod2_medSNR(:,10))
hold on
plot(PD_data.t,dod2_SG_medSNR(:,10))
plot(PD_data.t(find(tIncCh_baseline_dod2_medSNR(:,10)==0)),PD_data.dod2_medSNR(find(tIncCh_baseline_dod2_medSNR(:,10)==0),10),'m.','MarkerSize',1)
legend('dod2 medSNR','dod2 SG medSNR','BL shifts')
subplot(4,1,3)
plot(PD_data.t,PD_data.dod2_medSNR(:,13))
hold on
plot(PD_data.t,dod2_SG_medSNR(:,13))
plot(PD_data.t(find(tIncCh_baseline_dod2_medSNR(:,13)==0)),PD_data.dod2_medSNR(find(tIncCh_baseline_dod2_medSNR(:,13)==0),13),'m.','MarkerSize',1)
legend('dod2 medSNR','dod2 SG medSNR','BL shifts')
subplot(4,1,4)
plot(PD_data.t,PD_data.dod2_medSNR(:,20))
hold on
plot(PD_data.t,dod2_SG_medSNR(:,20))
plot(PD_data.t(find(tIncCh_baseline_dod2_medSNR(:,20)==0)),PD_data.dod2_medSNR(find(tIncCh_baseline_dod2_medSNR(:,20)==0),20),'m.','MarkerSize',1)
legend('dod2 medSNR','dod2 SG medSNR','BL shifts')

%% test PCA
nSV=13;
[dN,svs,nSV] = hmrMotionCorrectPCA( PD_data.SD, PD_data.dod2_medSNR, tInc_baseline_dod2_medSNR, nSV );

figure()
plot(PD_data.t,PD_data.dod2_medSNR(:,13))
hold on
plot(PD_data.t,dN(:,13))

%% Apply Motion Correction Techniques
for k=[5]
    %motion correction here
    % 5) Convert the semi-simulated data back to attenuation changes
    PD_data.dod2_medSNR = hmrConc2OD( PD_data.dc_1_sim_medSNR, PD_data.SD, PD_data.ppf );
    PD_data.dod2_lowSNR = hmrConc2OD( PD_data.dc_1_sim_lowSNR, PD_data.SD, PD_data.ppf );
    PD_data.dod2_highSNR = hmrConc2OD( PD_data.dc_1_sim_highSNR, PD_data.SD, PD_data.ppf );

    %IF NO Motion Correction Technique - run this
    PD_data.dod3noc_medSNR = PD_data.dod2_medSNR;
    PD_data.dod3noc_lowSNR = PD_data.dod2_lowSNR;
    PD_data.dod3noc_highSNR = PD_data.dod2_highSNR;
    
    PD_data = PD_data_motioncorr(PD_data,k); %either apply MA correction
    ch_idx = PD_data.goodch_idx;
    %0 = no corr.
    %1 = only BP
    %2 = only wavelet + bp

    %3 = only spline + bp(HmR MA dect.)
    %4 = Spline + Wavelet + bp(HmR MA dect.)

    %5 = only spline + bp (GVDT MA dect.)
    %6 = spline + wavelet + bp(GVDT MA dect.)
    %IF NO Motion Correction Technique - run this
    %PD_data.dod3noc_medSNR = PD_data.dod2_medSNR;
    %PD_data.dod3noc_lowSNR = PD_data.dod2_lowSNR;
    %PD_data.dod3noc_highSNR = PD_data.dod2_highSNR;

    %PD_data.dod3 = PD_data.inter_dod_new;
    %PD_data.dod3 = PD_data.dod2; %or don't
    % baseline correction do this after motion correction.
    %Use first 5 minutes as baseline to look at DOD across all time %USE THIS
    %baseline = [1 5]; %start and end of baseline
    %PD_data = block_avg_data(PD_data,baseline);
    % IF NO Motion Correction Technique - run this
    %PD_data.dod3noc_medSNR = PD_data.dod2_medSNR; %or don't
    %PD_data.dod3noc_lowSNR = PD_data.dod2_lowSNR; %or don't
    %PD_data.dod3noc_highSNR = PD_data.dod2_highSNR; %or don't
    % Recover HbO Hb
    %weekN = 30; %week age of infant
    %[HbO_c,Hb_c,HbT_c] = recon_HbO_Hb(PD_data,weekN,goodch_idx);
    dc_2_medSNR = hmrOD2Conc( PD_data.dod3_medSNR, PD_data.SD, PD_data.ppf );
    dc_2_lowSNR = hmrOD2Conc( PD_data.dod3_lowSNR, PD_data.SD, PD_data.ppf );
    dc_2_highSNR = hmrOD2Conc( PD_data.dod3_highSNR, PD_data.SD, PD_data.ppf );

    dc_2noc_medSNR = hmrOD2Conc( PD_data.dod3noc_medSNR, PD_data.SD, PD_data.ppf );
    dc_2noc_lowSNR = hmrOD2Conc( PD_data.dod3noc_lowSNR, PD_data.SD, PD_data.ppf );
    dc_2noc_highSNR = hmrOD2Conc( PD_data.dod3noc_highSNR, PD_data.SD, PD_data.ppf );

    [HbO_x,Hb_x,HbT_x] = dc2HbX(dc_2_medSNR);
    PD_data.HbO_2_medSNR = HbO_x;
    PD_data.Hb_2_medSNR = Hb_x;
    PD_data.HbT_2_medSNR = HbT_x;

    [HbO_x,Hb_x,HbT_x] = dc2HbX(dc_2_lowSNR);
    PD_data.HbO_2_lowSNR = HbO_x;
    PD_data.Hb_2_lowSNR = Hb_x;
    PD_data.HbT_2_lowSNR = HbT_x;

    [HbO_x,Hb_x,HbT_x] = dc2HbX(dc_2_highSNR);
    PD_data.HbO_2_highSNR = HbO_x;
    PD_data.Hb_2_highSNR = Hb_x;
    PD_data.HbT_2_highSNR = HbT_x;

    [HbO_x,Hb_x,HbT_x] = dc2HbX(dc_2noc_medSNR);
    PD_data.HbO_2noc_medSNR = HbO_x;
    PD_data.Hb_2noc_medSNR = Hb_x;
    PD_data.HbT_2noc_medSNR = HbT_x;

    [HbO_x,Hb_x,HbT_x] = dc2HbX(dc_2noc_lowSNR);
    PD_data.HbO_2noc_lowSNR = HbO_x;
    PD_data.Hb_2noc_lowSNR = Hb_x;
    PD_data.HbT_2noc_lowSNR = HbT_x;

    [HbO_x,Hb_x,HbT_x] = dc2HbX(dc_2noc_highSNR);
    PD_data.HbO_2noc_highSNR = HbO_x;
    PD_data.Hb_2noc_highSNR = Hb_x;
    PD_data.HbT_2noc_highSNR = HbT_x;

    clear A dc_1 dc_2_lowSNR dc_2_medSNR dc_2_highSNR dc_2noc_lowSNR dc_2noc_medSNR dc_2noc_highSNR dod dodConv dodConv_og Hb_x HbO_x HbT_x ;

    %
    % figure()
    % plot(PD_data.t, PD_data.HbO_1(:,1))
    % hold on
    % plot(PD_data.t, PD_data.HbO_sim1_medSNR(:,1))
    % plot(PD_data.t, PD_data.HbO_2_medSNR(:,1))
    % xlabel('Time / s');
    % ylabel('Sim HRF + HRF / M')
    % legend('HbO','HbO + Sim HRF','HbO + Sim HRF Recovered')
    %
    % figure()
    % plot(PD_data.t, PD_data.HbO_1(:,1))
    % hold on
    % plot(PD_data.t, PD_data.HbO_sim1_medSNR(:,1))
    % plot(PD_data.t, PD_data.HbO_2_medSNR(:,1))
    % xlabel('Time / s');
    % ylabel('Sim HRF + HRF / M')
    % legend('HbO','HbO + Sim HRF','HbO + Sim HRF MA Correction')

    % baseline correction for HbX2 (with no correction and correction))
    baseline = [1 5*10*60]; %start and end of baseline
    %PD_data = block_avg_data(PD_data,baseline);
    PD_data.HbO_2noc_medSNR = PD_data.HbO_2noc_medSNR(:,:) - mean(PD_data.HbO_2noc_medSNR(baseline(1):baseline(2),:));
    PD_data.HbO_2_medSNR = PD_data.HbO_2_medSNR(:,:) - mean(PD_data.HbO_2_medSNR(baseline(1):baseline(2),:));
    PD_data.HbO_2noc_lowSNR = PD_data.HbO_2noc_lowSNR(:,:) - mean(PD_data.HbO_2noc_lowSNR(baseline(1):baseline(2),:));
    PD_data.HbO_2_lowSNR = PD_data.HbO_2_lowSNR(:,:) - mean(PD_data.HbO_2_lowSNR(baseline(1):baseline(2),:));
    PD_data.HbO_2noc_highSNR = PD_data.HbO_2noc_highSNR(:,:) - mean(PD_data.HbO_2noc_highSNR(baseline(1):baseline(2),:));
    PD_data.HbO_2_highSNR = PD_data.HbO_2_highSNR(:,:) - mean(PD_data.HbO_2_highSNR(baseline(1):baseline(2),:));

    PD_data.Hb_2noc_medSNR = PD_data.Hb_2noc_medSNR(:,:) - mean(PD_data.Hb_2noc_medSNR(baseline(1):baseline(2),:));
    PD_data.Hb_2_medSNR = PD_data.Hb_2_medSNR(:,:) - mean(PD_data.Hb_2_medSNR(baseline(1):baseline(2),:));
    PD_data.Hb_2noc_lowSNR = PD_data.Hb_2noc_lowSNR(:,:) - mean(PD_data.Hb_2noc_lowSNR(baseline(1):baseline(2),:));
    PD_data.Hb_2_lowSNR = PD_data.Hb_2_lowSNR(:,:) - mean(PD_data.Hb_2_lowSNR(baseline(1):baseline(2),:));
    PD_data.Hb_2noc_highSNR = PD_data.Hb_2noc_highSNR(:,:) - mean(PD_data.Hb_2noc_highSNR(baseline(1):baseline(2),:));
    PD_data.Hb_2_highSNR = PD_data.Hb_2_highSNR(:,:) - mean(PD_data.Hb_2_highSNR(baseline(1):baseline(2),:));

    % figure()
    % plot(PD_data.t, PD_data.HbO_2noc_lowSNR(:,44))
    % hold on
    % plot(PD_data.t, PD_data.HbO_2_lowSNR(:,44),'LineWidth',2)
    % plot(PD_data.t, PD_data.simulationONLY(1).HRF.HbO(:,1),'LineWidth',2)
    % xlabel('Time / s')
    % ylabel('delta HbO / M')
    % xline([900],'b--')
    % xline([5100],'b--')
    % legend('No correction','Bandpass','HRF Sim')
    % title('Low SNR')

    % figure()
    % plot(PD_data.t, PD_data.HbO_2noc_medSNR(:,44))
    % hold on
    % plot(PD_data.t, PD_data.HbO_2_medSNR(:,44),'LineWidth',2)
    % plot(PD_data.t, PD_data.simulationONLY(2).HRF.HbO(:,1),'LineWidth',2)
    % xlabel('Time / s')
    % ylabel('delta HbO / M')
    % xline([900],'b--')
    % xline([5100],'b--')
    % legend('No correction','Bandpass','HRF Sim')
    % title('Med SNR')

    figure()
    plot(PD_data.t, PD_data.HbO_2noc_highSNR(:,44))
    hold on
    plot(PD_data.t, PD_data.HbO_2_highSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(3).HRF.HbO(:,1),'LineWidth',2)
    xlabel('Time / s')
    ylabel('delta HbO / M')
    xline([900],'b--')
    xline([5100],'b--')
    legend('No correction','Bandpass','HRF Sim')
    title('High SNR')

    figure()
    plot(PD_data.t, PD_data.HbO_2_lowSNR(:,44),'LineWidth',2)
    hold on
    plot(PD_data.t, PD_data.HbO_2_medSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.HbO_2_highSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(1).HRF.HbO(:,1),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(2).HRF.HbO(:,1),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(3).HRF.HbO(:,1),'LineWidth',2)
    legend('Bandpass lowSNR','Bandpass medSNR','Bandpass highSNR','HRFsim lowSNR','HRFsim medSNR','HRFsim highSNR')
    xlabel('Time / s')
    ylabel('delta HbO / M')
    xline([900],'b--')
    xline([5100],'b--')

    figure()
    plot(PD_data.t, PD_data.HbO_2noc_lowSNR(:,44),'LineWidth',2)
    hold on
    plot(PD_data.t, PD_data.HbO_2noc_medSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.HbO_2noc_highSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.HbO_2_lowSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.HbO_2_medSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.HbO_2_highSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(1).HRF.HbO(:,1),'--','LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(2).HRF.HbO(:,1),'--','LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(3).HRF.HbO(:,1),'--','LineWidth',2)
    legend('no corr lowSNR','no corr medSNR','no corr highSNR','bp lowSNR','bp medSNR','bp highSNR','sim lowSNR','sim medSNR','sim highSNR')
    xlabel('Time / s')
    ylabel('delta HbO / M')
    xline([900],'b--')
    xline([5100],'b--')

    figure()
    subplot(1,3,1)
    plot(PD_data.t, PD_data.HbO_2noc_lowSNR(:,44))
    hold on
    plot(PD_data.t, PD_data.HbO_2_lowSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(1).HRF.HbO(:,1),'LineWidth',2)
    xlabel('Time / s')
    ylabel('delta HbO / M')
    xline([900],'b--')
    xline([5100],'b--')
    subplot(1,3,2)
    plot(PD_data.t, PD_data.HbO_2noc_medSNR(:,44))
    hold on
    plot(PD_data.t, PD_data.HbO_2_medSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(2).HRF.HbO(:,1),'LineWidth',2)
    xlabel('Time / s')
    ylabel('delta HbO / M')
    xline([900],'b--')
    xline([5100],'b--')
    subplot(1,3,3)
    plot(PD_data.t, PD_data.HbO_2noc_highSNR(:,44))
    hold on
    plot(PD_data.t, PD_data.HbO_2_highSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(3).HRF.HbO(:,1),'LineWidth',2)
    xlabel('Time / s')
    ylabel('delta HbO / M')
    xline([900],'b--')
    xline([5100],'b--')

    figure()
    plot(PD_data.t, PD_data.HbO_2noc_medSNR(:,44))
    hold on
    plot(PD_data.t, PD_data.HbO_2_medSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(2).HRF.HbO(:,1),'LineWidth',2)
    plot(PD_data.t, PD_data.HbO_2noc_highSNR(:,44))
    plot(PD_data.t, PD_data.HbO_2_highSNR(:,44),'LineWidth',2)
    plot(PD_data.t, PD_data.simulationONLY(3).HRF.HbO(:,1),'LineWidth',2)
    xlabel('Time / s')
    ylabel('delta HbO / M')
    xline([900],'b--')
    xline([5100],'b--')
    legend('No corr med','Bandpass med','HRF Sim med','No corr high','Bandpass high','HRF Sim high')
    title('Med and High SNR')

    figure()
    sgtitle('Low SNR')
    for i=1:12
        subplot(3,4,i)
        plot(PD_data.t, PD_data.HbO_2noc_lowSNR(:,ch_idx(i)))
        hold on
        plot(PD_data.t, PD_data.HbO_2_lowSNR(:,ch_idx(i)),'LineWidth',2)
        plot(PD_data.t, PD_data.simulationONLY(1).HRF.HbO(:,1),'LineWidth',2)
        xline([900],'b--')
        xline([5100],'b--')
        title("ch "+num2str(ch_idx(i)))
        xlabel('T / s')
        ylabel('dHbO / M')
    end

    figure()
    sgtitle('Med SNR')
    for i=1:12
        subplot(3,4,i)
        plot(PD_data.t, PD_data.HbO_2noc_medSNR(:,ch_idx(i)))
        hold on
        plot(PD_data.t, PD_data.HbO_2_medSNR(:,ch_idx(i)),'LineWidth',2)
        plot(PD_data.t, PD_data.simulationONLY(2).HRF.HbO(:,1),'LineWidth',2)
        xline([900],'b--')
        xline([5100],'b--')
        title("ch "+num2str(ch_idx(i)))
        xlabel('T / s')
        ylabel('dHbO / M')
    end

    figure()
    sgtitle('High SNR')
    for i=1:12
        subplot(3,4,i)
        plot(PD_data.t, PD_data.HbO_2noc_highSNR(:,ch_idx(i)))
        hold on
        plot(PD_data.t, PD_data.HbO_2_highSNR(:,ch_idx(i)),'LineWidth',2)
        plot(PD_data.t, PD_data.simulationONLY(3).HRF.HbO(:,1),'LineWidth',2)
        xline([900],'b--')
        xline([5100],'b--')
        title("ch "+num2str(ch_idx(i)))
        xlabel('T / s')
        ylabel('dHbO / M')
    end
    % Evaluate changes
    % Getting Error metric between 2nd recovered HbX and the Simulated HbX without data.
    %PD_data.deltaHRF_medSNR = PD_data.HbO_2_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1);
    %PD_data.deltaHRF_rs_medSNR = (sqrt( PD_data.HbO_2_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1)).^2);
    PD_data.deltaHRF_rms.HbO_medSNR = sqrt(  mean( (PD_data.HbO_2_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1)).^2));
    PD_data.deltaHRF_rms.HbO_lowSNR = sqrt(  mean( (PD_data.HbO_2_lowSNR - PD_data.simulationONLY(1).HRF.HbO(:,1)).^2));
    PD_data.deltaHRF_rms.HbO_highSNR = sqrt(  mean( (PD_data.HbO_2_highSNR - PD_data.simulationONLY(3).HRF.HbO(:,1)).^2));

    %new RMS (right)
    %PD_data.deltaHRF_rms.HbO_medSNR = sqrt(  mean(  (PD_data.HbO_2_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1)).^2 ));
    %PD_data.deltaHRF_rms.HbO_noc_medSNR = sqrt(  mean( (PD_data.HbO_2noc_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1)).^2));

    %old RMS (wrong)
    %PD_data.deltaHRF_rms.HbO_noc_medSNR = (sqrt(  mean (PD_data.HbO_2noc_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1)).^2));

    PD_data.deltaHRF_rms.HbO_noc_medSNR = sqrt(  mean( (PD_data.HbO_2noc_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1)).^2));
    PD_data.deltaHRF_rms.HbO_noc_lowSNR = sqrt(  mean( (PD_data.HbO_2noc_lowSNR - PD_data.simulationONLY(1).HRF.HbO(:,1)).^2));
    PD_data.deltaHRF_rms.HbO_noc_highSNR = sqrt(  mean( (PD_data.HbO_2noc_highSNR - PD_data.simulationONLY(3).HRF.HbO(:,1)).^2));

    PD_data.deltaHRF_rs.HbO_medSNR = (sqrt(  (PD_data.HbO_2_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1)).^2));
    PD_data.deltaHRF_rs.HbO_noc_medSNR = (sqrt(  (PD_data.HbO_2noc_medSNR - PD_data.simulationONLY(2).HRF.HbO(:,1)).^2));

    figure()
    plot(PD_data.deltaHRF_rs.HbO_noc_medSNR(:,1))
    hold on
    plot(PD_data.deltaHRF_rs.HbO_medSNR(:,1),'LineWidth',2)
    xlabel('Time / S')
    ylabel('RS error')
    legend('RS no corr.','RS bp')
    % figure()
    % plot(PD_data.simulationONLY(3).HRF.HbO(:,1))
    % hold on
    % plot(PD_data.HbO_sim1_highSNR(:,3))
    % plot(PD_data.HbO_2_highSNR(:,4))

    PD_data.deltaHRF_rms.Hb_medSNR = sqrt(  mean( (PD_data.Hb_2_medSNR - PD_data.simulationONLY(2).HRF.HbR(:,1)).^2));
    PD_data.deltaHRF_rms.Hb_lowSNR = sqrt(  mean( (PD_data.Hb_2_lowSNR - PD_data.simulationONLY(1).HRF.HbR(:,1)).^2));
    PD_data.deltaHRF_rms.Hb_highSNR = sqrt(  mean( (PD_data.Hb_2_highSNR - PD_data.simulationONLY(3).HRF.HbR(:,1)).^2));

    PD_data.deltaHRF_rms.Hb_noc_medSNR = sqrt(  (mean(PD_data.Hb_2noc_medSNR - PD_data.simulationONLY(2).HRF.HbR(:,1)).^2));
    PD_data.deltaHRF_rms.Hb_noc_lowSNR = sqrt(  (mean(PD_data.Hb_2noc_lowSNR - PD_data.simulationONLY(1).HRF.HbR(:,1)).^2));
    PD_data.deltaHRF_rms.Hb_noc_highSNR = sqrt(  (mean(PD_data.Hb_2noc_highSNR - PD_data.simulationONLY(3).HRF.HbR(:,1)).^2));

    figure()
    semilogy(ch_idx(1:end/2),PD_data.deltaHRF_rms.HbO_medSNR(:,ch_idx(1:end/2)),'rx')
    hold on
    yline([2.41E-6],'b--')
    ylabel('RMS Error')
    xlabel('Channel N.')
    title('Med SNR 2.4E-6 M')

    max_rs = max ( [PD_data.deltaHRF_rms.HbO_noc_medSNR(:,ch_idx(1:end/2)) PD_data.deltaHRF_rms.HbO_medSNR(:,ch_idx(1:end/2))] );
    figure()
    plot([0 max_rs],[0 max_rs],'b')
    hold on
    plot(PD_data.deltaHRF_rms.HbO_noc_lowSNR(:,ch_idx(1:end/2)),PD_data.deltaHRF_rms.HbO_lowSNR(:,ch_idx(1:end/2)),'rx')
    plot(PD_data.deltaHRF_rms.HbO_noc_medSNR(:,ch_idx(1:end/2)),PD_data.deltaHRF_rms.HbO_medSNR(:,ch_idx(1:end/2)),'bx')
    plot(PD_data.deltaHRF_rms.HbO_noc_highSNR(:,ch_idx(1:end/2)),PD_data.deltaHRF_rms.HbO_highSNR(:,ch_idx(1:end/2)),'kx')
    ylabel('Bandpass filter')
    xlabel('No Correction')
    legend('x=y','lowSNR','medSNR','highSNR')
    %title('Med SNR 2.4E-6 M')

    figure()
    boxplot(PD_data.deltaHRF_rms.HbO_medSNR(:,ch_idx(1:end/2)))
    ylabel('RMS Error')
    xlabel('All Good Channels')
    title('Med SNR 2.4E-6 M')

    RMS_boxplot_data = [PD_data.deltaHRF_rms.HbO_lowSNR(:,ch_idx(1:end/2)) ;PD_data.deltaHRF_rms.HbO_medSNR(:,ch_idx(1:end/2)) ; PD_data.deltaHRF_rms.HbO_highSNR(:,ch_idx(1:end/2)) ]';
    RMS_boxplot_data_noc = [PD_data.deltaHRF_rms.HbO_noc_lowSNR(:,ch_idx(1:end/2)) ;PD_data.deltaHRF_rms.HbO_noc_medSNR(:,ch_idx(1:end/2)) ; PD_data.deltaHRF_rms.HbO_noc_highSNR(:,ch_idx(1:end/2)) ]';
    
    RMS_all.deltaHRF_rms.HbO_medSNR(k,:) = PD_data.deltaHRF_rms.HbO_medSNR;
    RMS_all.deltaHRF_rms.Hb_medSNR(k,:) = PD_data.deltaHRF_rms.Hb_medSNR;
    RMS_all.deltaHRF_rms.HbO_lowSNR(k,:) = PD_data.deltaHRF_rms.HbO_lowSNR;
    RMS_all.deltaHRF_rms.Hb_lowSNR(k,:) = PD_data.deltaHRF_rms.Hb_lowSNR;
    RMS_all.deltaHRF_rms.HbO_highSNR(k,:) = PD_data.deltaHRF_rms.HbO_highSNR;
    RMS_all.deltaHRF_rms.Hb_highSNR(k,:) = PD_data.deltaHRF_rms.Hb_highSNR;

    RMS_all.deltaHRF_rms.HbO_noc_medSNR(k,:) = PD_data.deltaHRF_rms.HbO_noc_medSNR;
    RMS_all.deltaHRF_rms.Hb_noc_medSNR(k,:) = PD_data.deltaHRF_rms.Hb_noc_medSNR;
    RMS_all.deltaHRF_rms.HbO_noc_lowSNR(k,:) = PD_data.deltaHRF_rms.HbO_noc_lowSNR;
    RMS_all.deltaHRF_rms.Hb_noc_lowSNR(k,:) = PD_data.deltaHRF_rms.Hb_noc_lowSNR;
    RMS_all.deltaHRF_rms.HbO_noc_highSNR(k,:) = PD_data.deltaHRF_rms.HbO_noc_highSNR;
    RMS_all.deltaHRF_rms.Hb_noc_highSNR(k,:) = PD_data.deltaHRF_rms.Hb_noc_highSNR;
    
    figure()
    boxplot(([RMS_boxplot_data_noc RMS_boxplot_data]),'Labels',{'Noc Low','Noc Med','Noc High','Low','Med','High'})
    ylabel('RMS Error')
    xlabel('All Good Channels')
    title("PD10 NIRS.20-21 HbO - Noc and method"+k+" ")
end

%%
figure()
boxplot(([RMS_boxplot_data_noc RMS_boxplot_data]),'Labels',{'Noc Low','Noc Med','Noc High','Low','Med','High'})
ylabel('RMS Error')
xlabel('All Good Channels')
title("PD10 NIRS.20-21 HbO - Noc and method"+k+" ")

%%
PD_data_bp_0_01 = PD_data;
PD_data_bp_0_001 = PD_data;
PD_data_bp_0_0001 = PD_data;
PD_data_bp_0_0005 = PD_data;

figure()
plot(PD_data.t, PD_data.HbO_2noc_medSNR(:,13))
hold on
plot(PD_data.t, PD_data.HbO_2_medSNR(:,13),'LineWidth',2)
plot(PD_data.t, PD_data.simulationONLY(2).HRF.HbO(:,1),'LineWidth',2)
xlabel('Time / s')
ylabel('delta HbO / M')
xline([900],'b--')
xline([5100],'b--')
legend('No correction','Bandpass','HRF Sim')
title('Med SNR')
%%

figure()
semilogx([0.01 0.001 0.0005 0.0001],[5.597E-6 4.769E-6 4.3145E-6 1.9072E-6],'r*-')
ylabel('RMS / A.U')
xlabel('Low pass filter / Hz')
title('Median RMS - Med SNR - HbO')
[5.597E-6 4.769E-6 4.3145E-6 1.9072E-6]


%%
y = ch_idx
x = [1 2]

PD_data.deltaHRF_rms.HbO_lowSNR
figure()
imagesc([1 ;2],ch_idx',[PD_data.deltaHRF_rms.HbO_lowSNR ;PD_data.deltaHRF_rms.HbO_noc_lowSNR])
colorbar

figure()
imagesc([PD_data.deltaHRF_rms.HbO_noc_lowSNR(:,ch_idx(1:end/2)) ; PD_data.deltaHRF_rms.HbO_lowSNR(:,ch_idx(1:end/2))] )
colorbar
xlabel('All Good Channels')
