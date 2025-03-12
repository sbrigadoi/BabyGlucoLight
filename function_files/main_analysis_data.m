%restoredefaultpath
%rehash toolboxcache

clear all;
close all;
clc;
%save('PD10_m_hypo_1.nirs','PD_data','-mat')
%load('PD10_m_hypo_1.nirs','-mat')
%% Step 1. Load NIRS data and Glucose Data

%Subject and event parameters
subject = 'PD10';
eventT = 'm_hypo';
eventN = 1;
storage = 'E'; %E for HDD with PC, F for HDD on Laptop;
fs = 10; %sampling rate (Hz)

%load data
PD_data = loadNIRS_glucose_data(subject,eventT,eventN,storage);
%% Step 2. Signal Qaulity check

dRange = [5E-4 3];
SNRrange = 1 %2;
[goodch_idx PD_data] = find_good_ch(PD_data,dRange,SNRrange);

%% Step 3. Find M.A Train (1=no MA Train, 0 =MA train) (1= no MA, 0=MA)
%(if ch is in MA train idx, it doesn't have a MA train)
meanValue = mean(PD_data.d);
dodConv = -log(abs(PD_data.d)./meanValue);
dodConv_og = -log(abs(PD_data.d)./meanValue);

[MA_train_idx  tIncCh] = find_MA_train(PD_data,dodConv,goodch_idx);
% Step 3.5 Show channels
ch_idx = goodch_idx;
show_data(PD_data,dodConv_og,ch_idx)

%% Step 4. Do Motion Corr. 
%A. Do spline on ch without MA train
p = 0.99; %0.99
dodSpline = hmrMotionCorrectSpline(dodConv,PD_data.t,PD_data.SD,tIncCh,p);
dodConv(:,MA_train_idx) = dodSpline(:,MA_train_idx); 
show_data(PD_data,dodConv,ch_idx)
%% Step 4.5 WAV Motion Corr.
%iqr = 1.1;
%dodWavelet = hmrMotionCorrectWavelet(dodConv,PD_data.SD,iqr); %run wavelet with no spline
%show_data(PD_data,dodWavelet,ch_idx)

%% Step 4.X Bandpass Filter

lowerCutOff = 0;
higherCutOff = 0.01;
fs; % compute sampling frequency
dodFilt = hmrBandpassFilt(dodConv,fs,lowerCutOff,higherCutOff);

%% Step 4.Y block avg/baseline
dod = dodFilt; %this should be the final DOD form of the processed data.
PD_data.dod = dod;
%% Step 5. Downsample
% Reshape and down-sample data to 0.1 Hz.
ds = fs*60; %1 min block avg (10frames per second * 60s = N frames in 1 min)
PD_data = downsample_data(PD_data,ds);

%% Step 6. Block Avg.
%Use first 5 minutes as baseline to look at DOD across all time %USE THIS
baseline = [1 5]; %start and end of baseline
PD_data = block_avg_data(PD_data,baseline);

%% Step 7. Recover HbO Hb HbT
weekN = 30; %week age of infant
[HbO_c,Hb_c,HbT_c] = recon_HbO_Hb(PD_data,weekN,goodch_idx);

%% View Tomography
view_HbO_Hb_tomography(HbO_c,Hb_c,HbT_c,goodch_idx,weekN)

%B. Do Wav on All

%C. Do we find MA
%if yes, do above 2 steps
%if no proceed








