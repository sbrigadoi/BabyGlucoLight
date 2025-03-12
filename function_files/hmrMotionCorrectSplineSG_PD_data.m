function [dod ,tIncCh_baseline_dod_SG,tInc_baseline_dod_SG SNR_Three_PD]= hmrMotionCorrectSplineSG_PD_data(dod, d, t, SD, p, FrameSize_sec, turnon,PD_data)
% Sahar Jahani, October 2017
% Added turn on/off option Meryem Nov 2017
if exist('turnon')
   if turnon==0
   return;
   end
end

[tIncCh, tInc SNR_Three_PD] = hmrtInc_baselineshift_Ch_PD_data(dod, t,PD_data); % finding the baseline shift motions

tIncCh_baseline_dod_SG = tIncCh;
tInc_baseline_dod_SG = tInc;

figure()
plot(t,dod(:,1))
hold on
plot(t(find(tIncCh(:,1)==0)),dod(find(tIncCh(:,1)==0),1),'m.','MarkerSize',6)
xlabel('Time / s')
ylabel('dOD / A.U')
title("BL changes found in SG")
%title("MedSNR Hmr MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR_medSNR+" ")

fs = abs(1/(t(2)-t(1)));
%% extending signal for motion detection purpose (12 sec from each edge)
extend = round(12*fs); 

tIncCh1=repmat(tIncCh(1,:),extend,1);
tIncCh2=repmat(tIncCh(end,:),extend,1);
tIncCh=[tIncCh1;tIncCh;tIncCh2];

d1=repmat(dod(1,:),extend,1);
d2=repmat(dod(end,:),extend,1);
dod=[d1;dod;d2];

t2=(0:(1/fs):(length(dod)/fs))';
t2=t2(1:length(dod),1);

[dodLP,ylpf] = hmrBandpassFilt( dod, fs, 0, 2 ); 

%% Spline Interpolation
dod = hmrMotionCorrectSpline(dodLP, t2, SD, tIncCh, p);
dod=dod(extend+1:end-extend,:); % removing the extention

%% Savitzky_Golay filter
K = 3; % polynomial order
FrameSize_sec = round(FrameSize_sec * fs);
if mod(FrameSize_sec,2)==0
    FrameSize_sec = FrameSize_sec  + 1;
end
dod=sgolayfilt(dod,K,FrameSize_sec);

