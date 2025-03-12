function PD_data = PD_data_motioncorr_15min(PD_data,methodN,figON)

        
fs = PD_data.fs; %Hz
    %0 = no corr.
    %1 = only BP
    %2 = only wavelet + bp

    %3 = only spline + bp(HmR MA dect.)
    %4 = Spline + Wavelet + bp(HmR MA dect.)

    %5 = only spline + bp (GVDT MA dect.)
    %6 = spline + wavelet + bp(GVDT MA dect.)
    

if methodN > 1

end

         %Do GVDT to find MA trains across entire data.
         PD_data = gvtd_10mincorr(PD_data,1);

       
         %Interpolate across MA trains from GVDT - GET dod_int
         PD_data = MA_inter_GVDT_10mins(PD_data,1);
         %consider factor of STD you add in noise for the interpolation
         

        

        %Do HMR to find spike/BL M.A
        tMotion = 1.0;%0.5; %0.8; % %0.5; %time range in seconds
        tMask = 2; %mark data *- time around m.a as m.a
        SDThresh = 7; %10; %12
        AmpThresh = 0.1; %0.35; %0.5;
        tIncMan = ones(length(PD_data.t ),1); % set it to vectors of ones (this is a vector used to remove manually parts of the data if needed)
        % Motion detection technique. tIncCh is a matrix number of samples x twice n of
        % channels which contains for each channel (column) the information about
        % whether an artifact was present (0s) or not (1s). tInc is a vector which
        % contains information on whether at that time sample in any of the channel
        % was present an artifact (0s) or not (1s). tInc can therefore be obtained
        % from tIncCh by setting to 0 every row that contains at least one 0.
        [tInc,tIncCh] = hmrMotionArtifactByChannel(PD_data.dod_int, fs, PD_data.SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);
        N_MA_HMR = size(find(tIncCh(:,PD_data.goodch_idx)==0),1)/size(PD_data.goodch_idx,1);
        figure()
        plot(PD_data.t,PD_data.dod_int(:,13))
        hold on
        plot(PD_data.t(find(tIncCh(:,13)==0)),PD_data.dod_int(find(tIncCh(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("Hmr YANG 22 for MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ")
        
        

%%  
        if figON == 0
            set(gcf,'Visible','off');              
            set(0,'DefaultFigureVisible','off');

        elseif figON == 1
            set(gcf,'Visible','on');              
            set(0,'DefaultFigureVisible','on');
        end

        p = 0.99; %0.99 - GET dod_spline
        PD_data.dod_spline = hmrMotionCorrectSpline(PD_data.dod_int,PD_data.t,PD_data.SD,tIncCh,p);

        figure()
        plot(PD_data.t,PD_data.dod_int(:,13))
        hold on
        plot(PD_data.t,PD_data.dod_spline(:,13))
        plot(PD_data.t(find(tIncCh(:,13)==0)),PD_data.dod_int(find(tIncCh(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("Use Spline on Hmr YANG 22 MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ")
    %%
        %test Spline SG - GET dod_SG
        p=0.99;

        %do SG on dodINT
        %FrameSize_sec = 1;%0.5 %1/fs;
        FrameSize_sec = 10;



        [dod_SG ,tIncCh_baseline_dod_SG,tInc_baseline_dod_SG] = hmrMotionCorrectSplineSG_PD_data(PD_data.dod_int, PD_data.dod_int, PD_data.t, PD_data.SD, p, FrameSize_sec,1,PD_data);
        PD_data.dod_SG = dod_SG;
        
        %do SG on dodSPline
        %FrameSize_sec = 1;%0.5 %1/fs;
        [dod_SG_s ,tIncCh_baseline_dod_SG_s,tInc_baseline_dod_SG_s] = hmrMotionCorrectSplineSG(PD_data.dod_spline, PD_data.dod_int, PD_data.t, PD_data.SD, p, FrameSize_sec,1);
        PD_data.dod_SG_s = dod_SG_s;
        
        %do SG on dodSG to make SG 2
        FrameSize_sec = 6;
        [dod_SG_2 ,tIncCh_baseline_dod_SG_2,tInc_baseline_dod_SG_2] = hmrMotionCorrectSplineSG(PD_data.dod_SG, PD_data.dod_int, PD_data.t, PD_data.SD, p, FrameSize_sec,1);
        PD_data.dod_SG_2 = dod_SG_2;
    
        figure()
        plot(PD_data.t,PD_data.dod_int(:,13))
        hold on
        plot(PD_data.t,PD_data.dod_SG_2(:,13))
        plot(PD_data.t(find(tIncCh_baseline_dod_SG_2(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG_2(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("Use Spline SG Hmr YANG 22 MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ")

        figure()
        plot(PD_data.t,PD_data.dod_int(:,13))
        hold on
        plot(PD_data.t,PD_data.dod_spline(:,13))
        plot(PD_data.t,PD_data.dod_SG_s(:,13))
        plot(PD_data.t,PD_data.dod_SG(:,13))
        plot(PD_data.t,PD_data.dod_SG_2(:,13))
        %plot(PD_data.t(find(tIncCh_baseline_dod_SG(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("dod int vs Spline vs vs Spline SG(s) vs SplineSG vs SGx2 Hmr YANG 22 MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ")
        legend('Int. dod','Spline','Spline SG (spline)','Spline SG','Spline SG x2')

        figure()
        subplot(3,1,1)
        plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(1) ))
        hold on
        plot(PD_data.t,PD_data.dod_spline(:,PD_data.goodch_idx(1)))
        plot(PD_data.t,PD_data.dod_SG_s(:,PD_data.goodch_idx(1)))
        plot(PD_data.t,PD_data.dod_SG(:,PD_data.goodch_idx(1)))
        plot(PD_data.t,PD_data.dod_SG_2(:,PD_data.goodch_idx(1)))
        %plot(PD_data.t(find(tIncCh_baseline_dod_SG(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("Ch "+num2str(PD_data.goodch_idx(1))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        %title("dod int vs Spline vs vs Spline SG(s) vs SplineSG vs SGx2 Hmr YANG 22 MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ch "+num2str(PD_data.goodch_idx(1))+"")
        legend('Int. dod','Spline','Spline SG (spline)','Spline SG','Spline SG x2')
        subplot(3,1,2)
        plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(2) ))
        hold on
        plot(PD_data.t,PD_data.dod_spline(:,PD_data.goodch_idx(2)))
        plot(PD_data.t,PD_data.dod_SG_s(:,PD_data.goodch_idx(2)))
        plot(PD_data.t,PD_data.dod_SG(:,PD_data.goodch_idx(2)))
        plot(PD_data.t,PD_data.dod_SG_2(:,PD_data.goodch_idx(2)))
        %plot(PD_data.t(find(tIncCh_baseline_dod_SG(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("Ch "+num2str(PD_data.goodch_idx(2))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        %title("dod int vs Spline vs vs Spline SG(s) vs SplineSG vs SGx2 Hmr YANG 22 MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ch "+num2str(PD_data.goodch_idx(2))+"")
        legend('Int. dod','Spline','Spline SG (spline)','Spline SG','Spline SG x2')
        subplot(3,1,3)
        plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(11) ))
        hold on
        plot(PD_data.t,PD_data.dod_spline(:,PD_data.goodch_idx(11)))
        plot(PD_data.t,PD_data.dod_SG_s(:,PD_data.goodch_idx(11)))
        plot(PD_data.t,PD_data.dod_SG(:,PD_data.goodch_idx(11)))
        plot(PD_data.t,PD_data.dod_SG_2(:,PD_data.goodch_idx(11)))
        %plot(PD_data.t(find(tIncCh_baseline_dod_SG(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("Ch "+num2str(PD_data.goodch_idx(11))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        %title("dod int vs Spline vs vs Spline SG(s) vs SplineSG vs SGx2 Hmr YANG 22 MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ch "+num2str(PD_data.goodch_idx(11))+"")
        legend('Int. dod','Spline','Spline SG (spline)','Spline SG','Spline SG x2')


        figure()
        plot(PD_data.t,PD_data.dod(:,13))
        hold on
        plot(PD_data.t,PD_data.dod_int(:,13))
        plot(PD_data.t,PD_data.dod_spline(:,13))
        plot(PD_data.t,PD_data.dod_SG(:,13))
        plot(PD_data.t,PD_data.dod_SG_2(:,13))
        %plot(PD_data.t(find(tIncCh_baseline_dod_SG(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("dod vs dod int vs Spline vs SplineSG SplineSG x2  Hmr YANG 22 MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ")
        legend('dod','Int. dod','Spline','Spline SG','Spline SG x2')

        %Compare different window sizes SG

        %%
        %do wavelet        
        iqr = 1.1;
        %PD_data.dod_wav_spline = hmrMotionCorrectWavelet(PD_data.dod_spline,PD_data.SD,iqr); %run wavelet
        %PD_data.dod_wav_SG_s = hmrMotionCorrectWavelet(PD_data.dod_SG_s,PD_data.SD,iqr); %run wavelet
        PD_data.dod_wav_SG = hmrMotionCorrectWavelet(PD_data.dod_SG,PD_data.SD,iqr); %run wavelet 
        %PD_data.dod_wav_SG2 = hmrMotionCorrectWavelet(PD_data.dod_SG_2,PD_data.SD,iqr); %run wavelet 


        % figure()
        % subplot(3,1,1)
        % plot(PD_data.t,PD_data.dod_wav_spline(:,PD_data.goodch_idx(1)))
        % hold on
        % plot(PD_data.t,PD_data.dod_wav_SG_s(:,PD_data.goodch_idx(1)))
        % plot(PD_data.t,PD_data.dod_wav_SG(:,PD_data.goodch_idx(1)))
        % plot(PD_data.t,PD_data.dod_wav_SG2(:,PD_data.goodch_idx(1)))
        % plot(PD_data.t,PD_data.dod_SG_2(:,PD_data.goodch_idx(1)))
        % plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(1)))
        % xlabel('Time / s')
        % ylabel('dOD / A.U')
        % legend('dod wav (spline)','dod wav (SG(spline))','dod wav(SG)','dod wav(SGx2)','dod (SGx2)','dod int')
        % title("Ch "+num2str(PD_data.goodch_idx(1))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        % subplot(3,1,2)
        % plot(PD_data.t,PD_data.dod_wav_spline(:,PD_data.goodch_idx(2)))
        % hold on
        % plot(PD_data.t,PD_data.dod_wav_SG_s(:,PD_data.goodch_idx(2)))
        % plot(PD_data.t,PD_data.dod_wav_SG(:,PD_data.goodch_idx(2)))
        % plot(PD_data.t,PD_data.dod_wav_SG2(:,PD_data.goodch_idx(2)))
        % plot(PD_data.t,PD_data.dod_SG_2(:,PD_data.goodch_idx(2)))
        % plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(2)))
        % xlabel('Time / s')
        % ylabel('dOD / A.U')
        % legend('dod wav (spline)','dod wav (SG(spline))','dod wav(SG)','dod wav(SGx2)','dod SGx2)','dod int')
        % title("Ch "+num2str(PD_data.goodch_idx(2))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        % subplot(3,1,3)
        % plot(PD_data.t,PD_data.dod_wav_spline(:,PD_data.goodch_idx(11)))
        % hold on
        % plot(PD_data.t,PD_data.dod_wav_SG_s(:,PD_data.goodch_idx(11)))
        % plot(PD_data.t,PD_data.dod_wav_SG(:,PD_data.goodch_idx(11)))
        % plot(PD_data.t,PD_data.dod_wav_SG2(:,PD_data.goodch_idx(11)))
        % plot(PD_data.t,PD_data.dod_SG_2(:,PD_data.goodch_idx(11)))
        % plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(11)))
        % xlabel('Time / s')
        % ylabel('dOD / A.U')
        % legend('dod wav (spline)','dod wav (SG(spline))','dod wav(SG)','dod wav(SGx2)','dod SGx2)','dod int')
        % title("Ch "+num2str(PD_data.goodch_idx(11))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        % sgtitle("wavlet comparisons CH")    

         %Correct these
         

    %%
        % Bandpass Filter
        lowerCutOff = 0;
        %higherCutOff = 0.01;
        %higherCutOff = 0.05; %pre 02 05 24
        higherCutOff = 0.0067;
        %0.001

        

        fs = 10; %sampling rate (Hz)


        PD_data.BP_SG2 = hmrBandpassFilt(PD_data.dod_SG_2,fs,lowerCutOff,higherCutOff);
        %PD_data.BP_wav_SG2 = hmrBandpassFilt(PD_data.dod_wav_SG2,fs,lowerCutOff,higherCutOff);
        PD_data.dod1_noc = hmrBandpassFilt(PD_data.dod,fs,lowerCutOff,higherCutOff);

        %
        %PD_data.BP_wav_spline = hmrBandpassFilt(PD_data.dod_wav_spline,fs,lowerCutOff,higherCutOff);
        %PD_data.BP_wav_SG_s = hmrBandpassFilt(PD_data.dod_wav_SG_s,fs,lowerCutOff,higherCutOff);
        PD_data.BP_wav_SG = hmrBandpassFilt(PD_data.dod_wav_SG,fs,lowerCutOff,higherCutOff);

        %choose dod1
        %PD_data.dod1 = hmrBandpassFilt(PD_data.dod_SG_2,fs,lowerCutOff,higherCutOff);
        PD_data.dod1 = hmrBandpassFilt(PD_data.dod_wav_SG,fs,lowerCutOff,higherCutOff);


        figure()
        subplot(3,1,1)
        plot(PD_data.t,PD_data.dod1_noc(:,PD_data.goodch_idx(1)))
        hold on
        plot(PD_data.t,PD_data.BP_SG2(:,PD_data.goodch_idx(1)))
        plot(PD_data.t,PD_data.BP_wav_SG(:,PD_data.goodch_idx(1)))
        %%%%plot(PD_data.t(find(tIncCh_baseline_dod_SG(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        legend('BP no corr','BP SGx2','BP Wav(SG)')
        title("Ch "+num2str(PD_data.goodch_idx(1))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        subplot(3,1,2)
        plot(PD_data.t,PD_data.dod1_noc(:,PD_data.goodch_idx(2)))
        hold on
        plot(PD_data.t,PD_data.BP_SG2(:,PD_data.goodch_idx(2)))
        plot(PD_data.t,PD_data.BP_wav_SG(:,PD_data.goodch_idx(2)))
        %%%%plot(PD_data.t(find(tIncCh_baseline_dod_SG(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        legend('BP no corr','BP SGx2','BP Wav(SG)')
        title("Ch "+num2str(PD_data.goodch_idx(2))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        subplot(3,1,3)
        plot(PD_data.t,PD_data.dod1_noc(:,PD_data.goodch_idx(11)))
        hold on
        plot(PD_data.t,PD_data.BP_SG2(:,PD_data.goodch_idx(11)))
        plot(PD_data.t,PD_data.BP_wav_SG(:,PD_data.goodch_idx(11)))
        %%%%%plot(PD_data.t(find(tIncCh_baseline_dod_SG(:,13)==0)),PD_data.dod_int(find(tIncCh_baseline_dod_SG(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        legend('BP no corr','BP SGx2','BP Wav(SG)')
        title("Ch "+num2str(PD_data.goodch_idx(11))+" Subject"+num2str(PD_data.subjectN)+ " Event "+PD_data.eventType+" Event N "+num2str(PD_data.eventN)+" Time W"+PD_data.time_window+"")
        sgtitle("BP comparisons CH")   
    
        figure()
        plot(PD_data.t,PD_data.dod_SG_2(:,13))
        hold on
        plot(PD_data.t,PD_data.dod1(:,13))
        legend('SGx2 dod','Bandpass on SGx2 dod')
        title("BP lower "+lowerCutOff+" Hz and Higher "+higherCutOff+" Hz");

        
        %%
        % PD_data.dod1 = hmrBandpassFilt(PD_data.dod_wav_SG2,fs,lowerCutOff,higherCutOff);
        % figure()
        % plot(PD_data.t,PD_data.dod_SG_2(:,13))
        % hold on
        % plot(PD_data.t,PD_data.dod1(:,13))
        % legend('SGx2 dod','Bandpass on wav(SGx2) dod')
        % title("BP lower "+lowerCutOff+" Hz and Higher "+higherCutOff+" Hz");


%PD_data.dod1_noc = PD_data.dod1;

end