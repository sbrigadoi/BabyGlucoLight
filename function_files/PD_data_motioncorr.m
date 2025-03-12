function PD_data = PD_data_motioncorr(PD_data,methodN)
    %0 = no corr.
    %1 = only BP
    %2 = only wavelet + bp

    %3 = only spline + bp(HmR MA dect.)
    %4 = Spline + Wavelet + bp(HmR MA dect.)

    %5 = only spline + bp (GVDT MA dect.)
    %6 = spline + wavelet + bp(GVDT MA dect.)
dRange = [5E-4 3];
SNRrange = 2; %2;
[goodch_idx, PD_data] = find_good_ch(PD_data,dRange,SNRrange);

fs = 10;
length_train = 600; %frames 
PD_data.dod2_og_lowSNR = PD_data.dod2_lowSNR;
PD_data.dod2_og_medSNR = PD_data.dod2_medSNR;
PD_data.dod2_og_highSNR = PD_data.dod2_highSNR;

if methodN > 1
    if methodN == 2
        iqr = 1.1;
        PD_data.dod2_medSNR = hmrMotionCorrectWavelet(PD_data.dod2_medSNR,PD_data.SD,iqr); %run wavelet with no spline
        PD_data.dod2_lowSNR = hmrMotionCorrectWavelet(PD_data.dod2_lowSNR,PD_data.SD,iqr); %run wavelet with no spline
        PD_data.dod2_highSNR = hmrMotionCorrectWavelet(PD_data.dod2_highSNR,PD_data.SD,iqr); %run wavelet with no spline
        
        %PD_data.dod2_medSNR = dodWavelet_medSNR;
        %PD_data.dod2_lowSNR = dodWavelet_lowSNR;
        %PD_data.dod2_highSNR = dodWavelet_highSNR;
        'wavelet'
    end
    
    if methodN == 3 | methodN == 4
        'Homer MA detect'
        % Detect motion artifacts in signal
        tMotion = 1;%0.5; %0.8; % %0.5; %time range in seconds
        tMask = 2; %mark data *- time around m.a as m.a
        SDThresh = 12; %10; %12
        AmpThresh = 0.15; %0.35; %0.5;
        tIncMan = ones(length(PD_data.t ),1); % set it to vectors of ones (this is a vector used to remove manually parts of the data if needed)
        % Motion detection technique. tIncCh is a matrix number of samples x twice n of
        % channels which contains for each channel (column) the information about
        % whether an artifact was present (0s) or not (1s). tInc is a vector which
        % contains information on whether at that time sample in any of the channel
        % was present an artifact (0s) or not (1s). tInc can therefore be obtained
        % from tIncCh by setting to 0 every row that contains at least one 0.
        [tInc,tIncCh_medSNR] = hmrMotionArtifactByChannel(PD_data.dod2_medSNR, fs, PD_data.SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);
        [tInc,tIncCh_lowSNR] = hmrMotionArtifactByChannel(PD_data.dod2_lowSNR, fs, PD_data.SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);
        [tInc,tIncCh_highSNR] = hmrMotionArtifactByChannel(PD_data.dod2_highSNR, fs, PD_data.SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);
        
        N_MA_HMR_medSNR = size(find(tIncCh_medSNR(:,goodch_idx)==0),1)/size(goodch_idx,1);
        figure()
        plot(PD_data.t,PD_data.dod2_medSNR(:,13))
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),PD_data.dod2_medSNR(find(tIncCh_medSNR(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("MedSNR Hmr MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR_medSNR+" ")

       
    

        %find M.A train
        %IncCh_train_medSNR = ones(1,size(tIncCh_medSNR,2));
        %IncCh_train_loc_medSNR = ones(size(tIncCh_medSNR,1),size(tIncCh_medSNR,2));
        %for i=1:size(tIncCh_medSNR,2)
        %    if size(unique(tIncCh_medSNR(:,i)),1) > 1
        %        for j=1: size(tIncCh_medSNR,1)-length_train
        %            if nnz(tIncCh_medSNR(j:j+(length_train-1),i)) < (length_train/2)
        %                IncCh_train_medSNR(1,i) = 0;
        %                IncCh_train_loc_medSNR(j,i) = 0;
        %            end
        %        end
        %    end
        %end
        %IncCh_train_medSNR = IncCh_train_medSNR';
        %MA_train_idx_medSNR = IncCh_train_medSNR;
        %for i=1:size(MA_train_idx_medSNR,1)
        %    if MA_train_idx_medSNR(i,1) == 1;
        %        MA_train_idx_medSNR(i,1) = i;
        %    end
        %end
        %MA_train_idx_medSNR = nonzeros(MA_train_idx_medSNR);

        %IncCh_train_lowSNR = ones(1,size(tIncCh_lowSNR,2));
        %IncCh_train_loc_lowSNR = ones(size(tIncCh_lowSNR,1),size(tIncCh_lowSNR,2));
        %for i=1:size(tIncCh_lowSNR,2)
        %    if size(unique(tIncCh_lowSNR(:,i)),1) > 1
        %        for j=1: size(tIncCh_lowSNR,1)-length_train
        %            if nnz(tIncCh_lowSNR(j:j+(length_train-1),i)) < (length_train/2)
        %                IncCh_train_lowSNR(1,i) = 0;
        %                IncCh_train_loc_lowSNR(j,i) = 0;
        %            end
        %        end
        %    end
        %end
        %IncCh_train_lowSNR = IncCh_train_lowSNR';
        %MA_train_idx_lowSNR = IncCh_train_lowSNR;
        %for i=1:size(MA_train_idx_lowSNR,1)
        %    if MA_train_idx_lowSNR(i,1) == 1;
        %        MA_train_idx_lowSNR(i,1) = i;
        %    end
        %end
        %MA_train_idx_lowSNR = nonzeros(MA_train_idx_lowSNR);

        %IncCh_train_highSNR = ones(1,size(tIncCh_highSNR,2));
        %IncCh_train_loc_highSNR = ones(size(tIncCh_highSNR,1),size(tIncCh_highSNR,2));
        %for i=1:size(tIncCh_highSNR,2)
        %    if size(unique(tIncCh_highSNR(:,i)),1) > 1
        %        for j=1: size(tIncCh_highSNR,1)-length_train
        %            if nnz(tIncCh_highSNR(j:j+(length_train-1),i)) < (length_train/2)
        %                IncCh_train_highSNR(1,i) = 0;
        %                IncCh_train_loc_highSNR(j,i) = 0;
        %            end
        %        end
        %    end
        %end
        %IncCh_train_highSNR = IncCh_train_highSNR';
        %MA_train_idx_highSNR = IncCh_train_highSNR;
        %for i=1:size(MA_train_idx_highSNR,1)
        %    if MA_train_idx_highSNR(i,1) == 1;
        %        MA_train_idx_highSNR(i,1) = i;
        %    end
        %end
        %MA_train_idx_highSNR = nonzeros(MA_train_idx_highSNR);
        
        %%%1 = No train
        %%%0 = MA Train present in channel
        %for i=1:size(PD_data.dod,2)
        %    if size(unique(IncCh_train_loc_medSNR(:,i)),1) ==2 %if channel has MA train
        %        MA_locations_true_medSNR = [find(tIncCh_medSNR(:,i)==0)];
        %        MA_train_locations_true_medSNR = find(IncCh_train_loc_medSNR(:,i)==0);
        %    end
        %    if size(unique(IncCh_train_loc_lowSNR(:,i)),1) ==2 %if channel has MA train
        %        MA_locations_true_lowSNR = [find(tIncCh_lowSNR(:,i)==0)];
        %        MA_train_locations_true_lowSNR = find(IncCh_train_loc_lowSNR(:,i)==0);
        %    end
        %    if size(unique(IncCh_train_loc_highSNR(:,i)),1) ==2 %if channel has MA train
        %        MA_locations_true_highSNR = [find(tIncCh_highSNR(:,i)==0)];
        %        MA_train_locations_true_highSNR = find(IncCh_train_loc_highSNR(:,i)==0);
        %    end
        %end
        %%% MA terain end

         %figure()
         %plot(PD_data.t,PD_data.dod2_medSNR(:,13))
         %hold on
         %plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),PD_data.dod2_medSNR(find(tIncCh_medSNR(:,13)==0),13),'m.','MarkerSize',6)
         %plot(PD_data.t(find(IncCh_train_loc_medSNR(:,13)==0)),PD_data.dod2_medSNR(find(IncCh_train_loc_medSNR(:,13)==0),13),'y.','MarkerSize',6)
         %xlabel('Time / s')
         %ylabel('dOD / A.U')
        

        
        % Spline
        %A. Do spline on ch without MA train
        p = 0.99; %0.99

        PD_data.dod2_medSNR = hmrMotionCorrectSpline(PD_data.dod2_medSNR,PD_data.t,PD_data.SD,tIncCh_medSNR,p);
        %PD_data.dod2_lowSNR = hmrMotionCorrectSpline(PD_data.dod2_lowSNR,PD_data.t,PD_data.SD,tIncCh_lowSNR,p);
        %PD_data.dod2_highSNR = hmrMotionCorrectSpline(PD_data.dod2_highSNR,PD_data.t,PD_data.SD,tIncCh_highSNR,p);

        figure()
        plot(PD_data.t,PD_data.dod2_og_medSNR(:,13))
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),PD_data.dod2_og_medSNR(find(tIncCh_medSNR(:,13)==0),13),'m.','MarkerSize',6)
        plot(PD_data.t,PD_data.dod2_medSNR(:,13))
        xlabel('Time / s')
        ylabel('dOD / A.U')
        %title("MedSNR Hmr MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR_medSNR+" ")

        'Spline'

        if methodN == 4 
            iqr = 1.1;
            PD_data.dod2_medSNR = hmrMotionCorrectWavelet(PD_data.dod2_medSNR,PD_data.SD,iqr); %run wavelet 
            PD_data.dod2_lowSNR = hmrMotionCorrectWavelet(PD_data.dod2_lowSNR,PD_data.SD,iqr); %run wavelet 
            PD_data.dod2_highSNR = hmrMotionCorrectWavelet(PD_data.dod2_highSNR,PD_data.SD,iqr); %run wavelet 
            'Wavelet'
        end
        % figure()
        % plot(PD_data.t,PD_data.dod2_medSNR(:,13))
        % hold on
        % plot(PD_data.t,dodSpline_medSNR(:,13))

        %PD_data.dod2Spline_lowSNR = dodSpline_lowSNR;
        %PD_data.dod2Spline_medSNR = dodSpline_medSNR;
        %PD_data.dod2Spline_highSNR = dodSpline_highSNR;

        %PD_data.dod2_lowSNR

        %PD_data.dod3_lowSNR(:,MA_train_idx) = dodSpline_lowSNR(:,MA_train_idx);
        %PD_data.dod3_medSNR(:,MA_train_idx) = dodSpline_medSNR(:,MA_train_idx);
        %PD_data.dod3_highSNR(:,MA_train_idx) = dodSpline_highSNR(:,MA_train_idx);

        %PD_data.dod3_lowSNR(:,MA_train_idx) = dodSpline_lowSNR(:,MA_train_idx);
        %PD_data.dod3_medSNR(:,MA_train_idx) = dodSpline_medSNR(:,MA_train_idx);
        %PD_data.dod3_highSNR(:,MA_train_idx) = dodSpline_highSNR(:,MA_train_idx);
        %use below %%%%%%%%%%%
        %PD_data.dod2_lowSNR(:,MA_train_idx_lowSNR) = dodSpline_lowSNR(:,MA_train_idx_lowSNR);
        %PD_data.dod2_medSNR(:,MA_train_idx_medSNR) = dodSpline_medSNR(:,MA_train_idx_medSNR);
        %PD_data.dod2_highSNR(:,MA_train_idx_highSNR) = dodSpline_highSNR(:,MA_train_idx_highSNR);

        %%show_data(PD_data,dodConv,ch_idx)
        % Step 4.5 WAV Motion Corr.
        %iqr = 1.1;
        %dodWavelet = hmrMotionCorrectWavelet(PD_data.dod2_medSNR,PD_data.SD,iqr); %run wavelet with no spline
        %show_data(PD_data,dodWavelet,ch_idx)
        %PD_data.dod2splinewavelet = dodWavelet;
        %PD_data.dod2 = dodWavelet;
    end

    if methodN == 5 | methodN == 6
        %do GVDT
        'GVDT'
        %GVDT
         %for i=2:size(PD_data.dod2_medSNR,1)
         %    g_medSNR(i,:)=sqrt( (PD_data.dod2_medSNR(i,:)-PD_data.dod2_medSNR(i-1,:)).^2 );
         %end

        g_medSNR = zeros(size(PD_data.dod2_medSNR,1),size(PD_data.dod2_medSNR,2));
         for i=2:size(PD_data.dod2_medSNR,1)
                for j=1:size(PD_data.dod2_medSNR,2)
                    g_medSNR(i,j)=sqrt( (PD_data.dod2_medSNR(i,j)-PD_data.dod2_medSNR(i-1,j)).^2 );
                end
         end

        %abr = mean(g);
        %gvtdTimeTrace = gvtd(PD_data.dod2_medSNR'); %input ch x t, output t x 1
        %use only 'good channels' idx for GVDT distribution
        gvtdTimeTrace_medSNR = gvtd(PD_data.dod2_medSNR(:,goodch_idx)'); %input ch x t, output t x 1
        gvtdTimeTrace_lowSNR = gvtd(PD_data.dod2_lowSNR(:,goodch_idx)'); %input ch x t, output t x 1
        gvtdTimeTrace_highSNR = gvtd(PD_data.dod2_highSNR(:,goodch_idx)'); %input ch x t, output t x 1

        %figure()
        %plot(PD_data.t,gvtdTimeTrace_medSNR)
        %ylabel('GVTD (global RMS)/ A.U')
        %xlabel('Time / s')
        %title('GVTD of good channels')

        statType = StatType.Histogram_Mode;
        nStd = 4;
        thresh_medSNR = find_gvtd_thresh(gvtdTimeTrace_medSNR, statType, nStd);
        thresh_lowSNR = find_gvtd_thresh(gvtdTimeTrace_lowSNR, statType, nStd);
        thresh_highSNR = find_gvtd_thresh(gvtdTimeTrace_highSNR, statType, nStd);


        gvdt_more_thresh_medSNR = find(gvtdTimeTrace_medSNR>thresh_medSNR);
        N_MA_GVDT_medSNR = size(gvdt_more_thresh_medSNR,1);
        
        gvdt_more_thresh_lowSNR = find(gvtdTimeTrace_lowSNR>thresh_lowSNR);
        N_MA_GVDT_lowSNR = size(gvdt_more_thresh_lowSNR,1);
        
        gvdt_more_thresh_highSNR = find(gvtdTimeTrace_highSNR>thresh_highSNR);
        N_MA_GVDT_highSNR = size(gvdt_more_thresh_highSNR,1);
    
        %thresh = make_gvtd_hist(gvtdTimeTrace, plotThresh, statType, nStd, binSize)
        
        figure()
        histogram(gvtdTimeTrace_medSNR)
        hold on
        xline([thresh_medSNR],LineWidth=4)
        xlabel('GVTD value / A.U')
        ylabel('Frequency / N. of timepoints')
        title("Histogram GVTD nSTD = "+nStd+" yields threshold of "+thresh_medSNR+" with "+N_MA_GVDT_medSNR+" M.As")

        % Spline
        %A. Do spline on ch without MA train
        p = 0.99; %0.99
        
        %set MA markers. 1= NO MA, 0 = MA
        tIncCh_medSNR = ones(size(PD_data.dod2_medSNR,1)  , size(PD_data.dod2_medSNR,2)  );
        tIncCh_medSNR(gvdt_more_thresh_medSNR,:) = 0;

        N_MA_HMR_medSNR = size(find(tIncCh_medSNR(:,goodch_idx)==0),1)/size(goodch_idx,1);
        figure()
        plot(PD_data.t,PD_data.dod2_medSNR(:,13))
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),PD_data.dod2_medSNR(find(tIncCh_medSNR(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("MedSNR Hmr MAD Using GVDT nStd = "+nStd+" |N M.A = "+N_MA_HMR_medSNR+" ")

         figure()
         sgtitle("Med SNR GVTD lambda "+nStd+" ")
         for i=1:12
            subplot(3,4,i)
            plot(PD_data.t, PD_data.dod2_medSNR(:,goodch_idx(i)))
            hold on
            plot(PD_data.t(find(tIncCh_medSNR(:,goodch_idx(i))==0)),PD_data.dod2_medSNR(find(tIncCh_medSNR(:,goodch_idx(i))==0),goodch_idx(i)),'m.','MarkerSize',6)
            xline([900],'b--')
            xline([5100],'b--')
            title("ch "+num2str(goodch_idx(i)))
            xlabel('T / s')
            ylabel('dod / AU')
          end

        %gvtd_medfilt_medSNR = medfilt1(gvtdTimeTrace_medSNR,10);
        figure()
        subplot(2,1,1)
        imagesc(g_medSNR(:,goodch_idx)')
        colorbar, xlabel('Time / s'), ylabel('Measurement #'),title("GVDT GoodCh medSNR")
        subplot(2,1,2)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR(:,13)==0),1),'m.','MarkerSize',6)
         
        avg_g_medSNR = mean(g_medSNR(:,goodch_idx));
        nStd_meanGVDT =1;
        threshold_meanGVDT = [mean(avg_g_medSNR) + (nStd_meanGVDT*std(avg_g_medSNR))];
        figure()
        plot(mean(g_medSNR(:,goodch_idx)),'r+');
        hold on
        yline([mean(avg_g_medSNR) + (nStd_meanGVDT*std(avg_g_medSNR))],'b--')
        ylabel('Mean GVDT Value')
        xlabel('Good Ch N')
        title("Mean GVDT good ch, threshold mean +"+nStd_meanGVDT+"x STD ("+threshold_meanGVDT+")")

        %NOW recalcute 'goodch idx' by rejecting channels above threshold
        %in GVDT
        goodch_idx = goodch_idx(find(avg_g_medSNR < threshold_meanGVDT));

        PD_data.goodch_idx = goodch_idx;

        %Perform new GVDT
        nStd = 4;
        gvtdTimeTrace_medSNR = gvtd(PD_data.dod2_medSNR(:,goodch_idx)'); %input ch x t, output t x 1
        thresh_medSNR = find_gvtd_thresh(gvtdTimeTrace_medSNR, statType, nStd);
        gvdt_more_thresh_medSNR = find(gvtdTimeTrace_medSNR>thresh_medSNR);
        N_MA_GVDT_medSNR = size(gvdt_more_thresh_medSNR,1);
        
        figure()
        histogram(gvtdTimeTrace_medSNR)
        hold on
        xline([thresh_medSNR],LineWidth=4)
        xlabel('GVTD value / A.U')
        ylabel('Frequency / N. of timepoints')
        title("Histogram GVTD nSTD = "+nStd+" yields threshold of "+thresh_medSNR+" with "+N_MA_GVDT_medSNR+" M.As")
        
        tIncCh_medSNR = ones(size(PD_data.dod2_medSNR,1)  , size(PD_data.dod2_medSNR,2)  );
        tIncCh_medSNR(gvdt_more_thresh_medSNR,:) = 0;

        N_MA_HMR_medSNR = size(find(tIncCh_medSNR(:,goodch_idx)==0),1)/size(goodch_idx,1);
        figure()
        plot(PD_data.t,PD_data.dod2_medSNR(:,13))
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),PD_data.dod2_medSNR(find(tIncCh_medSNR(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("MedSNR Hmr MAD Using GVDT new nStd = "+nStd+" |N M.A = "+N_MA_HMR_medSNR+" ")

         figure()
         sgtitle("Med SNR GVTD new lambda "+nStd+" ")
         for i=1:12
            subplot(3,4,i)
            plot(PD_data.t, PD_data.dod2_medSNR(:,goodch_idx(i)))
            hold on
            plot(PD_data.t(find(tIncCh_medSNR(:,goodch_idx(i))==0)),PD_data.dod2_medSNR(find(tIncCh_medSNR(:,goodch_idx(i))==0),goodch_idx(i)),'m.','MarkerSize',6)
            xline([900],'b--')
            xline([5100],'b--')
            title("ch "+num2str(goodch_idx(i)))
            xlabel('T / s')
            ylabel('dod / AU')
          end
        
        g_medSNR_new = zeros(size(PD_data.dod2_medSNR,1),size(PD_data.dod2_medSNR,2));
        for i=2:size(PD_data.dod2_medSNR,1)
             g_medSNR_new(i,:)=sqrt( (PD_data.dod2_medSNR(i,:)-PD_data.dod2_medSNR(i-1,:)).^2 );
        end

        figure()
        subplot(2,1,1)
        imagesc(g_medSNR_new(:,goodch_idx)')
        colorbar, xlabel('Time / s'), ylabel('Measurement #'),title("GVDT new GoodCh medSNR lambda "+nStd+" ")
        subplot(2,1,2)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR(:,13)==0),1),'m.','MarkerSize',6)
         
        figure()
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR(:,13)==0),1),'m.','MarkerSize',6)

        % find trains/start and end of MA trains, without adding mask.
        mean(gvtdTimeTrace_medSNR(gvtdTimeTrace_medSNR < thresh_medSNR));
        
        tIncCh_medSNR;
        fs = 10;
        length_train = 100; %frames 
        IncCh_train = ones(1,size(tIncCh_medSNR,2));
        IncCh_train_loc = ones(size(tIncCh_medSNR,1),size(tIncCh_medSNR,2));
        for i=1:size(tIncCh_medSNR,2)
            if size(unique(tIncCh_medSNR(:,i)),1) > 1
                for j=1: size(tIncCh_medSNR,1)-length_train
                    if nnz(tIncCh_medSNR(j:j+(length_train-1),i)) < (length_train*0.85)
                        IncCh_train(1,i) = 0;
                        IncCh_train_loc(j,i) = 0;
                    end
                end
            end
        end

        IncCh_train = IncCh_train';
        MA_train_idx = IncCh_train;

        for i=1:size(MA_train_idx,1)
            if MA_train_idx(i,1) == 1;
                MA_train_idx(i,1) = i;
            end
        end

        MA_train_idx = nonzeros(MA_train_idx);
        figure()
        subplot(2,1,1)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR(:,13)==0),1),'m.','MarkerSize',6)
        subplot(2,1,2)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR(:,13)==0),1),'m.','MarkerSize',6)
        plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),gvtdTimeTrace_medSNR(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
        
        %find delta's in data. Consectuvive frames have delta =1. 
        %so a MA train starts/ends when there is a jump of 100
        IncCh_train_loc_delta = find(IncCh_train_loc(:,1)==0);
        %find(diff(find(IncCh_train_loc(:,1)==0))>200);
        %end of MA train
        end_MA_train(:,1) = find(diff(find(IncCh_train_loc(:,1)==0))>200);
        start_MA_train(:,1) = [1 ; end_MA_train(1:end-1,1)+1];
        %location of start/end of MA train in data (frame N)
        end_MA_train_locations_frame = IncCh_train_loc_delta(end_MA_train,1);
        start_MA_train_locations_frame = IncCh_train_loc_delta(start_MA_train,1);

        figure()
        subplot(2,1,1)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR(:,13)==0),1),'m.','MarkerSize',6)
        subplot(2,1,2)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR(:,13)==0),1),'m.','MarkerSize',6)
        plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),gvtdTimeTrace_medSNR(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
        plot(PD_data.t(start_MA_train_locations_frame),thresh_medSNR,'g*','MarkerSize',8)
        plot(PD_data.t(end_MA_train_locations_frame),thresh_medSNR,'b*','MarkerSize',8)

        
        figure()
        plot(PD_data.t,PD_data.dod2_medSNR(:,1))
        hold on
        plot(PD_data.t(start_MA_train_locations_frame),PD_data.dod2_medSNR(start_MA_train_locations_frame,1),'g*','MarkerSize',8)
        plot(PD_data.t(end_MA_train_locations_frame),PD_data.dod2_medSNR(end_MA_train_locations_frame,1),'b*','MarkerSize',8)
        %start_MA_train_locations_frame
       



        %tIncCh_medSNR
        %ADDING MASK and seeing MA trains
        %for i=1:size(tIncCh_medSNR,1)
        tMask = 0.2*fs; %mark data *- time around m.a as m.a
        tIncCh_medSNR_new = ones(size(PD_data.dod2_medSNR,1)  , size(PD_data.dod2_medSNR,2)  );
        for i=1:size(gvdt_more_thresh_medSNR,1)
            tIncCh_medSNR_new(gvdt_more_thresh_medSNR(i,1)-tMask:gvdt_more_thresh_medSNR(i,1)+tMask,:) = 0;
        end
        length_train = 100; %frames 
        IncCh_train = ones(1,size(tIncCh_medSNR,2));
        IncCh_train_loc = ones(size(tIncCh_medSNR,1),size(tIncCh_medSNR,2));
        for i=1:size(tIncCh_medSNR,2)
            if size(unique(tIncCh_medSNR(:,i)),1) > 1
                for j=1: size(tIncCh_medSNR,1)-length_train
                    if nnz(tIncCh_medSNR_new(j:j+(length_train-1),i)) < (length_train*0.85)
                        IncCh_train(1,i) = 0;
                        IncCh_train_loc(j,i) = 0;
                    end
                end
            end
        end

        IncCh_train = IncCh_train';
        MA_train_idx = IncCh_train;

        for i=1:size(MA_train_idx,1)
            if MA_train_idx(i,1) == 1;
                MA_train_idx(i,1) = i;
            end
        end

        MA_train_idx = nonzeros(MA_train_idx);
        figure()
        subplot(2,1,1)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR_new(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR_new(:,13)==0),1),'m.','MarkerSize',6)
        subplot(2,1,2)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR_new(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR_new(:,13)==0),1),'m.','MarkerSize',6)
        plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),gvtdTimeTrace_medSNR(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)

        %find delta's in data. Consectuvive frames have delta =1. 
        %so a MA train starts/ends when there is a jump of 100
        IncCh_train_loc_delta = find(IncCh_train_loc(:,1)==0);
        %find(diff(find(IncCh_train_loc(:,1)==0))>200);
        %end of MA train
        end_MA_train(:,1) = find(diff(find(IncCh_train_loc(:,1)==0))>200);
        end_MA_train = [end_MA_train ; size(IncCh_train_loc_delta,1)];
        start_MA_train(:,1) = [1 ; end_MA_train(1:end-1,1)+1];
        %location of start/end of MA train in data (frame N)
        end_MA_train_locations_frame = IncCh_train_loc_delta(end_MA_train,1);
        start_MA_train_locations_frame = IncCh_train_loc_delta(start_MA_train,1);
  
        figure()
        subplot(2,1,1)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR_new(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR_new(:,13)==0),1),'m.','MarkerSize',6)
        subplot(2,1,2)
        plot(PD_data.t,gvtdTimeTrace_medSNR)
        hold on
        plot(PD_data.t(find(tIncCh_medSNR_new(:,13)==0)),gvtdTimeTrace_medSNR(find(tIncCh_medSNR_new(:,13)==0),1),'m.','MarkerSize',6)
        plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),gvtdTimeTrace_medSNR(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
        plot(PD_data.t(start_MA_train_locations_frame),thresh_medSNR,'g*','MarkerSize',8)
        plot(PD_data.t(end_MA_train_locations_frame),thresh_medSNR,'b*','MarkerSize',8)

        figure()
        plot(PD_data.t,PD_data.dod2_medSNR(:,1))
        hold on
        plot(PD_data.t(start_MA_train_locations_frame),thresh_medSNR,'g*','MarkerSize',8)
        plot(PD_data.t(end_MA_train_locations_frame),thresh_medSNR,'b*','MarkerSize',8)

        %%% Determine mean before and after MA train
        for j=1:size(PD_data.dod2_medSNR,2)
            for i=1:size(end_MA_train_locations_frame,1)
                delta_dod_MA_train(i,j) = abs( 100* (1-(abs( mean(PD_data.dod2_medSNR(start_MA_train_locations_frame(i)-20:start_MA_train_locations_frame(i),j)) / mean(PD_data.dod2_medSNR(end_MA_train_locations_frame(i):end_MA_train_locations_frame(i)+20,j)) )))); 
            end
        end
        
        %%% If change is more than 10% - baseline shift
        %dod_MA_train_BLshift
        %find(delta_dod_MA_train(:,1)>10);
        
        %%% Assign a BL shift MA dect. We will run spline only on this.
        tIncCh_medSNR_new_BLshift = ones(size(PD_data.dod2_medSNR,1),size(PD_data.dod2_medSNR,2) );
        
        for i=1:size(PD_data.dod2_medSNR,2)
            
            start_BLshift = start_MA_train_locations_frame(find(delta_dod_MA_train(:,i)>10));
            end_BLshift = end_MA_train_locations_frame(find(delta_dod_MA_train(:,i)>10));

            for j=1:size(end_BLshift,1)
                tIncCh_medSNR_new_BLshift(start_BLshift(j):end_BLshift(j),i) = 0;
            end
        end
        
       testt = hmrMotionCorrectSpline(PD_data.dod2_medSNR,PD_data.t,PD_data.SD,tIncCh_medSNR_new_BLshift,p);
        
        figure()
        plot(PD_data.t,PD_data.dod2_medSNR(:,1))
        hold on
        plot(PD_data.t,testt(:,1))
        plot(PD_data.t(find(tIncCh_medSNR_new_BLshift(:,1)==0)),PD_data.dod2_medSNR(find(tIncCh_medSNR_new_BLshift(:,1)==0),1),'m.','MarkerSize',1)
        
        
        N_MA_HMR_medSNR = size(find(tIncCh_medSNR_new(:,goodch_idx)==0),1)/size(goodch_idx,1);
        figure()
        plot(PD_data.t,PD_data.dod2_medSNR(:,13))
        hold on
        plot(PD_data.t(find(tIncCh_medSNR_new(:,13)==0)),PD_data.dod2_medSNR(find(tIncCh_medSNR_new(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("MedSNR Hmr MAD Using GVDT new nStd = "+nStd+" |N M.A = "+N_MA_HMR_medSNR+" ")


        %figure()
        %histogram(mean(g_medSNR(:,goodch_idx)))

        tIncCh_lowSNR = ones(size(PD_data.dod2_lowSNR,1)  , size(PD_data.dod2_lowSNR,2)  );
        for i=1:size(gvdt_more_thresh_lowSNR,1)
            tIncCh_lowSNR(gvdt_more_thresh_lowSNR(i,1)-tMask:gvdt_more_thresh_lowSNR(i,1)+tMask,:) = 0;
        end
        %tIncCh_lowSNR(gvdt_more_thresh_lowSNR,:) = 0;
        %tIncCh_lowSNR(gvdt_more_thresh_lowSNR+1,:) = 0;
        
        tIncCh_highSNR = ones(size(PD_data.dod2_highSNR,1)  , size(PD_data.dod2_highSNR,2)  );
        for i=1:size(gvdt_more_thresh_highSNR,1)
            tIncCh_highSNR(gvdt_more_thresh_highSNR(i,1)-tMask:gvdt_more_thresh_highSNR(i,1)+tMask,:) = 0;
        end
        %tIncCh_highSNR(gvdt_more_thresh_highSNR,:) = 0;
        %tIncCh_highSNR(gvdt_more_thresh_highSNR+1,:) = 0;


        PD_data.dod2_medSNR = hmrMotionCorrectSpline(PD_data.dod2_medSNR,PD_data.t,PD_data.SD,tIncCh_medSNR_new,p);
        
        figure()
        plot(PD_data.t,PD_data.dod2_og_medSNR(:,1))
        hold on
        plot(PD_data.t,PD_data.dod2_medSNR(:,1))
        title("Spline with GVTD-new nSTD = "+nStd+" ")


        PD_data.dod2_lowSNR = hmrMotionCorrectSpline(PD_data.dod2_lowSNR,PD_data.t,PD_data.SD,tIncCh_lowSNR,p);
        PD_data.dod2_highSNR = hmrMotionCorrectSpline(PD_data.dod2_highSNR,PD_data.t,PD_data.SD,tIncCh_highSNR,p);
        'Spline'
        
        if methodN == 6 
            iqr = 1.1;
            PD_data.dod2_medSNR = hmrMotionCorrectWavelet(PD_data.dod2_medSNR,PD_data.SD,iqr); %run wavelet 
            PD_data.dod2_lowSNR = hmrMotionCorrectWavelet(PD_data.dod2_lowSNR,PD_data.SD,iqr); %run wavelet 
            PD_data.dod2_highSNR = hmrMotionCorrectWavelet(PD_data.dod2_highSNR,PD_data.SD,iqr); %run wavelet 
            'Wavlet'
        end

    end
end

% Bandpass Filter
lowerCutOff = 0;
higherCutOff = 0.001;
fs = 10; %sampling rate (Hz)
PD_data.dod3_lowSNR = hmrBandpassFilt(PD_data.dod2_lowSNR,fs,lowerCutOff,higherCutOff);
PD_data.dod3_medSNR = hmrBandpassFilt(PD_data.dod2_medSNR,fs,lowerCutOff,higherCutOff);
PD_data.dod3_highSNR = hmrBandpassFilt(PD_data.dod2_highSNR,fs,lowerCutOff,higherCutOff);
'BP'

%PD_data.dod3 = PD_data.dod2;
end