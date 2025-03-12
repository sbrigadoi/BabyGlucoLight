function [goodch_idx PD_data] = find_good_ch_15mins(PD_data,dRange,SNRrange,figON)
      
    fs = 10; %sampling rate (Hz)

        % if figON == 0
        %     set(gcf,'Visible','off');              
        %     set(0,'DefaultFigureVisible','off');
        % 
        % elseif figON == 1
        %     set(gcf,'Visible','on');              
        %     set(0,'DefaultFigureVisible','on');
        % end

    gap_between_SNRcheck_min = 2;
    length_SNR_window_min = 2;
    length_data_min = size(PD_data.t,1)/(fs*60);
    N_15mins = floor((length_data_min-2)/gap_between_SNRcheck_min );
    %
    for i=0:N_15mins
        start_SNR(i+1,1) = ((fs*1*60*i*gap_between_SNRcheck_min)+1);
        end_SNR(i+1,1) = start_SNR(i+1,1) + (fs*60*length_SNR_window_min);
    end


    remCh_allT=zeros(size(PD_data.d,2),N_15mins);
    for i=1:N_15mins
        d = PD_data.d(start_SNR(i):end_SNR(i),:);
        remCh_allT(:,i) = removeNoisyChannelsPDdata(d,dRange,SNRrange);
    end

    %remCh_baseline_event = [remCh_baseline remCh_event];
    for i=1:size(remCh_allT,1)
        if sum(remCh_allT(i,:))==N_15mins
            remCh(i,1) =1;
        else
            remCh(i,1) =0;
        end
    end
    PD_data.SD.MeasListAct = remCh;

    goodch_idx = zeros(size(remCh,1),1);
    goodch_idx(find(PD_data.SD.MeasListAct==1),1) = find(PD_data.SD.MeasListAct==1);
    
    %goodch_idx = nonzeros(goodch_idx);
    
    %get rid of channels with jump > 0.15
    for i=1:size(remCh,1)
        if abs( max(PD_data.d(:,i))-min(PD_data.d(:,i)) ) > 0.15 %bad ch
            goodch_idx(i,1)=0;
            remCh(i,1) =0;
        end
    end
    PD_data.SD.MeasListAct = remCh;
    goodch_idx = nonzeros(goodch_idx);

    %get rid of channels with detector saturation high freq. oscillation
    %FFT
    %Fs = boxy_hdr.sample_rate; % Sampling frequency                    
    T = 1/fs;  % Sampling period
    L=size(PD_data.d,1);
    f = fs*(0:(L/2))/L;
    % Length of signal
    t = (0:L-1)*T;

    for k= 1:size(PD_data.d,2);
        Xiac(k,:) = PD_data.d(:,k);
        Yiac(k,:) = fft(Xiac(k,:));

        P2iac(k,:) = abs(Yiac(k,:)/L);
        P1iac(k,:) = P2iac(k,1:L/2+1);
        P1iac(k,2:end-1) = 2*P1iac(k,2:end-1);     
    end

    %figure()
    %plot(PD_data.t,PD_data.d(:,goodch_idx))
    %xlabel("T / s")

    %figure()
    %semilogx(f,P1iac(goodch_idx,:))
    %xlabel("f / Hz")

    %main freq. peak of detector oscilliations is 0.0516 Hz
    % col 32 of f is 0.0516 Hz
    % 0.04 Hz in f is col 25 and 0.06 Hz in f is col 37

    %mean_p_f0_0516 = mean(P1iac(goodch_idx,25:37)');
    %std_mean_p_f0_0516 = std(mean_p_f0_0516);
    %avg_mean_p_f0_0516 = mean(mean(P1iac(goodch_idx,25:37)'));
    %outl_p_f_0_0516 = find(mean_p_f0_0516> (avg_mean_p_f0_0516 + (2*std_mean_p_f0_0516)) )

    %max pvalue
    max_p_f0_0516 = max(P1iac(goodch_idx,25:41)'); %max values of p between 0.04 and 0.0667 Hz
    avg_max_p_f0_0516 = mean(max_p_f0_0516); %avg value of the max values
    std_max_p_f0_0516 = std(max_p_f0_0516); %std of the max values
    outl_p_f_0_0516 = find(max_p_f0_0516> (avg_max_p_f0_0516 + (1.5*std_max_p_f0_0516)) ); %outlier channels are 1.5 STD above the average

    max_p_f0_0067 = max(P1iac(goodch_idx,5:7)'); %max values of p between 0.04 and 0.0667 Hz
    avg_max_p_f0_0067 = mean(max_p_f0_0067); %avg value of the max values
    std_max_p_f0_0067 = std(max_p_f0_0067); %std of the max values
    outl_p_f_0_0067 = find(max_p_f0_0067> (avg_max_p_f0_0067 + (1.5*std_max_p_f0_0067)) ); %outlier channels are 1.5 STD above the average

    max_p_f0_0150 = max(P1iac(goodch_idx,9:12)'); %max values of p between 0.04 and 0.0667 Hz
    avg_max_p_f0_0150 = mean(max_p_f0_0150); %avg value of the max values
    std_max_p_f0_0150 = std(max_p_f0_0150); %std of the max values
    outl_p_f_0_0150 = find(max_p_f0_0150> (avg_max_p_f0_0150 + (1.5*std_max_p_f0_0150)) ); %outlier channels are 1.5 STD above the average

    outl_p_f = unique([outl_p_f_0_0516  outl_p_f_0_0067  outl_p_f_0_0150]);

    %figure()
    %semilogx(f,P1iac(goodch_idx(outl_p_f_0_0516),:))
    if size(outl_p_f,2)>1
        %set these outliers to 0 index in goodch idx
        remCh(goodch_idx(outl_p_f,1),1)=0;
        goodch_idx(outl_p_f,1) = 0;
        goodch_idx = nonzeros(goodch_idx);
    end

    % if size(outl_p_f_0_0067,2)>1
    %     %set these outliers to 0 index in goodch idx
    %     remCh(goodch_idx(outl_p_f_0_0067,1),1)=0;
    %     goodch_idx(outl_p_f_0_0067,1) = 0;
    %     goodch_idx = nonzeros(goodch_idx);
    % end
    % 
    % if size(outl_p_f_0_0150,2)>1
    %     %set these outliers to 0 index in goodch idx
    %     remCh(goodch_idx(outl_p_f_0_0150,1),1)=0;
    %     goodch_idx(outl_p_f_0_0150,1) = 0;
    %     goodch_idx = nonzeros(goodch_idx);
    % end

    %figure()
    %semilogx(f,(P1iac(goodch_idx,:)))

    %figure()
    %semilogx(f,P1iac(goodch_idx(outl_p_f_0_0516),:))



    %need to check that channels are on both WL's and update 
    %goodch idx
    %rem ch and PD data DS measlistAct.
    %size(goodch_idx,1)
    
    goodch_idx_check = goodch_idx;
    for i=1:size(goodch_idx,1)
        if (find(goodch_idx_check == goodch_idx_check(i)+64)) | (find(goodch_idx_check == goodch_idx_check(i)-64))
            %goodch_idx(i)
            goodch_idx(i) = goodch_idx(i);
            %remCh(i)=1;
        else
            goodch_idx(i) = 0;
        end
    end
    
    goodch_idx = nonzeros(goodch_idx);
    goodch_idx = [goodch_idx; goodch_idx+64];

    if size(goodch_idx,1) <2 %less than 2good chs, i.e 0 or 1
        goodch_idx = [1;65];
    end

    %remCh
    %size(goodch_idx)
    remCh(goodch_idx)=1;
    PD_data.SD.MeasListAct(goodch_idx) = 1;
    PD_data.goodch_idx = goodch_idx;
    
    if figON == 1
        figure()
        plot(PD_data.t,PD_data.d(:,PD_data.goodch_idx))
        xlabel('Time / S')
        ylabel('D / A.U')
        title("Raw D - goodch idx - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+" ")
    end
    

    %goodch_idx = zeros(size(remCh,1),1);
    %convert remCh to channel index
    %for i=1:size(remCh,1)
    %    if remCh(i,1)==1
    %        goodch_idx(i,1)=i;
    %    end
    %end
    %goodch_idx = nonzeros(goodch_idx);

    % set(gcf,'Visible','on');              
    % set(0,'DefaultFigureVisible','on');


end