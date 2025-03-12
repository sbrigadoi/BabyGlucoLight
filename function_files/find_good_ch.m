function [goodch_idx PD_data] = find_good_ch(PD_data,dRange,SNRrange)
        
    fs = 10; %sampling rate (Hz)

    gap_between_SNRcheck_min = 15;
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
    %convert remCh to channel index
    for i=1:size(remCh,1)
        if remCh(i,1)==1
            goodch_idx(i,1)=i;
        end
    end
    goodch_idx = nonzeros(goodch_idx);

end