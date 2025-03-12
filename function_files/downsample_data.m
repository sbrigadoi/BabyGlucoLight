function PD_data = downsample_data(PD_data,ds) 

    for i = 1 : size(PD_data.dod,2)
        Y_logI_1Hz(:,i) = mean(reshape(PD_data.dod(1:size(PD_data.dod,1) - mod(size(PD_data.dod,1),ds),i),ds,[])); %160 was 40
    end

    t_ds = 1:1:size(Y_logI_1Hz,1)';

    dod = Y_logI_1Hz;

    PD_data.dod = dod;
    PD_data.t_ds = t_ds;


end