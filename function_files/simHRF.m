function PD_data = simHRF(PD_data)
%check data_RS if it needs to be transposed by output of yout if it's frames x 1 or frames x frames
Fc = 10; %1 sample = 1 minute
%data_RS = mean(HbO_c(cortex_nodes,:));
duration_hrf = 30*60; 
nHRF = 1;
distance = 9000; %Shits the HRF back/forward in time %9000 is at 15 minutes. (15min * 60s * 10Hz)
nChrom = 1;
block =2;
blockDuration = 1000; %length of block being conv. with
data_RS = PD_data.HbO_1(:,1);
[yout,vett_hrf,u,t_hrf,hrf_avg]=addHRF_infant_version_bothChrom(PD_data.t,Fc,data_RS,duration_hrf,nHRF,distance,nChrom,block,blockDuration);

% figure()
% plot(PD_data.t, vett_hrf)
% xlabel('Time / s');
% ylabel('Sim HRF / M')

    for i=1:size(PD_data.dod1,2)/2
        %data_RS = HbO_1(:,ch_idx(i));
        data_RS = PD_data.HbO_1(:,i);

        [yout,vett_hrf,u,t_hrf,hrf_avg]=addHRF_infant_version_bothChrom(PD_data.t,Fc,data_RS,duration_hrf,nHRF,distance,nChrom,block,blockDuration);
        PD_data.HbO_sim1(:,i) = yout;
    
        %data_RS = Hb_1(:,ch_idx(i));
        data_RS = PD_data.Hb_1(:,i);

        [yout,vett_hrf,u,t_hrf,hrf_avg]=addHRF_infant_version_bothChrom(PD_data.t,Fc,data_RS,duration_hrf,nHRF,distance,nChrom,block,blockDuration);
        PD_data.Hb_sim1(:,i) = yout;

        %data_RS = HbT_1(:,ch_idx(i));
        data_RS = PD_data.HbT_1(:,i);

        [yout,vett_hrf,u,t_hrf,hrf_avg]=addHRF_infant_version_bothChrom(PD_data.t,Fc,data_RS,duration_hrf,nHRF,distance,nChrom,block,blockDuration);
        PD_data.HbT_sim1(:,i) = yout;
end

dc_1_sim(:,1,:) = PD_data.HbO_sim1;
dc_1_sim(:,2,:) = PD_data.Hb_sim1;
dc_1_sim(:,3,:) = PD_data.HbT_sim1;

PD_data.dc_1_sim = dc_1_sim;

PD_data.vett_hrf = vett_hrf;
PD_data.u = u;
PD.data.t_hrf = t_hrf;
PD.data.hrf_avg = hrf_avg;


end