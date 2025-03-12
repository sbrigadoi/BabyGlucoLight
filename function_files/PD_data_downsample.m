function [PD_data] = PD_data_downsample(PD_data,ds_rate)

PD_data.dod1 = downsample(PD_data.dod1,ds_rate);
PD_data.dod1_noc = downsample(PD_data.dod1_noc,ds_rate);
PD_data.t = downsample(PD_data.t,ds_rate);
PD_data.s = downsample(PD_data.s,ds_rate);
PD_data.d = downsample(PD_data.d,ds_rate);
PD_data.aux = downsample(PD_data.aux,ds_rate);

% figure()
% plot(PD_data.t,PD_data.dod1_noc(:,1))
% hold on
% plot(PD_data.t,PD_data.dod1(:,1))

figure()
plot(PD_data.t,PD_data.dod1(:,PD_data.goodch_idx))
xlabel('Time / S')
ylabel('\Delta OD')
title("Dod 1 (Mo Corr) - DS by factor "+ds_rate+" PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+"")

end