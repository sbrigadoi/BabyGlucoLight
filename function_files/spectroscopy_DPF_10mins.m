function [PD_data] = spectroscopy_DPF_10mins(PD_data)

% Step 5. Spectoscropy
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
dc_1_noc = hmrOD2Conc( PD_data.dod1_noc, PD_data.SD, PD_data.ppf );

[HbO_x,Hb_x,HbT_x] = dc2HbX(dc_1);
PD_data.HbO_1 = HbO_x;
PD_data.Hb_1 = Hb_x;
PD_data.HbT_1 = HbT_x;

[HbO_x,Hb_x,HbT_x] = dc2HbX(dc_1_noc);
PD_data.HbO_1_noc = HbO_x;
PD_data.Hb_1_noc = Hb_x;
PD_data.HbT_1_noc = HbT_x;

figure()
plot(PD_data.HbO_1_noc(:,1))
hold on
plot(PD_data.HbO_1(:,1))
legend('Only BP','Int,splineSG,BP')
ylabel('d HbO (M)')
xlabel('Time / s')
xline(900,'g--');

figure()
plot(PD_data.HbO_1(:,PD_data.goodch_idx(1:end/2)))
ylabel('d HbO (M)')
xlabel('Time / s')
xline(900,'g--');
title("Good Ch HbO from dod1")

figure()
plot(PD_data.Hb_1(:,PD_data.goodch_idx(1:end/2)))
ylabel('d HbO (M)')
xlabel('Time / s')
xline(900,'g--');
title("Good Ch Hb from dod1")

figure()
plot(PD_data.HbT_1(:,PD_data.goodch_idx(1:end/2)))
ylabel('d HbO (M)')
xlabel('Time / s')
xline(900,'g--');
title("Good Ch HbT from dod1")

end 