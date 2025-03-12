% T = readtable('CGM_NICU2020.xlsx');
% 
% a = T(1,6)
% %%
% T = readmatrix('CGM_NICU2020.xlsx');
% %%
% 
% A = readmatrix('PD1.xlsx');
% 
% A = readmatrix('PD1_1.xlsx');
% %%
% B = readtable('PD1.xlsx');
% 
% B = readcell('PD1.xlsx');
% 
% 
% B = readcell('PD1_1.xlsx');
% B = readtable('PD1_1.xlsx','VariableNamingRule','preserve');
% 
% %C = cell2mat( B(2:end,3) )
% 
% C = table2array(B(1:end,3))
% D = cell2mat(C)
% str2num(D)
% 
% 
% D =str2double(C);
%% Get Glucose value for format '188,00' i.e string with comma seperator
% copy and paste just glucose values into new excel sheet
%B = readcell('PD1_1.xlsx');
clear all;
close all;
clc;

%% work
name_spreadsheet = 'PD1.xlsx';
type_spreadsheet = 1; %For glucose saved as a string
PD1 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% (if you get the readcell error, only set path to the external hard drive ONLY BABYGLUCO )
% work
name_spreadsheet = 'PD2.xlsx';
type_spreadsheet = 1; %For glucose saved as a string
PD2 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
% work
name_spreadsheet = 'PD23.xlsx';
type_spreadsheet = 1; %For glucose saved as a string
PD23 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD3.xlsx';
type_spreadsheet = 2; %For glucose saved as a number
PD3 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
name_spreadsheet = 'PD4.xlsx';
PD4 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD5.xlsx';
PD5 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD6.xlsx';
PD6 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD7.xlsx';
PD7 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD8.xlsx';
PD8 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD9.xlsx';
PD9 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work (if you get the readcell error, only set path to the external hard drive)
type_spreadsheet = 2; %For glucose saved as a number
name_spreadsheet = 'PD10.xlsx';
PD10 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
%work
name_spreadsheet = 'PD12.xlsx';
PD12 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD13.xlsx';
PD13 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
%  work
name_spreadsheet = 'PD14.xlsx';
PD14 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD15.xlsx';
PD15 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%  work
name_spreadsheet = 'PD16.xlsx';
PD16 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD17.xlsx';
PD17 = readglucose_string_number(name_spreadsheet,type_spreadsheet)

% work
name_spreadsheet = 'PD19.xlsx';
PD19 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%  work
name_spreadsheet = 'PD20.xlsx';
PD20 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%  work
name_spreadsheet = 'PD24.xlsx';
PD24 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD25.xlsx';
PD25 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD26.xlsx';
PD26 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD28.xlsx';
PD28 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%work
name_spreadsheet = 'PD33.xlsx';
PD33 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD34.xlsx';
PD34 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%  work
name_spreadsheet = 'PD35.xlsx';
PD35 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work - date being read strange??
%
name_spreadsheet = 'PD41.xlsx';
PD41 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
%  work
name_spreadsheet = 'PD42.xlsx';
PD42 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD45.xlsx';
PD45 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
name_spreadsheet = 'PD48.xlsx';
PD48 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%work
name_spreadsheet = 'PD49.xlsx';
PD49 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD50.xlsx';
PD50 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%  work
name_spreadsheet = 'PD53.xlsx';
PD53 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
% work
name_spreadsheet = 'PD55.xlsx';
PD55 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
% making new ones work
name_spreadsheet = 'PD11.xlsx'; %type 2
PD11 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD36.xlsx'; %type 2
%PD36 = readglucose_number(name_spreadsheet)
PD36 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
name_spreadsheet = 'PD37.xlsx'; %type 2
%PD37 = readglucose_number(name_spreadsheet)
PD37 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %date needed to be converted to datetime from string

name_spreadsheet = 'PD39.xlsx'; %type 2
PD39 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %date needed to be converted to datetime from str
%
name_spreadsheet = 'PD40.xlsx'; %type 2
PD40 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %date needed to be converted to datetime from str
%
name_spreadsheet = 'PD43.xlsx'; %type 2
PD43 = readglucose_string_number(name_spreadsheet,type_spreadsheet)  %date needed to be converted to datetime from str
%
name_spreadsheet = 'PD47.xlsx'; %type 2
PD47 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %dates were in mixed US and UK format

name_spreadsheet = 'PD51.xlsx'; %type 2
PD51 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD52.xlsx'; %type 2
PD52 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD54.xlsx'; %type 2
PD54 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD22.xlsx'; %type 2
PD22 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD44.xlsx'; %type 2
PD44 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD56.xlsx'; %type 2
PD56 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD57.xlsx'; %type 2
PD57 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD58.xlsx'; %type 2
PD58 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD59.xlsx'; %type 2
PD59 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

name_spreadsheet = 'PD60.xlsx'; %type 2
PD60 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

%had to make a manual solution
%name_spreadsheet = 'PD38.xlsx'; %type 2
%PD38 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double
load('PD38.mat');
%% make ones remaining work!
type_spreadsheet = 2; %For glucose saved as a string
name_spreadsheet = 'PD60.xlsx'; %type 2
PD60 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double
%%

type_spreadsheet = 2; %For glucose saved as a string
name_spreadsheet = 'PD48.xlsx'; %type 2
PD48 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %last time entry was in datetime format, should of been double

%%

PD_0 = [0.1 0.1 0.1 0.1];

PD_hyperSm_hypoSm = [PD1.N_events.S_hyper PD1.N_events.m_hyper PD1.N_events.S_hypo PD1.N_events.m_hypo;
PD2.N_events.S_hyper PD2.N_events.m_hyper PD2.N_events.S_hypo PD2.N_events.m_hypo;
PD3.N_events.S_hyper PD3.N_events.m_hyper PD3.N_events.S_hypo PD3.N_events.m_hypo;
PD4.N_events.S_hyper PD4.N_events.m_hyper PD4.N_events.S_hypo PD4.N_events.m_hypo;
PD5.N_events.S_hyper PD5.N_events.m_hyper PD5.N_events.S_hypo PD5.N_events.m_hypo;
PD6.N_events.S_hyper PD6.N_events.m_hyper PD6.N_events.S_hypo PD6.N_events.m_hypo;
PD7.N_events.S_hyper PD7.N_events.m_hyper PD7.N_events.S_hypo PD7.N_events.m_hypo;
PD8.N_events.S_hyper PD8.N_events.m_hyper PD8.N_events.S_hypo PD8.N_events.m_hypo;
PD_0;
PD10.N_events.S_hyper PD10.N_events.m_hyper PD10.N_events.S_hypo PD10.N_events.m_hypo;
PD_0;
PD_0;
PD13.N_events.S_hyper PD13.N_events.m_hyper PD13.N_events.S_hypo PD13.N_events.m_hypo;
PD14.N_events.S_hyper PD14.N_events.m_hyper PD14.N_events.S_hypo PD14.N_events.m_hypo;
PD15.N_events.S_hyper PD15.N_events.m_hyper PD15.N_events.S_hypo PD15.N_events.m_hypo;
PD16.N_events.S_hyper PD16.N_events.m_hyper PD16.N_events.S_hypo PD16.N_events.m_hypo;
PD17.N_events.S_hyper PD17.N_events.m_hyper PD17.N_events.S_hypo PD17.N_events.m_hypo;
PD_0;
PD19.N_events.S_hyper PD19.N_events.m_hyper PD19.N_events.S_hypo PD19.N_events.m_hypo;
PD20.N_events.S_hyper PD20.N_events.m_hyper PD20.N_events.S_hypo PD20.N_events.m_hypo;
PD_0;
PD_0;
PD23.N_events.S_hyper PD23.N_events.m_hyper PD23.N_events.S_hypo PD23.N_events.m_hypo;
PD24.N_events.S_hyper PD24.N_events.m_hyper PD24.N_events.S_hypo PD24.N_events.m_hypo;
PD25.N_events.S_hyper PD25.N_events.m_hyper PD25.N_events.S_hypo PD25.N_events.m_hypo;
PD26.N_events.S_hyper PD26.N_events.m_hyper PD26.N_events.S_hypo PD26.N_events.m_hypo;
PD_0;
PD_0;
PD_0;
PD_0;
PD_0;
PD_0;
PD33.N_events.S_hyper PD33.N_events.m_hyper PD33.N_events.S_hypo PD33.N_events.m_hypo;
PD34.N_events.S_hyper PD34.N_events.m_hyper PD34.N_events.S_hypo PD34.N_events.m_hypo;
PD35.N_events.S_hyper PD35.N_events.m_hyper PD35.N_events.S_hypo PD35.N_events.m_hypo;
PD_0;
PD_0;
PD_0;
PD_0;
PD_0;
PD41.N_events.S_hyper PD41.N_events.m_hyper PD41.N_events.S_hypo PD41.N_events.m_hypo;
PD42.N_events.S_hyper PD42.N_events.m_hyper PD42.N_events.S_hypo PD42.N_events.m_hypo;
PD_0;
PD_0;
PD45.N_events.S_hyper PD45.N_events.m_hyper PD45.N_events.S_hypo PD45.N_events.m_hypo;
PD_0;
PD_0;
PD48.N_events.S_hyper PD48.N_events.m_hyper PD48.N_events.S_hypo PD48.N_events.m_hypo;
PD49.N_events.S_hyper PD49.N_events.m_hyper PD49.N_events.S_hypo PD49.N_events.m_hypo;
PD50.N_events.S_hyper PD50.N_events.m_hyper PD50.N_events.S_hypo PD50.N_events.m_hypo;
PD_0;
PD_0;
PD53.N_events.S_hyper PD53.N_events.m_hyper PD53.N_events.S_hypo PD53.N_events.m_hypo;
PD_0;
PD55.N_events.S_hyper PD55.N_events.m_hyper PD55.N_events.S_hypo PD55.N_events.m_hypo;]

PD_hyperSm_hypoSm(:,3:4) = -1* PD_hyperSm_hypoSm(:,3:4);

figure()
bar(PD_hyperSm_hypoSm,'stacked')
legend('S Hyper','m Hyper','S Hypo','m Hypo')
xlabel('PD Number');
ylabel('N events');

%%

PD_compare(:,2) = PD_hyperSm_hypoSm(:,4) ;
PD_compare(:,3) = PD_hyperSm_hypoSm(:,3) ;
PD_compare(:,4) = PD_hyperSm_hypoSm(:,2) ;
PD_compare(:,5) = PD_hyperSm_hypoSm(:,1) ;
PD_compare(:,1) = sum(abs(PD_hyperSm_hypoSm),2);

%% run to get bar chart of all classified events 18 04 24
%updated for all PD from 01 to 60, without 8 with missing CGM 16 10 24

PD_0 = [0.1 0.1 0.1 0.1];

PD_0 = [0 0 0 0];

PD_hyperSm_hypoSm = [PD1.N_events.S_hyper PD1.N_events.m_hyper PD1.N_events.S_hypo PD1.N_events.m_hypo;
PD2.N_events.S_hyper PD2.N_events.m_hyper PD2.N_events.S_hypo PD2.N_events.m_hypo;
PD3.N_events.S_hyper PD3.N_events.m_hyper PD3.N_events.S_hypo PD3.N_events.m_hypo;
PD4.N_events.S_hyper PD4.N_events.m_hyper PD4.N_events.S_hypo PD4.N_events.m_hypo;
PD5.N_events.S_hyper PD5.N_events.m_hyper PD5.N_events.S_hypo PD5.N_events.m_hypo;
PD6.N_events.S_hyper PD6.N_events.m_hyper PD6.N_events.S_hypo PD6.N_events.m_hypo;
PD7.N_events.S_hyper PD7.N_events.m_hyper PD7.N_events.S_hypo PD7.N_events.m_hypo;
PD8.N_events.S_hyper PD8.N_events.m_hyper PD8.N_events.S_hypo PD8.N_events.m_hypo;
PD9.N_events.S_hyper PD9.N_events.m_hyper PD9.N_events.S_hypo PD9.N_events.m_hypo;
PD10.N_events.S_hyper PD10.N_events.m_hyper PD10.N_events.S_hypo PD10.N_events.m_hypo;
PD11.N_events.S_hyper PD11.N_events.m_hyper PD11.N_events.S_hypo PD11.N_events.m_hypo;
PD12.N_events.S_hyper PD12.N_events.m_hyper PD12.N_events.S_hypo PD12.N_events.m_hypo;
PD13.N_events.S_hyper PD13.N_events.m_hyper PD13.N_events.S_hypo PD13.N_events.m_hypo;
PD14.N_events.S_hyper PD14.N_events.m_hyper PD14.N_events.S_hypo PD14.N_events.m_hypo;
PD15.N_events.S_hyper PD15.N_events.m_hyper PD15.N_events.S_hypo PD15.N_events.m_hypo;
PD16.N_events.S_hyper PD16.N_events.m_hyper PD16.N_events.S_hypo PD16.N_events.m_hypo;
PD17.N_events.S_hyper PD17.N_events.m_hyper PD17.N_events.S_hypo PD17.N_events.m_hypo;
PD_0;
PD19.N_events.S_hyper PD19.N_events.m_hyper PD19.N_events.S_hypo PD19.N_events.m_hypo;
PD20.N_events.S_hyper PD20.N_events.m_hyper PD20.N_events.S_hypo PD20.N_events.m_hypo;
PD_0;
PD22.N_events.S_hyper PD22.N_events.m_hyper PD22.N_events.S_hypo PD22.N_events.m_hypo;
PD23.N_events.S_hyper PD23.N_events.m_hyper PD23.N_events.S_hypo PD23.N_events.m_hypo;
PD24.N_events.S_hyper PD24.N_events.m_hyper PD24.N_events.S_hypo PD24.N_events.m_hypo;
PD25.N_events.S_hyper PD25.N_events.m_hyper PD25.N_events.S_hypo PD25.N_events.m_hypo;
PD26.N_events.S_hyper PD26.N_events.m_hyper PD26.N_events.S_hypo PD26.N_events.m_hypo;
PD_0;
PD28.N_events.S_hyper PD28.N_events.m_hyper PD28.N_events.S_hypo PD28.N_events.m_hypo;
PD_0;
PD_0;
PD_0;
PD_0;
PD33.N_events.S_hyper PD33.N_events.m_hyper PD33.N_events.S_hypo PD33.N_events.m_hypo;
PD34.N_events.S_hyper PD34.N_events.m_hyper PD34.N_events.S_hypo PD34.N_events.m_hypo;
PD35.N_events.S_hyper PD35.N_events.m_hyper PD35.N_events.S_hypo PD35.N_events.m_hypo;
PD36.N_events.S_hyper PD36.N_events.m_hyper PD36.N_events.S_hypo PD36.N_events.m_hypo;
PD37.N_events.S_hyper PD37.N_events.m_hyper PD37.N_events.S_hypo PD37.N_events.m_hypo;
PD38.N_events.S_hyper PD38.N_events.m_hyper PD38.N_events.S_hypo PD38.N_events.m_hypo;
PD39.N_events.S_hyper PD39.N_events.m_hyper PD39.N_events.S_hypo PD39.N_events.m_hypo;
PD40.N_events.S_hyper PD40.N_events.m_hyper PD40.N_events.S_hypo PD40.N_events.m_hypo;
PD41.N_events.S_hyper PD41.N_events.m_hyper PD41.N_events.S_hypo PD41.N_events.m_hypo;
PD42.N_events.S_hyper PD42.N_events.m_hyper PD42.N_events.S_hypo PD42.N_events.m_hypo;
PD43.N_events.S_hyper PD43.N_events.m_hyper PD43.N_events.S_hypo PD43.N_events.m_hypo;
PD44.N_events.S_hyper PD44.N_events.m_hyper PD44.N_events.S_hypo PD44.N_events.m_hypo;
PD45.N_events.S_hyper PD45.N_events.m_hyper PD45.N_events.S_hypo PD45.N_events.m_hypo;
PD_0;
PD47.N_events.S_hyper PD47.N_events.m_hyper PD47.N_events.S_hypo PD47.N_events.m_hypo;
PD48.N_events.S_hyper PD48.N_events.m_hyper PD48.N_events.S_hypo PD48.N_events.m_hypo;
PD49.N_events.S_hyper PD49.N_events.m_hyper PD49.N_events.S_hypo PD49.N_events.m_hypo;
PD50.N_events.S_hyper PD50.N_events.m_hyper PD50.N_events.S_hypo PD50.N_events.m_hypo;
PD51.N_events.S_hyper PD51.N_events.m_hyper PD51.N_events.S_hypo PD51.N_events.m_hypo;
PD52.N_events.S_hyper PD52.N_events.m_hyper PD52.N_events.S_hypo PD52.N_events.m_hypo;
PD53.N_events.S_hyper PD53.N_events.m_hyper PD53.N_events.S_hypo PD53.N_events.m_hypo;
PD54.N_events.S_hyper PD54.N_events.m_hyper PD54.N_events.S_hypo PD54.N_events.m_hypo;
PD55.N_events.S_hyper PD55.N_events.m_hyper PD55.N_events.S_hypo PD55.N_events.m_hypo;
PD56.N_events.S_hyper PD56.N_events.m_hyper PD56.N_events.S_hypo PD56.N_events.m_hypo;
PD57.N_events.S_hyper PD57.N_events.m_hyper PD57.N_events.S_hypo PD57.N_events.m_hypo;
PD58.N_events.S_hyper PD58.N_events.m_hyper PD58.N_events.S_hypo PD58.N_events.m_hypo;
PD59.N_events.S_hyper PD59.N_events.m_hyper PD59.N_events.S_hypo PD59.N_events.m_hypo;
PD60.N_events.S_hyper PD60.N_events.m_hyper PD60.N_events.S_hypo PD60.N_events.m_hypo;];

PD_hyperSm_hypoSm(:,3:4) = -1* PD_hyperSm_hypoSm(:,3:4);

PD_N = [1:1:60];

figure()
b = bar(PD_N,PD_hyperSm_hypoSm,0.9,'stacked','FaceColor','flat')
%legend('S Hyper','m Hyper','S Hypo','m Hypo')
xlabel('PD Number');
ylabel('N events');
grid on
b(1).CData = [255/255 86/255 7/255];
b(2).CData = [255/255 147/255 6/255];
b(3).CData = [137/255 5/255 209/255];
b(4).CData = [35/255 33/255 172/255];
 ax = gca;
    ax.FontSize = 26; %20
%b = bar(y,'FaceColor','flat');
% for k = 1:size(PD_hyperSm_hypoSm,2)
%     b(k).CData = k;
% end

%% %updated for all PD from 01 to 60, without 8 with missing CGM 17 10 24

PD_N = [1:1:17 19 20 22:1:26 28 33:1:45 47:1:60  ]';

PD_0 = [0.1 0.1 0.1 0.1];
PD_hyperSm_hypoSm = [PD1.N_events.S_hyper PD1.N_events.m_hyper PD1.N_events.S_hypo PD1.N_events.m_hypo;
PD2.N_events.S_hyper PD2.N_events.m_hyper PD2.N_events.S_hypo PD2.N_events.m_hypo;
PD3.N_events.S_hyper PD3.N_events.m_hyper PD3.N_events.S_hypo PD3.N_events.m_hypo;
PD4.N_events.S_hyper PD4.N_events.m_hyper PD4.N_events.S_hypo PD4.N_events.m_hypo;
PD5.N_events.S_hyper PD5.N_events.m_hyper PD5.N_events.S_hypo PD5.N_events.m_hypo;
PD6.N_events.S_hyper PD6.N_events.m_hyper PD6.N_events.S_hypo PD6.N_events.m_hypo;
PD7.N_events.S_hyper PD7.N_events.m_hyper PD7.N_events.S_hypo PD7.N_events.m_hypo;
PD8.N_events.S_hyper PD8.N_events.m_hyper PD8.N_events.S_hypo PD8.N_events.m_hypo;
PD9.N_events.S_hyper PD9.N_events.m_hyper PD9.N_events.S_hypo PD9.N_events.m_hypo;
PD10.N_events.S_hyper PD10.N_events.m_hyper PD10.N_events.S_hypo PD10.N_events.m_hypo;
PD11.N_events.S_hyper PD11.N_events.m_hyper PD11.N_events.S_hypo PD11.N_events.m_hypo;
PD12.N_events.S_hyper PD12.N_events.m_hyper PD12.N_events.S_hypo PD12.N_events.m_hypo;
PD13.N_events.S_hyper PD13.N_events.m_hyper PD13.N_events.S_hypo PD13.N_events.m_hypo;
PD14.N_events.S_hyper PD14.N_events.m_hyper PD14.N_events.S_hypo PD14.N_events.m_hypo;
PD15.N_events.S_hyper PD15.N_events.m_hyper PD15.N_events.S_hypo PD15.N_events.m_hypo;
PD16.N_events.S_hyper PD16.N_events.m_hyper PD16.N_events.S_hypo PD16.N_events.m_hypo;
PD17.N_events.S_hyper PD17.N_events.m_hyper PD17.N_events.S_hypo PD17.N_events.m_hypo;
PD19.N_events.S_hyper PD19.N_events.m_hyper PD19.N_events.S_hypo PD19.N_events.m_hypo;
PD20.N_events.S_hyper PD20.N_events.m_hyper PD20.N_events.S_hypo PD20.N_events.m_hypo;
PD22.N_events.S_hyper PD22.N_events.m_hyper PD22.N_events.S_hypo PD22.N_events.m_hypo;
PD23.N_events.S_hyper PD23.N_events.m_hyper PD23.N_events.S_hypo PD23.N_events.m_hypo;
PD24.N_events.S_hyper PD24.N_events.m_hyper PD24.N_events.S_hypo PD24.N_events.m_hypo;
PD25.N_events.S_hyper PD25.N_events.m_hyper PD25.N_events.S_hypo PD25.N_events.m_hypo;
PD26.N_events.S_hyper PD26.N_events.m_hyper PD26.N_events.S_hypo PD26.N_events.m_hypo;
PD28.N_events.S_hyper PD28.N_events.m_hyper PD28.N_events.S_hypo PD28.N_events.m_hypo;
PD33.N_events.S_hyper PD33.N_events.m_hyper PD33.N_events.S_hypo PD33.N_events.m_hypo;
PD34.N_events.S_hyper PD34.N_events.m_hyper PD34.N_events.S_hypo PD34.N_events.m_hypo;
PD35.N_events.S_hyper PD35.N_events.m_hyper PD35.N_events.S_hypo PD35.N_events.m_hypo;
PD36.N_events.S_hyper PD36.N_events.m_hyper PD36.N_events.S_hypo PD36.N_events.m_hypo;
PD37.N_events.S_hyper PD37.N_events.m_hyper PD37.N_events.S_hypo PD37.N_events.m_hypo;
PD38.N_events.S_hyper PD38.N_events.m_hyper PD38.N_events.S_hypo PD38.N_events.m_hypo;
PD39.N_events.S_hyper PD39.N_events.m_hyper PD39.N_events.S_hypo PD39.N_events.m_hypo;
PD40.N_events.S_hyper PD40.N_events.m_hyper PD40.N_events.S_hypo PD40.N_events.m_hypo;
PD41.N_events.S_hyper PD41.N_events.m_hyper PD41.N_events.S_hypo PD41.N_events.m_hypo;
PD42.N_events.S_hyper PD42.N_events.m_hyper PD42.N_events.S_hypo PD42.N_events.m_hypo;
PD43.N_events.S_hyper PD43.N_events.m_hyper PD43.N_events.S_hypo PD43.N_events.m_hypo;
PD44.N_events.S_hyper PD44.N_events.m_hyper PD44.N_events.S_hypo PD44.N_events.m_hypo;
PD45.N_events.S_hyper PD45.N_events.m_hyper PD45.N_events.S_hypo PD45.N_events.m_hypo;
PD47.N_events.S_hyper PD47.N_events.m_hyper PD47.N_events.S_hypo PD47.N_events.m_hypo;
PD48.N_events.S_hyper PD48.N_events.m_hyper PD48.N_events.S_hypo PD48.N_events.m_hypo;
PD49.N_events.S_hyper PD49.N_events.m_hyper PD49.N_events.S_hypo PD49.N_events.m_hypo;
PD50.N_events.S_hyper PD50.N_events.m_hyper PD50.N_events.S_hypo PD50.N_events.m_hypo;
PD51.N_events.S_hyper PD51.N_events.m_hyper PD51.N_events.S_hypo PD51.N_events.m_hypo;
PD52.N_events.S_hyper PD52.N_events.m_hyper PD52.N_events.S_hypo PD52.N_events.m_hypo;
PD53.N_events.S_hyper PD53.N_events.m_hyper PD53.N_events.S_hypo PD53.N_events.m_hypo;
PD54.N_events.S_hyper PD54.N_events.m_hyper PD54.N_events.S_hypo PD54.N_events.m_hypo;
PD55.N_events.S_hyper PD55.N_events.m_hyper PD55.N_events.S_hypo PD55.N_events.m_hypo;
PD56.N_events.S_hyper PD56.N_events.m_hyper PD56.N_events.S_hypo PD56.N_events.m_hypo;
PD57.N_events.S_hyper PD57.N_events.m_hyper PD57.N_events.S_hypo PD57.N_events.m_hypo;
PD58.N_events.S_hyper PD58.N_events.m_hyper PD58.N_events.S_hypo PD58.N_events.m_hypo;
PD59.N_events.S_hyper PD59.N_events.m_hyper PD59.N_events.S_hypo PD59.N_events.m_hypo;
PD60.N_events.S_hyper PD60.N_events.m_hyper PD60.N_events.S_hypo PD60.N_events.m_hypo;];

%% duration of events
PD_hyperSm_hypoSm_eventlength_sum = [sum(PD1.events_length.S_hyper) sum(PD1.events_length.m_hyper) sum(PD1.events_length.S_hypo) sum(PD1.events_length.m_hypo);
sum(PD2.events_length.S_hyper) sum(PD2.events_length.m_hyper) sum(PD2.events_length.S_hypo) sum(PD2.events_length.m_hypo);
sum(PD3.events_length.S_hyper) sum(PD3.events_length.m_hyper) sum(PD3.events_length.S_hypo) sum(PD3.events_length.m_hypo);
sum(PD4.events_length.S_hyper) sum(PD4.events_length.m_hyper) sum(PD4.events_length.S_hypo) sum(PD4.events_length.m_hypo);
sum(PD5.events_length.S_hyper) sum(PD5.events_length.m_hyper) sum(PD5.events_length.S_hypo) sum(PD5.events_length.m_hypo);
sum(PD6.events_length.S_hyper) sum(PD6.events_length.m_hyper) sum(PD6.events_length.S_hypo) sum(PD6.events_length.m_hypo);
sum(PD7.events_length.S_hyper) sum(PD7.events_length.m_hyper) sum(PD7.events_length.S_hypo) sum(PD7.events_length.m_hypo);
sum(PD8.events_length.S_hyper) sum(PD8.events_length.m_hyper) sum(PD8.events_length.S_hypo) sum(PD8.events_length.m_hypo);
sum(PD9.events_length.S_hyper) sum(PD9.events_length.m_hyper) sum(PD9.events_length.S_hypo) sum(PD9.events_length.m_hypo);
sum(PD10.events_length.S_hyper) sum(PD10.events_length.m_hyper) sum(PD10.events_length.S_hypo) sum(PD10.events_length.m_hypo);
sum(PD11.events_length.S_hyper) sum(PD11.events_length.m_hyper) sum(PD11.events_length.S_hypo) sum(PD11.events_length.m_hypo);
sum(PD12.events_length.S_hyper) sum(PD12.events_length.m_hyper) sum(PD12.events_length.S_hypo) sum(PD12.events_length.m_hypo);
sum(PD13.events_length.S_hyper) sum(PD13.events_length.m_hyper) sum(PD13.events_length.S_hypo) sum(PD13.events_length.m_hypo);
sum(PD14.events_length.S_hyper) sum(PD14.events_length.m_hyper) sum(PD14.events_length.S_hypo) sum(PD14.events_length.m_hypo);
sum(PD15.events_length.S_hyper) sum(PD15.events_length.m_hyper) sum(PD15.events_length.S_hypo) sum(PD15.events_length.m_hypo);
sum(PD16.events_length.S_hyper) sum(PD16.events_length.m_hyper) sum(PD16.events_length.S_hypo) sum(PD16.events_length.m_hypo);
sum(PD17.events_length.S_hyper) sum(PD17.events_length.m_hyper) sum(PD17.events_length.S_hypo) sum(PD17.events_length.m_hypo);
sum(PD19.events_length.S_hyper) sum(PD19.events_length.m_hyper) sum(PD19.events_length.S_hypo) sum(PD19.events_length.m_hypo);
sum(PD20.events_length.S_hyper) sum(PD20.events_length.m_hyper) sum(PD20.events_length.S_hypo) sum(PD20.events_length.m_hypo);
sum(PD22.events_length.S_hyper) sum(PD22.events_length.m_hyper) sum(PD22.events_length.S_hypo) sum(PD22.events_length.m_hypo);
sum(PD23.events_length.S_hyper) sum(PD23.events_length.m_hyper) sum(PD23.events_length.S_hypo) sum(PD23.events_length.m_hypo);
sum(PD24.events_length.S_hyper) sum(PD24.events_length.m_hyper) sum(PD24.events_length.S_hypo) sum(PD24.events_length.m_hypo);
sum(PD25.events_length.S_hyper) sum(PD25.events_length.m_hyper) sum(PD25.events_length.S_hypo) sum(PD25.events_length.m_hypo);
sum(PD26.events_length.S_hyper) sum(PD26.events_length.m_hyper) sum(PD26.events_length.S_hypo) sum(PD26.events_length.m_hypo);
sum(PD28.events_length.S_hyper) sum(PD28.events_length.m_hyper) sum(PD28.events_length.S_hypo) sum(PD28.events_length.m_hypo);
sum(PD33.events_length.S_hyper) sum(PD33.events_length.m_hyper) sum(PD33.events_length.S_hypo) sum(PD33.events_length.m_hypo);
sum(PD34.events_length.S_hyper) sum(PD34.events_length.m_hyper) sum(PD34.events_length.S_hypo) sum(PD34.events_length.m_hypo);
sum(PD35.events_length.S_hyper) sum(PD35.events_length.m_hyper) sum(PD35.events_length.S_hypo) sum(PD35.events_length.m_hypo);
sum(PD36.events_length.S_hyper) sum(PD36.events_length.m_hyper) sum(PD36.events_length.S_hypo) sum(PD36.events_length.m_hypo);
sum(PD37.events_length.S_hyper) sum(PD37.events_length.m_hyper) sum(PD37.events_length.S_hypo) sum(PD37.events_length.m_hypo);
sum(PD38.events_length.S_hyper) sum(PD38.events_length.m_hyper) sum(PD38.events_length.S_hypo) sum(PD38.events_length.m_hypo);
sum(PD39.events_length.S_hyper) sum(PD39.events_length.m_hyper) sum(PD39.events_length.S_hypo) sum(PD39.events_length.m_hypo);
sum(PD40.events_length.S_hyper) sum(PD40.events_length.m_hyper) sum(PD40.events_length.S_hypo) sum(PD40.events_length.m_hypo);
sum(PD41.events_length.S_hyper) sum(PD41.events_length.m_hyper) sum(PD41.events_length.S_hypo) sum(PD41.events_length.m_hypo);
sum(PD42.events_length.S_hyper) sum(PD42.events_length.m_hyper) sum(PD42.events_length.S_hypo) sum(PD42.events_length.m_hypo);
sum(PD43.events_length.S_hyper) sum(PD43.events_length.m_hyper) sum(PD43.events_length.S_hypo) sum(PD43.events_length.m_hypo);
sum(PD44.events_length.S_hyper) sum(PD44.events_length.m_hyper) sum(PD44.events_length.S_hypo) sum(PD44.events_length.m_hypo);
sum(PD45.events_length.S_hyper) sum(PD45.events_length.m_hyper) sum(PD45.events_length.S_hypo) sum(PD45.events_length.m_hypo);
sum(PD47.events_length.S_hyper) sum(PD47.events_length.m_hyper) sum(PD47.events_length.S_hypo) sum(PD47.events_length.m_hypo);
sum(PD48.events_length.S_hyper) sum(PD48.events_length.m_hyper) sum(PD48.events_length.S_hypo) sum(PD48.events_length.m_hypo);
sum(PD49.events_length.S_hyper) sum(PD49.events_length.m_hyper) sum(PD49.events_length.S_hypo) sum(PD49.events_length.m_hypo);
sum(PD50.events_length.S_hyper) sum(PD50.events_length.m_hyper) sum(PD50.events_length.S_hypo) sum(PD50.events_length.m_hypo);
sum(PD51.events_length.S_hyper) sum(PD51.events_length.m_hyper) sum(PD51.events_length.S_hypo) sum(PD51.events_length.m_hypo);
sum(PD52.events_length.S_hyper) sum(PD52.events_length.m_hyper) sum(PD52.events_length.S_hypo) sum(PD52.events_length.m_hypo);
sum(PD53.events_length.S_hyper) sum(PD53.events_length.m_hyper) sum(PD53.events_length.S_hypo) sum(PD53.events_length.m_hypo);
sum(PD54.events_length.S_hyper) sum(PD54.events_length.m_hyper) sum(PD54.events_length.S_hypo) sum(PD54.events_length.m_hypo);
sum(PD55.events_length.S_hyper) sum(PD55.events_length.m_hyper) sum(PD55.events_length.S_hypo) sum(PD55.events_length.m_hypo);
sum(PD56.events_length.S_hyper) sum(PD56.events_length.m_hyper) sum(PD56.events_length.S_hypo) sum(PD56.events_length.m_hypo);
sum(PD57.events_length.S_hyper) sum(PD57.events_length.m_hyper) sum(PD57.events_length.S_hypo) sum(PD57.events_length.m_hypo);
sum(PD58.events_length.S_hyper) sum(PD58.events_length.m_hyper) sum(PD58.events_length.S_hypo) sum(PD58.events_length.m_hypo);
sum(PD59.events_length.S_hyper) sum(PD59.events_length.m_hyper) sum(PD59.events_length.S_hypo) sum(PD59.events_length.m_hypo);
sum(PD60.events_length.S_hyper) sum(PD60.events_length.m_hyper) sum(PD60.events_length.S_hypo) sum(PD60.events_length.m_hypo);];

PD_hyperSm_hypoSm_eventlength_sum = PD_hyperSm_hypoSm_eventlength_sum*5;

%% k means cluster into 4 groups
rng(1); % For reproducibility
idx_K_4 = kmeans(PD_hyperSm_hypoSm,4);
idx_K_2 = kmeans(PD_hyperSm_hypoSm,2);
idx_K_3 = kmeans(PD_hyperSm_hypoSm,3);
idx_K_5 = kmeans(PD_hyperSm_hypoSm,5);

idx_K_4_eventlength_sum = kmeans(PD_hyperSm_hypoSm_eventlength_sum,4);
idx_K_2_eventlength_sum = kmeans(PD_hyperSm_hypoSm_eventlength_sum,2);
idx_K_3_eventlength_sum = kmeans(PD_hyperSm_hypoSm_eventlength_sum,3);
idx_K_5_eventlength_sum = kmeans(PD_hyperSm_hypoSm_eventlength_sum,5);

[idx_K_6, C6] = kmeans([sum(PD_hyperSm_hypoSm(:,1:2),2) sum(PD_hyperSm_hypoSm(:,3:4),2)],6);
[idx_K_5, C5] = kmeans([sum(PD_hyperSm_hypoSm(:,1:2),2) sum(PD_hyperSm_hypoSm(:,3:4),2)],5);
[idx_K_4, C4] = kmeans([sum(PD_hyperSm_hypoSm(:,1:2),2) sum(PD_hyperSm_hypoSm(:,3:4),2)],4);
[idx_K_2, C2] = kmeans([sum(PD_hyperSm_hypoSm(:,1:2),2) sum(PD_hyperSm_hypoSm(:,3:4),2)],2);
[idx_K_3, C3] = kmeans([sum(PD_hyperSm_hypoSm(:,1:2),2) sum(PD_hyperSm_hypoSm(:,3:4),2)],3);

[idx_K_6_eventlength_sum, C6_eventlength_sum] = kmeans([sum(PD_hyperSm_hypoSm_eventlength_sum(:,1:2),2) sum(PD_hyperSm_hypoSm_eventlength_sum(:,3:4),2)],6);
[idx_K_5_eventlength_sum, C5_eventlength_sum] = kmeans([sum(PD_hyperSm_hypoSm_eventlength_sum(:,1:2),2) sum(PD_hyperSm_hypoSm_eventlength_sum(:,3:4),2)],5);
[idx_K_4_eventlength_sum, C4_eventlength_sum] = kmeans([sum(PD_hyperSm_hypoSm_eventlength_sum(:,1:2),2) sum(PD_hyperSm_hypoSm_eventlength_sum(:,3:4),2)],4);
[idx_K_2_eventlength_sum, C2_eventlength_sum] = kmeans([sum(PD_hyperSm_hypoSm_eventlength_sum(:,1:2),2) sum(PD_hyperSm_hypoSm_eventlength_sum(:,3:4),2)],2);
[idx_K_3_eventlength_sum, C3_eventlength_sum] = kmeans([sum(PD_hyperSm_hypoSm_eventlength_sum(:,1:2),2) sum(PD_hyperSm_hypoSm_eventlength_sum(:,3:4),2)],3);

a = [1:1:17 19 20 22:1:26 28 33:1:45 47:1:60  ]'; b = num2str(a); c = cellstr(b);
dx = 0.3; dy = 0.3; % displacement so the text does not overlay the data points

figure()
subplot(1,2,1)
scatter(sum(PD_hyperSm_hypoSm(:,1:2),2),sum(PD_hyperSm_hypoSm(:,3:4),2),'blue')
ylabel('Sum N Hypo S m')
xlabel('Sum N Hyper S m')
title("Sum N events")
xlim([0 35]);ylim([0 35]);
text(sum(PD_hyperSm_hypoSm(:,1:2),2)+dx, sum(PD_hyperSm_hypoSm(:,3:4),2)+dy, c);
subplot(1,2,2)
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(:,1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(:,3:4),2),'red')
ylabel('Sum Length Hypo S m / minutes')
xlabel('Sum Length Hyper S m / minutes')
title("Sum Length events / minutes")
xlim([0 5000]);ylim([0 5000]);
text(sum(PD_hyperSm_hypoSm_eventlength_sum(:,1:2),2)+dx, sum(PD_hyperSm_hypoSm_eventlength_sum(:,3:4),2)+dy, c);


idx_no_Sm_hyper= find(sum(PD_hyperSm_hypoSm(:,1:2),2)==0);
b = num2str(a(idx_no_Sm_hyper)); c = cellstr(b);

figure()
scatter(sum(PD_hyperSm_hypoSm(idx_no_Sm_hyper,1:2),2),sum(PD_hyperSm_hypoSm(idx_no_Sm_hyper,3:4),2),'blue')
ylabel('Sum N Hypo S m')
xlabel('Sum N Hyper S m')
title("Sum N events - Only Sm Hypo")
text(sum(PD_hyperSm_hypoSm(idx_no_Sm_hyper,1:2),2)+dx, sum(PD_hyperSm_hypoSm(idx_no_Sm_hyper,3:4),2)+dy, c);


%% K means cluster subplots K2:6 sum N
figure()
sgtitle("Sum N events K means Clustering")
subplot(2,3,1)
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_2==1),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_2==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_2==2),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_2==2),3:4),2),'cyan')
plot(C2(:,1),C2(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Sum N Hypo S m')
xlabel('Sum N Hyper S m')
title('K = 2')
xlim([0 35]);ylim([0 35]);

subplot(2,3,2)
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_3==1),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_3==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_3==2),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_3==2),3:4),2),'cyan')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_3==3),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_3==3),3:4),2),'red')
plot(C3(:,1),C3(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Sum N Hypo S m')
xlabel('Sum N Hyper S m')
title('K = 3')
xlim([0 35]);ylim([0 35]);

subplot(2,3,3)
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_4==1),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_4==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_4==2),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_4==2),3:4),2),'cyan')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_4==3),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_4==3),3:4),2),'red')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_4==4),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_4==4),3:4),2),'black')
plot(C4(:,1),C4(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Sum N Hypo S m')
xlabel('Sum N Hyper S m')
title('K = 4')
xlim([0 35]);ylim([0 35]);


subplot(2,3,4)
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_5==1),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_5==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_5==2),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_5==2),3:4),2),'cyan')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_5==3),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_5==3),3:4),2),'red')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_5==4),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_5==4),3:4),2),'black')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_5==5),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_5==5),3:4),2),'green')
plot(C5(:,1),C5(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Sum N Hypo S m')
xlabel('Sum N Hyper S m')
title('K = 5')
xlim([0 35]);ylim([0 35]);

subplot(2,3,5)
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_6==1),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_6==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_6==2),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_6==2),3:4),2),'cyan')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_6==3),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_6==3),3:4),2),'red')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_6==4),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_6==4),3:4),2),'black')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_6==5),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_6==5),3:4),2),'green')
scatter(sum(PD_hyperSm_hypoSm(find(idx_K_6==6),1:2),2),sum(PD_hyperSm_hypoSm(find(idx_K_6==6),3:4),2),'magenta')
plot(C6(:,1),C6(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Sum N Hypo S m')
xlabel('Sum N Hyper S m')
title('K = 6')
xlim([0 35]);ylim([0 35]);

%% K means cluster subplots K2:6 event length
figure()
sgtitle("Sum Events Length K means Clustering")
subplot(2,3,1)
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_2_eventlength_sum==1),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_2_eventlength_sum==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_2_eventlength_sum==2),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_2_eventlength_sum==2),3:4),2),'cyan')
plot(C2_eventlength_sum(:,1),C2_eventlength_sum(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Eventlength sum Hypo S m / minutes')
xlabel('Eventlength sum Hyper S m / minutes')
title('K = 2')

subplot(2,3,2)
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_3_eventlength_sum==1),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_3_eventlength_sum==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_3_eventlength_sum==2),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_3_eventlength_sum==2),3:4),2),'cyan')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_3_eventlength_sum==3),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_3_eventlength_sum==3),3:4),2),'red')
plot(C3_eventlength_sum(:,1),C3_eventlength_sum(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Eventlength sum Hypo S m / minutes')
xlabel('Eventlength sum Hyper S m / minutes')
title('K = 3')

subplot(2,3,3)
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_4_eventlength_sum==1),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_4_eventlength_sum==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_4_eventlength_sum==2),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_4_eventlength_sum==2),3:4),2),'cyan')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_4_eventlength_sum==3),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_4_eventlength_sum==3),3:4),2),'red')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_4_eventlength_sum==4),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_4_eventlength_sum==4),3:4),2),'black')
plot(C4_eventlength_sum(:,1),C4_eventlength_sum(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Eventlength sum Hypo S m / minutes')
xlabel('Eventlength sum Hyper S m / minutes')
title('K = 4')


subplot(2,3,4)
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==1),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==2),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==2),3:4),2),'cyan')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==3),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==3),3:4),2),'red')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==4),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==4),3:4),2),'black')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==5),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_5_eventlength_sum==5),3:4),2),'green')
plot(C5_eventlength_sum(:,1),C5_eventlength_sum(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Eventlength sum Hypo S m / minutes')
xlabel('Eventlength sum Hyper S m / minutes')
title('K = 5')

subplot(2,3,5)
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==1),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==1),3:4),2),'blue')
hold on
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==2),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==2),3:4),2),'cyan')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==3),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==3),3:4),2),'red')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==4),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==4),3:4),2),'black')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==5),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==5),3:4),2),'green')
scatter(sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==6),1:2),2),sum(PD_hyperSm_hypoSm_eventlength_sum(find(idx_K_6_eventlength_sum==6),3:4),2),'magenta')
plot(C6_eventlength_sum(:,1),C6_eventlength_sum(:,2),'kx',MarkerSize=16,LineWidth=3)
ylabel('Eventlength sum Hypo S m / minutes')
xlabel('Eventlength sum Hyper S m / minutes')
title('K = 6')

%%
PD_hyperSm_hypoSm(:,3:4) = -1* PD_hyperSm_hypoSm(:,3:4);
figure()
bar(PD_N,PD_hyperSm_hypoSm,'stacked')
legend('S Hyper','m Hyper','S Hypo','m Hypo')
xlabel('PD Number');
ylabel('N events');
grid on

%% mhypo or shypo

%PD 10 13 15 25	33 8 28 34 39 41 45 48 49 55

length_m_hypo = [PD10.events_length.m_hypo;
PD13.events_length.m_hypo;
PD15.events_length.m_hypo;
PD25.events_length.m_hypo;
PD33.events_length.m_hypo;
PD8.events_length.m_hypo;
PD28.events_length.m_hypo;
PD34.events_length.m_hypo;
PD39.events_length.m_hypo;
PD41.events_length.m_hypo;
PD45.events_length.m_hypo;
PD48.events_length.m_hypo;
PD49.events_length.m_hypo;
PD55.events_length.m_hypo;];

length_S_hypo = [PD10.events_length.S_hypo;
PD13.events_length.S_hypo;
PD15.events_length.S_hypo;
PD25.events_length.S_hypo;
PD33.events_length.S_hypo;
PD8.events_length.S_hypo;
PD28.events_length.S_hypo;
PD34.events_length.S_hypo;
PD39.events_length.S_hypo;
PD41.events_length.S_hypo;
PD45.events_length.S_hypo;
PD48.events_length.S_hypo;
PD49.events_length.S_hypo;
PD55.events_length.S_hypo;];

mean(length_m_hypo)
std(length_m_hypo)

length_S_m_hypo = [length_S_hypo;length_m_hypo];
length_S_m_hypo = nonzeros(length_S_m_hypo);
mean(length_S_m_hypo)*5
std(length_S_m_hypo)*5


%%
% figure()
% b = bar(PD_N,PD_hyperSm_hypoSm,'stacked');
% xtips1 = PD_N;
% ytips1 = PD_hyperSm_hypoSm;
% labels1 = string(PD_N);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
%

% PD1.events_length.m_hyper
% PD2.events_length.m_hyper
% PD3.events_length.m_hyper
% PD4.events_length.m_hyper



%find(PD_hyperSm_hypoSm(:,1)


%%
PD_hyperSm_hypoSm = [PD1 ;PD2;PD3;PD4;PD5;PD6;PD7;PD8;PD9;PD10; PD12;PD13;PD14;PD15;PD16;PD17;PD19;PD20;PD23;PD24;PD25;PD26;PD28;PD33;PD34;PD35;PD41;PD42;PD45;PD48;PD49;PD50;PD53;PD55 ];

PD_hyperSm_hypoSm(:,3:4) = -1* PD_hyperSm_hypoSm(:,3:4);

figure()
bar(PD_hyperSm_hypoSm,'stacked')
legend('S Hyper','m Hyper','S Hypo','m Hypo')
xlabel('PD Number');
ylabel('N events');
%% bad
type_spreadsheet = 2;
name_spreadsheet = 'PD11.xlsx';
PD11 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%%
name_spreadsheet = 'PD22.xlsx';
PD22 = readglucose_number(name_spreadsheet)

name_spreadsheet = 'PD23.xlsx';
PD23 = readglucose_number(name_spreadsheet)
%% sort I think to do with the date being read as string and not date type

name_spreadsheet = 'PD36.xlsx'; %type 2
%PD36 = readglucose_number(name_spreadsheet)
PD36 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%
name_spreadsheet = 'PD37.xlsx'; %type 2
%PD37 = readglucose_number(name_spreadsheet)
PD37 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %date needed to be converted to datetime from string
%% not working
name_spreadsheet = 'PD38.xlsx';
PD38 = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%% working
name_spreadsheet = 'PD39.xlsx'; %type 2
PD39 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %date needed to be converted to datetime from str
%
name_spreadsheet = 'PD40.xlsx'; %type 2
PD40 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %date needed to be converted to datetime from str
%
name_spreadsheet = 'PD43.xlsx'; %type 2
PD43 = readglucose_string_number(name_spreadsheet,type_spreadsheet)  %date needed to be converted to datetime from str

name_spreadsheet = 'PD47.xlsx'; %type 2
PD47 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %dates were in mixed US and UK format
%% sort time reading
name_spreadsheet = 'PD44.xlsx';
PD44 = readglucose_number(name_spreadsheet)

%% sort date switches US UK done
name_spreadsheet = 'PD47.xlsx'; %type 2
PD47 = readglucose_string_number(name_spreadsheet,type_spreadsheet) %dates were in mixed US and UK format
%% sort 
name_spreadsheet = 'PD51.xlsx';
PD51 = readglucose_number(name_spreadsheet)

%% sort
name_spreadsheet = 'PD52.xlsx';
PD52 = readglucose_number(name_spreadsheet)

name_spreadsheet = 'PD54.xlsx';
PD54 = readglucose_number(name_spreadsheet)
%%
clear all
PD_num = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17  19 20  22 24 25 26 27 28 33 34 35 36 37 38 39 40 41 42 43 44 45 47 48 49 50 51 52 53 54 55];

for i=1:size(PD_num,2)
%PD_all(i,1:4) =
 readglucose_number("PD"+num2str(i)+".xlsx")
end

%%

aaa = testfun(name_spreadsheet)

function testindex = testfun(name_spreadsheet)
        testindex.field1 = 1;
        testindex.field2 = [4 ; 5 ; 6];

end

%%
function glucose_events_info = readglucose_string_number(name_spreadsheet,type_spreadsheet)
%%
S_hyperG = 180;
m_hyperG = 144;

S_hypoG =47;
m_hypoG =72;


%name_spreadsheet = convertCharsToStrings(name_spreadsheet);
spreadsheet_data = readcell(name_spreadsheet);

%spreadsheet_data = readcell("PD10.xlsx");

%PD1.xlsx
%spreadsheet_data = readcell('PD1.xlsx');

%convert PD37 date from string to datetime
if name_spreadsheet == "PD60.xlsx" | name_spreadsheet == "PD59.xlsx" | name_spreadsheet == "PD56.xlsx" | name_spreadsheet == "PD38.xlsx" | name_spreadsheet == "PD22.xlsx" |name_spreadsheet == "PD37.xlsx" | name_spreadsheet == "PD38.xlsx" | name_spreadsheet == "PD39.xlsx" | name_spreadsheet == "PD40.xlsx" | name_spreadsheet == "PD43.xlsx" | name_spreadsheet == "PD51.xlsx"
    for i=1:size(spreadsheet_data,1)
        spreadsheet_data{i,1} = datetime(spreadsheet_data{i,1},'InputFormat','dd/MM/yyyy');
    end
end

%for PD51, need to convert time from 0.04 to an actual 24HR time

%glucose saved as string
%spreadsheet_data = readcell('test_time.xlsx');
if type_spreadsheet ==1
    for i=1:size(spreadsheet_data,1)
        glucoseC(i)= (str2double (cell2mat(spreadsheet_data(i,3)))/100)' ; % div 100 as glucose values have 2 DP of 0.
        %glucoseC(i)= ( (cell2mat(spreadsheet_data(i,3))))' ; % if glucose
        %saved as number
        date(i)=spreadsheet_data{i,1};
        time(i)=spreadsheet_data{i,2};
    end
end


%datetime(spreadsheet_data{i,1},'InputFormat','dd/MM/yyyy')

% if glucose
%saved as number
if type_spreadsheet ==2
    for i=1:size(spreadsheet_data,1)
        %glucoseC(i)= (str2double (cell2mat(spreadsheet_data(i,3)))/100)' ; % div 100 as glucose values have 2 DP of 0.
        glucoseC(i)= ( (cell2mat(spreadsheet_data(i,3))))' ; % if glucose
        %saved as number
        date(i)=spreadsheet_data{i,1};
        time(i)=spreadsheet_data{i,2};
    end
end


glucoseC=glucoseC';
%date(1:size(spreadsheet_data,1),1)=spreadsheet_data{1:size(spreadsheet_data,1),1};
%time(1:size(spreadsheet_data,1),1)=spreadsheet_data{1:size(spreadsheet_data,1),2};
date=date';
time= time';


delta_t = 1.157407407407407e-05; %1/ s in a day
s_in_day = 86400;
scale_t = 1/delta_t;
time = time*scale_t;

days = unique(date);
for i=1:size(date,1)
    dayN(i) = find(date(i)==days);
end
dayN = dayN';

universal_time = time + (s_in_day*(dayN-1));

if universal_time(1) > universal_time(end) %this means time is going backwards so we need to flip it
    universal_time = flip(universal_time);
    glucoseC = flip(glucoseC);
    date = flip(date);
    time = flip(time);
end

% Eliminate glucose values that are 
%a. duplicate readings
%b. not part of continuous every 5 min sampling

%297s not 300s as there might be small rounding errors.
rows_dup = find(round(diff(universal_time))<297); %finding a. and b. since gap > 303 classed as gaps
if size(rows_dup,1) > 0 %if there are any duplicates at all
    %keep 1st of each cluster, remove 2nd-N elements
    dup_remove = [false ; diff(rows_dup)==1]; %0 = KEEP , 1=Remove (%keeping 1st of each cluster, removing the n=2 to n=max nth elements of cluster

    if size(dup_remove,1)==1 %only 1 duplicate meas. so get rid of element after it
        dup_remove = true ;
    else
        diff_rows_dup = diff(rows_dup);
        %set single elements in rows dup to false, to get rid of. 
        %i.e elements not in a cluster..

        %check 1st element
        %if diff of 1st element > 1, it should be removed, set to true (1)
        if diff_rows_dup(1) >1
            dup_remove(1) = true;
        end
        %check last element for same check
        if diff_rows_dup(end) >1
            dup_remove(end) = true;
        end

        %Check all middle elements
        diff_rows_dup_zerod = diff_rows_dup-1;
        for i=1:size(diff_rows_dup,1)-1
            if nnz(diff_rows_dup_zerod(i:i+1)) ==2

                dup_remove(i+1) = true; %remove this element
            end
        end

    end
    %%%%remove row indices

    %dup_remove contains FALSE (0) (KEEP) and TRUE (1) (REMOVE) idxs
    %rows_dup contains the row number from the glucose data to keep or remove
    remove_row = rows_dup(dup_remove);
    %keep_row = rows_dup(~dup_remove) ;

    all_rows = 1:1:size(glucoseC,1);

    all_rows(remove_row) = 0;
    all_rows_keep = nonzeros(all_rows);

    %edit these variables
    glucoseC = glucoseC(all_rows_keep,:);
    date = date(all_rows_keep,:);
    time = time(all_rows_keep,:);
    spreadsheet_data = spreadsheet_data(all_rows_keep,:);
    universal_time = universal_time (all_rows_keep,:);
    dayN = dayN(all_rows_keep,:);
end
%
%find jumps in time
% for i=1:size(universal_time,1)-2
%     if round(universal_time(i,1) - universal_time(i+1,1)) ~= -300  & round(universal_time(i,1) - universal_time(i+1,1)) ~= -240
%         jump_start(i) = i;    
%     end
%     %if universal_time(i,1) - universal_time(i+1,1) >  - 400
%     %    cont_reading(i) = i;    
%     %end
%     if round(universal_time(end-i,1) - universal_time(end-i-1,1)) ~= 300 & round(universal_time(end-i,1) - universal_time(end-i-1,1)) ~= 240
%         jump_end(i) = i;    
%     end
% end

% for i=1:size(universal_time,1)-2
%     %if (round(universal_time(i,1) - universal_time(i+1,1)) < -305 | round(universal_time(i,1) - universal_time(i+1,1)) > -295 ) & ( round(universal_time(i,1) - universal_time(i+1,1)) < -245 | round(universal_time(i,1) - universal_time(i+1,1)) > -235) 
%     if (round(universal_time(i,1) - universal_time(i+1,1)) < -305 | round(universal_time(i,1) - universal_time(i+1,1)) > -295 )  % & ( round(universal_time(i,1) - universal_time(i+1,1)) < -245 | round(universal_time(i,1) - universal_time(i+1,1)) > -235) 
% 
%         jump_start(i) = i;    
%     end
%     %if universal_time(i,1) - universal_time(i+1,1) >  - 400
%     %    cont_reading(i) = i;    
%     %end
%     %if (round(universal_time(end-i,1) - universal_time(end-i-1,1)) < 295 | round(universal_time(end-i,1) - universal_time(end-i-1,1)) > 305) &  (round(universal_time(end-i,1) - universal_time(end-i-1,1)) < 235 | round(universal_time(end-i,1) - universal_time(end-i-1,1)) > 245)
%     if (round(universal_time(end-i,1) - universal_time(end-i-1,1)) < 295 | round(universal_time(end-i,1) - universal_time(end-i-1,1)) > 305) %&  (round(universal_time(end-i,1) - universal_time(end-i-1,1)) < 235 | round(universal_time(end-i,1) - universal_time(end-i-1,1)) > 245)
% 
%         jump_end(i) = i;    
%     end
% end

% new gap finding
for i=1:size(universal_time,1)-2
    delta_t_gap(i,1) = universal_time(i+1,1)- universal_time(i,1);
end

%delta_t_gap_idx = find( (delta_t_gap < 297)| (delta_t_gap > 303) );
delta_t_gap_idx = find( (delta_t_gap > 303) );

size(delta_t_gap_idx,1)>0;

if size(delta_t_gap_idx,1)>0
    for i=1:size(delta_t_gap_idx,1)
        delta_t_gap_length(i,1) =  -( universal_time(delta_t_gap_idx(i)) - universal_time(delta_t_gap_idx(i)+1));
    end

    %delta_t_gap_Nsamples = round(delta_t_gap_length./300) -1;
    %delta_t_gap_Nsamples = floor(delta_t_gap_length./300) -1; %16 04 24
    %delta_t_gap_Nsamples = floor(round(delta_t_gap_length./300)) -1; %16 04 24

    %delta_t_gap_length./300
    %delta_t_gap_Nsamples = floor(ceil(delta_t_gap_length./300)) -1; %16 04 24

    %delta_t_gap_Nsamples = floor(ceil(delta_t_gap_length./300)) -1; %16 04 24

    %delta_t_gap_Nsamples = ceil((round(delta_t_gap_length)./300)-1);

    delta_t_gap_Nsamples = (floor(round(delta_t_gap_length)./300)) -1; %if gap <600, rtn 0 gap. gap 600 rtn 1 gap, 

    delta_t_gap_start=delta_t_gap_idx;
    delta_t_gap_end=delta_t_gap_idx+1;

    universal_time(delta_t_gap_start);
    universal_time(delta_t_gap_end);
end

if size(delta_t_gap_idx,1)>0
    jump_start = delta_t_gap_start;
    jump_end = delta_t_gap_end;

    %cont_reading = nonzeros(cont_reading);

    %for i=1:size(jump,1)
    %    length_jump(i,1) = universal_time(jump(i)) - universal_time(jump(i)+1);
    %end

    length_jump = delta_t_gap_Nsamples;
    start_gap = delta_t_gap_start;
    end_gap = delta_t_gap_end;
else
    jump = 0;
    length_jump =0;
    start_gap = 0;
    end_gap = 0;
    delta_t_gap_Nsamples = 0;
    delta_t_gap_length = 0;
end

%% Interpolate gaps 15 04 24

%number of time points needed to fill the gap. (i.e it doesn't include the
%time point of the start/end of the gap)

if start_gap(1)>0
for i=1:size(delta_t_gap_Nsamples,1)%%%%
if delta_t_gap_Nsamples(i) > 0
    %%
    %delta_t_gap_Nsamples(1);

    %interpolate glucoseC between the gaps by getting the mean
    glucoseC_mean_btwn_gaps = round(mean([glucoseC(start_gap(i)) glucoseC(end_gap(i))]));
    glucoseC_int_gaps = ones(delta_t_gap_Nsamples(i),1)*glucoseC_mean_btwn_gaps;

    %Universal time
    universal_time_int_gaps = [round(universal_time(start_gap(i))):300:round(universal_time(end_gap(i)))]';
    universal_time_int_gaps = universal_time_int_gaps(2:end-1);

    %Time

    %special case if the gap will just be 1 interpolated value
    

    %if time(start_gap(i)) < time(end_gap(i)) %time hasn't looped behind due to passing midnight i.e 5pm to 3pm = 22hrs passed
    %    time_int_gaps=[time(start_gap(i)):300:time(end_gap(i))]';
    %    time_int_gaps = time_int_gaps(2:end-1);
    %end

    %if time(start_gap(i)) > time(end_gap(i)) %time has looped behind due to passing midnight i.e 5pm to 3pm = 22hrs passed
        %time_int_gaps=[time(start_gap(i)):300:time(end_gap(i))]';
        %time_int_gaps = time_int_gaps(2:end-1);

        time_int_gaps = [round(time(start_gap(i))):300:round(time(start_gap(i)))+ (300*delta_t_gap_Nsamples(i))+300]'; %maybe this should be done for both cases, as it will be corrected in the next step anyway
        time_int_gaps = time_int_gaps(2:end-1);
    %end
    %delta_t_gap_Nsamples(1)
    %size of interpolated time points
    dayN_gap = zeros(size(time_int_gaps,1),1);


    %if the start/end of gap are on the same day
    if dayN(end_gap(i)) - dayN(start_gap(i)) == 0
        dayN_gap(:) = 0;
    end

    %if the start/end of gap are on different days
    if dayN(end_gap(i)) - dayN(start_gap(i)) > 0
        for j=1:dayN(end_gap(i)) - dayN(start_gap(i)) %for number of days gap , 1,2,3 etc.
            dayN_gap( find(time_int_gaps > s_in_day*j)) = j; %if the time_int_gaps crosses a multiple of 24hr (seconds) it is day J
        end
    end

    %correct dayN int
    dayN_int = ones(size(time_int_gaps,1),1)*dayN(start_gap(i)); %starting value for the dayN interpolation
    dayN_int = dayN_int + dayN_gap; %adding the delta dayN gap, for times that go past midnight

    time_int_gaps = time_int_gaps - (s_in_day*dayN_gap); %correcting for new days, by subtracting 24hrs from times from new days

    %correct date int
    date_int = zeros(size(time_int_gaps,1),1) + date(start_gap(i)); %starting value for the date interpolation
    date_int = date_int + hours(24)*dayN_gap;


% delta_t_gap_length(1)/s_in_day
%     for j=1:ceil(delta_t_gap_length(1)/s_in_day)
%         dayN_gap = find(time_int_gaps > s_in_day*j)
%     end


    %check time doesn't cross midnight 

    %Date - need to check that gap doesn't cross midnight
    %date(start_gap(i))
    %date(end_gap(i))

    %day N - '' make sure doesn't cross midnight

    %spreadsheet_data - '' make sure doesn't cross midnight

    % Append data %Old data up to gap + New inter data + old data after gap
    new_glucoseC = [glucoseC(1:start_gap(i)) ;  glucoseC_int_gaps ; glucoseC(end_gap(i):end) ];
    new_time = [time(1:start_gap(i)); time_int_gaps;  time(end_gap(i):end)];
    new_dayN = [dayN(1:start_gap(i)) ; dayN_int ; dayN(end_gap(i):end)  ];
    new_universal_time = [universal_time(1:start_gap(i)) ; universal_time_int_gaps ; universal_time(end_gap(i):end) ];
    new_date = [date(1:start_gap(i)) ; date_int ; date(end_gap(i):end) ];

    %update new spreadsheet data
    %create interpolated spreedsheet data
    %spreadsheet_int = zeros(size(new_time,1),3);
    spreadsheet_int = cell(size(time_int_gaps,1),3);

    for k=1:size(spreadsheet_int,1)
        spreadsheet_int{k, 1}  = date_int(k);
        spreadsheet_int{k,2} = time_int_gaps(k)/s_in_day;
        spreadsheet_int{k,3} = glucoseC_int_gaps(k); %glucoseC
    end

    %now append and make the new spreadsheets
    new_spreadsheet_data = cell(size(new_time,1),3);
    for k=1:start_gap(i)
        new_spreadsheet_data{k,1} = spreadsheet_data{k, 1};
        new_spreadsheet_data{k,2} = spreadsheet_data{k, 2};
        new_spreadsheet_data{k,3} = spreadsheet_data{k, 3};
    end
    %middle of spreadsheet data
    for k=1:delta_t_gap_Nsamples(i)
        new_spreadsheet_data{k+start_gap(i),1} = spreadsheet_int{k, 1};
        new_spreadsheet_data{k+start_gap(i),2} = spreadsheet_int{k, 2};
        new_spreadsheet_data{k+start_gap(i),3} = spreadsheet_int{k, 3};
    end
    for k=end_gap(i)+delta_t_gap_Nsamples(i):size(new_time,1)
        new_spreadsheet_data{k,1} = spreadsheet_data{k-delta_t_gap_Nsamples(i), 1};
        new_spreadsheet_data{k,2} = spreadsheet_data{k-delta_t_gap_Nsamples(i), 2};
        new_spreadsheet_data{k,3} = spreadsheet_data{k-delta_t_gap_Nsamples(i), 3};
    end

    %update gap end index pointer (and start pointers of next gaps)
    end_gap(i:end) = end_gap(i:end) + delta_t_gap_Nsamples(i);%+1; %+1;
    start_gap(i+1:end) = start_gap(i+1:end) + delta_t_gap_Nsamples(i);%+1; %+1;
    
    %update new values
    glucoseC = new_glucoseC;
    time = new_time;
    dayN = new_dayN;
    universal_time = new_universal_time;
    date = new_date;
    spreadsheet_data = new_spreadsheet_data;
    %%
end
end
    
end

glucose_events_info.gap_info.start_gap = start_gap;
glucose_events_info.gap_info.end_gap = end_gap;
glucose_events_info.gap_info.delta_t_gap_Nsamples = delta_t_gap_Nsamples;
glucose_events_info.gap_info.delta_t_gap_length = delta_t_gap_length;
%variables that need updating
%glucoseC = glucoseC(all_rows_keep,:);
%date = date(all_rows_keep,:);
%time = time(all_rows_keep,:);
%spreadsheet_data = spreadsheet_data(all_rows_keep,:);
%universal_time = universal_time (all_rows_keep,:);
%dayN = dayN(all_rows_keep,:);
%% Eliminate glucose values that are (12 04 24) again - after gap filling
%a. duplicate readings
%b. not part of continuous every 5 min sampling

clear vars rows_dup dup_remove diff_rows_dup diff_rows_dup_zerod remove_row all_rows all_rows_keep;

%297s not 300s as there might be small rounding errors.
rows_dup = find(round(diff(universal_time))<297); %finding a. and b. since gap > 303 classed as gaps
if size(rows_dup,1) > 0 %if there are any duplicates at all
    %keep 1st of each cluster, remove 2nd-N elements
    dup_remove = [false ; diff(rows_dup)==1]; %0 = KEEP , 1=Remove (%keeping 1st of each cluster, removing the n=2 to n=max nth elements of cluster
    
    if size(dup_remove,1)==1 %only 1 duplicate meas. so get rid of element after it
        dup_remove = true ;
    else
        diff_rows_dup = diff(rows_dup);
        %set single elements in rows dup to false, to get rid of. 
        %i.e elements not in a cluster..

        %check 1st element
        %if diff of 1st element > 1, it should be removed, set to true (1)
        if diff_rows_dup(1) >1
            dup_remove(1) = true;
        end
        %check last element for same check
        if diff_rows_dup(end) >1
            dup_remove(end) = true;
        end

        %Check all middle elements
        diff_rows_dup_zerod = diff_rows_dup-1;
        for i=1:size(diff_rows_dup,1)-1
            if nnz(diff_rows_dup_zerod(i:i+1)) ==2
            
                dup_remove(i+1) = true; %remove this element
            end
        end

    end
    %%%%remove row indices

    %dup_remove contains FALSE (0) (KEEP) and TRUE (1) (REMOVE) idxs
    %rows_dup contains the row number from the glucose data to keep or remove
    remove_row = rows_dup(dup_remove);
    %keep_row = rows_dup(~dup_remove) ;

    all_rows = 1:1:size(glucoseC,1);

    all_rows(remove_row) = 0;
    all_rows_keep = nonzeros(all_rows);

    %edit these variables
    glucoseC = glucoseC(all_rows_keep,:);
    date = date(all_rows_keep,:);
    time = time(all_rows_keep,:);
    spreadsheet_data = spreadsheet_data(all_rows_keep,:);
    universal_time = universal_time (all_rows_keep,:);
    dayN = dayN(all_rows_keep,:);
end
%%
%change above 
%new gap finding end

% find( (delta_t_gap < 235)| (delta_t_gap > 245) )
% 
% unique(delta_t_gap)
% 
% 
% 320 > 305 | 320 < 295 | 320 > 245 | 320 < 235
% 
% (306 > 305 | 306 < 295) & (306 > 245 | 306 < 235)
% 
% 295 < 306 < 305

% jump = 1;
% 
% if exist('jump_start','var')
%     jump_start = nonzeros(jump_start);
%     jump_end = nonzeros(jump_end);
% 
%     %cont_reading = nonzeros(cont_reading);
% 
%     %for i=1:size(jump,1)
%     %    length_jump(i,1) = universal_time(jump(i)) - universal_time(jump(i)+1);
%     %end
% 
%     length_jump = delta_t_gap_Nsamples.
%     %start_gap = jump(1:2:end-1);
%     %end_gap = jump(2:2:end-1);
%     start_gap = jump_start;
%     end_gap = jump_end;
% else
%     jump = 0;
%     length_jump =0;
%     start_gap = 0;
%     end_gap = 0;
% end



S_hyper_point = find(glucoseC>S_hyperG);
%m_hyper_point = find(glucoseC<S_hyperG & glucoseC>=m_hyperG);
m_hyper_point = find(glucoseC>m_hyperG);


S_hypo_point = find(glucoseC<S_hypoG);
%m_hypo_point = find(glucoseC>S_hypoG & glucoseC<=m_hypoG);
m_hypo_point = find(glucoseC<m_hypoG);

euglycemia_point = find(glucoseC<=m_hyperG & glucoseC>=m_hypoG);

for i=1:size(m_hyper_point,1)-2
    if m_hyper_point(i) == m_hyper_point(i+1)-1 && m_hyper_point(i) == m_hyper_point(i+2)-2
        keep_m_hyper_point(i,1) = m_hyper_point(i);
    end 
end
for i=1:size(m_hyper_point,1)-2
    if m_hyper_point(end-i+1) == m_hyper_point(end-i)+1 && m_hyper_point(end-i+1) == m_hyper_point(end-i-1)+2
        keep_m_hyper_point2(i,1) = m_hyper_point(end-i+1);
    end 
end
for i=2:size(m_hyper_point,1)-1
    if m_hyper_point(i) == m_hyper_point(i-1)+1 && m_hyper_point(i) == m_hyper_point(i+1)-1
        keep_m_hyper_point3(i,1) = m_hyper_point(i);
    end 
end
if exist('keep_m_hyper_point','var')
    if exist('keep_m_hyper_point3','var')
        keep_m_hyper_point = unique([keep_m_hyper_point ; keep_m_hyper_point2 ; keep_m_hyper_point3]);
        m_hyper_point = nonzeros(keep_m_hyper_point);
    end
end
if ~exist('keep_m_hyper_point','var') & ~exist('keep_m_hyper_point2','var') & ~exist('keep_m_hyper_point3','var')
    m_hyper_point = 0;
end

for i=1:size(S_hyper_point,1)-2
    if S_hyper_point(i) == S_hyper_point(i+1)-1 && S_hyper_point(i) == S_hyper_point(i+2)-2
        keep_S_hyper_point(i,1) = S_hyper_point(i);
    end 
end
for i=1:size(S_hyper_point,1)-2
    if S_hyper_point(end-i+1) == S_hyper_point(end-i)+1 && S_hyper_point(end-i+1) == S_hyper_point(end-i-1)+2
        keep_S_hyper_point2(i,1) = S_hyper_point(end-i+1);
    end 
end
for i=2:size(S_hyper_point,1)-1
    if S_hyper_point(i) == S_hyper_point(i-1)+1 && S_hyper_point(i) == S_hyper_point(i+1)-1
        keep_S_hyper_point3(i,1) = S_hyper_point(i);
    end 
end
if exist('keep_S_hyper_point','var')
    if exist('keep_S_hyper_point3','var')
        keep_S_hyper_point = unique([keep_S_hyper_point ; keep_S_hyper_point2 ; keep_S_hyper_point3]);
        S_hyper_point = nonzeros(keep_S_hyper_point);
    end
end
if ~exist('keep_S_hyper_point','var') & ~exist('keep_S_hyper_point2','var') & ~exist('keep_S_hyper_point3','var')
    S_hyper_point = 0;
end

for i=1:size(m_hypo_point,1)-2
    if m_hypo_point(i) == m_hypo_point(i+1)-1 && m_hypo_point(i) == m_hypo_point(i+2)-2
        keep_m_hypo_point(i,1) = m_hypo_point(i);
    end 
end
for i=1:size(m_hypo_point,1)-2
    if m_hypo_point(end-i+1) == m_hypo_point(end-i)+1 && m_hypo_point(end-i+1) == m_hypo_point(end-i-1)+2
        keep_m_hypo_point2(i,1) = m_hypo_point(end-i+1);
    end 
end
for i=2:size(m_hypo_point,1)-1
    if m_hypo_point(i) == m_hypo_point(i-1)+1 && m_hypo_point(i) == m_hypo_point(i+1)-1
        keep_m_hypo_point3(i,1) = m_hypo_point(i);
    end 
end
if exist('keep_m_hypo_point','var')
    if exist('keep_m_hypo_point3','var')
        keep_m_hypo_point = unique([keep_m_hypo_point ; keep_m_hypo_point2 ;keep_m_hypo_point3]);
        m_hypo_point = nonzeros(keep_m_hypo_point);
    end
end
if ~exist('keep_m_hypo_point','var') & ~exist('keep_m_hypo_point2','var') & ~exist('keep_m_hypo_point3','var')
    m_hypo_point = 0;
end
for i=1:size(S_hypo_point,1)-2
    if S_hypo_point(i) == S_hypo_point(i+1)-1 && S_hypo_point(i) == S_hypo_point(i+2)-2
        keep_S_hypo_point(i,1) = S_hypo_point(i);
    end 
end
for i=1:size(S_hypo_point,1)-2
    if S_hypo_point(end-i+1) == S_hypo_point(end-i)+1 && S_hypo_point(end-i+1) == S_hypo_point(end-i-1)+2
        keep_S_hypo_point2(i,1) = S_hypo_point(end-i+1);
    end 
end
for i=2:size(S_hypo_point,1)-1
    if S_hypo_point(i) == S_hypo_point(i-1)+1 && S_hypo_point(i) == S_hypo_point(i+1)-1
        keep_S_hypo_point3(i,1) = S_hypo_point(i);
    end 
end
if exist('keep_S_hypo_point','var')
    if exist('keep_S_hypo_point3','var')
        keep_S_hypo_point = unique([keep_S_hypo_point ; keep_S_hypo_point2 ;keep_S_hypo_point3]);
        S_hypo_point = nonzeros(keep_S_hypo_point);
    end
end
if ~exist('keep_S_hypo_point','var') & ~exist('keep_S_hypo_point2','var') & ~exist('keep_S_hypo_point3','var')
    S_hypo_point = 0;
end

%%
%%%% m hyper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (size(m_hyper_point,1)>0)
    for i=1:size(m_hyper_point,1)-2 %was -3
        if m_hyper_point(i)>3
            j=0;
            found = 0;
            while found == 0 & (m_hyper_point(i)-4-j >2)
                    if ( (glucoseC(m_hyper_point(i)-1-j)<m_hyperG & glucoseC(m_hyper_point(i)-1-j)>m_hypoG) ...
                        ...
                        ...
                        & (glucoseC(m_hyper_point(i)-2-j)<m_hyperG & glucoseC(m_hyper_point(i)-2-j)>m_hypoG) ...
                        ...
                        & (glucoseC(m_hyper_point(i)-3-j)<m_hyperG & glucoseC(m_hyper_point(i)-3-j)>m_hypoG)  ...
                        ...
                        & (glucoseC(m_hyper_point(i)-4-j)<m_hyperG & glucoseC(m_hyper_point(i)-4-j)>m_hypoG)  ...
                        ...
                        &   abs (glucoseC(m_hyper_point(i)-1-j) - glucoseC(m_hyper_point(i)-2-j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)-1-j) - glucoseC(m_hyper_point(i)-3-j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)-2-j) - glucoseC(m_hyper_point(i)-3-j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)-1-j) - glucoseC(m_hyper_point(i)-4-j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)-2-j) - glucoseC(m_hyper_point(i)-4-j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)-3-j) - glucoseC(m_hyper_point(i)-4-j) ) <= 15  )
                        %start_m_hyper(i) = m_hyper_point(i-3);
                        start_m_hyper(i) = m_hyper_point(i)-4-j;
                        found = 1; 
                    end
                    j=j+1;
             end
            j=0;
            found = 0;
            while found == 0 & (m_hyper_point(i)+4+j <size(glucoseC,1))

                    if ( (glucoseC(m_hyper_point(i)+1+j)<m_hyperG & glucoseC(m_hyper_point(i)+1+j)>m_hypoG) ...
                        ...
                         ...
                        & (glucoseC(m_hyper_point(i)+2+j)<m_hyperG & glucoseC(m_hyper_point(i)+2+j)>m_hypoG) ...
                        ...
                        & (glucoseC(m_hyper_point(i)+3+j)<m_hyperG & glucoseC(m_hyper_point(i)+3+j)>m_hypoG) ...
                        & (glucoseC(m_hyper_point(i)+4+j)<m_hyperG & glucoseC(m_hyper_point(i)+4+j)>m_hypoG) ...
                        &   abs (glucoseC(m_hyper_point(i)+1+j) - glucoseC(m_hyper_point(i)+2+j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)+1+j) - glucoseC(m_hyper_point(i)+3+j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)+2+j) - glucoseC(m_hyper_point(i)+3+j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)+1+j) - glucoseC(m_hyper_point(i)+4+j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)+2+j) - glucoseC(m_hyper_point(i)+4+j) ) <= 15 ...
                        &   abs (glucoseC(m_hyper_point(i)+3+j) - glucoseC(m_hyper_point(i)+4+j) ) <= 15 )
                        %end_m_hyper(i) = m_hyper_point(i+3);
                        end_m_hyper(i) = m_hyper_point(i)+4+j;
                        found =1;
                    end
                  j=j+1;
            end
         end
     end
end
%%
if exist('start_m_hyper','var') & exist('end_m_hyper','var')
start_m_hyper = nonzeros(start_m_hyper);
end_m_hyper = nonzeros(end_m_hyper);

start_m_hyper = unique(start_m_hyper);
end_m_hyper = unique(end_m_hyper);

% figure();plot(universal_time,glucoseC)
% hold on
% plot(universal_time(start_m_hyper),glucoseC(start_m_hyper),'ro')
% plot(universal_time(end_m_hyper),glucoseC(end_m_hyper),'rx')


 if universal_time(1) > universal_time(2)
     ordered_start_m_hyper = (end_m_hyper);
     ordered_end_m_hyper = (start_m_hyper);
 else
    ordered_start_m_hyper = (start_m_hyper);
    ordered_end_m_hyper = (end_m_hyper);
 end

 if ordered_start_m_hyper(1)>=ordered_end_m_hyper(1)
        ordered_end_m_hyper = ordered_end_m_hyper(2:end);
 end

if ordered_start_m_hyper(end)>=ordered_end_m_hyper(end)
        ordered_start_m_hyper = ordered_start_m_hyper(1:end-1);
end

%find(m_hyper_point == ordered_start_m_hyper(1,1):ordered_end_m_hyper(1,1))
%check if the 1st end point should be there, if there are m hyper pooints
%between the 1st start and 1st end
if (size(ordered_start_m_hyper,1)>1 & size(ordered_end_m_hyper,1)>1)
if size(find(m_hyper_point >ordered_start_m_hyper(1,1) & m_hyper_point <ordered_end_m_hyper(1,1)),1) < 3
    ordered_end_m_hyper = ordered_end_m_hyper(2:end,1);
end
%check if the last start point should be there, if there are m hyper pooints
%between the last start and last end
if size(find(m_hyper_point >ordered_start_m_hyper(end,1) & m_hyper_point <ordered_end_m_hyper(end,1)),1) < 3
    ordered_start_m_hyper = ordered_start_m_hyper(1:end-1,1);
end
end
%ordered_end_m_hyper

%
    for i=1:size(ordered_start_m_hyper,1)
        m_hyper_event(i,1)=ordered_start_m_hyper(i);
        m_hyper_event(i,2)=ordered_end_m_hyper(i);
    end
m_hyperN = size(m_hyper_event,1);

else
    m_hyperN =0;
end

if m_hyperN >0
for i=1:min( [size(ordered_start_m_hyper,1) size(ordered_end_m_hyper,1) ] )
    ordered_start_end_m_hyper(i,1) =ordered_start_m_hyper(i);
    ordered_start_end_m_hyper(i,2) =ordered_end_m_hyper(i);
end



%find(m_hyper_point == ordered_start_m_hyper(1,1):ordered_end_m_hyper(1,1))
%check if the 1st end point should be there, if there are m hyper pooints
%between the 1st start and 1st end
% if (size(ordered_start_m_hyper,1) & size(ordered_end_m_hyper,1)) > 1
% if size(find(m_hyper_point >ordered_start_m_hyper(1,1) & m_hyper_point <ordered_end_m_hyper(1,1)),1) < 3
%     ordered_end_m_hyper = ordered_end_m_hyper(2:end,1);
% end
%check if the last start point should be there, if there are m hyper pooints
%between the last start and last end
% if size(find(m_hyper_point >ordered_start_m_hyper(end,1) & m_hyper_point <ordered_end_m_hyper(end,1)),1) < 3
%     ordered_start_m_hyper = ordered_start_m_hyper(1:end-1,1);
% end
% end
%ordered_end_m_hyper

%m_hyper_point


if size(ordered_start_m_hyper,1) > size(ordered_end_m_hyper,1)
    ordered_start_end_m_hyper_post_start(:,1)=ordered_start_m_hyper(1:end-1);
    ordered_start_end_m_hyper_post_start(:,2)=ordered_end_m_hyper(1:end);
end

if size(ordered_start_m_hyper,1) < size(ordered_end_m_hyper,1)
    ordered_start_end_m_hyper_post_start(:,1)=ordered_start_m_hyper(1:end);
    ordered_start_end_m_hyper_post_start(:,2)=ordered_end_m_hyper(1:end-1);
end

if size(ordered_start_m_hyper,1) == size(ordered_end_m_hyper,1)
    ordered_start_end_m_hyper_post_start(:,1)=ordered_start_m_hyper(1:end);
    ordered_start_end_m_hyper_post_start(:,2)=ordered_end_m_hyper(1:end);
end

length_m_hyper = ordered_start_end_m_hyper_post_start(:,2) - ordered_start_end_m_hyper_post_start(:,1);



 date_ordered_start_end_m_hyper_post_start(:,1) =  date(ordered_start_end_m_hyper_post_start(1:end,1));
 date_ordered_start_end_m_hyper_post_start(:,2) =  date(ordered_start_end_m_hyper_post_start(1:end,2));

 time_ordered_start_end_m_hyper_post_start(:,1) =  time(ordered_start_end_m_hyper_post_start(1:end,1))./scale_t;
 time_ordered_start_end_m_hyper_post_start(:,2) =  time(ordered_start_end_m_hyper_post_start(1:end,2))./scale_t;


% inputdate = '17037.46688591';
% jday = inputdate(1:5);
% date = datetime(jday, 'InputFormat', 'yyDDD') + days(mod(str2double(inputdate), 1))

%time_ordered_start_end_m_hyper_post_start(1,1)*scale_t   ./3600
%datestr(time_ordered_start_end_m_hyper_post_start(1,1),'hh:mm:ss')

%time_ordered_start_end_m_hyper_post_start(1,1)

for i=1:size(time_ordered_start_end_m_hyper_post_start,1)
    N_hours = fix( time_ordered_start_end_m_hyper_post_start(i,1) ./ (delta_t*3600) );
    N_minutes = fix(mod( time_ordered_start_end_m_hyper_post_start(i,1), (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds = mod( mod( time_ordered_start_end_m_hyper_post_start(i,1), (delta_t*3600) ) ./ (delta_t*60) , 1)*60;

    time_HMS_ordered_start_end_m_hyper_post_start(i,1:3)= [N_hours N_minutes N_seconds];

    N_hours = fix( time_ordered_start_end_m_hyper_post_start(i,2) ./ (delta_t*3600) );
    N_minutes = fix(mod( time_ordered_start_end_m_hyper_post_start(i,2), (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds = mod( mod( time_ordered_start_end_m_hyper_post_start(i,2), (delta_t*3600) ) ./ (delta_t*60) , 1)*60;

    time_HMS_ordered_start_end_m_hyper_post_start(i,4:6)= [N_hours N_minutes N_seconds];
end
end

   
%%
%datestr(0.25,'hh:mm:ss')

%datestr(time_ordered_start_end_m_hyper_post_start(1,1),'hh:mm')
%datestr(T,'hh:mm:ss PM')    


%%%% S hyper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (size(S_hyper_point,1)>0)
for i=1:size(S_hyper_point,1)-2 %was -3
    if S_hyper_point(i)>3
    
            j=0;
            found = 0;
            while found == 0 & (S_hyper_point(i)-4-j >2)
                 if ( (glucoseC(S_hyper_point(i)-1-j)<m_hyperG & glucoseC(S_hyper_point(i)-1-j)>m_hypoG) ...
                        ...
                     ...
                     & (glucoseC(S_hyper_point(i)-2-j)<m_hyperG & glucoseC(S_hyper_point(i)-2-j)>m_hypoG) ...
                        ...
                     & (glucoseC(S_hyper_point(i)-3-j)<m_hyperG & glucoseC(S_hyper_point(i)-3-j)>m_hypoG)  ...
                     & (glucoseC(S_hyper_point(i)-4-j)<m_hyperG & glucoseC(S_hyper_point(i)-4-j)>m_hypoG)  ...
                     &   abs (glucoseC(S_hyper_point(i)-1-j) - glucoseC(S_hyper_point(i)-2-j) ) <= 15 ...
                     &   abs (glucoseC(S_hyper_point(i)-1-j) - glucoseC(S_hyper_point(i)-3-j) ) <= 15 ...
                     &   abs (glucoseC(S_hyper_point(i)-2-j) - glucoseC(S_hyper_point(i)-3-j) ) <= 15 ...
                     &   abs (glucoseC(S_hyper_point(i)-1-j) - glucoseC(S_hyper_point(i)-4-j) ) <= 15 ...
                     &   abs (glucoseC(S_hyper_point(i)-2-j) - glucoseC(S_hyper_point(i)-4-j) ) <= 15 ...
                     &   abs (glucoseC(S_hyper_point(i)-3-j) - glucoseC(S_hyper_point(i)-4-j) ) <= 15 )
                    start_S_hyper(i) = S_hyper_point(i)-4-j;
                     found=1;
                end
            j=j+1;
            end
            j=0;
            found = 0;
            while found == 0 & (S_hyper_point(i)+4+j <size(glucoseC,1))
                if ( (glucoseC(S_hyper_point(i)+1+j)<m_hyperG & glucoseC(S_hyper_point(i)+1+j)>m_hypoG) ...
                     ...
                     ...
                     & (glucoseC(S_hyper_point(i)+2+j)<m_hyperG & glucoseC(S_hyper_point(i)+2+j)>m_hypoG) ...
                       ...
                    & (glucoseC(S_hyper_point(i)+3+j)<m_hyperG & glucoseC(S_hyper_point(i)+3+j)>m_hypoG)  ...
                     & (glucoseC(S_hyper_point(i)+4+j)<m_hyperG & glucoseC(S_hyper_point(i)+4+j)>m_hypoG)  ...
                 &   abs (glucoseC(S_hyper_point(i)+1+j) - glucoseC(S_hyper_point(i)+2+j) ) <= 15 ...
                 &   abs (glucoseC(S_hyper_point(i)+1+j) - glucoseC(S_hyper_point(i)+3+j) ) <= 15 ...
                    &   abs (glucoseC(S_hyper_point(i)+2+j) - glucoseC(S_hyper_point(i)+3+j) ) <= 15 ...  
                &   abs (glucoseC(S_hyper_point(i)+1+j) - glucoseC(S_hyper_point(i)+4+j) ) <= 15 ...
                 &   abs (glucoseC(S_hyper_point(i)+2+j) - glucoseC(S_hyper_point(i)+4+j) ) <= 15 ...
                    &   abs (glucoseC(S_hyper_point(i)+3+j) - glucoseC(S_hyper_point(i)+4+j) ) <= 15 )
                end_S_hyper(i) = S_hyper_point(i)+4+j;
                found=1;
                end
                j=j+1;
             end
    end
end
end
%%
if exist('start_S_hyper','var') & exist('end_S_hyper','var')

    start_S_hyper = nonzeros(start_S_hyper);
    end_S_hyper = nonzeros(end_S_hyper);

    start_S_hyper = unique(start_S_hyper);
    end_S_hyper = unique(end_S_hyper);
    if universal_time(1) > universal_time(2)
        ordered_start_S_hyper = (end_S_hyper);
        ordered_end_S_hyper = (start_S_hyper);
    else
        ordered_start_S_hyper = (start_S_hyper);
        ordered_end_S_hyper = (end_S_hyper);
    end

     if ordered_start_S_hyper(1)>=ordered_end_S_hyper(1)
        ordered_end_S_hyper = ordered_end_S_hyper(2:end);
    end
    if ordered_start_S_hyper(end)>=ordered_end_S_hyper(end)
        ordered_start_S_hyper = ordered_start_S_hyper(1:end-1);
    end

    %check if the 1st end point should be there, if there are m hyper pooints
    %between the 1st start and 1st end
    if (size(ordered_start_S_hyper,1)>1 & size(ordered_end_S_hyper,1)>1)
    if size(find(S_hyper_point >ordered_start_S_hyper(1,1) & S_hyper_point <ordered_end_S_hyper(1,1)),1) < 3
        ordered_end_S_hyper = ordered_end_S_hyper(2:end,1);
    end
    %check if the last start point should be there, if there are m hyper pooints
    %between the last start and last end
    if size(find(S_hyper_point >ordered_start_S_hyper(end,1) & S_hyper_point <ordered_end_S_hyper(end,1)),1) < 3
        ordered_start_S_hyper = ordered_start_S_hyper(1:end-1,1);
    end
    end

        for i=1:size(ordered_start_S_hyper,1)
            S_hyperG_event(i,1)=ordered_start_S_hyper(i);
            S_hyperG_event(i,2)=ordered_end_S_hyper(i);
        end

    S_hyperN = size(S_hyperG_event,1);

else
    S_hyperN =0;
end

if S_hyperN >0
for i=1:min( [size(ordered_start_S_hyper,1) size(ordered_end_S_hyper,1) ] )
    ordered_start_end_S_hyper(i,1) =ordered_start_S_hyper(i);
    ordered_start_end_S_hyper(i,2) =ordered_end_S_hyper(i);
end

% %check if the 1st end point should be there, if there are m hyper pooints
% %between the 1st start and 1st end
% if (size(ordered_start_S_hyper,1) & size(ordered_end_S_hyper,1)) > 1
% if size(find(S_hyper_point >ordered_start_S_hyper(1,1) & S_hyper_point <ordered_end_S_hyper(1,1)),1) < 3
%     ordered_end_S_hyper = ordered_end_S_hyper(2:end,1);
% end
% %check if the last start point should be there, if there are m hyper pooints
% %between the last start and last end
% if size(find(S_hyper_point >ordered_start_S_hyper(end,1) & S_hyper_point <ordered_end_S_hyper(end,1)),1) < 3
%     ordered_start_S_hyper = ordered_start_S_hyper(1:end-1,1);
% end
% end
if size(ordered_start_S_hyper,1) > size(ordered_end_S_hyper,1)
    ordered_start_end_S_hyper_post_start(:,1)=ordered_start_S_hyper(1:end-1);
    ordered_start_end_S_hyper_post_start(:,2)=ordered_end_S_hyper(1:end);
end

if size(ordered_start_S_hyper,1) < size(ordered_end_S_hyper,1)
    ordered_start_end_S_hyper_post_start(:,1)=ordered_start_S_hyper(1:end);
    ordered_start_end_S_hyper_post_start(:,2)=ordered_end_S_hyper(1:end-1);
end

if size(ordered_start_S_hyper,1) == size(ordered_end_S_hyper,1)
    ordered_start_end_S_hyper_post_start(:,1)=ordered_start_S_hyper(1:end);
    ordered_start_end_S_hyper_post_start(:,2)=ordered_end_S_hyper(1:end);
end

length_S_hyper = ordered_start_end_S_hyper_post_start(:,2) - ordered_start_end_S_hyper_post_start(:,1);

 date_ordered_start_end_S_hyper_post_start(:,1) =  date(ordered_start_end_S_hyper_post_start(1:end,1));
 date_ordered_start_end_S_hyper_post_start(:,2) =  date(ordered_start_end_S_hyper_post_start(1:end,2));

 time_ordered_start_end_S_hyper_post_start(:,1) =  time(ordered_start_end_S_hyper_post_start(1:end,1))./scale_t;
 time_ordered_start_end_S_hyper_post_start(:,2) =  time(ordered_start_end_S_hyper_post_start(1:end,2))./scale_t;

% inputdate = '17037.46688591';
% jday = inputdate(1:5);
% date = datetime(jday, 'InputFormat', 'yyDDD') + days(mod(str2double(inputdate), 1))

%time_ordered_start_end_m_hyper_post_start(1,1)*scale_t   ./3600
%datestr(time_ordered_start_end_m_hyper_post_start(1,1),'hh:mm:ss')

%time_ordered_start_end_m_hyper_post_start(1,1)

for i=1:size(time_ordered_start_end_S_hyper_post_start,1)
    N_hours = fix( time_ordered_start_end_S_hyper_post_start(i,1) ./ (delta_t*3600) );
    N_minutes = fix(mod( time_ordered_start_end_S_hyper_post_start(i,1), (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds = mod( mod( time_ordered_start_end_S_hyper_post_start(i,1), (delta_t*3600) ) ./ (delta_t*60) , 1)*60;

    time_HMS_ordered_start_end_S_hyper_post_start(i,1:3)= [N_hours N_minutes N_seconds];

    N_hours = fix( time_ordered_start_end_S_hyper_post_start(i,2) ./ (delta_t*3600) );
    N_minutes = fix(mod( time_ordered_start_end_S_hyper_post_start(i,2), (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds = mod( mod( time_ordered_start_end_S_hyper_post_start(i,2), (delta_t*3600) ) ./ (delta_t*60) , 1)*60;

    time_HMS_ordered_start_end_S_hyper_post_start(i,4:6)= [N_hours N_minutes N_seconds];
end
end

%%%% m hypo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (size(m_hypo_point,1) >0)
% for i=1:size(m_hypo_point,1)-3
%     if m_hypo_point(i)>3
%     j=0;
%     found=0;
%         while found ==0 & (m_hypo_point(i)-3-j>2)
%             if ( (glucoseC(m_hypo_point(i)-1)<m_hyperG & glucoseC(m_hypo_point(i)-1)>m_hypoG) ...
%                     ...
%                     ...
%                     & (glucoseC(m_hypo_point(i)-2)<m_hyperG & glucoseC(m_hypo_point(i)-2)>m_hypoG) ...
%                     ...
%                     & (glucoseC(m_hypo_point(i)-3)<m_hyperG & glucoseC(m_hypo_point(i)-3)>m_hypoG)  ...
%                     &   abs (glucoseC(m_hypo_point(i)-1) - glucoseC(m_hypo_point(i)-2) ) <= 15 ...
%                     &   abs (glucoseC(m_hypo_point(i)-1) - glucoseC(m_hypo_point(i)-3) ) <= 15 ...
%                     &   abs (glucoseC(m_hypo_point(i)-2) - glucoseC(m_hypo_point(i)-3) ) <= 15 ) 
%                     start_m_hypo(i) = m_hypo_point(i) -3;
%                     found =1;
%             end
%             j=j+1;
%         end
%         j=0;
%         found =0;
%         while found ==0 & (m_hypo_point(i)+3+j< size(glucoseC,1))
%             if ( (glucoseC(m_hypo_point(i)+1)<m_hyperG & glucoseC(m_hypo_point(i)+1)>m_hypoG) ...
%                     ...
%                     ...
%                     & (glucoseC(m_hypo_point(i)+2)<m_hyperG & glucoseC(m_hypo_point(i)+2)>m_hypoG) ...
%                     ...
%                     & (glucoseC(m_hypo_point(i)+3)<m_hyperG & glucoseC(m_hypo_point(i)+3)>m_hypoG) ... 
%                     &   abs (glucoseC(m_hypo_point(i)+1) - glucoseC(m_hypo_point(i)+2) ) <= 15 ...
%                     &   abs (glucoseC(m_hypo_point(i)+1) - glucoseC(m_hypo_point(i)+3) ) <= 15 ...
%                     &   abs (glucoseC(m_hypo_point(i)+2) - glucoseC(m_hypo_point(i)+3) ) <= 15 ) 
%                 end_m_hypo(i) = m_hypo_point(i)+3;
%                 found=1;
%             end
%             j=j+1;
%         end
%     end
% 
% end
% end
% test j
if (size(m_hypo_point,1) >0)
for i=1:size(m_hypo_point,1)-2 %was -3
    if m_hypo_point(i)>3
    j=0;
    found=0;
        while found ==0 & (m_hypo_point(i)-4-j>2)
            if ( (glucoseC(m_hypo_point(i)-1-j)<m_hyperG & glucoseC(m_hypo_point(i)-1-j)>m_hypoG) ...
                    ...
                    ...
                    & (glucoseC(m_hypo_point(i)-2-j)<m_hyperG & glucoseC(m_hypo_point(i)-2-j)>m_hypoG) ...
                    ...
                    & (glucoseC(m_hypo_point(i)-3-j)<m_hyperG & glucoseC(m_hypo_point(i)-3-j)>m_hypoG)  ...
                    & (glucoseC(m_hypo_point(i)-4-j)<m_hyperG & glucoseC(m_hypo_point(i)-4-j)>m_hypoG)  ...
                    &   abs (glucoseC(m_hypo_point(i)-1-j) - glucoseC(m_hypo_point(i)-2-j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)-1-j) - glucoseC(m_hypo_point(i)-3-j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)-2-j) - glucoseC(m_hypo_point(i)-3-j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)-1-j) - glucoseC(m_hypo_point(i)-4-j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)-2-j) - glucoseC(m_hypo_point(i)-4-j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)-3-j) - glucoseC(m_hypo_point(i)-4-j) ) <= 15 ) 
                    start_m_hypo(i) = m_hypo_point(i) -4-j;
                    found =1;
            end
            j=j+1;
        end
        j=0;
        found =0;
        while found ==0 & (m_hypo_point(i)+4+j< size(glucoseC,1))
            if ( (glucoseC(m_hypo_point(i)+1+j)<m_hyperG & glucoseC(m_hypo_point(i)+1+j)>m_hypoG) ...
                    ...
                    ...
                    & (glucoseC(m_hypo_point(i)+2+j)<m_hyperG & glucoseC(m_hypo_point(i)+2+j)>m_hypoG) ...
                    ...
                    & (glucoseC(m_hypo_point(i)+3+j)<m_hyperG & glucoseC(m_hypo_point(i)+3+j)>m_hypoG) ...
                    & (glucoseC(m_hypo_point(i)+4+j)<m_hyperG & glucoseC(m_hypo_point(i)+4+j)>m_hypoG) ...
                    &   abs (glucoseC(m_hypo_point(i)+1+j) - glucoseC(m_hypo_point(i)+2+j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)+1+j) - glucoseC(m_hypo_point(i)+3+j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)+2+j) - glucoseC(m_hypo_point(i)+3+j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)+1+j) - glucoseC(m_hypo_point(i)+4+j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)+2+j) - glucoseC(m_hypo_point(i)+4+j) ) <= 15 ...
                    &   abs (glucoseC(m_hypo_point(i)+3+j) - glucoseC(m_hypo_point(i)+4+j) ) <= 15 ) 
                end_m_hypo(i) = m_hypo_point(i)+4+j;
                found=1;
            end
            j=j+1;
        end
    end
    
end
end

%%
if exist('start_m_hypo','var') & exist('end_m_hypo','var')
start_m_hypo = nonzeros(start_m_hypo);
end_m_hypo = nonzeros(end_m_hypo);

start_m_hypo = unique(start_m_hypo);
end_m_hypo = unique(end_m_hypo);

 if universal_time(1) > universal_time(2)
     ordered_start_m_hypo = (end_m_hypo);
     ordered_end_m_hypo = (start_m_hypo);
 else
    ordered_start_m_hypo = (start_m_hypo);
    ordered_end_m_hypo = (end_m_hypo);
 end
  if ordered_start_m_hypo(1)>=ordered_end_m_hypo(1)
        ordered_end_m_hypo = ordered_end_m_hypo(2:end);
  end

  if ordered_start_m_hypo(end)>=ordered_end_m_hypo(end)
        ordered_start_m_hypo = ordered_start_m_hypo(1:end-1);
  end

   %check if the 1st end point should be there, if there are m hyper pooints
    %between the 1st start and 1st end
    if (size(ordered_start_m_hypo,1)>1 & size(ordered_end_m_hypo,1)>1)
    if size(find(m_hypo_point >ordered_start_m_hypo(1,1) & m_hypo_point <ordered_end_m_hypo(1,1)),1) < 3
        ordered_end_m_hypo = ordered_end_m_hypo(2:end,1);
    end
    %check if the last start point should be there, if there are m hyper pooints
    %between the last start and last end
    if size(find(m_hypo_point >ordered_start_m_hypo(end,1) & m_hypo_point <ordered_end_m_hypo(end,1)),1) < 3
        ordered_start_m_hypo = ordered_start_m_hypo(1:end-1,1);
    end
    end

    for i=1:size(ordered_start_m_hypo,1)
        m_hypo_event(i,1)=ordered_start_m_hypo(i);
        m_hypo_event(i,2)=ordered_end_m_hypo(i);
    end
m_hypoN = size(m_hypo_event,1);

else
    m_hypoN =0;
end

%time stuff

if m_hypoN >0 
for i=1:min( [size(ordered_start_m_hypo,1) size(ordered_end_m_hypo,1) ] )
    ordered_start_end_m_hypo(i,1) =ordered_start_m_hypo(i);
    ordered_start_end_m_hypo(i,2) =ordered_end_m_hypo(i);
end

%check if the 1st end point should be there, if there are m hyper pooints
%between the 1st start and 1st end
% if (size(ordered_start_m_hypo,1) & size(ordered_end_m_hypo,1)) > 1
% if size(find(m_hypo_point >ordered_start_m_hypo(1,1) & m_hypo_point <ordered_end_m_hypo(1,1)),1) < 3
%     ordered_end_m_hypo = ordered_end_m_hypo(2:end,1);
% end
% %check if the last start point should be there, if there are m hyper pooints
% %between the last start and last end
% if size(find(m_hypo_point >ordered_start_m_hypo(end,1) & m_hypo_point <ordered_end_m_hypo(end,1)),1) < 3
%     ordered_start_m_hypo = ordered_start_m_hypo(1:end-1,1);
% end
% end
if size(ordered_start_m_hypo,1) > size(ordered_end_m_hypo,1)
    ordered_start_end_m_hypo_post_start(:,1)=ordered_start_m_hypo(1:end-1);
    ordered_start_end_m_hypo_post_start(:,2)=ordered_end_m_hypo(1:end);
end

if size(ordered_start_m_hypo,1) < size(ordered_end_m_hypo,1)
    ordered_start_end_m_hypo_post_start(:,1)=ordered_start_m_hypo(1:end);
    ordered_start_end_m_hypo_post_start(:,2)=ordered_end_m_hypo(1:end-1);
end

if size(ordered_start_m_hypo,1) == size(ordered_end_m_hypo,1)
    ordered_start_end_m_hypo_post_start(:,1)=ordered_start_m_hypo(1:end);
    ordered_start_end_m_hypo_post_start(:,2)=ordered_end_m_hypo(1:end);
end
length_m_hypo = ordered_start_end_m_hypo_post_start(:,2) - ordered_start_end_m_hypo_post_start(:,1);

 date_ordered_start_end_m_hypo_post_start(:,1) =  date(ordered_start_end_m_hypo_post_start(1:end,1));
 date_ordered_start_end_m_hypo_post_start(:,2) =  date(ordered_start_end_m_hypo_post_start(1:end,2));

 time_ordered_start_end_m_hypo_post_start(:,1) =  time(ordered_start_end_m_hypo_post_start(1:end,1))./scale_t;
 time_ordered_start_end_m_hypo_post_start(:,2) =  time(ordered_start_end_m_hypo_post_start(1:end,2))./scale_t;

for i=1:size(time_ordered_start_end_m_hypo_post_start,1)
    N_hours = fix( time_ordered_start_end_m_hypo_post_start(i,1) ./ (delta_t*3600) );
    N_minutes = fix(mod( time_ordered_start_end_m_hypo_post_start(i,1), (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds = mod( mod( time_ordered_start_end_m_hypo_post_start(i,1), (delta_t*3600) ) ./ (delta_t*60) , 1)*60;

    time_HMS_ordered_start_end_m_hypo_post_start(i,1:3)= [N_hours N_minutes N_seconds];

    N_hours = fix( time_ordered_start_end_m_hypo_post_start(i,2) ./ (delta_t*3600) );
    N_minutes = fix(mod( time_ordered_start_end_m_hypo_post_start(i,2), (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds = mod( mod( time_ordered_start_end_m_hypo_post_start(i,2), (delta_t*3600) ) ./ (delta_t*60) , 1)*60;

    time_HMS_ordered_start_end_m_hypo_post_start(i,4:6)= [N_hours N_minutes N_seconds];
end
end
%stop here Guy!!
%%
%%%% S hypo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (size(S_hypo_point,1) >0)
for i=1:size(S_hypo_point,1)-2 %was -3
    if S_hypo_point(i)>3
            j=0;
            found = 0;
            while found == 0 & (S_hypo_point(i)-4-j >2)
                     if ( (glucoseC(S_hypo_point(i)-1-j)<m_hyperG & glucoseC(S_hypo_point(i)-1-j)>m_hypoG) ...
                        ...
                        ...
                        & (glucoseC(S_hypo_point(i)-2-j)<m_hyperG & glucoseC(S_hypo_point(i)-2-j)>m_hypoG) ...
                        ...
                        & (glucoseC(S_hypo_point(i)-3-j)<m_hyperG & glucoseC(S_hypo_point(i)-3-j)>m_hypoG)  ...
                        & (glucoseC(S_hypo_point(i)-4-j)<m_hyperG & glucoseC(S_hypo_point(i)-4-j)>m_hypoG)  ...
                        &   abs (glucoseC(S_hypo_point(i)-1-j) - glucoseC(S_hypo_point(i)-2-j) ) <= 15 ...
                        &   abs (glucoseC(S_hypo_point(i)-1-j) - glucoseC(S_hypo_point(i)-3-j) ) <= 15 ...
                        &   abs (glucoseC(S_hypo_point(i)-2-j) - glucoseC(S_hypo_point(i)-3-j) ) <= 15 ...
                        &   abs (glucoseC(S_hypo_point(i)-1-j) - glucoseC(S_hypo_point(i)-4-j) ) <= 15 ...
                        &   abs (glucoseC(S_hypo_point(i)-2-j) - glucoseC(S_hypo_point(i)-4-j) ) <= 15 ...
                        &   abs (glucoseC(S_hypo_point(i)-3-j) - glucoseC(S_hypo_point(i)-4-j) ) <= 15 )
                        start_S_hypo(i) = S_hypo_point(i)-4-j;
                        found = 1;
                      end
                j=j+1;
            end
        j=0;
        found =0;
        while found ==0 & (S_hypo_point(i)+4+j< size(glucoseC,1))
            if ( (glucoseC(S_hypo_point(i)+1+j)<m_hyperG & glucoseC(S_hypo_point(i)+1+j)>m_hypoG) ...
                 ...
                    ...
                    & (glucoseC(S_hypo_point(i)+2+j)<m_hyperG & glucoseC(S_hypo_point(i)+2+j)>m_hypoG) ...
                    ...
                    & (glucoseC(S_hypo_point(i)+3+j)<m_hyperG & glucoseC(S_hypo_point(i)+3+j)>m_hypoG)  ...
                    & (glucoseC(S_hypo_point(i)+4+j)<m_hyperG & glucoseC(S_hypo_point(i)+4+j)>m_hypoG)  ...
                    &   abs (glucoseC(S_hypo_point(i)+1+j) - glucoseC(S_hypo_point(i)+2+j) ) <= 15 ...
                    &   abs (glucoseC(S_hypo_point(i)+1+j) - glucoseC(S_hypo_point(i)+3+j) ) <= 15 ...
                    &   abs (glucoseC(S_hypo_point(i)+2+j) - glucoseC(S_hypo_point(i)+3+j) ) <= 15 ...
                    &   abs (glucoseC(S_hypo_point(i)+1+j) - glucoseC(S_hypo_point(i)+4+j) ) <= 15 ...
                    &   abs (glucoseC(S_hypo_point(i)+2+j) - glucoseC(S_hypo_point(i)+4+j) ) <= 15 ...
                    &   abs (glucoseC(S_hypo_point(i)+3+j) - glucoseC(S_hypo_point(i)+4+j) ) <= 15 )
                end_S_hypo(i) = S_hypo_point(i)+4+j;
                found = 1;
            end
            j=j+1;
        end
    end
end
end

if exist('start_S_hypo','var') & exist('end_S_hypo','var')
start_S_hypo = nonzeros(start_S_hypo);
end_S_hypo = nonzeros(end_S_hypo);

start_S_hypo = unique(start_S_hypo);
end_S_hypo = unique(end_S_hypo);

 if universal_time(1) > universal_time(2)
     ordered_start_S_hypo = (end_S_hypo);
     ordered_end_S_hypo = (start_S_hypo);
 else
    ordered_start_S_hypo = (start_S_hypo);
    ordered_end_S_hypo = (end_S_hypo);
 end
 if ordered_start_S_hypo(1)>=ordered_end_S_hypo(1)
        ordered_end_S_hypo = ordered_end_S_hypo(2:end);
 end

  if ordered_start_S_hypo(end)>=ordered_end_S_hypo(end)
        ordered_start_S_hypo = ordered_start_S_hypo(1:end-1);
 end
    %check if the 1st end point should be there, if there are m hyper pooints
    %between the 1st start and 1st end
    if (size(ordered_start_S_hypo,1)>1 & size(ordered_end_S_hypo,1)>1)
    if size(find(S_hypo_point >ordered_start_S_hypo(1,1) & S_hypo_point <ordered_end_S_hypo(1,1)),1) < 3
        ordered_end_S_hypo = ordered_end_S_hypo(2:end,1);
    end
    %check if the last start point should be there, if there are m hyper pooints
    %between the last start and last end
    if size(find(S_hypo_point >ordered_start_S_hypo(end,1) & S_hypo_point <ordered_end_S_hypo(end,1)),1) < 3
        ordered_start_S_hypo = ordered_start_S_hypo(1:end-1,1);
    end
    end

    for i=1:size(ordered_start_S_hypo,1)
        S_hypoG_event(i,1)=ordered_start_S_hypo(i);
        S_hypoG_event(i,2)=ordered_end_S_hypo(i);
    end
S_hypoN = size(S_hypoG_event,1);


else
    S_hypoN =0;
end

%time stuff

if S_hypoN >0 
for i=1:min( [size(ordered_start_S_hypo,1) size(ordered_end_S_hypo,1) ] )
    ordered_start_end_S_hypoG(i,1) =ordered_start_S_hypo(i);
    ordered_start_end_S_hypoG(i,2) =ordered_end_S_hypo(i);
end

%check if the 1st end point should be there, if there are m hyper pooints
%between the 1st start and 1st end
% if (size(ordered_start_S_hypo,1) & size(ordered_end_S_hypo,1)) > 1
% if size(find(S_hypo_point >ordered_start_S_hypo(1,1) & S_hypo_point <ordered_end_S_hypo(1,1)),1) < 3
%     ordered_end_S_hypo = ordered_end_S_hypo(2:end,1);
% end
% %check if the last start point should be there, if there are m hyper pooints
% %between the last start and last end
% if size(find(S_hypo_point >ordered_start_S_hypo(end,1) & S_hypo_point <ordered_end_S_hypo(end,1)),1) < 3
%     ordered_start_S_hypo = ordered_start_S_hypo(1:end-1,1);
% end
% end
if size(ordered_start_S_hypo,1) > size(ordered_end_S_hypo,1)
    ordered_start_end_S_hypoG_post_start(:,1)=ordered_start_S_hypo(1:end-1);
    ordered_start_end_S_hypoG_post_start(:,2)=ordered_end_S_hypo(1:end);
end

if size(ordered_start_S_hypo,1) < size(ordered_end_S_hypo,1)
    ordered_start_end_S_hypoG_post_start(:,1)=ordered_start_S_hypo(1:end);
    ordered_start_end_S_hypoG_post_start(:,2)=ordered_end_S_hypo(1:end-1);
end

if size(ordered_start_S_hypo,1) == size(ordered_end_S_hypo,1)
    ordered_start_end_S_hypoG_post_start(:,1)=ordered_start_S_hypo(1:end);
    ordered_start_end_S_hypoG_post_start(:,2)=ordered_end_S_hypo(1:end);
end
length_S_hypoG = ordered_start_end_S_hypoG_post_start(:,2) - ordered_start_end_S_hypoG_post_start(:,1);

 date_ordered_start_end_S_hypoG_post_start(:,1) =  date(ordered_start_end_S_hypoG_post_start(1:end,1));
 date_ordered_start_end_S_hypoG_post_start(:,2) =  date(ordered_start_end_S_hypoG_post_start(1:end,2));

 time_ordered_start_end_S_hypoG_post_start(:,1) =  time(ordered_start_end_S_hypoG_post_start(1:end,1))./scale_t;
 time_ordered_start_end_S_hypoG_post_start(:,2) =  time(ordered_start_end_S_hypoG_post_start(1:end,2))./scale_t;

for i=1:size(time_ordered_start_end_S_hypoG_post_start,1)
    N_hours = fix( time_ordered_start_end_S_hypoG_post_start(i,1) ./ (delta_t*3600) );
    N_minutes = fix(mod( time_ordered_start_end_S_hypoG_post_start(i,1), (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds = mod( mod( time_ordered_start_end_S_hypoG_post_start(i,1), (delta_t*3600) ) ./ (delta_t*60) , 1)*60;

    time_HMS_ordered_start_end_S_hypoG_post_start(i,1:3)= [N_hours N_minutes N_seconds];

    N_hours = fix( time_ordered_start_end_S_hypoG_post_start(i,2) ./ (delta_t*3600) );
    N_minutes = fix(mod( time_ordered_start_end_S_hypoG_post_start(i,2), (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds = mod( mod( time_ordered_start_end_S_hypoG_post_start(i,2), (delta_t*3600) ) ./ (delta_t*60) , 1)*60;

    time_HMS_ordered_start_end_S_hypoG_post_start(i,4:6)= [N_hours N_minutes N_seconds];
end
end

%%%%%%%%%%%%%%%%%% Check if any of the mild Hyper/Hypo events cross into
%%%%%%%%%%%%%%%%%% Severe
%%

%ordered_start_end_m_hyper_post_start;

if exist('ordered_start_end_m_hyper_post_start','var')
    for i=1:size(ordered_start_end_m_hyper_post_start,1)
        for j=ordered_start_end_m_hyper_post_start(i,1):ordered_start_end_m_hyper_post_start(i,2)-2
            if glucoseC(j,1)>S_hyperG & glucoseC(j+1,1)>S_hyperG & glucoseC(j+2,1)>S_hyperG
                found_S_hyperG(i) =i;
            end
        end
    end
% need to append on to existing S hyper params if they exist
    if exist('found_S_hyperG','var')
        found_S_hyperG = nonzeros(found_S_hyperG);
        length_S_hyper = length_m_hyper(found_S_hyperG,:);
        S_hyperN = size(found_S_hyperG,1);
        time_HMS_ordered_start_end_S_hyper_post_start = time_HMS_ordered_start_end_m_hyper_post_start(found_S_hyperG,:);
        ordered_start_S_hyper = ordered_start_m_hyper(found_S_hyperG,:);
        ordered_end_S_hyper = ordered_end_m_hyper(found_S_hyperG,:);
        ordered_start_end_S_hyper_post_start = ordered_start_end_m_hyper_post_start(found_S_hyperG,:);

        found_m_hyper = (1:m_hyperN)';

        unique_m_hyper = setdiff(found_m_hyper, found_S_hyperG);
        time_HMS_ordered_start_end_m_hyper_post_start = time_HMS_ordered_start_end_m_hyper_post_start(unique_m_hyper,:);
        length_m_hyper = length_m_hyper(unique_m_hyper,:);
        m_hyperN = size(unique_m_hyper,1);

        ordered_start_m_hyper = ordered_start_m_hyper(unique_m_hyper,:);
        ordered_end_m_hyper = ordered_end_m_hyper(unique_m_hyper,:);
        ordered_start_end_m_hyper_post_start = ordered_start_end_m_hyper_post_start(unique_m_hyper,:);
    end
end
%%
if exist('ordered_start_end_m_hypo_post_start','var')
    for i=1:size(ordered_start_end_m_hypo_post_start,1)
        for j=ordered_start_end_m_hypo_post_start(i,1):ordered_start_end_m_hypo_post_start(i,2)-2
            if glucoseC(j,1)<S_hypoG & glucoseC(j+1,1)<S_hypoG & glucoseC(j+2,1)<S_hypoG
                found_S_hypoG(i) =i;
            end
        end
    end
% need to append on to existing S hypo params if they exist
    if exist('found_S_hypoG','var')
        found_S_hypoG = nonzeros(found_S_hypoG);
        length_S_hypo = length_m_hypo(found_S_hypoG,:);
        S_hypoN = size(found_S_hypoG,1);
        time_HMS_ordered_start_end_S_hypo_post_start = time_HMS_ordered_start_end_m_hypo_post_start(found_S_hypoG,:);
        ordered_start_S_hypo = ordered_start_m_hypo(found_S_hypoG,:);
        ordered_end_S_hypo = ordered_end_m_hypo(found_S_hypoG,:);
        ordered_start_end_S_hypo_post_start = ordered_start_end_m_hypo_post_start(found_S_hypoG,:);

        found_m_hypo = (1:m_hypoN)';

        unique_m_hypo = setdiff(found_m_hypo, found_S_hypoG);
        time_HMS_ordered_start_end_m_hypo_post_start = time_HMS_ordered_start_end_m_hypo_post_start(unique_m_hypo,:);
        %time_HMS_ordered_start_end_m_hypo_post_start = time_HMS_ordered_start_end_m_hypo_post_start(unique_m_hypo,:);
        length_m_hypo = length_m_hypo(unique_m_hypo,:);
        m_hypoN = size(unique_m_hypo,1);

        ordered_start_m_hypo = ordered_start_m_hypo(unique_m_hypo,:);
        ordered_end_m_hypo = ordered_end_m_hypo(unique_m_hypo,:);
        ordered_start_end_m_hypo_post_start = ordered_start_end_m_hypo_post_start(unique_m_hypo,:);

        


    end
end



%%
% check if any events have jumps in them
% 
% jump
% length_jump
% 
% ordered_start_end_S_hyper_post_start
% ordered_start_end_m_hyper_post_start
% ordered_start_end_S_hypo_post_start
% ordered_start_end_m_hypo_post_start
if size(delta_t_gap_idx,1)>0
    jump = delta_t_gap_start;
end
if exist('jump','var')
    if exist('ordered_start_end_m_hyper_post_start','var')
        for i=1:size(ordered_start_end_m_hyper_post_start,1)
           N_m_hyper_jumps_length(i,1) = size (find(jump > ordered_start_end_m_hyper_post_start(i,1) & jump < ordered_start_end_m_hyper_post_start(i,2) ) , 1);
           N_m_hyper_jumps_length(i,2) =  sum(length_jump(find(jump > ordered_start_end_m_hyper_post_start(i,1) & jump < ordered_start_end_m_hyper_post_start(i,2) )));
        end
    end
    
    if exist('ordered_start_end_S_hyper_post_start','var')
        for i=1:size(ordered_start_end_S_hyper_post_start,1)
           N_S_hyper_jumps_length(i,1) = size (find(jump > ordered_start_end_S_hyper_post_start(i,1) & jump < ordered_start_end_S_hyper_post_start(i,2) ) , 1);
           N_S_hyper_jumps_length(i,2) =  sum(length_jump(find(jump > ordered_start_end_S_hyper_post_start(i,1) & jump < ordered_start_end_S_hyper_post_start(i,2) )));
        end
    end

    if exist('ordered_start_end_m_hypo_post_start','var')
        for i=1:size(ordered_start_end_m_hypo_post_start,1)
           N_m_hypo_jumps_length(i,1) = size (find(jump > ordered_start_end_m_hypo_post_start(i,1) & jump < ordered_start_end_m_hypo_post_start(i,2) ) , 1);
           N_m_hypo_jumps_length(i,2) =  sum(length_jump(find(jump > ordered_start_end_m_hypo_post_start(i,1) & jump < ordered_start_end_m_hypo_post_start(i,2) )));
        end
    end
    
    if exist('ordered_start_end_S_hypo_post_start','var')
        for i=1:size(ordered_start_end_S_hypo_post_start,1)
           N_S_hypo_jumps_length(i,1) = size (find(jump > ordered_start_end_S_hypo_post_start(i,1) & jump < ordered_start_end_S_hypo_post_start(i,2) ) , 1);
           N_S_hypo_jumps_length(i,2) =  sum(length_jump(find(jump > ordered_start_end_S_hypo_post_start(i,1) & jump < ordered_start_end_S_hypo_post_start(i,2) )));
        end
    end
end

    


% S_hyperN =1;
% S_hypoN =1;
% 
% m_hyperN =1;
% m_hypoN =1;
% 
% if size(S_hyper_point,1) == 0
%     S_hyperN =0;
% else
%     for i=2:size(S_hyper_point,1)
% 
%         if S_hyper_point(i)  ~= S_hyper_point(i-1)+1
%             S_hyperN = S_hyperN + 1;
%         end
%     end
% end
% 
% if size(S_hypo_point,1) == 0
%     S_hypoN =0;
% else
% 
%     for i=2:size(S_hypo_point,1)
% 
%         if S_hypo_point(i)  ~= S_hypo_point(i-1)+1
%             S_hypoN = S_hypoN + 1;
%         end
%     end
% end
% 
% if size(m_hyper_point,1) == 0
%     m_hyperN =0;
% else
%     for i=2:size(m_hyper_point,1)
% 
%         if m_hyper_point(i)  ~= m_hyper_point(i-1)+1
%             m_hyperN = m_hyperN + 1;
%         end
%     end
% end
% 
% if size(m_hypo_point,1) == 0
%     m_hypoN =0;
% else
% 
%     for i=2:size(m_hypo_point,1)
% 
%         if m_hypo_point(i)  ~= m_hypo_point(i-1)+1
%             m_hypoN = m_hypoN + 1;
%         end
%     end
% end

% 'S Hyper N'
% S_hyperN
% 
% 'm Hyper N'
% m_hyperN
% 
% 'S Hypo N'
% S_hypoN
% 
% 'm Hypo N'
% m_hypoN

%parameters we want to return
hyperSm_hypoSm_N = [S_hyperN m_hyperN S_hypoN m_hypoN];

%glucose_events_info.N_events = hyperSm_hypoSm_N;

glucose_events_info.N_events.S_hyper = S_hyperN;
glucose_events_info.N_events.m_hyper = m_hyperN;
glucose_events_info.N_events.S_hypo = S_hypoN;
glucose_events_info.N_events.m_hypo = m_hypoN;


if not(exist('length_S_hyper','var'))
    length_S_hyper = 0;
end
if not(exist('length_m_hyper','var'))
    length_m_hyper = 0;
end
if not(exist('length_S_hypo','var'))
    length_S_hypo = 0;
end
if not(exist('length_m_hypo','var'))
    length_m_hypo = 0;
end

glucose_events_info.events_length.S_hyper = length_S_hyper;
glucose_events_info.events_length.m_hyper = length_m_hyper;
glucose_events_info.events_length.S_hypo = length_S_hypo;
glucose_events_info.events_length.m_hypo = length_m_hypo;

% length_S_hyper
% length_m_hyper
% length_S_hypo
% length_m_hypo
% 
% ordered_start_end_S_hyper_post_start
% ordered_start_end_m_hyper_post_start
% ordered_start_end_S_hypo_post_start
% ordered_start_end_m_hypo_post_start
if not(exist('ordered_start_end_S_hyper_post_start','var'))
    ordered_start_end_S_hyper_post_start = [0 0];
end
if not(exist('ordered_start_end_m_hyper_post_start','var'))
    ordered_start_end_m_hyper_post_start = [0 0];
end
if not(exist('ordered_start_end_S_hypo_post_start','var'))
    ordered_start_end_S_hypo_post_start = [0 0]; 
end
if not(exist('ordered_start_end_m_hypo_post_start','var'))
    ordered_start_end_m_hypo_post_start = [0 0];
end
glucose_events_info.events_start_end.S_hyper = ordered_start_end_S_hyper_post_start;
glucose_events_info.events_start_end.m_hyper = ordered_start_end_m_hyper_post_start;
glucose_events_info.events_start_end.S_hypo = ordered_start_end_S_hypo_post_start;
glucose_events_info.events_start_end.m_hypo = ordered_start_end_m_hypo_post_start;

% time_HMS_ordered_start_end_S_hyper_post_start
% time_HMS_ordered_start_end_m_hyper_post_start
% time_HMS_ordered_start_end_S_hypo_post_start
% time_HMS_ordered_start_end_m_hypo_post_start

if not(exist('time_HMS_ordered_start_end_S_hyper_post_start','var'))
    time_HMS_ordered_start_end_S_hyper_post_start = [0 0 0 0];
end
if not(exist('time_HMS_ordered_start_end_m_hyper_post_start','var'))
    time_HMS_ordered_start_end_m_hyper_post_start = [0 0 0 0];
end
if not(exist('time_HMS_ordered_start_end_S_hypo_post_start','var'))
    time_HMS_ordered_start_end_S_hypo_post_start = [0 0 0 0]; 
end
if not(exist('time_HMS_ordered_start_end_m_hypo_post_start','var'))
    time_HMS_ordered_start_end_m_hypo_post_start = [0 0 0 0];
end

glucose_events_info.time_of_events.S_hyper = time_HMS_ordered_start_end_S_hyper_post_start;
glucose_events_info.time_of_events.m_hyper = time_HMS_ordered_start_end_m_hyper_post_start;
glucose_events_info.time_of_events.S_hypo = time_HMS_ordered_start_end_S_hypo_post_start;
glucose_events_info.time_of_events.m_hypo = time_HMS_ordered_start_end_m_hypo_post_start;


% N_S_hyper_jumps_length
% N_m_hyper_jumps_length
% N_S_hypo_jumps_length
% N_m_hypo_jumps_length

if not(exist('N_S_hyper_jumps_length','var'))
    N_S_hyper_jumps_length = [0 0];
end
if not(exist('N_m_hyper_jumps_length','var'))
    N_m_hyper_jumps_length = [0 0];
end
if not(exist('N_S_hypo_jumps_length','var'))
    N_S_hypo_jumps_length = [0 0]; 
end
if not(exist('N_m_hypo_jumps_length','var'))
    N_m_hypo_jumps_length = [0 0];
end

glucose_events_info.glucose_data_gaps.S_hyper = N_S_hyper_jumps_length;
glucose_events_info.glucose_data_gaps.m_hyper = N_m_hyper_jumps_length;
glucose_events_info.glucose_data_gaps.S_hypo = N_S_hypo_jumps_length;
glucose_events_info.glucose_data_gaps.m_hypo = N_m_hypo_jumps_length;
%%
%remove events if they have jumps in them of length more than 8 samples. 

%glucose_events_info.events_start_end.S_hyper
if exist('N_m_hyper_jumps_length','var')
    for i = 1:size(N_m_hyper_jumps_length,1)
        if N_m_hyper_jumps_length(i,2) > 8
         m_hyper_jumps_idx(i,1)=i;
        end
    end
    if exist('m_hyper_jumps_idx','var')
    m_hyper_jumps_idx = nonzeros(m_hyper_jumps_idx);
    glucose_events_info.N_events.m_hyper = glucose_events_info.N_events.m_hyper - size(m_hyper_jumps_idx,1);
    glucose_events_info.events_length.m_hyper(m_hyper_jumps_idx)=0;  
    glucose_events_info.events_start_end.m_hyper(m_hyper_jumps_idx,:)=0;
    glucose_events_info.time_of_events.m_hyper(m_hyper_jumps_idx,:)=0;
    end
end

if exist('N_S_hyper_jumps_length','var')
    for i = 1:size(N_S_hyper_jumps_length,1)
        if N_S_hyper_jumps_length(i,2) > 8
         S_hyper_jumps_idx(i,1)=i;
        end
    end
    if exist('S_hyper_jumps_idx','var')
    S_hyper_jumps_idx = nonzeros(S_hyper_jumps_idx);
    glucose_events_info.N_events.S_hyper = glucose_events_info.N_events.S_hyper - size(S_hyper_jumps_idx,1);
    glucose_events_info.events_length.S_hyper(S_hyper_jumps_idx)=0;  
    glucose_events_info.events_start_end.S_hyper(S_hyper_jumps_idx,:)=0;
    glucose_events_info.time_of_events.S_hyper(S_hyper_jumps_idx,:)=0;
    end
end


if exist('N_m_hypo_jumps_length','var')
    for i = 1:size(N_m_hypo_jumps_length,1)
        if N_m_hypo_jumps_length(i,2) > 8
         m_hypo_jumps_idx(i,1)=i;
        end
    end
    if exist('m_hypo_jumps_idx','var')
    m_hypo_jumps_idx = nonzeros(m_hypo_jumps_idx);
    glucose_events_info.N_events.m_hypo = glucose_events_info.N_events.m_hypo - size(m_hypo_jumps_idx,1);
    glucose_events_info.events_length.m_hypo(m_hypo_jumps_idx)=0;  
    glucose_events_info.events_start_end.m_hypo(m_hypo_jumps_idx,:)=0;
    glucose_events_info.time_of_events.m_hypo(m_hypo_jumps_idx,:)=0;
    end
end

if exist('N_S_hypo_jumps_length','var')
    for i = 1:size(N_S_hypo_jumps_length,1)
        if N_S_hypo_jumps_length(i,2) > 8
         S_hypo_jumps_idx(i,1)=i;
        end
    end
    if exist('S_hypo_jumps_idx','var')
    S_hypo_jumps_idx = nonzeros(S_hypo_jumps_idx);
    glucose_events_info.N_events.S_hypo = glucose_events_info.N_events.S_hypo - size(S_hypo_jumps_idx,1);
    glucose_events_info.events_length.S_hypo(S_hypo_jumps_idx)=0;  
    glucose_events_info.events_start_end.S_hypo(S_hypo_jumps_idx,:)=0;
    glucose_events_info.time_of_events.S_hypo(S_hypo_jumps_idx,:)=0;
    end
end

%% remove non zeros?

if size(find(glucose_events_info.events_start_end.S_hyper(:,1)),1) > 0
    non_zero_idx = find(glucose_events_info.events_start_end.S_hyper(:,1));
    glucose_events_info.events_start_end.S_hyper = glucose_events_info.events_start_end.S_hyper(non_zero_idx,:);
    glucose_events_info.time_of_events.S_hyper = glucose_events_info.time_of_events.S_hyper(non_zero_idx,:);
    glucose_events_info.events_length.S_hyper = glucose_events_info.events_length.S_hyper(non_zero_idx,:);
    glucose_events_info.glucose_data_gaps.S_hyper = glucose_events_info.glucose_data_gaps.S_hyper(non_zero_idx,:);
end

if size(find(glucose_events_info.events_start_end.m_hyper(:,1)),1) > 0
    non_zero_idx = find(glucose_events_info.events_start_end.m_hyper(:,1));
    glucose_events_info.events_start_end.m_hyper = glucose_events_info.events_start_end.m_hyper(non_zero_idx,:);
    glucose_events_info.time_of_events.m_hyper = glucose_events_info.time_of_events.m_hyper(non_zero_idx,:);
    glucose_events_info.events_length.m_hyper = glucose_events_info.events_length.m_hyper(non_zero_idx,:);
    glucose_events_info.glucose_data_gaps.m_hyper = glucose_events_info.glucose_data_gaps.m_hyper(non_zero_idx,:);
end

if size(find(glucose_events_info.events_start_end.S_hypo(:,1)),1) > 0
    non_zero_idx = find(glucose_events_info.events_start_end.S_hypo(:,1));
    glucose_events_info.events_start_end.S_hypo = glucose_events_info.events_start_end.S_hypo(non_zero_idx,:);
    glucose_events_info.time_of_events.S_hypo = glucose_events_info.time_of_events.S_hypo(non_zero_idx,:);
    glucose_events_info.events_length.S_hypo = glucose_events_info.events_length.S_hypo(non_zero_idx,:);
    glucose_events_info.glucose_data_gaps.S_hypo = glucose_events_info.glucose_data_gaps.S_hypo(non_zero_idx,:);
end

if size(find(glucose_events_info.events_start_end.m_hypo(:,1)),1) > 0
    non_zero_idx = find(glucose_events_info.events_start_end.m_hypo(:,1));
    glucose_events_info.events_start_end.m_hypo = glucose_events_info.events_start_end.m_hypo(non_zero_idx,:);
    glucose_events_info.time_of_events.m_hypo = glucose_events_info.time_of_events.m_hypo(non_zero_idx,:);
    glucose_events_info.events_length.m_hypo = glucose_events_info.events_length.m_hypo(non_zero_idx,:);
    glucose_events_info.glucose_data_gaps.m_hypo = glucose_events_info.glucose_data_gaps.m_hypo(non_zero_idx,:);
end



%%
S_hyperN = glucose_events_info.N_events.S_hyper;
m_hyperN = glucose_events_info.N_events.m_hyper;
S_hypoN = glucose_events_info.N_events.S_hypo;
m_hypoN = glucose_events_info.N_events.m_hypo;

hyperSm_hypoSm_N = [S_hyperN m_hyperN S_hypoN m_hypoN];

% if size(glucose_events_info.events_start_end.S_hypo,1) > 1
%     glucose_events_info.events_start_end.S_hypo = nonzeros(glucose_events_info.events_start_end.S_hypo);
%     glucose_events_info.events_length.S_hypo = nonzeros(glucose_events_info.events_length.S_hypo);
%     %glucose_events_info.time_of_events.S_hypo  
% end

%%
% figure()
% plot(universal_time,glucoseC)
% xlabel('time / s')
% ylabel('Glucose C mg/Dl')
% hold on
% yline(S_hyperG);
% yline(m_hyperG);
% yline(S_hypoG);
% yline(m_hypoG);
% title(" "+name_spreadsheet+" hyper S m hypo S m "+num2str(hyperSm_hypoSm_N))

%%% if size(glucose_events_info.events_start_end.m_hyper,1) > 1
%%%     if glucose_events_info.events_start_end.m_hyper ~= 0
%%%     glucose_events_info.events_start_end.m_hyper(:,1) = nonzeros(glucose_events_info.events_start_end.m_hyper(:,1));
%%%     glucose_events_info.events_start_end.m_hyper(:,2) = nonzeros(glucose_events_info.events_start_end.m_hyper(:,2));
%%%     end
%%% end

% 
%%% if size(glucose_events_info.events_start_end.S_hyper,1) > 1
%%%     if glucose_events_info.events_start_end.S_hyper ~= 0
%%%     glucose_events_info.events_start_end.S_hyper(:,1) = nonzeros(glucose_events_info.events_start_end.S_hyper(:,1));
%%%     glucose_events_info.events_start_end.S_hyper(:,2) = nonzeros(glucose_events_info.events_start_end.S_hyper(:,2));
%%%     end
%%% end

% if size(glucose_events_info.events_start_end.m_hyper,1) >= 1
%     if size( find(glucose_events_info.events_start_end.m_hyper),2) > 0;
%     glucose_events_info.events_start_end.plot.m_hyper(:,1) = nonzeros(glucose_events_info.events_start_end.m_hyper(:,1));
%     glucose_events_info.events_start_end.plot.m_hyper(:,2) = nonzeros(glucose_events_info.events_start_end.m_hyper(:,2));
%     end
% end
% if size(glucose_events_info.events_start_end.S_hyper,1) >= 1
%     if size( find(glucose_events_info.events_start_end.S_hyper),2) > 0;
%     glucose_events_info.events_start_end.plot.S_hyper(:,1) = nonzeros(glucose_events_info.events_start_end.S_hyper(:,1));
%     glucose_events_info.events_start_end.plot.S_hyper(:,2) = nonzeros(glucose_events_info.events_start_end.S_hyper(:,2));
%     end
% end
% if size(glucose_events_info.events_start_end.m_hypo,1) >= 1
%     if size( find(glucose_events_info.events_start_end.m_hypo),2) > 0;
%     glucose_events_info.events_start_end.plot.m_hypo(:,1) = nonzeros(glucose_events_info.events_start_end.m_hypo(:,1));
%     glucose_events_info.events_start_end.plot.m_hypo(:,2) = nonzeros(glucose_events_info.events_start_end.m_hypo(:,2));
%     end
% end
% if size(glucose_events_info.events_start_end.S_hypo,1) >= 1
%     if size( find(glucose_events_info.events_start_end.S_hypo),2) > 0;
%     glucose_events_info.events_start_end.plot.S_hypo(:,1) = nonzeros(glucose_events_info.events_start_end.S_hypo(:,1));
%     glucose_events_info.events_start_end.plot.S_hypo(:,2) = nonzeros(glucose_events_info.events_start_end.S_hypo(:,2));
%     end
% end

%final overide of glucose events
if glucose_events_info.N_events.m_hyper == 0
    glucose_events_info.N_events.m_hyper = 0;
    glucose_events_info.events_length.m_hyper = 0;  
    glucose_events_info.events_start_end.m_hyper = 0;
    glucose_events_info.time_of_events.m_hyper = 0;
end
if glucose_events_info.N_events.S_hyper == 0
    glucose_events_info.N_events.S_hyper = 0;
    glucose_events_info.events_length.S_hyper = 0;  
    glucose_events_info.events_start_end.S_hyper = 0;
    glucose_events_info.time_of_events.S_hyper = 0;
end
if glucose_events_info.N_events.m_hypo == 0
    glucose_events_info.N_events.m_hypo = 0;
    glucose_events_info.events_length.m_hypo = 0;  
    glucose_events_info.events_start_end.m_hypo = 0;
    glucose_events_info.time_of_events.m_hypo = 0;
end
if glucose_events_info.N_events.S_hypo == 0
    glucose_events_info.N_events.S_hypo = 0;
    glucose_events_info.events_length.S_hypo = 0;  
    glucose_events_info.events_start_end.S_hypo = 0;
    glucose_events_info.time_of_events.S_hypo = 0;
end

if exist('unique_m_hypo','var')
    if size(unique_m_hypo,2) == 0;
        glucose_events_info.N_events.m_hypo = 0;
        glucose_events_info.events_length.m_hypo = 0;  
        glucose_events_info.events_start_end.m_hypo = 0;
        glucose_events_info.time_of_events.m_hypo = 0;
    end
end
if exist('unique_m_hyper','var')
    if size(unique_m_hyper,2) == 0;
        glucose_events_info.N_events.m_hyper = 0;
        glucose_events_info.events_length.m_hyper = 0;  
        glucose_events_info.events_start_end.m_hyper = 0;
        glucose_events_info.time_of_events.m_hyper = 0;
    end
end
S_hyperN = glucose_events_info.N_events.S_hyper;
m_hyperN = glucose_events_info.N_events.m_hyper;
S_hypoN = glucose_events_info.N_events.S_hypo;
m_hypoN = glucose_events_info.N_events.m_hypo;

hyperSm_hypoSm_N = [S_hyperN m_hyperN S_hypoN m_hypoN];
% glucose_events_info.N_events.S_hyper
% glucose_events_info.N_events.m_hypo
% glucose_events_info.N_events.S_hypo


% if size(glucose_events_info.events_start_end.S_hypo,1) > 1
%     if glucose_events_info.events_start_end.S_hypo ~= 0
%     glucose_events_info.events_start_end.S_hypo(:,1) = nonzeros(glucose_events_info.events_start_end.S_hypo(:,1));
%     glucose_events_info.events_start_end.S_hypo(:,2) = nonzeros(glucose_events_info.events_start_end.S_hypo(:,2));
%     end
% end

%% combine events if they are overlapping. 

if size(glucose_events_info.events_start_end.m_hyper,1) > 1;
    for i= 2:size(glucose_events_info.events_start_end.m_hyper,1) 
            if glucose_events_info.events_start_end.m_hyper(i,1) <= glucose_events_info.events_start_end.m_hyper(i-1,2)
                glucose_events_info.events_start_end.m_hyper(i,1) = glucose_events_info.events_start_end.m_hyper(i-1,1);
                glucose_events_info.events_start_end.m_hyper(i-1,2) = glucose_events_info.events_start_end.m_hyper(i,2);
            end
    end
    %unique_m_hyper2 = [ unique(glucose_events_info.events_start_end.m_hyper(:,1)) unique(glucose_events_info.events_start_end.m_hyper(:,2))]  ;
    %glucose_events_info.events_start_end.m_hyper = unique_m_hyper2;
end
if size(glucose_events_info.events_start_end.m_hyper,1) > 1;
    for i= 2:size(glucose_events_info.events_start_end.m_hyper,1) 
            if glucose_events_info.events_start_end.m_hyper(i,1) <= glucose_events_info.events_start_end.m_hyper(i-1,2)
                glucose_events_info.events_start_end.m_hyper(i,1) = glucose_events_info.events_start_end.m_hyper(i-1,1);
                glucose_events_info.events_start_end.m_hyper(i-1,2) = glucose_events_info.events_start_end.m_hyper(i,2);
            end
    end
    unique_m_hyper2 = [ unique(glucose_events_info.events_start_end.m_hyper(:,1)) unique(glucose_events_info.events_start_end.m_hyper(:,2))]  ;
    glucose_events_info.events_start_end.m_hyper = unique_m_hyper2;
end
%%
if size(glucose_events_info.events_start_end.S_hyper,1) > 1;
    for i= 2:size(glucose_events_info.events_start_end.S_hyper,1)
        if glucose_events_info.events_start_end.S_hyper(i,1)~=0
                if glucose_events_info.events_start_end.S_hyper(i,1) <= glucose_events_info.events_start_end.S_hyper(i-1,2)
                    glucose_events_info.events_start_end.S_hyper(i,1) = glucose_events_info.events_start_end.S_hyper(i-1,1);
                    glucose_events_info.events_start_end.S_hyper(i-1,2) = glucose_events_info.events_start_end.S_hyper(i,2);
                end
        end
    end
    %unique_S_hyper2 = [ unique(glucose_events_info.events_start_end.S_hyper(:,1)) unique(glucose_events_info.events_start_end.S_hyper(:,2))]  ;
    %glucose_events_info.events_start_end.S_hyper = unique_S_hyper2;
end
if size(glucose_events_info.events_start_end.S_hyper,1) > 1;
        for i= 2:size(glucose_events_info.events_start_end.S_hyper,1)
            if glucose_events_info.events_start_end.S_hyper(i,1)~=0
                if glucose_events_info.events_start_end.S_hyper(i,1) <= glucose_events_info.events_start_end.S_hyper(i-1,2)
                    glucose_events_info.events_start_end.S_hyper(i,1) = glucose_events_info.events_start_end.S_hyper(i-1,1);
                    glucose_events_info.events_start_end.S_hyper(i-1,2) = glucose_events_info.events_start_end.S_hyper(i,2);
                end
            end
        end
    unique_S_hyper2 = [ unique(glucose_events_info.events_start_end.S_hyper(:,1)) unique(glucose_events_info.events_start_end.S_hyper(:,2))]  ;
    glucose_events_info.events_start_end.S_hyper = unique_S_hyper2;
end
%%
if size(glucose_events_info.events_start_end.m_hypo,1) > 1;
    for i= 2:size(glucose_events_info.events_start_end.m_hypo,1) 
            if glucose_events_info.events_start_end.m_hypo(i,1) <= glucose_events_info.events_start_end.m_hypo(i-1,2)
                glucose_events_info.events_start_end.m_hypo(i,1) = glucose_events_info.events_start_end.m_hypo(i-1,1);
                glucose_events_info.events_start_end.m_hypo(i-1,2) = glucose_events_info.events_start_end.m_hypo(i,2);
                
                %glucose_events_info.time_of_events.m_hypo(i,1) = glucose_events_info.time_of_events.m_hypo(i-1,1);
                %glucose_events_info.time_of_events.m_hypo(i-1,2) = glucose_events_info.time_of_events.m_hypo(i,2);
            end
    end
    %unique_m_hypo2 = [ unique(glucose_events_info.events_start_end.m_hypo(:,1)) unique(glucose_events_info.events_start_end.m_hypo(:,2))];
    %glucose_events_info.events_start_end.m_hypo = unique_m_hypo2;
end
if size(glucose_events_info.events_start_end.m_hypo,1) > 1;
    for i= 2:size(glucose_events_info.events_start_end.m_hypo,1) 
            if glucose_events_info.events_start_end.m_hypo(i,1) <= glucose_events_info.events_start_end.m_hypo(i-1,2)
                glucose_events_info.events_start_end.m_hypo(i,1) = glucose_events_info.events_start_end.m_hypo(i-1,1);
                glucose_events_info.events_start_end.m_hypo(i-1,2) = glucose_events_info.events_start_end.m_hypo(i,2);
                
                %glucose_events_info.time_of_events.m_hypo(i,1) = glucose_events_info.time_of_events.m_hypo(i-1,1);
                %glucose_events_info.time_of_events.m_hypo(i-1,2) = glucose_events_info.time_of_events.m_hypo(i,2);
            end
    end
    %hard programme solution for PD38
    %undo % to run the loop for PD38
    %if name_spreadsheet == 'PD38.xlsx'
        %'ddd'
     %   unique_m_hypo2_end_PD38 = unique(glucose_events_info.events_start_end.m_hypo(:,2))
     %  unique_m_hypo2_end_PD38 = [unique_m_hypo2_end_PD38(1:6,1); unique_m_hypo2_end_PD38(8:11,1)];
     %   unique_m_hypo2 = [ unique(glucose_events_info.events_start_end.m_hypo(:,1)) unique_m_hypo2_end_PD38];


    %else
    unique_m_hypo2 = [ unique(glucose_events_info.events_start_end.m_hypo(:,1)) unique(glucose_events_info.events_start_end.m_hypo(:,2))];
    glucose_events_info.events_start_end.m_hypo = unique_m_hypo2;
    %end
end
%%
if size(glucose_events_info.events_start_end.S_hypo,1) > 1;
    for i= 2:size(glucose_events_info.events_start_end.S_hypo,1) 
            if glucose_events_info.events_start_end.S_hypo(i,1) <= glucose_events_info.events_start_end.S_hypo(i-1,2)
                glucose_events_info.events_start_end.S_hypo(i,1) = glucose_events_info.events_start_end.S_hypo(i-1,1);
                glucose_events_info.events_start_end.S_hypo(i-1,2) = glucose_events_info.events_start_end.S_hypo(i,2);
            end
    end
    %unique_S_hypo2 = [ unique(glucose_events_info.events_start_end.S_hypo(:,1)) unique(glucose_events_info.events_start_end.S_hypo(:,2))]  ;
    %glucose_events_info.events_start_end.S_hypo = unique_S_hypo2;
end
if size(glucose_events_info.events_start_end.S_hypo,1) > 1;
    for i= 2:size(glucose_events_info.events_start_end.S_hypo,1) 
            if glucose_events_info.events_start_end.S_hypo(i,1) <= glucose_events_info.events_start_end.S_hypo(i-1,2)
                glucose_events_info.events_start_end.S_hypo(i,1) = glucose_events_info.events_start_end.S_hypo(i-1,1);
                glucose_events_info.events_start_end.S_hypo(i-1,2) = glucose_events_info.events_start_end.S_hypo(i,2);
            end
    end
    unique_S_hypo2 = [ unique(glucose_events_info.events_start_end.S_hypo(:,1)) unique(glucose_events_info.events_start_end.S_hypo(:,2))]  ;
    glucose_events_info.events_start_end.S_hypo = unique_S_hypo2;
end

%% adjust new lengths
if exist('unique_m_hypo2','var')
for i=1:size(glucose_events_info.events_start_end.m_hypo,1)
    unique_length_m_hypo(i,1) = glucose_events_info.events_start_end.m_hypo(i,2) - glucose_events_info.events_start_end.m_hypo(i,1);
end
glucose_events_info.N_events.m_hypo = size(glucose_events_info.events_start_end.m_hypo,1);
glucose_events_info.events_length.m_hypo = unique_length_m_hypo;
glucose_events_info.N_events.m_hypo = size(unique_length_m_hypo,1);
glucose_events_info.N_events.m_hypo = size(nonzeros(unique_length_m_hypo(:,1)),1);
end

if exist('unique_S_hypo2','var')
for i=1:size(glucose_events_info.events_start_end.S_hypo,1)
    unique_length_S_hypo(i,1) = glucose_events_info.events_start_end.S_hypo(i,2) - glucose_events_info.events_start_end.S_hypo(i,1);
end
glucose_events_info.N_events.S_hypo = size(glucose_events_info.events_start_end.S_hypo,1);
glucose_events_info.events_length.S_hypo = unique_length_S_hypo;
glucose_events_info.N_events.S_hypo = size(unique_length_S_hypo,1);
glucose_events_info.N_events.S_hypo = size(nonzeros(unique_length_S_hypo(:,1)),1);
end

if exist('unique_m_hyper2','var')
for i=1:size(glucose_events_info.events_start_end.m_hyper,1)
    unique_length_m_hyper(i,1) = glucose_events_info.events_start_end.m_hyper(i,2) - glucose_events_info.events_start_end.m_hyper(i,1);
end
glucose_events_info.N_events.m_hyper = size(glucose_events_info.events_start_end.m_hyper,1);
glucose_events_info.events_length.m_hyper = unique_length_m_hyper;
glucose_events_info.N_events.m_hyper = size(unique_length_m_hyper,1);
glucose_events_info.N_events.m_hyper = size(nonzeros(unique_length_m_hyper(:,1)),1);
end
%%
if exist('unique_S_hyper2','var')
for i=1:size(glucose_events_info.events_start_end.S_hyper,1)
    unique_length_S_hyper(i,1) = glucose_events_info.events_start_end.S_hyper(i,2) - glucose_events_info.events_start_end.S_hyper(i,1);
end
glucose_events_info.N_events.S_hyper = size(glucose_events_info.events_start_end.S_hyper,1);
glucose_events_info.events_length.S_hyper = unique_length_S_hyper;
glucose_events_info.N_events.S_hyper = size(unique_length_S_hyper,1);
glucose_events_info.N_events.S_hyper = size(nonzeros(unique_length_S_hyper(:,1)),1);
end


S_hyperN = glucose_events_info.N_events.S_hyper;
m_hyperN = glucose_events_info.N_events.m_hyper;
S_hypoN = glucose_events_info.N_events.S_hypo;
m_hypoN = glucose_events_info.N_events.m_hypo;

hyperSm_hypoSm_N = [S_hyperN m_hyperN S_hypoN m_hypoN];

%% Add date to glucose events info!!

if size(find(glucose_events_info.events_start_end.m_hypo),2) > 0 
    glucose_events_info.date_of_events.m_hypo(:,1) = date(glucose_events_info.events_start_end.m_hypo(:,1));
    glucose_events_info.date_of_events.m_hypo(:,2) = date(glucose_events_info.events_start_end.m_hypo(:,2));
end

if size(find(glucose_events_info.events_start_end.S_hypo),2) > 0 
    glucose_events_info.date_of_events.S_hypo(:,1) = date(glucose_events_info.events_start_end.S_hypo(:,1));
    glucose_events_info.date_of_events.S_hypo(:,2) = date(glucose_events_info.events_start_end.S_hypo(:,2));
end

if size(find(glucose_events_info.events_start_end.m_hyper),2) > 0 
    glucose_events_info.date_of_events.m_hyper(:,1) = date(glucose_events_info.events_start_end.m_hyper(:,1));
    glucose_events_info.date_of_events.m_hyper(:,2) = date(glucose_events_info.events_start_end.m_hyper(:,2));
end

if size(find(glucose_events_info.events_start_end.S_hyper),2) > 0 
    glucose_events_info.date_of_events.S_hyper(:,1) = date(glucose_events_info.events_start_end.S_hyper(:,1));
    glucose_events_info.date_of_events.S_hyper(:,2) = date(glucose_events_info.events_start_end.S_hyper(:,2));
end

%% re add times, so it's just times updated
    % all time vector 
    all_time = time./scale_t;

    N_hours_all = fix( (all_time ./ (delta_t*3600) ));
    N_minutes_all = fix(mod( all_time, (delta_t*3600) ) ./ (delta_t*60) );
    N_seconds_all = floor(mod( mod( all_time, (delta_t*3600) ) ./ (delta_t*60) , 1)*60);

    

    time_HMS_ordered_all = [N_hours_all N_minutes_all N_seconds_all];

    time_rows = datetime(int2str(time_HMS_ordered_all), 'Format', 'HH mm ss');

%time_m_hyper = glucose_events_info.time_of_events.m_hyper ; 
glucose_events_info = rmfield(glucose_events_info,'time_of_events')

if size(find(glucose_events_info.events_start_end.m_hypo),2) > 0 
    %glucose_events_info.time_of_events.m_hypo(:,1:3) = time_HMS_ordered_all(glucose_events_info.events_start_end.m_hypo(:,1),:);
    %glucose_events_info.time_of_events.m_hypo(:,4:6) = time_HMS_ordered_all(glucose_events_info.events_start_end.m_hypo(:,2),:);
    glucose_events_info.time_of_events.m_hypo(:,1) = time_rows(glucose_events_info.events_start_end.m_hypo(:,1),:);
    glucose_events_info.time_of_events.m_hypo(:,2) = time_rows(glucose_events_info.events_start_end.m_hypo(:,2),:);
end

if size(find(glucose_events_info.events_start_end.S_hypo),2) > 0 
    %glucose_events_info.time_of_events.S_hypo(:,1:3) = time_HMS_ordered_all(glucose_events_info.events_start_end.S_hypo(:,1),:);
    %glucose_events_info.time_of_events.S_hypo(:,4:6) = time_HMS_ordered_all(glucose_events_info.events_start_end.S_hypo(:,2),:);
    glucose_events_info.time_of_events.S_hypo(:,1) = time_rows(glucose_events_info.events_start_end.S_hypo(:,1),:);
    glucose_events_info.time_of_events.S_hypo(:,2) = time_rows(glucose_events_info.events_start_end.S_hypo(:,2),:);
end

if size(find(glucose_events_info.events_start_end.m_hyper),2) > 0 
    %glucose_events_info.time_of_events.m_hyper(:,1:3) = time_HMS_ordered_all(glucose_events_info.events_start_end.m_hyper(:,1),:);
    %glucose_events_info.time_of_events.m_hyper(:,4:6) = time_HMS_ordered_all(glucose_events_info.events_start_end.m_hyper(:,2),:);
    glucose_events_info.time_of_events.m_hyper(:,1) = time_rows(glucose_events_info.events_start_end.m_hyper(:,1),:);
    glucose_events_info.time_of_events.m_hyper(:,2) = time_rows(glucose_events_info.events_start_end.m_hyper(:,2),:);
end

if size(find(glucose_events_info.events_start_end.S_hyper),2) > 0 
    %glucose_events_info.time_of_events.S_hyper(:,1:3) = time_HMS_ordered_all(glucose_events_info.events_start_end.S_hyper(:,1),:);
    %glucose_events_info.time_of_events.S_hyper(:,4:6) = time_HMS_ordered_all(glucose_events_info.events_start_end.S_hyper(:,2),:);
    glucose_events_info.time_of_events.S_hyper(:,1) = time_rows(glucose_events_info.events_start_end.S_hyper(:,1),:);
    glucose_events_info.time_of_events.S_hyper(:,2) = time_rows(glucose_events_info.events_start_end.S_hyper(:,2),:);
end



%% Adjust new plots
if size(glucose_events_info.events_start_end.m_hyper,1) >= 1
    if size( find(glucose_events_info.events_start_end.m_hyper),2) > 0;
    glucose_events_info.events_start_end.plot.m_hyper(:,1) = nonzeros(glucose_events_info.events_start_end.m_hyper(:,1));
    glucose_events_info.events_start_end.plot.m_hyper(:,2) = nonzeros(glucose_events_info.events_start_end.m_hyper(:,2));
    end
end
if size(glucose_events_info.events_start_end.S_hyper,1) >= 1
    if size( find(glucose_events_info.events_start_end.S_hyper),2) > 0;
    glucose_events_info.events_start_end.plot.S_hyper(:,1) = nonzeros(glucose_events_info.events_start_end.S_hyper(:,1));
    glucose_events_info.events_start_end.plot.S_hyper(:,2) = nonzeros(glucose_events_info.events_start_end.S_hyper(:,2));

    glucose_events_info.events_length.plot.S_hyper(:,1) = nonzeros(glucose_events_info.events_length.S_hyper(:,1));
    
    end
end


if size(glucose_events_info.events_start_end.m_hypo,1) >= 1
    if size( find(glucose_events_info.events_start_end.m_hypo),2) > 0;
    glucose_events_info.events_start_end.plot.m_hypo(:,1) = nonzeros(glucose_events_info.events_start_end.m_hypo(:,1));
    glucose_events_info.events_start_end.plot.m_hypo(:,2) = nonzeros(glucose_events_info.events_start_end.m_hypo(:,2));
    end
end
if size(glucose_events_info.events_start_end.S_hypo,1) >= 1
    if size( find(glucose_events_info.events_start_end.S_hypo),2) > 0;
    glucose_events_info.events_start_end.plot.S_hypo(:,1) = nonzeros(glucose_events_info.events_start_end.S_hypo(:,1));
    glucose_events_info.events_start_end.plot.S_hypo(:,2) = nonzeros(glucose_events_info.events_start_end.S_hypo(:,2));
    end
end
%% add glucose data

glucose_events_info.glucose = glucoseC;


%%
figure()
plot(universal_time,glucoseC,'k')
xlabel('time / s')
ylabel('Glucose C mg/Dl')
hold on
yline(S_hyperG,'--');
yline(m_hyperG,'--');
yline(S_hypoG,'--');
yline(m_hypoG,'--');
title(" "+name_spreadsheet+" hyper S m hypo S m "+num2str(hyperSm_hypoSm_N))
if size( find(glucose_events_info.events_start_end.S_hyper),2) > 0 %glucose_events_info.events_start_end.S_hyper ~= 0
    plot(universal_time(glucose_events_info.events_start_end.plot.S_hyper(:,1)) , glucoseC(glucose_events_info.events_start_end.plot.S_hyper(:,1)),'m.','MarkerSize',12)
    plot(universal_time(glucose_events_info.events_start_end.plot.S_hyper(:,2)) , glucoseC(glucose_events_info.events_start_end.plot.S_hyper(:,2)),'mx','MarkerSize',12)
end
if size( find(glucose_events_info.events_start_end.m_hyper),2) > 0%glucose_events_info.events_start_end.m_hyper ~= 0
    plot(universal_time(glucose_events_info.events_start_end.plot.m_hyper(:,1)) , glucoseC(glucose_events_info.events_start_end.plot.m_hyper(:,1)),'r.','MarkerSize',12)
    plot(universal_time(glucose_events_info.events_start_end.plot.m_hyper(:,2)) , glucoseC(glucose_events_info.events_start_end.plot.m_hyper(:,2)),'rx','MarkerSize',12)
end
if size( find(glucose_events_info.events_start_end.S_hypo),2) > 0%glucose_events_info.events_start_end.S_hypo ~= 0
    plot(universal_time(glucose_events_info.events_start_end.plot.S_hypo(:,1)) , glucoseC(glucose_events_info.events_start_end.plot.S_hypo(:,1)),'c.','MarkerSize',12)
    plot(universal_time(glucose_events_info.events_start_end.plot.S_hypo(:,2)) , glucoseC(glucose_events_info.events_start_end.plot.S_hypo(:,2)),'cx','MarkerSize',12)
end
if size( find(glucose_events_info.events_start_end.m_hypo),2) > 0%glucose_events_info.events_start_end.m_hypo ~= 0
    plot(universal_time(glucose_events_info.events_start_end.plot.m_hypo(:,1)) , glucoseC(glucose_events_info.events_start_end.plot.m_hypo(:,1)),'b.','MarkerSize',12)
    plot(universal_time(glucose_events_info.events_start_end.plot.m_hypo(:,2)) , glucoseC(glucose_events_info.events_start_end.plot.m_hypo(:,2)),'bx','MarkerSize',12)
end
if start_gap ~= 0
    plot(universal_time(start_gap) , glucoseC(start_gap),'<','MarkerSize',12)
    plot(universal_time(end_gap) , glucoseC(end_gap),'>','MarkerSize',12)
end
%legend()

% figure()
% plot(universal_time,glucoseC)
% xlabel('time / s')
% ylabel('Glucose C mg/Dl')
% hold on
% yline(S_hyperG);
% yline(m_hyperG);
% yline(S_hypoG);
% yline(m_hypoG);
% plot(universal_time(ordered_start_m_hypo),glucoseC(ordered_start_m_hypo),'rx')
% plot(universal_time(ordered_end_m_hypo),glucoseC(ordered_end_m_hypo),'bx')

% figure()
% plot(universal_time,glucoseC)
% xlabel('time / s')
% ylabel('Glucose C mg/Dl')
% hold on
% yline(S_hyperG);
% yline(m_hyperG);
% yline(S_hypoG);
% yline(m_hypoG);
% plot(universal_time(m_hyper_point),glucoseC(m_hyper_point),'b+')
% plot(universal_time(ordered_start_m_hyper),glucoseC(ordered_start_m_hyper),'r*')
% plot(universal_time(ordered_end_m_hyper),glucoseC(ordered_end_m_hyper),'go')


    
end
%%
function hyperSm_hypoSm_N = readglucose_number(name_spreadsheet)

S_hyperG = 180;
m_hyperG = 144;


S_hypoG =47;
m_hypoG =72;

spreadsheet_data = readcell(name_spreadsheet);
%spreadsheet_data = readcell('test_time.xlsx');

for i=1:size(spreadsheet_data,1)
    %glucoseC(i)= (str2double (cell2mat(spreadsheet_data(i,3)))/100)' ; % div 100 as glucose values have 2 DP of 0.
    glucoseC(i)= ( (cell2mat(spreadsheet_data(i,3))))' ; % if glucose
    %saved as number

    date(i)=spreadsheet_data{i,1};
    time(i)=spreadsheet_data{i,2};
end
glucoseC=glucoseC';
%date(1:size(spreadsheet_data,1),1)=spreadsheet_data{1:size(spreadsheet_data,1),1};
%time(1:size(spreadsheet_data,1),1)=spreadsheet_data{1:size(spreadsheet_data,1),2};
date=date';
time= time';

delta_t = 1.157407407407407e-05;
s_in_day = 86400;
scale_t = 1/delta_t;
time = time*scale_t;

days = unique(date);
for i=1:size(date,1)
    dayN(i) = find(date(i)==days);
end
dayN = dayN';

universal_time = time + (s_in_day*(dayN-1));


S_hyper_point = find(glucoseC>S_hyperG);
m_hyper_point = find(glucoseC<S_hyperG & glucoseC>m_hyperG);

S_hypo_point = find(glucoseC<S_hypoG);
m_hypo_point = find(glucoseC>S_hypoG & glucoseC<m_hypoG);

S_hyperN =1;
S_hypoN =1;

m_hyperN =1;
m_hypoN =1;

if size(S_hyper_point,1) == 0
    S_hyperN =0;
else
    for i=2:size(S_hyper_point,1)

        if S_hyper_point(i)  ~= S_hyper_point(i-1)+1
            S_hyperN = S_hyperN + 1;
        end
    end
end

if size(S_hypo_point,1) == 0
    S_hypoN =0;
else

    for i=2:size(S_hypo_point,1)

        if S_hypo_point(i)  ~= S_hypo_point(i-1)+1
            S_hypoN = S_hypoN + 1;
        end
    end
end

if size(m_hyper_point,1) == 0
    m_hyperN =0;
else
    for i=2:size(m_hyper_point,1)

        if m_hyper_point(i)  ~= m_hyper_point(i-1)+1
            m_hyperN = m_hyperN + 1;
        end
    end
end

if size(m_hypo_point,1) == 0
    m_hypoN =0;
else

    for i=2:size(m_hypo_point,1)

        if m_hypo_point(i)  ~= m_hypo_point(i-1)+1
            m_hypoN = m_hypoN + 1;
        end
    end
end

% 'S Hyper N'
% S_hyperN
% 
% 'm Hyper N'
% m_hyperN
% 
% 'S Hypo N'
% S_hypoN
% 
% 'm Hypo N'
% m_hypoN

hyperSm_hypoSm_N = [S_hyperN m_hyperN S_hypoN m_hypoN];

figure()
plot(universal_time,glucoseC)
xlabel('time / s')
ylabel('Glucose C mg/Dl')
hold on
yline(S_hyperG);
yline(m_hyperG);
yline(S_hypoG);
yline(m_hypoG);
title(" "+name_spreadsheet+" hyper S m hypo S m "+num2str(hyperSm_hypoSm_N))


end


%% 

% PD10.events_start_end.m_hypo  
% 
% PD10.events_start_end.m_hypo(1,1)
% PD10.events_start_end.m_hypo(1,2)
% 
% min_m_hypo = find(PD10.glucose(PD10.events_start_end.m_hypo(1,1):PD10.events_start_end.m_hypo(1,2),1) == min(PD10.glucose(PD10.events_start_end.m_hypo(1,1):PD10.events_start_end.m_hypo(1,2),1)));
% 
% min_m_hypo*5 %N minutes after the event when the maximum hypo occured. 


%PD10.events_start_end.m_hypo(1,1)+min_m_hypo(1)

%15*5

