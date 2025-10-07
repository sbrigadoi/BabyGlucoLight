close all
clear all
clc


%%
% % idx_PD = [4;5;8;9;11;14;15;19;47];
idx_PD = [4;5;8;9;11;15;19;47];
% idx_PD = 14;


name_first = 'PD_%d_FIRST_rsFC.mat';
name_last  = 'PD_%d_LAST_rsFC.mat';

rsFC_hbo_first = [];
rsFC_hbr_first = [];
rsFC_hbo_last = [];
rsFC_hbr_last = [];

for i = 1:1:length(idx_PD)
    
    curr_idx = idx_PD(i);

    % First euglycemia rsFC
    tmp_first_name = sprintf(name_first, curr_idx);
    load(tmp_first_name)
    disp(['LOAD : ',tmp_first_name])
    disp(' ')
    tmp_hbo_zFC_first = (FC_struct.hbo_FC);
    tmp_hbr_zFC_first = (FC_struct.hbr_FC);
    tmp_hbo_zFC_first = atanh(FC_struct.hbo_FC);
    tmp_hbr_zFC_first = atanh(FC_struct.hbr_FC);
    % tmp_hbo_zFC_first = zscore(FC_struct.hbo_FC);
    % tmp_hbr_zFC_first= zscore(FC_struct.hbr_FC);

    % Last Euglycemia interval
    tmp_last_name = sprintf(name_last, curr_idx);
    load(tmp_last_name)
    disp(['LOAD : ',tmp_last_name])
    disp(' ')
    tmp_hbo_zFC_last = (FC_struct.hbo_FC);
    tmp_hbr_zFC_last = (FC_struct.hbr_FC);
    tmp_hbo_zFC_last = atanh(FC_struct.hbo_FC);
    tmp_hbr_zFC_last = atanh(FC_struct.hbr_FC);
    % tmp_hbo_zFC_last = zscore(FC_struct.hbo_FC);
    % tmp_hbr_zFC_last = zscore(FC_struct.hbr_FC);


    tmp_up = triu(tmp_hbo_zFC_first);
    At = tmp_up .';
    m  = (1:size(At,1)).' >= (1:size(At,2));
    tmp_hbo_first_up  = At(m);
    tmp_hbo_first_up = tmp_hbo_first_up(~isinf(tmp_hbo_first_up))';
    % to_remove = find(tmp_hbo_first_up==1);
    % tmp_hbo_first_up(to_remove) = [];    
    % tmp_hbo_first_up = tmp_hbo_first_up';
    if length(tmp_hbo_first_up) < 36
        to_add = NaN*ones(1,36-length(tmp_hbo_first_up));
        tmp_hbo_first_up = [to_add tmp_hbo_first_up]
    end
    rsFC_hbo_first = [rsFC_hbo_first ; tmp_hbo_first_up];

    tmp_up = triu(tmp_hbr_zFC_first);
    At = tmp_up .';
    m  = (1:size(At,1)).' >= (1:size(At,2));
    tmp_hbr_first_up  = At(m);
    tmp_hbr_first_up = tmp_hbr_first_up(~isinf(tmp_hbr_first_up))';
    % to_remove = find(tmp_hbr_first_up==1);
    % tmp_hbr_first_up(to_remove) = [];    
    % tmp_hbr_first_up = tmp_hbr_first_up';    
    if length(tmp_hbr_first_up) < 36
        to_add = NaN*ones(1,36-length(tmp_hbr_first_up));
        tmp_hbr_first_up = [to_add tmp_hbr_first_up]
    end
    rsFC_hbr_first = [rsFC_hbr_first ; tmp_hbr_first_up];

    tmp_up = triu(tmp_hbo_zFC_last);
    At = tmp_up .';
    m  = (1:size(At,1)).' >= (1:size(At,2));
    tmp_hbo_last_up  = At(m);
    tmp_hbo_last_up = tmp_hbo_last_up(~isinf(tmp_hbo_last_up))';
    % to_remove = find(tmp_hbo_last_up==1);
    % tmp_hbo_last_up(to_remove) = [];    
    % tmp_hbo_last_up = tmp_hbo_last_up';

    if length(tmp_hbo_last_up) < 36
        to_add = NaN*ones(1,36-length(tmp_hbo_last_up));
        tmp_hbo_last_up = [to_add tmp_hbo_last_up]
    end
    rsFC_hbo_last = [rsFC_hbo_last ; tmp_hbo_last_up];

    tmp_up = triu(tmp_hbr_zFC_last);
    At = tmp_up .';
    m  = (1:size(At,1)).' >= (1:size(At,2));
    tmp_hbr_last_up  = At(m);
    tmp_hbr_last_up = tmp_hbr_last_up(~isinf(tmp_hbr_last_up))';
    % to_remove = find(tmp_hbr_last_up==1);
    % tmp_hbr_last_up(to_remove) = [];    
    % tmp_hbr_last_up = tmp_hbr_last_up';    
    if length(tmp_hbr_last_up) < 36
        to_add = NaN*ones(1,36-length(tmp_hbr_last_up));
        tmp_hbr_last_up = [to_add tmp_hbr_last_up]
    end
    rsFC_hbr_last = [rsFC_hbr_last ; tmp_hbr_last_up];

end

%%
rsFC_metrics_label = string();

for i = 1:1:length(FC_struct.labels_short)
    for j = i+1:1:length(FC_struct.labels_short)
        tmp_label = strcat(FC_struct.labels_short(i)," - ",FC_struct.labels_short(j));
        rsFC_metrics_label = [rsFC_metrics_label tmp_label  ];
    end
end

rsFC_metrics_label(1) = [];

%%
% rsFC_hbo_first(:,1:15) = NaN;
% rsFC_hbr_last(:,1:15) = NaN;
% rsFC_hbo_first(:,1:15) = NaN;
% rsFC_hbr_last(:,1:15) = NaN;
% %%
% save rr.mat rsFC_hbo_first rsFC_hbo_last rsFC_hbr_first rsFC_hbr_last
%%

% to_remove = [1,10,18,25,31,36,40,43,45];
% rsFC_hbo_first(:,to_remove) = [];
% rsFC_hbr_first(:,to_remove) = [];
% rsFC_hbo_last(:,to_remove) = [];
% rsFC_hbr_last(:,to_remove) = [];
%%
% a = load('rr.mat');
% 
% rsFC_hbo_first = [rsFC_hbo_first;a.rsFC_hbo_first];
% rsFC_hbo_last = [rsFC_hbo_last;a.rsFC_hbo_last];
% rsFC_hbr_first = [rsFC_hbr_first;a.rsFC_hbr_first];
% rsFC_hbr_last = [rsFC_hbr_last;a.rsFC_hbr_last];
%%

[h_hbo,p_hbo] = ttest2(rsFC_hbo_first,rsFC_hbo_last)
[h_hbr,p_hbr,cc,stats] = ttest2(rsFC_hbr_first,rsFC_hbr_last)

idx_sig_hbo = find(h_hbo==1);

idx_sig_hbr = find(h_hbr==1);

disp('======================================')
disp(['HbO significative = ',num2str(length(idx_sig_hbo))])
if length(idx_sig_hbo)>0
    for i = length(idx_sig_hbo)
        disp([rsFC_metrics_label(idx_sig_hbr(i))])
    end
end
disp(['HbR significative = ',num2str(length(idx_sig_hbr))])
if length(idx_sig_hbr)>0
    for i = length(idx_sig_hbr)
        disp([rsFC_metrics_label(idx_sig_hbr(i))])
        disp(p_hbr(idx_sig_hbr(i)))
        disp(stats.tstat(idx_sig_hbr(i)))
        disp(stats.df(idx_sig_hbr(i)))
    end
end
