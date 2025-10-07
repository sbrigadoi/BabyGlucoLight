% CORRELAZIONE TRA DIFFERENZA DI rsFC E METRICHE GLICEMICHE

close all
clear all
clc

%% LOAD OF THE DATA

% data.mat contains the cgm_metrics, the difference between rsFC for both
% HbO and HbR (already z-scored) and the label of the ROI couples

% load data.mat  % PD14 WITH 2 EXCLUDED ROI (NaN)
% load data_all_ROI.mat  % PD14 ALL ROI

%load data.mat

load('rsFC_hbo_metrics.mat')
load('rsFC_hbr_metrics.mat')
load('cgm_metrics.mat')
load('rsFC_metrics_label.mat')

%run just OD PDs
%cgm_metrics = cgm_metrics([2;3;4;5;7;8;9;10;12],:);
%rsFC_hbo_metrics = rsFC_hbo_metrics([2;3;4;5;7;8;9;10;12],:);
%rsFC_hbr_metrics = rsFC_hbr_metrics([2;3;4;5;7;8;9;10;12],:);

%just add PD3
%cgm_metrics = cgm_metrics([1;2;3;4;5;7;8;9;10;12],:);
%rsFC_hbo_metrics = rsFC_hbo_metrics([1;2;3;4;5;7;8;9;10;12],:);
%rsFC_hbr_metrics = rsFC_hbr_metrics([1;2;3;4;5;7;8;9;10;12],:);

%just add PD10
%cgm_metrics = cgm_metrics([2;3;4;5;6;7;8;9;10;12],:);
%rsFC_hbo_metrics = rsFC_hbo_metrics([2;3;4;5;6;7;8;9;10;12],:);
%rsFC_hbr_metrics = rsFC_hbr_metrics([2;3;4;5;6;7;8;9;10;12],:);

%just add PD25
%cgm_metrics = cgm_metrics([2;3;4;5;7;8;9;10;11;12],:);
%rsFC_hbo_metrics = rsFC_hbo_metrics([2;3;4;5;7;8;9;10;11;12],:);
%rsFC_hbr_metrics = rsFC_hbr_metrics([2;3;4;5;7;8;9;10;11;12],:);

%add PD10 and PD25
%cgm_metrics = cgm_metrics([2;3;4;5;6;7;8;9;10;11;12],:);
%rsFC_hbo_metrics = rsFC_hbo_metrics([2;3;4;5;6;7;8;9;10;11;12],:);
%rsFC_hbr_metrics = rsFC_hbr_metrics([2;3;4;5;6;7;8;9;10;11;12],:);


%run with new PDS 3 10 25
%idx_PD    = [4;5;8;9;11;14;15;19;47];
%idx_PD    = [3;4;5;8;9;10;11;14;15;19;25;47];
%idx_PD([2;3;4;5;7;8;9;10;12])


%Guy update 10 04 25 - just use 8 basic metrics
cgm_metrics= [cgm_metrics(:,1:2) cgm_metrics(:,16:21)];


% From table to matrix and z-score of the cgm metrics
cgm_metrics_doub = table2array(cgm_metrics);
cgm_metrics_doub_z = zscore(cgm_metrics_doub);
% Extraction of cgm metrics names
cgm_metrics_label = string(cgm_metrics.Properties.VariableNames);% Canonical Correlation and Spearman's correlation



%% SINGOLA COPPIA - TUTTA CGM

% alfa = 0.05/23;
alfa = 0.05;

% HbO

r_hbo = [];
p_hbo = [];
count_tot_hbo = 0;

% Controllo se le metriche (colonne) del glucosio sono gaussiane, se anche
% una non lo è devo usare Spearman
for jj = 1:1:size(cgm_metrics_doub_z,2)
    if lillietest(cgm_metrics_doub_z(:,jj))==1
        % disp('Spearman')
        % disp(num2str(jj))
    end
end

disp('===================================================================')
disp('HbO')
disp(' ')

for i = 1:size(rsFC_hbo_metrics,2)

    count_tot_hbo = count_tot_hbo + 1;

    tmp_rsFC_hbo_metrics = rsFC_hbo_metrics(:,i);
    tmp_cgm_metrics_doub_z = cgm_metrics_doub_z;

    if sum(isnan(tmp_rsFC_hbo_metrics))>0
        to_remove = find(isnan(tmp_rsFC_hbo_metrics));
        tmp_rsFC_hbo_metrics(to_remove) = [];
        tmp_cgm_metrics_doub_z(to_remove,:) = [];
    end

    [r_hbo(i,:),p_hbo(i,:)] = corr(tmp_rsFC_hbo_metrics,tmp_cgm_metrics_doub_z,'type','Spearman');

    idx = find(p_hbo(i,:)<alfa);
    for jj = 1:length(idx)
        text = '%s with %-7s \t --> r = %f \t p-value = %f   --> SPEARMAN \n';
        fprintf(text,rsFC_metrics_label(1,i),cgm_metrics_label(1,idx(jj)),r_hbo(i,idx(jj)),p_hbo(i,idx(jj)))
    end

end

[riga,colonna] = find(p_hbo<alfa);
count_hbo = length(riga);

% disp('===================================================================')
% disp('HbO')
% disp(' ')
% for i = 1:length(riga)
%     text = '%s with %-7s \t --> r = %f \t p-value = %f   --> SPEARMAN \n';
%     fprintf(text,rsFC_metrics_label(1,riga(i)),cgm_metrics_label(1,colonna(i)),r_hbo(riga(i),colonna(i)),p_hbo(riga(i),colonna(i)))
% end


% HbR
r_hbr = [];
p_hbr = [];
count_tot_hbr = 0;

% Controllo se le metriche (colonne) del glucosio sono gaussiane, se anche
% una non lo è devo usare Spearman
% for jj = 1:1:size(cgm_metrics_doub_z)
%     if lillietest(cgm_metrics_doub_z(:,jj))==1
%         disp('Spearman')
%     end
% end

disp('===================================================================')
disp('HbR')
disp(' ')

for i = 1:size(rsFC_hbr_metrics,2)

    count_tot_hbr = count_tot_hbr + 1;

    tmp_rsFC_hbr_metrics = rsFC_hbr_metrics(:,i);
    tmp_cgm_metrics_doub_z = cgm_metrics_doub_z;

    if sum(isnan(tmp_rsFC_hbr_metrics))>0
        to_remove = find(isnan(tmp_rsFC_hbr_metrics));
        tmp_rsFC_hbr_metrics(to_remove) = [];
        tmp_cgm_metrics_doub_z(to_remove,:) = [];
    end

    [r_hbr(i,:),p_hbr(i,:)] = corr(tmp_rsFC_hbr_metrics,tmp_cgm_metrics_doub_z,'type','Spearman');
    
    idx = find(p_hbr(i,:)<alfa);
    for jj = 1:length(idx)
        text = '%s with %-7s \t --> r = %f \t p-value = %f   --> SPEARMAN \n';
        fprintf(text,rsFC_metrics_label(1,i),cgm_metrics_label(1,idx(jj)),r_hbr(i,idx(jj)),p_hbr(i,idx(jj)))
    end


end

[riga,colonna] = find(p_hbr<alfa);
count_hbr = length(riga);

% disp('===================================================================')
% disp('HbR')
% disp(' ')
% for i = 1:length(riga)
%     text = '%s with %-7s \t --> r = %f \t p-value = %f   --> SPEARMAN \n';
%     fprintf(text,rsFC_metrics_label(1,riga(i)),cgm_metrics_label(1,colonna(i)),r_hbr(riga(i),colonna(i)),p_hbr(riga(i),colonna(i)))
% end

count_tot_hbo = size(rsFC_hbo_metrics,2)*size(cgm_metrics_doub_z,2);
count_tot_hbr = size(rsFC_hbo_metrics,2)*size(cgm_metrics_doub_z,2);

disp('===================================================================')
disp(['ALFA = ',num2str(alfa)])
disp(['HbO SIGNIFICATIVE --> ',num2str(count_hbo),'/',num2str(count_tot_hbo)])
disp(['HbR SIGNIFICATIVE --> ',num2str(count_hbr),'/',num2str(count_tot_hbr)])
disp(' ')
disp('DONE')

%% CORREZIONE CON FDR (Benjamini & Hochberg)

% If 'pdep,' the original Bejnamini & Hochberg
% FDR procedure is used, which is guaranteed to be accurate if
% the individual tests are independent or positively dependent
% (e.g., Gaussian variables that are positively correlated or
% independent).  If 'dep,' the FDR procedure
% described in Benjamini & Yekutieli (2001) that is guaranteed
% to be accurate for any test dependency structure (e.g.,
% Gaussian variables with any covariance matrix) is used. 'dep'
% is always appropriate to use but is less powerful than 'pdep.'
% {default: 'pdep'}

disp('===================================================================')

disp('FDR Benjamini & Hochberg')

h_hbo_by = [];
adj_p = [];

% for i = 1:size(p_hbo,1)
%     [h_hbo_by(i,:) , crit_p, adj_ci_cvrg, adj_p(i,:)]=fdr_bh(p_hbo(i,:),0.05,'pdep','no'); %fdr on the metrics
% end

 % for i = 1:size(p_hbo,2)
 %     [h_hbo_by(:,i)  , crit_p, adj_ci_cvrg, adj_p(:,i)]=fdr_bh(p_hbo(:,i),0.05,'pdep','no'); %fdr on the ROI
 % end


[h_hbo_by(:,:) , crit_p, adj_ci_cvrg, adj_p(:,:)]=fdr_bh(p_hbo(:,:),0.05,'pdep','no'); %fdr on the ROI and metrics

hbo_fdr_p= adj_p;

disp('HbO --> FDR Benjamini & Hochberg')
[riga_hbo,colonna_hbo] = find(h_hbo_by==1);
for jj = 1:length(riga_hbo)
    text = '%s with %-7s \t --> r = %f \t p-value_AGGIUSTATO = %f   --> SPEARMAN \n';
    fprintf(text,rsFC_metrics_label(1,riga_hbo(jj)),cgm_metrics_label(1,colonna_hbo(jj)),r_hbo(riga_hbo(jj),colonna_hbo(jj)),adj_p(riga_hbo(jj),colonna_hbo(jj)))
end

h_hbr_by = [];
adj_p = [];


% for i = 1:size(p_hbr,1)
%     [h_hbr_by(i,:), crit_p, adj_ci_cvrg, adj_p(i,:)]=fdr_bh(p_hbr(i,:),0.05,'pdep','no'); %fdr just on the metrics(8)
% end

% for i = 1:size(p_hbr,2)
%     [h_hbr_by(:,i) , crit_p, adj_ci_cvrg, adj_p(:,i) ]=fdr_bh(p_hbr(:,i),0.05,'pdep','no'); %fdr on the ROI(36)
% end

[h_hbr_by(:,:) , crit_p, adj_ci_cvrg, adj_p(:,:)]=fdr_bh(p_hbr(:,:),0.05,'pdep','no'); %fdr on the ROI and metrics

hbr_fdr_p= adj_p;


disp('HbR --> FDR Benjamini & Hochberg')
[riga_hbr,colonna_hbr] = find(h_hbr_by==1);
for jj = 1:length(riga_hbr)
    text = '%s with %-7s \t --> r = %f \t p-value_AGGIUSTATO = %f   --> SPEARMAN \n';
    fprintf(text,rsFC_metrics_label(1,riga_hbr(jj)),cgm_metrics_label(1,colonna_hbr(jj)),r_hbr(riga_hbr(jj),colonna_hbr(jj)),adj_p(riga_hbr(jj),colonna_hbr(jj)))
end

disp('===================================================================')
disp('FDR Benjamini & Hochberg')
disp(['HbO SIGNIFICATIVE --> ',num2str(length(riga_hbo)),'/',num2str(count_tot_hbo)])
disp(['HbR SIGNIFICATIVE --> ',num2str(length(riga_hbr)),'/',num2str(count_tot_hbr)])

%% BONFERRONI HOLM

disp('===================================================================')
disp('BONFERRONI HOLM')

%HbO
for i = 1:size(p_hbo,1)
    [corrected_p_hbo(i,:), h_hbo_bonfh(i,:)]=bonf_holm(p_hbo(i,:),0.05);
end
disp('HbO --> BONFERRONI HOLM')
[riga_hbo,colonna_hbo] = find(h_hbo_bonfh==1);
for jj = 1:length(riga_hbo)
    text = '%s with %-7s \t --> r = %f \t p-value_CORRETTO = %f   --> SPEARMAN \n';
    fprintf(text,rsFC_metrics_label(1,riga_hbo(jj)),cgm_metrics_label(1,colonna_hbo(jj)),r_hbo(riga_hbo(jj),colonna_hbo(jj)),corrected_p_hbO(riga_hbo(jj),colonna_hbo(jj)))
end

% HbR
for i = 1:size(p_hbo,1)
    [corrected_p_hbr(i,:), h_hbr_bonfh(i,:)]=bonf_holm(p_hbr(i,:),0.05);
end
disp('HbR --> BONFERRONI HOLM')
[riga_hbr,colonna_hbr] = find(h_hbr_bonfh==1);
for jj = 1:length(riga_hbr)
    text = '%s with %-7s \t --> r = %f \t p-value_CORRETTO = %f   --> SPEARMAN \n';
    fprintf(text,rsFC_metrics_label(1,riga_hbr(jj)),cgm_metrics_label(1,colonna_hbr(jj)),r_hbr(riga_hbr(jj),colonna_hbr(jj)),corrected_p_hbr(riga_hbr(jj),colonna_hbr(jj)))
end

disp('===================================================================')
disp('BONFERRONI HOLM')
disp(['HbO SIGNIFICATIVE --> ',num2str(length(riga_hbo)),'/',num2str(count_tot_hbo)])
disp(['HbR SIGNIFICATIVE --> ',num2str(length(riga_hbr)),'/',num2str(count_tot_hbr)])

%% FDR MATLAB

disp('===================================================================')
disp('FDR MATLAB')
disp(' ')

FDR_hbo = [];
q_hbo = [];
for i = 1:1:size(p_hbo,1)
    % [FDR_hbo(i,:),q_hbo(i,:)] = mafdr(p_hbo(i,:));
    [FDR_hbo(i,:)] = mafdr(p_hbo(i,:),'BHFDR','true');
end

disp('HbO --> FDR MATLAB')
[riga_hbo,colonna_hbo] = find(FDR_hbo<0.05);
% [riga_hbo,colonna_hbo] = find(q_hbo<0.05);
for jj = 1:length(riga_hbo)
    text = '%s with %-7s \t --> r = %f \t p-value_FDR = %f   --> SPEARMAN \n';
    fprintf(text,rsFC_metrics_label(1,riga_hbo(jj)),cgm_metrics_label(1,colonna_hbo(jj)),r_hbo(riga_hbo(jj),colonna_hbo(jj)),FDR_hbo(riga_hbo(jj),colonna_hbo(jj)))
end

FDR_hbr = [];
q_hbr = [];
for i = 1:1:size(p_hbr,1)
    % [FDR_hbr(i,:),q_hbr(i,:)] = mafdr(p_hbr(i,:));
    FDR_hbr(i,:) = mafdr(p_hbr(i,:),'BHFDR','true');
end

[riga_hbr,colonna_hbr] = find(FDR_hbr<0.05);
% [riga_hbr,colonna_hbr] = find(q_hbr<0.05);
disp('HbR --> FDR MATLAB')
for jj = 1:length(riga_hbr)
    text = '%s with %-7s \t --> r = %f \t p-value_FDR = %f   --> SPEARMAN \n';
    fprintf(text,rsFC_metrics_label(1,riga_hbr(jj)),cgm_metrics_label(1,colonna_hbr(jj)),r_hbr(riga_hbr(jj),colonna_hbr(jj)),FDR_hbr(riga_hbr(jj),colonna_hbr(jj)))
end

disp('===================================================================')
disp('FDR MATLAB')
disp(['HbO SIGNIFICATIVE --> ',num2str(length(riga_hbo)),'/',num2str(count_tot_hbo)])
disp(['HbR SIGNIFICATIVE --> ',num2str(length(riga_hbr)),'/',num2str(count_tot_hbr)])
disp(' ')
disp('DONE')

%% SINGOLA METRICA - SINGOLA rsFC (RAGGRUPPO SECONDO METRICHE)

disp('===================================================================')
disp('HbO')
disp(' ')

count_hbo = 0;
count_tot_hbo = 0;
for i = 1:size(cgm_metrics_doub_z,2)

    for j = 1:size(rsFC_hbo_metrics,2)

        count_tot_hbo = count_tot_hbo + 1;
    
        tmp_rsFC_hbo_metrics = rsFC_hbo_metrics(:,j);
        tmp_cgm_metrics_doub_z = cgm_metrics_doub_z(:,i);
    
        if sum(isnan(tmp_rsFC_hbo_metrics))>0
            to_remove = find(isnan(tmp_rsFC_hbo_metrics));
            tmp_rsFC_hbo_metrics(to_remove) = [];
            tmp_cgm_metrics_doub_z(to_remove,:) = [];
        end
    
        [r_hbo(i,j),p_hbo(i,j)] = corr(tmp_rsFC_hbo_metrics,tmp_cgm_metrics_doub_z,'type','Spearman');
    
        if p_hbo(i,j)<alfa;
            count_hbo = count_hbo+1;
            text = '%-10s with %s \t --> r = %f \t p-value = %f   --> SPEARMAN \n';
            fprintf(text,cgm_metrics_label(1,i),rsFC_metrics_label(1,j),r_hbo(i,j),p_hbo(i,j))
        end
    end
end


disp('===================================================================')
disp('HbR')
disp(' ')

count_hbr = 0;
count_tot_hbr = 0;

for i = 1:size(cgm_metrics_doub_z,2)

    for j = 1:size(rsFC_hbr_metrics,2)

        count_tot_hbr = count_tot_hbr + 1;
    
        tmp_rsFC_hbr_metrics = rsFC_hbr_metrics(:,j);
        tmp_cgm_metrics_doub_z = cgm_metrics_doub_z(:,i);
    
        if sum(isnan(tmp_rsFC_hbr_metrics))>0
            to_remove = find(isnan(tmp_rsFC_hbr_metrics));
            tmp_rsFC_hbr_metrics(to_remove) = [];
            tmp_cgm_metrics_doub_z(to_remove,:) = [];
        end
    
        [r_hbr(i,j),p_hbr(i,j)] = corr(tmp_rsFC_hbr_metrics,tmp_cgm_metrics_doub_z,'type','Spearman');
    
        if p_hbr(i,j)<alfa;
            count_hbr = count_hbr+1;
            text = '%-10s with %s \t --> r = %f \t p-value = %f   --> SPEARMAN \n';
            fprintf(text,cgm_metrics_label(1,i),rsFC_metrics_label(1,j),r_hbr(i,j),p_hbr(i,j))
        end
    end
end


count_tot_hbo = size(rsFC_hbo_metrics,2)*size(cgm_metrics_doub_z,2);
count_tot_hbr = size(rsFC_hbo_metrics,2)*size(cgm_metrics_doub_z,2);

disp('===================================================================')
disp(['ALFA = ',num2str(alfa)])
disp(['HbO SIGNIFICATIVE --> ',num2str(count_hbo),'/',num2str(count_tot_hbo)])
disp(['HbR SIGNIFICATIVE --> ',num2str(count_hbr),'/',num2str(count_tot_hbr)])
disp(' ')
disp('DONE')

%% plotting figures for paper


%section 1 - BGC metrics
cgm_metrics

figure()
spider_plot( ([cgm_metrics_doub(1,:) ; cgm_metrics_doub(2,:)] )  )


figure()
spider_plot( ([cgm_metrics_doub(:,:)] )  )


%sd=2 col
%f3-p3= col6
figure()
plot(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,6)),2),  rsFC_hbr_metrics(~isnan(rsFC_hbr_metrics(:,6)),6),'bx','LineWidth',8)
xlim([min(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,6)),2))  max(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,6)),2))]);
xlabel("SD (mg/dL)")
ylabel("\Delta rsFC")
title("HbR F3-P3 SD   R= "+num2str(r_hbr(6,2))+" fdr p = "+num2str(hbr_fdr_p(6,2))+"")
%r_hbr(6,2)
p = polyfit(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,6)),2),rsFC_hbr_metrics(~isnan(rsFC_hbr_metrics(:,6)),6)   ,1);
hold on
plot(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,6)),2), [p(2) +  (p(1)*cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,6)),2))],'k-','LineWidth',8)


%fz-f4 = 9
%max =4
figure()
plot(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,9)),4),  rsFC_hbr_metrics(~isnan(rsFC_hbr_metrics(:,9)),9),'bx','LineWidth',8)
xlim([min(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,9)),4))  max(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,9)),4))]);
xlabel("Max (mg/dL)")
ylabel("\Delta rsFC")
title("HbR Fz-F4 Max   R = "+num2str(r_hbr(9,4))+" fdr p = "+num2str(hbr_fdr_p(9,4))+"")
%r_hbr(9,4)
%hbr_fdr_p(9,4)
p = polyfit(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,9)),4),rsFC_hbr_metrics(~isnan(rsFC_hbr_metrics(:,9)),9)   ,1);
hold on
plot(cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,9)),4), [p(2) +  (p(1)*cgm_metrics_doub(~isnan(rsFC_hbr_metrics(:,9)),4))],'k-','LineWidth',8)


figure()
imagesc([0;r_hbr(6,2);r_hbr(9,4)])
colorbar();
clim([-1 1])
colormap cool



