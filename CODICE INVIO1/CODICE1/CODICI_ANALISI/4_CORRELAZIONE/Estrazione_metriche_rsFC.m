% ESTRAZIONE METRICHE rsFC

% From some CGM metrics and the difference in the rsFC it is computed the
% canonical correlation analysis in order to find possible correlation
% between glucose and brain activity

close all
clear all
clc

disp_f = 4;

% rng(100)
% rng("default")

%% SETTING OF FOLDER PATH 

base_path = pwd;
data_path_cgm = fullfile(base_path,'DATI','CGM_INTERP');
rsFC_results  = fullfile("RISULTATI DEFINITIVI/rsFC/");
cgm_results_path = fullfile(base_path,'RISULTATI DEFINITIVI/CGM_METRICS');
function_path = fullfile(base_path,'CODICE','FUNCTION');
dati_ext_disc = fullfile('E:\DATI');

addpath(genpath(pwd))
addpath(genpath(data_path_cgm))
addpath(genpath(rsFC_results))
addpath(genpath(function_path))
addpath(genpath(dati_ext_disc))
addpath(genpath(cgm_results_path))

%% Load of rsFC results

% PD : [4;5;8;9;11;14;15;19;47]
% For each PD load of the rsFC results (first and last euglycemia), z-score
% them, perform the difference on the upper triangular matrix and save them
% into the first CCA variable.

%idx_PD = [4;5;8;9;11;14;15;19;47];
%idx_PD = 8;

idx_PD = [3;4;5;8;9;10;11;14;15;19;25;47];

name_first = 'PD_%d_FIRST_rsFC.mat';
name_last  = 'PD_%d_LAST_rsFC.mat';

rsFC_hbo_metrics = [];
rsFC_hbr_metrics = [];

for i = 1:1:length(idx_PD)
    
    curr_idx = idx_PD(i);

    % First euglycemia rsFC
    tmp_first_name = sprintf(name_first, curr_idx);
    load(tmp_first_name)
    disp(['LOAD : ',tmp_first_name])
    disp(' ')
    tmp_hbo_zFC_first = atanh(FC_struct.hbo_FC);
    tmp_hbr_zFC_first = atanh(FC_struct.hbr_FC);
    % tmp_hbo_zFC_first = zscore(FC_struct.hbo_FC);
    % tmp_hbr_zFC_first = zscore(FC_struct.hbr_FC);

    % Last Euglycemia interval
    tmp_last_name = sprintf(name_last, curr_idx);
    load(tmp_last_name)
    disp(['LOAD : ',tmp_last_name])
    disp(' ')
    tmp_hbo_zFC_last = atanh(FC_struct.hbo_FC);
    tmp_hbr_zFC_last = atanh(FC_struct.hbr_FC);
    % tmp_hbo_zFC_last = zscore(FC_struct.hbo_FC);
    % tmp_hbr_zFC_last = zscore(FC_struct.hbr_FC);

    % Compute the difference between the last and first rsFC Z-Score matrix
    tmp_hbo_diff_rsFC = tmp_hbo_zFC_last-tmp_hbo_zFC_first;
    tmp_hbr_diff_rsFC = tmp_hbr_zFC_last-tmp_hbr_zFC_first;
    
    if disp_f == 3
        figure()
        
        imagesc(tmp_hbo_diff_rsFC)
        axis square
        colormap jet
        caxis([-1 1])
        colorbar
        title(['PD ',num2str(curr_idx),' - HbO rsFC MATRIX - Z-SCORE'])
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_comp;
        xticks = linspace(1, size(tmp_hbo_zFC_first, 2), numel(xticklabels));
        xtickangle(45)
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
        yticklabels = FC_struct.labels_comp;
        yticks = linspace(1, size(tmp_hbo_zFC_first, 2), numel(yticklabels));
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    end


    % Extraction of the upper triangolar matrix and creation of the rsFC
    % metrics [n_PD x n_conf]
    tmp_up = triu(tmp_hbo_diff_rsFC);
    At = tmp_up.';
    m  = (1:size(At,1)).' >= (1:size(At,2));
    tmp_hbo_diff_up  = At(m);
    tmp_hbo_diff_up = tmp_hbo_diff_up(~isnan(tmp_hbo_diff_up))';
    if length(tmp_hbo_diff_up) < 36
        to_add = NaN*ones(1,36-length(tmp_hbo_diff_up));
        tmp_hbo_diff_up = [to_add tmp_hbo_diff_up]
    end
    rsFC_hbo_metrics = [rsFC_hbo_metrics ; tmp_hbo_diff_up];

    tmp_up = triu(tmp_hbr_diff_rsFC);
    At = tmp_up.';
    m  = (1:size(At,1)).' >= (1:size(At,2));
    tmp_hbr_diff_up  = At(m);
    tmp_hbr_diff_up = tmp_hbr_diff_up(~isnan(tmp_hbr_diff_up))';
    if length(tmp_hbr_diff_up) < 36
        to_add = NaN*ones(1,36-length(tmp_hbr_diff_up));
        tmp_hbr_diff_up = [to_add tmp_hbr_diff_up]
    end
    rsFC_hbr_metrics = [rsFC_hbr_metrics ; tmp_hbr_diff_up];


    % DISPLAY OF THE rsFC Z-SCORE MATRIX
    if disp_f == 1
        % HbO
        figure()
    
        t_hbo = tiledlayout(1,2);   
        t_hbo.TileSpacing = 'compact';
        t_hbo.Padding = 'compact';
        
        nexttile
        imagesc(tmp_hbo_zFC_first)
        axis square
        colormap jet
        caxis([-1 1])
        colorbar
        title(['PD ',num2str(curr_idx),' - HbO rsFC MATRIX - Z-SCORE'])
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_comp;
        xticks = linspace(1, size(tmp_hbo_zFC_first, 2), numel(xticklabels));
        xtickangle(45)
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
        yticklabels = FC_struct.labels_comp;
        yticks = linspace(1, size(tmp_hbo_zFC_first, 2), numel(yticklabels));
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
        
        nexttile
    
        imagesc(FC_struct.hbo_p)
        colormap jet
        axis square
        caxis([0 0.0014])
        colorbar
        title('HbO P-VALUES')
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_short;
        xtickangle(45)
        xticks = linspace(1, size(tmp_hbr_zFC_first, 2), numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        yticklabels = FC_struct.labels_short;
        yticks = linspace(1, size(tmp_hbr_zFC_first, 2), numel(yticklabels));
        % ytickangle(45)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    
        % HbR
        figure()
    
        t_hbo = tiledlayout(1,2);
        t_hbo.TileSpacing = 'compact';
        t_hbo.Padding = 'compact';
        
        nexttile
        imagesc(tmp_hbr_zFC_first)
        axis square
        colormap jet
        caxis([-1 1])
        colorbar
        title(['PD ',num2str(curr_idx),' - HbR rsFC MATRIX - Z-SCORE'])
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_comp;
        % xticklabels = label_active_nodes_sort;
        xticks = linspace(1, size(FC_struct.hbo_p, 2), numel(xticklabels));
        xtickangle(45)
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
        yticklabels = FC_struct.labels_comp;
        % yticklabels = label_active_nodes_sort;
        yticks = linspace(1, size(FC_struct.hbo_FC, 2), numel(yticklabels));
        % ytickangle(45)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
        
        nexttile
    
        imagesc(FC_struct.hbr_p)
        colormap jet
        axis square
        caxis([0 0.0014])
        colorbar
        % title('HbR P-VALUES')
        title('HbR P-VALUES')
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_short;
        xtickangle(45)
        xticks = linspace(1, size(FC_struct.hbo_FC, 2), numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        yticklabels = FC_struct.labels_short;
        yticks = linspace(1, size(FC_struct.hbo_FC, 2), numel(yticklabels));
        % ytickangle(45)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

        pause()
        close all
    end

    if disp_f == 2
        % HbO
        figure()
    
        t_hbo = tiledlayout(1,2);   
        t_hbo.TileSpacing = 'compact';
        t_hbo.Padding = 'compact';
        
        nexttile
        imagesc(tmp_hbo_zFC_last)
        axis square
        colormap jet
        caxis([-1 1])
        colorbar
        title(['PD ',num2str(curr_idx),' - HbO rsFC MATRIX - Z-SCORE'])
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_comp;
        xticks = linspace(1, size(tmp_hbo_zFC_last, 2), numel(xticklabels));
        xtickangle(45)
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
        yticklabels = FC_struct.labels_comp;
        yticks = linspace(1, size(tmp_hbo_zFC_last, 2), numel(yticklabels));
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
        
        nexttile
    
        imagesc(FC_struct.hbo_p)
        colormap jet
        axis square
        caxis([0 0.0014])
        colorbar
        title('HbO P-VALUES')
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_short;
        xtickangle(45)
        xticks = linspace(1, size(tmp_hbr_zFC_last, 2), numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        yticklabels = FC_struct.labels_short;
        yticks = linspace(1, size(tmp_hbr_zFC_last, 2), numel(yticklabels));
        % ytickangle(45)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    
        % HbR
        figure()
    
        t_hbo = tiledlayout(1,2);
        t_hbo.TileSpacing = 'compact';
        t_hbo.Padding = 'compact';
        
        nexttile
        imagesc(tmp_hbr_zFC_last)
        axis square
        colormap jet
        caxis([-1 1])
        colorbar
        title(['PD ',num2str(curr_idx),' - HbR rsFC MATRIX - Z-SCORE'])
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_comp;
        % xticklabels = label_active_nodes_sort;
        xticks = linspace(1, size(FC_struct.hbo_p, 2), numel(xticklabels));
        xtickangle(45)
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
        yticklabels = FC_struct.labels_comp;
        % yticklabels = label_active_nodes_sort;
        yticks = linspace(1, size(FC_struct.hbo_FC, 2), numel(yticklabels));
        % ytickangle(45)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
        
        nexttile
    
        imagesc(FC_struct.hbr_p)
        colormap jet
        axis square
        caxis([0 0.0014])
        colorbar
        % title('HbR P-VALUES')
        title('HbR P-VALUES')
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gca,'fontsize',12,'fontweight','bold')
        
        xticklabels = FC_struct.labels_short;
        xtickangle(45)
        xticks = linspace(1, size(FC_struct.hbo_FC, 2), numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        yticklabels = FC_struct.labels_short;
        yticks = linspace(1, size(FC_struct.hbo_FC, 2), numel(yticklabels));
        % ytickangle(45)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

        pause()
        close all
    end
end

