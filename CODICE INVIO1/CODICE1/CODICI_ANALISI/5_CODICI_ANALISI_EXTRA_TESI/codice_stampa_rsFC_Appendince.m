%% DISPLAY rsFC MATRIX

% For better rappresentation of results use tildelayout (paper/thesis)
close all
clear all
clc

%%
% idx_PD = [4;5;8;9;11;14;15;19;47];
idx_PD = [8];

name_first = 'PD_%d_FIRST_rsFC.mat';
name_first_save = 'PD_%d_FIRST_rsFC.jpg';
name_last  = 'PD_%d_LAST_rsFC.mat';
name_last_save  = 'PD_%d_LAST_rsFC.jpg';

rsFC_hbo_metrics = [];
rsFC_hbr_metrics = [];

%%
for i = 1:1:length(idx_PD)

    curr_idx = idx_PD(i);

    % First euglycemia rsFC
    tmp_first_name = sprintf(name_first, curr_idx);
    tmp_name_disp = sprintf('PD %d - FIRST',curr_idx);
    tmp_first_name_save = sprintf(name_first_save,curr_idx);
    load(tmp_first_name)
    disp(['LOAD : ',tmp_first_name])
    disp(' ')

    hbo_FC_first = FC_struct.hbo_FC;
    hbr_FC_first = FC_struct.hbr_FC;
    label_active_nodes_sort = FC_struct.labels_comp;

    label_active_nodes_short = FC_struct.labels_short;

    % HbO
    fig = figure()
    set(gcf, 'Position', get(0, 'Screensize'));

    f.WindowState = 'maximized';


    t = tiledlayout(2,2);   
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    nexttile
    imagesc(hbo_FC_first)
    axis square
    colormap jet
    caxis([-1 1])
    colorbar
    title('HbO rsFC MATRIX')
    % set(gca,'fontsize',12,'fontweight','bold')


    xticklabels = label_active_nodes_short;
    % xticklabels = label_active_nodes_sort;
    xticks = linspace(1, size(hbo_FC_first, 2), numel(xticklabels));
    xtickangle(45)
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    yticklabels = label_active_nodes_sort;
    % yticklabels = label_active_nodes_sort;
    yticks = linspace(1, size(hbo_FC_first, 2), numel(yticklabels));
    % ytickangle(45)
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

    % HbR

    % figure()

    % t = tiledlayout(1,2);
    % t.TileSpacing = 'compact';
    % t.Padding = 'compact';

    nexttile
    imagesc(hbr_FC_first)
    axis square
    colormap jet
    caxis([-1 1])
    colorbar
    title('HbR rsFC MATRIX')

    % set(gca,'fontsize',12,'fontweight','bold')

    xticklabels = label_active_nodes_short;
    xtickangle(45)
    xticks = linspace(1, size(hbr_FC_first, 2), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    yticklabels = label_active_nodes_short;
    yticks = linspace(1, size(hbr_FC_first, 2), numel(yticklabels));
    % ytickangle(45)
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

   


    % 
    % sgtitle(tmp_name_disp,'fontweight','bold')

    % path = ['C:\Users\giaco\OneDrive\Desktop\TESI MAGISTRALE\RISULTATI DEFINITIVI\rsFC_APPENDICE\',tmp_first_name_save]
    % saveas(fig,path)
    % nome = ['PD ',num2str(curr_idx)];
    % print(nome,'-djpeg','-r600')
end

for i = 1:1:length(idx_PD)

    curr_idx = idx_PD(i);

    % First euglycemia rsFC
    tmp_last_name = sprintf(name_last, curr_idx);
    tmp_name_disp = sprintf('PD %d - LAST',curr_idx);
    tmp_last_name_save = sprintf(name_last_save,curr_idx);
    load(tmp_last_name)
    disp(['LOAD : ',tmp_last_name])
    disp(' ')

    hbo_FC_last = FC_struct.hbo_FC;
    hbo_p_sort = FC_struct.hbo_p;
    p_bonf = FC_struct.p_bonf;
    hbr_FC_last = FC_struct.hbr_FC;
    hbr_p_sort = FC_struct.hbr_p;
    label_active_nodes_sort = FC_struct.labels_comp;

    label_active_nodes_short = FC_struct.labels_short;

    % HbO
    nexttile
    imagesc(hbo_FC_last)
    axis square
    colormap jet
    caxis([-1 1])
    colorbar
    title('HbO rsFC MATRIX')
    % set(gca,'fontsize',12,'fontweight','bold')


    xticklabels = label_active_nodes_short;
    % xticklabels = label_active_nodes_sort;
    xticks = linspace(1, size(hbo_FC_last, 2), numel(xticklabels));
    xtickangle(45)
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    yticklabels = label_active_nodes_sort;
    % yticklabels = label_active_nodes_sort;
    yticks = linspace(1, size(hbo_FC_last, 2), numel(yticklabels));
    % ytickangle(45)
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

    

    % HbR

    % figure()

    % t = tiledlayout(1,2);
    % t.TileSpacing = 'compact';
    % t.Padding = 'compact';

    nexttile
    imagesc(hbr_FC_last)
    axis square
    colormap jet
    caxis([-1 1])
    colorbar
    title('HbR rsFC MATRIX')

    % set(gca,'fontsize',12,'fontweight','bold')

    xticklabels = label_active_nodes_short;
    xtickangle(45)
    xticks = linspace(1, size(hbr_FC_last, 2), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    yticklabels = label_active_nodes_short;
    yticks = linspace(1, size(hbr_FC_last, 2), numel(yticklabels));
    % ytickangle(45)
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

    % sgtitle('PRIMA ED ULTIMA FINESTRA DI EUGLICEMIA PD08')
    print('rsFC_PD8','-djpeg','-r600')

end