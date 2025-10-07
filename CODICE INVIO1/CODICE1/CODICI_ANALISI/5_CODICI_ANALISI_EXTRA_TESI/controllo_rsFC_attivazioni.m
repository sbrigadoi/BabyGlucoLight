% DISPLAY rsFC MATRIX

% For better rappresentation of results use tildelayout (paper/thesis)
close all
clear all
clc

%
idx_PD = [4;5;8;9;11;14;15;19;47];

name_first = 'PD_%d_FIRST_rsFC.mat';
name_last  = 'PD_%d_LAST_rsFC.mat';

soglia = 0.5;

conto_vect_HbO_z = zeros(1,8);
conto_vect_HbR_z = zeros(1,8);

conto_vect_HbO_n = zeros(1,7);
conto_vect_HbR_n = zeros(1,7);

% FIRST

for i = 1:1:length(idx_PD)

    curr_idx = idx_PD(i);

    % First euglycemia rsFC
    tmp_first_name = sprintf(name_first, curr_idx);
    tmp_name_disp = sprintf('PD %d - FIRST',curr_idx);
    load(tmp_first_name)
    disp(['LOAD : ',tmp_first_name])
    disp(' ')

    hbo_FC = FC_struct.hbo_FC;
    hbo_p_sort = FC_struct.hbo_p;
    p_bonf = FC_struct.p_bonf;
    hbr_FC = FC_struct.hbr_FC;
    hbr_p_sort = FC_struct.hbr_p;
    label_active_nodes_sort = FC_struct.labels_comp;

    sotto_diago_hbo = diag(hbo_FC,-1);
    sotto_diago_hbo_2 = diag(hbo_FC,-2);
    sotto_diago_hbr = diag(hbr_FC,-1);
    sotto_diago_hbr_2 = diag(hbo_FC,-2);

    for jj = 1:1:length(sotto_diago_hbo)
        if sotto_diago_hbo(jj)>soglia
            conto_vect_HbO_z(1,jj) = conto_vect_HbO_z(1,jj)+1;
        end

        if sotto_diago_hbr(jj)>soglia
            conto_vect_HbR_z(1,jj) = conto_vect_HbR_z(1,jj)+1;
        end

    end

    for jj = 1:1:length(sotto_diago_hbr_2)
        if sotto_diago_hbo_2(jj)>soglia
            conto_vect_HbO_n(1,jj) = conto_vect_HbO_n(1,jj)+1;
        end

        if sotto_diago_hbr_2(jj)>soglia
            conto_vect_HbR_n(1,jj) = conto_vect_HbR_n(1,jj)+1;
        end
    end    

end


for i = 1:1:length(idx_PD)

    curr_idx = idx_PD(i);

    % First euglycemia rsFC
    tmp_last_name = sprintf(name_last, curr_idx);
    tmp_name_disp = sprintf('PD %d - FIRST',curr_idx);
    load(tmp_last_name)
    disp(['LOAD : ',tmp_last_name])
    disp(' ')

    hbo_FC = FC_struct.hbo_FC;
    hbo_p_sort = FC_struct.hbo_p;
    p_bonf = FC_struct.p_bonf;
    hbr_FC = FC_struct.hbr_FC;
    hbr_p_sort = FC_struct.hbr_p;
    label_active_nodes_sort = FC_struct.labels_comp;

    sotto_diago_hbo = diag(hbo_FC,-1);
    sotto_diago_hbo_2 = diag(hbo_FC,-2);
    sotto_diago_hbr = diag(hbr_FC,-1);
    sotto_diago_hbr_2 = diag(hbr_FC,-2);

    for jj = 1:1:length(sotto_diago_hbo)
        if sotto_diago_hbo(jj)>soglia
            conto_vect_HbO_z(jj) = conto_vect_HbO_z(jj)+1;
        end

        if sotto_diago_hbr(jj)>soglia
            conto_vect_HbR_z(jj) = conto_vect_HbR_z(jj)+1;
        end

    end
    for jj = 1:1:length(sotto_diago_hbr_2)
        if sotto_diago_hbo_2(jj)> soglia
            conto_vect_HbO_n(1,jj) = conto_vect_HbO_n(1,jj)+1;
        end

        if sotto_diago_hbr_2(jj)>soglia
            conto_vect_HbR_n(1,jj) = conto_vect_HbR_n(1,jj)+1;
        end
    end
end


%%
conteggi_totali_centrali = [conto_vect_HbO_z;conto_vect_HbR_z];

conteggi_totali_centrali(:,[3 6]) = [];

conteggi_totali_lat = [conto_vect_HbO_n([1 4 7]);conto_vect_HbR_n([1 4 7])];



%%
label = {'HbO';'HbR'};
tab_finale_centrale = table(label,conteggi_totali_centrali(:,1),conteggi_totali_centrali(:,2),conteggi_totali_centrali(:,3),conteggi_totali_centrali(:,4),conteggi_totali_centrali(:,5),conteggi_totali_centrali(:,6));
tab_finale_centrale.Properties.VariableNames = {'Name','Fz-F3','Fz-F4','Cz-C3','Cz-C4','Pz-P3','Pz-P4'};

tab_finale_lat = table(label,conteggi_totali_lat(:,1),conteggi_totali_lat(:,2),conteggi_totali_lat(:,3));
tab_finale_lat.Properties.VariableNames = {'Name','F3-F4','C3-C4','P3-P4'};
 
tab_finale = table(label,conteggi_totali_centrali(:,1),conteggi_totali_centrali(:,2),conteggi_totali_lat(:,1),conteggi_totali_centrali(:,3),conteggi_totali_centrali(:,4),conteggi_totali_lat(:,2),conteggi_totali_centrali(:,5),conteggi_totali_centrali(:,6),conteggi_totali_lat(:,3));
tab_finale.Properties.VariableNames = {'Name','Fz-F3','Fz-F4','F3-F4','Cz-C3','Cz-C4','C3-C4','Pz-P3','Pz-P4','P3-P4'};

%%

mat = table2array(tab_finale(:,2:end));
Y = mat(1,:);
YY = mat(2,:);


figure()

set(gcf, 'Position', get(0, 'Screensize'));

subplot(2,1,1)
b = bar(mat(1,:),'FaceColor','flat')
xticklabels(tab_finale.Properties.VariableNames(2:end))
for i = 1:1:9
    b.CData(i,:) = [0.6350 0.0780 0.1840];
end
ylim([0 15])
title('CORRELAZIONI INTRALOBO CON r > 0.5 - HbO')
ylabel('Conteggio')
h=gca; h.XAxis.TickLength = [0 0];
set(gca,'FontWeight','bold')
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center','fontweight','bold');
legend('HbO')


subplot(2,1,2)
b = bar(mat(2,:))
xticklabels(tab_finale.Properties.VariableNames(2:end))
title('CORRELAZIONI INTRALOBO CON r > 0.5 - HbR')
ylabel('Conteggio')
h=gca; h.XAxis.TickLength = [0 0];
ylim([0 15])
set(gca,'FontWeight','bold')
text(1:length(YY),YY,num2str(YY'),'vert','bottom','horiz','center','fontweight','bold');
legend('HbR')
% 
% saveas(gcf,'CORR_INTRALOBO.jpg')
% print('FIG_49_BAR_conteggio_r_mag_','-djpeg','-r600')
%%

% figure()
% bar(mat')
% xticklabels(tab_finale.Properties.VariableNames(2:end))
% legend('HbO','HbR')

%%