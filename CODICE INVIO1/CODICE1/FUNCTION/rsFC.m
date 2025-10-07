function [ROI_table,ROI_active_table,tenTwenty,FC_struct] = rsFC(tenFive,gmSurfaceMesh,vol2gm,all_J,C_thresh,r_ROI,hbo,hbr,f_res)
% Perform the resting state functional connectivity

% INPUT: tenFive       --> structure with labels and positions of 10/5 EEG system
%        gmSurfaceMesh --> Age matched GM surface
%        vol2gm        --> Age matched transfer matrix (volumetric to GM)
%        all_J         --> Summation over all rows and coloumn of a J matrix
%        C_thresh      --> Coverage threshold (*)
%        r_ROI         --> radius of sferical ROIs
%        hbo           --> 5' hbo signal for each nodes (samples x nodes)
%        hbr           --> 5' hbr signal for each nodes (samples x nodes)
%        f_res         --> 1 = display results

% OUTPUT: ROI_table        --> informative table with center ROIs
%         ROI_table_active --> informative table with active center ROIs
%         tenTwenty        --> Labels and position of 10/20 EEG system
%         FC_struct        --> structure with rsFC, p_value, p_bonf e ROI labels
%                              for HbO and HbR

% (*) = 'Array Designer: automated optimized array design for functional
% near-infrared spectroscopy' - Brigadoi et. al.
%% IDENTIFICATION OF 10/20 EEG SYSTEM

% Extraction of the labels and positions of the 10/20 EEG system sources
% from the 10/5 EEG system.

tenTwenty = struct();
tenTwenty.labels = ["Cz", "Fpz", "Fz", "Pz", "Oz", "T7", "C3", "C4", "T8", "Fp1", "F7", "P7", "O1", "Fp2", "F8", "P8", "O2", "F3", "F4", "P3","P4"]';
tenTwenty.positions = [];

for i = 1:1:length(tenTwenty.labels)
    
    tmp_idx = find(tenTwenty.labels(i) == tenFive.labels);
    tenTwenty.positions(i,:) = tenFive.positions(tmp_idx,:);
    
end

%% PROJECTION OF SCALP SOURCES TO GM NODES 

% Firstly, we need to make sure the co-ordinates above match to a node in 
% the mesh. 

probes = tenTwenty.positions;
for i=1:size(probes,1)
    % Find the distance between the probe positions and every node in the mesh
    dist_probe_node = pdist2(probes(i,:),gmSurfaceMesh.node(:,1:3)); 
    % Find the minimum of these distances. 
    % I.e for each probe, we find the closest node in the mesh.
    find(dist_probe_node == min(dist_probe_node)); 
    % Set each probe position to the nearest mesh node
    probes(i,1:3) = gmSurfaceMesh.node(find(dist_probe_node == min(dist_probe_node)),1:3); 
end

%% ROIs CENTER SELECTION BASED ON COVERAGE THRESHOLD

% J matrix mapped over GM
J_gm_conv = vol2gm*all_J';

n_nodes = length(J_gm_conv);

% If a probe (node) is covered from the array-sensitivity 
% (sensitivity > coverage threshold) it is flagged as active and taken into
% the account fot the rsFC as ROI center.
f_active_probe = [];

for i = 1:1:size(probes,1)
    tmp_idx = find(probes(i,1) == gmSurfaceMesh.node(:,1) & probes(i,2) == gmSurfaceMesh.node(:,2) & probes(i,3) == gmSurfaceMesh.node(:,3));
    tmp_J = abs(J_gm_conv(tmp_idx));
    if tmp_J > C_thresh
        f_active_probe(i,1) = 1;
    else
        f_active_probe(i,1) = 0;
    end
end

%% Informative table for all nodes

ROI_table = table();
ROI_table.labels = tenTwenty.labels;
ROI_table.coordx  = probes(:,1);
ROI_table.coordy  = probes(:,2);
ROI_table.coordz  = probes(:,3);   
% ROI_table.active  = f_active_probe;
ROI_table.active = ones(21,1);

%% TRATTENGO SOLO LE ROI COMUNI (9)

standard_ROI_label = ["F3","Fz","F4","C3","Cz","C4","P3","Pz","P4"]';
n_roi_active = length(ROI_table.active);
to_delete = [];
 
for i = 1:1:n_roi_active
    tmp_name = ROI_table.labels(i);
    f_find = find(tmp_name == standard_ROI_label);
    if isempty(f_find) & ROI_table.active(i)==1
        ROI_table.active(i,:) = 0;
        disp(['ROI ',tmp_name,' REMOVED'])
    end
end


%% Informative table with only good nodes
ROI_active_table = ROI_table;
idx_to_remove = find(ROI_active_table.active==0);
ROI_active_table(idx_to_remove,:) = [];


%% IDENTIFICATION OF ROIs

% Distance between ROIs center

probes_active = probes(ROI_table.active==1,:);
idx_active_p = find(ROI_table.active==1);
label_active_nodes = ROI_table.labels(idx_active_p);             

% dist_active = pdist(probes_active,'euclidean');

% For each active node (ROIs center) it is computed the distance from each
% nodes inside the GM surface mesh. Each node with a distance below the ROI
% radius is assigned to that ROI.

% Initialization of support variables
n_active_probe = size(probes_active,1);
gmSurfaceMesh.node(:,4) = zeros(size(gmSurfaceMesh.node,1),1);

idx = find(ROI_table.active==1);

for i = 1:1:n_active_probe

    % Distance between ROIs center and each nodes in the GM surface mesh
    tmp_d = pdist2(probes_active(i,:),gmSurfaceMesh.node(:,1:3));
    % Nodes to ROIs
    tmp_idx = find(tmp_d<=r_ROI);
    gmSurfaceMesh.node(tmp_idx,4) = i; 
 
end

% n_nodes_ROI = [];
% 
% for i = 1:1:size(probes,1)
% 
%     aa = find(gmSurfaceMesh.node(:,4) == i);
%     n_nodes_ROI = [n_nodes_ROI ; length(aa)];
% 
% end
% 
% n_nodes_ROI
%% ROIs INFORMATION

% Count of active nodes for each ROIs
ROI_vect = gmSurfaceMesh.node(:,4);
ROI_count = [0:1:n_active_probe]';

pos = 0;
for i = 0:1:n_active_probe
    pos = pos+1;
    tmp_idx = find(ROI_vect==i);
    tmp_n = length(tmp_idx);
    ROI_count(pos,2) = tmp_n;
end

% Usefull information about ROIs
n_roi = length(find(ROI_count(:,2)>0))-1; % ROI 0 (inactive node) excluded
n_not_roi = ROI_count(1,2);
ROI_sort = unique(sort(ROI_count(:,2)));
min_ROI_active = ROI_sort(2);
idx_min_ROI_active = find(ROI_count(:,2)==min_ROI_active)-1;
max_ROI_active = ROI_sort(end-1);
idx_max_ROI_active = find(ROI_count(:,2)==max_ROI_active)-1;

ROI_count_table = table();
ROI_count_table.name = ROI_active_table.labels;
ROI_count_table.nodes_total = ROI_count(2:end,2);

% Numero di nodi attivi per ogni roi

gmSurfaceMesh.node(:,5) = 5;

for i = 1:1:n_roi
    tmp_idx_node = find(gmSurfaceMesh.node(:,4)==i);
    for j = 1:1:size(tmp_idx_node,1)
        curr_idx = tmp_idx_node(j);
        tmp_J = abs(J_gm_conv(curr_idx));
        if tmp_J > C_thresh
            gmSurfaceMesh.node(curr_idx,5) = 1;
        else
            gmSurfaceMesh.node(curr_idx,5) = 0;
        end
    end
end

% Conto nodi attivi per ogni ROI
ROI_count_table.node_active = zeros(n_roi,1);
for i = 1:n_roi
    tmp_n_active = find(gmSurfaceMesh.node(:,4) == i & gmSurfaceMesh.node(:,5)==1);
    ROI_count_table.node_active(i) = length(tmp_n_active);
end

% Se il numero dei nodi attivi è al di sotto di un 1/3 dei nodi della ROI,
% tale ROI non viene presa in considerazione per l'analisi rsFC poichè
% risulterebbe poco coperta dalla sensitività dell'array.
ROI_count_table.ROI_active = 5*ones(n_roi,1);
gmSurfaceProva = gmSurfaceMesh;
gmSurfaceProva.node(:,5) = [];
for i = 1:1:n_roi
    if ROI_count_table.node_active(i)<round(ROI_count_table.nodes_total(i)/3)
        ROI_count_table.ROI_active(i) = 0;
        tmp_idx = find(ROI_count_table.name(i)==ROI_table.labels);
        tmp_idx_2 = find(ROI_count_table.name(i) == ROI_active_table.labels);
        ROI_table.active(tmp_idx) = 0;
        ROI_active_table.active(i) = 0;
        % tmp_remove = find(gmSurfaceMesh.node(:,4)==i);
        % gmSurfaceProva.node(tmp_remove,4) = 0;
    else
        ROI_count_table.ROI_active(i) = 1;
    end
end



if length(find(ROI_count_table.ROI_active == 5))>0 
    disp('ERROR IN THE CODE')
end

if length(find(ROI_count_table.ROI_active == 0))>0 
    disp('ERROR: Some ROIs are not suffciently covered from the sensor array')
end


%% RESTING STATE FUNCTIONAL CONNECTIVITY

% Mean HbO and HbR signal for each ROIs (samples x ROIs)
hbo_mean_roi = [];
hbr_mean_roi = [];

for i = 1:1:n_roi
    % For each ROI, it is extraced the indexes of each nodes inside it
    tmp_idx_node_roi = find(gmSurfaceMesh.node(:,4)==i);
    % hbo.gm and hbr.gm contains the signal for each nodes(samples x nodes)
    tmp_hbo = hbo.gm(:,tmp_idx_node_roi);
    tmp_hbr = hbr.gm(:,tmp_idx_node_roi);
    % Mean over rows to obtain the mean HbO and HbR for each time points
    % for each ROIs --> (samples x ROI)
    hbo_mean_roi(:,i) = mean(tmp_hbo,2);
    hbr_mean_roi(:,i) = mean(tmp_hbr,2);
end


% rsFC Computation : Pearson correlation

% Compute the Pearson's correlation with relative p_value for each pairs of
% ROIs HbO and HbR time series. The command [corr,p_value] = corr(X)
% returns the Pearson's correlation coeficient computes over COLOUMNS. 

[FC_hbo, p_hbo] = corr(hbo_mean_roi);  
[FC_hbr, p_hbr] = corr(hbr_mean_roi);  

% Bonferroni correction (same for both HbO and HbR)
n_comparison = (size(FC_hbo,2)*size(FC_hbo,2)-size(FC_hbo,2))/2;
p_bonf = 0.05/n_comparison;

% Sort the ROIs in the following order: left center rigth, frontal to
% parietal

hbo_FC_sort = FC_hbo([6,2,7,4,1,5,8,3,9],[6,2,7,4,1,5,8,3,9]);
hbr_FC_sort = FC_hbr([6,2,7,4,1,5,8,3,9],[6,2,7,4,1,5,8,3,9]);

hbo_p_sort = p_hbo([6,2,7,4,1,5,8,3,9],[6,2,7,4,1,5,8,3,9]);
hbr_p_sort = p_hbr([6,2,7,4,1,5,8,3,9],[6,2,7,4,1,5,8,3,9]);


label_active_nodes_sort = label_active_nodes([6,2,7,4,1,5,8,3,9]);
label_active_nodes_sort_complete = ["F3 (Frontal left)", "Fz (Frontal center)", "F4 (Frontal right)", "C3 (Motor left)","Cz (Motor center)","C4 (Motor right)","P3 (Parietal left)","Pz (Parietal center)","P4 (Parietal right)"]';

%% METTO NAN PER ESCLUDERE ROI CON NUMERO NODI ATTIVI MINORI SOGLIA (PD14)
for i = 1:1:length(ROI_active_table.active)
    if ROI_active_table.active(i) == 0
        tmp_idx = find(ROI_active_table.labels(i)==label_active_nodes_sort);
        % Colonna
        hbo_FC_sort(:,tmp_idx) = NaN;
        hbr_FC_sort(:,tmp_idx) = NaN;
        hbo_p_sort(:,tmp_idx)  = NaN;
        hbr_p_sort(:,tmp_idx)  = NaN;
        % Riga
        hbo_FC_sort(tmp_idx,:) = NaN;
        hbr_FC_sort(tmp_idx,:) = NaN;
        hbo_p_sort(tmp_idx,:)  = NaN;
        hbr_p_sort(tmp_idx,:)  = NaN;
    end
end

%% FUNCTION OUTPUT

% Initialization of output variables

FC_struct        = struct();
FC_struct.hbo_FC = hbo_FC_sort;
FC_struct.hbo_p  = hbo_p_sort;
FC_struct.hbr_FC = hbr_FC_sort;
FC_struct.hbr_p  = hbr_p_sort;
FC_struct.p_bonf = p_bonf;

FC_struct.labels_comp  = label_active_nodes_sort_complete;
FC_struct.labels_short = label_active_nodes_sort;


%% DISPLAY OF RESULTS

if f_res == 1
    %% PLOT ALL ROIs CENTER 

    % Plot ROIs center in 3D spaces
    % figure()
    % plot3(probes(:,1),probes(:,2),probes(:,3),'ro')
    % xlabel('x / mm')
    % ylabel('y / mm')
    % zlabel('z / mm')
    % 
    % 
    % % Plot ROIs center in 2D
    % figure()
    % plot(probes(:,1),probes(:,2),'ro')
    % xlabel('x / mm')
    % ylabel('y / mm')


    %% ACTIVE ROIs CENTER DISTRIBUTION ON GM (red = atvice, black = removed)

    warning('off')
    rmpath(genpath('C:\Users\giaco\OneDrive\Desktop\TESI MAGISTRALE\CODICE\FUNCTION\NIRFAST-9.1'));
    rmpath(genpath('C:\Users\giaco\OneDrive\Desktop\TESI MAGISTRALE\CODICE\FUNCTION\NIRFASTer-2.0'));
    warning('on')

    figure()

    set(gcf, 'Position', get(0, 'Screensize'));
    t = tiledlayout(2,2);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    nexttile
    plotmeshISO(gmSurfaceMesh.node(:,1:3),gmSurfaceMesh.face)
    view(180,0)
    hold on
    set(gcf, 'Position', get(0, 'Screensize'));
    for i = 1:1:length(f_active_probe)
        if f_active_probe(i) == 1
            h(1)=plot3(probes(i,1),probes(i,2),probes(i,3),'.r','MarkerSize',30,'DisplayName','Active ROI')
        else
            h(2)=plot3(probes(i,1),probes(i,2),probes(i,3),'.k','MarkerSize',30,'DisplayName','Deactive ROI')
        end
    end
    legend(h(1:2),'Active ROI','Inactive ROI')
    title('Anterior')

    nexttile
    plotmeshISO(gmSurfaceMesh.node(:,1:3),gmSurfaceMesh.face)
    view(0,0)
    hold on
    set(gcf, 'Position', get(0, 'Screensize'));
    for i = 1:1:length(f_active_probe)
        if f_active_probe(i) == 1
            h(1)=plot3(probes(i,1),probes(i,2),probes(i,3),'.r','MarkerSize',30,'DisplayName','Active ROI')
        else
            h(2)=plot3(probes(i,1),probes(i,2),probes(i,3),'.k','MarkerSize',30,'DisplayName','Deactive ROI')
        end
    end
    legend(h(1:2),'Active ROI','Inactive ROI')
    title('Posterior')

    nexttile
    plotmeshISO(gmSurfaceMesh.node(:,1:3),gmSurfaceMesh.face)
    view(90,0)
    hold on
    set(gcf, 'Position', get(0, 'Screensize'));
    for i = 1:1:length(f_active_probe)
        if f_active_probe(i) == 1
            h(1)=plot3(probes(i,1),probes(i,2),probes(i,3),'.r','MarkerSize',30,'DisplayName','Active ROI')
        else
            h(2)=plot3(probes(i,1),probes(i,2),probes(i,3),'.k','MarkerSize',30,'DisplayName','Deactive ROI')
        end
    end
    legend(h(1:2),'Active ROI','Inactive ROI')  
    title('Lateral')

    nexttile()
    plotmeshISO(gmSurfaceMesh.node(:,1:3),gmSurfaceMesh.face)
    view(0,90)
    hold on
    set(gcf, 'Position', get(0, 'Screensize'));
    for i = 1:1:length(f_active_probe)
        if f_active_probe(i) == 1
            h(1)=plot3(probes(i,1),probes(i,2),probes(i,3),'.r','MarkerSize',30,'DisplayName','Active ROI')
        else
            h(2)=plot3(probes(i,1),probes(i,2),probes(i,3),'.k','MarkerSize',30,'DisplayName','Deactive ROI')
        end
    end
    legend(h(1:2),'Active ROI','Inactive ROI')  
    title('Superior')


    sgtitle('ROIs center','FontSize',15)


    % print('FIG_36_Roi_act_PD8','-djpeg','-r600')

    %% ROI RAPPRESENTATION OVER GM SURFACE
    % Use 'tildelayout' to have big figure to use for thesis/publication
    figure()
    gmSurfaceMesh.node(:,5) = [];
    % t = tiledlayout(1,2);
    t = tiledlayout(1,1);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    nexttile
    plotmeshISO(gmSurfaceMesh.node,gmSurfaceMesh.face);
    view(2)
    set(gcf, 'Position', get(0, 'Screensize'));
    title(['ROIs distribution (r = ',num2str(r_ROI),' mm)'])


    %ROI_vect
    figure()
    plotmeshISO([gmSurfaceMesh.node],gmSurfaceMesh.face);
    clim([0 9])
    colormap jet
    view([0 90])
    
    % title(['POSIZIONE ROI (r = ',num2str(r_ROI),' mm)'])
    % print('ROI_dist','-djpeg','-r600')


    % nexttile
    % gmSurfaceProva.node(find(gmSurfaceProva.node(:,4)==2),4) = 0;
    % gmSurfaceProva.node(find(gmSurfaceProva.node(:,4)==6),4) = 0;
    % plotmesh(gmSurfaceProva.node,gmSurfaceProva.face);
    % view(2)
    % set(gcf, 'Position', get(0, 'Screensize'));
    % title(['ROIs distribution (r = ',num2str(r_ROI),' mm)'])
    % set(gcf, 'Position', get(0, 'Screensize'));

    % print('FIG_37_ALL_vs_PD14','-djpeg','-r600')

    %% ROIs USEFUL INFORMATION
    disp('=================================================================')
    disp('ROI')
    disp(' ')
    disp(['Number of ROI = ',num2str(n_roi),'/',num2str(n_active_probe)])
    disp(['Number of nodes without ROI  = ',num2str(n_not_roi),'/',num2str(n_nodes)])
    disp(['ROI smallest number of nodes = ',num2str(min_ROI_active),' - Indexes = ',num2str(idx_min_ROI_active')])
    disp(['ROI largest number of nodes  = ',num2str(max_ROI_active),' - Indexes = ',num2str(idx_max_ROI_active')])
    disp(' ')
    disp('ACTIVE NODES FOR EACH ROIs')
    disp(ROI_count_table)
    disp('=================================================================')


    %% DISPLAY rsFC MATRIX
    
    % For better rappresentation of results use tildelayout (paper/thesis)
    
    % HbO
    figure()

    t = tiledlayout(1,2);   
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    nexttile
    imagesc(hbo_FC_sort)
    axis square
    colormap jet
    caxis([-1 1])
    colorbar
    title('HbO rsFC MATRIX')
    set(gca,'fontsize',12,'fontweight','bold')
    
    xticklabels = label_active_nodes_sort_complete;
    % xticklabels = label_active_nodes_sort;
    xticks = linspace(1, size(FC_hbo, 2), numel(xticklabels));
    xtickangle(45)
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    yticklabels = label_active_nodes_sort_complete;
    % yticklabels = label_active_nodes_sort;
    yticks = linspace(1, size(FC_hbo, 2), numel(yticklabels));
    % ytickangle(45)
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    
    nexttile

    imagesc(hbo_p_sort)
    colormap jet
    axis square
    caxis([0 p_bonf])
    colorbar
    % title('HbO P-VALUES')
    title('HbO P-VALUES')
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'fontsize',12,'fontweight','bold')
    
    xticklabels = label_active_nodes_sort;
    xtickangle(45)
    xticks = linspace(1, size(FC_hbo, 2), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    yticklabels = label_active_nodes_sort;
    yticks = linspace(1, size(FC_hbo, 2), numel(yticklabels));
    % ytickangle(45)
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    
    % HbR
    
    figure()

    t = tiledlayout(1,2);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    nexttile
    imagesc(hbr_FC_sort)
    axis square
    colormap jet
    caxis([-1 1])
    colorbar
    title('HbR rsFC MATRIX')
    set(gca,'fontsize',12,'fontweight','bold')
    
    xticklabels = label_active_nodes_sort_complete;
    % xticklabels = label_active_nodes_sort;
    xticks = linspace(1, size(FC_hbr, 2), numel(xticklabels));
    xtickangle(45)
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    yticklabels = label_active_nodes_sort_complete;
    % yticklabels = label_active_nodes_sort;
    yticks = linspace(1, size(FC_hbr, 2), numel(yticklabels));
    % ytickangle(45)
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    
    nexttile

    imagesc(hbr_p_sort)
    colormap jet
    axis square
    caxis([0 p_bonf])
    colorbar
    % title('HbR P-VALUES')
    title('HbR P-VALUES')
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'fontsize',12,'fontweight','bold')
    
    xticklabels = label_active_nodes_sort;
    xtickangle(45)
    xticks = linspace(1, size(FC_hbr, 2), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    yticklabels = label_active_nodes_sort;
    yticks = linspace(1, size(FC_hbr, 2), numel(yticklabels));
    % ytickangle(45)
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

end

end % FUNCTION END