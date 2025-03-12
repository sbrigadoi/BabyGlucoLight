%Guy Perkins 2025 (I hope my code is somewhat useable for you, please email
%me at guyantony.perkins@unipd.it if you need help!! Thank you.
%Follow each numbered section in order and read if you run it each time, or only
%once
%% SECTION 1 Check all events for given subject - make sure to change subject N and event type name do once for each subject/time window/event type
% This also shows the glucose plots for each event
clear all;close all; clc;
subjectN = 49;
eventType = "m_hypo"; %S_hypo
%eventType = "S_hypo"; %S_hypo
time_window = "A"; % A B C D E F (do A D F)
%time_window = "D"; % A B C D E F (do A D F)
%time_window = "F"; % A B C D E F (do A D F)
[PD_data] = view_all_events_subject(subjectN,eventType,time_window);
%% SECTION 2 Check the events to see how good/bad the signal is - do once for each subject/time window/event type
close all; clear all; clc;
subjectN = 49;
eventType = "m_hypo"; %S_hypo
%eventType = "S_hypo"; %S_hypo
time_window = "A"; % A B C D E F
%time_window = "D"; % A B C D E F
%time_window = "F"; % A B C D E F
[events_stats] = N_good_events_check_PD_10mins(subjectN,eventType,time_window,1)
save("event_stats_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","events_stats");
%% SECTION 3  Load data by subject N, event N, event type - do this everytime
%clear all; 
clc; close all; clearvars -except J mesh scalpSurfaceMesh gmSurfaceMesh headVolumeMesh landmarks vol2gm tenFive coverageThresh
w=0;
subjectN = 49;
%select one option
eventType = "m_hypo"; %for m_hypo only or m and S combined
%eventType = "S_hypo"; %S_hypo only

%select one option
eventType_all = "m-hypo"; %just m-hypo
%eventType_all = "S-hypo"; %just S-hypo
%eventType_all = "S-and-m-hypo"; %both m and s hypo events stacked
time_window = "A";

switch eventType
    case "m_hypo"
        mild_and_severe = 0; %1 YES , 0 NO

    switch eventType_all
        case "m-hypo"
             mild_and_severe = 0; %1 YES , 0 NO
        case "S-and-m-hypo"
             mild_and_severe = 1; %1 YES , 0 NO
    end
    case "S_hypo"
         mild_and_severe = 0; %1 YES , 0 NO
end

t_points = 601; %10 mins (s)
weekN_J = choose_weekN_mesh_jac_infant(subjectN);
%use this 26 07 23                  `               

load("Jacobian_infant"+num2str(weekN_J)+"week_850nm.mat")
load("mesh_infant"+num2str(weekN_J)+"week_850nm.mat")
load("AllMeshes_"+num2str(weekN_J)+"weeks.mshs","-mat")

% %%% TO DO - Set mesh to correct week for subject N

gmSurfaceMesh.vol2gm = vol2gm;
[coverageThresh] = get_Jac_threshold(headVolumeMesh,vol2gm);
N_nodes = size(headVolumeMesh.node,1);
% Based upon events stats - select events to pass through to gen. HbX
%returns event numbers that are good
%good meaning > 23 good chs and <50% MA train out of all T pts
%and meet Jacobian threshold ROI
[event_good] = choose_good_events_for_tomo(eventType,subjectN,time_window);
N_eventN = size(event_good,2);

PD_meta_data = get_PD_meta_data();

% Run this only once - Event Jac threshold
%[event_jac_threshold_idx] = event_jac_threshold_check(headVolumeMesh,gmSurfaceMesh,J_gmsurface,good_ch_idx_events,masks,coverageThresh);
%event_jac_threshold_idx = [event_jac_threshold_idx PD_m_hypo_events'];
%save("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
%% SECTION 4 - run this to gather HbT for all events - Only have to do once for each subject/time window/event type
%get_HbT_all_events_single_subject
clc;
switch eventType
    case "m_hypo"
        %HbT for all GM nodes in full mesh
        %HbT_all_m_hypo = zeros(size(GM_nodes,1),t_points,N_eventN);
        %HbT for all GM nodes in smooth surface mesh (save space)
        HbT_all_m_hypo_gmsurface = zeros(size(gmSurfaceMesh.node,1),t_points,N_eventN); 

        peak_glucose_all_m_hypo = zeros(N_eventN,1);
        length_glucose_all_m_hypo = zeros(N_eventN,1);
        good_ch_idx_events = zeros(128,N_eventN);
        
        for i=1:N_eventN
            %remove giacomo gdrive from path is HmrSG not working
            [PD_data] = Iteration10minDOT(subjectN,eventType,time_window,i,J,mesh,gmSurfaceMesh,t_points,vol2gm,coverageThresh,event_good);
            %HbO_all_m_hypo(:,:,i) = PD_data.HbO;
            %HbT_all_m_hypo(:,:,i) = PD_data.HbT(GM_nodes,:);
            HbT_all_m_hypo_gmsurface(:,:,i) = PD_data.HbT(:,:);
            peak_glucose_all_m_hypo(i) = PD_data.peak_glucose;
            length_glucose_all_m_hypo(i) = PD_data.length_glucose;
            good_ch_idx_events(1:size(PD_data.goodch_idx,1),i) = PD_data.goodch_idx;
            clc;
            i %keep track of what number we are on
        end
    
    PD_m_hypo_events = event_good;
    save("HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    save("peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    save("length_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","length_glucose_all_m_hypo");
    save("good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
    save("PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    %save("PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_good");


    case "S_hypo"
        %HbO_all_S_hypo = zeros(size(GM_nodes,1),t_points,N_eventN);
        %HbT_all_S_hypo = zeros(size(GM_nodes,1),t_points,N_eventN);
        HbT_all_S_hypo_gmsurface = zeros(size(gmSurfaceMesh.node,1),t_points,N_eventN); 
        peak_glucose_all_S_hypo = zeros(N_eventN,1);
        length_glucose_all_S_hypo = zeros(N_eventN,1);
        good_ch_idx_events = zeros(128,N_eventN);

        
        for i=1:N_eventN
            [PD_data] = Iteration10minDOT(subjectN,eventType,time_window,i,J,mesh,gmSurfaceMesh,t_points,vol2gm,coverageThresh,event_good);
            %HbO_all_S_hypo(:,:,i) = PD_data.HbO;
            HbT_all_S_hypo_gmsurface(:,:,i) = PD_data.HbT(:,:);
            peak_glucose_all_S_hypo(i) = PD_data.peak_glucose;
            length_glucose_all_S_hypo(i) = PD_data.length_glucose;
            good_ch_idx_events(1:size(PD_data.goodch_idx,1),i) = PD_data.goodch_idx;
            clc;
            i
        end
    
    PD_S_hypo_events = event_good;
    save("HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
    save("peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
    save("length_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","length_glucose_all_S_hypo");    
    save("good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
    save("PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
end
%% SECTION 5 - Load HbT data - do this EVERY TIME
HDD_text = "H"; %H PC , F laptop
%you can change your path to find the correct variables that you've
%generated, or can use mine.
switch eventType
    case "m_hypo"
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\length_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","length_glucose_all_m_hypo");
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");

    case "S_hypo"
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\length_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","length_glucose_all_S_hypo");
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
end
% Run each time
% Run this create mask for 95 and 99% and J thresh of Jacobian - do this evertime
PD.subjectN = subjectN;
PD_data.time_window= time_window;
PD.eventType = eventType;
PD.eventN = 0;
[masks] = get_J_mask(mesh,J,PD,PD_data,vol2gm,gmSurfaceMesh,coverageThresh);

% RUN THIS Jacobian mapping
%map J.complete (Intensity J) to gmsurface
J_gmsurface = vol2gm*J.complete';
%sum this across all 64 channels
J_all_gmsurface = sum(J_gmsurface,2);
%map meshsupport to gmsurface
mesh_support_gmsurface = vol2gm*mesh.support(:,1);
% spatially normalise J all gmsurface 
J_all_gmsurface_norm=J_all_gmsurface./mesh_support_gmsurface;
%plot J all gmsurface
figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface],gmSurfaceMesh.face);view(0,90);
cb = colorbar('horiz');
title("Jacobian - GM Smooth Week30")
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm");

%% SECTION 6 Run this only once - to get Event Jac threshold - single subject - single time window and single event type
[event_jac_threshold_idx] = event_jac_threshold_check(headVolumeMesh,gmSurfaceMesh,J_gmsurface,good_ch_idx_events,masks,coverageThresh);
switch eventType
    case "m_hypo"
        event_jac_threshold_idx = [event_jac_threshold_idx PD_m_hypo_events'];
    case "S_hypo"
        event_jac_threshold_idx = [event_jac_threshold_idx PD_S_hypo_events'];
end
save("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");

%% SECTION 7 For stacking all subjects only - S or m hypo ONLY (NO HYPER) 24 10 24 - RUN EACH TIME

switch eventType 
    case "m_hypo"
        if mild_and_severe == 1 %m and s hypo
            all_subject_N = [8 10 15 25 33 34 39 41 44 45 48 49 55 56 58 59 60 15 25 33 49];
                                %mhypo and then the shypo subjects
        else %just mhypo
            all_subject_N = [8 10 15 25 33 34 39 41 44 45 48 49 55 56 58 59 60];
        end
    case "S_hypo"
        all_subject_N = [15 25 33 49];

end

if eventType == "m_hypo"
    subjectN=8;
    HbT_all_m_hypo_gmsurface_all_8= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_8 = getfield(HbT_all_m_hypo_gmsurface_all_8,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_8 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_8 = getfield(event_jac_threshold_idx_8,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_8= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_8 = getfield(peak_glucose_all_m_hypo_8,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_8=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_8 = getfield(PD_m_hypo_events_8,'PD_m_hypo_events');

    subjectN=10;
    HbT_all_m_hypo_gmsurface_all_10= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_10 = getfield(HbT_all_m_hypo_gmsurface_all_10,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_10 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_10 = getfield(event_jac_threshold_idx_10,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_10= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_10 = getfield(peak_glucose_all_m_hypo_10,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_10=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_10 = getfield(PD_m_hypo_events_10,'PD_m_hypo_events');

    % subjectN=13;
    %         HbT_all_m_hypo_gmsurface_all_13= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    %         HbT_all_m_hypo_gmsurface_all_13 = getfield(HbT_all_m_hypo_gmsurface_all_13,'HbT_all_m_hypo_gmsurface');
    %         event_jac_threshold_idx_13 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    %         event_jac_threshold_idx_13 = getfield(event_jac_threshold_idx_13,'event_jac_threshold_idx');
    %         peak_glucose_all_m_hypo_13 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    %         peak_glucose_all_m_hypo_13 = getfield(peak_glucose_all_m_hypo_13,'peak_glucose_all_m_hypo');
    %         PD_m_hypo_events_13 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    %         PD_m_hypo_events_13 = getfield(PD_m_hypo_events_13,'PD_m_hypo_events');

    subjectN=15;
    HbT_all_m_hypo_gmsurface_all_15= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_15 = getfield(HbT_all_m_hypo_gmsurface_all_15,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_15 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_15 = getfield(event_jac_threshold_idx_15,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_15 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_15 = getfield(peak_glucose_all_m_hypo_15,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_15 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_15 = getfield(PD_m_hypo_events_15,'PD_m_hypo_events');

    subjectN=25;
    HbT_all_m_hypo_gmsurface_all_25= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_25 = getfield(HbT_all_m_hypo_gmsurface_all_25,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_25 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_25 = getfield(event_jac_threshold_idx_25,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_25 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_25 = getfield(peak_glucose_all_m_hypo_25,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_25 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_25 = getfield(PD_m_hypo_events_25,'PD_m_hypo_events');

    % %Yes A, Yes D, NOT F
    % subjectN=28;
    %         HbT_all_m_hypo_gmsurface_all_28= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    %         HbT_all_m_hypo_gmsurface_all_28 = getfield(HbT_all_m_hypo_gmsurface_all_28,'HbT_all_m_hypo_gmsurface');
    %         %event_jac_threshold_idx_15 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    %         %event_jac_threshold_idx_15 = getfield(event_jac_threshold_idx_15,'event_jac_threshold_idx');
    %         peak_glucose_all_m_hypo_28= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    %         peak_glucose_all_m_hypo_28 = getfield(peak_glucose_all_m_hypo_28,'peak_glucose_all_m_hypo');
    %         PD_m_hypo_events_28 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    %         PD_m_hypo_events_28 = getfield(PD_m_hypo_events_28,'PD_m_hypo_events');

    subjectN=33;
    HbT_all_m_hypo_gmsurface_all_33 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_33 = getfield(HbT_all_m_hypo_gmsurface_all_33,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_33 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_33 = getfield(event_jac_threshold_idx_33,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_33 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_33 = getfield(peak_glucose_all_m_hypo_33,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_33 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_33 = getfield(PD_m_hypo_events_33,'PD_m_hypo_events');

    %Not A, Yes D, Not F
    subjectN=34;
    %HbT_all_m_hypo_gmsurface_all_34= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    %HbT_all_m_hypo_gmsurface_all_34 = getfield(HbT_all_m_hypo_gmsurface_all_34,'HbT_all_m_hypo_gmsurface');
    %event_jac_threshold_idx_15 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    %event_jac_threshold_idx_15 = getfield(event_jac_threshold_idx_15,'event_jac_threshold_idx');
    %peak_glucose_all_m_hypo_34= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    %peak_glucose_all_m_hypo_34 = getfield(peak_glucose_all_m_hypo_34,'peak_glucose_all_m_hypo');
    %PD_m_hypo_events_34 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    %PD_m_hypo_events_34 = getfield(PD_m_hypo_events_34,'PD_m_hypo_events');

    subjectN=39;
    HbT_all_m_hypo_gmsurface_all_39 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_39 = getfield(HbT_all_m_hypo_gmsurface_all_39,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_39 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_39 = getfield(event_jac_threshold_idx_39,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_39 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_39 = getfield(peak_glucose_all_m_hypo_39,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_39 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_39 = getfield(PD_m_hypo_events_39,'PD_m_hypo_events');

    subjectN=41;
    HbT_all_m_hypo_gmsurface_all_41 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_41 = getfield(HbT_all_m_hypo_gmsurface_all_41,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_41 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_41 = getfield(event_jac_threshold_idx_41,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_41 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_41 = getfield(peak_glucose_all_m_hypo_41,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_41 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_41 = getfield(PD_m_hypo_events_41,'PD_m_hypo_events');

    subjectN=44;
    HbT_all_m_hypo_gmsurface_all_44 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_44 = getfield(HbT_all_m_hypo_gmsurface_all_44,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_44 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_44 = getfield(event_jac_threshold_idx_44,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_44 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_44 = getfield(peak_glucose_all_m_hypo_44,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_44 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_44 = getfield(PD_m_hypo_events_44,'PD_m_hypo_events');

    subjectN=45;
    HbT_all_m_hypo_gmsurface_all_45 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_45 = getfield(HbT_all_m_hypo_gmsurface_all_45,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_45 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_45 = getfield(event_jac_threshold_idx_45,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_45 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_45 = getfield(peak_glucose_all_m_hypo_45,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_45 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_45 = getfield(PD_m_hypo_events_45,'PD_m_hypo_events');


    subjectN=49;
    HbT_all_m_hypo_gmsurface_all_49 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_49 = getfield(HbT_all_m_hypo_gmsurface_all_49,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_49 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_49 = getfield(event_jac_threshold_idx_49,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_49 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_49 = getfield(peak_glucose_all_m_hypo_49,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_49 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_49 = getfield(PD_m_hypo_events_49,'PD_m_hypo_events');

    subjectN=55;
    HbT_all_m_hypo_gmsurface_all_55 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_55 = getfield(HbT_all_m_hypo_gmsurface_all_55,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_55 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_55 = getfield(event_jac_threshold_idx_55,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_55 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_55 = getfield(peak_glucose_all_m_hypo_55,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_55 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_55 = getfield(PD_m_hypo_events_55,'PD_m_hypo_events');

    subjectN=56;
    HbT_all_m_hypo_gmsurface_all_56 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_56 = getfield(HbT_all_m_hypo_gmsurface_all_56,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_56 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_56 = getfield(event_jac_threshold_idx_56,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_56 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_56 = getfield(peak_glucose_all_m_hypo_56,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_56 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_56 = getfield(PD_m_hypo_events_56,'PD_m_hypo_events');
    %Not A, Yes D , Yes F
    subjectN=58;
    HbT_all_m_hypo_gmsurface_all_58 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_58 = getfield(HbT_all_m_hypo_gmsurface_all_58,'HbT_all_m_hypo_gmsurface');
    %event_jac_threshold_idx_58 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    %event_jac_threshold_idx_58 = getfield(event_jac_threshold_idx_58,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_58 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_58 = getfield(peak_glucose_all_m_hypo_58,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_58 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_58 = getfield(PD_m_hypo_events_58,'PD_m_hypo_events');

    subjectN=59;
    HbT_all_m_hypo_gmsurface_all_59 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_59 = getfield(HbT_all_m_hypo_gmsurface_all_59,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_59 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_59 = getfield(event_jac_threshold_idx_59,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_59 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_59 = getfield(peak_glucose_all_m_hypo_59,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_59 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_59 = getfield(PD_m_hypo_events_59,'PD_m_hypo_events');

    subjectN=60;
    HbT_all_m_hypo_gmsurface_all_60 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
    HbT_all_m_hypo_gmsurface_all_60 = getfield(HbT_all_m_hypo_gmsurface_all_60,'HbT_all_m_hypo_gmsurface');
    event_jac_threshold_idx_60 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_60 = getfield(event_jac_threshold_idx_60,'event_jac_threshold_idx');
    peak_glucose_all_m_hypo_60 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
    peak_glucose_all_m_hypo_60 = getfield(peak_glucose_all_m_hypo_60,'peak_glucose_all_m_hypo');
    PD_m_hypo_events_60 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
    PD_m_hypo_events_60 = getfield(PD_m_hypo_events_60,'PD_m_hypo_events');


    %load  event jac  and HbT and peak from good events in subjects with
    %missing time windows
    switch time_window
        case "A"
            all_subject_N = [8 10 15 25 33 39 41 44 45 48 49 55 56 59 60];

            subjectN=48;
            HbT_all_m_hypo_gmsurface_all_48 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
            HbT_all_m_hypo_gmsurface_all_48 = getfield(HbT_all_m_hypo_gmsurface_all_48,'HbT_all_m_hypo_gmsurface');
            event_jac_threshold_idx_48 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
            event_jac_threshold_idx_48 = getfield(event_jac_threshold_idx_48,'event_jac_threshold_idx');
            peak_glucose_all_m_hypo_48 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
            peak_glucose_all_m_hypo_48 = getfield(peak_glucose_all_m_hypo_48,'peak_glucose_all_m_hypo');
            PD_m_hypo_events_48 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
            PD_m_hypo_events_48 = getfield(PD_m_hypo_events_48,'PD_m_hypo_events');
            good_m_hypo_48=find(event_jac_threshold_idx_48(:,1)==1);

            HbT_all_gmsurface_48 = HbT_all_m_hypo_gmsurface_all_48(:,:,good_m_hypo_48);
            peak_glucose_all_48 = peak_glucose_all_m_hypo_48(good_m_hypo_48);

            good_events_listidx_48 = good_m_hypo_48;
            good_PD_events_48 = PD_m_hypo_events_48(good_events_listidx_48);


            % event_jac_threshold_idx_28 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(28)+"_"+time_window+".mat","event_jac_threshold_idx");
            % event_jac_threshold_idx_28 = getfield(event_jac_threshold_idx_28,'event_jac_threshold_idx');
            %
            % good_m_hypo_28=find(event_jac_threshold_idx_28(:,1)==1);
            % HbT_all_gmsurface_28 = HbT_all_m_hypo_gmsurface_all_28(:,:,good_m_hypo_28);
            % peak_glucose_all_28 = peak_glucose_all_m_hypo_28(good_m_hypo_28);
            %
            % good_events_listidx_28 = good_m_hypo_28;
            % good_PD_events_28 = PD_m_hypo_events_28(good_events_listidx_28);

        case "D"
            all_subject_N = [8 10 15 25 33 34 39 41 44 45 49 55 56 58 59 60];

            event_jac_threshold_idx_58 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(58)+"_"+time_window+".mat","event_jac_threshold_idx");
            event_jac_threshold_idx_58 = getfield(event_jac_threshold_idx_58,'event_jac_threshold_idx');

            good_m_hypo_58=find(event_jac_threshold_idx_58(:,1)==1);

            HbT_all_gmsurface_58 = HbT_all_m_hypo_gmsurface_all_58(:,:,good_m_hypo_58);
            peak_glucose_all_58 = peak_glucose_all_m_hypo_58(good_m_hypo_58);

            good_events_listidx_58 = good_m_hypo_58;
            good_PD_events_58 = PD_m_hypo_events_58(good_events_listidx_58);

            event_jac_threshold_idx_34 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(34)+"_"+time_window+".mat","event_jac_threshold_idx");
            event_jac_threshold_idx_34 = getfield(event_jac_threshold_idx_34,'event_jac_threshold_idx');

            good_m_hypo_34=find(event_jac_threshold_idx_34(:,1)==1);
            HbT_all_m_hypo_gmsurface_all_34= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
            HbT_all_m_hypo_gmsurface_all_34 = getfield(HbT_all_m_hypo_gmsurface_all_34,'HbT_all_m_hypo_gmsurface');

            HbT_all_gmsurface_34 = HbT_all_m_hypo_gmsurface_all_34(:,:,good_m_hypo_34);

            peak_glucose_all_m_hypo_34= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
            peak_glucose_all_m_hypo_34 = getfield(peak_glucose_all_m_hypo_34,'peak_glucose_all_m_hypo');

            peak_glucose_all_34 = peak_glucose_all_m_hypo_34(good_m_hypo_34);

            PD_m_hypo_events_34 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
            PD_m_hypo_events_34 = getfield(PD_m_hypo_events_34,'PD_m_hypo_events');
            good_events_listidx_34 = good_m_hypo_34;
            good_PD_events_34 = PD_m_hypo_events_34(good_events_listidx_34);


        case "F"
            all_subject_N = [8 10 15 25 33 39 41 44 45 49 55 56 58 59 60];

            event_jac_threshold_idx_58 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(58)+"_"+time_window+".mat","event_jac_threshold_idx");
            event_jac_threshold_idx_58 = getfield(event_jac_threshold_idx_58,'event_jac_threshold_idx');

            good_m_hypo_58=find(event_jac_threshold_idx_58(:,1)==1);
            HbT_all_gmsurface_58 = HbT_all_m_hypo_gmsurface_all_58(:,:,good_m_hypo_58);
            peak_glucose_all_58 = peak_glucose_all_m_hypo_58(good_m_hypo_58);

            good_events_listidx_58 = good_m_hypo_58;
            good_PD_events_58 = PD_m_hypo_events_58(good_events_listidx_58);
    end

    %putting everything together %miss out34 and 58 since these r defined above
    good_m_hypo_8=find(event_jac_threshold_idx_8(:,1)==1);
    good_m_hypo_10=find(event_jac_threshold_idx_10(:,1)==1);
    %good_m_hypo_13=find(event_jac_threshold_idx_13(:,1)==1);
    good_m_hypo_15=find(event_jac_threshold_idx_15(:,1)==1);
    good_m_hypo_25=find(event_jac_threshold_idx_25(:,1)==1);
    good_m_hypo_33=find(event_jac_threshold_idx_33(:,1)==1);
    good_m_hypo_39=find(event_jac_threshold_idx_39(:,1)==1);
    good_m_hypo_41=find(event_jac_threshold_idx_41(:,1)==1);
    good_m_hypo_44=find(event_jac_threshold_idx_44(:,1)==1);
    good_m_hypo_45=find(event_jac_threshold_idx_45(:,1)==1);
    good_m_hypo_49=find(event_jac_threshold_idx_49(:,1)==1);
    good_m_hypo_55=find(event_jac_threshold_idx_55(:,1)==1);
    good_m_hypo_56=find(event_jac_threshold_idx_56(:,1)==1);
    %good_m_hypo_58=find(event_jac_threshold_idx_58(:,1)==1);
    good_m_hypo_59=find(event_jac_threshold_idx_59(:,1)==1);
    good_m_hypo_60=find(event_jac_threshold_idx_60(:,1)==1);

    HbT_all_gmsurface_8 = HbT_all_m_hypo_gmsurface_all_8(:,:,good_m_hypo_8);
    HbT_all_gmsurface_10 = HbT_all_m_hypo_gmsurface_all_10(:,:,good_m_hypo_10);
    %HbT_all_gmsurface_13 = HbT_all_m_hypo_gmsurface_all_13(:,:,good_m_hypo_13);
    HbT_all_gmsurface_15 = HbT_all_m_hypo_gmsurface_all_15(:,:,good_m_hypo_15);
    HbT_all_gmsurface_25 = HbT_all_m_hypo_gmsurface_all_25(:,:,good_m_hypo_25);
    HbT_all_gmsurface_33 = HbT_all_m_hypo_gmsurface_all_33(:,:,good_m_hypo_33);
    HbT_all_gmsurface_39 = HbT_all_m_hypo_gmsurface_all_39(:,:,good_m_hypo_39);
    HbT_all_gmsurface_41 = HbT_all_m_hypo_gmsurface_all_41(:,:,good_m_hypo_41);
    HbT_all_gmsurface_44 = HbT_all_m_hypo_gmsurface_all_44(:,:,good_m_hypo_44);
    HbT_all_gmsurface_45 = HbT_all_m_hypo_gmsurface_all_45(:,:,good_m_hypo_45);
    HbT_all_gmsurface_49 = HbT_all_m_hypo_gmsurface_all_49(:,:,good_m_hypo_49);
    HbT_all_gmsurface_55 = HbT_all_m_hypo_gmsurface_all_55(:,:,good_m_hypo_55);
    HbT_all_gmsurface_56 = HbT_all_m_hypo_gmsurface_all_56(:,:,good_m_hypo_56);
    %HbT_all_gmsurface_58 = HbT_all_m_hypo_gmsurface_all_58(:,:,good_m_hypo_58);
    HbT_all_gmsurface_59 = HbT_all_m_hypo_gmsurface_all_59(:,:,good_m_hypo_59);
    HbT_all_gmsurface_60 = HbT_all_m_hypo_gmsurface_all_60(:,:,good_m_hypo_60);

    peak_glucose_all_8 = peak_glucose_all_m_hypo_8(good_m_hypo_8);
    peak_glucose_all_10 = peak_glucose_all_m_hypo_10(good_m_hypo_10);
    %peak_glucose_all_13 = peak_glucose_all_m_hypo_13(good_m_hypo_13);
    peak_glucose_all_15 = peak_glucose_all_m_hypo_15(good_m_hypo_15);
    peak_glucose_all_25 = peak_glucose_all_m_hypo_25(good_m_hypo_25);
    peak_glucose_all_33 = peak_glucose_all_m_hypo_33(good_m_hypo_33);
    peak_glucose_all_39 = peak_glucose_all_m_hypo_39(good_m_hypo_39);
    peak_glucose_all_41 = peak_glucose_all_m_hypo_41(good_m_hypo_41);
    peak_glucose_all_44 = peak_glucose_all_m_hypo_44(good_m_hypo_44);
    peak_glucose_all_45 = peak_glucose_all_m_hypo_45(good_m_hypo_45);
    peak_glucose_all_49 = peak_glucose_all_m_hypo_49(good_m_hypo_49);
    peak_glucose_all_55 = peak_glucose_all_m_hypo_55(good_m_hypo_55);
    peak_glucose_all_56 = peak_glucose_all_m_hypo_56(good_m_hypo_56);
    %peak_glucose_all_58 = peak_glucose_all_m_hypo_58(good_m_hypo_58);
    peak_glucose_all_59 = peak_glucose_all_m_hypo_59(good_m_hypo_59);
    peak_glucose_all_60 = peak_glucose_all_m_hypo_60(good_m_hypo_60);



    good_events_listidx_8 = good_m_hypo_8;
    good_PD_events_8 = PD_m_hypo_events_8(good_events_listidx_8);
    good_events_listidx_10 = good_m_hypo_10;
    good_PD_events_10 = PD_m_hypo_events_10(good_events_listidx_10);
    %good_events_listidx_13 = good_m_hypo_13;
    %            good_PD_events_13 = PD_m_hypo_events_13(good_events_listidx_13);
    good_events_listidx_15 = good_m_hypo_15;
    good_PD_events_15 = PD_m_hypo_events_15(good_events_listidx_15);
    good_events_listidx_25 = good_m_hypo_25;
    good_PD_events_25 = PD_m_hypo_events_25(good_events_listidx_25);
    good_events_listidx_33 = good_m_hypo_33;
    good_PD_events_33 = PD_m_hypo_events_33(good_events_listidx_33);
    good_PD_events_39 = PD_m_hypo_events_39(good_m_hypo_39);
    good_PD_events_41 = PD_m_hypo_events_41(good_m_hypo_41);
    good_PD_events_44 = PD_m_hypo_events_44(good_m_hypo_44);
    good_PD_events_45 = PD_m_hypo_events_45(good_m_hypo_45);
    good_events_listidx_49 = good_m_hypo_49;
    good_PD_events_49 = PD_m_hypo_events_49(good_events_listidx_49);
    good_PD_events_55 = PD_m_hypo_events_55(good_m_hypo_55);
    good_PD_events_56 = PD_m_hypo_events_56(good_m_hypo_56);
    %good_PD_events_58 = PD_m_hypo_events_58(good_m_hypo_58);
    good_PD_events_59 = PD_m_hypo_events_59(good_m_hypo_59);
    good_PD_events_60 = PD_m_hypo_events_60(good_m_hypo_60);

    % stack variables
    %rm 13 from all windows
    switch time_window
        case "A" %NO 34 58 %yes 48
            HbT_all_gmsurface = cat(3,HbT_all_gmsurface_8,HbT_all_gmsurface_10,HbT_all_gmsurface_15,HbT_all_gmsurface_25,HbT_all_gmsurface_33,HbT_all_gmsurface_39,HbT_all_gmsurface_41,HbT_all_gmsurface_44,HbT_all_gmsurface_45,HbT_all_gmsurface_48,HbT_all_gmsurface_49,HbT_all_gmsurface_55,HbT_all_gmsurface_56,HbT_all_gmsurface_59,HbT_all_gmsurface_60);
            peak_glucose_all = [peak_glucose_all_8;peak_glucose_all_10;peak_glucose_all_15;peak_glucose_all_25;peak_glucose_all_33;peak_glucose_all_39;peak_glucose_all_41;peak_glucose_all_44;peak_glucose_all_45;peak_glucose_all_48;peak_glucose_all_49;peak_glucose_all_55;peak_glucose_all_56;peak_glucose_all_59;peak_glucose_all_60];
            good_PD_events = [good_PD_events_8 good_PD_events_10 good_PD_events_15 good_PD_events_25 good_PD_events_33 good_PD_events_39 good_PD_events_41 good_PD_events_44 good_PD_events_45 good_PD_events_48 good_PD_events_49 good_PD_events_55 good_PD_events_56 good_PD_events_59 good_PD_events_60];
            good_PD_events_subN = [ones(1,size(good_PD_events_8,2))*8 ones(1,size(good_PD_events_10,2))*10 ones(1,size(good_PD_events_15,2))*15 ones(1,size(good_PD_events_25,2))*25 ones(1,size(good_PD_events_33,2))*33 ones(1,size(good_PD_events_39,2))*39 ones(1,size(good_PD_events_41,2))*41 ones(1,size(good_PD_events_44,2))*44 ones(1,size(good_PD_events_45,2))*45 ones(1,size(good_PD_events_48,2))*48 ones(1,size(good_PD_events_49,2))*49 ones(1,size(good_PD_events_55,2))*55 ones(1,size(good_PD_events_56,2))*56 ones(1,size(good_PD_events_59,2))*59 ones(1,size(good_PD_events_60,2))*60];
            %good_events_listidx = [good_m_hypo_8;good_m_hypo_8];


        case "D" %yes 34 58 %no 48
            HbT_all_gmsurface = cat(3,HbT_all_gmsurface_8,HbT_all_gmsurface_10,HbT_all_gmsurface_15,HbT_all_gmsurface_25,HbT_all_gmsurface_33,HbT_all_gmsurface_34,HbT_all_gmsurface_39,HbT_all_gmsurface_41,HbT_all_gmsurface_44,HbT_all_gmsurface_45,HbT_all_gmsurface_49,HbT_all_gmsurface_55,HbT_all_gmsurface_56,HbT_all_gmsurface_58,HbT_all_gmsurface_59,HbT_all_gmsurface_60);
            peak_glucose_all = [peak_glucose_all_8;peak_glucose_all_10;peak_glucose_all_15;peak_glucose_all_25;peak_glucose_all_33;peak_glucose_all_34;peak_glucose_all_39;peak_glucose_all_41;peak_glucose_all_44;peak_glucose_all_45;peak_glucose_all_49;peak_glucose_all_55;peak_glucose_all_56;peak_glucose_all_58;peak_glucose_all_59;peak_glucose_all_60];
            good_PD_events = [good_PD_events_8 good_PD_events_10 good_PD_events_15 good_PD_events_25 good_PD_events_33 good_PD_events_34 good_PD_events_39 good_PD_events_41 good_PD_events_44 good_PD_events_45 good_PD_events_49 good_PD_events_55 good_PD_events_56 good_PD_events_58 good_PD_events_59 good_PD_events_60];
            good_PD_events_subN = [ones(1,size(good_PD_events_8,2))*8 ones(1,size(good_PD_events_10,2))*10 ones(1,size(good_PD_events_15,2))*15 ones(1,size(good_PD_events_25,2))*25 ones(1,size(good_PD_events_33,2))*33 ones(1,size(good_PD_events_34,2))*34 ones(1,size(good_PD_events_39,2))*39 ones(1,size(good_PD_events_41,2))*41 ones(1,size(good_PD_events_44,2))*44 ones(1,size(good_PD_events_45,2))*45 ones(1,size(good_PD_events_49,2))*49 ones(1,size(good_PD_events_55,2))*55 ones(1,size(good_PD_events_56,2))*56 ones(1,size(good_PD_events_58,2))*58 ones(1,size(good_PD_events_59,2))*59 ones(1,size(good_PD_events_60,2))*60];


        case "F" %no 34 48
            HbT_all_gmsurface = cat(3,HbT_all_gmsurface_8,HbT_all_gmsurface_10,HbT_all_gmsurface_15,HbT_all_gmsurface_25,HbT_all_gmsurface_33,HbT_all_gmsurface_39,HbT_all_gmsurface_41,HbT_all_gmsurface_44,HbT_all_gmsurface_45,HbT_all_gmsurface_49,HbT_all_gmsurface_55,HbT_all_gmsurface_56,HbT_all_gmsurface_58,HbT_all_gmsurface_59,HbT_all_gmsurface_60);
            peak_glucose_all = [peak_glucose_all_8;peak_glucose_all_10;peak_glucose_all_15;peak_glucose_all_25;peak_glucose_all_33;peak_glucose_all_39;peak_glucose_all_41;peak_glucose_all_44;peak_glucose_all_45;peak_glucose_all_49;peak_glucose_all_55;peak_glucose_all_56;peak_glucose_all_58;peak_glucose_all_59;peak_glucose_all_60];
            good_PD_events = [good_PD_events_8 good_PD_events_10 good_PD_events_15 good_PD_events_25 good_PD_events_33 good_PD_events_39 good_PD_events_41 good_PD_events_44 good_PD_events_45 good_PD_events_49 good_PD_events_55 good_PD_events_56 good_PD_events_58 good_PD_events_59 good_PD_events_60];
            good_PD_events_subN = [ones(1,size(good_PD_events_8,2))*8 ones(1,size(good_PD_events_10,2))*10 ones(1,size(good_PD_events_15,2))*15 ones(1,size(good_PD_events_25,2))*25 ones(1,size(good_PD_events_33,2))*33 ones(1,size(good_PD_events_39,2))*39 ones(1,size(good_PD_events_41,2))*41 ones(1,size(good_PD_events_44,2))*44 ones(1,size(good_PD_events_45,2))*45 ones(1,size(good_PD_events_49,2))*49 ones(1,size(good_PD_events_55,2))*55 ones(1,size(good_PD_events_56,2))*56 ones(1,size(good_PD_events_58,2))*58 ones(1,size(good_PD_events_59,2))*59 ones(1,size(good_PD_events_60,2))*60];

    end
    %update eventtype Sub N string %all mhypo
    eventType_subN = strings(size(good_PD_events,2),1);
    eventType_subN(:,1) = "m_hypo";

    eventType_subN_idx = ones(size(good_PD_events,2),1); %1=mhypo,2=Shypo,3=mhyper,4=Shyper


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mild_and_severe == 1 %=1 means we want both m and S hypo
        eventType = "S_hypo"; %for correct loading
        subjectN=15;
        HbT_all_S_hypo_gmsurface_all_15= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
        HbT_all_S_hypo_gmsurface_all_15 = getfield(HbT_all_S_hypo_gmsurface_all_15,'HbT_all_S_hypo_gmsurface');
        event_jac_threshold_idx_15 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_15 = getfield(event_jac_threshold_idx_15,'event_jac_threshold_idx');
        peak_glucose_all_S_hypo_15 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
        peak_glucose_all_S_hypo_15 = getfield(peak_glucose_all_S_hypo_15,'peak_glucose_all_S_hypo');
        PD_S_hypo_events_15 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
        PD_S_hypo_events_15 = getfield(PD_S_hypo_events_15,'PD_S_hypo_events');

        subjectN=25;
        HbT_all_S_hypo_gmsurface_all_25= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
        HbT_all_S_hypo_gmsurface_all_25 = getfield(HbT_all_S_hypo_gmsurface_all_25,'HbT_all_S_hypo_gmsurface');
        event_jac_threshold_idx_25 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_25 = getfield(event_jac_threshold_idx_25,'event_jac_threshold_idx');
        peak_glucose_all_S_hypo_25 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
        peak_glucose_all_S_hypo_25 = getfield(peak_glucose_all_S_hypo_25,'peak_glucose_all_S_hypo');
        PD_S_hypo_events_25 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
        PD_S_hypo_events_25 = getfield(PD_S_hypo_events_25,'PD_S_hypo_events');

        subjectN=33;
        HbT_all_S_hypo_gmsurface_all_33= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
        HbT_all_S_hypo_gmsurface_all_33 = getfield(HbT_all_S_hypo_gmsurface_all_33,'HbT_all_S_hypo_gmsurface');
        event_jac_threshold_idx_33 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_33 = getfield(event_jac_threshold_idx_33,'event_jac_threshold_idx');
        peak_glucose_all_S_hypo_33 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
        peak_glucose_all_S_hypo_33 = getfield(peak_glucose_all_S_hypo_33,'peak_glucose_all_S_hypo');
        PD_S_hypo_events_33 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
        PD_S_hypo_events_33 = getfield(PD_S_hypo_events_33,'PD_S_hypo_events');

        subjectN=49;
        HbT_all_S_hypo_gmsurface_all_49= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
        HbT_all_S_hypo_gmsurface_all_49 = getfield(HbT_all_S_hypo_gmsurface_all_49,'HbT_all_S_hypo_gmsurface');
        event_jac_threshold_idx_49 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_49 = getfield(event_jac_threshold_idx_49,'event_jac_threshold_idx');
        peak_glucose_all_S_hypo_49 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
        peak_glucose_all_S_hypo_49 = getfield(peak_glucose_all_S_hypo_49,'peak_glucose_all_S_hypo');
        PD_S_hypo_events_49 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
        PD_S_hypo_events_49 = getfield(PD_S_hypo_events_49,'PD_S_hypo_events');


        %putting everything together
        good_S_hypo_15=find(event_jac_threshold_idx_15(:,1)==1);
        good_S_hypo_25=find(event_jac_threshold_idx_25(:,1)==1);
        good_S_hypo_33=find(event_jac_threshold_idx_33(:,1)==1);
        good_S_hypo_49=find(event_jac_threshold_idx_49(:,1)==1);
   
        HbT_all_gmsurface_15 = HbT_all_S_hypo_gmsurface_all_15(:,:,good_S_hypo_15);
        HbT_all_gmsurface_25 = HbT_all_S_hypo_gmsurface_all_25(:,:,good_S_hypo_25);
        HbT_all_gmsurface_33 = HbT_all_S_hypo_gmsurface_all_33(:,:,good_S_hypo_33);
        HbT_all_gmsurface_49 = HbT_all_S_hypo_gmsurface_all_49(:,:,good_S_hypo_49);
    
        peak_glucose_all_15 = peak_glucose_all_S_hypo_15(good_S_hypo_15);
        peak_glucose_all_25 = peak_glucose_all_S_hypo_25(good_S_hypo_25);
        peak_glucose_all_33 = peak_glucose_all_S_hypo_33(good_S_hypo_33);
        peak_glucose_all_49 = peak_glucose_all_S_hypo_49(good_S_hypo_49);
   
        good_events_listidx_15 = good_S_hypo_15;
        good_PD_events_15 = PD_S_hypo_events_15(good_events_listidx_15);
        good_events_listidx_25 = good_S_hypo_25;
        good_PD_events_25 = PD_S_hypo_events_25(good_events_listidx_25);
        good_events_listidx_33 = good_S_hypo_33;
        good_PD_events_33 = PD_S_hypo_events_33(good_events_listidx_33);
        good_events_listidx_49 = good_S_hypo_49;
        good_PD_events_49 = PD_S_hypo_events_49(good_events_listidx_49);
    

        % stack variables on top of exisiting m hypo, mhypo at the start
        % and then we add on the Shypo info to the end
   
            HbT_all_gmsurface = cat(3,HbT_all_gmsurface,HbT_all_gmsurface_15,HbT_all_gmsurface_25,HbT_all_gmsurface_33,HbT_all_gmsurface_49);
            peak_glucose_all = [peak_glucose_all; peak_glucose_all_15;peak_glucose_all_25;peak_glucose_all_33;peak_glucose_all_49];
            good_PD_events_S_hypo_only = [good_PD_events_15 good_PD_events_25 good_PD_events_33 good_PD_events_49];

            good_PD_events = [good_PD_events good_PD_events_15 good_PD_events_25 good_PD_events_33 good_PD_events_49];
            good_PD_events_subN = [good_PD_events_subN ones(1,size(good_PD_events_15,2))*15 ones(1,size(good_PD_events_25,2))*25 ones(1,size(good_PD_events_33,2))*33 ones(1,size(good_PD_events_49,2))*49];
            
            %update eventtype Sub N string %all mhypo and add S hypo
            eventType_subN_S_hypo = strings(size(good_PD_events_S_hypo_only,2),1);
            eventType_subN_S_hypo(:,1) = "S_hypo";
            eventType_subN = [eventType_subN ; eventType_subN_S_hypo ];

            eventType = "m_hypo"; %correct change after loading Shypo data

            eventType_subN_idx_S_hypo_only = ones(size(good_PD_events_S_hypo_only,2),1)*2;
            eventType_subN_idx = [eventType_subN_idx;eventType_subN_idx_S_hypo_only]; %1=mhypo,2=Shypo,3=mhyper,4=Shyper
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif eventType == "S_hypo" %if eventtype = mhypo %JUST S Hypo


    subjectN=15;
    HbT_all_S_hypo_gmsurface_all_15= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
    HbT_all_S_hypo_gmsurface_all_15 = getfield(HbT_all_S_hypo_gmsurface_all_15,'HbT_all_S_hypo_gmsurface');
    event_jac_threshold_idx_15 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_15 = getfield(event_jac_threshold_idx_15,'event_jac_threshold_idx');
    peak_glucose_all_S_hypo_15 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
    peak_glucose_all_S_hypo_15 = getfield(peak_glucose_all_S_hypo_15,'peak_glucose_all_S_hypo');
    PD_S_hypo_events_15 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
    PD_S_hypo_events_15 = getfield(PD_S_hypo_events_15,'PD_S_hypo_events');

    subjectN=25;
    HbT_all_S_hypo_gmsurface_all_25= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
    HbT_all_S_hypo_gmsurface_all_25 = getfield(HbT_all_S_hypo_gmsurface_all_25,'HbT_all_S_hypo_gmsurface');
    event_jac_threshold_idx_25 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_25 = getfield(event_jac_threshold_idx_25,'event_jac_threshold_idx');
    peak_glucose_all_S_hypo_25 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
    peak_glucose_all_S_hypo_25 = getfield(peak_glucose_all_S_hypo_25,'peak_glucose_all_S_hypo');
    PD_S_hypo_events_25 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
    PD_S_hypo_events_25 = getfield(PD_S_hypo_events_25,'PD_S_hypo_events');

    subjectN=33;
    HbT_all_S_hypo_gmsurface_all_33= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
    HbT_all_S_hypo_gmsurface_all_33 = getfield(HbT_all_S_hypo_gmsurface_all_33,'HbT_all_S_hypo_gmsurface');
    event_jac_threshold_idx_33 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_33 = getfield(event_jac_threshold_idx_33,'event_jac_threshold_idx');
    peak_glucose_all_S_hypo_33 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
    peak_glucose_all_S_hypo_33 = getfield(peak_glucose_all_S_hypo_33,'peak_glucose_all_S_hypo');
    PD_S_hypo_events_33 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
    PD_S_hypo_events_33 = getfield(PD_S_hypo_events_33,'PD_S_hypo_events');

    subjectN=49;
    HbT_all_S_hypo_gmsurface_all_49= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
    HbT_all_S_hypo_gmsurface_all_49 = getfield(HbT_all_S_hypo_gmsurface_all_49,'HbT_all_S_hypo_gmsurface');
    event_jac_threshold_idx_49 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_49 = getfield(event_jac_threshold_idx_49,'event_jac_threshold_idx');
    peak_glucose_all_S_hypo_49 = load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
    peak_glucose_all_S_hypo_49 = getfield(peak_glucose_all_S_hypo_49,'peak_glucose_all_S_hypo');
    PD_S_hypo_events_49 =load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
    PD_S_hypo_events_49 = getfield(PD_S_hypo_events_49,'PD_S_hypo_events');


    %putting everything together
    good_S_hypo_15=find(event_jac_threshold_idx_15(:,1)==1);
    good_S_hypo_25=find(event_jac_threshold_idx_25(:,1)==1);
    good_S_hypo_33=find(event_jac_threshold_idx_33(:,1)==1);
    good_S_hypo_49=find(event_jac_threshold_idx_49(:,1)==1);
   
    HbT_all_gmsurface_15 = HbT_all_S_hypo_gmsurface_all_15(:,:,good_S_hypo_15);
    HbT_all_gmsurface_25 = HbT_all_S_hypo_gmsurface_all_25(:,:,good_S_hypo_25);
    HbT_all_gmsurface_33 = HbT_all_S_hypo_gmsurface_all_33(:,:,good_S_hypo_33);
    HbT_all_gmsurface_49 = HbT_all_S_hypo_gmsurface_all_49(:,:,good_S_hypo_49);
    
    peak_glucose_all_15 = peak_glucose_all_S_hypo_15(good_S_hypo_15);
    peak_glucose_all_25 = peak_glucose_all_S_hypo_25(good_S_hypo_25);
    peak_glucose_all_33 = peak_glucose_all_S_hypo_33(good_S_hypo_33);
    peak_glucose_all_49 = peak_glucose_all_S_hypo_49(good_S_hypo_49);
   
    good_events_listidx_15 = good_S_hypo_15;
    good_PD_events_15 = PD_S_hypo_events_15(good_events_listidx_15);
    good_events_listidx_25 = good_S_hypo_25;
    good_PD_events_25 = PD_S_hypo_events_25(good_events_listidx_25);
    good_events_listidx_33 = good_S_hypo_33;
    good_PD_events_33 = PD_S_hypo_events_33(good_events_listidx_33);
    good_events_listidx_49 = good_S_hypo_49;
    good_PD_events_49 = PD_S_hypo_events_49(good_events_listidx_49);
    

    % stack variables
   
            HbT_all_gmsurface = cat(3,HbT_all_gmsurface_15,HbT_all_gmsurface_25,HbT_all_gmsurface_33,HbT_all_gmsurface_49);
            peak_glucose_all = [peak_glucose_all_15;peak_glucose_all_25;peak_glucose_all_33;peak_glucose_all_49];
            good_PD_events = [good_PD_events_15 good_PD_events_25 good_PD_events_33 good_PD_events_49];
            good_PD_events_subN = [ones(1,size(good_PD_events_15,2))*15 ones(1,size(good_PD_events_25,2))*25 ones(1,size(good_PD_events_33,2))*33 ones(1,size(good_PD_events_49,2))*49];
            
        %update eventtype Sub N string %all Shypo
        eventType_subN = strings(size(good_PD_events,2),1);
        eventType_subN(:,1) = "S_hypo";
        eventType_subN_idx = ones(size(good_PD_events,2),1)*2; %1=mhypo,2=Shypo,3=mhyper,4=Shyper

end %if event type = m_hypo

% HbT_all_gmsurface = cat(3,HbT_all_gmsurface(:,:,1:6),HbT_all_gmsurface(:,:,8:9));
% peak_glucose_all = [peak_glucose_all(1:6) ; peak_glucose_all(8:9) ];
% good_PD_events = [good_PD_events(1:6) good_PD_events(8:9) ];
% good_PD_events_subN = [good_PD_events_subN(1:6) good_PD_events_subN(8:9) ];

%% SECTION 8 for single or stacked - run this each time
% run this - r and t analysis
savefig = 0; %1=save figures, 0 =don't save figures
stacked =1; %1 YES , 0 NO (stacked means are there multiple subjects stacked!)
mild_and_severe = 0; %1 YES , 0 NO

%get_R_T_values(subjectN,eventType,savefig,stacked,mild_and_severe,time_window,PD,vol2gm,masks,gmSurfaceMesh,HbT_all_gmsurface,peak_glucose_all,good_PD_events,eventType_all,PD_data,good_PD_events_subN);
%% SECTION 9 Get BGC and HbT metrics - run this each time
% eventType_subN , eventType_subN_idx 
[BGC_metrics_PD_events] = get_timewindow_10mins_BGCmetrics(good_PD_events_subN,good_PD_events,time_window,eventType,eventType_subN , eventType_subN_idx );
[HbT_metrics_PD_events] = get_timewindow_10mins_HbTmetrics(good_PD_events_subN,good_PD_events,time_window,eventType,HbT_all_gmsurface,gmSurfaceMesh,masks,eventType_subN , eventType_subN_idx );

%% SECTION 10 - Run this to produce R val graphs
R_val_2_metrics(subjectN,eventType,savefig,stacked,mild_and_severe,time_window,PD,vol2gm,masks,gmSurfaceMesh,HbT_all_gmsurface,peak_glucose_all,good_PD_events,eventType_all,PD_data,good_PD_events_subN,BGC_metrics_PD_events ,HbT_metrics_PD_events,eventType_subN , eventType_subN_idx )

%get some meta data about the L and R hemisphere
L_R_hemisphere_meta_data(HbT_metrics_PD_events,PD_meta_data,good_PD_events_subN,good_PD_events)

%% SECTION 11 - Run this to get T value data
t_val_metric(subjectN,eventType,savefig,stacked,mild_and_severe,time_window,PD,vol2gm,masks,gmSurfaceMesh,HbT_all_gmsurface,peak_glucose_all,good_PD_events,eventType_all,PD_data,good_PD_events_subN,BGC_metrics_PD_events ,HbT_metrics_PD_events,eventType_subN , eventType_subN_idx )



%%%%%%%% END
%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SKIP THIS For stacking all subjects only SEPT 2024

%stack_subjects_m_hypo
all_subject_N = [7 8 13 15 20 25 33 49];
%HbT_all_m_hypo_gmsurface_all_sub = zeros(size(vol2gm,1),t_points,100);

        subjectN=7;
        HbT_all_m_hypo_gmsurface_all_7= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        HbT_all_m_hypo_gmsurface_all_7 = getfield(HbT_all_m_hypo_gmsurface_all_7,'HbT_all_m_hypo_gmsurface');
        event_jac_threshold_idx_7 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_7 = getfield(event_jac_threshold_idx_7,'event_jac_threshold_idx');
        peak_glucose_all_m_hypo_7= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        peak_glucose_all_m_hypo_7 = getfield(peak_glucose_all_m_hypo_7,'peak_glucose_all_m_hypo');
        PD_m_hypo_events_7=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        PD_m_hypo_events_7 = getfield(PD_m_hypo_events_7,'PD_m_hypo_events');


        subjectN=8;
        HbT_all_m_hypo_gmsurface_all_8= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        HbT_all_m_hypo_gmsurface_all_8 = getfield(HbT_all_m_hypo_gmsurface_all_8,'HbT_all_m_hypo_gmsurface');        
        event_jac_threshold_idx_8 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_8 = getfield(event_jac_threshold_idx_8,'event_jac_threshold_idx');       
        peak_glucose_all_m_hypo_8= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        peak_glucose_all_m_hypo_8 = getfield(peak_glucose_all_m_hypo_8,'peak_glucose_all_m_hypo');
        PD_m_hypo_events_8=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        PD_m_hypo_events_8 = getfield(PD_m_hypo_events_8,'PD_m_hypo_events');

        subjectN=13;
        %HbT_all_m_hypo_gmsurface_all_13= load("F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        %HbT_all_m_hypo_gmsurface_all_13 = getfield(HbT_all_m_hypo_gmsurface_all_13,'HbT_all_m_hypo_gmsurface');        
        %event_jac_threshold_idx_13 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        %event_jac_threshold_idx_13 = getfield(event_jac_threshold_idx_13,'event_jac_threshold_idx');       
        %peak_glucose_all_m_hypo_13= load("F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        %peak_glucose_all_m_hypo_13 = getfield(peak_glucose_all_m_hypo_13,'peak_glucose_all_m_hypo');
        %PD_m_hypo_events_13=load("F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        %PD_m_hypo_events_13 = getfield(PD_m_hypo_events_13,'PD_m_hypo_events');   

        %Not A, Yes D, Yes F
        subjectN=15;
        HbT_all_m_hypo_gmsurface_all_15= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        HbT_all_m_hypo_gmsurface_all_15 = getfield(HbT_all_m_hypo_gmsurface_all_15,'HbT_all_m_hypo_gmsurface');        
        %event_jac_threshold_idx_15 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        %event_jac_threshold_idx_15 = getfield(event_jac_threshold_idx_15,'event_jac_threshold_idx');       
        peak_glucose_all_m_hypo_15= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        peak_glucose_all_m_hypo_15 = getfield(peak_glucose_all_m_hypo_15,'peak_glucose_all_m_hypo');
        PD_m_hypo_events_15=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        PD_m_hypo_events_15 = getfield(PD_m_hypo_events_15,'PD_m_hypo_events');

        subjectN=20;
        HbT_all_m_hypo_gmsurface_all_20= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        HbT_all_m_hypo_gmsurface_all_20 = getfield(HbT_all_m_hypo_gmsurface_all_20,'HbT_all_m_hypo_gmsurface');   
        event_jac_threshold_idx_20 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_20 = getfield(event_jac_threshold_idx_20,'event_jac_threshold_idx');
        peak_glucose_all_m_hypo_20= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        peak_glucose_all_m_hypo_20 = getfield(peak_glucose_all_m_hypo_20,'peak_glucose_all_m_hypo');
        PD_m_hypo_events_20=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        PD_m_hypo_events_20 = getfield(PD_m_hypo_events_20,'PD_m_hypo_events');

        

        subjectN=25;
        HbT_all_m_hypo_gmsurface_all_25= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        HbT_all_m_hypo_gmsurface_all_25 = getfield(HbT_all_m_hypo_gmsurface_all_25,'HbT_all_m_hypo_gmsurface');        
        event_jac_threshold_idx_25 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_25 = getfield(event_jac_threshold_idx_25,'event_jac_threshold_idx');       
        peak_glucose_all_m_hypo_25= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        peak_glucose_all_m_hypo_25 = getfield(peak_glucose_all_m_hypo_25,'peak_glucose_all_m_hypo');
        PD_m_hypo_events_25=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        PD_m_hypo_events_25 = getfield(PD_m_hypo_events_25,'PD_m_hypo_events');

        
        subjectN=33;
        HbT_all_m_hypo_gmsurface_all_33= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        HbT_all_m_hypo_gmsurface_all_33 = getfield(HbT_all_m_hypo_gmsurface_all_33,'HbT_all_m_hypo_gmsurface');       
        event_jac_threshold_idx_33 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_33 = getfield(event_jac_threshold_idx_33,'event_jac_threshold_idx'); 
        peak_glucose_all_m_hypo_33= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        peak_glucose_all_m_hypo_33 = getfield(peak_glucose_all_m_hypo_33,'peak_glucose_all_m_hypo');
        PD_m_hypo_events_33=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        PD_m_hypo_events_33 = getfield(PD_m_hypo_events_33,'PD_m_hypo_events');

       
        subjectN=49;
        HbT_all_m_hypo_gmsurface_all_49= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        HbT_all_m_hypo_gmsurface_all_49 = getfield(HbT_all_m_hypo_gmsurface_all_49,'HbT_all_m_hypo_gmsurface');
        event_jac_threshold_idx_49 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_49 = getfield(event_jac_threshold_idx_49,'event_jac_threshold_idx'); 
        peak_glucose_all_m_hypo_49= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        peak_glucose_all_m_hypo_49 = getfield(peak_glucose_all_m_hypo_49,'peak_glucose_all_m_hypo');       
        PD_m_hypo_events_49=load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        PD_m_hypo_events_49 = getfield(PD_m_hypo_events_49,'PD_m_hypo_events');


        %size(HbT_all_m_hypo_gmsurface_all_7)

        %load("F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        %load("F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\length_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","length_glucose_all_m_hypo");
        %load("F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
        %load("F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");

%load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
            
            switch time_window
                case "A"
                
                case "D"
                   event_jac_threshold_idx_15 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(15)+"_"+time_window+".mat","event_jac_threshold_idx");
                   event_jac_threshold_idx_15 = getfield(event_jac_threshold_idx_15,'event_jac_threshold_idx');       
       
                   good_m_hypo_15=find(event_jac_threshold_idx_15(:,1)==1);
                   HbT_all_gmsurface_15 = HbT_all_m_hypo_gmsurface_all_15(:,:,good_m_hypo_15);
                   peak_glucose_all_15 = peak_glucose_all_m_hypo_15(good_m_hypo_15);

                   good_events_listidx_15 = good_m_hypo_15;
                   good_PD_events_15 = PD_m_hypo_events_15(good_events_listidx_15);
            case "F"
                    event_jac_threshold_idx_15 = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(15)+"_"+time_window+".mat","event_jac_threshold_idx");
                    event_jac_threshold_idx_15 = getfield(event_jac_threshold_idx_15,'event_jac_threshold_idx');       
       
                    good_m_hypo_15=find(event_jac_threshold_idx_15(:,1)==1);
                   HbT_all_gmsurface_15 = HbT_all_m_hypo_gmsurface_all_15(:,:,good_m_hypo_15);
                   peak_glucose_all_15 = peak_glucose_all_m_hypo_15(good_m_hypo_15);

                   good_events_listidx_15 = good_m_hypo_15;
                   good_PD_events_15 = PD_m_hypo_events_15(good_events_listidx_15);


            end
            good_m_hypo_7=find(event_jac_threshold_idx_7(:,1)==1);
            good_m_hypo_8=find(event_jac_threshold_idx_8(:,1)==1);
            %good_m_hypo_13=find(event_jac_threshold_idx_13(:,1)==1);
            %good_m_hypo_15=find(event_jac_threshold_idx_15(:,1)==1);
            good_m_hypo_20=find(event_jac_threshold_idx_20(:,1)==1);
            good_m_hypo_25=find(event_jac_threshold_idx_25(:,1)==1);
            good_m_hypo_33=find(event_jac_threshold_idx_33(:,1)==1);
            good_m_hypo_49=find(event_jac_threshold_idx_49(:,1)==1);

    
            HbT_all_gmsurface_7 = HbT_all_m_hypo_gmsurface_all_7(:,:,good_m_hypo_7);
            HbT_all_gmsurface_8 = HbT_all_m_hypo_gmsurface_all_8(:,:,good_m_hypo_8);
            %HbT_all_gmsurface_13 = HbT_all_m_hypo_gmsurface_all_13(:,:,good_m_hypo_13);
            %HbT_all_gmsurface_15 = HbT_all_m_hypo_gmsurface_all_15(:,:,good_m_hypo_15);
            HbT_all_gmsurface_20 = HbT_all_m_hypo_gmsurface_all_20(:,:,good_m_hypo_20);
            HbT_all_gmsurface_25 = HbT_all_m_hypo_gmsurface_all_25(:,:,good_m_hypo_25);
            HbT_all_gmsurface_33 = HbT_all_m_hypo_gmsurface_all_33(:,:,good_m_hypo_33);
            HbT_all_gmsurface_49 = HbT_all_m_hypo_gmsurface_all_49(:,:,good_m_hypo_49);

            peak_glucose_all_7 = peak_glucose_all_m_hypo_7(good_m_hypo_7);
            peak_glucose_all_8 = peak_glucose_all_m_hypo_8(good_m_hypo_8);
            %peak_glucose_all_13 = peak_glucose_all_m_hypo_13(good_m_hypo_13);
            %peak_glucose_all_15 = peak_glucose_all_m_hypo_15(good_m_hypo_15);
            peak_glucose_all_20 = peak_glucose_all_m_hypo_20(good_m_hypo_20);
            peak_glucose_all_25 = peak_glucose_all_m_hypo_25(good_m_hypo_25);
            peak_glucose_all_33 = peak_glucose_all_m_hypo_33(good_m_hypo_33);
            peak_glucose_all_49 = peak_glucose_all_m_hypo_49(good_m_hypo_49);


good_events_listidx_7 = good_m_hypo_7;
            good_PD_events_7 = PD_m_hypo_events_7(good_events_listidx_7);
good_events_listidx_8 = good_m_hypo_8;
            good_PD_events_8 = PD_m_hypo_events_8(good_events_listidx_8);
%good_events_listidx_13 = good_m_hypo_13;
%            good_PD_events_13 = PD_m_hypo_events_13(good_events_listidx_13);
%good_events_listidx_15 = good_m_hypo_15;
%            good_PD_events_15 = PD_m_hypo_events_15(good_events_listidx_15);
good_events_listidx_20 = good_m_hypo_20;
            good_PD_events_20 = PD_m_hypo_events_20(good_events_listidx_20);
good_events_listidx_25 = good_m_hypo_25;
            good_PD_events_25 = PD_m_hypo_events_25(good_events_listidx_25);
good_events_listidx_33 = good_m_hypo_33;
            good_PD_events_33 = PD_m_hypo_events_33(good_events_listidx_33);
good_events_listidx_49 = good_m_hypo_49;
            good_PD_events_49 = PD_m_hypo_events_49(good_events_listidx_49);



% stack variable for
%add 13 and 15n(not forA)
switch time_window
    case "A"
    HbT_all_gmsurface = cat(3,HbT_all_gmsurface_7,HbT_all_gmsurface_8,HbT_all_gmsurface_20,HbT_all_gmsurface_25,HbT_all_gmsurface_33,HbT_all_gmsurface_49);
    peak_glucose_all = [peak_glucose_all_7;peak_glucose_all_8;peak_glucose_all_20;peak_glucose_all_25;peak_glucose_all_33;peak_glucose_all_49];
    good_PD_events = [good_PD_events_7 good_PD_events_8 good_PD_events_20 good_PD_events_25 good_PD_events_33 good_PD_events_49];
    good_PD_events_subN = [ones(1,size(good_PD_events_7,2))*7 ones(1,size(good_PD_events_8,2))*8 ones(1,size(good_PD_events_20,2))*20 ones(1,size(good_PD_events_25,2))*25 ones(1,size(good_PD_events_33,2))*33 ones(1,size(good_PD_events_49,2))*49];


    case "D"
    HbT_all_gmsurface = cat(3,HbT_all_gmsurface_7,HbT_all_gmsurface_8,HbT_all_gmsurface_15,HbT_all_gmsurface_20,HbT_all_gmsurface_25,HbT_all_gmsurface_33,HbT_all_gmsurface_49);
    peak_glucose_all = [peak_glucose_all_7;peak_glucose_all_8;peak_glucose_all_15;peak_glucose_all_20;peak_glucose_all_25;peak_glucose_all_33;peak_glucose_all_49];
    good_PD_events = [good_PD_events_7 good_PD_events_8 good_PD_events_15 good_PD_events_20 good_PD_events_25 good_PD_events_33 good_PD_events_49];
    good_PD_events_subN = [ones(1,size(good_PD_events_7,2))*7 ones(1,size(good_PD_events_8,2))*8 ones(1,size(good_PD_events_15,2))*15 ones(1,size(good_PD_events_20,2))*20 ones(1,size(good_PD_events_25,2))*25 ones(1,size(good_PD_events_33,2))*33 ones(1,size(good_PD_events_49,2))*49];


    case "F"
    HbT_all_gmsurface = cat(3,HbT_all_gmsurface_7,HbT_all_gmsurface_8,HbT_all_gmsurface_15,HbT_all_gmsurface_20,HbT_all_gmsurface_25,HbT_all_gmsurface_33,HbT_all_gmsurface_49);
    peak_glucose_all = [peak_glucose_all_7;peak_glucose_all_8;peak_glucose_all_15;peak_glucose_all_20;peak_glucose_all_25;peak_glucose_all_33;peak_glucose_all_49];
    good_PD_events = [good_PD_events_7 good_PD_events_8 good_PD_events_15 good_PD_events_20 good_PD_events_25 good_PD_events_33 good_PD_events_49];
    good_PD_events_subN = [ones(1,size(good_PD_events_7,2))*7 ones(1,size(good_PD_events_8,2))*8 ones(1,size(good_PD_events_15,2))*15 ones(1,size(good_PD_events_20,2))*20 ones(1,size(good_PD_events_25,2))*25 ones(1,size(good_PD_events_33,2))*33 ones(1,size(good_PD_events_49,2))*49];

end
%% SKIP THIS Run for single subject - mild and severe combined

        eventType = "m_hypo";
        HbT_all_m_hypo_gmsurface_all= load("F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_m_hypo_gmsurface");
        HbT_all_m_hypo_gmsurface_all = getfield(HbT_all_m_hypo_gmsurface_all,'HbT_all_m_hypo_gmsurface');        
        event_jac_threshold_idx_m_hypo = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_m_hypo = getfield(event_jac_threshold_idx_m_hypo,'event_jac_threshold_idx');       
        peak_glucose_all_m_hypo= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_m_hypo");
        peak_glucose_all_m_hypo = getfield(peak_glucose_all_m_hypo,'peak_glucose_all_m_hypo');
        PD_m_hypo_events= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_m_hypo_events");
        PD_m_hypo_events = getfield(PD_m_hypo_events,'PD_m_hypo_events');

        eventType = "S_hypo";
        HbT_all_S_hypo_gmsurface_all= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\HbT_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","HbT_all_S_hypo_gmsurface");
        HbT_all_S_hypo_gmsurface_all = getfield(HbT_all_S_hypo_gmsurface_all,'HbT_all_S_hypo_gmsurface');        
        event_jac_threshold_idx_S_hypo = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
        event_jac_threshold_idx_S_hypo = getfield(event_jac_threshold_idx_S_hypo,'event_jac_threshold_idx');       
        peak_glucose_all_S_hypo= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\peak_glucose_all_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","peak_glucose_all_S_hypo");
        peak_glucose_all_S_hypo = getfield(peak_glucose_all_S_hypo,'peak_glucose_all_S_hypo');
        PD_S_hypo_events= load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\PD_good_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","PD_S_hypo_events");
        PD_S_hypo_events = getfield(PD_S_hypo_events,'PD_S_hypo_events');

good_m_hypo=find(event_jac_threshold_idx_m_hypo(:,1)==1);
good_S_hypo=find(event_jac_threshold_idx_S_hypo(:,1)==1);

HbT_all_gmsurface_m_hypo = HbT_all_m_hypo_gmsurface_all(:,:,good_m_hypo);
HbT_all_gmsurface_S_hypo = HbT_all_S_hypo_gmsurface_all(:,:,good_S_hypo);

peak_glucose_all_m_hypo = peak_glucose_all_m_hypo(good_m_hypo);
peak_glucose_all_S_hypo = peak_glucose_all_S_hypo(good_S_hypo);

    
good_events_listidx_m_hypo = good_m_hypo;
good_events_listidx_S_hypo = good_S_hypo;

good_PD_events_m_hypo = PD_m_hypo_events(good_events_listidx_m_hypo);
good_PD_events_S_hypo = PD_S_hypo_events(good_events_listidx_S_hypo);

 HbT_all_gmsurface = cat(3,HbT_all_gmsurface_m_hypo,HbT_all_gmsurface_S_hypo);
    peak_glucose_all = [peak_glucose_all_m_hypo;peak_glucose_all_S_hypo]; 
    good_PD_events = [good_PD_events_m_hypo good_PD_events_S_hypo];

%% skip run this - everytime - for individual subjects - ind. event types
load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");

% %just using 'good' events
switch subjectN
    case 7
        switch eventType
            case "m_hypo"
            good_m_hypo=find(event_jac_threshold_idx(:,1)==1);
            HbT_all_gmsurface = HbT_all_m_hypo_gmsurface(:,:,good_m_hypo);
            peak_glucose_all = peak_glucose_all_m_hypo(good_m_hypo);
        end
    case 8
        switch eventType
            case "m_hypo"
            good_m_hypo=find(event_jac_threshold_idx(:,1)==1);
            HbT_all_gmsurface = HbT_all_m_hypo_gmsurface(:,:,good_m_hypo);
            peak_glucose_all = peak_glucose_all_m_hypo(good_m_hypo);
        end
    case 20
        switch eventType
            case "m_hypo"
            good_m_hypo=find(event_jac_threshold_idx(:,1)==1);
            HbT_all_gmsurface = HbT_all_m_hypo_gmsurface(:,:,good_m_hypo);
            peak_glucose_all = peak_glucose_all_m_hypo(good_m_hypo);
        end
    case 25
        switch eventType
            case "m_hypo"
            good_m_hypo=find(event_jac_threshold_idx(:,1)==1);
            HbT_all_gmsurface = HbT_all_m_hypo_gmsurface(:,:,good_m_hypo);
            peak_glucose_all = peak_glucose_all_m_hypo(good_m_hypo);
        end
    case 33
        switch eventType
            case "m_hypo"
            good_m_hypo=find(event_jac_threshold_idx(:,1)==1);
            HbT_all_gmsurface = HbT_all_m_hypo_gmsurface(:,:,good_m_hypo);
            peak_glucose_all = peak_glucose_all_m_hypo(good_m_hypo);
        end
    case 49
        switch eventType
            case "m_hypo"
            good_m_hypo=find(event_jac_threshold_idx(:,1)==1);
            HbT_all_gmsurface = HbT_all_m_hypo_gmsurface(:,:,good_m_hypo);
            peak_glucose_all = peak_glucose_all_m_hypo(good_m_hypo);
            case "S_hypo"
            good_S_hypo=find(event_jac_threshold_idx(:,1)==1);
            HbT_all_gmsurface = HbT_all_S_hypo_gmsurface(:,:,good_S_hypo);
            peak_glucose_all = peak_glucose_all_S_hypo(good_S_hypo);
        end
    case 60
        switch eventType
            case "m_hypo"
            good_m_hypo=find(event_jac_threshold_idx(:,1)==1);
            HbT_all_gmsurface = HbT_all_m_hypo_gmsurface(:,:,good_m_hypo);
            peak_glucose_all = peak_glucose_all_m_hypo(good_m_hypo);
        end
end


%set up like PD20

switch eventType
    case "m_hypo"
    good_events_listidx = good_m_hypo;
    good_PD_events = PD_m_hypo_events(good_events_listidx);
    eventType_subN_idx = ones(size(good_PD_events,2),1); %1=mhypo,2=Shypo,3=mhyper,4=Shyper


    case "S_hypo"
    good_events_listidx = good_S_hypo;
    good_PD_events = PD_S_hypo_events(good_events_listidx);
    eventType_subN_idx = ones(size(good_PD_events,2),1)*2; %1=mhypo,2=Shypo,3=mhyper,4=Shyper


end


good_PD_events_subN = ones(size(event_good,2),1)*subjectN;
%update eventtype Sub N string
eventType_subN = strings(size(good_PD_events,2),1);
eventType_subN(:,1) = eventType;

eventType_subN_idx = ones(size(good_PD_events,2),1); %1=mhypo,2=Shypo,3=mhyper,4=Shyper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rest of code - MISC.

%%%%%%%%%%%%%%%%%%% end
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
for i=1:size(HbT_metrics_PD_events.EEG_9_ROI_HbT_all_gmsurface,1)
    subplot(3,3,i)
    plot(BGC_metrics_PD_events.local_BCG_mean,HbT_metrics_PD_events.EEG_9_ROI_HbT_all_gmsurface(i,:),'rx')
    ylim([-10 5])
end

























%% Metrics SB 19 02 25 meeting
clc;
subjectN=45;
event_jac_threshold_idx = event_jac_threshold_idx_45;
%load("event_stats_m_hypo_PD"+num2str(subjectN)+"_A")
load("event_stats_"+eventType+"_PD"+num2str(subjectN)+"_F")
%Q1
x = good_PD_events(1,good_PD_events_subN==subjectN);
pc_Tpts_MA_train_allE = zeros(size(x,2),1);
N_MA_trains_allE = zeros(size(x,2),1);
for i=1:size(x,2)
    %events_stats.eventN_allE  == x(i)
    pc_Tpts_MA_train_allE(i) = events_stats.pc_Tpts_MA_train_allE(events_stats.eventN_allE  == x(i));
    N_MA_trains_allE(i) = events_stats.N_MA_trains_allE(events_stats.eventN_allE  == x(i));
end
pc_Tpts_MA_train_allE
N_MA_trains_allE
%Q2 
%change subjectN
%subjectN=8; 
%load("event_stats_m_hypo_PD"+num2str(subjectN)+"_A")
load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");

x = good_PD_events(1,good_PD_events_subN==subjectN);
goodch_signal_prune = zeros(size(x,2),1);
for i=1:size(good_PD_events(1,good_PD_events_subN==subjectN),2)
    goodch_signal_prune(i) = nnz(events_stats.goodch_idx_allE(:,events_stats.eventN_allE  == x(i)));
end
goodch_signal_prune = sum(goodch_signal_prune);
goodch_after_rm_outliers = nnz(good_ch_idx_events(:,event_jac_threshold_idx(:,1)==1)); %change eventhjac idx

[goodch_signal_prune goodch_after_rm_outliers]

%% Figures for Paper
% figure; plotmesh_iso2([scalpSurfaceMesh.node ones(size(scalpSurfaceMesh.node,1),1) ],scalpSurfaceMesh.face);view(0,90);
% hold on
% plot3(mesh.source.coord(:,1), mesh.source.coord(:,2),mesh.source.coord(:,3),'ro','LineWidth',16))
% cb = colorbar('horiz');
% title("Jacobian - GM Smooth Week30")
% xlabel("x / mm");
% ylabel("y / mm");
% zlabel("z / mm");
% 
% figure; plotmesh_iso2([mesh.nodes ones(size(mesh.nodes,1),1) ],mesh.elements  );view(0,90);
% cb = colorbar('horiz');
% hold on
% plot3(mesh.source.coord(:,1), mesh.source.coord(:,2),mesh.source.coord(:,3),'ro','LineWidth',16,'MarkerSize',16)

fontsize_n = 20 ; 
%Source and detectors over gmsmooth surface
figure; plotmesh_iso2([gmSurfaceMesh.node ones(size(gmSurfaceMesh.node,1),1)*1 ],gmSurfaceMesh.face);view(0,90);
hold on
plot3(mesh.source.coord(:,1), mesh.source.coord(:,2),mesh.source.coord(:,3),'ro','LineWidth',12)
plot3(mesh.meas.coord(:,1), mesh.meas.coord(:,2),mesh.meas.coord(:,3),'bx','LineWidth',12)
cb = colorbar('horiz');
clim([0 2])
title("Mesh- GM Smooth Week30")
xlabel("x / mm",FontSize=fontsize_n);
ylabel("y / mm",FontSize=fontsize_n);
zlabel("z / mm",FontSize=fontsize_n);
set(gca,'fontsize',fontsize_n)

figure; 
subplot(1,2,1)
plotmesh_iso2([gmSurfaceMesh.node ones(size(gmSurfaceMesh.node,1),1)*1 ],gmSurfaceMesh.face);view(-90,0);
hold on
plot3(mesh.source.coord(:,1), mesh.source.coord(:,2),mesh.source.coord(:,3),'ro','LineWidth',12)
plot3(mesh.meas.coord(:,1), mesh.meas.coord(:,2),mesh.meas.coord(:,3),'bx','LineWidth',12)
cb = colorbar('horiz');
clim([0 2])
xlabel("x / mm",FontSize=fontsize_n);
ylabel("y / mm",FontSize=fontsize_n);
zlabel("z / mm",FontSize=fontsize_n);
subplot(1,2,2)
plotmesh_iso2([gmSurfaceMesh.node ones(size(gmSurfaceMesh.node,1),1)*1 ],gmSurfaceMesh.face);view(90,0);
hold on
plot3(mesh.source.coord(:,1), mesh.source.coord(:,2),mesh.source.coord(:,3),'ro','LineWidth',12)
plot3(mesh.meas.coord(:,1), mesh.meas.coord(:,2),mesh.meas.coord(:,3),'bx','LineWidth',12)
cb = colorbar('horiz');
clim([0 2])
xlabel("x / mm",FontSize=fontsize_n);
ylabel("y / mm",FontSize=fontsize_n);
zlabel("z / mm",FontSize=fontsize_n);
set(gca,'fontsize',fontsize_n)


figure; 
subplot(1,2,1)
plotmesh_iso2([gmSurfaceMesh.node  masks.mask_j_thresh_gmsurface ],gmSurfaceMesh.face);view(-90,0);
hold on
plot3(mesh.source.coord(:,1), mesh.source.coord(:,2),mesh.source.coord(:,3),'ro','LineWidth',24)
plot3(mesh.meas.coord(:,1), mesh.meas.coord(:,2),mesh.meas.coord(:,3),'bx','LineWidth',24)
cb = colorbar('horiz');
clim([0 2])
xlabel("x / mm",FontSize=fontsize_n);
ylabel("y / mm",FontSize=fontsize_n);
zlabel("z / mm",FontSize=fontsize_n);
set(gca,'fontsize',fontsize_n)
ax = gca;
ax.FontSize = 20;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);
set(gca,'XColor', 'none','YColor','none','ZColor','none')
set(gca, 'color', 'none');
subplot(1,2,2)
plotmesh_iso2([gmSurfaceMesh.node  masks.mask_j_thresh_gmsurface ],gmSurfaceMesh.face);view(90,0);
hold on
plot3(mesh.source.coord(:,1), mesh.source.coord(:,2),mesh.source.coord(:,3),'ro','LineWidth',24)
plot3(mesh.meas.coord(:,1), mesh.meas.coord(:,2),mesh.meas.coord(:,3),'bx','LineWidth',24)
cb = colorbar('horiz');
clim([0 2])
xlabel("x / mm",FontSize=fontsize_n);
ylabel("y / mm",FontSize=fontsize_n);
zlabel("z / mm",FontSize=fontsize_n);
set(gca,'fontsize',fontsize_n)
ax = gca;
ax.FontSize = 20;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'ZTickLabel',[]);
set(gca,'XColor', 'none','YColor','none','ZColor','none')
set(gca, 'color', 'none');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%figure - channel map % of good channels for each channel - colour coded -
%2D?
switch eventType 
    case "m_hypo"
        if mild_and_severe == 1 %m and s hypo
            all_subject_N = [8 10 15 25 33 34 39 41 44 45 48 49 55 56 58 59 60 15 25 33 49];
                                %mhypo and then the shypo subjects
        else %just mhypo
            switch time_window
                case "A"
                    all_subject_N = [8 10 15 25 33 39 41 44 45 48 49 55 56 59 60];
                case "D"
                    all_subject_N = [8 10 15 25 33 34 39 41 44 45 49 55 56 58 59 60];
                case "F"
                    all_subject_N = [8 10 15 25 33 39 41 44 45 49 55 56 58 59 60];
            end
        end
    case "S_hypo"
        all_subject_N = [15 25 33 49];
end

%load good_ch_idx_events
good_ch_idx_events_all_sub_N = zeros(128,50);
col = 1;
for i=1:size(all_subject_N,2)
    subjectN = all_subject_N(i);
    load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
    event_jac_threshold_idx_N = load("event_jac_thresh_idx_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+".mat","event_jac_threshold_idx");
    event_jac_threshold_idx_N = getfield(event_jac_threshold_idx_N,'event_jac_threshold_idx');
    good_ch_idx_events_all_sub_N(:,col:col +sum(event_jac_threshold_idx_N(:,1)==1)-1  ) = good_ch_idx_events(:,event_jac_threshold_idx_N(:,1)==1) ;
    col = col + sum(event_jac_threshold_idx_N(:,1)==1);
end

N_good_ch_idx_events_all_sub_N = zeros(1,128);

for j=1:128
    N_good_ch_idx_events_all_sub_N(1,j) = size(find(good_ch_idx_events_all_sub_N==j),1);
end

%make sds
sds = zeros(1,64);
for i=1:64
    %(mesh.source.coord(mesh.link(i,1),1), mesh.source.coord(mesh.link(i,1),2),mesh.source.coord(mesh.link(i,1),3) )
    %(mesh.meas.coord(mesh.link(i,2),1), mesh.meas.coord(mesh.link(i,2),2),mesh.meas.coord(mesh.link(i,2),3) )
    sds(i) = sqrt(  (mesh.source.coord(mesh.link(i,1),1)-mesh.meas.coord(mesh.link(i,2),1))^2 + (mesh.source.coord(mesh.link(i,1),2)-mesh.meas.coord(mesh.link(i,2),2))^2  + (mesh.source.coord(mesh.link(i,1),3)-mesh.meas.coord(mesh.link(i,2),3))^2    );
end
%this is a generic SDS, irl SDSs may change due to different head sizes
%etc. also young babies will be smaller and have more transpatent heads

%do plot
% figure()
% hold on
% for i=1:8
%     scatter(mesh.source.coord(mesh.link(i,1),1),mesh.source.coord(mesh.link(i,1),2),'ro','linewidth',8)
%     scatter(mesh.meas.coord(mesh.link(i,1),1),mesh.meas.coord(mesh.link(i,1),2),'bx','linewidth',8)
%     plot([mesh.source.coord(mesh.link(i,1),1) mesh.meas.coord(mesh.link(i,1),1)], [mesh.source.coord(mesh.link(i,1),2) mesh.meas.coord(mesh.link(i,1),2)],'LineWidth',6)
% end

% figure()
% 
%     scatter(mesh.source.coord(mesh.link(:,1),1),mesh.source.coord(mesh.link(:,1),2),'ro','linewidth',8)
% hold on
%     scatter(mesh.meas.coord(mesh.link(:,1),1),mesh.meas.coord(mesh.link(:,1),2),'bx','linewidth',8)
%     plot([mesh.source.coord(mesh.link(:,1),1) mesh.meas.coord(mesh.link(:,2),1)], [mesh.source.coord(mesh.link(:,1),2) mesh.meas.coord(mesh.link(:,2),2)],'LineWidth',6)

greyJetSB = [0	0	0.562500000000000;
0	0.0208333339542151	0.708333313465118;
0	0.0416666679084301	0.854166686534882;
0	0.0625000000000000	1;
0	0.129464283585548	1;
0	0.196428567171097	1;
0	0.263392865657806	1;
0	0.330357134342194	1;
0	0.397321432828903	1;
0	0.464285701513290	1;
0	0.531250000000000	1;
0	0.598214268684387	1;
0	0.665178596973419	1;
0	0.732142865657806	1;
0	0.799107134342194	1;
0	0.866071403026581	1;
0	0.933035731315613	1;
0	1	1;
0.0593837536871433	0.986834764480591	0.984593868255615;
0.118767507374287	0.973669469356537	0.969187676906586;
0.178151264786720	0.960504174232483	0.953781485557556;
0.237535014748573	0.947338938713074	0.938375353813171;
0.296918779611588	0.934173703193665	0.922969222068787;
0.356302529573441	0.921008408069611	0.907563030719757;
0.415686279535294	0.907843112945557	0.892156839370728;
0.475070029497147	0.894677877426148	0.876750707626343;
0.534453809261322	0.881512641906738	0.861344575881958;
0.593837559223175	0.868347346782684	0.845938384532929;
0.653221309185028	0.855182051658630	0.830532193183899;
0.712605059146881	0.842016816139221	0.815126061439514;
0.771988809108734	0.828851580619812	0.799719929695129;
0.831372559070587	0.815686285495758	0.784313738346100;
0.831372559070587	0.815686285495758	0.784313738346100;
0.844343900680542	0.829864263534546	0.723981916904450;
0.857315242290497	0.844042241573334	0.663650095462799;
0.870286583900452	0.858220219612122	0.603318274021149;
0.883257925510407	0.872398197650909	0.542986452579498;
0.896229267120361	0.886576175689697	0.482654601335526;
0.909200608730316	0.900754153728485	0.422322779893875;
0.922171950340271	0.914932131767273	0.361990958452225;
0.935143291950226	0.929110109806061	0.301659137010574;
0.948114633560181	0.943288087844849	0.241327300667763;
0.961085975170136	0.957466065883637	0.180995479226112;
0.974057316780090	0.971644043922424	0.120663650333881;
0.987028658390045	0.985822021961212	0.0603318251669407;
1	1	0;
1	0.933333337306976	0;
1	0.866666674613953	0;
1	0.800000011920929	0;
1	0.733333349227905	0;
1	0.666666686534882	0;
1	0.600000023841858	0;
1	0.533333361148834	0;
1	0.466666668653488	0;
1	0.400000005960465	0;
1	0.333333343267441	0;
1	0.266666680574417	0;
1	0.200000002980232	0;
1	0.133333340287209	0;
1	0.0666666701436043	0;
1	0	0;
0.833333313465118	0	0;
0.666666686534882	0	0;
0.500000000000000	0	0];

greyJetSB_48col = greyJetSB(end-47:end,:); %48 colours from something to red
greyJetSB_50col = greyJetSB(end-49:end,:); %48 colours from something to red

greyJetSBJac = greyJetSB(end/2:end,:); %just grey to red


N_good_ch_idx_events_all_sub_N; %values can be max of 48 (as there are 48 events for time D - all m hypo, no hyper)
%first 1-64 are wvl 1
%second 65-128 are wv 2
%and they repeat the same

landmark_EEG1020 = [
	-23.57	14.32	24.47;
	-1.09	15.9	37.04;
 	21.63	12.93	26.08;
 	-32.89	-16.63	35.48;
 	-2.01	-17.37	51.07;
 	29.67	-16.27	36.75;
 	-24.76	-50.57	27.24;
 	-1.33	-52.76	39.53;
 	22.94	-49.58	29.2;
    -1.59 	 0.86 	 46.78;
    -2.05 	 -36.67 	 48.65];


% %FCz
% -1.59 	 0.86 	 46.78 
% 
% %Cpz
% -2.05 	 -36.67 	 48.65 
%%

figure()
plot([mesh.source.coord(mesh.link(1,1),1) mesh.meas.coord(mesh.link(1,2),1)], [mesh.source.coord(mesh.link(1,1),2) mesh.meas.coord(mesh.link(1,2),2)],'LineWidth',4,'Color', greyJetSBJac(N_good_ch_idx_events_all_sub_N(1),:)   )
hold on
for i=2:64
    plot([mesh.source.coord(mesh.link(i,1),1) mesh.meas.coord(mesh.link(i,2),1)], [mesh.source.coord(mesh.link(i,1),2) mesh.meas.coord(mesh.link(i,2),2)],'LineWidth',4,'Color',greyJetSBJac(N_good_ch_idx_events_all_sub_N(i)-1,:) )
end
scatter(mesh.source.coord(2,1),mesh.source.coord(2,2),'o','linewidth',34,'MarkerEdgeColor',[253/256 152/256 64/256])
scatter(mesh.source.coord(5,1),mesh.source.coord(5,2),'o','linewidth',34,'MarkerEdgeColor',[253/256 152/256 64/256])
%scatter(mesh.source.coord(mesh.link(:,1),1),mesh.source.coord(mesh.link(:,1),2),'ro','linewidth',24)
%hold on
%scatter(mesh.meas.coord(mesh.link(:,1),1),mesh.meas.coord(mesh.link(:,1),2),'bx','linewidth',24)
%text(landmark_EEG1020(end-1,1)-1,landmark_EEG1020(end-1,2)+7,'FCz','FontSize',20) %text
%text(landmark_EEG1020(end,1)-1,landmark_EEG1020(end,2)-4,'CPz','FontSize',20) %text
%scatter(mesh.source.coord(mesh.link(5,1),1),mesh.source.coord(mesh.link(5,1),2),'ro','linewidth',18)
scatter(landmark_EEG1020(4:6,1),landmark_EEG1020(4:6,2),'o','linewidth',14,'MarkerEdgeColor',[253/256 152/256 64/256])
scatter(mesh.meas.coord(mesh.link(:,1),1),mesh.meas.coord(mesh.link(:,1),2),'bx','linewidth',24) %detectors
scatter(mesh.source.coord(mesh.link(:,1),1),mesh.source.coord(mesh.link(:,1),2),'ro','linewidth',24) %sources
%scatter(landmark_EEG1020(6,1),landmark_EEG1020(6,2),'o','linewidth',12,'MarkerEdgeColor',[117/256 84/256 174/256])
%xlabel("x / mm");
%ylabel("y / mm");
colormap(greyJetSBJac)
cb = colorbar('vertic');
caxis([1 34])
ax = gca;
ax.FontSize = 20;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
%%

%[253/256 152/256 64/256] orange
%[117/256 84/256 174/256] purple

figure()
scatter3(mesh.source.coord(mesh.link(:,1),1),mesh.source.coord(mesh.link(:,1),2),mesh.source.coord(mesh.link(:,1),3),'ro','linewidth',18)
hold on
scatter3(mesh.meas.coord(mesh.link(:,1),1),mesh.meas.coord(mesh.link(:,1),2),mesh.meas.coord(mesh.link(:,1),3),'bx','linewidth',18)
plot3([mesh.source.coord(mesh.link(1,1),1) mesh.meas.coord(mesh.link(1,2),1)], [mesh.source.coord(mesh.link(1,1),2) mesh.meas.coord(mesh.link(1,2),2)],[mesh.source.coord(mesh.link(1,1),3) mesh.meas.coord(mesh.link(1,2),3)],'LineWidth',2,'Color', greyJetSBJac(N_good_ch_idx_events_all_sub_N(1),:)   )
for i=2:64
    plot3([mesh.source.coord(mesh.link(i,1),1) mesh.meas.coord(mesh.link(i,2),1)], [mesh.source.coord(mesh.link(i,1),2) mesh.meas.coord(mesh.link(i,2),2)],[mesh.source.coord(mesh.link(i,1),3) mesh.meas.coord(mesh.link(i,2),3)],'LineWidth',2,'Color',greyJetSBJac(N_good_ch_idx_events_all_sub_N(i),:) )
end
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm")
colormap(greyJetSBJac)
cb = colorbar('vertic');
caxis([0 33])
ax = gca;
ax.FontSize = 20;


median(N_good_ch_idx_events_all_sub_N);


ch_above_median_use = find(N_good_ch_idx_events_all_sub_N>=median(N_good_ch_idx_events_all_sub_N));
ch_above_median_use = ch_above_median_use(1:end/2);


%find(N_good_ch_idx_events_all_sub_N>24)
%%
% RUN THIS Jacobian mapping
%map J.complete (Intensity J) to gmsurface
J_gmsurface = vol2gm*J.complete';
%sum this across all 64 channels
J_median_gmsurface = sum(J_gmsurface(:,ch_above_median_use)')';
%map meshsupport to gmsurface
mesh_support_gmsurface = vol2gm*mesh.support(:,1);
% spatially normalise J all gmsurface 
J_median_gmsurface_norm=J_median_gmsurface./mesh_support_gmsurface;
%plot J all gmsurface

%%
figure();plotmesh_JAC([gmSurfaceMesh.node J_median_gmsurface],gmSurfaceMesh.face);view(0,90);
hold on
scatter3(mesh.source.coord(mesh.link(:,1),1),mesh.source.coord(mesh.link(:,1),2),mesh.source.coord(mesh.link(:,1),3),'ro','linewidth',24)
scatter3(mesh.meas.coord(mesh.link(:,1),1),mesh.meas.coord(mesh.link(:,1),2),mesh.meas.coord(mesh.link(:,1),3),'bx','linewidth',24)
%plot3([mesh.source.coord(mesh.link(ch_above_median_use(1),1),1) mesh.meas.coord(mesh.link(ch_above_median_use(1),2),1)], [mesh.source.coord(mesh.link(ch_above_median_use(1),1),2) mesh.meas.coord(mesh.link(ch_above_median_use(1),2),2)],[mesh.source.coord(mesh.link(ch_above_median_use(1),1),3) mesh.meas.coord(mesh.link(ch_above_median_use(1),2),3)],'LineWidth',2  )
%for i=2:size(ch_above_median_use,2)
%    plot3([mesh.source.coord(mesh.link(ch_above_median_use(i),1),1) mesh.meas.coord(mesh.link(ch_above_median_use(i),2),1)], [mesh.source.coord(mesh.link(ch_above_median_use(i),1),2) mesh.meas.coord(mesh.link(ch_above_median_use(i),2),2)],[mesh.source.coord(mesh.link(ch_above_median_use(i),1),3) mesh.meas.coord(mesh.link(ch_above_median_use(i),2),3)],'LineWidth',2 )
%end
cb = colorbar('horiz');
title("Median Ch Jacobian - GM Smooth Week30 TW D all m hypo")
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm");
ax = gca;
ax.FontSize = 20;
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');


figure();sgtitle("Median Ch Jacobian - GM Smooth Week30 TW D all m hypo");
subplot(1,2,1);plotmesh_JAC([gmSurfaceMesh.node J_median_gmsurface],gmSurfaceMesh.face);view(0,90);
hold on
scatter3(mesh.source.coord(mesh.link(:,1),1),mesh.source.coord(mesh.link(:,1),2),mesh.source.coord(mesh.link(:,1),3),'ro','linewidth',24)
scatter3(mesh.meas.coord(mesh.link(:,1),1),mesh.meas.coord(mesh.link(:,1),2),mesh.meas.coord(mesh.link(:,1),3),'bx','linewidth',24)
%cb = colorbar('vertic');
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm");
ax = gca;
ax.FontSize = 20;
view([-90 0])
set(gca,'XColor', 'none','YColor','none','ZColor','none')
set(gca, 'color', 'none');
subplot(1,2,2);plotmesh_JAC([gmSurfaceMesh.node J_median_gmsurface],gmSurfaceMesh.face);view(0,90);
hold on
scatter3(mesh.source.coord(mesh.link(:,1),1),mesh.source.coord(mesh.link(:,1),2),mesh.source.coord(mesh.link(:,1),3),'ro','linewidth',24)
scatter3(mesh.meas.coord(mesh.link(:,1),1),mesh.meas.coord(mesh.link(:,1),2),mesh.meas.coord(mesh.link(:,1),3),'bx','linewidth',24)
%plot([mesh.source.coord(mesh.link(ch_above_median_use(1),1),1) mesh.meas.coord(mesh.link(ch_above_median_use(1),2),1)], [mesh.source.coord(mesh.link(ch_above_median_use(1),1),2) mesh.meas.coord(mesh.link(ch_above_median_use(1),2),2)],'LineWidth',2  )
%for i=2:size(ch_above_median_use,2)
%    plot([mesh.source.coord(mesh.link(ch_above_median_use(i),1),1) mesh.meas.coord(mesh.link(ch_above_median_use(i),2),1)], [mesh.source.coord(mesh.link(ch_above_median_use(i),1),2) mesh.meas.coord(mesh.link(ch_above_median_use(i),2),2)],'LineWidth',2 )
%end
%cb = colorbar('vertic');
xlabel("x / mm");
ylabel("y / mm");
zlabel("z / mm");
ax = gca;
ax.FontSize = 20;
view([90 0])
set(gca,'XColor', 'none','YColor','none','ZColor','none')
set(gca, 'color', 'none');




%%
figure()
hist(N_good_ch_idx_events_all_sub_N(1,1:end/2))
xlabel("N of good ch's (max 64)");
ylabel("N")

nnz(good_ch_idx_events_all_sub_N(:,1:48))/(48*2)


for i=1:48
    nnz(good_ch_idx_events_all_sub_N(:,i))
end




    %scatter(BCG_metric,max_HbT_p_val_0_05_R_all_pos(:,:),130,col_scat,'filled')

    %scatter(BCG_metric,max_HbT_p_val_0_05_R_all_pos(:,:),130,col_scat,'filled')







%% Glucose event - window A B C D E F G

%for i=1:size(good_PD_events_subN,2)
i=1;

    %load glucose data from subject N
    PD_N = load("F:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(good_PD_events_subN(i))+"");
    PD_N = getfield(PD_N,"PD"+num2str(good_PD_events_subN(i))+"");

    %get start and end time points for glucose event
    %PD_N.events_start_end.m_hypo(good_PD_events(i),:); %eventN = good_PD_events(i)

    %BCG values during entire event (GLOBAL BCG metrics)
    global_BCG = PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1) : PD_N.events_start_end.m_hypo(good_PD_events(i),2));
%end

time_window_t_BCG = zeros(7,2);

%"A"
            %time_window_t=[0 10];
            time_window_t_BCG(1,:)=[1 3]; %sample 1 to 3 covers first 10 mins 0,5,10 min
            %time_window_t_BCG=[2 4]; %sample 1 to 3 covers first 10 mins 0,5,10 min
%"B"
            time_window_t=[13 23]; %sample 1 to 3 covers first 10 mins 0,5,10 min
            time_window_t_BCG(2,:) =[3 5]; %sample 3 to 5 covers 10 to 20 mins 10,15,20 min
%"C" % 5Min before min glucose
% "m_hypo"
min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)))  );

    time_window_t_BCG(3,:) =[min_glucose_pos(end)-1 min_glucose_pos(end)+1]; %time sample of min glucose pos -1 (5mins)
                
% "D" % min glucose
% "m_hypo"
                    %find min glucose pos (index will corrospond to inddex in PD_time)
                    min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)))  );
                    %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                    %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                    %make it -2mins before point and 8mins after =10mins total
                    %time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];
                    time_window_t_BCG(4,:)=[min_glucose_pos(end) min_glucose_pos(end)+2];
                    %can change the (end) to (1) . (end) means it will always look at
                    %the LATEST glucoce MINIMUM, (1) means it looks at the FIRST
                    %glucose minimum
 %"E" % 5 min post min glucose
 %eventType
 %"m_hypo"
                    %find min glucose pos (index will corrospond to inddex in PD_time)
                    min_glucose_pos = find(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)) == min(PD_N.glucose(PD_N.events_start_end.m_hypo(good_PD_events(i),1):PD_N.events_start_end.m_hypo(good_PD_events(i),2)))  );
                    %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                    %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                    %make it -2mins before point and 8mins after =10mins total
                    %time_window_t=[PD_time(min_glucose_pos(end)+1 )-2 PD_time(min_glucose_pos(end)+1 )+8];
                    time_window_t_BCG(5,:) =[min_glucose_pos(end)+1 min_glucose_pos(end)+3]; %time sample of min glucose pos -1 (5mins)n

 %"F" %start of end baseline
            %time_window_t=[PD_time(end)-17 PD_time(end)-7];
            %incorrectF, similar to G
            %time_window_t_BCG =[size(global_BCG,1)-2 size(global_BCG,1)]; %last-2 lastlast time sample %time win F on 30 10 24 was this by mistake
            
            %correct F
            time_window_t_BCG(6,:) =[size(global_BCG,1)-3 size(global_BCG,1)-1]; %last-3 last-1 last time sample
            
            %F to end
            %time_window_t_BCG =[size(global_BCG,1)-3 size(global_BCG,1)]; %last-3 last last time sample

%case "G" %last point pre hypo glucose
%"G"
time_window_t_BCG(7,:) =[size(global_BCG,1)-2 size(global_BCG,1)];
    figure()
    for i=1:size(global_BCG,1)
    subplot(2,4,i)
    plot([0:1:size(global_BCG,1)-1]*5,global_BCG,'b-o')
    title("Time W"+time_window+" PD "+num2str(good_PD_events_subN(i))+" "+eventType+" "+num2str(good_PD_events(i))+". nirs")
    xlabel("Time / mins")
    ylabel("BGC / mg/DL")
    xline(15,'g--');xline((size(global_BCG,1)-4)*5,'r--');
    xline((time_window_t_BCG(i,1)-1)*5,'m--');xline((time_window_t_BCG(i,2)-1)*5,'m--');
    yline(72,'k--');
    end

%%
%Get HbT metrics
%Correlate between any two metrics
%info on stacked Sub events
edges = unique(good_PD_events_subN);
counts = histc(good_PD_events_subN(:), edges);[edges' counts]
size(edges,2)

%get_ROI_9_R_values
%% 
if mild_and_severe ==  1
    eventType= "S_m_hypo";
end

%for stacking subjects
if stacked == 1; %yes we are stacking
switch time_window
    case "A"
        subjectN = 820253349;
    case "D"
        subjectN = 81520253349;
    case "F"
        subjectN = 781520253349;
end
end
%subjectN = 820253349;
%subjectN = 81520253349;
PD.eventType = eventType;
PD.subjectN = subjectN;

%Mapping mask to gmSurface
% mask_95_allnodes_allch = zeros(size(vol2gm,2),1);
% mask_95_allnodes_allch(cortex_nodes(masks.mask_95_cortexnodes_allch)) = 1;
% mask_95_allnodes_allch_gmsurface = vol2gm*mask_95_allnodes_allch;
% 
% mask_95_allnodes_allch_gmsurface(find(mask_95_allnodes_allch_gmsurface~=0))=1;

%mask_j_thresh_cortexnodes_allch
%masks.mask_nodes_bottom_j_thresh_cortexnodes_allch

mask_j_thresh_allnodes_allch_gmsurface = zeros(size(vol2gm,1),1);
mask_j_thresh_allnodes_allch_gmsurface(masks.mask_j_thresh_cortexnodes_allch,1)=1;

%top 95% JAC All Ch = 1, bottom 5% = 0;
% figure; plotmesh_iso2([gmSurfaceMesh.node mask_95_allnodes_allch_gmsurface],gmSurfaceMesh.face)
% cb = colorbar('horiz');

figure; plotmesh_iso2([gmSurfaceMesh.node mask_j_thresh_allnodes_allch_gmsurface],gmSurfaceMesh.face)
cb = colorbar('horiz');

% nnz(mask_95_allnodes_allch_gmsurface)
% 
% cortex_nodes(masks.mask_95_cortexnodes_allch);
% cortex_nodes;

% run this
rho = zeros(size(HbT_all_gmsurface,1),1);
pval = zeros(size(HbT_all_gmsurface,1),1);

t_val = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
h_val_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
ci_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));
pval_t = zeros(size(HbT_all_gmsurface,1),size(HbT_all_gmsurface,3));


%for i=1:size(cortex_nodes,1)
for i=1:size(HbT_all_gmsurface,1)

    %HbT_all_S_hypo_max = (max(HbT_all_S_hypo(i,121:end,:)));
    
    %HbT_all_max = (max(HbT_all(i,121:end,:))); %get max HbT from T=2:10mins
    
    %used pre 31 05 24
    %HbT_all_max = (max(HbT_all_m_hypo_gmsurface(i,121:end,:))); %get max HbT from T=2:10mins
    %used 31 05 24
    %HbT_all_max = (max(abs(HbT_all_gmsurface(i,121:end,:)))); %get max HbT from T=2:10mins


    [HbT_all_max HbT_all_max_idx] = (max(abs(HbT_all_gmsurface(i,121:end,:)))); %get max HbT from T=2:10mins
    HbT_all_max_idx = reshape(HbT_all_max_idx(1,1,:), 1, []);
    
    HbT_all_max = HbT_all_gmsurface(i,121+HbT_all_max_idx-1,:); %get max HbT from T=2:10mins
    HbT_all_max = diag(reshape(HbT_all_max,[],size(HbT_all_max,3))) ;%rehsape into sqaure, and get diag

    

    %[a b]= max(.....) %one output will be index, one output will be the
    %value (from input i.e 121:end)
    %find index of abs max
    %but keep information of + or - value
    % updated 26 08 24


    %get b and rho values for each node
    %[b, stats] = robustfit(peak_glucose_all_S_hypo,HbT_all_S_hypo_max(1,:)'   );
    %[rho(i),pval(i)] = corr(HbT_all_S_hypo_max(1,:)',   b(1) + b(2)*peak_glucose_all_S_hypo     );
    
    %peak glucose
    if size(peak_glucose_all,1) > 2 %b fit only works with > 2 events 

        %[b, stats] = robustfit(peak_glucose_all,HbT_all_max(1,:)'   );
        [b, stats] = robustfit(peak_glucose_all,HbT_all_max   );

        %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   b(1) + b(2)*peak_glucose_all     );
        [rho(i),pval(i)] = corr(HbT_all_max,   b(1) + b(2)*peak_glucose_all     );

    else
        %[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   peak_glucose_all     );
        [rho(i),pval(i)] = corr(HbT_all_max,   peak_glucose_all     );

    end

    %length of glucose
    %[b, stats] = robustfit(length_glucose_all_m_S_hypo,HbT_all_m_S_hypo_max(1,:)'   );
    %[rho(i),pval(i)] = corr(HbT_all_m_S_hypo_max(1,:)',   b(1) + b(2)*length_glucose_all_m_S_hypo     );
    for j=1:size(HbT_all_gmsurface,3)
        [t_test(i,j),pval_t(i,j)] = ttest(HbT_all_gmsurface(i,:,j));
        [h_val_t(i,j),pval_t(i,j),ci_t,stats] = ttest(HbT_all_gmsurface(i,:,j));
        t_val(i,j) = stats.tstat;
    end
end 


%[h,p,ci,stats] = ttest(HbT_all_m_hypo_gmsurface(1,:,1));
%t_val = stats.tstat
%
% run this
%pval(logical(mask_95_allnodes_allch_gmsurface))
% set up adj p matrix of ones
%mask_j_thresh_allnodes_allch_gmsurface

adj_p = ones(1,size(mask_j_thresh_allnodes_allch_gmsurface,1));
adj_p_t = ones(size(HbT_all_gmsurface,3),size(mask_j_thresh_allnodes_allch_gmsurface,1));

%adj_p_mask_t = ones(size(mask_95_allnodes_allch_gmsurface,1),size(HbT_all_m_hypo_gmsurface,3));

% FDR correction
[h, crit_p, adj_ci_cvrg, adj_p_mask]=fdr_bh(pval(logical(mask_j_thresh_allnodes_allch_gmsurface)));


%fix this
for i=1:size(HbT_all_gmsurface,3)
    [h_t, crit_p_t, adj_ci_cvrg_t, adj_p_mask_t(:,i)]=fdr_bh(pval_t(logical(mask_j_thresh_allnodes_allch_gmsurface)));
end

%set nodes in top 95% mask to be corrected p vals.
adj_p(logical(mask_j_thresh_allnodes_allch_gmsurface)) = adj_p_mask;

for i=1:size(adj_p_mask_t,2)
    adj_p_t(i,logical(mask_j_thresh_allnodes_allch_gmsurface)) = adj_p_mask_t(:,i);
end

%max(rho)
%set R P to be zero for nodes outside of mask.
%rho(masks.mask_nodes_bottom_5_cortexnodes_allch)=0;
%pval(masks.mask_nodes_bottom_5_cortexnodes_allch)=-1;

rho(~logical(mask_j_thresh_allnodes_allch_gmsurface))=0;
pval(~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;
adj_p(~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;

t_val(~logical(mask_j_thresh_allnodes_allch_gmsurface),:)=0;
pval_t(~logical(mask_j_thresh_allnodes_allch_gmsurface),:)=-1;
adj_p_t(:,~logical(mask_j_thresh_allnodes_allch_gmsurface))=-1;

if size(peak_glucose_all,1) > 1 %max rho > 1 events 
    max_rho = max([max(rho) abs(min(rho))]);
else
    max_rho = 1;
end
max_pval = max([max(pval) abs(min(pval))]);
max_adj_p = max([max(adj_p) abs(min(adj_p))]);

max_t_val = max([max(t_val) abs(min(t_val))]);
max_pval_t = max([max(pval_t) abs(min(pval_t))]);
max_adj_p_t = max([max(adj_p_t) abs(min(adj_p_t))]);
% run this plot mesh (iso2mesh)
%rho=rho';
%pval=pval';
figure()
subplot(1,3,1)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(0,90); %TOP VIEW % used for infant week 30
clim([-max_rho max_rho]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
subplot(1,3,2)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(0,90);
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(1,3,3)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(0,90);
clim([-max_adj_p max_adj_p]);
cb = colorbar('horiz');
ylabel(cb,"FDR cor. P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")

figure()
subplot(1,3,1)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(0,90); %TOP VIEW % used for infant week 30
%clim([-max_rho max_rho]);
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
subplot(1,3,2)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(0,90);
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(1,3,3)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(0,90);
clim([-max_adj_p max_adj_p]);
cb = colorbar('horiz');
ylabel(cb,"FDR cor. P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_val.png")
end

figure()
subplot(1,3,1)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(0,90); %TOP VIEW % used for infant week 30
%clim([-max_rho max_rho]);
clim([-0.05 0.05]);
cb = colorbar('vertical');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
ax = gca;
ax.FontSize = 16;
subplot(1,3,2)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(0,90);
clim([-0.05 0.05]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
ax = gca;
ax.FontSize = 16;
subplot(1,3,3)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(0,90);
clim([-max_adj_p max_adj_p]);
cb = colorbar('horiz');
ylabel(cb,"FDR cor. P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_val_p0_05.png")
end
ax = gca;
ax.FontSize = 16;

%average across p<0.05 nodes to get scatter plot

p_val_0_05_R = intersect(find(pval<=0.05),find(pval>=0));
max_HbT_p_val_0_05_R_all = zeros(size(p_val_0_05_R,1),size(good_PD_events,2));
%max_HbT_p_val_0_05_R;
%peak_glucose_all
for i=1:size(p_val_0_05_R,1)
    %Abs HbT
    %max_HbT_p_val_0_05_R(i,:) = max(abs(HbT_all_gmsurface(p_val_0_05_R(i),121:end,:)));

    [max_HbT_p_val_0_05_R max_HbT_p_val_0_05_R_idx] = (max(abs(HbT_all_gmsurface(p_val_0_05_R(i),121:end,:)))); %get max HbT from T=2:10mins
    max_HbT_p_val_0_05_R_idx = reshape(max_HbT_p_val_0_05_R_idx(1,1,:), 1, []);
    
    max_HbT_p_val_0_05_R = HbT_all_gmsurface(p_val_0_05_R(i),121+max_HbT_p_val_0_05_R_idx-1,:); %get max HbT from T=2:10mins
    max_HbT_p_val_0_05_R = diag(reshape(max_HbT_p_val_0_05_R,[],size(max_HbT_p_val_0_05_R,3)))';%rehsape into sqaure, and get diag
    max_HbT_p_val_0_05_R_all(i,:) = max_HbT_p_val_0_05_R;

end

c_max_HbT_p_val_0_05_R = polyfit(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:)),1);
fit_y_c_p_glucose_max_r_HbT = c_max_HbT_p_val_0_05_R(1)*peak_glucose_all + c_max_HbT_p_val_0_05_R(2);


figure()
plot(peak_glucose_all,max_HbT_p_val_0_05_R_all(:,:),'bo')
hold on
plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'r--')
xlabel("Peak Glucose Value / mg/dL")
ylabel("Max HbT / \mu M")
title("Combined plot R of P<0.05- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" N CORTEX NODES "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_p_0_05nodes_scatter.png")
end

%c_max_HbT_p_val_0_05_R = polyfit(peak_glucose_all,mean(max_HbT_p_val_0_05_R(:,:)),1);
%fit_y_c_p_glucose_max_r_HbT = c_max_HbT_p_val_0_05_R(1)*peak_glucose_all + c_max_HbT_p_val_0_05_R(2);


figure()
plot(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:)),'bo')
hold on
plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'r--')
xlabel("Peak Glucose Value / mg/dL")
ylabel("Max HbT / \mu M")
title("Mean Combined plot R of P<0.05- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" N CORTEX NODES "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_p_0_05nodes_scatter_mean.png")
end
% testing coloured plot/scatter
%figure()
%plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'r--','LineWidth',1)
%hold on
%figure()
%hold on
%for i=1:size(good_PD_events_subN,2)
%get colour
%good_PD_events_subN(i);  %'MarkerFaceColor',[0.5,0.5,0.5]
%plot(peak_glucose_all(i),mean(max_HbT_p_val_0_05_R_all(:,i)),'x','LineWidth',7,'MarkerFaceColor',[good_PD_events_subN(i)/100,good_PD_events_subN(i)/100, good_PD_events_subN(i)/100]);
%s1 = plot(peak_glucose_all(i),mean(max_HbT_p_val_0_05_R_all(:,i)),'x','LineWidth',7);
%s1.MarkerFaceColor = [good_PD_events_subN(i)/100 good_PD_events_subN(i)/100 good_PD_events_subN(i)/100];
%end
%xlabel("Peak Glucose Value / mg/dL")
%ylabel("Max HbT / \mu M")
%title("Mean Combined plot R of P<0.05- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" N CORTEX NODES "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
%testing colour scatter options
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%legend('fit',split(num2str(good_PD_events_subN)))
%legend('x','c')
%legend(num2str(good_PD_events_subN))
%newStr = split(num2str(good_PD_events_subN))
%scatter plot colours

%figure()
%plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'r--','LineWidth',1)
%hold on

%1:1:size(unique(good_PD_events_subN'),1)
%good_PD_events_subN
%rng(2,"twister")
%seed = rng;

%RGB_col = [102/255 197/255 204/255;
%246/255 207/255 113/255;
%248/255 156/255 116/255;
%220/255 176/255 242/255;
%135/255 197/255 95/255;
%179/255 179/255 179/255];

RGB_col = [241/255 71/255 37/255;
241/255 245/255 33/255;
17/255 251/255 56/255;
255/255 144/255 37/255;
236/255 137/255 132/255;
202/255 125/255 187/255;
%136/255 8/255 165/255;
63/255 73/255 254/255];

%col_scat = [1.6*good_PD_events_subN'/100, 1.5*good_PD_events_subN'/100, 0.8*good_PD_events_subN'/100];

col_scat_idx = zeros(size(good_PD_events_subN,2),1);
col_scat = zeros(size(good_PD_events_subN,2),1);
for i=1:size(good_PD_events_subN,2)
    col_scat_idx(i,1) = find(good_PD_events_subN(i) == unique(good_PD_events_subN));
end

[b, stats] = robustfit(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:))   );

%[rho(i),pval(i)] = corr(HbT_all_max(1,:)',   b(1) + b(2)*peak_glucose_all     );
[rho_HbT_pval_0_05_R_all,pvalHbT_pval_0_05_R_all] = corr(mean(max_HbT_p_val_0_05_R_all(:,:))',   b(1) + b(2)*peak_glucose_all     );

%plot(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:)),'bo')

%% 
%timeWA psoter fNIRS 04 09 24
%time w A
col_scat_idx(1:3,1) = 2;
col_scat_idx(4:8,1) = 4;
col_scat_idx(9:11,1) = 5;
col_scat_idx(12:22,1) = 6;
col_scat_idx(23:31,1) = 7;
%time W D
col_scat_idx(1:3,1) = 2;
col_scat_idx(4:5,1) = 3;
col_scat_idx(6:13,1) = 4;
col_scat_idx(14,1) = 5;
col_scat_idx(15:27,1) = 6;
col_scat_idx(28:38,1) = 7;
%time W F
%7
col_scat_idx(1,1) = 1;
col_scat_idx(2:6,1) = 2;
col_scat_idx(7:8,1) = 3;
col_scat_idx(9:13,1) = 4;
col_scat_idx(14:15,1) = 5;
col_scat_idx(16:25,1) = 6;
col_scat_idx(26:37,1) = 7;

col_scat = RGB_col(col_scat_idx(:,1),:);

%savefig=1;

figure()
plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'k-','LineWidth',8)
hold on
scatter(peak_glucose_all,mean(max_HbT_p_val_0_05_R_all(:,:)),130,col_scat,'filled')
xlim([40 72])
ylim([-13 13])
xlabel("Peak Glucose Value / mg/dL")
ylabel("Max HbT / \mu M")
title("Mean Combined plot R of P<0.05- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" N CORTEX NODES "+num2str(size(p_val_0_05_R,1))+" TW "+time_window+" p1 "+num2str(c_max_HbT_p_val_0_05_R(1))+" p0 "+num2str(c_max_HbT_p_val_0_05_R(2))+" ")
ax = gca;
ax.FontSize = 16;

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_p_0_05nodes_scatter_mean_colour.png")
end

savefig=0;

figure()
subplot(2,3,1)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(90,0); %TOP VIEW % used for infant week 30
%clim([-max_rho max_rho]);
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
subplot(2,3,4)
plotmesh_iso2([gmSurfaceMesh.node rho],gmSurfaceMesh.face)
view(-90,0);
%clim([-max_rho max_rho]);
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"R value")
title("R val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+"")
subplot(2,3,2)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(90,0); %TOP VIEW % used for infant week 30
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(2,3,5)
plotmesh_iso2([gmSurfaceMesh.node pval],gmSurfaceMesh.face)
view(-90,0);
clim([-0.1 0.1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(2,3,3)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(90,0); %TOP VIEW % used for infant week 30
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
subplot(2,3,6)
plotmesh_iso2([gmSurfaceMesh.node adj_p'],gmSurfaceMesh.face)
view(-90,0);
clim([-1 1]);
cb = colorbar('horiz');
ylabel(cb,"P value")
title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_R_val_side.png")
end

% run this
find(rho==max(rho));
pval(find(rho==max(rho)));

%max_rho_HbT = max(HbT_all(find(rho==max(rho)),121:end,:));

%max_rho_HbT = max(HbT_all_m_hypo_gmsurface(find(rho==max(rho)),121:end,:));

%max rho HbT using max abs HbT all pre 26 08 24
%max_rho_HbT = max(abs(HbT_all_gmsurface(find(rho==max(rho)),121:end,:)));

%max rho HbT using + and - polarity
%max_rho_HbT = HbT_all_gmsurface(find(rho==max(rho)),121:end,:);

[max_rho_HbT max_rho_HbT_idx] = (max(abs(HbT_all_gmsurface(find(rho==max(rho)),121:end,:)))); %get max HbT from T=2:10mins
max_rho_HbT_idx = reshape(max_rho_HbT_idx(1,1,:), 1, []);
    
max_rho_HbT = HbT_all_gmsurface(find(rho==max(rho)),121+max_rho_HbT_idx-1,:); %get max HbT from T=2:10mins
max_rho_HbT = diag(reshape(max_rho_HbT,[],size(max_rho_HbT,3)));%rehsape into sqaure, and get diag



%max_rho_HbT = max(HbT_all(33044,121:end,:));
%max_rho_HbT = max_rho_HbT(1,:)'; %pre 26 08 24 
peak_glucose_all;

switch subjectN
     case 7
        figure()
        plot(peak_glucose_all,max_rho_HbT,'x')
        legend('m hypo')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+"")
 
     case 8
        figure()
        plot(peak_glucose_all,max_rho_HbT,'x')
        legend('m hypo')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+"")
 
    case 20
        figure()
        plot(peak_glucose_all,max_rho_HbT,'x')
        legend('m hypo')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+"")
    case 33
        figure()
        plot(peak_glucose_all(1:end-2),max_rho_HbT(1:end-2),'x')
        hold on
        plot(peak_glucose_all(end-1:end),max_rho_HbT(end-1:end),'rx')
        legend('m hypo','S hypo')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+"")
end


c_p_glucose_max_r_HbT = polyfit(peak_glucose_all,max_rho_HbT,1);
fit_y_c_p_glucose_max_r_HbT = c_p_glucose_max_r_HbT(1)*peak_glucose_all + c_p_glucose_max_r_HbT(2);

figure()
plot(peak_glucose_all,max_rho_HbT,'x')
hold on
        plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'o--')
        legend('m hypo','p1 fit','Location','northwest')
        xlabel("Peak Glucose Value / mg/dL")
        ylabel("Max HbT / \mu M")
        title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+" p1 "+num2str(c_p_glucose_max_r_HbT(1))+" p0 "+num2str(c_p_glucose_max_r_HbT(2))+" ")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_maxR_node.png")
end

%savefig=1;
%colour
figure()
scatter(peak_glucose_all,max_rho_HbT,65,col_scat,'filled')
hold on
plot(peak_glucose_all,fit_y_c_p_glucose_max_r_HbT,'k-','LineWidth',1)
xlim([40 72])
ylim([-13 13])
%legend('m hypo','p1 fit','Location','northwest')
xlabel("Peak Glucose Value / mg/dL")
ylabel("Max HbT / \mu M")
title("Max R "+num2str(max(rho))+" P("+num2str(pval(find(rho==max(rho))))+")- Glc vs Max HbT. Sub. "+num2str(subjectN)+" "+eventType+" CORTEX NODE "+num2str(find(rho==max(rho)))+" TW "+time_window+" p1 "+num2str(c_p_glucose_max_r_HbT(1))+" p0 "+num2str(c_p_glucose_max_r_HbT(2))+" ")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
    saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_maxR_node_colour.png")
end

savefig=0;
%%
%t test, t val.
% switch subjectN
%     case 7
%         eventN_t_test = good_PD_events;
%     case 8
%         eventN_t_test = good_PD_events;
%     case 20
%         eventN_t_test = good_PD_events;
%     case 25
%         eventN_t_test = good_PD_events;
%     case 33
%         eventN_t_test = good_PD_events;
% end


% t val stuff below
eventN_t_test = good_PD_events;

for i=1:size(t_val,2)
    figure()
    subplot(1,3,1)
    plotmesh_iso2([gmSurfaceMesh.node t_val(:,i)],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    clim([-max_t_val max_t_val]);
    cb = colorbar('horiz');
    ylabel(cb,"T value")
    title("T val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" "+num2str(eventN_t_test(i))+" TW "+time_window+"")
    subplot(1,3,2)
    plotmesh_iso2([gmSurfaceMesh.node pval_t(:,i)],gmSurfaceMesh.face)
    view(0,90);
    clim([-0.1 0.1]);
    cb = colorbar('horiz');
    ylabel(cb,"P value")
    title("P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" "+num2str(eventN_t_test(i))+" TW "+time_window+"")
    subplot(1,3,3)
    plotmesh_iso2([gmSurfaceMesh.node adj_p_t(i,:)'],gmSurfaceMesh.face)
    view(0,90);
    %clim([-max_adj_p_t max_adj_p_t]);
    clim([-0.1 0.1]);
    cb = colorbar('horiz');
    ylabel(cb,"FDR cor. P value")
    title("FDR adj P val. Peak Glc vs Max HbT, PD"+num2str(subjectN)+" "+eventType+" "+num2str(eventN_t_test(i))+" TW "+time_window+"")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val"+num2str(eventN_t_test(i))+".png")
    end
end


% t test DICE

%t_val
%mask_95_allnodes_allch_gmsurface

t_val_pos_mask = zeros(size(t_val,1),size(t_val,2));
t_val_neg_mask = zeros(size(t_val,1),size(t_val,2));
t_val_zero_mask = zeros(size(t_val,1),size(t_val,2));

t_val_pos_neg_mask_net = zeros(size(t_val,1),1);

for i=1:size(t_val,2)
    t_val_pos = find(t_val(:,i)>0);
    t_val_neg = find(t_val(:,i)<0);
    t_val_zero = find(t_val(:,i)==0);


    t_val_pos_mask(t_val_pos,i) = 1;
    t_val_neg_mask(t_val_neg,i) = 1;
    t_val_zero_mask(t_val_zero,i) = 1;
end

figure()
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_mask(:,1)],gmSurfaceMesh.face)
cb = colorbar('horiz');

%for each node, sum if mask is pos or neg, 
%pos value t_val_pos_neg_mask_net means more events HbT increase
%neg value t_val_pos_neg_mask_net means more events HbT decrease

%only consider nodes with p<0.05 for the tvalue


%find(pval_t(:,1)>0.05);

%set nodes which have P val (t) > 0.05 to zero
for i=1:size(t_val,2)
    t_val_pos_mask(find(pval_t(:,i)>0.05),i) = 0;
    t_val_neg_mask(find(pval_t(:,i)>0.05),i) = 0;
end

t_val_pos_mask_net = sum(t_val_pos_mask,2);
t_val_neg_mask_net = sum(t_val_neg_mask,2);

t_val_pos_neg_mask_net = t_val_pos_mask_net - t_val_neg_mask_net;

%t_val_pos_neg_mask_net_norm = t_val_pos_neg_mask_net/max(abs(t_val_pos_neg_mask_net));

t_val_pos_neg_mask_net_norm = t_val_pos_neg_mask_net/size(t_val,2);

%size(t_val,2)

figure()
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(0,90); %TOP VIEW % used for infant week 30
title("T val NET (p<0.05) (N Pos - N Neg), PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_net.png")
end

figure()
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(0,90); %TOP VIEW % used for infant week 30
title("T val NET Norm (p<0.05) (N Pos - N Neg)/N Events, PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_net_norm.png")
end

figure()
subplot(2,2,1)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(-90,0); %TOP VIEW % used for infant week 30
subplot(2,2,2)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(90,0); %TOP VIEW % used for infant week 30
subplot(2,2,3)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(0,0); %TOP VIEW % used for infant week 30
subplot(2,2,4)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net],gmSurfaceMesh.face)
clim([-size(t_val,2) size(t_val,2)])
cb = colorbar('horiz');
view(180,0); %TOP VIEW % used for infant week 30
sgtitle("T val NET (p<0.05) (N Pos - N Neg), PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_net_side.png")
end

figure()
subplot(2,2,1)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(-90,0); %TOP VIEW % used for infant week 30
subplot(2,2,2)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(90,0); %TOP VIEW % used for infant week 30
subplot(2,2,3)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(0,0); %TOP VIEW % used for infant week 30
subplot(2,2,4)
plotmesh_iso2([gmSurfaceMesh.node t_val_pos_neg_mask_net_norm],gmSurfaceMesh.face)
clim([-1 1])
cb = colorbar('horiz');
view(180,0); %TOP VIEW % used for infant week 30
sgtitle("T val NET Norm (p<0.05) (N Pos - N Neg)/N Events, PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_net_norm_side.png")
end


for i=1:size(t_val,2)
    for j=1:size(t_val,2)
        t_val_dice_pos(i,j) = dice(t_val_pos_mask(:,i),t_val_pos_mask(:,j));
        t_val_dice_neg(i,j) = dice(t_val_neg_mask(:,i),t_val_neg_mask(:,j));
    end
end

figure()
subplot(1,2,1)
imagesc(t_val_dice_pos); 
xlabel("Event N")
xlabel("Event N")
title("T val DICE (Pos), PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
clim([0 1]);
colorbar;
subplot(1,2,2)
imagesc(t_val_dice_neg); 
xlabel("Event N")
xlabel("Event N")
title("T val DICE (Neg), PD"+num2str(subjectN)+" "+eventType_all+" TW "+time_window+"")
clim([0 1]);
colorbar;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_matrix.png")
end

% Average HbT across all events

    t_val_AVG = mean(t_val')';
    t_val_pc_pos = size(find(t_val_AVG>0),1)/nnz(mask_j_thresh_allnodes_allch_gmsurface)*100;
    t_val_pc_neg = size(find(t_val_AVG<0),1)/nnz(mask_j_thresh_allnodes_allch_gmsurface)*100;

    %max(abs(t_val_AVG))

    figure()
    plotmesh_iso2([gmSurfaceMesh.node t_val_AVG],gmSurfaceMesh.face)
    view(0,90); %TOP VIEW % used for infant week 30
    %clim([-max_t_val max_t_val]);
    clim([-max(abs(t_val_AVG)) max(abs(t_val_AVG))]);
    cb = colorbar('horiz');
    ylabel(cb,"T value")
    title("T val. AVG HbT, PD"+num2str(subjectN)+" "+eventType_all+" "+num2str(eventN_t_test(i))+" TW "+time_window+" "+num2str(t_val_pc_pos)+"% + "+num2str(t_val_pc_neg)+" % -")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_T_val_avg.png")
    end


    %HbT_all_m_hypo_gmsurface_AVG = mean(HbT_all_m_hypo_gmsurface);

%create overlap between active R value nodes and Net T value nodes

%savefig = 1;
%pick three colours
%red - R val active only 
%blue - Net T active only
%green - Both R and Net V active
%active could be defined by >0.5 or <-0.5 in each metric

%r 
rho;
r_thresh_active = 0.5;
r_thresh_active_pval = 0.05;
t_thresh_active = 0.25;

%rho_active_mask = find(abs(rho)>=r_thresh_active);
rho_active_mask = find(abs(pval)<=r_thresh_active_pval);

%t net
t_val_pos_neg_mask_net_norm;

t_val_active_mask = find(abs(t_val_pos_neg_mask_net_norm)>=t_thresh_active);

r_t_val_active_mask = zeros(size(HbT_all_gmsurface,1),1);


r_t_val_active_mask(rho_active_mask,1) = r_t_val_active_mask(rho_active_mask,1)+1;
r_t_val_active_mask(t_val_active_mask,1) = r_t_val_active_mask(t_val_active_mask,1)+2;

figure()
plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
view(0,90)
clim([0 3])
cb = colorbar('horiz');
ylabel(cb,"Active R(1) T(2) RT(3)")
    title("Active R T val HbT, PD"+num2str(subjectN)+" "+eventType_all+" N events = "+num2str(size(eventN_t_test,2))+" TW "+time_window+" R thresh pval< "+num2str(r_thresh_active_pval)+" T thresh "+num2str(t_thresh_active)+"")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_active_R_T_val.png")
    end

% figure()
% plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
% view(0,90)
% clim([0 3])
% cb = colorbar('horiz');
% ylabel(cb,"Active R(1) T(2) RT(3)")
%     title("Active R T val HbT, PD"+num2str(subjectN)+" "+eventType_all+" N events = "+num2str(size(eventN_t_test,2))+" TW "+time_window+" R thresh pval "+num2str(r_thresh_active_pval)+" T thresh "+num2str(t_thresh_active)+"")
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     if savefig == 1
%         saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_active_R_T_val.png")
%     end

figure()
subplot(1,3,1)
plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
view(0,90)
clim([0 3])
cb = colorbar('horiz');
ylabel(cb,"Active R(1) T(2) RT(3)")
    title("Active R T val HbT, PD"+num2str(subjectN)+" "+eventType_all+" N events = "+num2str(size(eventN_t_test,2))+" TW "+time_window+" R thresh pval "+num2str(r_thresh_active_pval)+" T thresh "+num2str(t_thresh_active)+"")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(1,3,2)
plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
view(-90,0)
clim([0 3])
cb = colorbar('horiz');
ylabel(cb,"Active R(1) T(2) RT(3)")
title("Left Hemishpere")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(1,3,3)
plotmesh_iso2([gmSurfaceMesh.node r_t_val_active_mask],gmSurfaceMesh.face)
view(90,0)
clim([0 3])
cb = colorbar('horiz');
ylabel(cb,"Active R(1) T(2) RT(3)")
title("Right Hemishpere")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_active_R_T_val_CLRview.png")
    end


%% R value in 9 ROI
% Get 6 regions of cortex
landmark_EEG1020 = [
	-23.57	14.32	24.47;
	-1.09	15.9	37.04;
 	21.63	12.93	26.08;
 	-32.89	-16.63	35.48;
 	-2.01	-17.37	51.07;
 	29.67	-16.27	36.75;
 	-24.76	-50.57	27.24;
 	-1.33	-52.76	39.53;
 	22.94	-49.58	29.2 ];

ROI_data = ones(size(gmSurfaceMesh.node,1),1)*4 ;%  abs(randn(size(gmSurfaceMesh.node,1),1))+100; %RANDOMISE AND SET SO IT'S AT MAX CBAR
%GETTTING CORTEX ZONES 
%F3  FZ  F4
%C3  CZ  C4
%P3  PZ  P4
radius_LM = 20;
landmark_EEG1020_cortexF3 = find(pdist2(landmark_EEG1020(1,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexF3) = -5;
landmark_EEG1020_cortexFZ = find(pdist2(landmark_EEG1020(2,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexFZ) = -4;
landmark_EEG1020_cortexF4 = find(pdist2(landmark_EEG1020(3,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexF4) = -3;
landmark_EEG1020_cortexC3 = find(pdist2(landmark_EEG1020(4,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexC3) = -2;
landmark_EEG1020_cortexCZ = find(pdist2(landmark_EEG1020(5,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexCZ) = -1;
landmark_EEG1020_cortexC4 = find(pdist2(landmark_EEG1020(6,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexC4) = -0;
landmark_EEG1020_cortexP3 = find(pdist2(landmark_EEG1020(7,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexP3) = 1;
landmark_EEG1020_cortexPZ = find(pdist2(landmark_EEG1020(8,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexPZ) = 2;
landmark_EEG1020_cortexP4 = find(pdist2(landmark_EEG1020(9,:),gmSurfaceMesh.node)<radius_LM);
ROI_data(landmark_EEG1020_cortexP4) = 3;

EEG_lm_9ROI = [];

figure()
plotmesh_iso2([gmSurfaceMesh.node ROI_data],gmSurfaceMesh.face)
clim([-5 4])
cb = colorbar('horiz');
ylabel(cb,'EEG LM Zone','FontSize',10,'Rotation',0)
title("9 EEG zones on cortex F3z4 C3z4 P3z4")

HbT_all_max_ROI = zeros(9,size(peak_glucose_all,1));

%HbT_all_max_ROI(1,:) = max(max(abs(HbT_all_m_hypo_gmsurface(landmark_EEG1020_cortexF3,121:end,:))));
%abs(mean(mean(HbT_all_m_hypo_gmsurface(landmark_EEG1020_cortexF3,121:end,:))));
%HbT_all_m_hypo_gmsurface(landmark_EEG1020_cortexF3,121:end,:)

%HbT_all_m_hypo_gmsurface

%max_HbT_ROI_9_all = zeros(size(p_val_0_05_R,1),size(good_PD_events,2));
%need to change to max, so that it considers + and - max
%[max_HbT_ROI_9 max_HbT_ROI_9_idx] = max(HbT_all_gmsurface(landmark_EEG1020_cortexF3,121:end,:)); %get max HbT from T=2:10mins
%HbT_all_max_ROI = zeros(9,size(peak_glucose_all,1));
[HbT_all_max_ROI(1,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexF3);
[HbT_all_max_ROI(2,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexFZ);
[HbT_all_max_ROI(3,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexF4);

[HbT_all_max_ROI(4,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexC3);
[HbT_all_max_ROI(5,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexCZ);
[HbT_all_max_ROI(6,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexC4);

[HbT_all_max_ROI(7,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexP3);
[HbT_all_max_ROI(8,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexPZ);
[HbT_all_max_ROI(9,:)] = get_max_min_HbT_ROI_9(HbT_all_gmsurface,landmark_EEG1020_cortexP4);

%OLD JUST USING MAX ABS
% HbT_all_max_ROI(1,:) = max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexF3,121:end,:))));
% HbT_all_max_ROI(2,:) =max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexFZ,121:end,:))));
% HbT_all_max_ROI(3,:) =max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexF4,121:end,:))));
% 
% HbT_all_max_ROI(4,:) = max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexC3,121:end,:))));
% HbT_all_max_ROI(5,:) =max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexCZ,121:end,:))));
% HbT_all_max_ROI(6,:) =max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexC4,121:end,:))));
% 
% HbT_all_max_ROI(7,:) = max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexP3,121:end,:))));
% HbT_all_max_ROI(8,:) =max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexPZ,121:end,:))));
% HbT_all_max_ROI(9,:) =max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexP4,121:end,:))));

%HbT_all_max_ROI(1,:) =  max(abs(mean(HbT_all_m_hypo_gmsurface(landmark_EEG1020_cortexF3,121:end,:))))

%avg
HbT_all_avg_ROI(1,:) = max(max(abs(HbT_all_gmsurface(landmark_EEG1020_cortexF3,121:end,:))));
%mean(HbT_all_gmsurface(landmark_EEG1020_cortexF3,121:end,:));

%get R value in each ROI
rho = zeros(size(HbT_all_max_ROI,1),1);
pval = zeros(size(HbT_all_max_ROI,1),1);

c_HbT_all_max_ROI = zeros(size(HbT_all_max_ROI,1),2);

fit_y_c_HbT_all_max_ROI = zeros(size(HbT_all_max_ROI,1),size(peak_glucose_all,1));

for i=1:9
    %HbT_all_S_hypo_max = (max(HbT_all_S_hypo(i,121:end,:)));
    %HbT_all_max = (max(HbT_all(i,121:end,:))); %get max HbT from T=2:10mins
    %used pre 31 05 24
    %HbT_all_max = (max(HbT_all_m_hypo_gmsurface(i,121:end,:))); %get max HbT from T=2:10mins
    %used 31 05 24
    %HbT_all_max = (max(abs(HbT_all_m_hypo_gmsurface(i,121:end,:)))); %get max HbT from T=2:10mins
    %get b and rho values for each node
    %[b, stats] = robustfit(peak_glucose_all_S_hypo,HbT_all_S_hypo_max(1,:)'   );
    %[rho(i),pval(i)] = corr(HbT_all_S_hypo_max(1,:)',   b(1) + b(2)*peak_glucose_all_S_hypo     );
 
    %peak glucose
    [b, stats] = robustfit(peak_glucose_all,HbT_all_max_ROI(i,:)'   );
    [rho(i),pval(i)] = corr(HbT_all_max_ROI(i,:)',   b(1) + b(2)*peak_glucose_all     );

    %length of glucose
    %[b, stats] = robustfit(length_glucose_all_m_S_hypo,HbT_all_m_S_hypo_max(1,:)'   );
    %[rho(i),pval(i)] = corr(HbT_all_m_S_hypo_max(1,:)',   b(1) + b(2)*length_glucose_all_m_S_hypo     );
    %for j=1:size(HbT_all_m_hypo_gmsurface,3)
    %    [t_test(i,j),pval_t(i,j)] = ttest(HbT_all_m_hypo_gmsurface(i,:,j));
    %    [h_val_t(i,j),pval_t(i,j),ci_t,stats] = ttest(HbT_all_m_hypo_gmsurface(i,:,j));
    %    t_val(i,j) = stats.tstat;
    %end
    c_HbT_all_max_ROI(i,:) = polyfit(peak_glucose_all,HbT_all_max_ROI(i,:)',1);
    fit_y_c_HbT_all_max_ROI(i,:) = c_HbT_all_max_ROI(i,1)*peak_glucose_all + c_HbT_all_max_ROI(i,2);

end 

%get correlation line of best fit
%c_HbT_all_max_ROI(1,:) = polyfit(peak_glucose_all,HbT_all_max_ROI(i,:)',1);
%fit_y_c_HbT_all_max_ROI = c_HbT_all_max_ROI(:,1)*peak_glucose_all + c_HbT_all_max_ROI(:,2);

%fit_y_c_HbT_all_max_ROI(:,i) = c_HbT_all_max_ROI(i,1)*peak_glucose_all + c_HbT_all_max_ROI(i,2);

%PLOT THIS
%%
figure()
sgtitle("EEG ROI HbT v BCG PD"+num2str(subjectN)+" "+eventType_all+" N events = "+num2str(size(HbT_all_gmsurface,3))+" TW "+time_window+"")
subplot(3,3,1)
plot(peak_glucose_all,HbT_all_max_ROI(1,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(1,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("F3 R="+num2str(rho(1))+" P="+num2str(pval(1))+"");
subplot(3,3,2)
plot(peak_glucose_all,HbT_all_max_ROI(2,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(2,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("FZ R="+num2str(rho(2))+" P="+num2str(pval(2))+"");
subplot(3,3,3)
plot(peak_glucose_all,HbT_all_max_ROI(3,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(3,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("F4 R="+num2str(rho(3))+" P="+num2str(pval(3))+"");

subplot(3,3,4)
plot(peak_glucose_all,HbT_all_max_ROI(4,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(4,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("C3 R="+num2str(rho(4))+" P="+num2str(pval(4))+"");
subplot(3,3,5)
plot(peak_glucose_all,HbT_all_max_ROI(5,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(5,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("CZ R="+num2str(rho(5))+" P="+num2str(pval(5))+"");
subplot(3,3,6)
plot(peak_glucose_all,HbT_all_max_ROI(6,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(6,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("C4 R="+num2str(rho(6))+" P="+num2str(pval(6))+"");

subplot(3,3,7)
plot(peak_glucose_all,HbT_all_max_ROI(7,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(7,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("P3 R="+num2str(rho(7))+" P="+num2str(pval(7))+"");
subplot(3,3,8)
plot(peak_glucose_all,HbT_all_max_ROI(8,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(8,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("PZ R="+num2str(rho(8))+" P="+num2str(pval(8))+"");
subplot(3,3,9)
plot(peak_glucose_all,HbT_all_max_ROI(9,:)','bo')
hold on
plot(peak_glucose_all,fit_y_c_HbT_all_max_ROI(9,:)','r--')
ylabel("Max \Delta HbT / \muM");
xlabel("Min. BCG / mg/DL");
title("P4 R="+num2str(rho(9))+" P="+num2str(pval(9))+"");

 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD.subjectN)+"_"+PD.eventType+"_"+num2str(PD.eventN)+"_TimeW_"+PD_data.time_window+"_ROI_9_EEG_R_P_val.png")
    end


%% 

%Intersect between the nodes in the EEG LM SPHERES and Cortex nodes
%but dont need to use as gmsurfacemesh is already just gm surface nodes
% landmark_EEG1020_cortexF3 = intersect(landmark_EEG1020_cortexF3,cortex_nodes);
% landmark_EEG1020_cortexFZ = intersect(landmark_EEG1020_cortexFZ,cortex_nodes);
% landmark_EEG1020_cortexF4 = intersect(landmark_EEG1020_cortexF4,cortex_nodes);
% landmark_EEG1020_cortexC3 = intersect(landmark_EEG1020_cortexC3,cortex_nodes);
% landmark_EEG1020_cortexCZ = intersect(landmark_EEG1020_cortexCZ,cortex_nodes);
% landmark_EEG1020_cortexC4 = intersect(landmark_EEG1020_cortexC4,cortex_nodes);
% landmark_EEG1020_cortexP3 = intersect(landmark_EEG1020_cortexP3,cortex_nodes);
% landmark_EEG1020_cortexPZ = intersect(landmark_EEG1020_cortexPZ,cortex_nodes);
% landmark_EEG1020_cortexP4 = intersect(landmark_EEG1020_cortexP4,cortex_nodes);
% 
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexF3 = landmark_EEG1020_cortexF3;
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexFZ = landmark_EEG1020_cortexFZ;
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexF4 = landmark_EEG1020_cortexF4;
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexC3 = landmark_EEG1020_cortexC3;
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexCZ = landmark_EEG1020_cortexCZ;
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexC4 = landmark_EEG1020_cortexC4;
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexP3 = landmark_EEG1020_cortexP3;
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexPZ = landmark_EEG1020_cortexPZ;
% PD.data.landmarkEEGnodes.landmark_EEG1020_cortexP4 = landmark_EEG1020_cortexP4;



figure()
plotniceimages_1_greyJET(mesh3,mesh_recon);
cb = colorbar('horiz'); 
ylabel(cb,'EEG LM Zone','FontSize',10,'Rotation',0)
max_HbO_cortex = max(mesh3.data(cortex_nodes,1));
min_HbO_cortex = min(mesh3.data(cortex_nodes,1));
caxis([caxis_lim(1) caxis_lim(2)]);
title("9 EEG zones on cortex F3z4 C3z4 P3z4")
  

%%
%FOR XXXXX
% peak_glucose
% HbO_all
% Hb_all
% HbT_all
clear all; clc; close all;

subjectN = 15;
eventN =1; %4 5 8 10 11 12 13 14
%eventType = "m_hypo"; %S_hypo

w=0;
%use this 26 07 23
load('Jacobian_infant30week_850nm.mat')
load('mesh_infant30week_850nm.mat')
cortex_nodes = find(mesh.region==3);

eventType = "m_hypo"; %S_hypo
time_window = "D";



PD_data = load("E:\PadovaPostDoc\BabyGluCo\NIRS_data\PD"+num2str(subjectN)+"\formatted\PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(eventN)+".nirs",'-mat');
PD_data = getfield(PD_data,"PD"+num2str(subjectN)+"_"+eventType+"_X");
%PD_data = PD_data.PD33_m_hypo_X; %change this
PD = load("E:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(subjectN)+".mat");
%PD = load("E:\PadovaPostDoc\BabyGluCo\formatted_PD_GuyPerkins\PD"+num2str(subjectN)+".mat");

PD = getfield(PD,"PD"+num2str(subjectN)+"");
%PD = PD.PD33; %change this
switch eventType
        case "m_hypo"
            PD_time = 0:5:5*PD.events_length.m_hypo(eventN);
        case "S_hypo"
            PD_time = 0:5:5*PD.events_length.S_hypo(eventN);
end
PD_data.time_window = time_window;

PD.subjectN = subjectN;
PD.eventN = eventN;
PD.eventType = eventType;

PD_data.subjectN = subjectN;
PD_data.eventN = eventN;
PD_data.eventType = eventType;

% Crop data. Just want 2mins before, 8 mins after (so total 10)
[PD_data time_window_t] = get_timewindow_10mins(PD_data,PD,PD_time,time_window);

fs=10;
PD_data.fs = fs;

if eventType == "m_hypo"; %S_hypo
    figure()
    plot(PD_time,PD.glucose(PD.events_start_end.m_hypo(eventN,1):PD.events_start_end.m_hypo(eventN,2)),'b-o')
    title("Time W"+time_window+" PD "+num2str(subjectN)+" "+eventType+" "+num2str(eventN)+". nirs")
    xlabel("Time / mins")
    ylabel("BGC / mg/DL")
    xline(15,'g--');xline(PD_time(end)-15,'r--');
    xline(time_window_t(1),'m--');xline(time_window_t(2),'m--');
    yline(72,'k--');
end
if eventType == "S_hypo"; %S_hypo
    figure()
    plot(PD_time,PD.glucose(PD.events_start_end.S_hypo(eventN,1):PD.events_start_end.S_hypo(eventN,2)),'b-o')
    title("Time W"+time_window+" PD "+num2str(subjectN)+" "+eventType+" "+num2str(eventN)+". nirs")
    xlabel("Time / mins")
    ylabel("BGC / mg/DL")
    xline(15,'g--');xline(PD_time(end)-15,'r--');
    xline(time_window_t(1),'m--');xline(time_window_t(2),'m--');
    yline(72,'k--');
end

% Step 2. Signal Qaulity check
dRange = [5E-4 3];
SNRrange = 2;
[goodch_idx, PD_data] = find_good_ch_15mins(PD_data,dRange,SNRrange,1);

% Step 3. Get DoD
meanValue = mean(PD_data.d);
dodConv = -log(abs(PD_data.d)./meanValue);
dodConv_og = -log(abs(PD_data.d)./meanValue);
PD_data.dod = dodConv;

figure()
plot(PD_data.t,PD_data.dod(:,PD_data.goodch_idx))
xlabel('Time / S')
ylabel('\Delta OD')
title("Dod goodch idx - PD "+num2str(subjectN)+" event "+eventType+" N "+num2str(eventN)+" ")


%%%% Step 4. Motion corr. OUTPUT dod1
methodN =1;
PD_data = PD_data_motioncorr_15min(PD_data,methodN,1);

figure()
plot(PD_data.t,PD_data.dod(:,PD_data.goodch_idx(1)))
hold on
plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(1)))
xlabel('Time / S')
ylabel('\Delta OD')
title("Dod1 (Mo corr.) goodch idx - PD "+num2str(subjectN)+" event "+eventType+" N "+num2str(eventN)+" ")

figure()
plot(PD_data.t,PD_data.dod1(:,PD_data.goodch_idx))
xlabel('Time / S')
ylabel('\Delta OD')
title("Dod1 (Mo corr.) goodch idx - PD "+num2str(subjectN)+" event "+eventType+" N "+num2str(eventN)+" ")

%%% Baseline subtraction
PD_data.dod1(:,:) = PD_data.dod1(:,:) - mean(PD_data.dod1(1:2*fs*60,:));
PD_data.dod1_noc(:,:) = PD_data.dod1_noc(:,:) - mean(PD_data.dod1_noc(1:2*fs*60,:));

figure()
plot(PD_data.t,PD_data.dod1(:,PD_data.goodch_idx))
xlabel('Time / S')
ylabel('\Delta OD')
title("Baseline subtraction goodch idx - Subject"+num2str(subjectN)+" event "+eventType+" N "+num2str(eventN)+" ")
% Remove outliers from dod1

std_factor = 2;%1.5; %boundary of N STD's from mean
t_pointsN = 0.25;%0.4; %fraction of timepoints outside of boundary to be an outlier channel
[PD_data] = rem_OD_outliers_10mins(PD_data,std_factor,t_pointsN);

% Step 5. Spectoscropy
close all;
[PD_data] = spectroscopy_DPF_10mins(PD_data);

% Remove outliers from HbO (don't have to do this step)
[PD_data] = rem_HbX_outliers_10mins(PD_data,std_factor,t_pointsN);
%currently only uses HbO to find outliers. 
%%% Downsample before tomography (by a factor of fs, fs is 10Hz, so will DS to 1Hz)
[PD_data] = PD_data_downsample(PD_data,fs);

% %%% Step 6. Image Recon
% w=0;
% %use this 26 07 23
% load('Jacobian_infant30week_850nm.mat')
% load('mesh_infant30week_850nm.mat')
[mesh] = plotjac_infant_mesh(mesh,PD_data,J,PD);
[masks] = get_J_mask(mesh,J,PD,PD_data);
%%% Tomography
[PD_data] = tomography_PD_data_10mins(mesh,J,PD_data,PD_time,PD,0);

%% Calculate R values
        switch PD_data.eventType
            case "m_hypo"
                %find min glucose pos (index will corrospond to inddex in PD_time)
                %min_glucose_pos = find(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)))  );
                
                peak_glucose(i) = min(PD.glucose(PD.events_start_end.m_hypo(PD_data.eventN,1):PD.events_start_end.m_hypo(PD_data.eventN,2)));

                %PD_time(min_glucose_pos(1)) (1) incase there are >1 min points
                %time window in mins. the -1 is for 5mins PRE MIN, the -2 and plus 8
                %make it -2mins before point and 8mins after =10mins total
                %%time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];
                %can change the (end) to (1) . (end) means it will always look at
                %the LATEST glucoce MINIMUM, (1) means it looks at the FIRST
                %glucose minimum
            case "S_hypo"
                %min_glucose_pos = find(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)) == min(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)))  );
                %%time_window_t=[PD_time(min_glucose_pos(end) )-2 PD_time(min_glucose_pos(end) )+8];

                peak_glucose(i) = min(PD.glucose(PD.events_start_end.S_hypo(PD_data.eventN,1):PD.events_start_end.S_hypo(PD_data.eventN,2)));
         end

HbO_all(:,:,i) = PD_data.HbO; 
Hb_all(:,:,i) = PD_data.Hb; 
HbT_all(:,:,i) = PD_data.HbT; 

%% ONLY LAPTOP R values for PD33 on laptop
%myFile = matfile('F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\Saved_HbT_data\PD33_BP_0_0_0067\HbT_all_m_hypo_PD33_D.mat');
%% ONLY LAPTOP m Hypo PD33
% good_m_hypo = [1 3 4 5 6 8 9 10 11 12 14 18];
% good_m_hypo = [1 3 4 5 6 8 9 10 11 12 14 18];
% loadedData = myFile.HbT_all_m_hypo(:,:,1);
% loadedData2 = myFile.HbT_all_m_hypo(:,:,[3 4 5 6]);
% loadedData3 = myFile.HbT_all_m_hypo(:,:,[8 9 10 11 12]);
% loadedData4 = myFile.HbT_all_m_hypo(:,:,14);
% loadedData5 = myFile.HbT_all_m_hypo(:,:,18);
% % S Hypo PD33
% myFile = matfile('F:\PadovaPostDoc\BabyGluCo\Study_10min_DOT\Saved_HbT_data\PD33_BP_0_0_0067\HbT_all_S_hypo_PD33_D.mat');
% loadedData6 = myFile.HbT_all_S_hypo(:,:,[2 3]);

%vol2gm*HbT_all_m_hypo

%threshold 95% max J all ch's


%%%mesh_recon = mesh;
%ind = reshape(mesh.region(mesh.elements),[],1);
%%%ind = reshape(ind>=3,[],4); %find index that are in region 3 4 5 or 6 %was
%%%using this before 06 05 24
%ind = reshape(ind==3,[],4); %find index that are in region 3 06 05 24
%ind = sum(ind,2);
%ind = find(ind==4); %sum = 4, mean all element co-ors are within 3 4 5 or 6
%[mesh3.elements,mesh3.nodes]=boundfaces(mesh.nodes,mesh.elements(ind,:),0);

%[mesh_all.elements,mesh_all.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%[meshZ.elements,meshZ.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);

%cortex_nodes = find(mesh.region>=3);

%GM_WM_nodes = find(mesh.region==3 |mesh.region==4 );
%all_J = sum(J_amp(:,:)); %this sums the jacobian for all measurements and nodes
%so we are left with the total sensitivity for EACH node

% sample_data=randn(size(mesh.nodes,1),1);
% mesh3.data = sample_data; %set mesh3.data as the total sensitivty
% 
% figure()
% plotniceimages_1_Jacobian(mesh3,mesh); %plotting the jacobian (requires function plotniceimages_1)
% 
% mesh3.data(cortex_nodes) = 5; %set mesh3.data as the total sensitivty
% mesh3.data = randn(size(cortex_nodes,1),1);
% 
% mesh3.data(GM_WM_nodes) = 5; %set mesh3.data as the total sensitivty
% mesh3.data(GM_nodes) = -5; %set mesh3.data as the total sensitivty
% mesh3.data(WM_nodes) = 0; %set mesh3.data as the total sensitivty
% mesh3.data(BS_nodes) = 5; %set mesh3.data as the total sensitivty
% mesh3.data(CB_nodes) = 10; %set mesh3.data as the total sensitivty
% figure()
% plotniceimages_1_Jacobian(mesh3,mesh); %plotting the jacobian (requires function plotniceimages_1)
% colorbar('horiz');
%


%% (optional) setting HbT data to size of mesh and then converting to gmsurface mesh
%using vol2gm transformation
HbT_all_m_hypo_allNodes = zeros(size(mesh.nodes,1),1);
HbT_all_m_hypo_allNodes(cortex_nodes) = HbT_all_m_hypo(:,1,1);
HbT_all_m_hypo_gmsurface= vol2gm*HbT_all_m_hypo_allNodes;

%plotting the results from this to see the transformation
figure; plotmesh_iso2([headVolumeMesh.node(:,1:3) HbT_all_m_hypo_allNodes],headVolumeMesh.elem(headVolumeMesh.elem(:,5) == 3,1:4))
cb = colorbar('horiz');
figure; plotmesh_iso2([gmSurfaceMesh.node HbT_all_m_hypo_gmsurface],gmSurfaceMesh.face)
cb = colorbar('horiz');
%HbT_all_m_hypo_allNodes = zeros(size(mesh.nodes,1),size(HbT_all_m_hypo,2));
%% dont run this is using gm smooth surface mesh
rho = zeros(size(cortex_nodes,1),1);
pval = zeros(size(cortex_nodes,1),1);

%% laptop map gmsurface 
%PD33
HbT_all_m_hypo_gmsurface = zeros(size(vol2gm,1),601,14);
HbT_all_allNodes = zeros(size(mesh.nodes,1),601,1);
% laptop only
HbT_all_allNodes(cortex_nodes,:) = loadedData;
HbT_all_m_hypo_gmsurface(:,:,1)= vol2gm*HbT_all_allNodes;

HbT_all_allNodes(cortex_nodes,:) = loadedData2(:,:,1);
HbT_all_m_hypo_gmsurface(:,:,2)= vol2gm*HbT_all_allNodes;
HbT_all_allNodes(cortex_nodes,:) = loadedData2(:,:,2);
HbT_all_m_hypo_gmsurface(:,:,3)= vol2gm*HbT_all_allNodes;
HbT_all_allNodes(cortex_nodes,:) = loadedData2(:,:,3);
HbT_all_m_hypo_gmsurface(:,:,4)= vol2gm*HbT_all_allNodes;
HbT_all_allNodes(cortex_nodes,:) = loadedData2(:,:,4);
HbT_all_m_hypo_gmsurface(:,:,5)= vol2gm*HbT_all_allNodes;

HbT_all_allNodes(cortex_nodes,:) = loadedData3(:,:,1);
HbT_all_m_hypo_gmsurface(:,:,6)= vol2gm*HbT_all_allNodes;
HbT_all_allNodes(cortex_nodes,:) = loadedData3(:,:,2);
HbT_all_m_hypo_gmsurface(:,:,7)= vol2gm*HbT_all_allNodes;
HbT_all_allNodes(cortex_nodes,:) = loadedData3(:,:,3);
HbT_all_m_hypo_gmsurface(:,:,8)= vol2gm*HbT_all_allNodes;
HbT_all_allNodes(cortex_nodes,:) = loadedData3(:,:,4);
HbT_all_m_hypo_gmsurface(:,:,9)= vol2gm*HbT_all_allNodes;
HbT_all_allNodes(cortex_nodes,:) = loadedData3(:,:,5);
HbT_all_m_hypo_gmsurface(:,:,10)= vol2gm*HbT_all_allNodes;

HbT_all_allNodes(cortex_nodes,:) = loadedData4;
HbT_all_m_hypo_gmsurface(:,:,11)= vol2gm*HbT_all_allNodes;

HbT_all_allNodes(cortex_nodes,:) = loadedData5;
HbT_all_m_hypo_gmsurface(:,:,12)= vol2gm*HbT_all_allNodes;

HbT_all_allNodes(cortex_nodes,:) = loadedData6(:,:,1);
HbT_all_m_hypo_gmsurface(:,:,13)= vol2gm*HbT_all_allNodes;
HbT_all_allNodes(cortex_nodes,:) = loadedData6(:,:,2);
HbT_all_m_hypo_gmsurface(:,:,14)= vol2gm*HbT_all_allNodes;

%cortex_nodes = find(mesh.region==3);
%GM_nodes = find(mesh.region==3 );
%WM_nodes = find(mesh.region==4 );
%BS_nodes = find(mesh.region==5 );
%CB_nodes = find(mesh.region==6 );
%N_nodes = size(mesh.nodes,1);

% %mesh_recon = mesh;
% ind = reshape(mesh.region(mesh.elements),[],1);
% %ind = reshape(ind>=3,[],4); %find index that are in region 3 4 5 or 6 %was
% %using this before 06 05 24
% ind = reshape(ind==3,[],4); %find index that are in region 3 06 05 24
% ind = sum(ind,2);
% ind = find(ind==4); %sum = 4, mean all element co-ors are within 3 4 5 or 6
% [meshR.elements,meshR.nodes]=boundfaces(mesh.nodes,mesh.elements(ind,:),0);
% 
% max_rho = max([max(rho) abs(min(rho))]);
% max_pval = max([max(pval) abs(min(pval))]);

%max(HbT_all_m_S_hypo(:,:,:));

% meshR.data(GM_nodes) = rho; %set mesh3.data as the total sensitivty
% figure()
% plotniceimages_1_greyJET(meshR,mesh); %plotting the jacobian (requires function plotniceimages_1)
% cb = colorbar('horiz');
% clim([-max_rho max_rho])
% ylabel(cb,"R value")
% title("R value Peak (Max or Min) Glucose vs Max HbT, Subject "+num2str(subjectN)+" "+eventType+"")
% 
% meshR.data(GM_nodes) = pval; %set mesh3.data as the total sensitivty
% figure()
% plotniceimages_1_greyJET(meshR,mesh); %plotting the jacobian (requires function plotniceimages_1)
% cb = colorbar('horiz');
% clim([-max_pval max_pval])
% ylabel(cb,"P value")
% title("P value Peak (Max or Min) Glucose vs Max HbT, Subject "+num2str(subjectN)+" "+eventType+"")

% figure()
% subplot(1,3,1)
% meshR.data(GM_nodes) = rho; %set mesh3.data as the total sensitivty
% plotniceimages_1_greyJET(meshR,mesh); %plotting the jacobian (requires function plotniceimages_1)
% cb = colorbar('horiz');
% clim([-max_rho max_rho])
% ylabel(cb,"R value")
% title("R value Peak (Max or Min) Glucose vs Max HbT, Subject "+num2str(subjectN)+" "+eventType+"")
% subplot(1,3,2)
% meshR.data(GM_nodes) = pval; %set mesh3.data as the total sensitivty
% plotniceimages_1_greyJET(meshR,mesh); %plotting the jacobian (requires function plotniceimages_1)
% cb = colorbar('horiz');
% clim([-max_pval max_pval])
% ylabel(cb,"P value")
% title("P value Peak (Max or Min) Glucose vs Max HbT, Subject "+num2str(subjectN)+" "+eventType+" TW "+time_window+"")
% subplot(1,3,3)
% meshR.data(GM_nodes) = adj_p; %set mesh3.data as the total sensitivty
% plotniceimages_1_greyJET(meshR,mesh); %plotting the jacobian (requires function plotniceimages_1)
% cb = colorbar('horiz');
% clim([-max_pval max_pval])
% ylabel(cb,"FDR cor. P value")
% title("FDR corrected P value Peak (Max or Min) Glucose vs Max HbT, Subject "+num2str(subjectN)+" "+eventType+" TW "+time_window+"")

% RUN THIS
%loading mshs file from 4D mesh's and plotting them
%load('AllMeshes_30weeks.mshs','-mat')
% figure; plotmesh_iso2(headVolumeMesh.node,headVolumeMesh.elem(headVolumeMesh.elem(:,5) == 3,1:4))
% figure; plotmesh_iso2(headVolumeMesh.node,headVolumeMesh.elem(:,1:4),'x>0');cb = colorbar('horiz'); view([-90 0]);
% figure; plotmesh_iso2(headVolumeMesh.node,headVolumeMesh.elem(:,1:4),'z>20');cb = colorbar('horiz'); view([-90 0]);
%figure; plotmesh_JAC([headVolumeMesh.node(:,1:3) sum(J.complete)'; ] ,headVolumeMesh.elem(:,1:4),'x>0');cb = colorbar('horiz');
% figure()
% plotmesh_iso2(gmSurfaceMesh.node,gmSurfaceMesh.face)
% xlabel('x / mm')
% ylabel('y / mm')
% zlabel('z / mm')
% figure()
% plotmesh_iso2(gmSurfaceMesh.node,gmSurfaceMesh.face)
% figure()
% subplot(2,2,1)
% plotmesh_iso2(gmSurfaceMesh.node,gmSurfaceMesh.face)
% view([0 90])
% xlabel('x / mm')
% ylabel('y / mm')
% zlabel('z / mm')
% subplot(2,2,2)
% plotmesh_iso2(gmSurfaceMesh.node,gmSurfaceMesh.face)
% view([0 90])
% xlabel('x / mm')
% ylabel('y / mm')
% zlabel('z / mm')
% subplot(2,2,3)
% plotmesh_iso2(gmSurfaceMesh.node,gmSurfaceMesh.face)
% view([0 90])
% xlabel('x / mm')
% ylabel('y / mm')
% zlabel('z / mm')
% subplot(2,2,4)
% plotmesh_iso2(gmSurfaceMesh.node,gmSurfaceMesh.face)
% view([0 90])
% xlabel('x / mm')
% ylabel('y / mm')
% zlabel('z / mm')
% CHANGE THIS FOR EACH SUBJECT - loading HbT data
%load('HbT_all_m_hypo_PD33_D.mat','-mat');
% 
% switch subjectN
% 
%     case 7 %add other subject Ns
%     load('HbT_all_m_hypo_PD7_D_gmsurface.mat','-mat');
%     load('peak_glucose_all_m_hypo_PD7_D.mat','-mat');
% 
%     case 8 %add other subject Ns
%     load('HbT_all_m_hypo_PD8_D_gmsurface.mat','-mat');
%     load('peak_glucose_all_m_hypo_PD8_D.mat','-mat');
% 
%     case 20 %add other subject Ns
%     load('HbT_all_m_hypo_PD20_D_gmsurface.mat','-mat');
%     load('peak_glucose_all_m_hypo_PD20_D.mat','-mat');
% 
%     case 33 %add other subject Ns
%     load('HbT_all_m_S_hypo_PD33_D_gmsurface.mat','-mat');
%     load('peak_glucose_all_m_hypo_PD33_D.mat','-mat');
%     load('peak_glucose_all_S_hypo_PD33_D.mat','-mat');
% end

% %plot J all gmsurface spatially normalized
% figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface_norm],gmSurfaceMesh.face);view(0,90);
% cb = colorbar('horiz');
% title("Spatially Normalised Jacobian - GM Smooth Week30")
% xlabel("x / mm");
% ylabel("y / mm");
% zlabel("z / mm");
% figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface_norm],gmSurfaceMesh.face,'x>0');view(0,90);
% cb = colorbar('horiz');
% figure; plotmesh_JAC([gmSurfaceMesh.node J_all_gmsurface_norm],gmSurfaceMesh.face);view(0,90);
% cb = colorbar('horiz');
%J_all_gmsurface

%% Run Mapping to gmSurfaceMesh
% HbT_all_m_hypo_gmsurface = zeros(size(vol2gm,1),size(HbT_all,2),size(HbT_all,3));
% HbT_all_allNodes = zeros(size(mesh.nodes,1),size(HbT_all,2),1);
% 
% for i=1:size(HbT_all,3)
%     HbT_all_allNodes(cortex_nodes,:) = HbT_all(:,:,i);
%     HbT_all_m_hypo_gmsurface(:,:,i)= vol2gm*HbT_all_allNodes;
% end


%% paper 21 01 25 - good channels of all subjects in m-s hypo only subgroup
HDD_text = "H";
subjectN = 49;
time_window = "A";
eventType = "m_hypo";

time_window_all = ["A","D","F"];
%mhypo
all_subject_N = [8 10 15 25 33 39 41 44 49 55 56 59 60]; %add 34 and 58

N_good_chs = 0;
N_events = 0;

min_N_good_chs = 128;
max_N_good_chs = 0;

for i=1:size(all_subject_N,2)
    subjectN = all_subject_N(i);
    for j=1:3
        time_window = time_window_all(j);
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
        N_good_chs = N_good_chs + nnz(good_ch_idx_events);
        N_events = N_events + size(good_ch_idx_events,2);
        
        for k = 1:size(good_ch_idx_events,2)
            if nnz(good_ch_idx_events(:,k)) <min_N_good_chs 
                min_N_good_chs = nnz(good_ch_idx_events(:,k));
            end

            if nnz(good_ch_idx_events(:,k)) > max_N_good_chs 
                max_N_good_chs = nnz(good_ch_idx_events(:,k));
            end
        end
    end
end
% leftover mhypo D 34 58, F 58
N_good_chs = 0;
N_events = 0;
min_N_good_chs = 128;
max_N_good_chs = 0;
eventType = "m_hypo";

        subjectN = 34;
        time_window_all = "D";
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
        N_good_chs = N_good_chs + nnz(good_ch_idx_events);
        N_events = N_events + size(good_ch_idx_events,2);
        for k = 1:size(good_ch_idx_events,2)
            if nnz(good_ch_idx_events(:,k)) <min_N_good_chs 
                min_N_good_chs = nnz(good_ch_idx_events(:,k));
            end

            if nnz(good_ch_idx_events(:,k)) > max_N_good_chs 
                max_N_good_chs = nnz(good_ch_idx_events(:,k));
            end
        end

        subjectN = 58;
        time_window_all = "D";
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
        N_good_chs = N_good_chs + nnz(good_ch_idx_events);
        N_events = N_events + size(good_ch_idx_events,2);
        for k = 1:size(good_ch_idx_events,2)
            if nnz(good_ch_idx_events(:,k)) <min_N_good_chs 
                min_N_good_chs = nnz(good_ch_idx_events(:,k));
            end

            if nnz(good_ch_idx_events(:,k)) > max_N_good_chs 
                max_N_good_chs = nnz(good_ch_idx_events(:,k));
            end
        end

        time_window_all = "F";
        load(""+HDD_text+":\PadovaPostDoc\BabyGluCo\Study_10min_DOT\All_analysis_all_subjects\Gmsurface_analysis_HbX\PD"+num2str(subjectN)+"\time_window_"+time_window+"\good_ch_idx_events_"+eventType+"_PD"+num2str(subjectN)+"_"+time_window+"_gmsurface.mat","good_ch_idx_events");
        N_good_chs = N_good_chs + nnz(good_ch_idx_events);
        N_events = N_events + size(good_ch_idx_events,2);
        for k = 1:size(good_ch_idx_events,2)
            if nnz(good_ch_idx_events(:,k)) <min_N_good_chs 
                min_N_good_chs = nnz(good_ch_idx_events(:,k));
            end

            if nnz(good_ch_idx_events(:,k)) > max_N_good_chs 
                max_N_good_chs = nnz(good_ch_idx_events(:,k));
            end
        end

% Shypo
eventType = "S_hypo";
all_subject_N = [15 25 33 49];
N_good_chs = 0;
N_events = 0;

min_N_good_chs = 128;
max_N_good_chs = 0;
for i=1:size(all_subject_N,2)
    subjectN = all_subject_N(i);
    for k = 1:size(good_ch_idx_events,2)
            if nnz(good_ch_idx_events(:,k)) <min_N_good_chs 
                min_N_good_chs = nnz(good_ch_idx_events(:,k));
            end

            if nnz(good_ch_idx_events(:,k)) > max_N_good_chs 
                max_N_good_chs = nnz(good_ch_idx_events(:,k));
            end
      end
end

%% Paper 21 05 2025 %find SNRThree
clc; close all; clearvars -except J mesh scalpSurfaceMesh gmSurfaceMesh headVolumeMesh landmarks vol2gm tenFive coverageThresh
w=0;
subjectN = 49;
eventType = "m_hypo"; %for m_hypo only or m and S combined
%eventType = "S_hypo"; %S_hypo only
eventType_all = "m-hypo";
%eventType_all = "S-hypo";
%eventType_all = "S-and-m-hypo";
time_window = "F";

switch eventType
    case "m_hypo"
        mild_and_severe = 0; %1 YES , 0 NO

    switch eventType_all
        case "m-hypo"
             mild_and_severe = 0; %1 YES , 0 NO
        case "S-and-m-hypo"
             mild_and_severe = 1; %1 YES , 0 NO
    end
    case "S_hypo"
         mild_and_severe = 0; %1 YES , 0 NO
end

t_points = 601; %10 mins (s)
weekN_J = choose_weekN_mesh_jac_infant(subjectN);
%use this 26 07 23                  `               

load("Jacobian_infant"+num2str(weekN_J)+"week_850nm.mat")
load("mesh_infant"+num2str(weekN_J)+"week_850nm.mat")
load("AllMeshes_"+num2str(weekN_J)+"weeks.mshs","-mat")

% %%% TO DO - Set mesh to correct week for subject N

gmSurfaceMesh.vol2gm = vol2gm;
[coverageThresh] = get_Jac_threshold(headVolumeMesh,vol2gm);
N_nodes = size(headVolumeMesh.node,1);
% Based upon events stats - select events to pass through to gen. HbX
%returns event numbers that are good
%good meaning > 23 good chs and <50% MA train out of all T pts
%and meet Jacobian threshold ROI
[event_good] = choose_good_events_for_tomo(eventType,subjectN,time_window);
N_eventN = size(event_good,2);

PD_meta_data = get_PD_meta_data();
clc;
switch eventType
    case "m_hypo"
        %HbT for all GM nodes in full mesh
        %HbT_all_m_hypo = zeros(size(GM_nodes,1),t_points,N_eventN);
        %HbT for all GM nodes in smooth surface mesh (save space)
        HbT_all_m_hypo_gmsurface = zeros(size(gmSurfaceMesh.node,1),t_points,N_eventN); 

        peak_glucose_all_m_hypo = zeros(N_eventN,1);
        length_glucose_all_m_hypo = zeros(N_eventN,1);
        good_ch_idx_events = zeros(128,N_eventN);
         SNR_Three_PD_all = zeros(128,N_eventN)';

        
        %for i=1:N_eventN
            %remove giacomo gdrive from path is HmrSG not working
            %[PD_data] = Iteration10minDOT(subjectN,eventType,time_window,i,J,mesh,gmSurfaceMesh,t_points,vol2gm,coverageThresh,event_good);
            %%HbO_all_m_hypo(:,:,i) = PD_data.HbO;
            %%HbT_all_m_hypo(:,:,i) = PD_data.HbT(GM_nodes,:);
            %HbT_all_m_hypo_gmsurface(:,:,i) = PD_data.HbT(:,:);
            %peak_glucose_all_m_hypo(i) = PD_data.peak_glucose;
            %length_glucose_all_m_hypo(i) = PD_data.length_glucose;
            %good_ch_idx_events(1:size(PD_data.goodch_idx,1),i) = PD_data.goodch_idx;
         %   clc;
         %   i
        %end
    
    %PD_m_hypo_events = event_good;
   
    %case "S_hypo"
        %HbO_all_S_hypo = zeros(size(GM_nodes,1),t_points,N_eventN);
        %HbT_all_S_hypo = zeros(size(GM_nodes,1),t_points,N_eventN);
        HbT_all_S_hypo_gmsurface = zeros(size(gmSurfaceMesh.node,1),t_points,N_eventN); 
        peak_glucose_all_S_hypo = zeros(N_eventN,1);
        length_glucose_all_S_hypo = zeros(N_eventN,1);
        good_ch_idx_events = zeros(128,N_eventN);

        
        for i=1:N_eventN
            %[PD_data] = Iteration10minDOT(subjectN,eventType,time_window,i,J,mesh,gmSurfaceMesh,t_points,vol2gm,coverageThresh,event_good);
            %%%%HbO_all_S_hypo(:,:,i) = PD_data.HbO;
            %HbT_all_S_hypo_gmsurface(:,:,i) = PD_data.HbT(:,:);
            %peak_glucose_all_S_hypo(i) = PD_data.peak_glucose;
            %length_glucose_all_S_hypo(i) = PD_data.length_glucose;
            %good_ch_idx_events(1:size(PD_data.goodch_idx,1),i) = PD_data.goodch_idx;
            clc;
            %i
        end
    
    %PD_S_hypo_events = event_good;
end


%
for i=1:N_eventN
    switch eventType
        case "m_hypo"
            eventN = event_good(i);
        case "S_hypo"
            eventN = event_good(i);
    end

    PD_data = load("H:\PadovaPostDoc\BabyGluCo\NIRS_data\PD"+num2str(subjectN)+"\formatted\PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(eventN)+".nirs",'-mat');
    PD_data = getfield(PD_data,"PD"+num2str(subjectN)+"_"+eventType+"_X");
    %PD_data = PD_data.PD33_m_hypo_X; %change this
    PD = load("H:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(subjectN)+".mat");
    %PD = load("E:\PadovaPostDoc\BabyGluCo\formatted_PD_GuyPerkins\PD"+num2str(subjectN)+".mat");

    PD = getfield(PD,"PD"+num2str(subjectN)+"");
    %PD = PD.PD33; %change this
    switch eventType
        case "m_hypo"
            PD_time = 0:5:5*PD.events_length.m_hypo(eventN);
        case "S_hypo"
            PD_time = 0:5:5*PD.events_length.S_hypo(eventN);
    end

    PD_data.time_window = time_window;

    PD.subjectN = subjectN;
    PD.eventN = eventN;
    PD.eventType = eventType;

    PD_data.subjectN = subjectN;
    PD_data.eventN = eventN;
    PD_data.eventType = eventType;

    % Crop data. Just want 2mins before, 8 mins after (so total 10)
    [PD_data time_window_t] = get_timewindow_10mins(PD_data,PD,PD_time,time_window);


    fs=10;
    PD_data.fs = fs;


    % Step 2. Signal Qaulity check
    dRange = [5E-4 3];
    SNRrange = 2;
    [goodch_idx, PD_data] = find_good_ch_15mins(PD_data,dRange,SNRrange,1);

    goodch_idx_step2 = PD_data.goodch_idx;

    % Step 3. Get DoD
    meanValue = mean(PD_data.d);
    dodConv = -log(abs(PD_data.d)./meanValue);
    dodConv_og = -log(abs(PD_data.d)./meanValue);
    PD_data.dod = dodConv;
    
    if size(goodch_idx,1) < 12
        PD_data.HbO = zeros(size(vol2gm,1),t_points);
        PD_data.Hb = zeros(size(vol2gm,1),t_points);
        PD_data.HbT = zeros(size(vol2gm,1),t_points);
    else
        %%%% Step 4. Motion corr. OUTPUT dod1
        methodN =1;
        %PD_data = PD_data_motioncorr_15min(PD_data,methodN,1);
        fs = PD_data.fs; %Hz

        %Do GVDT to find MA trains across entire data.
        PD_data = gvtd_10mincorr(PD_data,1);

        %Interpolate across MA trains from GVDT - GET dod_int
        PD_data = MA_inter_GVDT_10mins(PD_data,1);
        %consider factor of STD you add in noise for the interpolation

        %Do HMR to find spike/BL M.A
        tMotion = 1.0;%0.5; %0.8; % %0.5; %time range in seconds
        tMask = 2; %mark data *- time around m.a as m.a
        SDThresh = 7; %10; %12
        AmpThresh = 0.1; %0.35; %0.5;
        tIncMan = ones(length(PD_data.t ),1); % set it to vectors of ones (this is a vector used to remove manually parts of the data if needed)
        % Motion detection technique. tIncCh is a matrix number of samples x twice n of
        % channels which contains for each channel (column) the information about
        % whether an artifact was present (0s) or not (1s). tInc is a vector which
        % contains information on whether at that time sample in any of the channel
        % was present an artifact (0s) or not (1s). tInc can therefore be obtained
        % from tIncCh by setting to 0 every row that contains at least one 0.
        [tInc,tIncCh] = hmrMotionArtifactByChannel(PD_data.dod_int, fs, PD_data.SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);
        N_MA_HMR = size(find(tIncCh(:,PD_data.goodch_idx)==0),1)/size(PD_data.goodch_idx,1);

        %test Spline SG - GET dod_SG
        p=0.99;

        %do SG on dodINT
        %FrameSize_sec = 1;%0.5 %1/fs;
        FrameSize_sec = 10;

        [dod_SG ,tIncCh_baseline_dod_SG,tInc_baseline_dod_SG SNR_Three_PD] = hmrMotionCorrectSplineSG_PD_data(PD_data.dod_int, PD_data.dod_int, PD_data.t, PD_data.SD, p, FrameSize_sec,1,PD_data);
        PD_data.dod_SG = dod_SG;
        SNR_Three_PD_all(i,:) = SNR_Three_PD;


    end
    close all;
    clc;
end


save("SNR_Thre_PD_"+num2str(PD_data.subjectN)+"_"+PD_data.eventType+"_TimeW"+PD_data.time_window+"_all","SNR_Three_PD_all");

SNR_Three_PD_all_logical = zeros(size(SNR_Three_PD_all,1),128);
for i=1:size(SNR_Three_PD_all,1)
    SNR_Three_PD_all_logical(i,:) = SNR_Three_PD_all(i,:) > 3;
end
save("SNR_Thre_PD_"+num2str(PD_data.subjectN)+"_"+PD_data.eventType+"_TimeW"+PD_data.time_window+"_all_logical","SNR_Three_PD_all_logical");
