function [events_stats] = N_good_events_check_PD_10mins(subjectN,eventType,time_window,figON)

switch subjectN
    case 7
        switch eventType
            case "m_hypo"
                switch time_window
                    case "A"
                        %PD_events = [8 9 10 11 12 13 14 15];
                        PD_events = [8 9 10 11 12 13 14 15 16 17 18 19];
                    case "D"
                        %PD_events = [8 9 10 11];
                        PD_events = [8 9 10 11 12 13 14 15 16 17 18 19];
                    case "F"
                        %PD_events = [8 10 11 12 14 15 16];
                        PD_events = [8 9 10 11 12 13 14 15 16 17 18 19];
                end
        end
    case 8 %add other subject Ns
        switch eventType
            case "m_hypo"
                switch time_window
                    case "A"
                        %PD_events = [1 8 9 10];
                        PD_events = [1 2 3 4 5 6 7 8 9 10 11];
                    case "D"
                        %PD_events = [7 8 9 10];
                        PD_events = [1 2 3 4 5 6 7 8 9 10 11];
                    case "F"
                        %PD_events = [1 4 5 7 8 9 10 11];
                        PD_events = [1 2 3 4 5 6 7 8 9 10 11];
                end
                %PD_events = [1 2 3 4 5 6 7 8 9 10 11];
        end
    case 10
        switch eventType
            case "m_hypo"
                PD_events = [1];
            end
        
    case 13 %add other subject Ns
        switch eventType
            case "m_hypo"
                switch time_window
                    case "A"
                        PD_events = [5];
                    case "D"
                        PD_events = [5];
                    case "F"
                        PD_events = [5];
                end
                %PD_events = [1 2 3 4 5 6 7 8 9 10 11];
        end
    case 15
        switch eventType
            case "m_hypo"
                switch time_window
                    case "A"
                        PD_events = [1 2 3 4 5 6 7 8];
                    case "D"
                        PD_events = [1 2 3 4 5 6 7 8];
                    case "F"
                        PD_events = [1 2 3 4 5 6 7 8];
                end
                %PD_events = [1 2 3 4 5 6 7 8];
            case "S_hypo"
                switch time_window
                    case "A"
                        PD_events = [1 2];
                    case "D"
                        PD_events = [1 2];
                    case "F"
                        PD_events = [1 2];
                end
                %PD_events = [1 2];
        end
    case 20 %add other subject Ns
        switch eventType
            case "m_hypo"   
                switch time_window  
                    case "A"
                        %PD_events = [3 4 5 6 7 8 9 11];
                        PD_events = [3 4 5 6 7 8 9 10 11 12];
                    case "D"
                        PD_events = [3 4 5 6 7 8 9 10 11 12];
                        %PD_events = [3 4 5 6 7 9 10 11 12];
                    case "F"
                        PD_events = [3 4 5 6 7 8 9 10 11 12];
                        %PD_events = [3 4 5 6 7 9 10 11];
                end

          end

    case 25 %add other subject Ns
        switch eventType
            case "m_hypo"   
                switch time_window  
                    case "A"
                        PD_events = [1 2 3 4 5];
                    case "D"
                        %PD_events = [3 4 5 6 7 8 9 10 11 12];
                        %PD_events = [2 3 5];
                        PD_events = [1 2 3 4 5];
                    case "F"
                        %PD_events = [3 4 5 6 7 8 9 10 11 12];
                        PD_events = [1 2 3 4 5];
                        %PD_events = [1 2 3 5];
                end
            
             case "S_hypo"   
                switch time_window  
                    case "A"
                        PD_events = [1];
                    case "D"
                        PD_events = [1];
                    case "F"
                        PD_events = [1];
                end
        end
     case 28
        switch eventType
            case "m_hypo"
                PD_events = [1 2];
            end
        
     case 33
        switch eventType
            case "m_hypo"
                switch time_window  
                    case "A"
                        PD_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
                    case "D"
                        %PD_events = [3 4 5 6 7 8 9 10 11 12];
                        PD_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];

                    case "F"
                        %PD_events = [3 4 5 6 7 8 9 10 11 12];
                        PD_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
                end
                %PD_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
            case "S_hypo"
                switch time_window  
                    case "A"
                        PD_events = [1 2 3];
                    case "D"
                        %PD_events = [3 4 5 6 7 8 9 10 11 12];
                        PD_events = [1 2 3];
                    case "F"
                        %PD_events = [3 4 5 6 7 8 9 10 11 12];
                        PD_events = [1 2 3];
                end
                %PD_events = [1 2 3];
        end

        case 34
        switch eventType
            case "m_hypo"
                PD_events = [1];
            end
        
     case 39 %add other subject Ns
        switch eventType
            case "m_hypo"   
                switch time_window  
                    case "A"
                        PD_events = [1 2];
                    case "D"
                        PD_events = [1 2];
                    case "F"
                        PD_events = [1 2];
                end

          end
     case 41 %add other subject Ns
        switch eventType
            case "m_hypo"   
                switch time_window  
                    case "A"
                        PD_events = [2 3 4];
                    case "D"
                        PD_events = [2 3 4];
                    case "F"
                        PD_events = [2 3 4];
                end

          end

    case 44 %add other subject Ns
        switch eventType
            case "m_hypo"   
                switch time_window  
                    case "A"
                        PD_events = [1];
                    case "D"
                        PD_events = [1];
                    case "F"
                        PD_events = [1];
                end

          end

     case 45 %add other subject Ns
        switch eventType
            case "m_hypo"   
                switch time_window  
                    case "A"
                        PD_events = [1 2];
                    case "D"
                        PD_events = [1 2];
                    case "F"
                        PD_events = [1 2];
                end

          end

    case 48 %add other subject Ns
        switch eventType
            case "m_hypo"   
                switch time_window  
                    case "A"
                        PD_events = [1 2];
                    case "D"
                        PD_events = [1 2];
                    case "F"
                        PD_events = [1 2];
                end

          end

    case 49
        switch eventType
            case "m_hypo"
                switch time_window  
                    case "A"
                        %PD_events = [1 2 3 4 7 8 9 10 11 12 13 14 15 16 18 22 23 24 25 26 27 28];
                        PD_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];
                    case "D"
                        PD_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];
                    case "F"
                        %PD_events = [1 2 3 4 6 7 8 9 10 11 12 14 15 17 18 19 20 22 26 27 28 29];
                        PD_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];
                end
            case "S_hypo"
                switch time_window  
                    case "A"
                        PD_events = [1 2 3 4 5];
                    case "D"
                        %PD_events = [3 4 5 6 7 8 9 10 11 12];
                        PD_events = [1 2 3 4 5];
                    case "F"
                        %PD_events = [3 4 5 6 7 8 9 10 11 12];
                        PD_events = [1 2 3 4 5];
                end
                %PD_events = [1 2 3];
         end


        case 55 %add other subject Ns
            switch eventType
                case "m_hypo"   
                switch time_window  
                    case "A"
                        PD_events = [1 2 3 4];
                    case "D"
                        PD_events = [1 2 3 4];
                    case "F"
                        PD_events = [1 2 3 4];
                end

                end

        case 56
        switch eventType
            case "m_hypo"
                PD_events = [1 2 3 4 5];
            end

        case 58
        switch eventType
            case "m_hypo"
                PD_events = [1];
            end

        case 59
        switch eventType
            case "m_hypo"
                PD_events = [1 2];
            end

        case 60
        switch eventType
            case "m_hypo"
                PD_events = [1 2 3 4 5 6 7 8 9 10 11];
            end
        




end




        if figON == 0
            set(gcf,'Visible','off');              
            set(0,'DefaultFigureVisible','off');

        elseif figON == 1
            set(gcf,'Visible','on');              
            set(0,'DefaultFigureVisible','on');
        end

 savefig = 1;

 PD = load("H:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(subjectN)+".mat");
 PD = getfield(PD,"PD"+num2str(subjectN)+"");
 %n_events = getfield(PD.N_events,""+eventType+"");
n_events = size(PD_events,2);

fs = 10;
PD_data.fs = fs;

eventN = PD_events(1);
PD_data = load("PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(eventN)+".nirs",'-mat');
PD_data = getfield(PD_data,"PD"+num2str(subjectN)+"_"+eventType+"_X");

N_MA_trains = [];
pc_Tpts_MA_train = [];
N_MA_hmr_ch = zeros(size(PD_data.t,2),1);
pc_Tpts_MA_hmr_ch = zeros(size(PD_data.t,2),1);
goodch_idx_allE = zeros(size(PD_data.d,2), n_events);
eventN_allE = zeros(n_events,1);
eventN_allE = PD_events;

N_frames_MA_trains_allE = zeros(size(PD_data.d,2),n_events);

            % N_MA_trains_allE = zeros(N_MA_trains;
            % pc_Tpts_MA_train_allE = pc_Tpts_MA_train;
            % N_MA_hmr_ch_allE = N_MA_hmr_ch;
            % pc_Tpts_MA_hmr_ch_allE = pc_Tpts_MA_hmr_ch;
%%

    for i=1:n_events 
            i
            eventN = PD_events(i);
            %"event number"
            % change the names below for diff subject/event type
            PD_data = load("PD"+num2str(subjectN)+"_"+eventType+"_"+num2str(eventN)+".nirs",'-mat');
            PD_data = getfield(PD_data,"PD"+num2str(subjectN)+"_"+eventType+"_X");
            PD_data.subjectN = subjectN;
            PD_data.eventN = eventN;
            PD_data.eventType = eventType;
        
            PD = load("H:\PadovaPostDoc\BabyGluCo\formatted_PD_V2_GuyPerkins170424\PD"+num2str(subjectN)+".mat");
            PD = getfield(PD,"PD"+num2str(subjectN)+"");
            PD.subjectN = subjectN;
            PD.eventN = eventN;
            PD.eventType = eventType;
            
            switch PD_data.eventType
                case "m_hypo"
                    PD_time = 0:5:5*PD.events_length.m_hypo(eventN);
                case "S_hypo"
                    PD_time = 0:5:5*PD.events_length.S_hypo(eventN);
            end
            %get time window
            % Crop data. Just want 15mins before, 15 mins after (so total 30)
            [PD_data time_window_t] = get_timewindow_10mins(PD_data,PD,PD_time,time_window);
            PD_data.time_window = time_window;

            PD_data.subjectN = subjectN;
            PD_data.eventN = eventN;
            PD_data.eventType = eventType;

            % Step 2. Signal Qaulity check
            dRange = [5E-4 3];%3];
            SNRrange = 2; %2;
            [goodch_idx, PD_data] = find_good_ch_15mins(PD_data,dRange,SNRrange,0);
            PD_data.goodch_idx = goodch_idx;

            % Step 3. Get DoD
            meanValue = mean(PD_data.d);
            dodConv = -log(abs(PD_data.d)./meanValue);
            dodConv_og = -log(abs(PD_data.d)./meanValue);
            PD_data.dod = dodConv;
            
            %if number goodchs < 4, skip the MA train
            if size(goodch_idx,1)>3
                %Do GVDT to find MA trains across entire data.
                PD_data = gvtd_10mincorr(PD_data,0);
   
                %Interpolate across MA trains from GVDT - GET dod_int
                PD_data = MA_inter_GVDT_10mins(PD_data,0);
                %consider factor of STD you add in noise for the interpolation
            end

            if size(goodch_idx,1)<4
                PD_data.dod_int = PD_data.dod;
            end
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
            % figure()
            % plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(1)))
            % hold on
            % plot(PD_data.t(find(tIncCh(:,PD_data.goodch_idx(1))==0)),PD_data.dod_int(find(tIncCh(:,PD_data.goodch_idx(1))==0),PD_data.goodch_idx(1)),'m.','MarkerSize',6)
            % xlabel('Time / s')
            % ylabel('dOD / A.U')
            % title("Use Hmr YANG 22 for MAD "+tMotion+" |"+tMask+" |"+SDThresh+" |"+AmpThresh+" |N M.A = "+N_MA_HMR+" ")

            % Get stats for each event
            
            %if size(PD_data.MA_train_GVDT_new.start_train_frame,1) == 0 %if no MA trains, set values to zero and don't do calcs.
            if ~isfield(PD_data,"MA_train_GVDT_new") %Field doesn't exist = TRUE
                N_MA_trains = 0;
                pc_Tpts_MA_train = 0;
                N_MA_hmr_intrains = 0;
                N_Tpts_trains = 0;
                pc_Tpts_MA_hmr_ch = 0;
                N_frames_MA_trains = 0;

                %3. N M.A ch OUTSIDE OF Trains
                N_MA_hmr_ch = size(PD_data.t,1)- sum(tIncCh); %N MA tpts ch in total
                
                N_MA_hmr_ch = N_MA_hmr_ch - N_MA_hmr_intrains; %N MA tpts outside of trains
                
                %4. % of Tpts that are MA Ch outside of trains
                N_tpts_nottrains = size(PD_data.t,1);
                pc_Tpts_MA_hmr_ch = N_MA_hmr_ch ./N_tpts_nottrains;
            else
                %1. Number of MA trains
                N_MA_trains = size(PD_data.MA_train_GVDT_new.start_train_frame,1);
                %2. % of TPts that are MA trains
                pc_Tpts_MA_train = sum( PD_data.MA_train_GVDT_new.end_train_frame(:) - PD_data.MA_train_GVDT_new.start_train_frame(:) +1 )   /size(PD_data.t,1);
                %3. N M.A ch OUTSIDE OF Trains
                N_MA_hmr_ch = size(PD_data.t,1)- sum(tIncCh); %N MA tpts ch in total
                N_MA_hmr_intrains = size(PD_data.MA_train_GVDT_new.start_train_frame(:):PD_data.MA_train_GVDT_new.end_train_frame(:),2) - sum( tIncCh(PD_data.MA_train_GVDT_new.start_train_frame(:):PD_data.MA_train_GVDT_new.end_train_frame(:),:)); %N MA Tpts in trains only
                N_MA_hmr_ch = N_MA_hmr_ch - N_MA_hmr_intrains; %N MA tpts outside of trains
                %figure()
                %plot(PD_data.goodch_idx,N_MA_hmr_ch(1,PD_data.goodch_idx),'b*-')

                %4. % of Tpts that are MA Ch outside of trains
                N_Tpts_trains = size(PD_data.MA_train_GVDT_new.start_train_frame(:):PD_data.MA_train_GVDT_new.end_train_frame(:),2);
                N_tpts_nottrains = size(PD_data.t,1) - N_Tpts_trains;
                pc_Tpts_MA_hmr_ch = N_MA_hmr_ch ./N_tpts_nottrains; 

                %5. length of each train
                N_frames_MA_trains = PD_data.MA_train_GVDT_new.end_train_frame-PD_data.MA_train_GVDT_new.start_train_frame;

             end
            N_MA_trains_allE(i,:) = N_MA_trains;
            pc_Tpts_MA_train_allE(i,:) = pc_Tpts_MA_train;
            N_MA_hmr_ch_allE(i,:) = N_MA_hmr_ch;
            pc_Tpts_MA_hmr_ch_allE(i,:) = pc_Tpts_MA_hmr_ch;
            goodch_idx_allE(1:size(PD_data.goodch_idx,1),i) = PD_data.goodch_idx;

            N_frames_MA_trains_allE(1:size(N_frames_MA_trains,1),i) = N_frames_MA_trains; %each row is the N of frames for a MA train % each col is the event N
       
            % figure()
            % subplot(2,1,1)
            % imagesc(N_MA_hmr_ch(PD_data.goodch_idx)')
            % colorbar, xlabel('Time / samples'), ylabel('Measurement #'),title("GVDT new GoodCh medSNR lambda "+nStd+" ")
            % 
            % x = 1:13;
            % data = [37.6 24.5 14.6 18.1 19.5 8.1 28.5 7.9 3.3 4.1 7.9 1.9 4.3]';
            % errhigh = [2.1 4.4 0.4 3.3 2.5 0.4 1.6 0.8 0.6 0.8 2.2 0.9 1.5];
            % errlow  = [4.4 2.4 2.3 0.5 1.6 1.5 4.5 1.5 0.4 1.2 1.3 0.8 1.9];
            % figure()
            % barchart_idx = [1; 2; 3; 4];
            % barchart_data = [N_MA_trains pc_Tpts_MA_train mean(N_MA_hmr_ch(PD_data.goodch_idx)) mean(pc_Tpts_MA_hmr_ch(PD_data.goodch_idx))];
            % bar(barchart_idx,barchart_data)                
            % 
            % hold on
            % 
            % er = errorbar(x,data,errlow,errhigh);    
            % er.Color = [0 0 0];                            
            % er.LineStyle = 'none';  
            % 
            % hold off
            

    end
%%
    %create figures
    %1. N MA Trains
    figure()
    %barchart_idx =  [1:1:size(N_MA_trains_allE,1)]';
    barchart_idx = PD_events;
    barchart_data = [N_MA_trains_allE];
    bar(barchart_idx,barchart_data)
    yline(5,'r--')
    xlabel("Event N")
    ylabel("N MA Trains")
    title("N MA. Trains - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    
    %2. %Tpts in In MA Trains and %Tpts MA out of trains
    %Need to get mean for each event, since each event has different
    %goodch idx
    for i=1:size(pc_Tpts_MA_hmr_ch_allE,1)
        pc_Tpts_MA_hmr_ch_allE_mean(i) = mean(pc_Tpts_MA_hmr_ch_allE (i, nonzeros(goodch_idx_allE(:,i)) ));
        pc_Tpts_MA_hmr_ch_allE_std(i) = std(pc_Tpts_MA_hmr_ch_allE (i, nonzeros(goodch_idx_allE(:,i)) ));
        %N_MA_hmr_ch_allE_mean(i)
    end
    pc_Tpts_MA_hmr_ch_allE_std=pc_Tpts_MA_hmr_ch_allE_std*100;

    figure()
    %barchart_idx =  [1:1:size(N_MA_trains_allE,1)]';
    barchart_idx = PD_events;
    barchart_data = [pc_Tpts_MA_train_allE pc_Tpts_MA_hmr_ch_allE_mean' ]*100;
    bar(barchart_idx,barchart_data)%,%'stacked')
    yline(30,'r--')
    ylim([0 100])
    xlabel("Event N")
    ylabel("% Tpts")
    title("%Tpts in MA Trains & %Tpts MA outside of Trains - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    legend('In Train w.r.t all t','Out train w.r.t out train t')
    %add error bar to this to show STD of the mean of pc_Tpts_MA_hmr_ch_allE
    % hold on
    % er = errorbar(barchart_idx,barchart_data,barchart_data-pc_Tpts_MA_hmr_ch_allE_std,barchart_data+pc_Tpts_MA_hmr_ch_allE_std);    
    % er.Color = [0 0 0];                            
    % er.LineStyle = 'none';  

    %3. N. MA's outside of trains per channel
    %N_MA_hmr_ch_allE(1,nonzeros(goodch_idx_allE(:,1)))
    
    figure()
    boxplot(N_MA_hmr_ch_allE(1,nonzeros(goodch_idx_allE(:,1))   )./nnz(goodch_idx_allE(:,1)) ,1,'Positions',1)
    if n_events > 1;
        for i=2:n_events
        hold on
        boxplot(N_MA_hmr_ch_allE(i,nonzeros(goodch_idx_allE(:,i)))./nnz(goodch_idx_allE(:,i)) ,i,'Positions',i)
        end
    end
    title("N MA HmR in Good Ch idx OUTSIDE TRAINS DIV BY N GOOD CHs - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    xlabel('Event N'), ylabel('N. MA DIV N Good Chs')

    %4. Stacked bar graphs - length of MA trains and number
    N_frames_MA_trains_allE;
    figure()
    bar(barchart_idx,N_frames_MA_trains_allE','stacked')
    yline(0.4*size(PD_data.t,1),'r--')
    title("N & Length MA Trains (frames) Good Ch idx - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    xlabel('Event N')
    ylabel('Length of MA trains (frames)')
    ylim([0 size(PD_data.t,1)])

    %show 2 and 4 on subplot
    figure()
    subplot(2,1,1)
    barchart_data = [pc_Tpts_MA_train_allE pc_Tpts_MA_hmr_ch_allE_mean' ]*100;
    bar(barchart_idx,barchart_data)%,%'stacked')
    yline(30,'r--')
    ylim([0 100])
    xlabel("Event N")
    ylabel("% Tpts")
    title("%Tpts in MA Trains & %Tpts MA outside of Trains - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    legend('In Train w.r.t all t','Out train w.r.t out train t')
    subplot(2,1,2)
    bar(barchart_idx,N_frames_MA_trains_allE','stacked')
    yline(0.4*size(PD_data.t,1),'r--')
    title("N & Length MA Trains (frames) Good Ch idx - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    xlabel('Event N')
    ylabel('Length of MA trains (frames)')
    ylim([0 size(PD_data.t,1)])

    %show 1 and 4 on subplot
    figure()
    subplot(2,1,1)
    %barchart_idx =  [1:1:size(N_MA_trains_allE,1)]';
    barchart_idx = PD_events;
    barchart_data = [N_MA_trains_allE];
    bar(barchart_idx,barchart_data)
    yline(5,'r--')
    xlabel("Event N")
    ylabel("N MA Trains")
    title("N MA. Trains - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    subplot(2,1,2)
    bar(barchart_idx,N_frames_MA_trains_allE','stacked')
    yline(0.4*size(PD_data.t,1),'r--')
    title("N & Length MA Trains (frames) Good Ch idx - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    xlabel('Event N')
    ylabel('Length of MA trains (frames)')
    ylim([0 size(PD_data.t,1)])

    %show 2 4 1 on same plot
    figure()
    subplot(3,1,1)
    barchart_data = [pc_Tpts_MA_train_allE pc_Tpts_MA_hmr_ch_allE_mean' ]*100;
    bar(barchart_idx,barchart_data)%,%'stacked')
    yline(50,'r--')
    ylim([0 100])
    xlabel("Event N")
    ylabel("% Tpts")
    title("%Tpts in MA Trains & %Tpts MA outside of Trains - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    legend('In Train w.r.t all t','Out train w.r.t out train t')
    subplot(3,1,2)
    bar(barchart_idx,N_frames_MA_trains_allE','stacked')
    yline(0.5*size(PD_data.t,1),'r--')
    title("N & Length MA Trains (frames) Good Ch idx - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    xlabel('Event N')
    ylabel('Length of MA trains (frames)')
    ylim([0 size(PD_data.t,1)])
    subplot(3,1,3)
    %barchart_idx =  [1:1:size(N_MA_trains_allE,1)]';
    barchart_idx = PD_events;
    barchart_data = [N_MA_trains_allE];
    bar(barchart_idx,barchart_data)
    yline(5,'r--')
    xlabel("Event N")
    ylabel("N MA Trains")
    title("N MA. Trains - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    set(gcf,'Visible','on');              
    set(0,'DefaultFigureVisible','on');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD_data.subjectN)+"_"+PD_data.eventType+"_"+num2str(PD_data.eventN)+"_TimeW_"+PD_data.time_window+"_event_MA_check.png")
    end

    n_good_chs=zeros(1,n_events);
    for k=1:n_events
        n_good_chs(1,k) = nnz(goodch_idx_allE(:,k));
    end

    figure()
    bar(PD_events,n_good_chs)
    ylim([0 128])
    yline(25,'r--')
    xlabel('Event N')
    ylabel('N Good Chs')
    title("N Good Chs - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if savefig == 1
        saveas(gcf,"PD"+num2str(PD_data.subjectN)+"_"+PD_data.eventType+"_"+num2str(PD_data.eventN)+"_TimeW_"+PD_data.time_window+"_N_good_chs.png")
    end

    events_stats.N_MA_trains_allE = N_MA_trains_allE;
    events_stats.pc_Tpts_MA_train_allE = pc_Tpts_MA_train_allE; 
    events_stats.pc_Tpts_MA_hmr_ch_allE_mean = pc_Tpts_MA_hmr_ch_allE_mean;
    events_stats.N_MA_hmr_ch_allE = N_MA_hmr_ch_allE;
    events_stats.N_frames_MA_trains_allE = N_frames_MA_trains_allE;
    events_stats.goodch_idx_allE=goodch_idx_allE;
    events_stats.eventN_allE = eventN_allE;
    
    % figure()
    % imagesc(N_MA_hmr_ch_allE(:,PD_data.goodch_idx))
    % colorbar, xlabel('Good Ch Idx N'), ylabel('Event N'),title("N MA HmR in Good Ch idx OUTSIDE TRAINS - Time W "+time_window+" PD "+num2str(subjectN)+" "+eventType+"")

%%%
% %%
% [PD_m_hypo_events] = get_event_numbers(subjectN,eventType)
% switch subjectN
%     case 7 %add other subject Ns
%         switch eventType
%             case "m_hypo"
%                 %PD_m_hypo_events = [8 9 10 11 12 13 14 15 16 17 18 19];                
%                 %PD_m_hypo_events = [8 9 10 11 13 14 17];%17 too long ma
%                 %train - int not work
%                 %PD_m_hypo_events = [8 9 10 11 13 14];
%                 %all below
%                 PD_m_hypo_events = [8 9 10 11 12 13 14 15 16 17 18 19];
%                 %8 9 10 11 13 14 %only good events
%                 N_eventN = size(PD_m_hypo_events,2);
%         end
%     case 8 %add other subject Ns
%         switch eventType
%             case "m_hypo"
%                 %PD_m_hypo_events = [1 2 3 4 5 6 7 8 9 10 11];
%                 PD_m_hypo_events = [1 2 3 4 5 6 7 8 9 10 11];
%                 N_eventN = size(PD_m_hypo_events,2);
%         end
%     case 13
%         switch eventType
%             case "m_hypo"
%                 PD_m_hypo_events = [5];
%                 N_eventN = size(PD_m_hypo_events,2);
% 
%         end
%     case 15
%         switch eventType
%             case "m_hypo"
%                 %PD_m_hypo_events = [1 2 3 4 5 6 7 8];
%                 PD_m_hypo_events = [1 2 3 4 5 6 7 8];
%                 N_eventN = size(PD_m_hypo_events,2);
%             case "S_hypo"
%                 PD_S_hypo_events = [1 2];
%                 N_eventN = size(PD_S_hypo_events,2);
%         end
%     case 20 %add other subject Ns
%         switch eventType
%             case "m_hypo"
%                 %PD_m_hypo_events = [3 4 5 6 7 8 9 10 11 12];
%                 %PD_m_hypo_events = [3 4 5 6 7 9 11 12];
%                 %N_eventN = size(PD_m_hypo_events,2);
%                 PD_m_hypo_events = [3 4 5 6 7 8 9 10 11 12];
% 
%                 %switch time_window  
%                     %case "A"
%                     %PD_m_hypo_events = [3 4 5 6 8 9 11];
% 
%                     %case "D"
%                     %PD_m_hypo_events = [3 4 5 6 7 8 9 10 11 12];
%                     %PD_m_hypo_events = [3 4 5 6 7 9 11 12];
%                 %end
%             N_eventN = size(PD_m_hypo_events,2);
%         end
%     case 25
%         switch eventType
%             case "m_hypo"
%                 PD_m_hypo_events = [1 2 3 4 5];
%                 N_eventN = size(PD_m_hypo_events,2);
%             case "S_hypo"
%                 PD_S_hypo_events = [1];
%                 N_eventN = size(PD_S_hypo_events,2);
%         end
%     case 33 %add other subject Ns
%         switch eventType
%             case "m_hypo"
%                 %switch time_window
%                 PD_m_hypo_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
%                 N_eventN = size(PD_m_hypo_events,2);
% 
%                     %case "A"
%                     %    PD_m_hypo_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
% 
%                     %case "D"
%                     %    PD_m_hypo_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
%                     %    PD_m_hypo_events = [1 3 4 6 7 9 10 11 12 13 15 19]; %based on good idx 15 05 24
% 
% 
%                 %N_eventN = size(PD_m_hypo_events,2);
%             case "S_hypo"
%                 PD_S_hypo_events = [1 2 3]; %S hypo
%                 %PD_S_hypo_events = [2 3]; %S hypo  %based on good idx 15 05 24
%                 N_eventN = size(PD_S_hypo_events,2);
%         end
%     case 49
%         switch event_type
%             case "m_hypo"
%              PD_m_hypo_events = [1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];
%              N_eventN = size(PD_m_hypo_events,2);
% 
%             case "S_hypo"
%             PD_S_hypo_events = [1 2 3 4 5]; %S hypo
%             N_eventN = size(PD_S_hypo_events,2);
%         end
% end
    
end