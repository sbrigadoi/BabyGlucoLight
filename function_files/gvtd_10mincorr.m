function [PD_data] = gvtd_10mincorr(PD_data,figON)

        if figON == 0
            set(gcf,'Visible','off');              
            set(0,'DefaultFigureVisible','off');

        elseif figON == 1
            set(gcf,'Visible','on');              
            set(0,'DefaultFigureVisible','on');
        end

        g = zeros(size(PD_data.dod,1),size(PD_data.dod,2));
         for i=2:size(PD_data.dod,1)
                for j=1:size(PD_data.dod,2)
                    g(i,j)=sqrt( (PD_data.dod(i,j)-PD_data.dod(i-1,j)).^2 );
                end
         end

        %abr = mean(g);
        %gvtdTimeTrace = gvtd(PD_data.dod2_medSNR'); %input ch x t, output t x 1
        %use only 'good channels' idx for GVDT distribution
        gvtdTimeTrace = gvtd(PD_data.dod(:,PD_data.goodch_idx)'); %input ch x t, output t x 1
        

        %figure()
        %plot(PD_data.t,gvtdTimeTrace_medSNR)
        %ylabel('GVTD (global RMS)/ A.U')
        %xlabel('Time / s')
        %title('GVTD of good channels')

        statType = StatType.Histogram_Mode;
        nStd = 4;
        thresh = find_gvtd_thresh(gvtdTimeTrace, statType, nStd);
        
        gvdt_more_thresh = find(gvtdTimeTrace>thresh);
        N_MA_GVDT = size(gvdt_more_thresh,1);
        
        figure()
        histogram(gvtdTimeTrace)
        hold on
        xline([thresh],LineWidth=4)
        xlabel('GVTD value / A.U')
        ylabel('Frequency / N. of timepoints')
        title("Histogram GVTD nSTD = "+nStd+" yields threshold of "+thresh+" with "+N_MA_GVDT+" M.As")

        %set MA markers. 1= NO MA, 0 = MA
        tIncCh = ones(size(PD_data.dod,1)  , size(PD_data.dod,2)  );
        tIncCh(gvdt_more_thresh,:) = 0;

        N_MA_HMR = size(find(tIncCh(:,PD_data.goodch_idx)==0),1)/size(PD_data.goodch_idx,1);
        figure()
        plot(PD_data.t,PD_data.dod(:,PD_data.goodch_idx(1) ))
        hold on
        plot(PD_data.t(find(tIncCh(:,PD_data.goodch_idx(1) )==0)),PD_data.dod(find(tIncCh(:,PD_data.goodch_idx(1) )==0), PD_data.goodch_idx(1)),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("MAD Using GVDT nStd = "+nStd+" |N M.A = "+N_MA_HMR+" Ch"+num2str(PD_data.goodch_idx(1))+" ")

        % n_goodchs = size(PD_data.goodch_idx,1);
        % if  n_goodchs<4
        %     n_goodchs = 2
        % 
        % end

         figure()
         sgtitle("GVTD lambda "+nStd+" ")
         for i=1:4
            subplot(3,4,i)
            plot(PD_data.t, PD_data.dod(:,PD_data.goodch_idx(i)))
            hold on
            plot(PD_data.t(find(tIncCh(:,PD_data.goodch_idx(i))==0)),PD_data.dod(find(tIncCh(:,PD_data.goodch_idx(i))==0),PD_data.goodch_idx(i)),'m.','MarkerSize',6)
            %xline([900],'b--')
            %xline([5100],'b--')
            title("ch "+num2str(PD_data.goodch_idx(i)))
            xlabel('T / s')
            ylabel('dod / AU')
          end

        %gvtd_medfilt_medSNR = medfilt1(gvtdTimeTrace_medSNR,10);
        figure()
        subplot(2,1,1)
        imagesc(g(:,PD_data.goodch_idx)')
        colorbar, xlabel('Time / samples'), ylabel('Measurement #'),title("GVDT GoodChIdx")
        subplot(2,1,2)
        plot(PD_data.t,gvtdTimeTrace)
        hold on
        plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
        xlim([PD_data.t(1) PD_data.t(end)])
        xlabel('Time / S')
        ylabel('GVTD / A.U')

        %%% find MA trains using first GVDT run.
        %tIncCh = ones(size(PD_data.dod,1)  , size(PD_data.dod,2)  );
        %tIncCh(gvdt_more_thresh,:) = 0;

        length_train = 200; %frames 
        IncCh_train = ones(1,size(tIncCh,2));
        IncCh_train_loc = ones(size(tIncCh,1),size(tIncCh,2));
        for i=1:size(tIncCh,2)
            if size(unique(tIncCh(:,i)),1) > 1
                for j=1: size(tIncCh,1)-length_train %looking forward
                    if nnz(tIncCh(j:j+(length_train-1),i)) < (length_train*0.9) % IF >20f/200f are MA, train=true
                        IncCh_train(1,i) = 0; %train TRUE (0)
                        IncCh_train_loc(j,i) = 0; %train TRUE (0)
                    end
                end
                for j=1: size(tIncCh,1)-length_train %looking backward
                    if nnz(tIncCh(end-j-(length_train-1):end-j,i)) < (length_train*0.9) % IF >15f/100f are MA, train=true
                        IncCh_train(1,i) = 0; %train TRUE (0)
                        IncCh_train_loc(end-j,i) = 0; %train TRUE (0)
                    end
                end
            end
        end

        IncCh_train = IncCh_train';
        MA_train_idx = IncCh_train;

        for i=1:size(MA_train_idx,1)
            if MA_train_idx(i,1) == 1;
                MA_train_idx(i,1) = i;
            end
        end

        MA_train_idx = nonzeros(MA_train_idx);
        
        IncCh_train_loc_delta = find(IncCh_train_loc(:,1)==0);
        
        %if there are no MA, skip MA analysis
        if size(IncCh_train_loc_delta,1) == 0 

        else
            %if there are many MA trains, else only 1 MA train
            %MA train the same train if gap < 100 frames between them.
            if size(find(diff(find(IncCh_train_loc(:,1)==0))>100),1)>0
                %end_MA_train(:,1) = find(diff(find(IncCh_train_loc(:,1)==0))>100);
                end_MA_train(:,1) = [find(diff(IncCh_train_loc_delta)>100); size(IncCh_train_loc_delta,1)]; %so it gets end of the last train   
            else
                end_MA_train(1,1) = size(IncCh_train_loc_delta,1);
            end
            start_MA_train(:,1) = [1 ; end_MA_train(1:end-1,1)+1];
            %location of start/end of MA train in data (frame N)
            end_MA_train_locations_frame = IncCh_train_loc_delta(end_MA_train,1);
            start_MA_train_locations_frame = IncCh_train_loc_delta(start_MA_train,1);

            figure()
            sgtitle("SNR GVTD lambda "+nStd+" ")
            for i=1:4
                subplot(3,4,i)
                plot(PD_data.t, PD_data.dod(:,PD_data.goodch_idx(i)))
                hold on
                plot(PD_data.t(find(tIncCh(:,PD_data.goodch_idx(i))==0)),PD_data.dod(find(tIncCh(:,PD_data.goodch_idx(i))==0),PD_data.goodch_idx(i)),'m.','MarkerSize',6)
                plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),PD_data.dod(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
                plot(PD_data.t(start_MA_train_locations_frame),thresh,'g*','MarkerSize',8)
                plot(PD_data.t(end_MA_train_locations_frame),thresh,'b*','MarkerSize',8)
                %xline([900],'b--')
                %xline([5100],'b--')
                title("ch "+num2str(PD_data.goodch_idx(i)))
                xlabel('T / s')
                ylabel('dod / AU')
            end

            figure()
            subplot(2,1,1)
            imagesc(g(:,PD_data.goodch_idx)')
            colorbar, xlabel('Time / samples'), ylabel('Measurement #'),title("GVDT GoodCh medSNR lambda "+nStd+" ")
            subplot(2,1,2)
            plot(PD_data.t,gvtdTimeTrace)
            hold on
            plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
            xlim([PD_data.t(1) PD_data.t(end)]);
            xlabel('Time / s')
            ylabel('GVDT / A.U')
            hold on
            plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
            plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),gvtdTimeTrace(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
            plot(PD_data.t(start_MA_train_locations_frame),thresh,'g*','MarkerSize',8)
            plot(PD_data.t(end_MA_train_locations_frame),thresh,'b*','MarkerSize',8)

        
            PD_data.MA_train_GVDT.start_train_frame = start_MA_train_locations_frame;
            PD_data.MA_train_GVDT.end_train_frame = end_MA_train_locations_frame;

        end
        clear start_MA_train_locations_frame end_MA_train_locations_frame IncCh_train_loc_delta end_MA_train start_MA_train;
        clear MA_train_idx ;
        %%%%%% finding new good ch from threshold of GVDT
         
        avg_g = mean(g(:,PD_data.goodch_idx));
        nStd_meanGVDT =1;
        threshold_meanGVDT = [mean(avg_g) + (nStd_meanGVDT*std(avg_g))];
        figure()
        plot(mean(g(:,PD_data.goodch_idx)),'r+');
        hold on
        yline([mean(avg_g) + (nStd_meanGVDT*std(avg_g))],'b--')
        ylabel('Mean GVDT Value')
        xlabel('Good Ch N')
        title("Mean GVDT good ch, threshold mean +"+nStd_meanGVDT+"x STD ("+threshold_meanGVDT+")")

        %NOW recalcute 'goodch idx' by rejecting channels above threshold
        %in GVDT
        goodch_idx = PD_data.goodch_idx(find(avg_g < threshold_meanGVDT));
        
        for i=1:size(goodch_idx,1)
            if (find(goodch_idx == goodch_idx(i)+64)) | (find(goodch_idx == goodch_idx(i)-64))
                %goodch_idx(i)
                goodch_idx(i) = goodch_idx(i);
                %remCh(i)=1;
            else
                goodch_idx(i) = 0;
            end
        end
    
        goodch_idx = nonzeros(goodch_idx);
        goodch_idx = [goodch_idx; goodch_idx+64];

        %PD_data.goodch_idx = goodch_idx;

        %Perform new GVDT
        nStd = 4;
        gvtdTimeTrace = gvtd(PD_data.dod(:,goodch_idx)'); %input ch x t, output t x 1
        thresh = find_gvtd_thresh(gvtdTimeTrace, statType, nStd);
        gvdt_more_thresh = find(gvtdTimeTrace>thresh);
        N_MA_GVDT = size(gvdt_more_thresh,1);
        
        figure()
        histogram(gvtdTimeTrace)
        hold on
        xline([thresh],LineWidth=4)
        xlabel('GVTD value / A.U')
        ylabel('Frequency / N. of timepoints')
        title("New Histogram GVTD nSTD = "+nStd+" yields threshold of "+thresh+" with "+N_MA_GVDT+" M.As")
        
        tIncCh = ones(size(PD_data.dod,1)  , size(PD_data.dod,2)  );
        tIncCh(gvdt_more_thresh,:) = 0;

        N_MA_HMR = size(find(tIncCh(:,goodch_idx)==0),1)/size(goodch_idx,1);
        figure()
        plot(PD_data.t,PD_data.dod(:,13))
        hold on
        plot(PD_data.t(find(tIncCh(:,13)==0)),PD_data.dod(find(tIncCh(:,13)==0),13),'m.','MarkerSize',6)
        xlabel('Time / s')
        ylabel('dOD / A.U')
        title("Ch 13 Hmr MAD Using GVDT new nStd = "+nStd+" |N M.A = "+N_MA_HMR+" ")

         if size(goodch_idx,1)<12
            n_subfigs_goodchix = size(goodch_idx,1);
         else n_subfigs_goodchix = 12;

         end
        
         figure()
         sgtitle("SNR GVTD new lambda "+nStd+" ")
         for i=1:n_subfigs_goodchix
            subplot(3,4,i)
            plot(PD_data.t, PD_data.dod(:,goodch_idx(i)))
            hold on
            plot(PD_data.t(find(tIncCh(:,goodch_idx(i))==0)),PD_data.dod(find(tIncCh(:,goodch_idx(i))==0),goodch_idx(i)),'m.','MarkerSize',6)
            %xline([900],'b--')
            %xline([5100],'b--')
            title("ch "+num2str(goodch_idx(i)))
            xlabel('T / s')
            ylabel('dod / AU')
          end
        
        g_new = zeros(size(PD_data.dod,1),size(PD_data.dod,2));
        for i=2:size(PD_data.dod,1)
             g_new(i,:)=sqrt( (PD_data.dod(i,:)-PD_data.dod(i-1,:)).^2 );
        end

        figure()
        subplot(2,1,1)
        imagesc(g_new(:,goodch_idx)')
        colorbar, xlabel('Time / samples'), ylabel('Measurement #'),title("GVDT new GoodCh medSNR lambda "+nStd+" ")
        subplot(2,1,2)
        plot(PD_data.t,gvtdTimeTrace)
        hold on
        plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
        xlim([PD_data.t(1) PD_data.t(end)]);
        xlabel('Time / s')
        ylabel('GVDT / A.U')
         
        %figure()
        %plot(PD_data.t,gvtdTimeTrace)
        %hold on
        %plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)

        % find trains/start and end of MA trains, without adding mask.
        mean(gvtdTimeTrace(gvtdTimeTrace < thresh));
        
        tIncCh;
        fs = 10;
        length_train = 200; %frames 
        IncCh_train = ones(1,size(tIncCh,2));
        IncCh_train_loc = ones(size(tIncCh,1),size(tIncCh,2));
        for i=1:size(tIncCh,2)
            if size(unique(tIncCh(:,i)),1) > 1
                for j=1: size(tIncCh,1)-length_train %looking forward
                    if nnz(tIncCh(j:j+(length_train-1),i)) < (length_train*0.9) % IF >20f/200f are MA, train=true
                        IncCh_train(1,i) = 0; %train TRUE (0)
                        IncCh_train_loc(j,i) = 0; %train TRUE (0)
                    end
                end
                for j=1: size(tIncCh,1)-length_train %looking backward
                    if nnz(tIncCh(end-j-(length_train-1):end-j,i)) < (length_train*0.9) % IF >20f/200f are MA, train=true
                        IncCh_train(1,i) = 0; %train TRUE (0)
                        IncCh_train_loc(end-j,i) = 0; %train TRUE (0)
                    end
                end
            end
        end

        IncCh_train = IncCh_train';
        MA_train_idx = IncCh_train;

        for i=1:size(MA_train_idx,1)
            if MA_train_idx(i,1) == 1;
                MA_train_idx(i,1) = i;
            end
        end

        MA_train_idx = nonzeros(MA_train_idx);
        % figure()
        % subplot(2,1,1)
        % plot(PD_data.t,gvtdTimeTrace)
        % hold on
        % plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
        % subplot(2,1,2)
        % plot(PD_data.t,gvtdTimeTrace)
        % hold on
        % plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
        % plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),gvtdTimeTrace(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
         
        %find delta's in data. Consectuvive frames have delta =1. 
        %so a MA train starts/ends when there is a jump of 100
        IncCh_train_loc_delta = find(IncCh_train_loc(:,1)==0);

        %if there are no MA, skip MA analysis
        if size(IncCh_train_loc_delta,1) == 0 
            PD_data.goodch_idx_gvdt = goodch_idx; %added 24 02 25 
        else
        
            %find(diff(find(IncCh_train_loc(:,1)==0))>200);
            %end of MA train
            %if there are many MA trains, else only 1 MA train 
            if size(find(diff(find(IncCh_train_loc(:,1)==0))>100),1)>0
                %end_MA_train(:,1) = find(diff(find(IncCh_train_loc(:,1)==0))>100);
                end_MA_train(:,1) = [find(diff(IncCh_train_loc_delta)>100); size(IncCh_train_loc_delta,1)]; %so it gets end of the last train
            else
                end_MA_train(1,1) = size(IncCh_train_loc_delta,1);
            end
            %end_MA_train(:,1) = find(diff(find(IncCh_train_loc(:,1)==0))>100);
            start_MA_train(:,1) = [1 ; end_MA_train(1:end-1,1)+1];
            %location of start/end of MA train in data (frame N)
            end_MA_train_locations_frame = IncCh_train_loc_delta(end_MA_train,1);
            start_MA_train_locations_frame = IncCh_train_loc_delta(start_MA_train,1);

            PD_data.MA_train_GVDT_new.start_train_frame = start_MA_train_locations_frame;
            PD_data.MA_train_GVDT_new.end_train_frame = end_MA_train_locations_frame;

       
            figure()
            subplot(2,1,1)
            plot(PD_data.t,gvtdTimeTrace)
            hold on
            plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
            subplot(2,1,2)
            plot(PD_data.t,gvtdTimeTrace)
            hold on
            plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
            plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),gvtdTimeTrace(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
            plot(PD_data.t(start_MA_train_locations_frame),thresh,'g*','MarkerSize',8)
            plot(PD_data.t(end_MA_train_locations_frame),thresh,'b*','MarkerSize',8)

        
            figure()
            plot(PD_data.t,PD_data.dod(:,1))
            hold on
            plot(PD_data.t(start_MA_train_locations_frame),PD_data.dod(start_MA_train_locations_frame,1),'g*','MarkerSize',8)
            plot(PD_data.t(end_MA_train_locations_frame),PD_data.dod(end_MA_train_locations_frame,1),'b*','MarkerSize',8)
            %start_MA_train_locations_frame

            figure()
            subplot(2,1,1)
            imagesc(g_new(:,goodch_idx)')
            colorbar, xlabel('Time / samples'), ylabel('Measurement #'),title("GVDT new GoodCh medSNR lambda "+nStd+" ")
            subplot(2,1,2)
            plot(PD_data.t,gvtdTimeTrace)
            hold on
            plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
            xlim([PD_data.t(1) PD_data.t(end)]);
            xlabel('Time / s')
            ylabel('GVDT / A.U')
            hold on
            plot(PD_data.t(find(tIncCh(:,13)==0)),gvtdTimeTrace(find(tIncCh(:,13)==0),1),'m.','MarkerSize',6)
            plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),gvtdTimeTrace(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
            plot(PD_data.t(start_MA_train_locations_frame),thresh,'g*','MarkerSize',8)
            plot(PD_data.t(end_MA_train_locations_frame),thresh,'b*','MarkerSize',8)

            figure()
            sgtitle("SNR GVTD new lambda "+nStd+" ")
            for i=1:n_subfigs_goodchix
                subplot(3,4,i)
                plot(PD_data.t, PD_data.dod(:,goodch_idx(i)))
                hold on
                plot(PD_data.t(find(tIncCh(:,goodch_idx(i))==0)),PD_data.dod(find(tIncCh(:,goodch_idx(i))==0),goodch_idx(i)),'m.','MarkerSize',6)
                plot(PD_data.t(find(IncCh_train_loc(:,13)==0)),PD_data.dod(find(IncCh_train_loc(:,13)==0),1),'y.','MarkerSize',6)
                plot(PD_data.t(start_MA_train_locations_frame),thresh,'g*','MarkerSize',8)
                plot(PD_data.t(end_MA_train_locations_frame),thresh,'b*','MarkerSize',8)
                %xline([900],'b--')
                %xline([5100],'b--')
                title("ch "+num2str(goodch_idx(i)))
                xlabel('T / s')
                ylabel('dod / AU')
            end

            PD_data.goodch_idx_gvdt = goodch_idx;

        end


        if isfield(PD_data,"MA_train_GVDT") %Field exist = TRUE
            %get rid of MA trains with length of 0 or 1 frames
            train_length_GVDT = PD_data.MA_train_GVDT.end_train_frame-PD_data.MA_train_GVDT.start_train_frame; %train length
            %train_length_GVDT_new = PD_data.MA_train_GVDT_new.end_train_frame-PD_data.MA_train_GVDT_new.start_train_frame; %train length
        
            PD_data.MA_train_GVDT.start_train_frame(find(train_length_GVDT<2)) = 0; %set trains of length <2 to be zero
            PD_data.MA_train_GVDT.end_train_frame(find(train_length_GVDT<2)) = 0; %set trains of length <2 to be zero

            PD_data.MA_train_GVDT.start_train_frame = nonzeros( PD_data.MA_train_GVDT.start_train_frame); %remove zero elements
            PD_data.MA_train_GVDT.end_train_frame = nonzeros( PD_data.MA_train_GVDT.end_train_frame); %remove zero elements

            %PD_data.MA_train_GVDT_new.start_train_frame(find(train_length_GVDT_new<2)) = 0; %set trains of length <2 to be zero
            %PD_data.MA_train_GVDT_new.end_train_frame(find(train_length_GVDT_new<2)) = 0; %set trains of length <2 to be zero

            %PD_data.MA_train_GVDT_new.start_train_frame = nonzeros( PD_data.MA_train_GVDT_new.start_train_frame); %remove zero elements
            %PD_data.MA_train_GVDT_new.end_train_frame = nonzeros( PD_data.MA_train_GVDT_new.end_train_frame); %remove zero elements
        end

        if isfield(PD_data,"MA_train_GVDT_new") %Field exist = TRUE

            train_length_GVDT_new = PD_data.MA_train_GVDT_new.end_train_frame-PD_data.MA_train_GVDT_new.start_train_frame; %train length

            PD_data.MA_train_GVDT_new.start_train_frame(find(train_length_GVDT_new<2)) = 0; %set trains of length <2 to be zero
            PD_data.MA_train_GVDT_new.end_train_frame(find(train_length_GVDT_new<2)) = 0; %set trains of length <2 to be zero

            PD_data.MA_train_GVDT_new.start_train_frame = nonzeros( PD_data.MA_train_GVDT_new.start_train_frame); %remove zero elements
            PD_data.MA_train_GVDT_new.end_train_frame = nonzeros( PD_data.MA_train_GVDT_new.end_train_frame); %remove zero elements
        end

        set(gcf,'Visible','on');              
        set(0,'DefaultFigureVisible','on');

end
