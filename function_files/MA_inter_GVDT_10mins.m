function PD_data = MA_inter_GVDT_10mins(PD_data,figON)
        
        rng('default');

        if figON == 0
            set(gcf,'Visible','off');              
            set(0,'DefaultFigureVisible','off');

        elseif figON == 1
            set(gcf,'Visible','on');              
            set(0,'DefaultFigureVisible','on');
        end
        %sample time points
        sample_length_f = 600;
        %do interpolation for each MA train and each channel
        dod_int = PD_data.dod;

if ~isfield(PD_data,"MA_train_GVDT_new") %TRUE = FIELD DOESNT EXIST
            dod_int = PD_data.dod;% NO INT SINCE THERE ARE NO MA TRAINS
else
    for j=1:size(PD_data.MA_train_GVDT_new.end_train_frame,1) %N of MA trains
        
        int_dod_noise= [];

            %skip_MA_train = 0; %0 = dont skip train, 1=skip this MA train
            %check length of MA train 
            if PD_data.MA_train_GVDT_new.end_train_frame(j)-PD_data.MA_train_GVDT_new.start_train_frame(j) < 21
                sample_gap = PD_data.MA_train_GVDT_new.end_train_frame(j)-PD_data.MA_train_GVDT_new.start_train_frame(j);
                if sample_gap > 1
                    sample_gap = sample_gap -1; %reduce size of gap by one so no overlapping start:end of MA train sample pts
                end

                %if sample_gap < 2 %sample gap is 0 or 1, so skip this MA train!
                %    skip_MA_train = 1;
                %end
            else 
                sample_gap = 20; %default length of sample gap if the length of MA train is > 20frames
            end

        %if the start/end of the MA train is within sample_length_f to
        %the start/end of data, we need to reduce sample_length_f
        sample_length_f = 600;

        %if the start/end of the MA train are within 20f of frame 1 or
        %frame end, then just use the start/end for the sample

        if PD_data.MA_train_GVDT_new.start_train_frame(j) < 21 
        % start train is near start - so only use end %%%

            if PD_data.MA_train_GVDT_new.end_train_frame(j)+sample_length_f> size(PD_data.t,1)
                sample_length_f = sample_length_f - (PD_data.MA_train_GVDT_new.end_train_frame(j)+sample_length_f -size(PD_data.t,1) )  ; 
            end
            for i= 1:size(PD_data.dod,2) %channels
        
                sample_t = PD_data.t(PD_data.MA_train_GVDT_new.end_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j)+sample_length_f)';
                sample_t = [PD_data.t(PD_data.MA_train_GVDT_new.start_train_frame(j):PD_data.MA_train_GVDT_new.start_train_frame(j)+sample_gap )' sample_t]; %new - just using last 20frames of end   
                %sample gap was 20

                sample_dod = PD_data.dod(PD_data.MA_train_GVDT_new.end_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j)+sample_length_f,i)';
                sample_dod = [PD_data.dod(PD_data.MA_train_GVDT_new.start_train_frame(j):PD_data.MA_train_GVDT_new.start_train_frame(j)+sample_gap,i)' sample_dod ];%new - just using last 20frames of end
                %sample gap was 20

                %query time points (the MA train)
                query_t = PD_data.t(PD_data.MA_train_GVDT_new.start_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j),1)';
                %Interpolating across MA train
                int_dod = interpn(sample_t,sample_dod,query_t,'pchip');
                %trying to fix sample points must be sorted in ascending order
                %12 07 24 23:08 BST
                %int_dod = interpn(unique(sample_t)',unique(sample_dod)',unique(query_t)','pchip');

                %add random noise, using mean of STD 20f before and after MA train
                rng(0); %Fix seed of RNG, so it is always the same
                %just using start of sample, since it is a end sample
                std_int = mean( [std(sample_dod(1:20))] )*0.8; %0.8 of STD
                int_dod_noise(:,i) = int_dod + randn(1,size(int_dod,2))*std_int;

                dod_int(PD_data.MA_train_GVDT_new.start_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j),i) = int_dod_noise(:,i);
            end
        

        elseif PD_data.MA_train_GVDT_new.end_train_frame(j) > (size(PD_data.t,1)-21)
        % end train is near end - so only use start train
             if PD_data.MA_train_GVDT_new.start_train_frame(j)-sample_length_f< 0
                sample_length_f = sample_length_f + PD_data.MA_train_GVDT_new.start_train_frame(j)-sample_length_f  -1  ; %-1 avoids 0 values 
            end
            for i= 1:size(PD_data.dod,2) %channels

            sample_t = PD_data.t(PD_data.MA_train_GVDT_new.start_train_frame(j)-sample_length_f:PD_data.MA_train_GVDT_new.start_train_frame(j))';
            sample_t = [sample_t PD_data.t(PD_data.MA_train_GVDT_new.end_train_frame(j)-sample_gap:PD_data.MA_train_GVDT_new.end_train_frame(j) )' ]; %new - just using last 20frames of end   
            %                                                  this 20 causes overlap if the gap between start and end is less than 20frames
            %sample gap was 20

            sample_dod = PD_data.dod(PD_data.MA_train_GVDT_new.start_train_frame(j)-sample_length_f:PD_data.MA_train_GVDT_new.start_train_frame(j),i)';
            sample_dod = [sample_dod PD_data.dod(PD_data.MA_train_GVDT_new.end_train_frame(j)-sample_gap:PD_data.MA_train_GVDT_new.end_train_frame(j),i)' ];%new - just using last 20frames of end
            %sample gap was 20

            %query time points (the MA train)
            query_t = PD_data.t(PD_data.MA_train_GVDT_new.start_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j),1)';
            %Interpolating across MA train

            int_dod = interpn(sample_t,sample_dod,query_t,'pchip');
            %trying to fix sample points must be sorted in ascending order
            %12 07 24 23:08 BST
            %int_dod = interpn(unique(sample_t)',unique(sample_dod)',unique(query_t)','pchip');


            
            % hold on
            

            %add random noise, using mean of STD 20f before and after MA train
            rng(0); %Fix seed of RNG, so it is always the same
            %just using end of sample, since it is a start sample
            std_int = mean( [std(sample_dod(end-20:end)) ])*0.8; %0.8 of STD
            int_dod_noise(:,i) = int_dod + randn(1,size(int_dod,2))*std_int;
            dod_int(PD_data.MA_train_GVDT_new.start_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j),i) = int_dod_noise(:,i);

            %figure()
            %plot(dod_int(:,1))
            end

        else
            %if the start/end of the MA train is within sample_length_f to
            %the start/end of data, we need to reduce sample_length_f
            %(600f)
            if PD_data.MA_train_GVDT_new.start_train_frame(j)-sample_length_f< 0
                sample_length_f = sample_length_f + PD_data.MA_train_GVDT_new.start_train_frame(j)-sample_length_f  -1  ; %-1 avoids 0 values 
            end

            if PD_data.MA_train_GVDT_new.end_train_frame(j)+sample_length_f> size(PD_data.t,1)
                sample_length_f = sample_length_f - (PD_data.MA_train_GVDT_new.end_train_frame(j)+sample_length_f -size(PD_data.t,1) )  ; 
            end

            for i= 1:size(PD_data.dod,2) %channels
        
                sample_t = PD_data.t(PD_data.MA_train_GVDT_new.start_train_frame(j)-sample_length_f:PD_data.MA_train_GVDT_new.start_train_frame(j))';
                sample_t = [sample_t PD_data.t(PD_data.MA_train_GVDT_new.end_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j)+sample_length_f)' ];
        
                %sample dod
                sample_dod = PD_data.dod(PD_data.MA_train_GVDT_new.start_train_frame(j)-sample_length_f:PD_data.MA_train_GVDT_new.start_train_frame(j),i)';
                sample_dod = [sample_dod PD_data.dod(PD_data.MA_train_GVDT_new.end_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j)+sample_length_f,i)' ];
                
                %need to take into account MA trains that are small in
                %length, or even same start/end time pt
                
                %query time points (the MA train)
                query_t = PD_data.t(PD_data.MA_train_GVDT_new.start_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j),1)';
                %Interpolating across MA train
                int_dod = interpn(sample_t,sample_dod,query_t,'pchip');
                %trying to fix sample points must be sorted in ascending order
                %12 07 24 23:08 BST
                %int_dod = interpn(unique(sample_t)',unique(sample_dod)',unique(query_t)','pchip');
                
                %add random noise, using mean of STD 20f before and after MA train
                rng(0); %Fix seed of RNG, so it is always the same
            
                %if the sample DOD is smaller than the X frames we want to
                %sample for the STD calculation
                %if size(sample_dod,2)<20
                %    std_int = mean( [std(sample_dod((end/2)- size(sample_dod,2)*0.25  :end/2)) std(sample_dod((end/2):(end/2)+20))])*0.8; %0.8 of STD
                %end
        
                std_int = mean( [std(sample_dod((end/2)-20:end/2)) std(sample_dod((end/2):(end/2)+20))])*0.8; %0.8 of STD
                int_dod_noise(:,i) = int_dod + randn(1,size(int_dod,2))*std_int;

                dod_int(PD_data.MA_train_GVDT_new.start_train_frame(j):PD_data.MA_train_GVDT_new.end_train_frame(j),i) = int_dod_noise(:,i);

                %int_dod_noise_joined
            end


        


        end

    end %for N trains
 end %if exists trains
    PD_data.dod_int = dod_int;

        %figure()
        %plot(sample_t,sample_dod)
        %hold on
        %plot(query_t,int_dod_noise)
        % 
         figure()
         plot(PD_data.t,PD_data.dod(:, PD_data.goodch_idx(1) ))
         hold on
         plot(PD_data.t,PD_data.dod_int(:,PD_data.goodch_idx(1 )) )
         title("Dod int goodch idx - PD "+num2str(PD_data.subjectN)+" event "+PD_data.eventType+" N "+num2str(PD_data.eventN)+" goodch "+num2str(PD_data.goodch_idx(1))+" ")
         legend("Dod","Dod Int")
         xlabel("Time / s")
         ylabel(" \Delta OD / AU")
% 
        % figure()
        % plot(sample_t,sample_dod)
        % hold on
        % plot(query_t,int_dod_noise)


%%
        % interdod = zeros(size(PD_data.dod,1), size(end_MA_train_locations_true,1));
        % %do interpolation
        %     for j=1:size(end_MA_train_locations_true,1)
        % 
        %         sample_points = ([start_MA_train_locations_true(j) ; end_MA_train_locations_true(j)]);
        %         sample_values = (PD_data.dod(sample_points,1));
        %         q_points = [start_MA_train_locations_true(j):1:end_MA_train_locations_true(j) ]';
        % 
        %         interdod(1:size(q_points,1) ,j) = interpn(sample_points,sample_values,q_points,'linear');
        %     end
        % PD_data.inter_dod = interdod;
        % interdod_new = PD_data.dod(:,i);
        %     for k=1:size(interdod,2)
        %         %interdod_new(start_MA_train_locations_true(k):end_MA_train_locations_true(k)) = interdod(start_MA_train_locations_true(k):end_MA_train_locations_true(k),k);
        %         interdod_new(start_MA_train_locations_true(k):end_MA_train_locations_true(k)) = interdod(1: size( start_MA_train_locations_true(k):end_MA_train_locations_true(k) ,2)  ,k);
        % 
        %     end
        % PD_data.inter_dod_new(:,i) = interdod_new;


        set(gcf,'Visible','on');              
        set(0,'DefaultFigureVisible','on');

end