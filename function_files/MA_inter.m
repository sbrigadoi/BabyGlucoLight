function PD_data = MA_inter(PD_data)
PD_data.inter_dod_new = PD_data.dod;
%interdod_new = PD_data.dod;

for i=1:size(PD_data.dod,2)
    if size(unique(PD_data.IncCh_train_loc(:,i)),1) ==2 %if channel has MA train
        %get location of MA in data (frame N)
        MA_locations_true = [find(PD_data.tIncCh(:,i)==0)];
        MA_train_locations_true= find(PD_data.IncCh_train_loc(:,i)==0);

        %find positions of end of MA train
        end_MA_train_locations_true = find(  (diff(MA_train_locations_true)==1) ==0);
        end_MA_train_locations_true = [end_MA_train_locations_true; size(MA_train_locations_true,1)];

        %find positions of start of MA train
        start_MA_train_locations_true = 1;
        start_MA_train_locations_true = [ 1; end_MA_train_locations_true(1:end-1)+1];

        %location of start/end of MA train in data (frame N)
        end_MA_train_locations_true = MA_train_locations_true(end_MA_train_locations_true);
        start_MA_train_locations_true = MA_train_locations_true(start_MA_train_locations_true);

        interdod = zeros(size(PD_data.dod,1), size(end_MA_train_locations_true,1));
        %do interpolation
            for j=1:size(end_MA_train_locations_true,1)
        
                sample_points = ([start_MA_train_locations_true(j) ; end_MA_train_locations_true(j)]);
                sample_values = (PD_data.dod(sample_points,1));
                q_points = [start_MA_train_locations_true(j):1:end_MA_train_locations_true(j) ]';

                interdod(1:size(q_points,1) ,j) = interpn(sample_points,sample_values,q_points,'linear');
            end
        PD_data.inter_dod = interdod;
        interdod_new = PD_data.dod(:,i);
            for k=1:size(interdod,2)
                %interdod_new(start_MA_train_locations_true(k):end_MA_train_locations_true(k)) = interdod(start_MA_train_locations_true(k):end_MA_train_locations_true(k),k);
                interdod_new(start_MA_train_locations_true(k):end_MA_train_locations_true(k)) = interdod(1: size( start_MA_train_locations_true(k):end_MA_train_locations_true(k) ,2)  ,k);

            end
        PD_data.inter_dod_new(:,i) = interdod_new;
    end

end
figure()
plot(PD_data.t(start_MA_train_locations_true(1):end_MA_train_locations_true(1)),PD_data.dod(start_MA_train_locations_true(1):end_MA_train_locations_true(1),1))
hold on
plot(PD_data.t(start_MA_train_locations_true(1):end_MA_train_locations_true(1)),PD_data.inter_dod_new(start_MA_train_locations_true(1):end_MA_train_locations_true(1),1))
xlabel('Time / s')
ylabel('dOD / A.U')

figure()
plot(PD_data.t,PD_data.dod(:,1))
hold on
plot(PD_data.t,PD_data.inter_dod_new(:,1))
xlabel('Time / s')
ylabel('dOD / A.U')


end