function [MA_train_idx tIncCh IncCh_train IncCh_train_loc] = find_MA_train(PD_data,dodConv,goodch_idx)

fs = 10;
length_train = 600; %frames 

% Detect motion artifacts in signal
tMotion = 0.5;%1.0%0.5; %0.8; % %0.5; %time range in seconds
tMask = 2; %mark data *- time around m.a as m.a
SDThresh = 10%6.5 %10; %12
AmpThresh = 0.2 %0.35; %0.5; 
tIncMan = ones(length(PD_data.t ),1); % set it to vectors of ones (this is a vector used to remove manually parts of the data if needed)
% Motion detection technique. tIncCh is a matrix number of samples x twice n of
% channels which contains for each channel (column) the information about
% whether an artifact was present (0s) or not (1s). tInc is a vector which
% contains information on whether at that time sample in any of the channel
% was present an artifact (0s) or not (1s). tInc can therefore be obtained
% from tIncCh by setting to 0 every row that contains at least one 0. 

%spline wavelet
[tInc,tIncCh] = hmrMotionArtifactByChannel(dodConv, fs, PD_data.SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);
remCh = PD_data.SD.MeasListAct;

%find M.A train
IncCh_train = ones(1,size(tIncCh,2));
IncCh_train_loc = ones(size(tIncCh,1),size(tIncCh,2));
for i=1:size(tIncCh,2)
    if size(unique(tIncCh(:,i)),1) > 1
        for j=1: size(tIncCh,1)-length_train
            if nnz(tIncCh(j:j+(length_train-1),i)) < (length_train/2)
                IncCh_train(1,i) = 0;
                IncCh_train_loc(j,i) = 0;
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
%1 = No train
%0 = MA Train present in channel


% goodch_nontrain = zeros(size(IncCh_train,1),1);
% 
% for i=1:size(IncCh_train,1);
%     if sum([IncCh_train(i,1) remCh(i,1)]) == 2;
%         goodch_nontrain(i,1) = 1; 
%     end
% end
% %nnz(goodch_nontrain);
% 
% goodch_nontrain_idx = zeros(size(IncCh_train,1),1);
% for i=1:size(IncCh_train,1)
%     if goodch_nontrain(i,1) == 1;
%         goodch_nontrain_idx(i,1) = i;
%     end
% end
% goodch_nontrain_idx = nonzeros(goodch_nontrain_idx);




end