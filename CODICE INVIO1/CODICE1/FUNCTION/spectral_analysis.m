function [chRem_def] = spectral_analysis(d,Ts,dRange,SNRrange,th_power_1,th_power_2,f_inf_1,f_sup_1,f_inf_2,f_sup_2,disp_single,disp_result)
% OUTPUT
% chRem_def --> 1 good ch 0 bad ch

% INPUT
% d = Raw Nirs Signal
% Ts = 1/fs
% dRange and SNRrange : RemoveNoisyChannel function parameters
% th_power_1 and th_power_2 : power threshold for each frequency band
% f_inf_1 and f_sup 1 : first frequency band in Hz
% f_inf_2 and f_sup_2 : second frequency band in Hz
% disp_single = 1  --> display each channel both wavelength
% disp_result = 1  --> display final result both wavelength

%% Initialization of support variables
 
n_ch = size(d,2)/2;
t_p = [1:1:size(d,1)]';
t_rapp = Ts*[1:1:size(d,1)]';
%% REMOVE NOISY CHANNEL

% dRange = [5e-4 3];
% SNRrange = 1;
% 1 = good channel ; 0 = bad channel 
removeCh= removeNoisyChannel(d,dRange,SNRrange);

%% SPECTRAL ANALYSIS
% FFT

% BOTH wavelength
mean_d = mean(d);
    
for i = 1:size(d,2)
    d_spect(:,i) = fft((d(:,i)-mean_d(i)));
end 

% POWER
fs = 1/Ts;
n = size(d,1);                % number of samples
f = (0:n-1)*(fs/n);           % frequency range
power = abs(d_spect).^2/n;    % power of the DFT

%% FIRST WAVELENGTH

f_r = round(f,4);

idxStart = find(f_r >= f_inf_1,1,'first');
idxEnd   = find(f_r >= f_sup_1,1,'first');

idxStart2 = find(f_r >= f_inf_2,1,'first');
idxEnd2   = find(f_r >= f_sup_2,1,'first');

chRemSpec_1 = ones(n_ch,1);

% 0 REMOVED
% 1 KEPT

for iCh = 1:n_ch
    
    % Find all the frequencies for each channel under analysis which have
    % the power above th_power in the first frequency band or in the second
    % frequency band. If at least one frequency is above the threshold the
    % channel is flagged as to remove
    tmp  = find(power(idxStart:idxEnd,iCh)>th_power_1);
    tmp2 = find(power(idxStart2:idxEnd2,iCh)>th_power_2);
    
    % If tmp is not empty, the channel must be removed
    if ~isempty(tmp) | ~isempty(tmp2) 
        chRemSpec_1(iCh) = 0;
    end
end

%% DISPLAY RESULTS FIRST WAVELENGTH

% n_ch are the channels for the single wavelength
% if disp_single == 1
%     a = n_ch;                                                                  
%     for iCh = 1:n_ch
%         figure()
%         subplot(2,1,1)
%         plot(t_p,d_clean_first(:,iCh))
%         title(['RAW NIRS DATA ONLY GOOD CHANNEL FIRST WAVELENGTH - ',num2str(iCh)])
%         axis tight
%         subplot(2,1,2)
%         plot(f(1:n/2),power(1:n/2,iCh))
%         title('Power Spectra in [0 - 1] Hz')
%         hold on
%         xline(0.3732,'r--')
%         pause;
%         close;
%     end
% end

% if disp_single == 1
%     a = n_ch;                                                               %%%%%%%%%%%%%%%%%%%%%
%     for iCh = 1:a
%         figure()
%         subplot(2,1,1)
%         plot(t_p,d(:,iCh))
%         title('RAW NIRS DATA ONLY GOOD CHANNEL FIRST WALENGTH')
%         axis tight
%         subplot(2,1,2)
%         plot(f(1:n/2),power(1:n/2,iCh))
%         title('Power Spectra in [0 - 1] Hz')
%         hold on
%         xline(0.05,'r--')
%         xline(2,'r--')
%         if chRemSpec_1(iCh) == 0
%             title('Removed')
%         else
%             title('Kept')
%         end
%         pause;
%         close;
%     end
% end

% if disp_result == 1
%     d_clean_first_w = d(:,1:n_ch)
%     figure;
%     plot(t_p,d_clean_first_w)
%     axis tight
%     title('ALL CHANNEL FIRST WAVELENGTH')
% 
%     figure()
%     plot(t_p,d_clean_first_w(:,chRemSpec_1==1))
%     axis tight
%     title('CHANNEL KEPT FIRST WAVELENGTH')
% 
%     figure()
%     plot(t_p,d_clean_first_w(:,chRemSpec_1==0))
%     axis tight
%     title('CHANNEL REMOVED FIRST WAVELENGTH')
% end

%% SECOND WAVELENGTH


idxStart = find(f_r >= f_inf_1,1,'first');
idxEnd   = find(f_r >= f_sup_1,1,'first');

idxStart2 = find(f_r >= f_inf_2,1,'first');
idxEnd2   = find(f_r >= f_sup_2,1,'first');
chRemSpec_2 = ones(n_ch,1);

% 0 REMOVED
% 1 KEPT

pos = 0;

for iCh = n_ch+1:size(d,2)
    
    pos = pos + 1;

    % Find all the frequencies for each channel under analysis which have
    % the power above th_power in the first frequency band or in the second
    % frequency band. If at least one frequency is above the threshold the
    % channel is flagged as to remove

    tmp  = find(power(idxStart:idxEnd,iCh)>th_power_1);
    tmp2 = find(power(idxStart2:idxEnd2,iCh)>th_power_2);
    
    % If tmp is not empty, the channel must be removed
    if ~isempty(tmp) | ~isempty(tmp2)
        chRemSpec_2(pos) = 0;
    end
end

%% DISPLAY RESULTS SECOND WAVELENGTH

% if disp_single == 2
%     aa = size(d,2);                                                 
%     for iCh = n_ch+1:aa
%         figure()
%         subplot(2,1,1)
%         plot(t_p,d(:,iCh))
%         title(['RAW NIRS DATA ONLY GOOD CHANNEL SECOND WAVELENGTH - ',num2str(iCh)])
%         axis tight
%         subplot(2,1,2)
%         plot(f(1:n/2),power(1:n/2,iCh))
%         title('Power Spectra in [0 - 1] Hz')
%         hold on
%         xline(0.3732,'r--')
%         pause;
%         close;
%     end
% end

% if disp_single == 2
%     aa = size(d,2)                                                 
%     pos = 0;
%     for iCh = n_ch+1:aa
%         pos = pos + 1;
%         figure()
%         subplot(2,1,1)
%         plot(t_p,d(:,iCh))
%         title('RAW NIRS DATA ONLY GOOD CHANNEL SECOND WALENGTH')
%         axis tight
%         subplot(2,1,2)
%         plot(f(1:n/2),power(1:n/2,iCh))
%         hold on
%         xline(0.05,'r--')
%         xline(2,'r--')
%         if chRemSpec_2(pos) == 0
%             title('Removed')
%         else
%             title('Kept')
%         end
%         pause;
%         close;
%     end
% end

% if disp_result == 2
%     d_clean_second_w = d(:,n_ch+1:end);
% 
%     figure;
%     plot(t_p,d_clean_second_w)
%     title('ALL CHANNEL SECOND WAVELENGTH')
%     axis tight
% 
%     figure;
%     plot(t_p,d_clean_second_w(:,chRemSpec_2==1))
%     title('CHANNEL KEPT SECOND WAVELENGTH')
%     axis tight
% 
%     figure;
%     plot(t_p,d_clean_second_w(:,chRemSpec_2==0))
%     title('CHANNEL REMOVED SECOND WAVELENGTH')
%     axis tight
% 
% end

%% REMOVE CHANNEL

% If a channel is removed from the first wavelength, it must be removed
% also for the second one (and viceversa)
chRem_def = ones(size(d,2),1);

for i = 1:1:n_ch
    if chRemSpec_1(i) == 0 | chRemSpec_2(i) == 0 | removeCh(i) == 0
        chRem_def(i) = 0;
        chRem_def(i+n_ch) = 0;
    end
end

%%
if disp_single == 1

    for iCh = 1:1:size(d,2)
        figure()
        set(gcf, 'Position', get(0, 'Screensize'));

        subplot(2,1,1)
        plot(t_rapp,d(:,iCh))
        xlabel('Time [s]')
        ylabel('Amplitude [A.U.]')
        if iCh<=n_ch
            title(['RAW NIRS DATA FIRST WALENGTH : [',num2str(f_inf_1),'-',num2str(f_sup_1),'] - [',num2str(f_inf_2),'-',num2str(f_sup_2),'] Hz , Power th: ',num2str(th_power_1),' ',num2str(th_power_2)])
            title(['RAW NIRS DATA \lambda = 780 nm - Channel ',num2str(iCh)])
        else
            title(['RAW NIRS DATA SECOND WALENGTH : [',num2str(f_inf_1),'-',num2str(f_sup_1),'] - [',num2str(f_inf_2),'-',num2str(f_sup_2),'] Hz , Power th: ',num2str(th_power_1),' ',num2str(th_power_2)])
        end
        axis tight
        subplot(2,1,2)
        plot(f(1:n/2),power(1:n/2,iCh))
        xlabel('Frequency [Hz]')
        ylabel('PSD [A^2/Hz]')
        hold on
        xline(f_inf_1,'r--')
        xline(f_sup_1,'r--')
        xline(f_inf_2,'r--')
        xline(f_sup_2,'r--')
        if chRem_def(iCh) == 0
            title(['SPECTRAL ANALYSIS, Channel ',num2str(iCh),' - REMOVED'])
        else
            title(['SPECTRAL ANALYSIS, Channel ',num2str(iCh),' - KEPT'])
        end
        if iCh == 13
            print('FIG_28_SPETT_CH13_PD8','-djpeg','-r600')
        end
        pause;
        close;
    end
end

if disp_result == 1
    
    d_first  = d(:,1:n_ch) ;
    d_second = d(:,n_ch+1:end);

    figure;
    plot(t_p,d_first)
    title('ALL CHANNEL FIRST WAVELENGTH')
    axis tight

    figure;
    plot(t_p,d_second)
    title('ALL CHANNEL SECOND WAVELENGTH')
    axis tight

    if sum(removeCh==0)>0
        figure;
        plot(t_p,d_first(:,removeCh(1:n_ch)==0))
        title('CHANNEL REMOVED SNR AND AMP FIRST WAVELENGTH')
        axis tight
    else 
        disp('NO CHANNELS REMOVED FROM PRUNING AT THE FIRST WAVELENGTH')
    end
    
    if sum(removeCh==0)>0
        figure;
        plot(t_p,d_second(:,removeCh(n_ch+1:end)==0))
        title('CHANNEL REMOVED SNR AND AMP SECOND WAVELENGTH')
        axis tight
    else
        disp('NO CHANNELS REMOVED FROM PRUNING AT THE FIRST WAVELENGTH')
    end
    
    if sum(chRemSpec_1==0)>0
        figure;
        plot(t_p,d_first(:,chRemSpec_1==0))
        title('CHANNEL REMOVED SPECTRAL ANALYSIS FIRST WAVELEGTH')
        axis tight
    else
        disp('NO CHANNELS REMOVED FROM SPECTRAL ANALYSIS AT THE FIRST WAVELENGTH')
    end
    
    if sum(chRemSpec_2==0)>0
        figure;
        plot(t_p,d_second(:,chRemSpec_2==0))
        title('CHANNEL REMOVED SPECTRAL ANALYSIS SECOND WAVELEGTH')
        axis tight 
    else
        disp('NO CHANNELS REMOVED FROM SPECTRAL ANALYSIS IN THE SECOND WAVELENGTH')
    end

    figure;
    plot(t_p,d_first(:,chRem_def(1:n_ch)==1))
    title('GOOD CHANNEL FIRST WAVELENGTH')
    axis tight
    
    figure;
    plot(t_p,d_second(:,chRem_def(n_ch+1:end)==1))
    title('GOOD CHANNEL SECOND WAVELENGTH')
    axis tight
end

end  % end function