%% Assignment 8
%% 1.a) Generate a spike sereies that is 1024 points long with a spike at 512.
% filter the spike and then take the fft to see the amplitude and phase
% spectrum 
sp_sr = zeros(1,1024); %creating a serie of series
sp_sr(512)= 1; %creating a spike at n = 512 
NC = length(sp_sr); %length of spike series
dt = 1; %s 
t= (0:(NC-1)*dt); % creating a time array 
fc = 0.05 ; %corner frequency 
N = [1,2,4,8]; % Different poles 
fn = fftshift([0:NC/2,-NC/2+1:-1].*(1/dt)); %shifted 

%For each pole I am filtering the spike series, taking the fft, splitting
%into the amplitude and phase spectrum and then take the ifft. Each row in
%the sp_filt array has a row for each different pole. 
for n = 1:length(N)  
sp_filt(n,:) = lpfilt(sp_sr,dt,fc,N(n),'TYPE =1'); %filtering the spike series 
sp_fft(n,:) = fftshift(fft(sp_filt(n,:))); %taking the fft then shifting the spike series
sp_amp(n,:) = abs(sp_fft(n,:)); %geting the ampltiude spectrum of the fft of the spike series 
sp_phy(n,:) = atan2(imag(sp_fft(n,:)),real(sp_fft(n,:))); %getting the phase spectrum  
sp_inv(n,:) = ifft(fftshift(sp_fft(n,:))); %taking the ifft of the filtered spectrum 
end

%% Plotting 
% plotting the original time series and the time series filtered by 4
% different poles 
figure(3) 
ax1 = subplot(2,1,1)
plot(ax1,t,sp_sr(1,:),'r');
hold on; 
xlabel(ax1,'Time(s)'); 
ylabel(ax1,'Signal Strength'); 
title(ax1,'Original Spike Series and Filtered Spike Series for 4 Different Poles') 
subplot(2,1,2)
plot(t,sp_inv(1,:),'r');
hold on;
plot(t,sp_inv(2,:),'b');
hold on;
plot(t,sp_inv(3,:),'g');
hold on; 
plot(t,sp_inv(4,:),'k');
hold on; 
legend('N=1','N=2','N=4', 'N=8')
xlabel('Time(s)'); 
ylabel('Signal Strength'); 
xlim([450, 700])

figure(1) 
plot(fn,sp_amp(1,:),'r');
hold on; 
plot(fn,sp_amp(2,:),'b');
hold on;
plot(fn,sp_amp(3,:),'g');
hold on; 
plot(fn,sp_amp(4,:),'k');
hold on; 
xlabel('Frequency (Hz)'); 
ylabel('Amplitude'); 
legend('N=1','N=2','N=4', 'N=8')
title('Amplitude Diagram of Spike Series With a Low Pass Filter at fc=0.05Hz')

figure (2) 
plot(fn,sp_phy(4,:),'k');
hold on;
plot(fn,sp_phy(1,:),'r');
hold on; 
plot(fn,sp_phy(2,:),'b');
hold on;
plot(fn,sp_phy(3,:),'g');
hold on; 
xlabel('Frequency (Hz)'); 
ylabel('Phase (Radians)'); 
legend('N=8', 'N=1','N=2','N=4')
title('Phase Diagram of Spike Series With a Low Pass Filter at fc=0.05Hz')


%%  Comments 
%  Firstly, when looking at the amplitude spectrum of the spike series with a 
% low pass filter at fc = 0.05 Hz (Figure 1), we see that the steepness 
% increases with the number of poles. Somethine else interesting is that we
% see that starting point decrease closer to zero with increasing poles,
% where when N=1 (red curve) the amplitude spectrum starts around 0.1. 
% Looking at the phase spectrum for the spike series filtered through a 
% minimum phase low pass filter (Figure 2), we see that the phase spectrum 
% is the same for each number of poles. When plotted on top of each other,
% we can only see one of the curves at a time since they're the same. 
% In Figure 3, when we inverse the filtered spike series back into the time
% domain, instead of a large spipke at t = 512, we see a larger width of 
% the spike as N increases. Also with increasing N, we see more side lobes, 
% with N=8 having 2 side lobes on each side of the central spike. I 
% couldn't seem to normalize the filtered series to have the amplitude 
% reach 1, like in the spike series, especially since as you increase N,
% you descrease the overall height of the spike. 

%% c. now vary the cut-off frequency FC (pick 7 values, 3 above and 3 below 0.05, and 0.05)

fc_v = [0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2]; 
pole = 8;
colour = ['r','b','c','g', 'y', 'm', 'k'];


for l = 1:length(fc_v)
    spv_sr(l,:) = lpfilt(sp_sr,dt,fc_v(l),pole,'TYPE=1');
    spv_fft(l,:) = fftshift(fft(spv_sr(l,:)));
    spv_amp(l,:) = abs(spv_fft(l,:)); 
    spv_phy(l,:) = atan2(imag(spv_fft(l,:)),real(spv_fft(l,:))); 
    spv_ift(l,:) = ifft(fftshift(spv_fft(l,:))); 
    
    figure(4) 
    plot(t,spv_ift(l,:), colour(l))
    hold on; 
    xlabel('Time(s)'); 
    xlim([0,1024]);
    ylabel('Signal Strength'); 
    title('Filtered Spike Time Series with Various Cut-Off Frequency'); 
    legend(strcat('fc =',num2str(fc_v')))
    
    figure(5) 
    plot(fn, spv_amp(l,:), colour(l))
    hold on; 
    xlabel('Frequency (Hz)'); 
    ylabel('Signal Strength');
    title('FFT of Filtered Spike Time Series with Various Cut-Off Frequency'); 
    legend(strcat('fc =',num2str(fc_v')))   
end 
% legend(fg5, strcat('fc =',num2str(fc_v')))
% legend(fg4, strcat('fc =',num2str(fc_v')))

%% Comments
% Looking at the filtered spike time series with various cut-off frequencies 
% from fc =0.005 Hz to fc= 0.2 Hz (Figure 4), we see a trend in steepness 
% of the spike time series with increasing fc values. When fc = 0.005 Hz, 
% the spike series becomes bround around the spike value of t = 512 with a 
% signal strength of less than 0.1, decreasing the signal strength by over 
% a factor of 10. When looking at fc = 0.2 Hz, there is an increased 
% steepness with a large central lobe with an amplitude of approximately
% 0.45. With increasing fc, we also have increased number of side lobes. 
% Looking at the FFT of the filtered time series spike (Figure 5), we see 
% a trend of increasing width of the central lobes with increasing fc.
% Meaning with a higher cut-off frequency, the signal becomes made up 
% of more frequencies to make up the signal. 


%% d. 
ND = [4,8]; 

for d = 1:length(ND)
    sp_min(d,:) = lpfilt(sp_sr,dt,fc,ND(d),'TYPE=1');
    sp_zero(d,:) = lpfilt(sp_sr,dt,fc,ND(d),'TYPE=0');
    spm_fft(d,:) = fftshift(fft(sp_min(d,:)));
    spm_amp(d,:) = abs(sp_min(d,:)); 
    spz_fft(d,:) = fftshift(fft(sp_zero(d,:)));
    spz_amp(d,:) = abs(sp_zero(d,:)); 
    spm_ift(d,:) = ifft(fftshift(spm_fft(d,:))); 
    spz_ift(d,:) = ifft(fftshift(spz_fft(d,:)));
end 

figure(6) 
plot(t,spm_ift(1,:),'r.-'); hold on; 
plot(t,spm_ift(2,:), 'b.-'); hold on; 
plot(t, spz_ift(1,:), 'g-'); hold on; 
plot(t, spz_ift(2,:), 'c-'); hold on; 
xlabel('Time(s)'); 
xlim([400,600])
ylabel('Signal Strength'); 
title('Filtered Spike Time Series with Zero Phase or Minimum Phase Filter');
legend('Minimum Phase N=4','Minimum Phase N=8','Zero Phase N=4','Zero Phase N=8');

figure(7) 
plot(fn,spm_amp(1,:),'r.-'); hold on; 
plot(fn,spm_amp(2,:),'b.-'); hold on; 
plot(fn,spz_amp(1,:),'g-'); hold on; 
plot(fn,spz_amp(2,:),'c-'); hold on; 
xlabel('Frequency (Hz)'); 
ylabel('Signal Strength');
xlim([-200,200])
title('FFT of Filtered Spike Time Series with Zero Phase or Minimum Phase Filter');
legend('Minimum Phase N=4','Minimum Phase N=8','Zero Phase N=4','Zero Phase N=8');

%% Comments 
% When looking at the filtered spike time series for a minimum or zero 
% phase filter (Figure 6), we see no difference between the minimum or zero
% phase filter for a certain value of N. when looking at N=4, we see at the
% minimum (red curve with dots) and the zero phase curve (green) are 
% overlaid on each other. In the frequency domain (Figure 7), we see a 
% similar trend where the zero phase and the minimum phase filter produce 
% the same curve. Looking at N=8 in the frequency domain, we see that the 
% minimum phase (blue with dots) and the zero phase (cyan) are the same. 

