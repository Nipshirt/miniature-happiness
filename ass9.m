%% Assignment 9 

%% Compute Autocorrelation of Sweep by taking the ifft of the Power Spectrum - Zero Pad to 2N and plot the autocorrelation 
n_s = length(s);
time = [1:n_s].*dt;
freq = [0:n_s/2,-n_s/2+1:-1].*(1/dt); 

s_fft = fftshift(fft(s,2*n_s));
s_ps = s_fft.^2;

s_ifft = ifft(fftshift(s_ps),n_s)/n_s; 

figure(1) 
subplot(2,1,1)
plot(time,s)
xlabel('Time (s)') 
title('Time Series S')
ylabel('Amplitude')

subplot(2,1,2) 
plot(time, s_ifft)
xlabel('Time (s)') 
ylabel('Amplitude')
title('Autocorrelation of s')
%% Comments 
% The autocorrelation is also a sinusodial wave, but with an increased
% frequency and decreased amplitude. After one second, some irregularities
% in the sinusodial function start to occur with varying amplitudes that
% eventually decrease to 0 at time two. 
% Taking compactness to mean short duration, where the main energy in the
% waveform is localized in time. Both of these time series, seem to be
% similar in duration, both ending at two seconds. Perhaps the
% autocorrelation has more energy in the waveform since it has a higher
% frequency. 
%% b. Calculate the Cross Correlation of s with a in the frequency domain - plot both the cross correlation and original time series of a

a_fft= fftshift(fft(a, 2*n_s)); 
s_fft_conj = conj(s_fft);
auto_sa = s_fft_conj.*a_fft; 
auto_dt_sa = fftshift(ifft(auto_sa,n_s)/n_s); 
x_r = [1.549899193548388;1.713205645161291;2.051915322580645;2.112399193548387;2.269657258064516;2.448084677419355;2.569052419354839;3.043850806451613]; 
y_r = [-0.006173143979991;-0.006173143979991;0.227844545686244;0.113435897404974;-0.094579826742790;0.155039042234527;0.056231573264339;-0.058177075016932]; 
figure(2) 
subplot(2,1,1)
plot(time,a)
xlabel('Time (s)') 
title('Time Series A')
ylabel('Amplitude')
subplot(2,1,2) 
plot(time,auto_dt_sa); hold on; 
scatter(x_r, y_r, 'r+')
xlabel('Time (s)') 
title('Autocorrelation of S with A')
ylabel('Amplitude')
%[x_r,y_r] = ginput(8)
%% Comments
% The seismogram has the original reflectivity hidden by other convolved
% noise. By cross correlating with the sweep (signal driven into the
% ground, we can removed the the sweep and be left with the Earth's
% reponse. 
% Even with the autocorrelation, I still found it hard to find the 8
% reflectivity peaks since there's still a lot of other noise from the
% ground around. Picking 6 peaks was easy, but I choose two smaller peaks
% around where noise seems to increase. 
%% c. Compute Water Lebel Deconvulution for our Time series S - plot for range of delta against original seismo

delta = [0.001,0.01, 0.1, 1, 10] ; 

wlds_s = zeros(5, 2*n_s); 

for l = 1:length(delta)
    for m = 1:2*n_s
        wld_s(l,m) = a_fft(m)*s_fft_conj(m)/(s_fft(m)*s_fft_conj(m)+delta(l));
        ifft_wlds(l,:) = fftshift(ifft(wld_s(l,:), n_s)/n_s); 
    end
end

color = ['b','k','g','c','m'];
%% PLotting 
figure(3) 
plot(time, abs(ifft_wlds(1,:)), 'b'); hold on; 
plot(time, abs(ifft_wlds(2,:)), 'g'); hold on; 
plot(time, abs(ifft_wlds(3,:)), 'c'); hold on; 
plot(time, abs(ifft_wlds(4,:)), 'm'); hold on; 
plot(time, abs(ifft_wlds(5,:)), 'k'); hold on;
hold on

%% COmments 
% I think something is wrong with my water level deconvolution :(
% 

