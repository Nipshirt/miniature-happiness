%% Assignment 7 
%
%% a. The surface waves were so large, the instrument was saturated and took about 2 hours
% to recover. There is some non-linear behaviour before the signal settles. First remove this
% non-linear and saturated response from the beginning of the trace. Include time-domain plots
% (plotted versus time, not index number) of the raw data, the raw data with saturation removed,
% and the amplitude spectrum of the raw data with clipping removed (plotted versus frequency
% using fftshift from the negative Nyquist to positive Nyquist frequencies)
cmo_d = importdata('cmo.mat') ;
esk_d = importdata('esk.mat') ;
cmo_cut = 375; 
esk_cut = 375;
cmo_l = cmo_d(cmo_cut:3600);
esk_l = esk_d(esk_cut:3600); 
samp_in = 20; %s 
nq = 1/(2*samp_in); 
t = samp_in.*[1:length(cmo_d)];
t_h = t./3600 ; 
NC = length(cmo_l);
NE = length(esk_l);
N = 3600; 
cmo_fft = fft(cmo_l);
esk_fft = fft(esk_l);
cmo_ffts = fftshift(cmo_fft);
esk_ffts = fftshift(esk_fft);

fc =[0:NC/2,-NC/2+1:-1].*(1/samp_in);
fe = [0:NE/2,-NE/2+1:-1].*(1/samp_in);

figure(1) 
subplot(2,1,1) 
plot(t_h(cmo_cut:3600), cmo_l) 
xlabel('Time (hr)') 
title('Altered CMO Data')

subplot(2,1,2) 
plot(fftshift(fc), abs(cmo_ffts))
xlabel('Frequency (Hz)') 
ylabel('Amplitude')
title('Amplitude Spectrum of CMO Data')

figure(2) 
subplot(2,1,1) 
plot(t_h(esk_cut:3600),esk_l)
xlabel('Time (hr)') 
title('Altered ESK Data')

subplot(2,1,2) 
plot(fftshift(fe), abs(esk_ffts))
xlabel('Frequency (Hz)') 
ylabel('Amplitude')
title('Amplitude Spectrum of ESK Data')

%% b. Can you see the Earth is Flat!!!? tides now in this clipping-free data? To filter out the tides, es-
% timate their dominant frequency (ie periods of roughly 6-12 hrs) and consider the range of
% frequencies represented in the first few points of the Fourier transform of the data (ie fn =
% [0:N/2,-N/2+1:-1]/(T) for even N). Don’t forget that you need to delete both negative and
% positive frequencies in the Fourier spectrum. As suggested, filter out the tides and examine your
% “clean” seismogram. Include both time domain (versus time) and frequency domain (versus
% frequency) plots of the “cleaned” data. Comment on the difference between the tide-saturated
% and tide-free time series and amplitude spectrum plots. 

tide_freq = 1./([6:12].*3600);

cmo_fft(1) = 0; 
esk_fft(1) = 0;
cmo_fft(2) = 0; 
esk_fft(2) = 0;
cmo_fft(NC) = 0; 
esk_fft(NE) = 0;

cmo_ut = ifft(cmo_fft); 
esk_ut = ifft(esk_fft); 

figure(3) 
subplot(2,1,1) 
plot(t_h(cmo_cut:3600), cmo_ut) 
xlabel('Time (hr)') 
title('Altered CMO Data  Without Tidal Noise')

subplot(2,1,2) 
plot(fftshift(fc), abs(fftshift(cmo_fft)))
xlabel('Frequency (Hz)') 
ylabel('Amplitude')
title('Amplitude Spectrum of CMO Data Without Tidal Noise')

figure(4) 
subplot(2,1,1) 
plot(t_h(esk_cut:3600),esk_ut)
xlabel('Time (hr)') 
title('Altered ESK Data  Without Tidal Noise')

subplot(2,1,2) 
plot(fftshift(fe), abs(fftshift(esk_fft)))
xlabel('Frequency (Hz)') 
ylabel('Amplitude')
title('Amplitude Spectrum of ESK Data Without Tidal Noise')

%% c) The bursts of energy are surface waves that have travelled repeatedly round the great circle
% path through Oaxaca, Mexico and College, Alaska. Remember that every other burst of energy
% has travelled in opposite directions; if cR is the Rayleigh wave speed, rE the Earth’s radius,
% and  the angular distance between event and station, the first wavetrain arrives after a time
% rE/cR, the second after time rE(2 ? )/cR, etc. Although you have lost the first few
% arriving wavetrains due to instrument saturation, you can still estimate the distance  and,
% given rE = 6371 km, Rayleigh wave speed cR. Find these two values (Hint: You can use
% the ginput function to determine the timing of the surface wave “bursts” on the seismograms
% (type help ginput for more information). Choose the first 4 bursts of energy and look at
% the differences in timing between the bursts to figure out the distance  and Rayleigh wave
% velocity. You will have to determine which are the more important relative timing bursts to
% use.) Please hand in a time-domain plot indicating where you took your times (T1, T2, T3, T4)
% and the work required to compute  and cR.


figure(5)
plot(t_h(cmo_cut:3600), cmo_ut); hold on;
scatter(t_c, a_c, 22, 'r')
xlabel('Time (hr)') 
title('Altered CMO Data')
% [t_c, a_c] = ginput(4) ;

figure(6)
plot(t_h(esk_cut:3600),esk_ut); hold on; 
scatter(t_e, a_e, 22, 'r')
xlabel('Time (hr)') 
title('Altered ESK Data')

% [t_e,a_e] = ginput(4); 
% save('t_e'); save('a_e');save('t_c');save('a_c');

%% d) As in (d) on page 64 of TSAITG, identify the peaks in the frequency range 1.9-3.0 mHz and
% measure their frequencies on your tide-free amplitude spectrum plot. Compare the frequencies
% 2 you observe with those predicted in the text book, and identify the modes. Before you take
% the fft of the time series you will want to first window it with an appropriate taper window to
% avoid spectral leakage, and then pad with zeros. Use the Hann window which for a time series
% with N points as defined by 0.5*(1-cos(2*pi*[0:N-1]/(N-1)));. Compare the spectrum
% computed in this way with the spectrum you get with no taper (e.g. using a boxcar window
% function). Then compare your results with the predicted frequencies that are provided in the
% book. Include both the tapered and non-tapered amplitude spectrums plotted versus frequency,
% with comments on the difference between the two. Also include the normal mode frequencies
% that you determined from your amplitude spectrum (include a detailed, zoomed in plot of the
% frequency range to demonstrate that the normal modes are present)

