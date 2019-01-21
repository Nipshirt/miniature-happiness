%% Assignment 6 
% Generate a set of time series of a sinusodial function with varying
%  natural frequencies 
samp_rate = 20;%sps
dt  = 1/20 ; %s
f_n = 1/(2*dt); %hz 
N = 1024 ; 
t = [0: dt: (N-1)*dt];
freq = [1:1:19]; 

% double for loops in order to make a matrix with each sample frequency
% in a row with 1024 points corresponding to different times in our series 

for n = 1:length(freq)
    for m = 1:length(t)
        y(n,m) = sin(2*pi*freq(n).*t(m)) ;
    end
    figure(n)
    plot(t,y(n,:))
    title(sprintf( 'Time Series for Sin Function with f = %g', n))
    xlim([0, 2*pi]) 
    xlabel('Time (s)') 
end

%% Fourier Tranform  & plot? 
% TAKING THE FOURIER TRANSFORM THEN SHIFTING IT SO IT'LL BE IN THE RIGHT
% ORDER
fn = [[0:((N/2)-1)], [-N/2: -1]]; %natural frequencies 
for p = 1:length(freq) 
    fft_y(p,:) = fft(y(p,:), N);%N^2 to pad to increase resolution 
    sfft_y(p,:) = fftshift(fft_y(p,:));
    amp_y(p,:) = abs(sfft_y(p,:)); % amplitude
    phas_y(p,:)= atan2(imag(sfft_y(p,:)),real(sfft_y(p,:))); %phase spectrum 
end 
%% plotting amplitude spectrum 

for P = 1:length(freq) 
    %plotting amp 
    figure(P)
    plot(fftshift(fn), amp_y(P,:)/N)
    title( sprintf( 'Amplitude Spectum for Sin Function with f = %g', P))
    xlabel('Frequencies (Hz)') 
    ylabel(' Amplitude')
end 

%% Plotting Phase Spectrum 

for L = 1:length(freq)
    %plotting phase spectrum 
    figure(L) 
    plot(fftshift(fn), phas_y(L,:))
    title( sprintf( 'Phase Spectum for Sin Function with f = %g', L))
    xlabel('Frequencies (Hz)') 
    ylabel('Phase (Radians)')
end 