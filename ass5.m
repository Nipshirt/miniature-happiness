%% Assignment 5 
%% Question 1 - use fft to take the dft of a box car function and plot amplitude and phase spectrum 
N = 64; % fft should go from -N/2 to N/2 
M = [1, 2, 4, 8, 16, 32, 64] ; 
wn=[-32:31]*2*pi/N; %creates the frequency from -N/2 to N/2
%creates a matrix with the box step function for each M, for each 
for h = 1:length(M) 
for k = 1:N+1 % plus one due to matlab not indexing starting at one 
    if k >= 1 && k < M(h)+1 
        bk(h,k) = 1;
    elseif k >= M(h)+1 && k < N+1 
        bk(h,k) = 0 ;
    end 
end
% fourier transform for each box step function 
Bn(h,:)=fft(bk(h,:),64)./N;
Rn(h,:)=fftshift(abs(Bn(h,:)));
%getting the phase of Bn
for z = 1:length(Bn)
    phase(h,z) = atan2(imag(Bn(h,z)),real(Bn(h,z)));
end
figure (h)
subplot(3,1,1) 
plot(bk(h,:), 'r')

xlabel('Points (k)')
ylabel('Amplitude')
title('Box Car Function')

subplot(3,1,2) 
plot(wn, Rn(h,:))

title('Amplitude Plot') 
xlabel('Frequency (Hz)')
ylabel('Amplitude')

subplot(3,1,3) 
plot(phase(h,:))

title('Phase Spectrum')
xlabel('Frequency (Hz))')
ylabel('Phase (radians)')

end
%% Discussion Of Different Box Car Functions with Increasing Point of Change 

% What do you notice about the width of the central peak of the amplitude
% spectrum relative to the width of the box car function ?
% 
% As you increase M, which created the width for the box car function with
% respect to zero, the central peak of the amplitude spectrum decreases
% until M = 64, which is when we get the largest amplitude of 1. With
% increasing width, there is more frequencies needed at different
% amplitudes to make up the signal bk. 

%What do you notice about the phase spectrum? 
% As M increases, the phase spectrum includes more oscillations of the
% phase at each frequency. Each phase spectrum has a M number of peaks. 
% Also with an increase in M, we see a smaller
% period of oscillations. The period of the oscillations are approximately 
% M/62 

 %% Question 2 - properties of DFT 
 N2= 64;
 dt = 1; 
 index = 1:64; 
 % creating the discrete function 
 bj(1:8) = 1;
 bj(9:17) = -1;
 bj(18:64) = 0;

 %fourier transform 
 BJ = fftshift(fft(bj, N2)); 
omegan=[0:N2/2-1,-N2/2:-1]*2*pi/(N2*dt) ; 
figure (11)
subplot(2,1,1) 
plot(bj)
subplot(2,1,2)
plot(omegan(1:64),abs(BJ(1:64)))
%% a)shift property - shift by 10s and 58s 
 
aj10 = exp(-i.*omegan.*10);
aj58 = exp(-i.*omegan.*58); 
AJ10 = BJ.*aj10; 
AJ58 = BJ.*aj58;
i_aj10 = ifft(ifftshift(AJ10), 'symmetric');
i_aj58 = ifft(ifftshift(AJ58), 'symmetric');

figure(22) 
subplot(2,1,1) 
plot(i_aj10)
xlabel('Time (s)')
ylabel('Amplitude')
title('Discrete Function Shifted by 10s')
subplot(2,1,2) 
plot(i_aj58)
xlabel('Time (s)')
ylabel('Amplitude')
title('Discrete Function Shifted by 58s')

%% b) Time Reversal - Show that time reversal is the same as CONJUGATE 

BJ_con = conj(BJ); 
bj_con = ifft(ifftshift(BJ_con)); 

figure(33) 
subplot(2,1,1) 
plot(omegan,BJ_con)
subplot(2,1,2) 
plot(bj_con)
%% c) convolution 

conv_BJ = convolve(ifftshift(BJ),ifftshift(BJ)); %convolving bk in frequency domainf
i_conv_BJ = ifft(conv_BJ, N2, 'symmetric');
multi_bj = bj.*bj;

figure(44) 
subplot(2,1,1)
plot(i_conv_BJ)
subplot(2,1,2)
plot(multi_bj)

%% 2.c
% There is no need to pad with zeros when convolving the BJ since they are
% the same length so the length of the convolved function is just twice the
% length of bk 
%%  Differentiation of bk 
dif_bj = diff(bj); 
d_BJ= 1i*omegan.*BJ; 
i_d_BJ = ifft(ifftshift(d_BJ), N2); 

figure(55) 
subplot(2,1,1)
plot(dif_bj) 
subplot(2,1,2) 
plot(real(i_d_BJ))
%% Integration of bk 

int_bj = cumtrapz(bj); 
int_BJ = BJ./(1i*omegan); 
int_BJ(1) = 0; 
i_int_BJ = ifft(ifftshift(int_BJ), N2);

figure(66) 
subplot(2,1,1) 
plot(int_bj)
subplot(2,1,2) 
plot(abs(i_int_BJ))

%% Parsevals theorem

par_bj = sum(bj.^2)
par_BJ = (1/N2)*sum(BJ.^2)


