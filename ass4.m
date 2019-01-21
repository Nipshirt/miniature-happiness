%% Assignment 4 
close all; clear all;
%% Question 2 - use convolve.m to check hand worked answer 

a2 = [6, -5, 1, 0, -1, 3, -2]; 
b2 = [1/2, 1/2]; 

c2 = convolve(a2,b2); 

%% Question 3 
c3 = [1,5,14,27,34,32,17,4]; 
a3 = [1,2,3,1,2];

b3 = deconvolve(c3,a3); 

%% Question 4 
a1 = [1 -1];
A1 = fft(a1);
subplot(2,3,1);
plot(abs(A1));
title('a = [1,?1]');

a2 = [1 -1 0];
A2 = fft(a2);
subplot(2,3,2)
plot(abs(A2));
title('a = [1,?1,0]');


a3 = [1 -1 0 0];
A3 = fft(a3);
subplot(2,3,3)
plot(abs(A3));
title('a = [1,?1,0,0]');

a4 = [1 -1 0 0 0];
A4 = fft(a4);
subplot(2,3,4)
plot(abs(A4));
title('a = [1,?1,0,0,0]');

a5 = [1 -1 0 0 0 0];
A5 = fft(a5);
subplot(2,3,5)
plot(abs(A5));
title('a = [1,?1,0,0,0,0]');

% The first plot is a linear fit between 0 and two, which is the amplitude
% of the sequence [1, -1]
% As you continue to add on zeros to the sequence you are changing the N,
% amount of points, which has an affect to the over all amplitude by (1/N)
% and in the exponential term which controls the period of the FFT. 
% The addition of the zeros into the sequence, doesn't change the over all
% max amplitude but changes the change in amplitude and how fast it does or
% doesn't change from n-point to n-point. Odd number of zeros creates an
% non-change between the max amplitude and the next n where as a even
% number of zeros creates an negative change of the amplitude. 