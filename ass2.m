%% Assignment 2 
close all 
clear all 
t = 0:1023; % time series 
aa = rand(size(t)); %funcrion for 2a
%%
% First plot shows a function of random points between zero and one for each point in the time
% series.

ab = exp(-max(abs(t),500)/100*50); % function for 2b 
%%
% $$ e^((-1*max(abs(t),500)/100)*50) $$
% 
%Second plot shows a natural exponential with the maximum value of the time series divided by
% one half till time is 500 when the max function picks the time series
% over the the scaler 500. 
ac = 10*exp(-((t-500)/50).^2) ; %function for 2c
%%
% $$ 10*e^(-((t-500/50)^2) $$ 
% 
%Plot 3 exponential starts at a very small number then
% has increases rapidly at towards the peak at t=500, then decreases in the
% same function till the end of the time series. 
ad = exp(-abs(t-500)/100).*cos(2*pi*t/100); %function for 2d 
%%
% 
% $$ e^(-abs(t-500)/100)*cos(2*pi*t/100) $$
% 
%The 4th figure shows a function that seems to be a gaussian wave packet 
% with a period of 100 seconds and a peak at 500 seconds. 

ae = exp(-abs(t)/100).*sqrt(sin(2*pi*t/40)+1); %function for 2e 

%%
% $$ e^(-abs(t)/100)*sqrt(sin(2*pi*t/40)+1) $$
% 
%%This function is a periodic function with a constant decay factor of
% -1/100. The periodicity comes from the sin function with a period of 40
% seconds, which is positive always due to being a real function under a
% square root. 

%% Spike Train  
b_r = zeros(1,2048); 
ib = [50 150 250 1250]; 
nb = [-1 -2 0.5 -1];
b_r(ib) = nb ; 

figure(11)
plot(b_r); 
xlabel('Time (s)');
title('Plot of Spike Train'); 
%% Convolution of two functions 

aa_c = convolve(aa,b_r); 
ab_c = convolve(ab,b_r); 
ac_c = convolve(ac,b_r); 
ad_c = convolve(ad,b_r); 
ae_c = convolve(ae,b_r); 
ea_c = convolve(b_r,ae); % Here showing the communtativity principle of convolutions 

% Plotting 
figure(1) 
AA = subplot(2,1,1);
plot(aa); 
title('Plot of Original Function (top) and Plot (bottom) of Function Convolved with Spike Train')
xlabel('Time');
AA_C = subplot(2,1,2); 
plot(aa_c);
xlabel('Time'); 

figure(2) 
AB = subplot(2,1,1);
plot(ab); 
xlabel('Time');
title('Plot of Original Function (top) and Plot (bottom) of Function Convolved with Spike Train')
AB_C= subplot(2,1,2); 
plot(ab_c); 
xlabel('Time');

figure(3)
AC = subplot(2,1,1); 
plot(ac); 
xlabel('Time');
title('Plot of Original Function (top) and Plot (bottom) of Function Convolved with Spike Train')
AC_C = subplot(2,1,2); 
plot(ac_c); 
xlabel('Time');

figure(4) 
AD = subplot(2,1,1); 
plot(ad); 
xlabel('Time');
title('Plot of Original Function (top) and Plot (bottom) of Function Convolved with Spike Train')
AD_C = subplot(2,1,2); 
plot(ad_c); 
xlabel('Time');

figure(5)
AE = subplot(2,1,1); 
plot(ae); 
xlabel('Time');
title('Plot of Original Function (top) and Plot (bottom) of Function Convolved with Spike Train')
AE_C = subplot(2,1,2); 
plot(ae_c);
xlabel('Time') ;

% 
% %% Communtativity Principle 
figure(6) 
plot(ae_c(t+1),'b+');
hold on;
plot(ea_c(t+1),'k-'); 
title('Convolution of Two Identical Functions to Show Communtativity Principle'); 
xlabel('Time'); 



