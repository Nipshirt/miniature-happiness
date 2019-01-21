%% Assignment 8 part 2 - e. 
samp_f = 50; 
dt = (1/samp_f); 
N = length(nov1); 
t = [0:N-1].*dt; 

figure(1) 
plot(t, nov1, 'b'); hold on; 
scatter(x,y,40, 'r');
xlabel('Time (s)'); 
ylabel('Signal Strength'); 
title('November 1st Seismogram')
%[x,y] = ginput(3); 
x = [1.023041474654378e+02;2.276497695852534e+02;7.198156682027650e+02];
y = [0;0;0];
swave1_t = x(2); 
pwave1_t = x(1); 
%% f 
d_hp = hpfilt(nov1,dt,2.0, 2, 'TYPE=1');
xt = [5.852534562211982e+02;7.124423963133640e+02]; 
yt = [0;0;0];
x2 = [3.695852534562213e+02;5.852534562211982e+02];
y2 = [0;0;0]; 

figure(2) 
plot(t, d_hp, 'b'); hold on; 
plot(x2(1),y2(1), 'r+'); hold on; 
plot(x2(2), y2(2), 'ro'); 
%xlim([500, 900]);
xlabel('Time (s)'); 
ylabel('Signal Strength'); 
title('Filtered November 1st Seismogram With a Fc = 2.0Hz and N=2')
%[x2,y2] = ginput(2); 
%[xt, yt] = ginput(2); 

swave2_t = xt(2); 
pwave2_t = xt(1); 

%% g. getting the times for the swave and pwave - Estimate distance to earth quake 
vs = 4.6 ;%km/s
vp = 8.0; %km/s
dif_t_as = abs(pwave2_t-swave2_t); 
dif_t_1 = abs( pwave1_t-swave1_t);
ds_as = vs*swave2_t; %in km AFTERSHOCK
dp_as = vp*pwave2_t; % in km AFTERSHOCK
ds = vs*swave1_t; %km
dp = vp*pwave1_t; %km
dif_d_as = ds_as-dp_as; 
dif_d_1 = ds-dp; 
% Time difference between p-wave arrival and s-wave arrivial are about the
% same for both the primary event and the aftershock. 

%% h. 
d_bp = bpfilt(nov1,dt,0.8,1.5,2,'TYPE = 1'); 
xp1 = 1.065668202764977e+02; 
yp1 = 0; 
figure(3) 
subplot(2,1,1) 
plot(t,d_bp);hold on; 
plot(xp1,yp1,'ro')
xlabel('Time (s)'); 
ylabel('Signal Strength'); 

subplot(2,1,2) 
plot(t, d_bp);hold on; 
plot(xp1,yp1,'ro')
xlim([50, 200])
xlabel('Time (s)'); 
ylabel('Signal Strength'); 

%[xp1, yp1] = ginput(1);
%% i. 
d_hpi = hpfilt(nov1,dt,2.0, 2, 'TYPE=1');
xp2 = 97.357910906298000; 
figure(4) 
plot(t, d_hpi); hold on;
plot(xp2,0,'ro')
xlabel('Time (s)'); 
ylabel('Signal Strength'); 
xlim([80, 120]) ; 
%[xp2, yp2] = ginput(1); 

%% Comments 
% From part d. in this question. We saw that a spike spectrum filtered with
% a minimum phase or zero phase filter both had the same result, if we assume 
% this seismogram is made up of various spike series then we won't have a problem
% with either a minimum phase or a zero phase. I choose a minimum phase 
% filter since I'm assuming this time series is casual. 