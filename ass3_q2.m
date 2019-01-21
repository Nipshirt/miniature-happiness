%% Assignment 3 Question 2 Fourier Series 
clear all; close all; 
a0 = 1;
w = 2*pi;  
sm = 0; %This is an even function so there won't be any contribution from the odd function extention
t= linspace(0,1);

n = [1 2 4 16 128 2048]; %n values for separate trunacated series 
a1 = zeros(1, numel(t)); %makes matrix to preallocate for speed 
a2 = zeros(1, numel(t));
a3 = zeros(1, numel(t));
a4 = zeros(1, numel(t));
a5 = zeros(1, numel(t));
a6 = zeros(1, numel(t));
dc = a0/2; %this is the dc term pulled out 
 
for v=1:n(1)%for loop creates a vector from 1:n and steps through each value
    cn(v) = -2/(v*w)*((-1).^((v+1)/2));
    a1 = a1 + cn(v).*cos(v.*w*t); %ask how to make a second for loop work to go through the n values for the series so you don't gotta do this shit yourself 
end

for v=1:n(2)%for loop creates a vector from 1:n and steps through each value
    cn(v) = -2/(v*w)*((-1).^((v+1)/2));
    a2 = a2 + cn(v).*cos(v.*w*t); %ask how to make a second for loop work to go through the n values for the series so you don't gotta do this shit yourself 
end

for v=1:n(3)%for loop creates a vector from 1:n and steps through each value
    cn(v) = -2/(v*w)*((-1).^((v+1)/2));
    a3 = a3 + cn(v).*cos(v.*w*t); %ask how to make a second for loop work to go through the n values for the series so you don't gotta do this shit yourself 
end

for v=1:n(4)%for loop creates a vector from 1:n and steps through each value
    cn(v) = -2/(v*w)*((-1).^((v+1)/2));
    a4 = a4 + cn(v).*cos(v.*w*t); %ask how to make a second for loop work to go through the n values for the series so you don't gotta do this shit yourself 
end

for v=1:n(5)%for loop creates a vector from 1:n and steps through each value
    cn(v) = -2/(v*w)*((-1).^((v+1)/2));
    a5 = a5 + cn(v).*cos(v.*w*t); %ask how to make a second for loop work to go through the n values for the series so you don't gotta do this shit yourself 
end

for v=1:n(6)%for loop creates a vector from 1:n and steps through each value
    cn(v) = -2/(v*w)*((-1).^((v+1)/2));
    a6 = a6 + cn(v).*cos(v.*w*t); %ask how to make a second for loop work to go through the n values for the series so you don't gotta do this shit yourself 
end

A1 = dc+a1;
A2 = dc+a2; 
A3 = dc+a3; 
A4 = dc+a4; 
A5 = dc+a5; 
A6 = dc+a6; 

figure(1) 
plot(t,abs(A1))
hold on;
plot(t,abs(A2))
hold on;
plot(t,abs(A3))
hold on;
plot(t,abs(A4))
hold on;
plot(t,abs(A5))
hold on; 
plot(t,abs(A6))
hold on;
title('Fourier Series of A Even Box Function for Various Values of N') 
xlabel('Time') 
legend('n=1','n=2','n=4','n=16','n=128','n=2048')

%% Comments 
% The Fourier Series returns functions similar to the even box car function
% given to us. For n=1, we have a single sin wave to make up our function.
% This single sin wave cannot resolve the singularity at t=1/4 and t=3/4.
% As we increase n, we converge towards our original function. One issue we 
% run into is the overshoot both before and after the singularity. With 
% increasing n, that issue minimizes. 
% More of a comment on my code. I originally tried to use double for loop
% in order to just have it plot for each value of n, in the nested loops,
% but I found that when I tried to use the nested loop it created "bumps"
% along the function and didn't look quite right, hence the hard coding of
% it all. 

