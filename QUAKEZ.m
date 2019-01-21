clear 
A=[9.000 24.000 14.189 20.950
24.000 13.200 13.679 21.718
33.000 4.800 13.491 21.467
45.000 10.800 14.406 21.713
39.000 27.000 13.075 20.034
54.000 30.000 15.234 20.153
15.000 39.000 13.270 18.188
36.000 42.000 12.239 16.008
27.000 48.000 12.835 15.197
48.000 48.000 14.574 16.280
15.000 42.000 12.624 16.907
18.000 15.000 13.496 21.312
30.000 36.000 10.578 16.664];

% Station coordinates.
xcrd=A(:,1);
ycrd=A(:,2);

% Arrival times for earthquakes 1 and 2.
tq1=A(:,3);
tq2=A(:,4);
vel = 6.0;%km/h;
%% PART A
% Finding the smallest residual 


for n = 0:100 %x
    for i = 0:100 %y
        for k=1:13 
            t(k) = sqrt((n - xcrd(k))^2+(i-ycrd(k))^2)/vel ;
            resid1(k) = tq1(k)-t(k);
            resid2(k) = tq2(k)-t(k);   
        end
        otime1(i+1,n+1) = mean(resid1);
        otime2(i+1,n+1) = mean(resid2);
        resid1 = resid1 - mean(resid1); %bulk shift - gives the origin time
        resid2 = resid2 - mean(resid2); % need to figure out how to get o time 
        resid12 = resid1.^2;
        resid22 = resid2.^2;
        epsil1(i+1,n+1) = sum(resid12);
        epsil2(i+1,n+1) = sum(resid22);
    end 
end
% figure (1)
% imagesc(epsil1)
% figure (2)
% imagesc(epsil2)
%% PART B 
% Finding the best location and origin time for both earthquakes and the variance at that point  

[value1,ind1]= min(epsil1(:));
[value2, ind2]= min(epsil2(:));%ind will give you the index for the value of the matrix as if all of the columns are stacked on each other
ep1 = epsil1(:);
ep2 = epsil2(:);
mine1 = ep1(2861) 
mine2 = ep2(3293) 
Variance1 = mine1/13.0 %Variance of the best time and location
Variance2 = mine2/13.0 
ot1 = otime1(33,29) % Best Origin Times 
ot2 = otime2(61,33)
%% PART C 
% Finding the Uncertainty at the best location 

Uncer1 = mine1/10.0 %Uncertaintities  
Uncer2 = mine2/10.0
[row1,col1] = find(epsil1 == mine1)
[row2,col2] = find(epsil2 == mine2) %locations on 100 by 100 grid

%Need to find the time - mean of the residuals          
%find the chi_squared 
%% PART D
% Finding the Chi Squared for the model space 

for m = 0:100 %x
    for l = 0:100 %y
        for s=1:13 
            tf(s) = sqrt((m - xcrd(s))^2+(l-ycrd(s))^2)/vel ;
            red1(s) = tq1(s)-tf(s);
            red2(s) = tq2(s)-tf(s);   
        end;
        red1 = red1 - ot1; %bulk shift witht he best times from above 
        red2 = red2 - ot2;  
        red12 = red1.^2;
        red22 = red2.^2;
        epsilf1(l+1,m+1) = sum(red12);
        epsilf2(l+1,m+1) = sum(red22);
    end 
end
chi1 = epsilf1./Uncer1;
chi2 = epsilf2./Uncer2; % need to find the chi squared for every point  
bestchi1 = chi1(33,29) 
bestchi2 = chi2(61,33)

%%
%
% In the case where we use the individual standard deviation, we are using
% the expected standard deivation due to errors that occur where as we're
% just using the best standard deviation from where we have the smallest
% variance. 
% For both earthquakes we get a chi squared, which show's that this is the most
% accurate prediction for where the earth quakes are. 

%% PART E
% Plotting the 95% confidence interval for the Chi squared values 

%95% confidence interval 3.94 - 18.31
figure (3) 
%imagesc(epsilf1);hold on;
plot(xcrd,ycrd, '^'); hold on; 
[C1,h1] = contour(chi1, [0,18.31], 'r'); hold on; 
plot(29,33, 'o')
%clabel(C1,h1);
title('Earth Quake 1 95% Confidence Interval');
xlabel('Distance(Km)');
ylabel('Distance(Km)');
%%
% 
% <<EarthQuake1.png>>
% 


figure(4) 
%imagesc(epsilf2); hold on; 
plot(xcrd,ycrd, '^'); hold on;
[C2,h2] = contour(chi2, [0,18.31], 'r'); hold on; 
plot(33,61, 'o');
%clabel(C2,h2);
title('Earth Quake 2 95% Confidence Interval');
xlabel('Distance(Km)');
ylabel('Distance(Km)');

%%
%
%<<EarthQuake2.png>>
%

            