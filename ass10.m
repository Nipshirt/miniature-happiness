%% Assignment 10! 
%% Question 1. 

d1 = [0,0.1]; 
d2 = [1,1.2]; 
d3 = [2,2.7]; 
d4 = [3,6.9]; 
d5 = [4,15.0]; 
B = [ d1' d2' d3' d4' d5']; 
%Get x-vals from data
x_vals = B(1,:);
x_1 = zeros(5,1);
x_2 = x_vals;
x_3 = x_vals.^2;
A = [x_1  x_2' x_3'];
y = B(2,:)'; 

% Form the Normal Matrix 
norm  = A'*A; 

% Form the Vector 
y_vec = A'*y; 

%Obtain the least square solution m =(norm)^(-1)(y_vec) 
m = y_vec\norm; 

%% Plotting answer above 
x = linspace(0,4); 
y = m(1) + m(2)*x +m(3)*x.^2; 

figure(1)
plot(x,y) ;hold on;
scatter(B(1,:),B(2,:))
ylabel('Y') 
xlabel('X') 
title('Least Square fit for a Quadratic of y=c+bx+ax^2')


%%

A2 = [0 , -1; 0,0 ; 0,1] 
k = A2'*A2
d = [0.1;1.0;1.9]; 

