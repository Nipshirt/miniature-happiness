%% Assignment 11 
%% a. least squares fit of the data y = ax+b 
A1 = ones(24,1); 
% Matrix A 
A = [A1 x'];
% least squares soln 
m = A\y'; 


%% b. take out error by eye 
N=length(x);
error_e = [4,17] ;%errors index taken by eye 
y_e = [y(1:error_e(1)-1) y(error_e(1)+1:error_e(2)) y(error_e(2)+2:end) ];
x_e = [x(1:error_e(1)-1) x(error_e(1)+1:error_e(2)) x(error_e(2)+2:end)];
A_e = [A(1:error_e(1)-1,:); A(error_e(1)+1:error_e(2),:) ;A(error_e(2)+2:end,:)]; 

m_e = A_e\y_e'; 

%% d. Construct Ce then invert for a straight line again and plot it against the data 
sig_d = [sig(1:error_e(1)-1) sig(error_e(1)+1:error_e(2)) sig(error_e(2)+2:end)]';
% don't think the below process is right- tell Marek 
C = sig_d*sig_d'; 
Ce = C.^2; 


Ce_d = diag(sig_d); 
m_d = A_e'*inv(Ce_d)*A_e;
m_f = A_e'*inv(Ce_d)*y_e'; 
m_ls = inv(m_d)*m_f;
%% e. Weighted Residual 

e_a =y_e'-A_e*m_ls; 
e_b = inv(Ce_d)^(1/2); 

e_vec = abs(e_b*e_a); 
ind = find(e_vec>4); 
sig_dd = [sig_d(1:ind-1); sig_d(ind+1:end)]; 
Ce_dd = diag(sig_dd); 

y_d = [y_e(1:ind-1) y_e(ind+1:end)]; 
x_d = [x_e(1:ind-1) x_e(ind+1:end)]; 
A_d = [A_e(1:ind-1,:); A_e(ind+1:end,:)]; 

m_dd = A_d'*inv(Ce_dd)*A_d;
m_fd = A_d'*inv(Ce_dd)*y_d'; 
m_lsd =inv(m_dd)*m_fd;

%% my try -- weight solution witghout taking out any points 
Ce_m = diag(sig); 
e_ma =y'-A*m; 
e_mb = inv(Ce_m)^(1/2); 
error_m = abs(e_mb\e_ma);
ind_m  = find(error_m>10, 3);
sig_m = [sig(1:ind_m(1)-1) sig(ind_m(1)+1:ind_m(2)-1) sig(error_m(2)+1:error_m(3)-1)]; 
Ce_my = diag(sig_m); 

y_m = [y(1:ind_m(1)-1) y(ind_m(1)+1:ind_m(2)-1) y(error_m(2)+1:error_m(3)-1)];
x_m = [x(1:ind_m(1)-1) x(ind_m(1)+1:ind_m(2)-1) x(error_m(2)+1:error_m(3)-1)];
A_m = [A(1:ind_m(1)-1,:); A(ind_m(1)+1:ind_m(2)-1,:) ;A(error_m(2)+1:error_m(3)-1,:)]; 
 
m_dm = A_m'*inv(Ce_my)*A_m;
m_fm = A_m'*inv(Ce_my)*y_m'; 
m_lsm = inv(m_dm)*m_fm;

%% Plotting Least Squares 

subplot(3,1,1)
y_model = m(1) + m(2).*x; 
figure(1) 
scatter(x,y); hold on; 
plot(x,y_model);
xlabel('X Values') 
ylabel('Y Values') 
legend('Data','Least Squares Fit')
title('A.Least Squares Fit to y=mx+b')  

y_e_model = m_e(1) + m(2).*x;
subplot(3,1,2)
scatter(x_e,y_e); hold on; 
plot(x,y_e_model);
xlabel('X Values') 
ylabel('Y Values') 
legend('Data','Least Squares Fit')
title('B.Least Squares Fit to y=mx+b with Two Data Points Removed')  

y_weighted = m_ls(1)+ m_ls(2).*x; 
subplot(3,1,3) 
scatter(x_e,y_e); hold on; 
plot(x,y_weighted);
xlabel('X Values') 
ylabel('Y Values') 
legend('Data','Least Squares Fit')
title('d. Weighted Least Squares Fit to y=mx+b without 2 Data Points')  
figure(2)
y_final = m_lsd(1)+ m_lsd(2).*x; 
subplot(1,1,1) 
scatter(x_d,y_d); hold on; 
plot(x,y_final);
xlabel('X Values') 
ylabel('Y Values') 
legend('Data','Least Squares Fit y')
title('e. Weighted Least Squares Fit to y=mx+b without 3 Data Points')


%% Comments 
act = [5;2];
dif_b = abs(act - m_e) %difference between the actual slope and intercepts with the fit where I took out two points
dif_e = abs(act- m_ls) % same difference but with the weighted least squares answer 
dif_mt = abs(act-m_lsd) % same difference as dif_b but including my try of weighting all data points are first then taking out the 3 largest 
% Looking at the differences we see that the procedure done in part e 
% produces the best results, even when you look at the weighted errors at
% the beginning with all the data points(assuming I did all the procedures
% correctly). Taking out two points by eye created a decent fit, but
% picking a third point would have been very difficult, as it would have
% skewed the fit 