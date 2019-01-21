function [b] = deconvolve(c,a)
P = length(c); % number of points in c
N = length(a); % number of points in a
M = P+1-N; % number of points that should be in b
%?? create b

b = zeros(1,M);
for p = 1:M % a for loop for every point that should be in b 
    if p == 1 % for the first b so we can continue with the recursion 
        b(1) = a(1)/c(1) 
    else
        for k = 2:p %from k =1 to p-1 to create the substracted  from C 
            sum_b(k) = a(k)*b(p-k+1); %creates the sum to be substracted from the c
            b(p) = (1/a(1))*(c(p) - sum(sum_b)); % subtracts from the convolution C 
            end 
        end
    end
end

