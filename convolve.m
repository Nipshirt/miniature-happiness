function [c] = convolve(a,b)
%Convulution of two vectors
N=length(a);
M=length(b); 
c=zeros(1,N+M-1);
a = [a, zeros(1,M-1)]; % padding the time sequence 
b = [b, zeros(1,N-1)]; % padding the time sequence 


 for p = 1:N+M-1 %for the length of C 
    for k = 1:p 
        c(p) = c(p) + a(k)*b(p-k+1); 
    end
    end
 end 
    

