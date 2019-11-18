function [I, num_zero] = svdinv(X,threshold)
%Syntax: [I, num_zero] = svdinv(X,threshold)
%
%Calculate inv using svd, skip svd's with ratio to the maximum-svd less than the given threshold
if nargin<2
    threshold = 1e-10;
end;

[U,S,V] = svd(X);
 
d = diag(S);
d = abs(d/max(d))>threshold;

num_zero = 0;
for k=1:length(d)
    if(d(k)>0)
        S(k,k) = 1.00/S(k,k);
    else
        S(k,k) = 0;
        num_zero = num_zero+1;
    end;
end;

I = V*S*U';
