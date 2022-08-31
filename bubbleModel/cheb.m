function [ch]=cheb(N,x)

% subroutine used by waterprops
%computes chebyshev polynominal to order N,
% with first value in output matrix corresponding to zero order value
% written by S. Meers, modified by Gary Robb (20/09/05)
ch=zeros(1,N+1);
ch(1)=1;
ch(2)=x;
for n=2:N
    ch(n+1)=2*x*ch(n)-ch(n-1);
end