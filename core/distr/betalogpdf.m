function [logpdf] = betalogpdf(x, a, b )
%BETALOGPDF Calculate log-pdf (unnormalized) for beta distribution
%
%   [logpdf] = betalogpdf(x, a, b)
%   
%   where: B(x|a,b) = x^(a-1)*(1-x)^(b-1)/B(a,b)*I(0<=x<=1)
%    and   log(B(x|a,b)) ~ (a-1)*log(x) + (b-1)*log(1-x) 
%
%   this function is vectorized
%
%   NOTE: if the logpdf evaluates to Inf, then we return -Inf so the 
%   downstream applications won't go nuts. This will occur at the boundaries,
%   so it's not a serious problem.

if (b==1)
    logpdf = (a-1).*log(x);
    logpdf(find(or(x<=0,1<x))) = -Inf;
elseif (a==1)
    logpdf = (b-1).*log(1-x);
    logpdf(find(or(x<0,1<=x))) = -Inf;
else
    logpdf = (a-1).*log(x) + (b-1).*log(1-x);
    logpdf(find(or(x<=0,1<=x))) = -Inf;
end


