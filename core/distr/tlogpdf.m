function [logpdf] = tlogpdf(x, v)
%TLOGPDF Calculate log-pdf (unnormalized) for t-distribution
%
%   [logpdf] = tlogpdf(x, v)
%   
%   where: t(x|v) = Gam((v+1)/2)/Gam(v/2)/sqrt(v*pi)/(1+x^2/v)^((v+1)/2)
%    and   log(t(x|v) ~ -(v+1)/2*log(1+x^2/v)
%
%   this function is vectorized

logpdf = -(v+1)/2 .* log(1 + x.^2./v);

