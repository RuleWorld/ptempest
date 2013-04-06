function [logpdf] = explogpdf(x, mu)
%EXPLOGPDF Calculate log-pdf (unnormalized) for exponential distribution
%
%   [logpdf] = explogpdf(x, mu)
%   
%   where: exponential(x|mu) = 1/mu*exp( -x/mu )*I(x>0)
%    and   log(exponential(x|mu) ~ -x/mu
%
%   this function is vectorized

logpdf = -(x) ./ mu;
logpdf(find(x<=0)) = -Inf;

