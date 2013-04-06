function [logpdf] = lognlogpdf(x, mu, sigma)
%LOGNLOGPDF Calculate log-pdf (unnormalized) for lognormal distribution
%
%   [logpdf] = lognlogpdf(x, mu, sigma)
%   
%   where: lognormal(x|mu,sigma) = 1/(x*sqrt(2*pi*sigma^2) * exp( -(log(x)-mu)^2/(2*sigma^2) )*I(x>0)
%    and   log(lognormal(x|mu,sigma) ~ (-log(x) -(log(x)-mu)^2/(2*sigma^2))*I(x>0)
%
%   this function is vectorized

logpdf = -log(x) - (log(x) - mu).^2 ./ (2*sigma.^2);
logpdf(find(x<=0)) = -Inf;
