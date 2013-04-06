function [logpdf] = normlogpdf(x, mu, sigma)
%NORMLOGPDF Calculate log-pdf (unnormalized) for normal distribution
%
%   [logpdf] = normlogpdf(x, mu, sigma)
%   
%   where: normal(x|mu,sigma) = 1/sqrt(2*pi*sigma^2)*exp( -(x-mu)^2/(2*sigma^2) )
%    and   log(normal(x|mu,sigma) ~ -(x-mu)^2/(2*sigma^2)
%
%   this function is vectorized

logpdf = -(x - mu).^2 ./ (2*sigma.^2);
