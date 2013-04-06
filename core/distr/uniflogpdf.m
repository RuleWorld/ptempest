function [logpdf] = uniflogpdf(x, a, b)
%UNIFLOGPDF Calculate log-pdf (unnormalized) for uniform distribution
%
%   [logpdf] = uniflogpdf(x, a, b)
%   
%   where: uniform(x|a,b) = 1/(b-a)*I(a<=x<=b)
%    and   log(uniform(x|a,b) ~ 0 if I(a<=x<=b); -Inf otherwise
%
%   this function is vectorized

logpdf = zeros(size(x));
logpdf(find(or(x<a, b<x))) = -Inf;
