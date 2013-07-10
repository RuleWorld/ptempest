function [logpdf] = pointlogpdf(x, a)
%POINTLOGPDF Calculate log-pdf (unnormalized) for point distribution
%
%   [logpdf] = pointlogpdf(x, a)
%   
%   where: point(x|a) = 1*I(x==a)
%    and   log(point(x|a) ~ 0 if I(a==a); -Inf otherwise
%
%   this function is vectorized

logpdf = zeros(size(x));
logpdf(find(find(x~=a))) = -Inf;

