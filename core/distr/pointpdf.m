function [pdf] = pointpdf(x, a)
%POINTPDF Calculate pdf for point distribution
%
%   [pdf] = pointpdf(x, a)
%   
%   where: point(x|a) = 1*I(x==a)
%
%   this function is vectorized

pdf = ones(size(x));
pdf(find(x~=a)) = 0;

