function [pdf] = laplacepdf(x, mu, b)
%LAPLACEPDF Calculate pdf for laplace distribution
%
%   [pdf] = laplacepdf(x, mu, b)
%   
%   where: laplace(x|mu,b) = 1/(2*b)*exp( -abs(x-mu)/b )
%
%   this function is vectorized

pdf = 1./(2*b) .* exp(-abs(x-mu)./b );

