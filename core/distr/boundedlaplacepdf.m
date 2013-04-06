function [pdf] = boundedlaplacepdf(x, mu, b, l, u)
%BOUNDEDLAPLACEPDF Calculate pdf for bounded laplace distribution
%
%   [pdf] = laplacepdf(x, mu, b, l, u)
%   
%   where: laplace(x|mu,b,l,u) = 1/Z(b,l,u)*exp( -abs(x-mu)/b )*I(l<=x<=u)
%
%          with Z = Int[l,u]( exp(-abs(x-mu)/b )dx
%
%   this function is vectorized


% Note that CDF(x|mu,b) = 1/2*exp( (x-mu)/b ),    if x < mu
%                         1 - 1/2*exp(-(x-mu)/b), if x >= mu


% compute probability that a laplace random variable is between l and u.
Plu =  (1/2*exp((u-mu)./b).*(u<mu) + (1 - 1/2*exp(-(u-mu)./b)).*(u>=mu)) ...
      -(1/2*exp((l-mu)./b).*(l<mu) + (1 - 1/2*exp(-(l-mu)./b)).*(l>=mu));
% compute PDF for bounded laplacian (only normalization factor changes)
pdf = 1./(2*b) .* exp(-abs(x-mu)./b );
% renormalize to account for bounds
pdf = pdf./Plu;
% set pdf to zero if x is out of bounds
pdf(find(or(x<l,u<x))) = 0;

