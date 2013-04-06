function [logpdf] = boundedlaplacelogpdf(x, mu, b, l, u)
%BOUNDEDLAPLACELOGPDF Calculate log-pdf (unnormalized) for bounded laplace distribution
%
%   [logpdf] = laplacelogpdf(x, mu, b, l, u)
%   
%   where: laplace(x|mu,b,l,u) = 1/Z(b,l,u)*exp( -abs(x-mu)/b )*I(l<=x<=u)
%
%          with Z = Int[l,u]( exp(-abs(x-mu)/b )dx
%
%   this function is vectorized

logpdf = -abs(x-mu)./b;
logpdf(find(or(x<l,u<x))) = -Inf;
