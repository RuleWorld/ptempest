function [x] = laplacernd( mu, b )
%LAPLACERND Sample from the Laplace distribution
%
%   [x] = laplacernd(mu,b)
%
%   where: laplace(x|mu,b) = 1/(2*b)*exp( -abs(x-mu)/b )
%
%   this function is vectorized

U = rand(size(mu)) - 0.5;
x = mu - b.*sign(U).*log(1 - 2*abs(U));

