function [x] = boundedlaplacernd( mu, b, l, u )
%BOUNDEDLAPLACERND Sample from the bounded Laplace distribution
%
%   [x] = laplacernd(mu,b,l,u)
%
%   where: laplace(x|mu,b,l,u) = 1/Z(b,l,u)*exp( -abs(x-mu)/b )*I(l<=x<=u)
%
%          with Z = Int[l,u]( exp(-abs(x-mu)/b )dx
%
%   this function is vectorized

% initialize shape of x
x = zeros(size(mu));

% sample x from laplace and reject anything out of bounds
idx = 1:numel(mu);
while ~isempty(idx)

    U = rand(size(idx)) - 0.5;
    x(idx) = mu(idx) - b(idx).*sign(U).*log(1 - 2*abs(U));
    % find x that are out of range
    idx = find(or(x<l,u<x));

end

