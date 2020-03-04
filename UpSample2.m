function y = UpSample2(x,s)
%
% This is modified version of UpSample from Wavelab for a faster wavelet filter
% implementation
%
% (c) 2013 modified by Behnood Rasti
% behnood.rasti@gmail.com
%
% UpSample2 -- Upsampling operator
%  Usage
%    y = UpSample2(x,s)
%  Inputs
%    x   2-d matrix, of length n
%    s   upsampling scale, default = 2
%  Outputs
%    y   2-d matrix, of column length s*n with zeros
%        interpolating alternate samples columnwise
%        y(:,s*i-1) = x(:,i), i=1,...,n
%

if nargin == 1, s = 2; end
n = size(x,2)*s;
y = zeros(size(x,1),n);
y(:,1:s:(n-s+1) )=x;



%
% Copyright (c) 2006. David Donoho
%

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
