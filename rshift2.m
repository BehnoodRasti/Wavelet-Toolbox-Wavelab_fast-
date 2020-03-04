function y = rshift2(x)
%
%
% This is modified version of rshift from Wavelab for a faster wavelet filter
% implementation
%
% (c) 2013 modified by Behnood Rasti
% behnood.rasti@gmail.com
%
%
% rshift2 -- Circular right shift applied on columns of matrix x
%  Usage
%    y = rshift2(x)
%  Inputs
%    x   2-d matrix
%  Outputs
%    y   2-d matrix
%        y(:,i) = x(:,i-1) except y(:,1) = x(:,n)
%
% see also rshift, lshift, lshift2


n = size(x,2);
y = [ x(:,n) x(:, 1: (n-1) )];


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
