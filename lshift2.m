function y = lshift2(x)
%
%
% This is modified version of lshift from Wavelab for a faster wavelet filter
% implementation
%
% (c) 2013 modified by Behnood Rasti
% behnood.rasti@gmail.com
%
%
% lshift2 -- Circular left shift applied on columns of matrix x
%  Usage
%    y = lshift2(x)
%  Inputs
%    x   2-d matrix
%  Outputs
%    y   2-d matrix
%        y(:,i) = x(:,i+1) except y(:,n) = x(1)
%
% see also rshift, lshift, rshift2


y = [ x(:, 2:size(x,2) ) x(:,1) ];


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
