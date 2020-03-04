function y = iconv_2(f,x)
%
% This is a modified version of iconv from Wavelab for a faster wavelet filter
% implementation
%
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com
%
% iconv_2 -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = iconv_2(f,x)
%  Inputs
%    f   filter
%    x   2-d matrix
%  Outputs
%    y   filtered result
%
%  Description
%    Filtering by periodic convolution of x with f columnwise.
%
%  See Also
%    aconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%
n = size(x,2);
nr = size(x,1);
p = length(f);
if p <= n,
    xpadded = [x(:,(n+1-p):n) x];
else
    z = zeros(nr,p);
    for i=1:p,
        imod = 1 + rem(p*n -p + i-1,n);
        z(:,i) = x(:,imod);
    end
    xpadded = [z x];
end
ypadded = filter(f,1,xpadded')';
y = ypadded(:,(p+1):(n+p));

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
