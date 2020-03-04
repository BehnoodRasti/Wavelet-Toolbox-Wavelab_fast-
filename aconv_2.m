function y = aconv_2(f,x)
%
%
% This is a modified version of aconv from Wavelab for a faster wavelet filter
% implementation
%
% (c) 2013 modified by Behnood Rasti
% behnood.rasti@gmail.com
%
%
% aconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = aconv_2(f,x)
%  Inputs
%    f    filter
%    x    2-d matrix
%  Outputs
%    y    filtered result
%
%  Description
%    Filtering by periodic convolution of x with the
%    time-reverse of f columnwise.
%
%  See Also
%    iconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%

n = size(x,2);
p = length(f);
if p < n,
    xpadded = [x x(:,1:p)];
else
    z = zeros(size(x,1),p);
    for i=1:p,
        imod = 1 + rem(i-1,n);
        z(:,i) = x(:,imod);
    end
    xpadded = [x z];
end
fflip = reverse(f);
ypadded = filter(fflip,1,xpadded')';
y = ypadded(:,p:(n+p-1));

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
