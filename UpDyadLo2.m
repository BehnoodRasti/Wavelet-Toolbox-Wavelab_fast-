function y = UpDyadLo2(x,qmf)
%
% This is modified version of UpDyadLo from Wavelab, for a faster wavelet filter
% implementation
%
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com
%
% UpDyadLo2 -- Low-Pass Upsampling operator; periodized
%  Usage
%    y = UpDyadLo2(x,qmf)
%  Inputs
%    x    2-d matrix at coarser scale
%    qmf  filter
%  Outputs
%    y    2-d matrix at finer scale
%
%  See Also
%    DownDyadLo2, DownDyadLo, DownDyadHi, UpDyadHi, IWT_PO, iconv_2
%
y =  iconv_2(qmf, UpSample2(x) );


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
