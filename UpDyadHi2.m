function y = UpDyadHi2(x,qmf)
%
%
% This is modified version of UpDyadHi from Wavelab for a faster wavelet filter
% implementation
%
% (c) 2013 Modified by Behnood Rasti
% behnood.rasti@gmail.com
%
%
% UpDyadHi2 -- High-Pass Upsampling operator; periodized
%  Usage
%    y = UpDyadHi2(x,qmf)
%  Inputs
%    x    2-d matrix at coarser scale
%    qmf  filter
%  Outputs
%    y    2-d matrix at finer scale
%
%  See Also
%    DownDyadLo2,DownDyadLo, DownDyadHi, UpDyadLo, IWT_PO, aconv_2,UpSample2
%

y = aconv_2( MirrorFilt(qmf), rshift2( UpSample2(x) ) );


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
