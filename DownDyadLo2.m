function d = DownDyadLo2(x,qmf)
%
% This is modified version of DownDyadLo for a faster wavelet filter
% implementation
%
% Modified by Behnood Rasti
% behnood.rasti@gmail.com
%
% DownDyadLo -- Lo-Pass Downsampling operator (periodized)
%  Usage
%    d = DownDyadLo(x,qmf)
%  Inputs
%    x    1-d signal at fine scale
%    qmf    filter
%  Outputs
%    y    1-d signal at coarse scale
%
%  See Also
%    DownDyadHi, UpDyadHi, UpDyadLo, FWT_PO, aconv
%
	d = aconv_2(qmf,x);
	n = size(d,2);
	d = d(:,1:2:(n-1));

    
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
