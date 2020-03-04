function d = DownDyadHi2(x,qmf)
%
%
% This is modified version of DownDyadHi for a faster wavelet filter
% implementation
%
% Modified by Behnood Rasti
% behnood.rasti@gmail.com
%
%
% DownDyadHi -- Hi-Pass Downsampling operator (periodized)
%  Usage
%    d = DownDyadHi(x,f)
%  Inputs
%    x    1-d signal at fine scale
%    f    filter
%  Outputs
%    y    1-d signal at coarse scale
%
%  See Also
%    DownDyadLo, UpDyadHi, UpDyadLo, FWT_PO, iconv
%
	d = iconv_2( MirrorFilt(qmf),lshift2(x));
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
