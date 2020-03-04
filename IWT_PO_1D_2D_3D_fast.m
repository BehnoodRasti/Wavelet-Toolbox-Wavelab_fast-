function [x]=IWT_PO_1D_2D_3D_fast(wc,qmf,L,options)

% IWT_PO_1D_2D_3D_fast -- Inverse orthogonal wavelet transform (periodized)
%  Usage
%     [x]=IWT_PO_1D_2D_3D_fast(wc,qmf,L,options)
%
%  Inputs
%    wc    wavelet coefficients 
%    L     level of decomposition
%    qmf   quadrature mirror filter
%  Outputs
%    x     signal reconstructed from wc 
%
%  Description
%    If wc is the result of a forward wavelet transform, with
%    wc=FWT_PO_1D_2D_3D_fast(x,qmf,L), then [x]=IWT_PO_1D_2D_3D_fast(wc,qmf,L)
%    reconstructs x exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FWT_PO_1D_2D_3D_fast, MakeONFilter
%
% This is a fast implementation for 1D, 2D and 3D orthogonal wavelet transform
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com
if nargin<4 
    options=0;
end
[n1,n2,n3]=size(wc);
if n1==1 && n2>1 && n3==1 
    x = IWT_PO_fast(wc,L,qmf);
elseif n1==1 && n2>1 && n3==1 || options
    x = IWT_PO_fast(wc,L,qmf);
elseif n1>1 && n2>1 && n3==1
    x = IWT2_PO_fast(wc,L,qmf);
elseif n1>1 && n2>1 && n3>1
    x = IWT3_PO_fast(wc,L,qmf);
else
    error('wrong input-input has to be a signal (a row vector) or an image (a two dimentional matrix) or 3d matrix.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 