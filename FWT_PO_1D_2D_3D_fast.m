function wc=FWT_PO_1D_2D_3D_fast(x,qmf,L,options)

% FWT_PO_1D_2D_3D_fast -- Applies 1-d 2-d 3-d Forward orthogonal Wavelet Transform (periodized)
%                         corresponding to the size of the input data.   
%  Usage
%     wc = FWT_PO_1D_2D_3D_fast(x,qmf,L,options)
%  Inputs
%    x    1-d or 2-d or 3-d signal;
%    L    decomposition level
%    qmf  quadrature mirror filter (orthonormal)
%  Outputs
%    wc   1-d or 2-d or 3-d wavelet transform.
%
%  Description
%    1. qmf filter may be obtained from MakeONFilter   
%    2. usually, length(qmf) < 2^(L+1)
%    3. To reconstruct use IWT_PO_1D_2D_3D_fast
%
%  See Also
%    IWT_PO_1D_2D_3D_fast, MakeONFilter
%
% This is a fast implementation for 1D, 2D and 3D orthogonal wavelet transform
% (c) 2013 Written by Behnood Rasti
% behnood.rasti@gmail.com

if nargin<4 
    options=0;
end

[n1,n2,n3]=size(x);
if n1==1 && n2>1 && n3==1 || options
     wc = FWT_PO_fast(x,L,qmf);
elseif n1>1 && n2>1 && n3==1
     wc = FWT2_PO_fast(x,L,qmf);
elseif n1>1 && n2>1 && n3>1
     wc = FWT3_PO_fast(x,L,qmf);
else
    error('wrong input: input has to be a signal (a row vector) or an image (a two dimentional matrix) or 3d matrix.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


