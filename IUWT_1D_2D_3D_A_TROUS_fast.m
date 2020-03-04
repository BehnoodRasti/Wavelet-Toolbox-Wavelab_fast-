function [x,L]=IUWT_1D_2D_3D_A_TROUS_fast(wcL,wc,qmf,L,options)

% IUWT_1D_2D_3D_A_TROUS_fast -- Inverse MRA undecimated wavelet transform (periodized, a trous Algorithm)
%  Usage
%     [x,L]=IUWT_1D_2D_3D_A_TROUS_fast(wcL,wc,qmf,L,options)
%
%  Inputs
%    wcL   undecimated wavelet transform(LLL-coefficients)
%    wc    undecimated wavelet transform(LLH...HHH-coefficients)
%    L     level of decomposition
%    qmf   quadrature mirror filter
%    options: 1 applies 1D-IUWT on each row of matrices wcL and wc (this is fast implementation instead of using loop...
%             when you have many pixel vectors and you want to apply 1D-UWT)    
%  Outputs
%    x     signal reconstructed from wc and wcL
%
%  Description
%    If wcL and wcLL is the result of the forward undecimated wavelet transform, with
%    [wcL,wc]=FUWT_1D_2D_3D_A_TROUS_fast(x,qmf,L), then x=IUWT_1D_2D_3D_A_TROUS_fast(wcL,wc,qmf,L)
%    reconstructs x exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FUWT_1D_2D_3D_A_TROUS_fast, MakeONFilter
%
% (c) 2011, Written by Behnood Rasti
% behnood.rasti@gmail.com

if nargin<5 
    options=0;
end

[n1,n2,n3]=size(wcL);
if n1==1 && n2>1 && n3==1 || options
    x = IUWT_A_TROUS_fast(wcL,wc,L,qmf);
elseif n1>1 && n2>1 && n3==1
    x = IUWT2_A_TROUS_fast(wcL,wc,L,qmf);
elseif n1>1 && n2>1 && n3>1
    x = IUWT3_A_TROUS_fast(wcL,wc,L,qmf);
else
    error('wrong input--input has to be a signal (a row vector) or an image (a two dimentional matrix) or 3d matrix.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

