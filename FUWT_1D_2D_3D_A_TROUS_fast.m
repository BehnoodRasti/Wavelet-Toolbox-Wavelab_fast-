function [wcL,wc,L]=FUWT_1D_2D_3D_A_TROUS_fast(x,qmf,L,options)

% FUWT_1D_2D_3D_A_TROUS_fast -- Forward Undecimated Wavelet Transform (periodized, a trous Algorithm)
%  Usage
%     [wcL,wc,L] = FUWT_1D_2D_3D_A_TROUS_fast(x,qmf,L,options)
%  Inputs
%    x    1-d or 2-d or 3-d signal;
%    L    decomposition level
%    qmf  quadrature mirror filter (orthonormal)
%    options: 1 applies 1D-FUWT on each row of input matrix (this is fast implementation instead of using loop...    
%             when you have many pixel vectors and you want to apply 1D-UWT)  
%  Outputs
%    wcL   undecimated wavelet transform(LLL-coefficients)
%    wc    undecimated wavelet transform(LLH...HHH-coefficients)
%    L     level of decomposition
%
%  Description
%    1. qmf filter may be obtained from MakeONFilter   
%    2. usually, length(qmf) < 2^(L+1)
%    3. To reconstruct use  IUWT_1D_2D_3D_A_TROUS_fast
%
%  See Also
%    IUWT_1D_2D_3D_A_TROUS_fast, MakeONFilter
%
% (c) 2011, Written by Behnood Rasti
% behnood.rasti@gmail.com

if nargin<4 
    options=0;
end

[n1,n2,n3]=size(x);
if n1==1 && n2>1 && n3==1 || options
     [wcL,wc,L] = FUWT_A_TROUS_fast(x,L,qmf);
elseif n1>1 && n2>1 && n3==1
     [wcL,wc,L] = FUWT2_A_TROUS_fast(x,L,qmf);
elseif n1>1 && n2>1 && n3>1
     [wcL,wc,L] = FUWT3_A_TROUS_fast(x,L,qmf);
else
    error('wrong input--input has to be a signal (a row vector) or an image (a two dimentional matrix) or 3d matrix.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


