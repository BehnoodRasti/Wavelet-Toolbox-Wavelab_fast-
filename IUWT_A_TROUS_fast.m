function x = IUWT_A_TROUS_fast(wcL,wcoeff,L,qmf)
% IUWT_A_TROUS_fast -- Inverse 1-d MRA undecimated wavelet transform (periodized, a trous Algorithm)
%  Usage
%    x = IUWT_A_TROUS_fast(wcL,wcoeff,L,qmf)
%  Inputs
%    wcLL  undecimated wavelet transform(L-coefficients)
%    wcoeff  1-d undecimated wavelet transform(H-coefficients)
%    L     level of decomposition
%    qmf   quadrature mirror filter
%  Outputs
%    x     2-d signal reconstructed from wcoeff and wcLL
%
%  Description
%    If wcoeff and wcL is the result of a forward 1-d undecimated wavelet transform, with
%    [wcoeff,wcL] = FUWT_A_TROUS_fast(x,L,qmf), then x = IUWT_A_TROUS_fast(wcL,wcoeff,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FUWT_A_TROUS_fast, MakeONFilter
%
% This is a fast implementation for 2D undecimated wavelet transform
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com

[nr,nc] = size(wcL);
wcoeff=mat2cell(wcoeff,nr,(nc)*ones(1,L));
x = wcL;
%x1=x;x2=x;
qmf1=MirrorFilt(qmf);
for jl=1:1:L,
    qmf2=UpSampleN(qmf,2^(L-jl));
    qmf3=UpSampleN(qmf1,2^(L-jl));
    wcH = cell2mat(wcoeff(jl));
    ix=1:nr;
    x(ix,1:nc) = (iconv_2(qmf2, x(ix,:))  ...
        + aconv_2((qmf3), rshift2(wcH(ix,:))))/2;
    
end