function x = IUWT2_A_TROUS_fast(wcLL,wcoeff,L,qmf)
% IUWT2_A_TROUS_fast -- Inverse 2-d MRA undecimated wavelet transform (periodized, a trous Algorithm)
%  Usage
%    x = IUWT2_A_TROUS_fast(wcoeff,wcLL,L,qmf)
%  Inputs
%    wcLL  undecimated wavelet transform(LL-coefficients)
%    wcoeff  2-d undecimated wavelet transform(LH-HL-HH-coefficients)
%    L     level of decomposition
%    qmf   quadrature mirror filter
%  Outputs
%    x     2-d signal reconstructed from wcoeff and wcLL
%
%  Description
%    If wcoeff and wcLL is the result of a forward 2d undecimated wavelet transform, with
%    [wcoeff,wcLL] = FUWT2_A_TROUS_fast(x,L,qmf), then x = IUWT2_A_TROUS_fast(wcoeff,wcLL,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FUWT2_A_TROUS_fast, MakeONFilter
%
% This is a fast implementation for 2D undecimated wavelet transform
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com

[nr,nc] = size(wcLL);
wcoeff=mat2cell(wcoeff,nr,(nc)*ones(1,3*L));
x = wcLL;
x1=x;x2=x;
qmf1=MirrorFilt(qmf);
for jl=1:1:L,
    qmf2=UpSampleN(qmf,2^(L-jl));
    qmf3=UpSampleN(qmf1,2^(L-jl));
    wcHH = cell2mat(wcoeff(3*jl));
    wcHL = cell2mat(wcoeff(3*jl-1));
    wcLH = cell2mat(wcoeff(3*jl-2));
    iy=1:nc;
    x1(1:nr,iy) =  (iconv_2(qmf2, x(:,iy)')'  ...
        + aconv_2((qmf3), rshift2(wcLH(:,iy)'))')/2;
    x2(1:nr,iy) =  (iconv_2(qmf2, wcHL(:,iy)')'  ...
        + aconv_2((qmf3), rshift2(wcHH(:,iy)'))')/2;
    
    ix=1:nr;
    x(ix,1:nc) = (iconv_2(qmf2, x1(ix,:))  ...
        + aconv_2((qmf3), rshift2(x2(ix,:))))/2;
    
end