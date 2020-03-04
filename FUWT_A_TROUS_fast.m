function [wcL,wcoeff,L] = FUWT_A_TROUS_fast(x,L,qmf)
% FUWT_A_TROUS_fast -- 1-d MRA undecimated wavelet transform (periodized, a trous algorithm)
%  Usage
%     [wcoeff,wcLL] = FUWT_A_TROUS_fast(x,L,qmf)
%  Inputs
%    x     A vector or a 2-d matrix (nr by nc array)
%    L     level of deomposition
%    qmf   quadrature mirror filter
%  Outputs
%    wcLL    1-d undecimated wavelet transform(L-coefficients)
%    wcoeff  1-d undecimated wavelet transform(H-coefficients)
%    L       level of deomposition
%  Description
%    A one-dimensional Undecimated Wavelet Transform is computed for each row of the
%    array x. To reconstruct, use IUWT_A_TROUS_fast.
%
%  See Also
%    IUWT_A_TROUS_fast, MakeONFilter
%
% (c) 2011, Written by Behnood Rasti
% behnood.rasti@gmail.com

[nr,nc] = size(x);
wcH=x;wcL=x;
wcoef = cell(1,L);
qmf1=MirrorFilt(qmf);
for jl=1:1:L,
    ix=1:nr;
    row = x(ix,1:nc);
    wcL(ix, :) = aconv_2(qmf, row);
    wcH(ix, :) = iconv_2((qmf1),lshift2(row));
    A=wcH;
    wcoef(1,L+1-jl)=mat2cell(A,nr,nc);
    x=wcL;
    qmf=UpSampleN(qmf);
    qmf1=UpSampleN(qmf1);
end
wcoeff=cell2mat(wcoef);