function [wcLL,wcoeff,L] = FUWT2_A_TROUS_fast(x,L,qmf)
% FUWT2_A_TROUS_fast -- 2-d MRA undecimated wavelet transform (periodized, a trous algorithm)
%  Usage
%     [wcoeff,wcLL] = FUWT2_A_TROUS_fast(x,L,qmf)
%  Inputs
%    x     2-d image (n by n array)
%    L     level of deomposition
%    qmf   quadrature mirror filter
%  Outputs
%    wcLL    2-d undecimated wavelet transform(LL-coefficients)
%    wcoeff  2-d undecimated wavelet transform(LH-HL-HH-coefficients)
%    L       level of deomposition
%  Description
%    A two-dimensional undecimated Wavelet Transform is computed for the
%    array x.  To reconstruct, use IUWT2_A_TROUS_fast.
%
%  See Also
%    IUWT2_A_TROUS_fast, MakeONFilter
%
% (c) 2011, Written by Behnood Rasti
% behnood.rasti@gmail.com

[nr,nc] = size(x);
wcH=x;wcL=x;wcHH=x;wcHL=x;wcLH=x;wcLL=x;
wcoef = cell(1,L);
qmf1=MirrorFilt(qmf);
for jl=1:1:L,
    ix=1:nr;
    row = x(ix,1:nc);
    wcL(ix, :) = aconv_2(qmf, row);
    wcH(ix, :) = iconv_2((qmf1),lshift2(row));
    
    iy=1:nc;
    rowH = wcH(1:nr,iy)';
    wcHH(:, iy) = iconv_2((qmf1), lshift2(rowH))';
    wcHL(:, iy) = aconv_2(qmf, rowH)';
    rowL = wcL(1:nr,iy)';
    wcLH(:, iy) = iconv_2((qmf1), lshift2(rowL))';
    wcLL(:, iy) = aconv_2(qmf, rowL)';
    
    A=[wcLH,wcHL,wcHH];
    wcoef(1,L+1-jl)=mat2cell(A,nr,3*nc);
    x=wcLL;
    qmf=UpSampleN(qmf);
    qmf1=UpSampleN(qmf1);
end
wcoeff=cell2mat(wcoef);