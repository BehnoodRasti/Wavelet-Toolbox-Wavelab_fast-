function x = IWT_PO_fast(wc,L,qmf)
% IWT_PO_fast -- Inverse 1-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    x = IWT_PO_fast(wc,L,qmf)
%  Inputs
%    wc    1-d inverse wavelet transform of each row of wc [nx by ny array, ny dyadic]
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    x     1-d signal reconstructed from wc
%
%  Description
%    If wc is the result of a forward 1d wavelet transform, with
%    wc = FWT_PO_fast(x,L,qmf), then x = IWT_PO_fast(wc,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FWT_PO_fast, MakeONFilter
%
% This is a fast implementation for 1D orthogonal wavelet transform when
% there is more than one vector (signal) to be applied on. IWT_PO_fast applies 
% inverse 1d wavelet transform on row vectors.
% (c) 2013 Written by Behnood Rasti
% behnood.rasti@gmail.com
[nr,nc] = size(wc);
x = wc;
nc = nc/2^(L-1);

for jscal=1:L,
    ix=1:nr;
    topx = (nc/2+1):nc; botx = 1:(nc/2); allx = 1:nc;
    %---------------------------------------------------------------------
        x(ix,allx) = UpDyadLo2(x(ix,botx),qmf)  ...
            + UpDyadHi2(x(ix,topx),qmf);
    %---------------------------------------------------------------------
    nc = 2*nc;
end