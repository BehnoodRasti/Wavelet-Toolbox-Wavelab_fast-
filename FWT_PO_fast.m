function wc = FWT_PO_fast(x,L,qmf)
% FWT_PO_fast -- 1-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    wc = FWT_PO_fast(x,L,qmf)
%  Inputs
%    x     1-d or 2-d matrix (nx by ny array, ny dyadic)
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    wc    1-d wavelet transform applied on row of x
%
%  Description
%    A one-dimensional Wavelet Transform is computed for the
%    array x.  To reconstruct, use IWT_PO_fast.
%
%  See Also
%    IWT_PO_fast, MakeONFilter
%
% This is a fast implementation for 1D orthogonal wavelet transform when
% there is more than one vector (signal) to be applied on. FWT_PO_fast applies 
% forward 1d wavelet transform on row vectors.
% (c) 2013 Written by Behnood Rasti
% behnood.rasti@gmail.com
  [nr,nc] = size(x) ;
  wc = x ;
  beta = x;  %take samples at finest scale as beta-coeffts % applying 1D on each row of the matrix
  for j=1:1:L
       ix=1:nr;
       topx = (nc/2+1):nc;
       alfa = DownDyadHi2( beta, qmf);
       wc(ix,topx) = alfa;
       beta = DownDyadLo2(beta, qmf);
       nc=nc/2;
  end
  wc(ix,1:nc)= beta;
