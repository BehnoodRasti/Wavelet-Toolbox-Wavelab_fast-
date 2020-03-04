function wc = FWT2_PO_fast(x,L,qmf)
% FWT2_PO_fast -- 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    wc = FWT2_PO_fast(x,L,qmf)
%  Inputs
%    x     2-d object (nx by ny by nz array, nx,ny,nz dyadic)
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    wc    2-d wavelet transform
%
%  Description
%    A three-dimensional Wavelet Transform is computed for the
%    array x.  To reconstruct, use IWT2_PO_fast.
%
%  See Also
%    IWT2_PO_fast, MakeONFilter
%
% This is a fast implementation for FWT2_PO from Wavelab
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com
[nx,ny] = size(x);
wc = x;
for jscal=1:1:L,
    topx = (nx/2+1):nx;    botx = 1:(nx/2);
    topy = (ny/2+1):ny;    boty = 1:(ny/2);
    %----------------------------------------
    ix=1:nx;
            row = wc(ix,1:ny);
            wc(ix,topy) = DownDyadHi2(row,qmf);
            wc(ix,boty) = DownDyadLo2(row,qmf);
    
    %----------------------------------------
    iy=1:ny;
            row = wc(1:nx,iy)';
            wc(topx,iy) = DownDyadHi2(row,qmf)';
            wc(botx,iy) = DownDyadLo2(row,qmf)';
    
    %-----------------------------------------
    nx = nx/2;
    ny = ny/2;
end