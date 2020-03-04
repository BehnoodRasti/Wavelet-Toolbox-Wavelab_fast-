function x = IWT2_PO_fast(wc,L,qmf)
% IWT2_PO_fast -- Inverse 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    x = IWT2_PO_fast(wc,L,qmf)
%  Inputs
%    wc    2-d wavelet transform [nx by ny by nz array, nx,ny,nz dyadic]
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    x     2-d signal reconstructed from wc
%
%  Description
%    If wc is the result of a forward 2d wavelet transform, with
%    wc = FWT2_PO_fast(x,L,qmf), then x = IWT2_PO_fast(wc,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FWT2_PO_fast, MakeONFilter
%
% This is a fast implementation for IWT2_PO from Wavelab
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com

[nx,ny] = size(wc);
x = wc;
nx = nx/2^(L-1);
ny = ny/2^(L-1);

for jscal=1:L,
    topx = (nx/2+1):nx; botx = 1:(nx/2); allx = 1:nx;
    topy = (ny/2+1):ny; boty = 1:(ny/2); ally = 1:ny;
     iy=1:ny;
        x(allx,iy) = UpDyadLo2(x(botx,iy)',qmf)'...
            + UpDyadHi2(x(topx,iy)',qmf)';
   
    %---------------------------------------------------------------------
     ix=1:nx;
        x(ix,ally) = UpDyadLo2(x(ix,boty),qmf)  ...
            + UpDyadHi2(x(ix,topy),qmf);
    
    %---------------------------------------------------------------------
    nx = 2*nx;
    ny = 2*ny;
end