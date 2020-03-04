function x = IWT3_PO_fast(wc,L,qmf)
% IWT3_PO_fast -- Inverse 3-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    x = IWT3_PO_fast(wc,L,qmf)
%  Inputs
%    wc    3-d wavelet transform [nx by ny by nz array, nx,ny,nz dyadic]
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    x     3-d signal reconstructed from wc
%
%  Description
%    If wc is the result of a forward 3d wavelet transform, with
%    wc = FWT3_PO_fast(x,L,qmf), then x = IWT3_PO_fast(wc,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FWT3_PO_fast, MakeONFilter
%
% This is a fast implementation for IWT3_PO from Wavelab
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com
[nx,ny,nz] = size(wc);
x = wc;
%nc = 2^(L+1);
nx = nx/2^(L-1);
ny = ny/2^(L-1);
nz = nz/2^(L-1);

for jscal=1:L,
    topx = (nx/2+1):nx; botx = 1:(nx/2); allx = 1:nx;
    topy = (ny/2+1):ny; boty = 1:(ny/2); ally = 1:ny;
    topz = (nz/2+1):nz; botz = 1:(nz/2); allz = 1:nz;
    
    iy=1:ny;
    iz=1:nz;
    rowper1=permute(x(botx,iy,iz),[2,1,3]);
    rowper2=permute(rowper1,[1,3,2]);
    rowper3=permute(x(topx,iy,iz),[2,1,3]);
    rowper4=permute(rowper3,[1,3,2]);
    x(allx,iy,iz) =  permute(permute(reshape( UpDyadLo2(reshape(rowper2,ny*nz,nx/2),qmf)...
        + UpDyadHi2(reshape(rowper4,ny*nz,nx/2),qmf),ny,nz,nx),[1,3,2]),[2,1,3]);
    
    %---------------------------------------------------------------------
    ix=1:nx;
    iy=1:ny;
    x(ix,iy,allz) = reshape(UpDyadLo2(reshape(x(ix,iy,botz),nx*ny,nz/2),qmf)  ...
        + UpDyadHi2(reshape(x(ix,iy,topz),nx*ny,nz/2),qmf),nx,ny,nz);

    %----------------------------------------------------------------------

    ix=1:nx;
    iz=1:nz;
    rowper1=permute(x(ix,boty,iz),[1,3,2]);
    rowper2=permute(x(ix,topy,iz),[1,3,2]);
    x(ix,ally,iz) = permute(reshape(UpDyadLo2(reshape(rowper1,nx*nz,ny/2),qmf)  ...
        + UpDyadHi2(reshape(rowper2,nx*nz,ny/2),qmf),nx,nz,ny),[1,3,2]);
    
    %-------------------------------------------------------------------------
    nx = 2*nx;
    ny = 2*ny;
    nz = 2*nz;
end



