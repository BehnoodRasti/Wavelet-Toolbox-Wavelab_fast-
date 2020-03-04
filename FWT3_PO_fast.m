function wc = FWT3_PO_fast(x,L,qmf)
% FWT3_PO_fast -- 3-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    wc = FWT3_PO_fast(x,L,qmf)
%  Inputs
%    x     3-d object (nx by ny by nz array, nx,ny,nz dyadic)
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    wc    3-d wavelet transform
%
%  Description
%    A three-dimensional Wavelet Transform is computed for the
%    array x.  To reconstruct, use IWT3_PO_fast.
%
%  See Also
%    IWT3_PO_fast, MakeONFilter
%
% This is a fast implementation for FWT3_PO from Wavelab
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com


[nx,ny,nz] = size(x);
wc = x;
for jscal=1:1:L,
    topx = (nx/2+1):nx;    botx = 1:(nx/2);
    topy = (ny/2+1):ny;    boty = 1:(ny/2);
    topz = (nz/2+1):nz;    botz = 1:(nz/2);

     ix=1:nx;
         iy = 1:ny;
         ixy=1:nx*ny;
            row = wc(ix,iy,1:nz);
            wc1=zeros(nx*ny,nz);
            rowresh=reshape(row,nx*ny,nz);
            wc1(ixy,topz) = DownDyadHi2(rowresh,qmf);
            wc1(ixy,botz) = DownDyadLo2(rowresh,qmf);
            wc(ix,iy,1:nz)= reshape(wc1(ixy,[botz,topz]),nx,ny,nz);

    %----------------------------------------
         ix=1:nx;
         iz = 1:nz;
         ixz=1:nx*nz;
            row = wc(ix,1:ny,iz);
            rowper=permute(row,[1,3,2]);
            wc1=zeros(nx*nz,ny);
            rowpersh=reshape(rowper,nx*nz,ny);
            wc1(ixz,topy) = DownDyadHi2(rowpersh,qmf);
            wc1(ixz,boty) = DownDyadLo2(rowpersh,qmf);
            wc(ix,1:ny,iz)= permute(reshape(wc1(ixz,[boty,topy]),nx,nz,ny),[1,3,2]);

    %----------------------------------------
     iy=1:ny;
         iz = 1:nz;
         iyz=1:ny*nz;
            row = wc(1:nx,iy,iz);
            rowper1=permute(row,[2,1,3]);
            rowper=permute(rowper1,[1,3,2]);
            wc1=zeros(nz*ny,nx);
            rowpersh=reshape(rowper,ny*nz,nx);
            wc1(iyz,topx) = DownDyadHi2(rowpersh,qmf);
            wc1(iyz,botx) = DownDyadLo2(rowpersh,qmf);
            wc(1:nx,iy,iz)= permute(permute(reshape(wc1(iyz,[botx,topx]),ny,nz,nx),[1,3,2]),[2,1,3]);

    %-----------------------------------------
    nx = nx/2;
    ny = ny/2;
    nz = nz/2;
end

