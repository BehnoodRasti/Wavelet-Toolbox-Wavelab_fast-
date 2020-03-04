function [wcLLL,wcoeff,L] = FUWT3_A_TROUS_fast(x,L,qmf)
% FUWT3_A_TROUS_fast -- 3-d MRA undecimated wavelet transform (periodized, a trous algorithm)
%  Usage
%     [wcoeff,wcLL,L] = FUWT3_A_TROUS_fast(x,L,qmf)
%  Inputs
%    x     3-d image (nr x nc x nb array)
%    L     level of deomposition
%    qmf   quadrature mirror filter
%  Outputs
%    wcLL    3-d undecimated wavelet transform(LLL-coefficients)
%    wcoeff  3-d undecimated wavelet transform(LLH-LHL-LHH-HLL-HLH-HHL-HHH-coefficients)
%    L       level of deomposition
%  Description
%    A three-dimensional Undecimated Wavelet Transform is computed for the
%    array x.  To reconstruct, use IUWT3_A_TROUS_fast.
%
%  See Also
%    IUWT3_A_TROUS_fast, MakeONFilter
%
% (c) 2011, Written by Behnood Rasti
% behnood.rasti@gmail.com

[nr,nc,nb] = size(x);
wcH=x;wcL=x;wcHH=x;wcHL=x;wcLH=x;wcLL=x;wcHHH=x;wcHHL=x;wcHLH=x;wcHLL=x;
wcLHH=x;wcLHL=x;wcLLH=x;wcLLL=x;
wcoef = cell(1,L);
qmf1=MirrorFilt(qmf);
nx=nr;ny=nc;nz=nb;
for jl=1:1:L,
    iz=1:nb;
    ix=1:nr;
    row = x(ix,1:nc,iz);
    rowper=permute(row,[1,3,2]);
    rowpersh=reshape(rowper,nx*nz,ny);
    wcL(ix,:,iz) = permute(reshape(aconv_2(qmf, rowpersh),nx,nz,ny),[1,3,2]);
    wcH(ix,:,iz) = permute(reshape(iconv_2((qmf1), lshift2(rowpersh)),nx,nz,ny),[1,3,2]);
    %---------------------------------------------------------------------------------------------------
    iz=1:nb;
    iy=1:nc;
    rowH = wcH(1:nr,iy,iz);
    rowper1=permute(rowH,[2,1,3]);
    rowper=permute(rowper1,[1,3,2]);
    rowpersh=reshape(rowper,ny*nz,nx);
    wcHH(:,iy,iz) = permute(permute(reshape(iconv_2((qmf1), lshift2(rowpersh)),ny,nz,nx),[1,3,2]),[2,1,3]);
    wcHL(:,iy,iz) = permute(permute(reshape(aconv_2(qmf, rowpersh),ny,nz,nx),[1,3,2]),[2,1,3]);
    
    rowL = wcL(1:nr,iy,iz);
    rowper1=permute(rowL,[2,1,3]);
    rowper=permute(rowper1,[1,3,2]);
    rowpersh=reshape(rowper,ny*nz,nx);
    wcLH(:,iy,iz) = permute(permute(reshape(iconv_2((qmf1), lshift2(rowpersh)),ny,nz,nx),[1,3,2]),[2,1,3]);
    wcLL(:,iy,iz) = permute(permute(reshape(aconv_2(qmf, rowpersh),ny,nz,nx),[1,3,2]),[2,1,3]);
    
    %---------------------------------------------------------------------------------------------------
    iy=1:nc;
    ix=1:nr;
    rowHH = wcHH(ix,iy,1:nb);
    rowresh=reshape(rowHH,nx*ny,nz);
    wcHHH(ix,iy,:) = reshape(iconv_2((qmf1), lshift2(rowresh)),nx,ny,nz);
    wcHHL(ix,iy,:) = reshape(aconv_2(qmf, rowresh),nx,ny,nz);
    
    rowHL = wcHL(ix,iy,1:nb);
    rowresh=reshape(rowHL,nx*ny,nz);
    wcHLH(ix,iy,:) = reshape(iconv_2((qmf1), lshift2(rowresh)),nx,ny,nz);
    wcHLL(ix,iy,:) = reshape(aconv_2(qmf, rowresh),nx,ny,nz);
    
    
    rowLH = wcLH(ix,iy,1:nb);
    rowresh=reshape(rowLH,nx*ny,nz);
    wcLHH(ix,iy,:) = reshape(iconv_2((qmf1), lshift2(rowresh)),nx,ny,nz);
    wcLHL(ix,iy,:) = reshape(aconv_2(qmf, rowresh),nx,ny,nz);
    
    
    rowLL = wcLL(ix,iy,1:nb);
    rowresh=reshape(rowLL,nx*ny,nz);
    wcLLH(ix,iy,:) = reshape(iconv_2((qmf1), lshift2(rowresh)),nx,ny,nz);
    wcLLL(ix,iy,:) = reshape(aconv_2(qmf, rowresh),nx,ny,nz);
    %---------------------------------------------------------------------------------------------------
    
    A=[wcLLH,wcLHL,wcLHH,wcHLL,wcHLH,wcHHL,wcHHH];
    wcoef(1,L+1-jl)=mat2cell(A,nr,7*nc,nb);
    x=wcLLL;
    qmf=UpSampleN(qmf);
    qmf1=UpSampleN(qmf1);
end
wcoeff=cell2mat(wcoef);