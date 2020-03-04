function x = IUWT3_A_TROUS_fast(wcLLL,wcoeff,L,qmf)
% IUWT3_A_TROUS_fast -- Inverse 3-d MRA undecimated wavelet transform (periodized, a trous Algorithm)
%  Usage
%    x = IUWT3_A_TROUS_fast(wcoeff,wcLLL,L,qmf)
%  Inputs
%    wcLLL  undecimated wavelet transform(LLL-coefficients)
%    wcoeff 3-d undecimated wavelet transform(wcLLH,wcLHL,wcLHH,wcHLL,wcHLH,wcHHL,wcHHH-coefficients)
%    L      level of decomposition
%    qmf    quadrature mirror filter
%  Outputs
%    x     3-d signal reconstructed from wcoeff and wcLLL
%
%  Description
%    If wcoeff and wcLL is the result of a forward 3d undecimated wavelet transform, with
%    [wcoeff,wcLL] = FUWT3_A_TROUS_fast(x,L,qmf), then x = IUWT3_A_TROUS_fast(wcoeff,wcLLL,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FUWT3_A_TROUS_fast, MakeONFilter
%
% This is a fast implementation for 3D undecimated wavelet transform
% (c) 2013 Behnood Rasti
% behnood.rasti@gmail.com

[nr,nc,nb] = size(wcLLL);
wcoeff=mat2cell(wcoeff,nr,(nc)*ones(1,7*L),nb);
x = wcLLL;
x1=x;x2=x;x3=x;x4=x;
qmf1=MirrorFilt(qmf);
nx=nr;ny=nc;nz=nb;
for jl=1:1:L,
    qmf2=UpSampleN(qmf,2^(L-jl));
    qmf3=UpSampleN(qmf1,2^(L-jl));
    wcHHH = cell2mat(wcoeff(7*jl));
    wcHHL = cell2mat(wcoeff(7*jl-1));
    wcHLH = cell2mat(wcoeff(7*jl-2));
    wcHLL = cell2mat(wcoeff(7*jl-3));
    wcLHH = cell2mat(wcoeff(7*jl-4));
    wcLHL = cell2mat(wcoeff(7*jl-5));
    wcLLH = cell2mat(wcoeff(7*jl-6));
    %-------------------------------------------------------------------------------------------------
    ix=1:nr;
    iy=1:nc;
    
    x1(ix,iy,1:nb) =  reshape((iconv_2(qmf2, reshape(x(ix,iy,:),nx*ny,nz))  ...
        + aconv_2((qmf3), rshift2(reshape(wcLLH(ix,iy,:),nx*ny,nz)))),nx,ny,nz)/2;
    
    
    x2(ix,iy,1:nb) = reshape((iconv_2(qmf2,reshape(wcLHL(ix,iy,:),nx*ny,nz))  ...
        + aconv_2((qmf3), rshift2(reshape(wcLHH(ix,iy,:),nx*ny,nz)))),nx,ny,nz)/2;
    
    
    x3(ix,iy,1:nb) =  reshape((iconv_2(qmf2, reshape(wcHLL(ix,iy,:),nx*ny,nz))  ...
        + aconv_2((qmf3), rshift2(reshape(wcHLH(ix,iy,:),nx*ny,nz)))),nx,ny,nz)/2;
    
    
    x4(ix,iy,1:nb) =  reshape((iconv_2(qmf2, reshape(wcHHL(ix,iy,:),nx*ny,nz))  ...
        + aconv_2((qmf3), rshift2(reshape(wcHHH(ix,iy,:),nx*ny,nz)))),nx,ny,nz)/2;
    
    %-------------------------------------------------------------------------------------------------
    iz=1:nb;
    iy=1:nc;
    rowper1=permute(x1(:,iy,iz),[2,1,3]);
    rowper2=permute(rowper1,[1,3,2]);
    rowper3=permute(x2(:,iy,iz),[2,1,3]);
    rowper4=permute(rowper3,[1,3,2]);
    
    
    x1(1:nr,iy,iz) = permute(permute(reshape( (iconv_2(qmf2, reshape(rowper2,ny*nz,nx))  ...
        + aconv_2((qmf3), rshift2(reshape(rowper4,ny*nz,nx)))),ny,nz,nx),[1,3,2]),[2,1,3])/2;
    
    
    rowper1=permute(x3(:,iy,iz),[2,1,3]);
    rowper2=permute(rowper1,[1,3,2]);
    rowper3=permute(x4(:,iy,iz),[2,1,3]);
    rowper4=permute(rowper3,[1,3,2]);
    
    x2(1:nr,iy,iz) =  permute(permute(reshape((iconv_2(qmf2,reshape(rowper2,ny*nz,nx))  ...
        + aconv_2((qmf3), rshift2(reshape(rowper4,ny*nz,nx)))),ny,nz,nx),[1,3,2]),[2,1,3])/2;
    
    %-------------------------------------------------------------------------------------------------
    iz=1:nb;
    ix=1:nr;
    rowper1=permute(x1(ix,:,iz),[1,3,2]);
    rowper2=permute(x2(ix,:,iz),[1,3,2]);
    
    x(ix,1:nc,iz) = permute(reshape((iconv_2(qmf2, reshape(rowper1,nx*nz,ny))  ...
        + aconv_2((qmf3), rshift2(reshape(rowper2,nx*nz,ny)))),nx,nz,ny),[1,3,2])/2;
    
end