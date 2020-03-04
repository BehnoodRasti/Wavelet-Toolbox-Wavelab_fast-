%-------------------------------------------Contents--------------------------------------------------
% Wavelab_fast is a fast wavelet toolbox for one, two, and three dimensional signals.
% Wavelab_fast contains wavelet and undecimated wavelet transforms.
% Wavelet filters must be selected from wavelab toolbox by using MakeONFilter command or 
% from Rice wavelet toolbox by using daubcqf command. 
% Wavelab_fast is written based on modifying codes from Wavelab and adding some codes
% for higher dimensional signals and implementing undecimated wavelet transform using 
% algorithm a trous. The codes provided are much faster than the ones from wavelab for
% 2D and 3D signals. That has been done by skipping loops on pixels. 
% The toolbox is recommended for applying on large 2D and 3D datasets. 
% A brief description is given here.
%
%---------------------------------------Filter banks Tools-----------------------------------------------
% aconv_2, iconv_2-- Convolution Tools for Two-Scale Transform
% UpSample2 -- Upsampling operator
% rshift2 -- Circular right shift
% lshift2 -- Circular left shift 
% UpDyadHi2 -- High-Pass Upsampling operator; periodized
% UpDyadLo2 -- Low-Pass Upsampling operator; periodized
% DownDyadLo2-- Low-Pass Downsampling operator (periodized)
% DownDyadHi2-- Hi-Pass Downsampling operator (periodized)
%------------------------------------------Orthogonal Wavelets--------------------------------------------
% FWT_PO_fast -- 1-d MRA wavelet transform (periodized, orthogonal)
% IWT_PO_fast -- Inverse 1-d MRA wavelet transform (periodized, orthogonal)
% FWT2_PO_fast -- 2-d MRA wavelet transform (periodized, orthogonal)
% IWT2_PO_fast -- Inverse 2-d MRA wavelet transform (periodized, orthogonal)
% FWT3_PO_fast -- 3-d MRA wavelet transform (periodized, orthogonal)
% IWT3_PO_fast -- Inverse 3-d MRA wavelet transform (periodized, orthogonal)
%-------------------------------------------Undecimated Wavelets-------------------------------------------
% FUWT_A_TROUS_fast -- 1-d MRA undecimated wavelet transform (periodized, a trous algorithm)
% FUWT2_A_TROUS_fast -- 2-d MRA undecimated wavelet transform (periodized, a trous algorithm)
% FUWT3_A_TROUS_fast -- 3-d MRA undecimated wavelet transform (periodized, a trous algorithm)
% IUWT_A_TROUS_fast -- Inverse 1-d MRA undecimated wavelet transform (periodized, a trous Algorithm)
% IUWT2_A_TROUS_fast -- Inverse 2-d MRA undecimated wavelet transform (periodized, a trous Algorithm)
% IUWT3_A_TROUS_fast -- Inverse 3-d MRA undecimated wavelet transform (periodized, a trous Algorithm)
%-------------------------------------------------Using Wavelab_Fast Toolbox-------------------------------------
% The toolbox can be used by the following commands (recommended)
% FWT_PO_1D_2D_3D_fast -- Applies 1-d, 2-d, and 3-d Forward orthogonal Wavelet Transform (periodized)
%                         corresponding to the size of the input data.   
% IWT_PO_1D_2D_3D_fast -- Inverse 1-d, 2-d, and 3-d orthogonal wavelet transform (periodized)
%
% FUWT_1D_2D_3D_A_TROUS_fast -- 1-d, 2-d, and 3-d  Forward Undecimated Wavelet Transform (periodized, a trous Algorithm)
% IUWT_1D_2D_3D_A_TROUS_fast -- 1-d, 2-d, and 3-d Inverse Undecimated Wavelet Transform (periodized, a trous Algorithm)
%---------------------------------------------------------------------------------------------------------------------
% 