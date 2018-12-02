function E = EssentialMatrixFromFundamentalMatrix(F,K)
%% EssentialMatrixFromFundamentalMatrix
% Use the camera calibration matrix to esimate the Essential matrix
% Inputs:
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     F - size (3 x 3) fundamental matrix from EstimateFundamentalMatrix
% Outputs:
%     E - size (3 x 3) Essential matrix with singular values (1,1,0)

E=K'*F*K;
[U,D,V]=svd(E);
E_CleanUp=U*diag([1,1,0])*V';

E_Normlized=E_CleanUp/norm(E_CleanUp);

E=E_Normlized;  % Because the variable name in the definition of this function is F, not F_Normlized


