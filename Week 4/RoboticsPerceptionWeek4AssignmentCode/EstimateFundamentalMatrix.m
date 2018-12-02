function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2
A=[];
for i=1:size(x1,1)
    A(i,:)=[x1(i,1)*x2(i,1) x1(i,1)*x2(i,2) x1(i,1) x1(i,2)*x2(i,1) x1(i,2)*x2(i,2) x1(i,2) x2(i,1) x2(i,2) 1];
end
[U1,D1,V1]=svd(A);

f=V1(:,9);
F=[];
F(:,1)=f(1:3);
F(:,2)=f(4:6);
F(:,3)=f(7:9);

[U2,D2,V2]=svd(F);
D_CleanUp=zeros(3,3);
D_CleanUp(:,1:2)=D2(:,1:2);

F_CleanUp=U2*D_CleanUp*V2';

F_Normlized=F_CleanUp/norm(F_CleanUp);

F=F_Normlized;   % Because the variable name in the definition of this function is F, not F_Normlized




