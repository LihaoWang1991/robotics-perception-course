function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly


n=size(x,1);
A=[];
for i=1:n
    xci=inv(K)*[x(i,:),1]';     %extend x(i)=[x y] to [x y 1]  
    u=xci(1);
    v=xci(2);
    Xi=[X(i,:),1]';            %extend X(i)=[X Y Z] to [X Y Z 1]
    Ai=[zeros(1,4), -Xi', v*Xi';Xi' zeros(1,4) -u*Xi';-v*Xi' u*Xi' zeros(1,4)];
    A(3*i-2:3*i,1:12)=Ai;
end

% Calculate P=[R t]
[UA,DA,VA]=svd(A);
P_vector=VA(:,end);     
P1=P_vector(1:4,1);
P2=P_vector(5:8,1);
P3=P_vector(9:12,1);
P=[P1';P2';P3'];

% Calculate R and t
R=P(:,1:3);
t=P(:,4);

% Clean up R and t
[UR,DR,VR]=svd(R);
if det(UR*VR')>0
    R_CleanUp=UR*VR';
    t_CleanUp=t/DR(1,1);
else
    R_CleanUp=-UR*VR';
    t_CleanUp=-t/DR(1,1);
end

% Calculated t is camera 3 translation vector with respect to camera 3
% as described in the week 4 task pdf, the desired output C is the translation of camera 3
% with respect to camera 1, therefore C=-R'*t
R=R_CleanUp;
C=-R_CleanUp'*t_CleanUp;


 

