function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 



for i=1:size(x1,1)
    X(i,:)=Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), X0(i,:));

end

end


function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
for i=1:3
J=zeros(6,3);
J(1:2,:)=Jacobian_Triangulation(C1, R1, K, X0);
J(3:4,:)=Jacobian_Triangulation(C2, R2, K, X0);
J(5:6,:)=Jacobian_Triangulation(C3, R3, K, X0);
b=[x1(1),x1(2),x2(1),x2(2),x3(1),x3(2)]';
x1_reprojected=K*R1*(X0'-C1);
x2_reprojected=K*R2*(X0'-C2);
x3_reprojected=K*R3*(X0'-C3);
u1=x1_reprojected(1);
v1=x1_reprojected(2);
w1=x1_reprojected(3);
u2=x2_reprojected(1);
v2=x2_reprojected(2);
w2=x2_reprojected(3);
u3=x3_reprojected(1);
v3=x3_reprojected(2);
w3=x3_reprojected(3);
f_X=[u1/w1 v1/w1 u2/w2 v2/w2 u3/w3 v3/w3]';
delta_X=inv(J'*J)*J'*(b-f_X);
X0=X0+delta_X';
end
X=X0;

end



function J = Jacobian_Triangulation(C, R, K, X)
x_reprojected=K*R*(X'-C);
u=x_reprojected(1);
v=x_reprojected(2);
w=x_reprojected(3);
f=K(1,1);
px=K(1,3);
py=K(2,3);
r11=R(1,1);
r12=R(1,2);
r13=R(1,3);
r21=R(2,1);
r22=R(2,2);
r23=R(2,3);
r31=R(3,1);
r32=R(3,2);
r33=R(3,3);
du_dX=[f*r11+px*r31 f*r12+px*r32 f*r13+px*r33];
dv_dX=[f*r21+py*r31 f*r22+py*r32 f*r23+py*r33];
dw_dX=[r31 r32 r33];
J=[(w*du_dX-u*dw_dX)/w^2;(w*dv_dX-v*dw_dX)/w^2];
end



