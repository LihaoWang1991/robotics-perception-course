function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE
A=[];
for i=1:4
    x1=video_pts(i,1);
    x2=video_pts(i,2);
    x1prime=logo_pts(i,1);
    x2prime=logo_pts(i,2);
    A(2*i-1,:)=[-x1,-x2,-1,0,0,0,x1*x1prime,x2*x1prime,x1prime];
    A(2*i,:)=[0,0,0,-x1,-x2,-1,x1*x2prime,x2*x2prime,x2prime];
end

[U, S, V] = svd(A);

H(1,:)= V(1:3,9);
H(2,:)= V(4:6,9);
H(3,:)= V(7:9,9);
end

