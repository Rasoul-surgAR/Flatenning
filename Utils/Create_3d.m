% In this version we
% 1 - compute the 3d mesh based on a given step, so we have a downsampled version of mesh.

%function [X,F] = create_3d_ply_V1(ply_path , step)


clc
clear 
close all
% Read .ply file captured by 3d camera and reshape it to image size
ply_path = '/media/rasoul-surgar/830873c3-8201-4a0c-8abe-39d71bdf67d7/Flatenning/Flatenning-main/data/3.ply'
ptCloud = pcread(ply_path);
figure , pcshow(ptCloud)
X = ptCloud.Location(:,1);
Y = ptCloud.Location(:,2);
Z = ptCloud.Location(:,3);

width = 512;
height = 384;

XReshape = reshape(X , [width , height]);
YReshape = reshape(Y , [width , height]);
ZReshape = reshape(Z , [width , height]);

% triangulate it with a given step size

% creat a .obj mesh from that