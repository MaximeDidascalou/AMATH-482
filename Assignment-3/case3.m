clear all; close all; clc;


% -- Load Movie
load cam1_3.mat;
load cam2_3.mat;
load cam3_3.mat;
vidFrames1_3 = vidFrames1_3(:,:,:,14:end);
vidFrames3_3 = vidFrames3_3(:,:,:,5:end);
[height1 width1 rgb1 num_frames1] = size(vidFrames1_3);
[height2 width2 rgb2 num_frames2] = size(vidFrames2_3);
[height3 width3 rgb3 num_frames3] = size(vidFrames3_3);
num_frames = min([num_frames1 num_frames2 num_frames3]);


%%

% -- Watch Movie 1
for j=1:num_frames
  X=vidFrames1_3(:,:,:,j);
  X=rgb2gray(X);
  imshow(X); drawnow
end

%%

% -- Watch Movie 2
for j=1:num_frames
  X=vidFrames2_3(:,:,:,j);
  X=rgb2gray(X);
  imshow(X); drawnow
end

%%
% -- Watch Movie 3
for j=1:num_frames
  X=vidFrames3_3(:,:,:,j);
  X=rgb2gray(X);
  imshow(X); drawnow
end

%%

% -- Get coordinates for 1st vid
coords1 = zeros([2 num_frames]);
for j=1:num_frames
    X=rgb2gray(vidFrames1_3(:,:,:,j));
    X=im2double(X);
    X_temp=X(200:end,285:380);
    [M,I] = max(X_temp(:));
    dims = size(X_temp)
    [y,x] = ind2sub([dims(1) dims(2)],I);
%     y = y+191;
%     x = x+299;
    coords1(:,j) = [x y];
end
coords1(1,:) = coords1(1,:) - mean(coords1(1,:));
coords1(2,:) = coords1(2,:) - mean(coords1(2,:));
figure(1)
plot(1:num_frames,coords1(2,:))
plot3(coords1(2,:),coords1(1,:),1:num_frames)
xlim([-100,100])
ylim([-100,100])
xlabel('Y')
ylabel('X')
zlabel('time')
grid on

%%

% -- Get coordinates for second vid
coords2 = zeros([2 num_frames]);
for j=1:num_frames
    X=rgb2gray(vidFrames2_3(:,:,:,j));
    X=im2double(X);
    X_temp=X(160:end,100:400);
    [M,I] = max(X_temp(:));
    dims = size(X_temp)
    [y,x] = ind2sub([dims(1) dims(2)],I);
%     y = y+191;
%     x = x+255;
    coords2(:,j) = [x y];
end
coords2(1,:) = coords2(1,:) - mean(coords2(1,:));
coords2(2,:) = coords2(2,:) - mean(coords2(2,:));
figure(2)
plot(1:num_frames,coords2(2,:))
plot3(coords2(2,:),coords2(1,:),1:num_frames)
xlim([-150,150])
ylim([-150,150])
xlabel('Y')
ylabel('X')
zlabel('time')
grid on


%%

% -- Get coordinates for third vid
coords3 = zeros([2 num_frames]);
for j=1:num_frames
    X=rgb2gray(vidFrames3_3(:,:,:,j));
    X=im2double(X);
    X_temp=X(160:320,213:512);
    [M,I] = max(X_temp(:));
    dims = size(X_temp)
    [y,x] = ind2sub([dims(1) dims(2)],I);
%     y = y+191;
%     x = x+255;
    coords3(:,j) = [x y];
end
coords3(1,:) = coords3(1,:) - mean(coords3(1,:));
coords3(2,:) = coords3(2,:) - mean(coords3(2,:));
figure(3)
plot(1:num_frames,coords3(1,:))
plot3(coords3(2,:),coords3(1,:),1:num_frames)
xlim([-100,100])
ylim([-100,100])
xlabel('Y')
ylabel('X')
zlabel('time')
grid on


%%
A = zeros([6 num_frames]);
A(1:2,:) = coords1;
A(3:4,:) = coords2;
A(5:6,:) = coords3;
% A = A';

[U,S,V]=svd(A/sqrt(num_frames-1),'econ'); % perform the SVD
lambda=diag(S).^2; % produce diagonal variances
Y=U'*A; % produce the principal components projection




subplot(1,2,1)
hold on
plot(1:num_frames,Y(1,:))
plot(1:num_frames,Y(2,:))
plot(1:num_frames,Y(3,:))
xlabel('time')
ylabel('Principal Components')
legend('PC1','PC2','PC3')
pbaspect([1 1 1])
subplot(1,2,2)
plot(1:6,diag(S),'ko','Linewidth', 2)
ylabel('Singular Values')
xlim([1,6])
pbaspect([1 1 1])


% plot(1:num_frames,Y(3,:))
% plot(1:num_frames,Y(:,4))
% plot(1:num_frames,Y(:,5))
% plot(1:num_frames,Y(:,6))