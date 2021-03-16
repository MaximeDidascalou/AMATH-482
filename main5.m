clear all; close all; clc;

%%
%ski = VideoReader('monte_carlo_low.mp4');
ski = VideoReader('ski_drop_low.mp4');
video = read(ski);
[n m uint frames] = size(video);
dt = 1/ski.Framerate;
t = 0:dt:ski.Duration;

for j = 1:frames
    %ski_reshape = double(reshape(video(:,:,1,j), n*m, 1));
    ski_reshape = reshape(im2double(rgb2gray(video(:,:,:,j))), n*m, 1);
    v_ski(:,j) = ski_reshape;
end

X1 = v_ski(:,1:end-1);
X2 = v_ski(:,2:end);

%% SVD of X1 and Computation of ~S

[U, Sigma, V] = svd(X1,'econ');


%%
subplot(2,1,1)
plot(diag(Sigma),'k.','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 100])
ylabel("\sigma")
subplot(2,1,2)
semilogy(diag(Sigma),'k.','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 frames-1])
ylabel("\sigma (log scale)")
%%
r = 5;
U = U(:,1:r); 
Sigma = Sigma(1:r,1:r); 
V = V(:,1:r);

S = U'*X2*V*diag(1./diag(Sigma));
[eV, D] = eig(S); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues
omega = log(mu)/dt;
Phi = U*eV;
%%
subplot(1,2,1)
hold on
plot([-30 10], [0 0], 'k')
plot([0 0],[-30 30],'k')
plot(real(omega), imag(omega),'o')
ylabel("Imaginary")
xlabel("Real")
subplot(1,2,2)
hold on
plot([-30 0.5], [0 0], 'k')
plot([0 0],[-30 30],'k')
plot(real(omega), imag(omega),'o')
ylabel("Imaginary")
xlabel("Real")
xlim([-0.5 0.5])
ylim([-0.5 0.5])

%%

thresh = 0.1;
bg = find(abs(omega) < thresh);
omega_bg = omega(bg);
phi_bg = Phi(:,bg);

b = phi_bg\X1(:,1);
t = t(1:end-1);

u_modes = zeros([length(omega_bg) length(t)]);

for j = 1:length(t)
    u_modes(:,j) = b.*exp(omega_bg*t(j)); 
end

X_bg = phi_bg*u_modes;

%%
for j = 1:length(t)
    X_bg_vid(:,:,:,j) = uint8(reshape(X_bg(:,j),[],960)*255);
end

for j=1:length(t)
  frame=X_bg_vid(:,:,:,j);
  imshow(frame); drawnow
end

%% Subtract background to get foreground

X_fg = X1 - abs(X_bg);
R = X_fg - abs(X_fg);
R = R./2;
X_bg = R + abs(X_bg);
X_fg = X_fg - R;

%%

for j = 1:length(t)
    X_fg_vid(:,:,:,j) = uint8(reshape(X_fg(:,j),[],960)*255);
end

%%
for j=1:length(t)
  frame=X_fg_vid(:,:,:,j);
  imshow(frame*1.5); drawnow
end


%%
f1 = 115;
f2 = 290;
f3 = 420;

subplot(3,3,1)
imshow(rgb2gray(video(:,:,:,f1)))
subplot(3,3,2)
imshow(rgb2gray(video(:,:,:,f2)))
subplot(3,3,3)
imshow(rgb2gray(video(:,:,:,f3)))
subplot(3,3,4)
imshow(X_bg_vid(:,:,:,f1))
subplot(3,3,5)
imshow(X_bg_vid(:,:,:,f2))
subplot(3,3,6)
imshow(X_bg_vid(:,:,:,f3))
subplot(3,3,7)
imshow(X_fg_vid(:,:,:,f1)*2)
subplot(3,3,8)
imshow(X_fg_vid(:,:,:,f2)*2)
subplot(3,3,9)
imshow(X_fg_vid(:,:,:,f3)*2)











