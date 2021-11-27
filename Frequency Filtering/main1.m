clear all; close all; clc

load subdata.mat

L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); 
x = x2(1:n); 
y = x; 
z = x;

k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; 
ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

for j=1:49
    U_noisy(:,:,:,j)=reshape(subdata(:,j),n,n,n);
    M = max(abs(U_noisy),[],'all');
%     close all, isosurface(X,Y,Z,abs(U_noisy)/M,0.5)
%     axis([-20 20 -20 20 -20 20]), grid on, drawnow
%     pause(0.1)
end


%averaging in frequency space 
U_average_freq = zeros(n,n,n);
for j=1:49
    U_freq_noisy(:,:,:,j) = fftshift(fftn(U_noisy(:,:,:,j)));
    U_average_freq = U_average_freq + U_freq_noisy(:,:,:,j);
end
U_average_freq = U_average_freq/49;
   

% finding the highest frequency 
[m,ind] = max(U_average_freq(:));
[ind_kx,ind_ky,ind_kz] = ind2sub([n,n,n],ind);
center_Kx = Kx(ind_kx,ind_ky,ind_kz);
center_Ky = Ky(ind_kx,ind_ky,ind_kz);
center_Kz = Kz(ind_kx,ind_ky,ind_kz);


%setting up the filter
tau = 0.3;
filter = exp(-tau.*((Kx - center_Kx).^2+(Ky - center_Ky).^2 + (Kz - center_Kz).^2));

coordinates = zeros([49 3]);

for j=1:49
   % filtering the signal
    U_freq_filtered = ifftshift(U_freq_noisy(:,:,:,j).*filter);
    U_filtered = ifftn(U_freq_filtered);
    
    % finding the x and y coordinates for the plane to follow 
    [m,ind] = max(U_filtered(:));
    [ind_x,ind_y,ind_z] = ind2sub([n,n,n],ind);
    coordinates(j,1) = X(ind_x,ind_y,ind_z);
    coordinates(j,2) = Y(ind_x,ind_y,ind_z);
    coordinates(j,3) = Z(ind_x,ind_y,ind_z);
    
    % plotting
%     M = max(abs(U_filtered),[],'all');
%     close all, isosurface(X,Y,Z,abs(U_filtered)/M,0.5)
%     axis([-20 20 -20 20 -20 20]), grid on, drawnow
%     pause(1)
end
plot3(coordinates(:,1), coordinates(:,2),coordinates(:,3)) %plotting path
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
stringstart = ['start (',num2str(coordinates(1,1)),',',num2str(coordinates(1,2)),',',num2str(coordinates(1,3)),')']
stringend = ['end (',num2str(coordinates(end,1)),',',num2str(coordinates(end,2)),',',num2str(coordinates(end,3)),')']
text(coordinates(1,1), coordinates(1,2),coordinates(1,3),stringstart)
text(coordinates(end,1), coordinates(end,2),coordinates(end,3),stringend)
grid on
print -deps epsFig



















