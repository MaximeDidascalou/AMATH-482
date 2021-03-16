clear all; close all; clc

%%
figure(1)
[y_guns, Fs_guns] = audioread('GNR.m4a');
y_guns = y_guns(1:floor(length(y_guns)/4))';
%y_guns = y_guns';
n = length(y_guns);
tr_guns = n/Fs_guns; % record time in seconds
plot((1:n)/Fs_guns,y_guns);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Sweet Child O Mine');

t2 = linspace(0,tr_guns,n+1);
t = t2(1:n);
k = (1/tr_guns)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);
% p8 = audioplayer(y_guns,Fs_guns); playblocking(p8);

figure(2)
[y_floyd, Fs_floyd] = audioread('Floyd.m4a');
%y_floyd = y_floyd';
y_floyd = y_floyd(1:floor(length(y_floyd)/4))';
n2 = length(y_floyd);
tr_floyd = n2/Fs_floyd; % record time in seconds
plot((1:n2)/Fs_floyd,y_floyd);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Comfortably Numb');

t2_f = linspace(0,tr_floyd,n2+1);
t_f = t2_f(1:n2);
ks_f = (1/tr_floyd)*[-n2/2:n2/2-1];
k_f = ifftshift(ks_f);
% p8_f = audioplayer(y_floyd,Fs_floyd); playblocking(p8_f);


%%
% figure(3)
% a = 1;
% tau = 0:1:tr_guns;
% 
% for j = 1:length(tau)
%    g = exp(-a*(t - tau(j)).^2); % Window function
%    Sg = g.*y_guns;
%    Sgt = fft(Sg);
%    
%    ylm = max(abs(y_guns));
%    subplot(3,1,1) % Time domain
%     plot(t,y_guns,'k','Linewidth',1)
%     hold on
%     plot(t,g,'m','Linewidth',1)
%     set(gca,'Fontsize',16), xlabel('time (t)'), ylabel('S(t)')
% 
%     subplot(3,1,2) % Time domain
%     plot(t,Sg,'k','Linewidth',1)
%     set(gca,'Fontsize',16,'ylim',[-1 1]), xlabel('time (t)'), ylabel('S(t)*g(t-\tau)')
%     ylim([-ylm ylm])
%     
%     subplot(3,1,3) % Fourier domain
%     plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)),'r','Linewidth',2); 
%     set(gca,'Fontsize',16), xlabel('frequency (k)'), ylabel('fft(S(t)*g(t-\tau))')
%     xlim([0,3500])
%     drawnow
%      pause(0.1)
%     clf
% end

%%

figure(3)
a = 0.2;
tau = tr_guns/8:tr_guns/4:tr_guns

for j = 1:length(tau)
   Sg = zeros([1 n]);
   Sg((floor(-n/4+j*n/4)+1):floor(j*n/4)) = y_guns((floor(-n/4+j*n/4)+1):floor(j*n/4));
   Sgt = fft(Sg);
   Sgt_spec(:,j) = fftshift(abs(Sgt)); % We don't want to scale it
   
    subplot(2,2,j) % Fourier domain
    plot(ks,abs(Sgt_spec(:,j))/max(abs(Sgt)),'r','Linewidth',2); 
    set(gca,'Fontsize',16), xlabel('frequency (k)'), ylabel('fft(S(t)*g(t-\tau))')
    xlim([200,800])
    drawnow
end

% Same thing but centered on each bar to find all the frequencies in the
% bar
figure(6)


pcolor(tau,ks,Sgt_spec)
shading interp
set(gca,'ylim',[200 800],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')



%%

a = 40;
tau = 0:0.03:tr_guns;

for j = 1:length(tau)
   g = exp(-a*(t - tau(j)).^2); % Window function
   Sg = g.*y_guns;
   Sgt = fft(Sg);
   Sgt_spec(:,j) = fftshift(abs(Sgt)); % We don't want to scale it
end

figure(6)
pcolor(tau(1:70),ks,Sgt_spec(:,1:70))
shading interp
set(gca,'ylim',[200 800],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (Hz)')
text(1,276,"C#",'Color','white','FontSize',14)
text(1,553,"C#",'Color','white','FontSize',14)
text(0,415,"G#",'Color','white','FontSize',14)
text(0,369,"F#",'Color','white','FontSize',14)
text(0,741,"F#",'Color','white','FontSize',14)
text(1,701,"F",'Color','white','FontSize',14)
title("Spectrogram for GNR")



%%


a = 5;
tau = 0:0.1:tr_floyd;
a_fbig = 0.001;
tau_fbig = 100;
a_overtone = 0.05;

for j = 1:length(tau)
   g = exp(-a*(t_f - tau(j)).^2); % Window function
   % g_f = exp(-a_fbig*(ks_f-tau_fbig).^2); % For algorithm 3
   
   Sg = g.*y_floyd;
   Sgt = fft(Sg);
   Sgtshift = fftshift(abs(Sgt));
   [m, ind] = max(Sgtshift(330406:331464));
   freq = ks_f(ind + 330406)*2;
   
   filter1 = 1 - exp(-a_overtone*(ks_f-freq*1).^2);
   filter2 = 1 - exp(-a_overtone*(ks_f-freq*2).^2);
   filter3 = 1 - exp(-a_overtone*(ks_f-freq*3).^2);
   filter4 = 1 - exp(-a_overtone*(ks_f-freq*4).^2);
   filter5 = 1 - exp(-a_overtone*(ks_f-freq*5).^2);

   Sgt_spec(:,j) = Sgtshift.*filter1.*filter2.*filter3.*filter4.*filter5;
end

figure(6)
pcolor(tau,ks_f,Sgt_spec)
shading interp
set(gca,'ylim',[220 600],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (Hz)')
text(6,375,"F#",'Color','white','FontSize',14)
text(1.5,332,"E",'Color','white','FontSize',14)
text(11,587,"D",'Color','white','FontSize',14)
text(10,250,"B",'Color','white','FontSize',14)
title("Spectrogram for Pink Floyd Guitar")


%%


% a_f = 0.00000000000001;
% tau_f = 600;
% g_f = exp(-a_f*(ks_f-tau_f).^2);




freqs = fftshift(fft(y_floyd));
z = zeros([1 length(ks_f)]);
for i = 1:length(ks_f)
    if (-70 > ks_f(i) && ks_f(i) > -130) || (70 < ks_f(i) && ks_f(i) < 130)
        z(i) = freqs(i);
    end
end
z2 = ifft(ifftshift(z));
p8_filtered = audioplayer(abs(z2),Fs_floyd); playblocking(p8_filtered);




