% Clean workspace
clear all; close all; clc
% Import and convert audio pieces to vectors (Floyd)
[y2, Fs2] = audioread('Floyd.m4a');
tr_floyd = length(y2)/Fs2; 
y2s = y2(1:length(y2)-1);
% Set up time and Fourier domain (Floyd)
L2 = tr_floyd; n2 = length(y2s);
t2s = linspace(0,L2,n2+1); ts = t2s(1:n2);
k2 = (2*pi/L2)*[0:n2/2-1 -n2/2:-1]; ks2 = fftshift(k2);
% Apply Gabor filter on Floyd using Gaussian function
tau2 = 0:2:L2; a2 = 1000;
ygt2_spec = zeros(n2,length(tau2));
yg2 = zeros(n2,length(tau2));
for j = 1:length(tau2)
   g = exp(-a2*(ts-tau2(j)).^2); % Window function
   yg2(:,j) = g.*y2s';
   ygt2 = fft(yg2(:,j));
   ygt2_spec(:,j) = fftshift(abs(ygt2));
end
figure(1)
pcolor(tau2,ks2/(2*pi),ygt2_spec/(2*pi))
shading interp
set(gca,'ylim',[0 800],'Fontsize',16)
colormap(hot)
xlabel('time (s)'), ylabel('frequency (Hz)') title('Floyd Music Score')
%% Isolate guitar
ygtf2_spec = zeros(n2,length(tau2));
ksfilter = abs(ks2/(2*pi)-250);
[kmin,kind] = min(ksfilter);
arrayfilter = zeros(n2,1);
ksfilter2 = abs(ks2/(2*pi)-1000);
[kmin2,kind2] = min(ksfilter2);
arrayfilter(kind+1:kind2,1) = 1;
for j = 1:length(tau2)
    guitarrange = arrayfilter.*ygt2_spec(:,j);
    [m,ind] = max(guitarrange);
    k0 = abs(ks2(ind));
    taug2 = 0.001;
    filter = exp(-taug2*(ks2-k0).^2);
    ygtf2_spec(:,j) = filter'.*guitarrange;
end
% plot the guitar notes
figure(3)
pcolor(tau2,ks2/(2*pi),ygtf2_spec/(2*pi))
shading interp
set(gca,'ylim',[300 1000],'Fontsize',16)
colormap(hot)
xlabel('time (s)'), ylabel('frequency (Hz)') title('Notes of Guitar in Comfortably Numb')

yyaxis right
pcolor(tau2,ks2/(2*pi),ygtf2_spec/(2*pi))
shading interp
set(gca,'ylim',[300 1000],'Fontsize',16)
colormap(hot)
ylabel('Notes of guitar')
set(gca,'ytick',[330,370,440,494,587,740,880])
set(gca,'Yticklabel',{'E','F#','A','B','D','F#','A'})