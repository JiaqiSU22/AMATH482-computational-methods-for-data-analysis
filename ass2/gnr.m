% Clean workspace
clear all; close all; clc

% Import and convert audio pieces to vectors (GNR)
[y, Fs] = audioread('GNR.m4a');
tr_gnr = length(y)/Fs; % record time in seconds


% Set up time and Fourier domain (GNR)
L = tr_gnr;
n = length(y);
t2 = linspace(0,L,n+1);
t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Apply Gabor filter on GNR using Gaussian function
tau = 0:0.3:L;
a = 800;
ygt_spec = zeros(n,length(tau));
for j = 1:length(tau)
   g = exp(-a*(t - tau(j)).^2); % Window function
   yg = g.*y';
   ygt = fft(yg);
   ygt_spec(:, j) = fftshift(abs(ygt))/(2*pi);
end

figure(1)
yyaxis left
pcolor(tau, ks/(2*pi),ygt_spec)
shading interp
set(gca,'ylim',[0 1000],'Fontsize',16)
colormap(hot)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('GNR Guitar Music Score')

yyaxis right
pcolor(tau,ks/(2*pi),ygt_spec)
shading interp
set(gca,'ylim',[0 1000],'Fontsize',16)
colormap(hot)
ylabel('Notes of guitar')
set(gca,'ytick',[277,370,415,554,698,740])
set(gca,'Yticklabel',{'C#','F#','G#','C#','F','F#'})