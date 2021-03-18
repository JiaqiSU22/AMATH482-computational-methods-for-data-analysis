% Clean workspace
clear all; close all; clc

% Import the data as the 262144x49 (space by time) matrix called subdata
load subdata.mat 

% Set up spatial and Fourier domain
L=10; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1);
x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1];
ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% Average through spectrum
Unsum=zeros(n,n,n);
for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    Unsum=Unsum+fftn(Un);
end
Unave=abs(fftshift(Unsum))/49;
[M,I]=max(abs(Unave),[],'all','linear');
ix=Kx(I); iy=Ky(I); iz=Kz(I);

% Plot the surface
figure(1)
isosurface(Kx,Ky,Kz,abs(Unave)/M,0.7), grid on
title('3D Frequency Isourface Graph')
xlabel('X frequency(kx)'); ylabel('Y frequency(ky)'); zlabel('Z frequency(kz)');

% Create filter
tau=0.5;
xk0=ix; yk0=iy; zk0=iz;
filter=exp(-tau*(Kx-xk0).^2).*exp(-tau*(Ky-yk0).^2).*exp(-tau*(Kz-zk0).^2); 

% Find coordinate at each time
cor=zeros(3,49);
for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    Unt=filter.*fftshift(fftn(Un));
    Untf=ifftn(Unt);
    [Mt,idx]=max(abs(Untf),[],'all','linear');
    cor(1,j)=X(idx);
    cor(2,j)=Y(idx);
    cor(3,j)=Z(idx);
end

% Plot the path
figure(2)
plot3(cor(1,:),cor(2,:),cor(3,:)), grid on, hold on
plot3(cor(1,end),cor(2,end),cor(3,end),'r*')
text(cor(1,end),cor(2,end),cor(3,end),...
[' (' num2str(cor(1,end)) ', ' num2str(cor(2,end)) ', ' num2str(cor(3,end)) ')'])
title('24-hour 3D Path Graph of the Submarine')
xlabel('X'); ylabel('Y'); zlabel('Z');



 