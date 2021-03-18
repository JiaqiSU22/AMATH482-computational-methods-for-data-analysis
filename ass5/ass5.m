clear all; close all; clc

%% load data
v = VideoReader('monte_carlo_low.mp4');
frames = read(v);
[m,n,rgb,num_frames] = size(frames);
t = 1:10:num_frames;
dt = 10;
X = zeros(m*n,length(t));
for i = 1:length(t)
    single_frame = frames(:,:,:,(i-1)*dt+1);
    F = rgb2gray(single_frame);
    F = im2double(F);
    X(:,i) = reshape(F,[m*n,1]);
    imshow(F)
end
%% Create DMD matrices
X1 = X(:,1:end-1);
X2 = X(:,2:end);
%% SVD of X1 and Computation of ~S
[U, Sigma, V] = svd(X1,'econ');
U = U(:,1:2);
Sigma = Sigma(1:2,1:2);
V = V(:,1:2);
S = U'*X2*V*diag(1./diag(Sigma));
[eV, D] = eig(S); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues
omega = log(mu)/dt;
Phi = U*eV;
%% Create DMD Solution
y0 = Phi\X1(:,1); % pseudoinverse to get initial conditions
ubg_modes = zeros(length(y0),length(t));
for i = 1:length(t)
   ubg_modes(:,i) = y0.*exp(omega*t(i)); 
end
ubg_dmd = Phi*ubg_modes;
ufg_dmd = abs(X-ubg_dmd);
res = ufg_dmd.* (ufg_dmd < 0);
ubg_dmd = abs(ubg_dmd)+res;
ufg_dmd = ufg_dmd-res;
%%
for i = 1:length(t)
    frame_img = reshape(ubg_dmd(:,i),[m,n]);
    frame_img = mat2gray(frame_img);
    imshow(frame_img)
end
%%

for i = 1:length(t)
    frame_img = reshape(ufg_dmd(:,i),[m,n]);
    frame_img = mat2gray(frame_img);
    imshow(frame_img)
end



% %% Plotting Eigenvalues (omega)
% 
% % make axis lines
% line = -15:15;
% 
% plot(zeros(length(line),1),line,'k','Linewidth',2) % imaginary axis
% hold on
% plot(line,zeros(length(line),1),'k','Linewidth',2) % real axis
% plot(real(omega)*dt,imag(omega)*dt,'r.','Markersize',15)
% xlabel('Re(\omega)')
% ylabel('Im(\omega)')
% set(gca,'FontSize',16,'Xlim',[-1.5 0.5],'Ylim',[-1 1])