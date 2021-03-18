% Clean workspace
clear all; close all; clc

% -- Load Movie
load cam1_3.mat; [h1,w1,rgb1,num1] = size(vidFrames1_3);
load cam2_3.mat; [h2,w2,rgb2,num2] = size(vidFrames2_3);
load cam3_3.mat; [h3,w3,rgb3,num3] = size(vidFrames3_3);
num = min([num1,num2,num3]); t = 1:num;

% trim movie clips
vidF1 = vidFrames1_3;
vidF2 = vidFrames2_3(:,:,:,1:num);
vidF3 = vidFrames3_3(:,:,:,1:num);

% case1
% range1 = [300 200 50 250]; range2 = [220 60 70 320]; range3 = [120 300 150 200];
%case2
% range1 = [300 200 100 250]; range2 = [160 50 250 350]; range3 = [120 300 180 220];
%case3 & case4
range1 = [300 200 150 220]; range2 = [160 50 300 350]; range3 = [120 150 250 420];

vids = {vidF1,vidF2,vidF3};
ranges = {range1,range2,range3};
X = zeros(6,num);
%% plot position
for i=1:3
    for j=1:num
        M = rgb2gray(vids{i}(:,:,:,j));
        if i==3
            M = M';
        end
        M = imcrop(M,ranges{i});
        light = max(M,[],'all');
        [lightx, lighty] = ind2sub(size(M), find(M == light));
        X(2*i-1,j) = mean(lightx);
        X(2*i,j) = mean(lighty);
    end 
    X(2*i-1,:) = (X(2*i-1,:)-mean(X(2*i-1,:)))/sqrt(num-1);
    X(2*i,:) = (X(2*i,:)-mean(X(2*i,:)))/sqrt(num-1);
    
    figure(1)
    subplot(1,3,1)
    title('The Position of Mass in Z-direction(Case3)', 'fontsize', 14)
    plot(t, X(2*i-1,:)), hold on
    xlabel('Frames Numbers'), ylabel('Relative Scaled Position')
    axis([0 num -15 15])
    
    subplot(1,3,2) 
    title('The Position of Mass in XY-plane(Case3)', 'fontsize', 14)
    plot(t, X(2*i,:)), hold on
    xlabel('Frames Numbers'), ylabel('Relative Scaled Position')
    axis([0 num -15 15])
end
%% Case PCA analysis
[U,S,V] = svd(X, 'econ');
subplot(1,3,3) % PCA
plot(t, V(:,1:3)'),legend('PC1', 'PC2', 'PC3')
title('Principal Component Oscillation(Case3)', 'fontsize', 14)
xlabel('Frames Numbers'), ylabel('Relative Scaled Position')
axis([0 num -0.3 0.3])