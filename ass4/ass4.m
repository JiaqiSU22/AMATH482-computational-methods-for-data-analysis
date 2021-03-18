clear all; close all; clc

%% load data
[trainingdata, traingnd] = mnist_parse('train-images-idx3-ubyte', 'train-labels-idx1-ubyte');
trainingdata = double(reshape(trainingdata, size(trainingdata,1)*size(trainingdata,2), []).');
tr_data = trainingdata';
[m, n] = size(tr_data);
pxl = sqrt(m);
traingnd = double(traingnd);

[testdata, testgnd] = mnist_parse('t10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte');
testdata = double(reshape(testdata, size(testdata,1)*size(testdata,2), []).');
ts_data = testdata';
testgnd = double(testgnd);
%% SVD Analysis
[U, S, V] = svd(tr_data, 'econ');

figure(1) % singular value spectrum

subplot(2,1,1)
plot(diag(S),'o','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 100])
title('Singular Value Spectrum'), xlabel('modes'), ylabel('value')

subplot(2,1,2)
semilogy(diag(S),'o','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 100])
title('Singular Value Spectrum (log)'), xlabel('modes'), ylabel('value')

%% Projection onto 3 V-modes
figure(2)
for i=1:9
    label_ind = find(traingnd == i); % locate image labeled number i
    plot3(V(label_ind, 2), V(label_ind, 3), V(label_ind, 5),...
        'o', 'DisplayName', sprintf('%i',i), 'Linewidth', 2)
    hold on
end
xlabel('2nd V-Mode'), ylabel('3rd V-Mode'), zlabel('5th V-Mode')
title('3D Projection onto V-modes 2, 3, 5')
legend
set(gca,'Fontsize', 16)

%% LDA for 2 digits
feature = 30;
sucRate2 = zeros(10,10);
for i = 0:9
    for j = i+1:9
        dg1 = i; dg2 = j;
        dg1_ind = find(traingnd == dg1); dg1_data = tr_data(:,dg1_ind);
        dg2_ind = find(traingnd == dg2); dg2_data = tr_data(:,dg2_ind);

        [U2,S2,V2,threshold,w,sortdog,sortcat] = dc_trainer(dg1_data, dg2_data,feature);

        dg1_indts = find(testgnd == dg1);
        dg2_indts = find(testgnd == dg2);
        dg1_datats = ts_data(:,dg1_indts);
        dg2_datats = ts_data(:,dg2_indts);
        TestSet = [dg1_datats dg2_datats];
        TestNum = size(TestSet,2);
        TestLabel = zeros(1, TestNum);
        TestLabel(1, size(dg1_datats, 2)+1:end) = 1;
        TestMat = U2'*TestSet; % PCA projection
        pval = w'*TestMat;

        %  Check pval against threshold
        ResVec = (pval > threshold);
        err = abs(ResVec - TestLabel);
        errNum = sum(err);
        sucRate2(i+1,j+1) = 1 - errNum/TestNum;
    end
end
max(max(sucRate2))
min(min(sucRate2))
%% LDA for 3 digits
sucRate3 = zeros(1,120);
times = 0;
feature = 30;
for i = 0:9
    for j = i+1:9
        for k = j+1:9
            times = times+1;
            dg1 = i; dg2 = j; dg3 = k;
            dg1_ind = find(traingnd == dg1); dg1_data = tr_data(:,dg1_ind);
            dg2_ind = find(traingnd == dg2); dg2_data = tr_data(:,dg2_ind);
            dg3_ind = find(traingnd == dg3); dg3_data = tr_data(:,dg3_ind);

            [U3,S3,V3,thd1,thd2,w3,max_ind,min_ind] = dg3_trainer(dg1_data, dg2_data, dg3_data,feature);

            dg1_indts = find(testgnd == dg1);
            dg2_indts = find(testgnd == dg2);
            dg3_indts = find(testgnd == dg3);
            dg1_datats = ts_data(:,dg1_indts); n1_ts = size(dg1_datats, 2);
            dg2_datats = ts_data(:,dg2_indts); n2_ts = size(dg2_datats, 2);
            dg3_datats = ts_data(:,dg3_indts); n3_ts = size(dg3_datats, 2);
            TestSet = [dg1_datats dg2_datats dg3_datats];
            TestNum = size(TestSet,2);
            TestLabel = ones(1, TestNum);
            TestLabel(1, n1_ts+1:n2_ts) = 2;
            TestLabel(1, n1_ts+n2_ts+1:end) = 3;
            TestMat = U3'*TestSet; % PCA projection
            pval = w3'*TestMat;

            %  Check pval against threshold
            ResVec = zeros(1, feature);
            for t = 1:TestNum
                if pval(t) > thd1
                    ResVec(t) = max_ind;
                elseif pval(t) > thd2
                    ResVec(t) = 6-max_ind-min_ind;
                else
                    ResVec(t) = min_ind;
                end
            end
            err = (ResVec == TestLabel);
            errNum = sum(abs(err));
            sucRate3(times) = 1 - errNum/TestNum;
        end
    end
end
max(sucRate3)
min(sucRate3)

%% decision tree
% classification tree
tree = fitctree(tr_data',traingnd,'CrossVal','on');
view(tree.Trained{1},'Mode','graph');
classError = kfoldLoss(tree);


%% svm
% SVM classifier with training data, labels and test set
train = S*V'./(max(max(S))*V');
Mdl = fitcecoc(train',traingnd);
[Uts, Sts, Vts] = svd(ts_data, 'econ');
test = Sts*Vts'./(max(max(S))*V');
test_labels = predict(Mdl,test);




