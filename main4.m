clear all; close all; clc

[images_t, labels_t] = mnist_parse('train-images.idx3-ubyte', 'train-labels.idx1-ubyte');
[images, labels] = mnist_parse('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte');

[m,n,N] = size(images_t);
A = zeros(m*n,N);
for i = 1:N
    A(:,i) = reshape(images_t(:,:,i),m*n,1);
end
num = sum(labels_t == 0:9);
M = mean(A,2);
A_mean0 = A - M*ones([1 N]);

[m2,n2,N2] = size(images);
A2 = zeros(m2*n2,N2);
for i = 1:N2
    A2(:,i) = reshape(images(:,:,i),m2*n2,1);
end
M2 = mean(A2,2);
A2_mean0 = A2 - M2*ones([1 N2]);

%% svd decomposition
[U,S,V] = svd(A_mean0,'econ');
digits = S*V';

%% approximations of images: 
[U2,S2,V2] = svd(A,'econ');
r = [1 5 10 25 50 100 300 500];
figure(1)
for k = 1:8
   subplot(3,3,k)
   new = U2(:,1:r(k))*S2(1:r(k),1:r(k))*V2(:,1:r(k))';
   im = reshape(new(:,1),28,28);
   im = rescale(im);
   imshow(im)
   title(r(k))
end
subplot(3,3,9)
imshow(images_t(:,:,1))
title('Original')
sgtitle('Image reconstruction with different ranks')

%% plotting first 4 principal components
figure(2)
for k = 1:4
   subplot(2,2,k)
   ut1 = reshape(U(:,k),28,28);
   ut2 = rescale(ut1);
   imshow(ut2)
end

%% plotting singular values
figure(3)
subplot(2,1,1)
plot(diag(S),'k','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 m*n])
ylabel("\sigma")
subplot(2,1,2)
semilogy(diag(S),'k','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 m*n])
ylabel("\sigma (log scale)")

%% plotting projected data
figure(4)
col = [0 0 0;0 0 1;0 1 0;1 0 0;0 1 1;1 0 1;1 1 0;
    0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
for k = 0:9
    plot3(digits(1,find(labels_t == k)),digits(2,find(labels_t == k)),digits(3,find(labels_t == k)),'k.','color',col(k+1,:))
    hold on
end
[h,icons] = legend('0','1','2','3','4','5','6','7','8','9');
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',20);
xlabel("PC1")
ylabel("PC2")
zlabel("PC3")
title("Data projected on first 3 modes")

%% max and min succes rate
feature = 40;
min = 1;
max = 0;
min3 = 1;
max3 = 0;
sumrate = 0;
sumrate3 = 0;
for i = 0:9
   for j = i+1:9
        sucRate = LDA(U,digits,feature,labels_t,labels,i,j,A2_mean0,num);
        sumrate = sumrate + sucRate;
        if sucRate < min
            min = sucRate;
            mindigs = [i j];
        end
        if sucRate > max
            max = sucRate;
            maxdigs = [i j];
        end
        for k = j+1:9
            sucRate3 = LDA3(U,digits,feature,labels_t,labels,i,j,k,A2_mean0,num);
            sumrate3 = sumrate3 + sucRate3;
            if sucRate3 < min3
                min3 = sucRate3;
                mindigs3 = [i j k];
            end
            if sucRate3 > max3
                max3 = sucRate3;
                maxdigs3 = [i j k];
            end
        end
   end
end
sumrate = sumrate/45;
sumrate3 = sumrate3/120;

%% max and min succes rate
feature = 40;
min = 1;
max = 0;
min3 = 1;
max3 = 0;
sumrate = 0;
sumrate3 = 0;
for i = 0:9
   for j = i+1:9
        sucRate = LDA(U,digits,feature,labels_t,labels_t,i,j,A_mean0,num);
        sumrate = sumrate + sucRate;
        if sucRate < min
            min = sucRate;
            mindigs = [i j];
        end
        if sucRate > max
            max = sucRate;
            maxdigs = [i j];
        end
        for k = j+1:9
            sucRate3 = LDA3(U,digits,feature,labels_t,labels_t,i,j,k,A_mean0,num);
            sumrate3 = sumrate3 + sucRate3;
            if sucRate3 < min3
                min3 = sucRate3;
                mindigs3 = [i j k];
            end
            if sucRate3 > max3
                max3 = sucRate3;
                maxdigs3 = [i j k];
            end
        end
   end
end
sumrate = sumrate/45;
sumrate3 = sumrate3/120;
%% set up data
feature = 40;
X = digits(1:feature,:)';
Y = U'*A2_mean0;
Y = Y(1:feature,:)';

%% SVM for all 10:
Mdl = fitcecoc(X./max(X(:)),labels_t);
test_labels = predict(Mdl,Y./max(Y(:)));

% test accracy
err = sum(test_labels ~= labels);
succesrate = 1-err/N2

%% tree method for all 10:
tree=fitctree(X,labels_t,'CrossVal','on');
test_labels = predict(tree.Trained{1},Y);

% test it
err = sum(test_labels ~= labels);
succesrate = 1 - err/N2

%% set up data for 2 digits
dig1 = 4;
dig2 = 9;
Y2 = Y(labels == dig1 | labels == dig2,:);
X2 = X(labels_t == dig1 | labels_t == dig2,:);
newlab_t = labels_t(labels_t == dig1 | labels_t == dig2);
newlab = labels(labels ==dig1 | labels == dig2);

%% SVM for 2 digits
Mdl = fitcsvm(X2./max(X2(:)),newlab_t);
test_labels = predict(Mdl,Y2./max(Y2(:)));

% test accuracy
err = sum(test_labels ~= newlab);
succesrate = 1-err/length(test_labels);

%% tree for 2 digits:
tree=fitctree(X2,newlab_t,'CrossVal','on');
test_labels = predict(tree.Trained{1},Y2);

% test it
err = sum(test_labels ~= newlab);
succesrate = 1 - err/length(test_labels);
