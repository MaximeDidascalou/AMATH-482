function sucRate = LDA(U,digits,feature,labels_t,labels,dig1,dig2,A,num)
    digit1 = digits(1:feature,find(labels_t == dig1));
    digit2 = digits(1:feature,find(labels_t == dig2));
    %scatter matrices:
    m1 = mean(digit1,2);
    m2 = mean(digit2,2);
    Sw = 0; % within class variances
    for k = 1:num(dig1+1)
        Sw = Sw + (digit1(:,k) - m1)*(digit1(:,k) - m1)';
    end
    for k = 1:num(dig2+1)
        Sw =  Sw + (digit2(:,k) - m2)*(digit2(:,k) - m2)';
    end
    Sb = (m1-m2)*(m1-m2)'; % between class
    % Find the best projection line
    [V2, D] = eig(Sb,Sw); % linear disciminant analysis
    [lambda, ind] = max(abs(diag(D)));
    w = V2(:,ind);
    w = w/norm(w,2);
    % Project onto w
    v1 = w'*digit1;
    v2 = w'*digit2;
    % Make 0 below the threshold
    if mean(v1) > mean(v2)
        w = -w;
        v1 = -v1;
        v2 = -v2;
    end
    % Find the threshold value
    sort1 = sort(v1);
    sort2 = sort(v2);
    t1 = length(sort1);
    t2 = 1;
    while sort1(t1) > sort2(t2)
        t1 = t1 - 1;
        t2 = t2 + 1;
    end
    threshold = (sort1(t1) + sort2(t2))/2;
    % classifying unknown digits
    labelsnew = labels(labels == dig1 | labels == dig2);
    labelsnew = labelsnew == dig2;
    Y = U'*A(:,labels == dig1 | labels == dig2); % PCA projection
    pval = w'*Y(1:feature,:);
    ResVec = (pval > threshold);
    err = abs(ResVec - labelsnew');
    errNum = sum(err);
    sucRate = 1 - errNum/size(labelsnew,1);
end
