function sucRate = LDA(U,digits,feature,labels_t,labels,dig1,dig2,dig3,A,num)
    % Project onto PCA modes
    digit1 = digits(1:feature,find(labels_t == dig1));
    digit2 = digits(1:feature,find(labels_t == dig2));
    digit3 = digits(1:feature,find(labels_t == dig3));
    %scatter matrices:
    m1 = mean(digit1,2);
    m2 = mean(digit2,2);
    m3 = mean(digit3,2);
    mo = (m1+m2+m3)/3;
    Sw = 0; % within class variances
    for k = 1:num(dig1+1)
        Sw = Sw + (digit1(:,k) - m1)*(digit1(:,k) - m1)';
    end
    for k = 1:num(dig2+1)
        Sw =  Sw + (digit2(:,k) - m2)*(digit2(:,k) - m2)';
    end
    for k = 1:num(dig3+1)
        Sw =  Sw + (digit3(:,k) - m3)*(digit3(:,k) - m3)';
    end
    Sb1 = (m1-mo)*(m1-mo)'; Sb2 = (m2-mo)*(m2-mo)'; Sb3 = (m3-mo)*(m3-mo)';
    Sb = Sb1 + Sb2 + Sb3;
    % Find the best projection line
    [V2, D] = eig(Sb,Sw); % linear disciminant analysis
    [lambda, ind] = max(abs(diag(D)));
    w = V2(:,ind);
    w = w/norm(w,2);
    % Project onto w
    v1 = w'*digit1; v2 = w'*digit2; v3 = w'*digit3;
    
    % Find the threshold values
    [M,Imin] = min([sum(v1(1,:)) sum(v2(1,:)) sum(v3(1,:))]);
    [M,Imax] = max([sum(v1(1,:)) sum(v2(1,:)) sum(v3(1,:))]);
    Imid = 6 - Imin - Imax;
    sortt{1} = sort(v1(1,:)); sortt{2} = sort(v2(1,:)); sortt{3} = sort(v3(1,:));
    t1 = length(sortt{Imin}); t2 = 1;
    while sortt{Imin}(t1) > sortt{Imid}(t2)
        t1 = t1 - 1;
        t2 = t2 + 1;
    end
    threshold1 = (sortt{Imin}(t1) + sortt{Imin}(t2))/2;
    t2 = length(sortt{Imid}); t3 = 1;
    while sortt{Imid}(t2) > sortt{Imax}(t3)
        t2 = t2 - 1;
        t3 = t3 + 1;
    end
    threshold2 = (sortt{Imid}(t2) + sortt{Imax}(t3))/2;
    % classifying unknown digits
    digs = [dig1 dig2 dig3];
    Y = U'*A(:,labels == dig1 | labels == dig2| labels == dig3); % PCA projection
    pval = w'*Y(1:feature,:);
    labelsnew = labels(labels == dig1 | labels == dig2 | labels == dig3);
    err = zeros([1 length(pval)]);
    for i = 1:length(pval)
        if pval(i) < threshold1
            if digs(Imin) ~= labelsnew(i)
                err(i) = 1;
            end
        elseif pval(i) > threshold2
            if digs(Imax) ~= labelsnew(i)
                err(i) = 1;
            end
        else
            if digs(Imid) ~= labelsnew(i)
                err(i) = 1;
            end
        end
    end
    
    errNum = sum(err);
    sucRate = 1 - errNum/size(labelsnew,1);
end