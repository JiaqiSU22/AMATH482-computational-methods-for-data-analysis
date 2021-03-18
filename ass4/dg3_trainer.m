function [U,S,V,threshold1,threshold2,w,max_ind,min_ind] = dg3_trainer(wave1,wave2,wave3,feature)
    
    n1 = size(wave1,2); n2 = size(wave2,2); n3 = size(wave3,2);
    [U,S,V] = svd([wave1 wave2 wave3],'econ'); 
    total_wave = S*V';
    U = U(:,1:feature); % Add this in
    group1 = total_wave(1:feature,1:n1);
    group2 = total_wave(1:feature,n1+1:n1+n2);
    group3 = total_wave(1:feature,n1+n2+1:n1+n2+n3);
    m1 = mean(group1,2); m2 = mean(group2,2); m3 = mean(group3,2);
    tol_mean = mean(total_wave(1:feature,:),2);

    Sw = 0;
    for k=1:n1
        Sw = Sw + (group1(:,k)-m1)*(group1(:,k)-m1)';
    end
    for k=1:n2
        Sw = Sw + (group2(:,k)-m2)*(group2(:,k)-m2)';
    end
    for k=1:n3
        Sw = Sw + (group3(:,k)-m3)*(group3(:,k)-m3)';
    end
    Sb = (m1-tol_mean)*(m1-tol_mean)'+(m2-tol_mean)*(m2-tol_mean)'+(m3-tol_mean)*(m3-tol_mean)';

    [V3_eig,D3] = eig(Sb,Sw);
    [lambda,ind] = max(abs(diag(D3)));
    w = V3_eig(:,ind); w = w/norm(w,2);
    v1 = w'*group1; v2 = w'*group2; v3 = w'*group3;
    mean_ind = [mean(v1) mean(v2) mean(v3)];
    [max_mean, max_ind] = max(mean_ind);
    [min_mean, min_ind] = min(mean_ind);
    proj = {v1, v2, v3};

    sort1 = sort(proj{max_ind});
    sort2 = sort(proj{6-max_ind-min_ind});
    sort3 = sort(proj{min_ind});
    t1 = length(sort1); t2 = 1;
    while (sort1(t1)>sort2(t2) && (t1>1 || t2<length(sort2)))
        t1 = max(t1-1,1);
        t2 = min(t2+1,length(sort2));
    end
    threshold1 = (sort1(t1)+sort2(t2))/2;

    t1 = length(sort2);
    t2 = 1;
    while (sort2(t1)>sort3(t2) && (t1>1 || t2<length(sort3)))
        t1 = max(t1-1,1);
        t2 = min(t2+1,length(sort3));
    end
    threshold2 = (sort2(t1)+sort3(t2))/2;
end

