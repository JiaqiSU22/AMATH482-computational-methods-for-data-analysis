function dgData = digits_wavelet(dgfile)

    [m,n] = size(dgfile); 
    pxl = sqrt(m);
    nw = m/4; % wavelet resolution
    dgData = zeros(nw,n);
    
    for k = 1:n
        X = im2double(reshape(dgfile(:,k),pxl,pxl));
        [~,cH,cV,~]=dwt2(X,'haar');
        cod_cH1 = rescale(abs(cH));
        cod_cV1 = rescale(abs(cV));
        cod_edge = cod_cH1+cod_cV1;
        dgData(:,k) = reshape(cod_edge,nw,1);
    end
end
