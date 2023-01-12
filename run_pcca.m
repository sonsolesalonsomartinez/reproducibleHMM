function [coeff, scores, e, Q] = run_pcca(bnameHmm,index,K,P,ndim1,X)
    
    R = length(index);
    
    % Concatenate Gammas across states
    Midx  = repelem(1:R,K);
    M = zeros(ndim1,R*K);
    for r = 1:R 
        M(:,Midx==r) = cell2mat(struct2cell(load([bnameHmm num2str(index(r)) '.mat'],'Gamma')));
    end 
    
    % Run PCA
    M=single(M);
    [coeff,scores,e] = pca(M,'NumComponents',P);

    % Get Covariance Matrix
    if nargout ==4
        Q = getCovariancePC(X, bnameHmm, coeff,index, P);
    end


end