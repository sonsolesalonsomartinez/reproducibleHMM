function Q = getCovariancePC(X, bnameHmm, coeff, index, P)
   
    % Static FC (mC)
    sFC=[];
    for s  = 1:length(X)
        sFC(s,:,:) = corrcoef(X{s,:});    
    end
    mC=squeeze(mean(sFC));
    
    % Coeff PCA (A)
    A = coeff(:,1:P);
    R = length(index);
    
    % States Covariance Matrix (C)
    m=0;% counter for each R * K = M states
    C=[];
    for r = 1:R
        hmm = cell2mat(struct2cell(load(sprintf('%s%d.mat',bnameHmm,index(r)), 'hmm')));
        for k = 1:hmm.train.K
            m = m+1;
            C(:,:,m) = getFuncConn(hmm,k,1);
        end 
    end
    
    % Compute Covariance
    Q = PCCA(A,C,mC);

end

function PC = PCCA(A,C,mC)
    % C are the covariance matrices: p x p x K
    % A are the PCA weights (first output argument of PCA), states x components  
    % mC is the static FC, which woudl be computed as the state average if not provided
    % PC are p x p x components 
    
    if nargin < 3
        mC = riemann_mean(C);
    end
    X = Tangent_space(C,mC); % features x states 
    Y = X * A; % features x components
    PC = UnTangent_space(Y,mC);

end