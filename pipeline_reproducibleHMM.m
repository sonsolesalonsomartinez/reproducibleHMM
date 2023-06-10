% Code used in Alonso and Vidaurre (2023) "Stability of dynamic FC
% estimates in neuroimaging and electrophysiology: solutions and limits"
% 
% This script requires the HMM-MAR toolbox (available at
% https://github.com/OHBA-analysis/HMM-MAR)
%
% This script assumes that all the preprocessing has been done.
%
% To access the (preproceseds and parcellated) time series for all subjects,
% please follow this link:
% https://ora.ox.ac.uk/objects/uuid:2770bfd4-6ab8-4f1e-b5e7-06185e8e2ae1.
% Note that the analysis was conducted on the top 10 subjects 
%
% INPUT
% X    observations; 1xN cell containing the time series per subject,
% T    length of time series; 1xN cell containg the lengths of the timeseries
%
% This pipeline must be adapted to your particular configuration of files.
%
% Sonsoles Alonso
% Aarhus University, June 2023
%%%%%%%%%%%%%%%%%%%%%%%%%

% USER
R = 1000;     % Number of HMM runs
K = 6;        % Number of HMM states per run
maxclust = K; % Number of clusters in HC-HMM       

% Setup directories and filenames
mypath  = pwd;
codedir = '/path/2/toolboxes';
addpath(genpath([codedir '/ohba_analysis/HMM-MAR']));
hmmfile  = [mypath '/hmm_run'];     % Filename for HMM runs 
Rhmmfile= [mypath '/robusthmm'];    % Filename for BR-HMM and HC-HMM results 

% Read data
datafile = [mypath '/dataMEG_sample.mat']; % File for data timeseries
load(datafile,'X','T','Hz');
nRois = size(X{1},2);

% HMM model parameters
opt = struct();
opt.K = K; % number of states 
opt.order = 0; % no autoregressive components
opt.zeromean = 1; % do not model the mean
opt.covtype = 'full'; % full covariance matrix
opt.useParallel = 1;
opt.cyc = 100;
opt.standardise = 1;
opt.verbose = 1;
opt.inittype = 'HMM-MAR';
opt.initcyc = 10;
opt.initrep = 1;
%HMM model for MEG 
%comment the following 4 lines if analysing fMRI data
opt.onpower = 1;     % on power rather than raw time series
opt.filter = [8 12]; % alpha band
opt.Fs = Hz;
opt.downsample=100;  % downsample to speed computation


%% Run HMM R times
for r = 1:R
    [hmm, Gamma,~,~,~,~,fehist] = hmmmar(X,T,opt);
    f = fehist(end);
    save([hmmfile num2str(r) '.mat'], 'hmm', 'Gamma','f'); 
end

%% BR-HMM
F=nan(R,1);
for r=1:R
    load([hmmfile num2str(r) '.mat'],'f');
    F(r) = f;
end

% Get BR run
[~, br_run] = min(F);

% Save
save([Rhmmfile '_rep1.mat'],'br_run');


%% HC-HMM

% Concatenate gammas across HMM runs
gammaM=cell(R,1);
for r = 1:R
    load([hmmfile num2str(r) '.mat'],'Gamma');
    gammaM{r} = single(Gamma);
end
M = cat(2, gammaM{:});
Midx = [repelem(1:R,K); repmat(1:K,1,R)]; % [ runs ; states ]

% compute between-state correlations
C = corrcoef(M);

% Convert similarity matrix to dissimilarity matrix
D = 1-C;

% Perform hierarchical clustering
Z = linkage(D,'Ward');

% Cut the tree to form clusters
clusters = cluster(Z,'maxclust',maxclust);

% Calculate aggregated state timeseries and covariances
hc_gamma = zeros(size(M,1),maxclust);
hc_fc = nan(nRois,nRois,maxclust);

for i = 1:maxclust

    % States within cluster i
    idx = clusters==i;

    % States timeseries within cluster i
    cluster_states = M(:, idx); % get the state time series for cluster i

    % Aggregate state timeseries (i.e., hc_gamma) within cluster i
    hc_gamma(:,i) = mean(cluster_states,2);   % compute the mean across states
    
    % Aggregate covariances (i.e., hc_fc) within cluster i
    rr = Midx(1,idx); % index run
    kk = Midx(2,idx); % index state
    fc = nan(length(rr),nRois,nRois);

    for ii  = 1:length(rr)
        
        % Load hmm of run ii
        hmm_ii = cell2mat(struct2cell(load([hmmfile num2str(rr(ii)) '.mat'],'hmm')));
        
        % compute covariance matrix of orignal state
        covMatrix = getFuncConn(hmm_ii,kk(ii));
        
        % Weight of each state (i.e., Fractional Occupancy, FO)
        FO = mean(cluster_states(:,ii));
        
        % Store each covMatrix weighted by its FO
        fc(ii,:,:) = covMatrix * FO;
    end

    % Compute the weighted average of the covariances across states within cluster i
    hc_fc(:,:,i) = mean(fc);

end

% Normalize the rows of a matrix so that each row represents a probability distribution
hc_gamma = hc_gamma ./ sum(hc_gamma,2);

%Save
save([Rhmmfile '_rep1.mat'], 'clusters','Z','C','hc_gamma','hc_fc','-append')

%% BR-HMM stability

% By running the HMM R times again and selecting the one with minimum free
% energy; the stability of BR-HMM can be computed by comparing the
% similarity acros BR-HMM repetitions:

% load([Rhmmfile '_rep1.mat'], 'br_run');
% load([hmmfile num2str(br_run) '.mat'],'Gamma');
% br_gamma1 = Gamma;

% load([Rhmmfile '_rep2.mat'], 'br_run');
% load([hmmfile num2str(br_run) '.mat'],'Gamma');
% br_gamma2 = Gamma;

% S = [S,getGammaSimilarity(br_gamma1,br_gamma2)];


%% HC-HMM stability

% By repeating the HC-HMM with a different set of HMM runs, the similarity
% across hc_gamma repetitions can be computed as:

% load([Rhmmfile '_rep1.mat'], 'hc_gamma');
% hc_gamma1 = hc_gamma;

% load([Rhmmfile '_rep2.mat'], 'hc_gamma');
% hc_gamma2 = hc_gamma;

% S = [S,getGammaSimilarity(hc_gamma1,hc_gamma2)];

%% FIGURE

nrows = 3;ncols = 2; figure;

% Between-run similarity
SS = ones(R,R);
for r1=1:R-1    
    gamma1 = gammaM{r1};  
    for r2=r1+1:R       
        gamma2 = gammaM{r2};        
        S = getGammaSimilarity(gamma1,gamma2); 
        SS(r1, r2) = S;
        SS(r2, r1) = S;
    end   
end

% Pannel A: Histogram similarity
subplot(nrows,ncols,1:2);
isubd = find(triu(ones(size(SS)),1));
histogram(SS(isubd));
hold on;xline(min(SS(isubd)),'k-','minimum');
hold on;xline(mean(SS(isubd)),'k-','mean');
hold on;xline(max(SS(isubd)),'k-','maximum');
xlabel ('between-run similarity')
ylabel({['no. pairs from'] ;[ num2str(R) ' HMM runs']})

% Pannel B: Between-run similarity sorted by FE
subplot(nrows,ncols,ncols+1);
Fnorm = (F - min(F)) / (max(F) - min(F));
[~,orderfe]=sort(Fnorm,'ascend');
imagesc(SS(orderfe,orderfe))
xlabel({[ num2str(R) ' runs'] ;'(ranked order)'});
set(gca,'XTick',[],'YTick',[])
cb=colorbar();cb.Label.String ={'between-run similarity'};

% Pannel B: FE bar plot
subplot(nrows,ncols,ncols+2);
barh(Fnorm(orderfe(end:-1:1)))
xlabel('free energy')
ylabel({[ num2str(R) ' runs'] ;'(ranked order)'})

% Pannel C: HC-HMM Dendrogram
load([Rhmmfile '_rep1.mat'], 'clusters','C','Z')

subplot(nrows,ncols,2*ncols+2)
cutoff = Z(end-maxclust+2,3)-eps; % Find cutoff for K clusters
[~,~, ord] = dendrogram(Z, 0,'ColorThreshold',cutoff,'Orientation','right');
set(gca,'YTick',[],'XTick',[])
xlabel('distance')
ylabel([ num2str(R) ' x ' num2str(K) ' states'])

% Pannel C: Between-state correlation matrix
subplot(nrows,ncols,2*ncols+1);
imagesc(C(ord(end:-1:1),ord(end:-1:1)))
set(gca,'XTick',mytick,'YTick',mytick,...
    'XTickLabel',[],'YTickLabel',[]),grid on
xlabel([ num2str(R) ' x ' num2str(K) ' states'])
cb=colorbar();cb.Label.String = {'between-state correlation'};
