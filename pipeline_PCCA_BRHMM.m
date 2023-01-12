% Code used in Alonso Martinez et al. (2023) "Stability of dynamic FC
% estimates in neuroimaging and electrophysiology: solutions and limits"
% 
% This script requires the HMM-MAR toolbox (available at
% https://github.com/OHBA-analysis/HMM-MAR), the osl toolbox (available at
% https://github.com/OHBA-analysis/osl-core) and the covariance toolbox
% (available at https://github.com/alexandrebarachant/covariancetoolbox)
%
% This script assumes that all the preprocessing has been done 
%
% INPUT
% data          observations; cell containing the time series,
% T             length of time series

% This pipeline must be adapted to your particular configuration of files.
%
% Sonsoles Alonso Martinez
% Aarhus University, October 2023
%%%%%%%%%%%%%%%%%%%%%%%%%

mypath = pwd;

addpath(genpath('../ohba_analysis/HMM-MAR'));
addpath(genpath('../ohba_analysis/osl-core'));
addpath(genpath('../covariancetoolbox'));

% Read data
load([mypath 'dataMEG.mat'],'X','T');
Hz = 250;
ndim1=sum(cell2mat(T)); % total number of time points

% Setup filenames and data variables
bnameHmm   = [mypath '/hmm_run'];         % Basename for HMM results 
bnamePcca  = [mypath '/pcca_rep'];        % Basename for PCCA results 
fnamePcca  = [bnamePcca '1_200runs.mat']; % Filename of 1st PCCA rep on 200 runs
fnameFig   = [mypath '/figsData.mat'];    % Filename of data to plot


%% 1. - Run the HMM 1000 times
% HMM model parameters
K=6;

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
%comment the following lines if analysing fMRI data
opt.onpower = 1;     % on power rather than raw time series
opt.filter = [8 12]; % alpha band
opt.Fs = Hz;


nTotalRuns = 1000; 
for n = 1:nTotalRuns
    [hmm, Gamma,~,~,~,~,fehist] = hmmmar(X,T,opt);
    f = fehist(end);
    save([bnameHmm num2str(n) '.mat'], 'hmm', 'Gamma','f'); 
end

% Compute Between-run similarity for all runs
SSall=ones(nTotalRuns,nTotalRuns);
for ri=1:nTotalRuns-1
    gamma1 = cell2mat(struct2cell(load([bnameHmm num2str(ri) '.mat'],'Gamma')));
    S=[];
    for rj=ri+1:nTotalRuns
        gamma2 = cell2mat(struct2cell(load([bnameHmm num2str(rj) '.mat'],'Gamma')));
        S = [S,getGammaSimilarity(gamma1,gamma2)];
    end
    SSall(ri,ri+1:end) = S;
    SSall(ri+1:end,ri) = S;
end
save(fnameFig,'SSall');


%% 2. - BR-HMM

% Store Free Energy for all runs
Fall=[];
for r=1:nTotalRuns
    Fall(r)=cell2mat(struct2cell(load([bnameHmm num2str(r) '.mat'],'f')));
end

% Identify BR runs across varying number of runs (R) for Nr repetitions of BR-HMM procedure
% USER
Nr=10; % Number of BR-HMM repetitions
nRuns2test2 = 5:5:300; % BR-HMM to be run on R runs

idxRR=nan(max(nRuns2test2),Nr);
for R = nRuns2test2
    for n=1:Nr
        idx = randperm(nTotalRuns,R);
        [~, i] = min(Fall(idx));
        idxRR(R,n) = idx(i);
    end
end

% Compute Between-state-correlation for BR runs
idxRR_unique = unique(idxRR(nRuns2test2,:));
M = nan(ndim1,length(idxRR_unique)*K);
CCidx = repelem(idxRR_unique,K);
for i = 1:length(idxRR_unique)
    r=idxRR_unique(i);
    gamma = cell2mat(struct2cell(load([bnameHmm num2str(r) '.mat'],'Gamma')));
    M(:,CCidx == r) = gamma;
end 
CC = corrcoef(M);
save(fnameFig, 'Fall','idxRR','CC','CCidx','-append');

%% 3. - Run PCCA 30 times for varying number of runs

% USER
nRuns2test = [20:10:100 125 150 200];% PCA to be run on R runs
Np = 30; % number of PCA repetitions
P = 20;  % Number of principal components

for R = nRuns2test
    for n = 1:Np
        sprintf('Running PCA (repetition %d) on %d hmm runs',n,R)        
        % select R runs from the nTotalRuns (=1000) collected in step-1
        rng shuffle
        index = run_index(randperm(nTotalRuns,R));
        % run PCA on the R x K states and compute PC covariance matrix (Q)
        [coeff, scores, e, Q] = run_pcca(bnameHmm,index,K,P,ndim1,X);
        % save results
        save([bnamePcca num2str(n) '_' num2str(R) 'runs.mat'],...
            'coeff','scores','e','Q','index');
    end 
end

% PCCA stability
npairs = (Np*(Np-1))/2;
CS = [];
CSidx = repelem(nRuns2test,npairs)';
for R = nRuns2test
    sprintf('Running correlations between pairs of PCs for PCA on %d runs',n)
    scoresn=zeros(Np,ndim1,P);  
    for n = 1:Np
        load([bnamePcca num2str(n) '_' num2str(R) 'runs.mat'],'scores');
        scoresn(n,:,:) = scores(:,1:P);
    end
    % compute corelations
    for n1 = 1:Np-1
        scores_n1 = squeeze(scoresn(n1,:,:));
        for n2 = n1+1:Np    
            scores_n2 = squeeze(scoresn(n2,:,:));
            S = corr(scores_n1, scores_n2);
            CS = [CS; diag(S)'];
        end
    end  
end
save(fnameFig,'CS','CSidx','nRuns2test','-append')










%% Plot figures

% Load results all runs
load(fnameFig);

% load results for PCA1
load(fnamePcca);

% Similarity and free energy for runs of PCCA1
SS=SSall(index,index);
F=Fall(index);

% Compute Between state correlation 
R = length(index);
Midx  = repelem(1:R,K);
M=nan(ndim1,R*K);
for r = 1:R
    M(:,Midx==r) = cell2mat(struct2cell(load([bnameHmm num2str(index(r)) '.mat'],'Gamma')));
end
CO = corrcoef(M);

%% Pannel HMM stability
figure(1);
% A) histogram HMM stability
subplot(5,1,1)
isubd = find(triu(ones(R,R),1));
h=histogram(SS(isubd));
xlabel('between-run similarity');ylabel('count')
% A) dendro HMM stability north
subplot(5,1,2)
Z=linkage(1-SS,'ward');
[~,~,orderR]=dendrogram(Z,0);
% A) matrix HMM stability
subplot(5,1,3)
imagesc(SS(orderR,orderR)), colorbar('east')
subtitle('between-run similarity');
% B) barplot Free Energy
subplot(5,1,4)
b=bar(F(orderR));hold on
[~,x]=min(b.YData);hold on;xline(x,'-','BR run');
ylabel('free energy'),xlabel([num2str(R) ' runs'])
% C) matrix between-state correlation
subplot(5,1,5)
Z=linkage(1-CO,'ward');
[~,~,orderK]=dendrogram(Z,0);
imagesc(CO(orderK,orderK));
ylabel([num2str(R) ' x ' num2str(K) ' states']);
xlabel([num2str(R) ' x ' num2str(K) ' states']);
subtitle({'between-state correlation'})
colormap(s8,col0);
s8.Position([1,3,4])=s3.Position([1,3,4]);
s8.Position(2)=diff(s5.Position([3,2]))*0.9;
cb=colorbar('east');cb.Position =  [0.2276    0.5165    0.0092    0.0317];

%% Pannel BR-HMM stability
figure(2);
colR = parula(length(nRuns2test));
% D1) dendro BR run similarities
subplot(2,3,1)
idxrr = unique(idxRR(R,:));
Nr=length(idxrr);
srr = SSall(idxrr,idxrr);
Z=linkage(1-srr,'ward');
[~,~, ordersrr]=dendrogram(Z,0);
% D2) matrix between-run similarity
subplot(2,3,4)
imagesc(srr(ordersrr,ordersrr)),colorbar('east')
xlabel([num2str(Nr) ' BR runs']),ylabel([num2str(Nr) ' BR runs'])
subtitle('between-run similarity')
% E) Matrix between-state similarity
subplot(2,3,2)
i = find(ismember(CCidx,idxrr));
crr = CC(i,i);
Z=linkage(1-crr,'ward');
[~,~,ordercrr]=dendrogram(Z,0);
subplot(2,3,5)
imagesc(crr(ordercrr,ordercrr));,colorbar('east')
xlabel([num2str(Nr) ' x ' num2str(K) ' states']);
ylabel([num2str(Nr) ' x ' num2str(K) ' states']);
subtitle('between-state correlation')

% F-G) BR-HMM stability (SS) as a function of free energy (FE)
for n = nRuns2test2 
    i= find(nRuns2test<=n);if isempty(i),i=1;else,i = i(end);end
    col = colR(i,:);
    idx  = unique(idxRR(n,:));
    nr = length(idx);isubd = find(triu(ones(nr,nr),1));    
    % SS
    srr = SSall(idx,idx);srr=srr(isubd);
    subplot(2,3,3);
    mean_srr = mean(srr);
    err_cc = std(srr)/sqrt(nr);
    errorbar(n,mean_srr,err_cc,'-','MarkerSize',2,'MarkerFaceColor',col,...
            'MarkerEdgeColor',col,'Color',col,'LineWidth',1);hold on
    ylabel('between-run similarity');

    % FE
    FE2  = Fall(idx);
    subplot(2,3,6);
    mean_FE2 = mean(FE2);
    err_FE2 = std(FE2)/sqrt(nr);
    errorbar(n,mean_FE2,err_FE2,'-','MarkerSize',2,'MarkerFaceColor',col,...
            'MarkerEdgeColor',col,'Color',col,'LineWidth',1);hold on
    ylabel('free energy'); 

end
xlabel('number of runs (R)');

%% Pannel PCCA stability
figure(3);
% G) maxtriX coeff PCA
subplot(1,3,1)
imagesc(coeff(orderK,:)),colorbar('east')
ylabel([num2str(R) ' x ' num2str(K) ' states']);
xlabel('principal components');
subtitle('PCA coefficients')
% H) barplot explained variance
subplot(1,3,2)
e_indiv = (e/sum(e))*100;
e_cummu = cumsum(e)/sum(e);e_cummu = e_cummu*100;
yyaxis right
plot(e_cummu(1:P),'.-');
ylabel('cummulative');
yyaxis left
bar(e_indiv(1:P))
ylabel('explained variance %'); 
xlabel('principal components')
% I) errorbar PCCA_stability
subplot(1,3,3)
x = linspace(1,P,P);
count=0;
for n = nRuns2test
    count=count+1;
    ind = find(CSidx(:,1)==n);
    mean_cc = mean(abs(CS(ind,1:P)));
    err_cc = std(abs(CS(ind,1:P)))/sqrt(Np);
    errorbar(x,mean_cc,err_cc,...
            '-o','MarkerSize',0.5,'MarkerFaceColor',colR(count,:),...
            'MarkerEdgeColor',colR(count,:),'Color',colR(count,:));hold on
end
ylabel({'between-PCCA-rep';'correlations'});
xlabel('principal components');
cb=colorbar();cb.Ticks=0.05:1/(length(nRuns2test)):1;
cb.TickLabels=nRuns2test;title(cb,'R');


%% Brain maps
wbdir    = 'usr/bin';
addpath(genpath(wbdir))
osl_startup

% Brain map PCCA
load(fnamePcca ,'Q');
mapfile = [fnamePcca '_map'];
makeMap(Q,[],parcfile,maskfile,2,1,0,mapfile,wbdir);

% Brain map BR-HMM
[~, br_run] = min(Fall);
fnameBrhmm = [bname0 num2str(br_run)];
mapfile = [fnameBrhmm  '_map'];
hmm = cell2mat(struct2cell(load([fnameBrhmm '.mat'],'hmm')));
makeMap(hmm,[],parcfile,maskfile,2,1,0,mapfile,wbdir);




