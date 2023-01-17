% Code used in Alonso et al. (2023) "Towards stability of dynamic FC
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
% Aarhus University, January 2023
%%%%%%%%%%%%%%%%%%%%%%%%%

% USER
R = 1000; % Number of HMM runs
K = 6;    % Number of HMM states per run
P = 20;   % Number of principal components

% Setup directories and filenames
mypath  = pwd;
wbdir   = '/usr/bin';
codedir = '/path/2/toolboxes';
addpath(genpath(wbdir))
addpath(genpath([codedir '/ohba_analysis/HMM-MAR']));
addpath(genpath([codedir '/ohba_analysis/osl-core']));osl_startup
addpath(genpath([codedir '/covariancetoolbox']));
hmmfile  = [mypath '/hmm_run'];     % Basename for HMM results 
pccafile = [mypath '/pcca'];   % Filename for PCCA results 
datafile = [mypath '/dataMEG.mat']; % File for data timeseries
parcfile = [codedir '/ohba_analysis/parcellations/fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm.nii.gz'];
maskfile = [codedir '/ohba_analysis/parcellations/MNI152_T1_8mm_brain'];;

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
%comment the following lines if analysing fMRI data
opt.onpower = 1;     % on power rather than raw time series
opt.filter = [8 12]; % alpha band
opt.Fs = Hz;

% Read data
load(datafile,'X','T','Hz');
ndim1=sum(cell2mat(T)); % total number of time points

%% Run HMM R times
for n = 1:R
    [hmm, Gamma,~,~,~,~,fehist] = hmmmar(X,T,opt);
    f = fehist(end);
    save([hmmfile num2str(n) '.mat'], 'hmm', 'Gamma','f'); 
end

%% HMM stability

% Compute Between-run similarity [S]
S=ones(R,R);
for ri=1:R-1
    gamma1 = cell2mat(struct2cell(load([hmmfile num2str(ri) '.mat'],'Gamma')));
    s=[];
    for rj=ri+1:R
        gamma2 = cell2mat(struct2cell(load([hmmfile num2str(rj) '.mat'],'Gamma')));
        s = [s,getGammaSimilarity(gamma1,gamma2)]; %#ok<AGROW> 
    end
    S(ri,ri+1:end) = s;
    S(ri+1:end,ri) = s;
end

% Compute Between-state correlation [CC]
M = nan(ndim1,R*K);
CCidx = repelem(1:R,K);
for r = 1:R
    gamma = cell2mat(struct2cell(load([hmmfile num2str(r) '.mat'],'Gamma')));
    M(:,CCidx == r) = gamma;
end 
CC = corrcoef(M);

% Figure
figure(1);
% histogram between-run similarities
subplot(411)
isubd = find(triu(ones(R,R),1));
h=histogram(S(isubd));
xlabel('between-run similarity');ylabel('count')
% dendrogram between-run similarities
subplot(412)
Z=linkage(1-S,'ward');
[~,~,order]=dendrogram(Z,0);
% matrix between-run similarities
subplot(413)
imagesc(S(order,order)), colorbar('east')
subtitle('between-run similarity');
% matrix between-state correlations
subplot(414)
Z=linkage(1-CC,'ward');
[~,~,orderM]=dendrogram(Z,0);
imagesc(CC(orderM,orderM)),colorbar('east')
ylabel([num2str(R) ' x ' num2str(K) ' states']);
xlabel([num2str(R) ' x ' num2str(K) ' states']);
subtitle({'between-state correlation'})

%% BR-HMM

% Load Free Energy [FE] of R runs
FE=nan(R,1);
for r=1:R
    FE(r)=cell2mat(struct2cell(load([hmmfile num2str(r) '.mat'],'f')));
end

% Figure
figure(2);
% Stability as a function of FE
[~,order]=sort(FE,'ascend');
subplot(211)
bar(FE(order))
ylabel('free energy');xlabel('sorted BR runs')
subplot(212)
Stmp = S; Stmp(eye(size(Stmp,1)))=nan;
imagesc(Stmp(order,order))
xlabel('sorted BR runs');ylabel('sorted BR runs')
set(gca,'XTick',[],'YTick',[])
cb=colorbar('east');title(cb,'S')

% Brain map BR-HMM
[~, br_run] = min(FE);
brhmmfile = [bname0 num2str(br_run)];
mapfile = [brhmmfile  '_map'];
hmm = cell2mat(struct2cell(load([brhmmfile '.mat'],'hmm')));
makeMap(hmm,[],parcfile,maskfile,2,1,0,mapfile,wbdir);


%% PCCA

% run PCA on the R x K states and compute PC covariance matrix (Q)
index=1:R;
[coeff, scores, e, Q] = run_pcca(hmmfile,index,K,P,ndim1,X);
save(pccafile,'coeff','scores','e','Q','index');

% figure
figure(3);
% matriX PCA coeff
subplot(1,3,1)
imagesc(coeff(orderM,:)),colorbar('east')
ylabel([num2str(R) ' x ' num2str(K) ' states']);
xlabel('principal components');
subtitle('PCA coefficients')
% barplot explained variance
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

% Brain map PCCA
mapfile = [pccafile '_map'];
makeMap(Q,[],parcfile,maskfile,2,1,0,mapfile,wbdir);





