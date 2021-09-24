function varargout=pp1_ana(what,varargin)
%% details

% analyses for passive patterns project
% published in:

% S.A.Arbuckle 2021

% requires the following toolboxes for out-of-the-box functionality:
%   dataframe toolbox: https://github.com/jdiedrichsen/dataframe
%   pattern component modeling toolbox: https://github.com/jdiedrichsen/pcm_toolbox
%   plotlib toolbox: https://github.com/nejaz1/plotlib

% Directories for data
dataDir = '/Users/sarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_analysis/pp1_encodingAnalyses/data';
%dataDir = '/home/saarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_analysis/pp1_encodingAnalyses/data';

subj_name = {'pd01','pd02','s01','s02','s03','s04','s05','s06','s07','s08','s09'}; 
%subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10'}; 

glm = 4;
sn  = 2:11;
roi = [1:6];

% TO DO:
% - re-number rois according to list below
% - re-name and re-number subjs according to list above
% - files that should be edited for this:
%       - behaviour
%       - region betas
%       - region G
%       - selectivity analysis
%       - re-run model analyses


% roi details:
%   #   |   name
% ---------------
%   1   | BA 4a (rostral M1) (left hemi)
%   2   | BA 4p (caudal M1) (left hemi)
%   3   | BA 3a (hand area) (left hemi)
%   4   | BA 3b (hand area) (left hemi)
%   5   | BA 1 (hand area) (left hemi)
%   6   | BA 2 (hand area) (left hemi)
%   7   | OP1 (left hemi)
%   8   | S1 (3a+3b+1+2)
%   9   | M1 (4a+4p)
roiPlotClrs= {[69 114 180]./255,...
    [116 173 209]./255,...
    [254 224 144]./255,...
    [253 174 97]./255,...
    [244 109 67]./255,...
    [215 48 39]./255,...
    [1 1 1]};
modelClrs={[0.12,0.3,0.58],...
    [0.8,0.15,0],...
    [0.88,0.51,0.36],...
    [0.96,0.80,0.32],...
    [0.12,0.3, 0.58],...
    [0 0 0],...
    [0 0 0]};
roiNames = {'3a','3b','1','2','4a','4p'};
roiPlotNum = [3 4 5 6 1 2 7];
roiPlotNames = {'4a','4p','3a','3b','1','2'};

switch what
    case '0' % behaviour analysis:
    case 'beha:performance'
        % calculate error rates
        T = load('/Users/sarbuckle/DATA/passivePatterns1/fmri/data/pp1_fmri_simpleAna.mat');%load(fullfile(dataDir,'all_behaviour.mat'));
        T.numTrials = ones(size(T.sn));
        T.mismatch = T.trialType==2;
        D=tapply(T,{'sn','trialType'},{'resp_CR','sum'},{'resp_FA','sum'},...
            {'resp_hit','sum'},{'resp_miss','sum'},{'numTrials','sum'});
        % calc error rates:
        D.prop_FA = D.resp_FA./D.numTrials;
        D.prop_CR = D.resp_CR./D.numTrials;
        D.prop_hit = D.resp_hit./D.numTrials;
        D.prop_miss = D.resp_miss./D.numTrials;
        
        % display error rates:
        prop_FA = D.prop_FA(D.trialType==1); % not mismatch
        prop_CR = D.prop_CR(D.trialType==1); % not mismatch
        prop_hit = D.prop_hit(D.trialType==2); % is mismatch
        prop_miss= D.prop_miss(D.trialType==2); % is mismatch
        prop_error = (D.resp_miss(D.trialType==2) + D.resp_FA(D.trialType==1)) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        
        tt_FA = D.resp_FA(D.trialType==1) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        tt_CR = D.resp_CR(D.trialType==1) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        tt_hit = D.resp_hit(D.trialType==2) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        tt_miss = D.resp_miss(D.trialType==2) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        
        tt_thumb = (D.resp_FA(D.trialType==1) + D.resp_hit(D.trialType==2)) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        
        % display to user
        fprintf('mean FALSE ALARM:       %1.2f +- %1.2f (of not mismatch trials) [%1.2f +- %1.2f of total trials]\n',mean(prop_FA)*100,stderr(prop_FA)*100, mean(tt_FA)*100,stderr(tt_FA)*100);
        fprintf('mean CORRECT REJECTION: %1.2f +- %1.2f (of not mismatch trials) [%1.2f +- %1.2f of total trials]\n',mean(prop_CR)*100,stderr(prop_FA)*100, mean(tt_CR)*100,stderr(tt_CR)*100);
        fprintf('mean HIT RATE:          %1.2f +- %1.2f (of     mismatch trials) [%1.2f +- %1.2f of total trials]\n',mean(prop_hit)*100,stderr(prop_hit)*100, mean(tt_hit)*100,stderr(tt_hit)*100);
        fprintf('mean MISS RATE:         %1.2f +- %1.2f (of     mismatch trials) [%1.2f +- %1.2f of total trials]\n',mean(prop_miss)*100,stderr(prop_miss)*100, mean(tt_miss)*100,stderr(tt_miss)*100);
        fprintf('mean ERROR rate (FA & misses): %1.2f +- %1.2f [of total trials]\n',mean(prop_error)*100,stderr(prop_error)*100);
        fprintf('mean THUMB PRESS rate: %1.2f +- %1.2f [of total trials]\n',mean(tt_thumb)*100,stderr(tt_thumb)*100);
        
        varargout = {D};
    
        
    case '0' % selectivity analysis:
    case 'selectivity:do_singfleFingers' 
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        sn  = 2:11;
        glm = 4;
        roi = [1:6];
        fthres = 0.95; % % cutoff value for f-crit (take top 5%)
        vararginoptions(varargin,{'glm','roi','sn','fthres'});
        
        % set random seed for reproducability
        rng(99)        
        
        numSim = 1000; % # simulated datasets per model per participant (random and selective tuning models)
        conds  = 1:5;  % conditions to analyze (single finger conditions)
        numConds = numel(conds);
        % load data
        T=[];
        for rr=roi
            t = load(fullfile(dataDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T = addstruct(T,t);
        end
        T=getrow(T,ismember(T.sn,sn));
        % do analysis
        fprintf('\nsubj\troi\tsig voxels (%%)\tsvar\tevar\tselectivity (normed)');
        fprintf('\n----\t---\t--------------\t----\t----\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b      = [];
            b.beta = T.betaUW{ii}; % univariately pre-whitened data 
            b.tt   = T.tt{ii};
            b.run  = T.run{ii};
            b      = getrow(b,ismember(b.tt,conds));
            
            % zero-centre the voxel tuning curves in each run
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            
            % do voxel sub-selection
            % restrict analyses to voxels with significant F-vals:
            [F,Fcrit]  = pp1_ana('selectivity:Ftest',b.beta,b.tt,b.run,fthres);
            numVoxOrig = size(b.beta,2);
            numVoxSig  = sum(F>=Fcrit);
            if numVoxSig==0 % if no voxels meet significance threshold, boot
                d.passive   = v;
                d.sn        = v.*T.sn(ii);
                d.roi       = v.*T.roi(ii);
                d.numVoxF   = v.*numVoxSig;
                d.numVoxTot = v.*size(T.raw_beta{ii},2);
                d.avgF      = v.*nan;
                d.glm       = v.*glm;
                d.sft       = [nan;nan;nan];
                d.sftProb   = [nan;nan;nan];
                d.isEV      = [0;1;2];
                D = addstruct(D,d);
                % display to user:c
                fprintf('%s\t%02d\t%2.3f\t\t%2.4f\t%2.4f\t%1.5f\n',subj_name{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,nan,nan,nan);
                continue
            end
            b.beta      = b.beta(:,F>=Fcrit); % drop non-sig. voxels
            
            % calculate signal and noise strengths:
            [evar,svar] = pp1_ana('selectivity:estVariances',b.beta,b.tt,b.run);
            if svar<0; svar = 0; end % negative signal is zero signal (not a common occurance)
            
            % calc selectivity of actual voxel data
            b = tapply(b,{'tt'},{'beta','mean'}); % avg. voxel tuning curves across runs (mean beta per voxel = 0)
            sftBeta = pp1_ana('selectivity:estSelectivity',b.beta);
            sftBeta = mean(sftBeta);
            
            % some simulation params:
            numVoxSim = ceil(numVoxSig/numConds)*numConds; % round up so equal # of voxels per condition (for sparse patterns)
            numRun    = numel(unique(T.run{ii}));
            Gtt       = corr(b.beta');%eye(numConds); % results for random tuning don't depend on finger covariance structure
            
            % calc expected selectivity under random tuning
            sftRand   = pp1_ana('selectivity:tuningRandom',evar,svar,Gtt,numVoxSim,numRun,numSim,fthres);
            sftRandEV = nanmean(sftRand);
            
            % calc expected selectivity under perfectly selective tuning (given noise) 
            sftSelect   = pp1_ana('selectivity:tuningSelective',evar,svar,numConds,numVoxSim,numRun,numSim,fthres);
            sftSelectEV = nanmean(sftSelect);
            
            % normalize selectivity of the data by the expected values of
            % random (0) and selective (1) tuning
            sft_norm = (sftBeta - sftRandEV) / (sftSelectEV - sftRandEV);
            
            % add to output structure
            d.sn        = T.sn(ii);
            d.roi       = T.roi(ii);
            d.glm       = glm;
            
            d.fthres    = fthres;
            d.fcrit     = Fcrit;
            d.numVoxF   = numVoxSig;
            d.numVoxTot = size(T.raw_beta{ii},2);
            d.avgF      = mean(F(F>=Fcrit));
            d.svar_est  = svar;
            d.evar_est  = evar;
            d.snr       = d.svar_est ./ d.evar_est;
            
            d.sft       = sftBeta;
            d.sft_Rand  = sftRandEV;
            d.sft_Select= sftSelectEV;
            d.sft_norm  = sft_norm;
            
            D = addstruct(D,d);
            % display to user:
            fprintf('%s\t%02d\t%2.3f\t\t%2.4f\t%2.4f\t%1.5f\n',subj_name{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,svar,evar,sft_norm);
        end
        save(fullfile(dataDir,sprintf('glm%d_selectivity_0604_fthres%1.2f.mat',glm,fthres)),'-struct','D');
        varargout = {D};   
    case 'selectivity:do_featurePatterns' 
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        sn  = 2:11;
        glm = 4;
        roi = [1:6];
        fthres = 0.95; % % cutoff value for f-crit (take top 5%)
        vararginoptions(varargin,{'glm','roi','sn','fthres'});
        
        % set random seed for reproducability
        %rng(99)        
        
        numSim = 1000; % # simulated datasets per model per participant (random and selective tuning models)
        conds  = 6:30; % interaction conditions
        numConds = numel(conds);
        % load data
        T=[];
        T=load(fullfile(dataDir,sprintf('glm%d_featurePatterns.mat',glm)));

        % do analysis
        fprintf('\nsubj\troi\tsig voxels (%%)\tsvar\tevar\tselectivity (normed)');
        fprintf('\n----\t---\t--------------\t----\t----\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            idx     = false(size(T.sn,1),1);
            idx(ii) = true;
            t = getrow(T,idx);
            U = t.U{1}(conds,:);
            
            % zero-centre the voxel tuning curves in each run
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            
            % do voxel sub-selection
            % restrict analyses to voxels with significant F-vals:
            [F,Fcrit]  = pp1_ana('selectivity:Ftest',b.beta,b.tt,b.run,fthres);
            numVoxOrig = size(b.beta,2);
            numVoxSig  = sum(F>=Fcrit);
            if numVoxSig==0 % if no voxels meet significance threshold, boot
                d.passive   = v;
                d.sn        = v.*T.sn(ii);
                d.roi       = v.*T.roi(ii);
                d.numVoxF   = v.*numVoxSig;
                d.numVoxTot = v.*size(T.raw_beta{ii},2);
                d.avgF      = v.*nan;
                d.glm       = v.*glm;
                d.sft       = [nan;nan;nan];
                d.sftProb   = [nan;nan;nan];
                d.isEV      = [0;1;2];
                D = addstruct(D,d);
                % display to user:c
                fprintf('%s\t%02d\t%2.3f\t\t%2.4f\t%2.4f\t%1.5f\n',subj_name{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,nan,nan,nan);
                continue
            end
            b.beta      = b.beta(:,F>=Fcrit); % drop non-sig. voxels
            
            % calculate signal and noise strengths:
            [evar,svar] = pp1_ana('selectivity:estVariances',b.beta,b.tt,b.run);
            if svar<0; svar = 0; end % negative signal is zero signal (not a common occurance)
            
            % calc selectivity of actual voxel data
            b = tapply(b,{'tt'},{'beta','mean'}); % avg. voxel tuning curves across runs (mean beta per voxel = 0)
            sftBeta = pp1_ana('selectivity:estSelectivity',b.beta);
            sftBeta = mean(sftBeta);
            
            % some simulation params:
            numVoxSim = ceil(numVoxSig/numConds)*numConds; % round up so equal # of voxels per condition (for sparse patterns)
            numRun    = numel(unique(T.run{ii}));
            Gtt       = corr(b.beta');%eye(numConds); % results for random tuning don't depend on finger covariance structure
            
            % calc expected selectivity under random tuning
            sftRand   = pp1_ana('selectivity:tuningRandom',evar,svar,Gtt,numVoxSim,numRun,numSim,fthres);
            sftRandEV = nanmean(sftRand);
            
            % calc expected selectivity under perfectly selective tuning (given noise) 
            sftSelect   = pp1_ana('selectivity:tuningSelective',evar,svar,numConds,numVoxSim,numRun,numSim,fthres);
            sftSelectEV = nanmean(sftSelect);
            
            % normalize selectivity of the data by the expected values of
            % random (0) and selective (1) tuning
            sft_norm = (sftBeta - sftRandEV) / (sftSelectEV - sftRandEV);
            
            % add to output structure
            d.sn        = T.sn(ii);
            d.roi       = T.roi(ii);
            d.glm       = glm;
            
            d.fthres    = fthres;
            d.fcrit     = Fcrit;
            d.numVoxF   = numVoxSig;
            d.numVoxTot = size(T.raw_beta{ii},2);
            d.avgF      = mean(F(F>=Fcrit));
            d.svar_est  = svar;
            d.evar_est  = evar;
            d.snr       = d.svar_est ./ d.evar_est;
            
            d.sft       = sftBeta;
            d.sft_Rand  = sftRandEV;
            d.sft_Select= sftSelectEV;
            d.sft_norm  = sft_norm;
            
            D = addstruct(D,d);
            % display to user:
            fprintf('%s\t%02d\t%2.3f\t\t%2.4f\t%2.4f\t%1.5f\n',subj_name{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,svar,evar,sft_norm);
        end
        save(fullfile(dataDir,sprintf('glm%d_selectivity_0604_fthres%1.2f.mat',glm,fthres)),'-struct','D');
        varargout = {D};   
    case 'selectivity:getTuningCurves'
        % Get the significant voxel tuning curves for plotting
        
        % separated the tuning curve harvesting to avoid bloat in the 'do'
        % case.
        
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = [1:6];
        vararginoptions(varargin,{'glm','roi','sn'});
        
        conds  = 1:5; % conditions to analyze (single finger conditions)
        % load data
        T=[];
        for rr=roi
            t = load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T=addstruct(T,t);
        end
        T=getrow(T,ismember(T.sn,sn));
        
        D = []; % output structure
        Q = [];
        v = ones(5,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.betaUW{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            
            % estimate signal and noise strengths
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            % restrict analyses to voxels with significant F-stats:
            [F,Fcrit]  = pp1_imana('SFT:calcFstat',b.beta,b.tt,b.run);
            sigIdx  = F>=Fcrit;
            if sum(sigIdx)==0
                keyboard
            end
            b.beta  = b.beta(:,sigIdx); % drop non-sig. voxels
            b = tapply(b,{'tt'},{'beta','mean'}); % avg. betas across runs (mean beta per voxel = 0)
            
            sftBeta = pp1_imana('SFT:estimateSFT',b.beta);
            
            % normalize values by lowest and max response
            [maxVal,maxIdx] = max(b.beta);
            [minVal,minIdx] = min(b.beta);
            b.betaNorm = b.beta-minVal;
            maxNormVal = max(b.betaNorm);
            b.betaNorm = b.betaNorm./maxNormVal;
            Cm = indicatorMatrix('identity',maxIdx');
            
            % check if we are missing any fingers:
            isMissing = ~ismember([1:5],unique(maxIdx));
            if sum(isMissing)>0
                Cm_old = Cm;
                Cm = zeros(size(Cm_old,1),5);
                Cm(:,~isMissing) = Cm_old;
            end
            d=[];
            d.tuningNorm = pinv(Cm)*b.betaNorm'; % each row is one finger
            d.tuningRaw  = pinv(Cm)*b.beta'; % each row is one finger
            d.avgOtherFingers = pinv(Cm) * ((sum(b.betaNorm)-1)./3)';
            d.avgSelectivity = pinv(Cm)*sftBeta';
            d.tuningNorm(isMissing,:) = nan;
            d.tuningRaw(isMissing,:) = nan;
            d.avgOtherFingers(isMissing,:) = nan;
            d.numVox     = sum(Cm)'; % # voxels averaged into each finger
            d.numVoxRoi  = v.*size(b.beta,2);
            d.svar  = v.*svar;
            d.evar  = v.*evar;
            d.digitMax = [1:5]';
            d.roi   = v.*T.roi(ii);
            d.sn    = v.*T.sn(ii);
            d.glm   = v.*glm;
            D = addstruct(D,d);
            
            % simpler structure:
            q.roi   = T.roi(ii);
            q.sn    = T.sn(ii);
            q.glm   = glm;
            q.numVoxF = sum(sigIdx);
            q.numVoxTot = size(b.beta,2);
            q.svar  = svar;
            q.evar  = evar;
            q.avgSFT = nanmean(pp1_imana('SFT:estimateSFT',b.betaNorm)); % avg. response to "less-preferred" fingers
            Q=addstruct(Q,q);
            
            
        end
        
        % make into more universal plotting structure:
        T=[];
        for ii=1:size(D.sn)
            t=[];
            t.tuningNorm = D.tuningNorm(ii,:)';
            t.tuningRaw  = D.tuningRaw(ii,:)';
            t.avgSelectivity = v.*D.avgSelectivity(ii);
            t.numVox     = v.*D.numVox(ii);
            t.numVoxRoi  = v.*D.numVoxRoi(ii);
            t.svar  = v.*D.svar(ii);
            t.evar  = v.*D.evar(ii);
            t.digitMax = v.*D.digitMax(ii);
            t.digit    = [1:5]';
            t.roi      = v.*D.roi(ii);
            t.glm      = v.*D.glm(ii);
            t.sn       = v.*D.sn(ii);
            T=addstruct(T,t);
        end

        varargout = {T,D,Q};      
    case 'selectivity:Ftest'
        % calculates F-statistic per voxel to determine if voxel is
        % significantly modulated by finger(s)
        Y  = varargin{1}; % N x P matrix of data. (N=numCond*numRun, P=numVox)
        cV = varargin{2}; % N x 1 vector of condition assignments
        pV = varargin{3}; % N x 1 vector of run assignments
        fthres = 0.95;
        if numel(varargin)>3
            fthres = varargin{4}; % percent cutoff for F-stat
        end
        
        % housekeeping
        numVox   = size(Y,2);
        conds    = unique(cV)';
        numCond  = numel(conds);
        runs     = unique(pV)';
        numRun   = numel(runs);
        df1 = numCond-1;
        df2 = numCond*numRun - numCond - numRun;
        % remove run means (to be safe- often the data will have run means
        % removed before calling this func)
        C0  = indicatorMatrix('identity',pV);
        Y   = Y - C0*pinv(C0)*Y; 
        % compute mean and covariance matrices
        muK = zeros(numCond,numVox); % condition means
        SSR = zeros(1,numVox);       % ssr vector (common across conditions)
        n   = zeros(1,numCond);      % # observations per condition
        for ii=1:numCond
            c = conds(ii);
            idx  = find(cV==c);
            n(ii)= numel(idx);
            muK(ii,:) = sum(Y(idx,:),1) ./ n(ii); % condition means
            res = bsxfun(@minus,Y(idx,:),muK(ii,:)); % residuals from the grand mean across all observations of this condition
            SSR = SSR + sum(res.^2,1) ./ n(ii); % SSR (across observations) scaled by number of observations (in case # obs differ per condition)
        end
        SSB = sum(muK.^2,1); 
        F   = (SSB./df1) ./ (SSR./df2);
        Fcrit = finv(fthres,df1,df2); % 95% cutoff for F-stat
        % calculate significance of F-stats:
        p=fcdf(F,df1,df2);
        eps=0.000001; 
        p(p>1-eps)=1-eps;
        p(p<eps)=eps; 

        varargout = {F,Fcrit,p};
    case 'selectivity:estVariances'
        % empirically estimates error variance in activity patterns across
        % runs. 
        % Estimate is more accurate when run means have been removed.

        % % NOTE: we integrate across conditions & voxels

        Y  = varargin{1};    % patterns   [regressors x voxels]
        cV = varargin{2};    % conditions [regressors x 1]
        pV = varargin{3};    % partitions [regressors x 1]

        nV = size(Y,2);         % # voxels
        nP = numel(unique(pV));  % # partitions
        nC = numel(unique(cV));  % # conditions
        Ya = zeros(nP,nV*nC);   % pre-allocate

        for pp = 1:nP
            y = Y(pV==pp,:);     % vectorize patterns across conditions per run
            Ya(pp,:) = y(:)';
        end

        take = logical(tril(ones(nP),-1)); % lower-triangular index
        G    = cov(Ya');  % covariances between runs (each row = one run, has zero mean)
        R    = corr(Ya'); % correlations between runs

        var_S = sum(G(take))/sum(sum(take)); % signal variance (avg. across runs)
        r     = sum(R(take))/sum(sum(take)); % signal correlation (avg. across runs)
        var_E = var_S/r - var_S;             % error variance (avg. across runs)

        varargout = {var_E,var_S,r};
    case 'selectivity:estSelectivity'
        % calculate single finger tuning using normalzied distance approach
        Y         = varargin{1}; % CxN matrix of data.
        numC      = size(Y,1);
        maxY      = max(Y,[],1);
        avgDistsN = (sum(maxY-Y)./(numC-1)) ./ (maxY-min(Y,[],1));
        varargout = {avgDistsN};
    case 'selectivity:tuningRandom'
        % calculates expected value (avg. sft across voxels) for voxels generated under specfified G

        % simulation params
        errvar  = varargin{1};
        sigvar  = varargin{2};
        G       = varargin{3}; % here we override this option and use group-average finger-by-finger correlation matrix from Ejaz (2015)
        numVox  = varargin{4};
        numRun  = varargin{5};
        numSim  = varargin{6};
        fthres  = varargin{7};
        
        % group-average finger-by-finger correlation matrix:
        %G = rsa_squareIPM([1,0.797,0.789,0.785,0.771,1,0.929,0.877,0.837,1,0.939,0.883,1,0.952,1]);

        % generate data
        D = pp1_ana('selectivity:simulateRandom',G,numVox,numRun,numSim,sigvar,errvar);
        % calc expected tuning on simulated datasets
        sft = nan(1,numSim);
        for s = 1:numSim
            d  = getrow(D,D.sn==s);
            % do voxel sub-selection
            C0  = indicatorMatrix('identity',d.run); 
            d.Y = d.Y -C0*pinv(C0)*d.Y; % remove run means
            [F,Fcrit] = pp1_ana('selectivity:Ftest',d.Y,d.cond,d.run,fthres);
            sigIdx = F>=Fcrit;
            %[evar,svar] = pp1_ana('selectivity:estVariances',d.Y,d.cond,d.run); % double-check simulated snr approximates perscribed snr
            % calc selectivity
            Cd = indicatorMatrix('identity',d.cond);
            Ysig = pinv(Cd)*d.Y(:,sigIdx); % avg. simulated tuning curves across runs for significant voxels
            tmp_sft = pp1_ana('selectivity:estSelectivity',Ysig);
            sft(s)  = mean(tmp_sft); % mean selectivity across significant voxels for this simulated dataset
        end
        varargout = {sft};
    case 'selectivity:tuningSelective'    
        % calculates expected value (avg. sft across voxels) for voxels
        % with sparse tuning but corrupted by some degree of noise

        % simulation params
        errvar  = varargin{1};
        sigvar  = varargin{2};
        numCond = varargin{3};
        numVox  = varargin{4};
        numRun  = varargin{5};
        numSim  = varargin{6};
        fthres  = varargin{7};

        % generate data
        D = pp1_ana('selectivity:simulateSelective',numCond,numVox,numRun,numSim,sigvar,errvar); % selective tuning patterns with noise
        % calc expected tuning on simulated datasets
        sft = nan(1,numSim);
        for s = 1:numSim
            d      = getrow(D,D.sn==s);
            % do voxel sub-selection
            C0     = indicatorMatrix('identity',d.run); 
            d.Y    = d.Y -C0*pinv(C0)*d.Y; % remove run means
            [F,Fcrit] = pp1_ana('selectivity:Ftest',d.Y,d.cond,d.run,fthres);
            sigIdx = F>=Fcrit;
            %[evar,svar] = pp1_ana('selectivity:estVariances',d.Y,d.cond,d.run); % double-check simulated snr approximates perscribed snr
            % calc selectivity
            Cd = indicatorMatrix('identity',d.cond);
            Ysig = pinv(Cd)*d.Y(:,sigIdx); % avg. simulated tuning curves across runs for significant voxels
            tmp_sft = pp1_ana('selectivity:estSelectivity',Ysig);
            sft(s)  = mean(tmp_sft); % mean selectivity across significant voxels for this simulated dataset
        end
        varargout = {sft};
    case 'selectivity:simulateRandom'
        % true betas drawn from normal distribution, with added i.i.d. noise
        G       = varargin{1}; % G is already double centred b/c run means removed
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5};
        noise   = varargin{6};

        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation

        numCond = size(G,1);
        % need to scale all elements of G by mean of the diagonal elements
        % (variances) to ensure appropriate signal scaling:
        G = G./ (sum(diag(G)) / (numCond-1));

        D = []; % output structure
        for s = 1:numSubj
            %U = mvnrnd_exact(G,numVox);

            pSignal    = unifrnd(0,1,numCond,numVox); 
            pNoise     = unifrnd(0,1,numCond*numRun,numVox); 
            % Generate true sparse patterns
            U = signalDist(pSignal); 
            E = (U*U'); 
            Z = E^(-0.5)*U;   % Make random orthonormal vectors
            A = pcm_diagonalize(G); 
            if (size(A,2)>numVox)
                error('not enough voxels to represent G'); 
            end
            trueU = A*Z(1:size(A,2),:)*sqrt(numVox); 
            trueU = bsxfun(@times,trueU,sqrt(signal));   % Multiply by signal scaling factor

            X = kron(ones(numRun,1),eye(numCond)); % design matrix

            % Now add the random noise 
            Noise  = noiseDist(pNoise); 
            Noise  = bsxfun(@times,Noise,sqrt(noise)); % Multiply by noise scaling factor
            d.Y    = X*trueU + Noise; % pull through condition specific patterns and add i.i.d. noise
            % indexing fields
            d.run   = kron([1:numRun],ones(1,numCond))';
            d.cond  = kron(ones(1,numRun),[1:numCond])';
            d.sn    = ones(size(d.run)).*s;
            D = addstruct(D,d);
        end
        varargout = {D};
    case 'selectivity:simulateSelective'
        % true patterns are sparse (0's and 1's), with iid noise
        numCond = varargin{1};% sparsity level (1=totally sparse, 2-two conditions tuned, etc.)
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5}; % signal variance
        noise   = varargin{6}; % noise variance
        
        % get # conditions based on sparsity level
        numVoxPerCond = ceil(numVox/numCond); % voxels tuned per chord
        numVox = numVoxPerCond*numCond;
        signal = signal*(numCond/1); % rescale the signal by # conditions (each condition contributes independent amount to signal)

        % define signal generation
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        D = []; % output structure
        for s = 1:numSubj % per simulated dataset
            % draw uniform values for signal and noise (to be inverted
            % through any arbitrary function later)
            pNoise = unifrnd(0,1,numCond*numRun,numVox); 
            % Generate true sparse patterns
            U = kron(eye(numCond),ones(1,numVoxPerCond)); % true patterns are 1s for voxels tuned to fingers, 0s for non-tuned
            % scale patterns to match specified signal strength
            U = bsxfun(@times,U,sqrt(signal));
            %U = U-mean(U);
            X = kron(ones(numRun,1),eye(numCond)); % design matrix
            trueU = X*U;
            % Now add the random noise 
            Noise  = noiseDist(pNoise); 
            Noise  = bsxfun(@times,Noise,sqrt(noise)); % sacle noise
            d.Y    = trueU + Noise;
            % indexing fields
            d.run   = kron([1:numRun],ones(1,numCond))';
            d.cond  = kron(ones(1,numRun),[1:numCond])';
            d.sn    = ones(size(d.run)).*s;
            D = addstruct(D,d);
        end
        varargout = {D};
        
        
    case '0' % representational model analysis:
    case 'model:do_regions'
        % fit representational encoding models to data from BA regions
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = [1:6];
        glm = 4;
        D = []; % output structure
        for rr=roi
            % load subject data (load in once for all subjs to save time):
            [Y,pV,cV] = pp1_ana('region:loadPatterns','sn',sn,'roi',rr,'glm',glm);
            Gm        = pp1_ana('region:loadG','roi',rr,'glm',glm);
            d         = pp1_ana('model:fitModels',Y,pV,cV,Gm);
            % add indexing fields to datastructure:
            v  = ones(size(d.dataset));
            d.roi = v.*rr;
            Cs = pcm_indicatorMatrix('identity',d.dataset);
            d  = rmfield(d,{'dataset'});
            d.sn = Cs*sn';
            D=addstruct(D,d);
        end
        % save output
        %save(fullfile(dataDir,sprintf('glm%d_repModelFits',glm)),'-struct','D');
        
        varargout = {D};
    case 'model:do_regions_active'
        % fit representational encoding models to data from BA regions
        sn  = 1:8; 
        roi = [1:6];
        glm = 1;
        D = []; % output structure
        for rr=roi
            % load subject data (load in once for all subjs to save time):
            [Y,pV,cV] = pp1_ana('region:loadPatterns_active','sn',sn,'roi',rr,'glm',glm);
            Gm        = pp1_ana('region:loadG_active','roi',rr,'glm',glm);
            d         = pp1_ana('model:fitModels',Y,pV,cV,Gm,'models',...
                {'null','linear','2finger','3finger','4finger','5finger','flexible','noiseCeiling'});
            % add indexing fields to datastructure:
            v  = ones(size(d.dataset));
            d.roi = v.*rr;
            Cs = pcm_indicatorMatrix('identity',d.dataset);
            d  = rmfield(d,{'dataset'});
            d.sn = Cs*sn';
            D=addstruct(D,d);
        end
        % save output
        %save(fullfile(dataDir,sprintf('glm%d_repModelFits_active',glm)),'-struct','D');
        
        varargout = {D};
    case 'model:do_tessels'    
    case 'model:do_featurePatterns'
        % fit representational encoding models to data from BA regions
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = [2];
        glm = 4;
        D = []; % output structure
        for rr=roi
            % load subject data (load in once for all subjs to save time):
            [Y,pV,cV,~,xyz] = pp1_ana('region:loadPatterns','sn',sn,'roi',rr,'glm',glm);
            Gm = pp1_ana('region:loadG','roi',rr,'glm',glm);
            for ii=1:numel(sn)
                fp = pp1_ana('model:estU_tikhonov',Y{ii},pV{ii},cV{ii},Gm,'singleFinger');
                d.model = 2;
                d.modelName = 'linear';
                d.U   = {fp};
                d.cV  = {[1:30]'};
                d.roi = rr;
                d.sn  = sn(ii);
                d.xyzcoord = {xyz{ii}};
    
                D=addstruct(D,d);
                
                d=[];
                fp = pp1_ana('model:estU_tikhonov',Y{ii},pV{ii},cV{ii},Gm,'condition');
                d.model = 1;
                d.modelName = 'noiseCeiling';
                d.U   = {fp};
                d.cV  = {[1:30]'};
                d.roi = rr;
                d.sn  = sn(ii);
                d.xyzcoord = {xyz{ii}};
    
                D=addstruct(D,d);
                
            end
%             df        = pp1_ana('model:getFeaturePatterns',Y,pV,cV,Gm,'models',{'4finger'});  
%             % add indexing fields to datastructure:
%             v  = ones(size(df.dataset));
%             df.roi = v.*rr;
%             Cs = pcm_indicatorMatrix('identity',df.dataset);
%             df  = rmfield(df,{'dataset'});
%             df.sn = Cs*sn';
%             df.xyzcoord = cell(size(df.model));
%             for ii=1:numel(unique(df.sn))
%                 for jj=unique(df.model)'
%                     idx = df.sn==sn(ii) & df.model==jj;
%                     df.xyzcoord{idx} = xyz{ii};
%                 end
%             end
%             
        end
        varargout = {D};
    case 'model:do_modelPatterns'
        % fit representational encoding models to data from BA regions
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = [4];
        glm = 4;
        D = []; % output structure
        for rr=roi
            % load subject data (load in once for all subjs to save time):
            [Y,pV,cV,~,xyz] = pp1_ana('region:loadPatterns','sn',sn,'roi',rr,'glm',glm);
            Gm = pp1_ana('region:loadG','roi',rr,'glm',glm);
            for ii=1:numel(sn)
                fp = pp1_ana('model:estU_tikhonov',Y{ii},pV{ii},cV{ii},Gm,'singleFinger');
                U_hat = pp1_ana('model:predictPatterns',fp,'linear',[]);
                d.model = 2;
                d.modelName = {'linear'};
                d.U   = {U_hat};
                d.cV  = {[1:31]'};
                d.roi = rr;
                d.sn  = sn(ii);
                d.xyzcoord = {xyz{ii}};
    
                D=addstruct(D,d);
                
                d=[];
                U_true = pp1_ana('model:estU_tikhonov',Y{ii},pV{ii},cV{ii},Gm,'condition');
                d.model = 1;
                d.modelName = {'noiseCeiling'};
                d.U   = {U_true};
                d.cV  = {[1:31]'};
                d.roi = rr;
                d.sn  = sn(ii);
                d.xyzcoord = {xyz{ii}};
    
                D=addstruct(D,d);
                
            end
%             df        = pp1_ana('model:getFeaturePatterns',Y,pV,cV,Gm,'models',{'4finger'});  
%             % add indexing fields to datastructure:
%             v  = ones(size(df.dataset));
%             df.roi = v.*rr;
%             Cs = pcm_indicatorMatrix('identity',df.dataset);
%             df  = rmfield(df,{'dataset'});
%             df.sn = Cs*sn';
%             df.xyzcoord = cell(size(df.model));
%             for ii=1:numel(unique(df.sn))
%                 for jj=unique(df.model)'
%                     idx = df.sn==sn(ii) & df.model==jj;
%                     df.xyzcoord{idx} = xyz{ii};
%                 end
%             end
%             
        end
        varargout = {D};
    case 'model:do_region_stats'
        % case to do statistical tests of model fits
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiNames = {'3a','3b','1','2','4a','4p','S2'};
        roiPlotNum = [3 4 5 6 1 2 7];
        roiPlotNames = {'4a','4p','3a','3b','1','2','S2'};
        
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_repModelFits.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % find the null and noise ceiling models:
        numModels  = numel(unique(T.model));
        modelNull  = unique(T.model(strcmp(T.modelName,'null')));
        modelNCeil = unique(T.model(strcmp(T.modelName,'noise ceiling')));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==modelNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==modelNCeil),ones(numModels,1));
        T.roiPlotNum = roiPlotNum(T.roi)';
        
        % find the other model ids:
        modelLinear  = unique(T.model(strcmp(T.modelName,'linear')));
        model2f      = unique(T.model(strcmp(T.modelName,'2finger')));
        model3f      = unique(T.model(strcmp(T.modelName,'3finger')));
        model4f      = unique(T.model(strcmp(T.modelName,'4finger')));
        modelFlex    = unique(T.model(strcmp(T.modelName,'flexible linear')));
        model2f_noAdjacent   = unique(T.model(strcmp(T.modelName,'2finger no adjacent')));
        model2f_onlyAdjacent = unique(T.model(strcmp(T.modelName,'2finger no non-adjacent')));
        
        % do stats:
        % make pivottables:
        fprintf('MEAN model fits:\n');
        pivottable(T.model,T.roiPlotNum,T.r_norm,'mean','subset',~ismember(T.model,[modelNull,modelNCeil]));
        fprintf('SEM model fits:\n');
        pivottable(T.model,T.roiPlotNum,T.r_norm,'stderr','subset',~ismember(T.model,[modelNull,modelNCeil]));

        % LINEAR MODEL RESULTS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- LINEAR MODEL FITS ------------------\n');
        fprintf('\n----------------------------------------\n');
        % compare normalized null vs. linear model fits:
        fprintf('normalized null vs. linear model fits:\n');
        idx = ismember(T.model,[modelNull,modelLinear]);
        stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',idx);
        stats.eff.p
        fprintf('\n----------------------------------------\n');
        % compare normalized linear model fit across regions:
        fprintf('normalized linear model fit across regions:\n');
        stats = anovaMixed(T.r_norm,T.sn,'within',[T.roiPlotNum],{'roi'},'subset',T.model==modelLinear);
        stats.eff.p
        fprintf('\n----------------------------------------\n');
        % compare linear model fit in BA 3b vs. other SMc regions:
        % (two sided paired t-test, bonferroni corrected alpha=0.01)
        fprintf('linear model fit in BA 3b vs. other SMc regions:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==4 & T.model==modelLinear),T.r_norm(T.roiPlotNum==rr & T.model==modelLinear),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % compare noise ceiling vs. linear model fits:
        fprintf('noise ceiling vs. linear model fit (2-sided paired t-tests):\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelNCeil),T.r_norm(T.roiPlotNum==rr & T.model==modelLinear),2,'paired');
            fprintf('\n')
        end


        % 2FINGER MODEL RESULTS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- 2FINGER MODEL FITS -----------------\n');
        fprintf('\n----------------------------------------\n');
        % compare normalized linear and 2f int model fits:
        fprintf('normalized linear vs. 2f int model fits:\n');
        idx = ismember(T.model,[modelLinear,model2f]);
        stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',idx);
        stats.eff.p
        fprintf('\n----------------------------------------\n');
        % compare noise ceiling vs. 2f model fits:
        fprintf('noise ceiling vs. 2f model fit (2-sided paired t-tests):\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelNCeil),T.r_norm(T.roiPlotNum==rr & T.model==model2f),2,'paired');
            fprintf('\n')
        end


        % ADJACENT vs NON-ADJACENT FINGERS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- ADJACENT VS. NON-ADJACENT FITS -----\n');
        fprintf('\n----------------------------------------\n');
        % ttest all vs. neighbouring interactions:
        fprintf('two-sided paired ttest NO ADJACENT PAIRS vs. ALL PAIRS model fits:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model2f_noAdjacent),T.r_norm(T.roiPlotNum==rr & T.model==model2f),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % ttest all vs. non-neighbouring interactions:
        fprintf('two-sided paired ttest NO NON-ADJACENT PAIRS vs. ALL PAIRS model fits:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model2f_onlyAdjacent),T.r_norm(T.roiPlotNum==rr & T.model==model2f),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % compare neighbouring vs. non-neighbouring model fits:
        fprintf('NO ADJACENT vs. NO NON-ADJACENT model fits:\n');
        idx = ismember(T.model,[model2f_noAdjacent,model2f_onlyAdjacent]);
        stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',idx);
        stats.eff.p
        fprintf('\n----------------------------------------\n');
        % ttest neighbouring vs. non-neighbouring model fits:
        fprintf('two-sided paired ttest NO ADJACENT vs. NO NON-ADJACENT model fits:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model2f_noAdjacent),T.r_norm(T.roiPlotNum==rr & T.model==model2f_onlyAdjacent),2,'paired');
            fprintf('\n')
        end


        % 3 & 4FINGER MODEL RESULTS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- 3FINGER & 4FINGER MODEL FITS -------\n');
        fprintf('\n----------------------------------------\n');
        % compare noise ceiling vs. 3f model fits:
        fprintf('noise ceiling vs. 3f model fit (2-sided paired t-tests):\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelNCeil),T.r_norm(T.roiPlotNum==rr & T.model==model3f),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % compare noise ceiling vs. 4f model fits:
        fprintf('noise ceiling vs. 4f model fit (2-sided paired t-tests):\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelNCeil),T.r_norm(T.roiPlotNum==rr & T.model==model4f),2,'paired');
            fprintf('\n')
        end


        % FLEXIBLE MODEL RESULTS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- FLEXIBLE MODEL FITS ----------------\n');
        fprintf('\n----------------------------------------\n');
        % compare flexible vs. linear model fits in each region:
        % (two sided paired t-test, bonferroni corrected alpha=0.0083)
        fprintf('flexible vs. linear model fits in each region:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelFlex),T.r_norm(T.roiPlotNum==rr & T.model==modelLinear),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % compare 2finger vs. flexible model fits in each region:
        % (two sided paired t-test, bonferroni corrected alpha=0.0083)
        fprintf('2finger vs. flexible model fits in each region:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model2f),T.r_norm(T.roiPlotNum==rr & T.model==modelFlex),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        D1=getrow(T,T.model==model2f);
        D2=getrow(T,T.model==modelFlex);
        D1.diff=D1.r_norm-D2.r_norm;
        fprintf('MEAN 2Finger-Flex fits:\n');
        pivottable([],D1.roiPlotNum,D1.diff,'mean')
        fprintf('SEM 2Finger-Flex fits:\n');
        pivottable([],D1.roiPlotNum,D1.diff,'stderr')
        fprintf('\n----------------------------------------\n');
        % compare 3finger vs. flexible model fits in each region:
        % (two sided paired t-test, bonferroni corrected alpha=0.0083)
        fprintf('3finger vs. flexible model fits in each region:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model3f),T.r_norm(T.roiPlotNum==rr & T.model==modelFlex),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        fprintf('4finger vs. flexible model fits in each region:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model4f),T.r_norm(T.roiPlotNum==rr & T.model==modelFlex),2,'paired');
            fprintf('\n')
        end
            
        varargout = {T};
    
    case 'model:fitModels'
        % case to fit models to participant data
        Y      = varargin{1}; % cell array of activity patterns
        pV     = varargin{2}; % cell array of partition assignment (for each row of Y)
        cV     = varargin{3}; % cell array of condition assignment (for each row of Y)
        Groi   = varargin{4}; % [31x31x #subj] model Gs (if 3rd dim==1, apply same G for all subjs)
        % define which models we are fitting:
        models = {'null','linear','2finger','3finger','4finger','flexible',...
            '2finger_noAdjacent','2finger_noNonAdjacent','noiseCeiling'};
        vararginoptions(varargin(5:end),{'models'});
        
        numModels = numel(models);
        
        % get chord matrix
        chords   = pp1_ana('misc:chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        
        % loop through tessels and fit each individually:
        D=[]; % output
        for s=1:numel(Y)
            if isempty(Y{s})
                % rare instance of no data in tesselated rois for
                % one participant
                fprintf('\n********** NO DATA **********');
                continue
            end
            % define model G:
            if size(Groi,3)==numel(Y)
                modelG = Groi(:,:,s); % different model prior for each dataset (tessellation analysis)
            else
                modelG = Groi; % same model prior for each dataset (ROI analysis)
            end
            % split data into partitions (for estimation of patterns under each model)
            % assign runs to each partition
            part = unique(pV{s});
            numPart = numel(part);
            partI = {};
            for ip=1:numPart
                partI{ip}=part(ip);
            end  
            % allocate space for predicted patterns, evaluation metrics:
            numVox = size(Y{s},2);
            G_pred = zeros(31,31,numModels);
            Y_avg  = nan([numModels,5,numPart]); % avg. activity per # digits
            modelTheta = {};
            % pre-allocate space for model fit evaluations
            SS1 = zeros(numPart,numModels);
            SS2 = SS1; SSC = SS1; RSS = SS1; TSS=SS1;
            SS1_train = SS1; SS2_train = SS1; SSC_train = SS1; 
            RSS_train=SS1; TSS_train=SS1;
            % loop through partitions and estimate patterns
            for ii=1:numel(partI)
                % estimate the true training condition activity patterns:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii}); 
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'condition');
                            
                % predict patterns under each model:
                for mm=1:numModels
                    modelName = models{mm};
                    fprintf('DATASET: %02d  |  CVFold: %d  |  Model: %s\n',s,ii,modelName); 
                    switch modelName
                        case 'null'
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_ana('model:predictPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = nan(1,4);  
                        case 'linear'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'singleFinger');
                            Ypred = pp1_ana('model:predictPatterns',U,'linear',[]);
                            thetaEst = nan(1,4);
                        case '2finger'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'2finger');
                            Ypred = pp1_ana('model:predictPatterns',U,'2finger',[]);
                            thetaEst = nan(1,4);
                        case '3finger'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'3finger');
                            Ypred = pp1_ana('model:predictPatterns',U,'3finger',[]);
                            thetaEst = nan(1,4);
                        case '4finger'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'4finger');
                            Ypred = pp1_ana('model:predictPatterns',U,'4finger',[]);
                            thetaEst = nan(1,4);
                        case '5finger'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'5finger');
                            Ypred = pp1_ana('model:predictPatterns',U,'5finger',[]);
                            thetaEst = nan(1,4);
                        case 'flexible'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'singleFinger');
                            theta0 = [log(0.9) log(0.8) log(0.7) log(0.6)];
                            thetaFcn = @(x) modelLossRSS(x,U,Utrain,'flexible'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));                        
                            Ypred = pp1_ana('model:predictPatterns',U,'flexible',thetaEst);
                        case 'flexiblePCM'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'singleFinger');
                            
                            thetaEst = pp1_ana('est_FlexibleParams',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            Ypred = pp1_ana('model:predictPatterns',U,'flexible',thetaEst);
                        case '2finger_noAdjacent'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'2finger_noAdjacent');
                            Ypred = pp1_ana('model:predictPatterns',U,'2finger_noAdjacent',[]);
                            thetaEst = nan(1,4);
                        case '2finger_noNonAdjacent'
                            [U,lambdaReg,thetaReg] = pp1_ana('model:estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG,'2finger_noNonAdjacent');
                            Ypred = pp1_ana('model:predictPatterns',U,'2finger_noNonAdjacent',[]);
                            thetaEst = nan(1,4);   
                        case 'noiseCeiling'
                            Ypred     = Utrain;
                            thetaReg  = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst  = nan(1,4);      
                        otherwise
                            error('no model named: %s',modelName)
                    end
                    modNames{mm,1} = modelName;
                    % calculate model predicted second moment
                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
                    % calculate model predicted avg. activity
                    Y_avg_cent(mm,:,ii)  = mean(numD_inv*(Ypred-mean(Ypred,1)),2)'; % avg. activity per # digits, centred at mean 0
                    Y_avg(mm,:,ii)       = mean(numD_inv*Ypred,2)'; % avg. activity per # digits
                    % calculate metrics for R and R2 of prediction against TRAINING data:
                    [SS1_train(ii,mm),SS2_train(ii,mm),SSC_train(ii,mm),RSS_train(ii,mm),TSS_train(ii,mm)] = pp1_ana('model:evaluate',Ypred,Utrain); % corr b/t pred and TRAINing patterns
                    % calculate metrics for R and R2 of prediction against TEST data:
                    [SS1(ii,mm), SS2(ii,mm), SSC(ii,mm), RSS(ii,mm), TSS(ii,mm)] = pp1_ana('model:evaluate',Ypred,Ytest); % corr b/t pred and TEST patterns
                end
            end
            G_pred = G_pred./numel(partI); % avg. model predicted Gs across folds
            d = [];
            % for each model, avg. thetas across folds:
            d.modelName  = modNames;
            d.modelTheta = cell2mat(cellfun(@(x) mean(x,1),modelTheta,'uni',0)');
            d.regTheta   = cellfun(@(x) mean(x,1),regTheta,'uni',0)';
            d.regLambda  = nanmean(regLambda,2);
            % Pearson's R:
            d.r_train = [mean(SSC_train./sqrt(SS1_train.*SS2_train))]'; % each row is one cv-fold, each column is a model
            d.r_test  = [mean(SSC./sqrt(SS1.*SS2))]';
            % R2:
            d.r2_train = [1-sum(RSS_train)./sum(TSS_train)]'; % each row is one cv-fold, each column is a model
            d.r2_test  = [1-sum(RSS)./sum(TSS)]';
            % arrange data into output structure:
            d.avgAct      = mean(Y_avg,3); % avg. activity per # digits
            d.avgAct_cent = mean(Y_avg_cent,3); % avg. activity per # digits
            for mm=1:numModels
                d.gpred(mm,:) = rsa_vectorizeIPM(G_pred(:,:,mm));
            end
            d.model   = [1:numModels]';
            d.dataset = ones(numModels,1).*s;
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    case 'model:evaluate'
        % Ypred and Ytest are assumed to be the same size [31xP]
        % Ypred and Ytest are assumed to be the same condition arrangement
        Ypred = varargin{1}; % predicted patterns [31xP]. Assume that condition 1 is row 1, etc..
        Ytest = varargin{2}; % test patterns. Should be [31xP]- leave-one-out evaluation
        
        Ypred = Ypred-mean(Ypred,1); % rmv voxel means
        Ytest = Ytest-mean(Ytest,1);
        
        % get metrics for correlation
        SS1 = sum(sum(Ytest.*Ytest)); % test SS
        SS2 = sum(sum(Ypred.*Ypred)); % pred SS
        SSC = sum(sum(Ypred.*Ytest)); % cov
        % get metrics for R2        
        RSS = sum(sum((Ytest-Ypred).^2));
        TSS = SS1;
        
        varargout = {SS1,SS2,SSC,RSS,TSS};        
    case 'model:estU_tikhonov'
        % Regularized regression estimate of feature patterns
        % Use pcm with fixed model G to estimate signal and noise parameters

        % inputs
        Y    = varargin{1}; % matrix of activity patterns [#conds*#runs x #vox]
        pV   = varargin{2}; % partition vector (assume cV and pV are same across subjs)
        cV   = varargin{3};
        G    = varargin{4}; % condition vector (chord #s)
        type = varargin{5}; % which features are we estimation
        
        % create feature matrix and estimate model prior G according to
        % 'type':
        switch type
            case 'condition'
                Z  = pcm_indicatorMatrix('identity',cV); % feature design matrix for activity patterns
                Z0 = eye(size(G,1));
                Gprior = pinv(Z0)*G*pinv(Z0)';
            case 'singleFinger'
                % create single finger feature matrix
                Z0 = pp1_ana('misc:chords');
                Gprior = pinv(Z0)*G*pinv(Z0)';
                Z  = kron(ones(numel(unique(pV)),1),Z0); % feature design matrix for activity patterns
            case '2finger'
                % create 2finger feature matrix
                c1 = pp1_ana('misc:chords'); % single finger features
                c2 = pp1_ana('misc:chord_pairs'); % finger pair features
                Z0 = [c1 c2];
                Gprior = pinv(Z0)*G*pinv(Z0)';
                Z = kron(ones(numel(unique(pV)),1),Z0); % feature design matrix for activity patterns  
            case '3finger'
                % create 3finger feature matrix
                c1 = pp1_ana('misc:chords'); % single finger features
                c2 = pp1_ana('misc:chord_pairs'); % finger pair features
                c3 = pp1_ana('misc:chord_triplets'); % finger triplet features
                Z0 = [c1 c2 c3];
                Gprior = pinv(Z0)*G*pinv(Z0)';
                Z = kron(ones(numel(unique(pV)),1),Z0); % feature design matrix for activity patterns
            case '4finger'
                % create 4finger feature matrix
                c1 = pp1_ana('misc:chords'); % single finger features
                c2 = pp1_ana('misc:chord_pairs'); % finger pair features
                c3 = pp1_ana('misc:chord_triplets'); % finger triplet features
                c4 = pp1_ana('misc:chord_quads'); % etc...
                Z0 = [c1 c2 c3 c4];
                Gprior = pinv(Z0)*G*pinv(Z0)';
                Z = kron(ones(numel(unique(pV)),1),Z0); % feature design matrix for activity patterns
            case '5finger'
                % create 4finger feature matrix
                c1 = pp1_ana('misc:chords'); % single finger features
                c2 = pp1_ana('misc:chord_pairs'); % finger pair features
                c3 = pp1_ana('misc:chord_triplets'); % finger triplet features
                c4 = pp1_ana('misc:chord_quads'); % etc...
                Z0 = [c1 c2 c3 c4 [zeros(30,1);1]]; % add 5-finger interaction
                Gprior = pinv(Z0)*G*pinv(Z0)';
                Z = kron(ones(numel(unique(pV)),1),Z0); % feature design matrix for activity patterns
            case '2finger_noAdjacent'
                % create finger feature matrix
                c1 = pp1_ana('misc:chords'); % single finger features
                c2 = pp1_ana('misc:chord_pairs_noAdjacent'); % finger pair features (for non-neighbouring chords)
                Z0 = [c1 c2];
                Gprior = pinv(Z0)*G*pinv(Z0)';
                Z = kron(ones(numel(unique(pV)),1),Z0); % feature design matrix for activity patterns
            case '2finger_noNonAdjacent'
                % create finger feature matrix
                c1 = pp1_ana('misc:chords'); % single finger features
                c2 = pp1_ana('misc:chord_pairs_noNonAdjacent'); % finger pair features (for non-neighbouring chords)
                Z0 = [c1 c2];
                Gprior = pinv(Z0)*G*pinv(Z0)';
                Z = kron(ones(numel(unique(pV)),1),Z0); % feature design matrix for activity patterns
            otherwise
                error('no feature model of this type')
        end
        
        M{1}.type = 'component';
        M{1}.Gc = Gprior;
        M{1}.numGparams = 1;
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,Z,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        U = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);

        lambda = exp(theta_hat{1}(2))/exp(theta_hat{1}(1)); %lambda is noise/scale
        
        varargout = {U,lambda,theta_hat{1}};    
    case 'model:predictPatterns'
        % factorization of encoding models:
        U     = varargin{1}; % feature patterns
        model = varargin{2};
        theta = varargin{3}; % model params (only needed for flexible model)
        
        switch model
            case 'null'
                % Model predicts overall scaling of avg. activity with #
                % fingers. Scaling matches true mean scaling in training
                % data.
                % Ysf here are all 31 conditions from the training data
                % Set each condition pattern to the be mean pattern for all
                % chords with the same number of fingers.
                chords = pp1_ana('misc:chords'); 
                X = pcm_indicatorMatrix('identity',sum(chords,2)); % which patterns have the same # of fingers?
                Ymf_hat = X*pinv(X)*U;
            case 'linear'
                % get multi-finger design:
                X = pp1_ana('misc:chords');
                Ymf_hat = X*U;
            case 'flexible'
                % theta(1:4) = finger combination param (per # fingers in
                % chords for 2:5 digits)
                
                X = pp1_ana('misc:chords');
                numD = sum(X,2);
                X = X.*[ones(1,5) exp(theta(numD(numD>1)-1))]'; % flexible scaling per # fingers stimulated (force positive values with exp)
                Ymf_hat = X*U; 
            case '2finger'
                % model that includes 2-finger interaction components
                
                X1 = pp1_ana('misc:chords');
                X2 = pp1_ana('misc:chord_pairs');
                X  = [X1 X2];
                Ymf_hat = X*U; 
            case '3finger'
                % model that includes 2-finger interaction components
                
                X1 = pp1_ana('misc:chords');
                X2 = pp1_ana('misc:chord_pairs');
                X3 = pp1_ana('misc:chord_triplets');
                X  = [X1 X2 X3];
                Ymf_hat = X*U; 
            case '4finger'
                X1 = pp1_ana('misc:chords');
                X2 = pp1_ana('misc:chord_pairs');
                X3 = pp1_ana('misc:chord_triplets');
                X4 = pp1_ana('misc:chord_quads');
                X  = [X1 X2 X3 X4];
                Ymf_hat = X*U; 
            case '5finger'
                X1 = pp1_ana('misc:chords');
                X2 = pp1_ana('misc:chord_pairs');
                X3 = pp1_ana('misc:chord_triplets');
                X4 = pp1_ana('misc:chord_quads');
                X  = [X1 X2 X3 X4 [zeros(30,1);1]];
                Ymf_hat = X*U; 
            case '2finger_noAdjacent'
                X1 = pp1_ana('misc:chords');
                X2 = pp1_ana('misc:chord_pairs_noAdjacent');
                X  = [X1 X2];
                Ymf_hat = X*U; 
            case '2finger_noNonAdjacent'
                X1 = pp1_ana('misc:chords');
                X2 = pp1_ana('misc:chord_pairs_noNonAdjacent');
                X  = [X1 X2];
                Ymf_hat = X*U;                 
        end
        varargout = {Ymf_hat};
    
        
    case '0' % plotting:   
    case 'plot_encodingFits'
        % plots the normalized representational model fits     
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_repModelFits.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % find the null and noise ceiling models:
        numModels  = numel(unique(T.model));
        modelNull  = unique(T.model(strcmp(T.modelName,'null')));
        modelNCeil = unique(T.model(strcmp(T.modelName,'noise ceiling')));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==modelNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==modelNCeil),ones(numModels,1));
        T.roiPlotNum = roiPlotNum(T.roi)';

        % nice plot
        sty = style.custom(modelClrs); 
        sty.general.linestyle = {'-','-','-','-','-.','-','-'};
        sty.general.markerfill = [{modelClrs{1:6}},{[1 1 1]}];
        
        plt.line(T.roiPlotNum,T.r_norm,'split',T.model,'plotfcn','mean','style',sty,'subset',~ismember(T.model,[modelNull,modelNCeil]));
        ylabel(sprintf('normalized model fits\n(Pearson''s R)'));
        xlabel('Brodmann area');
        title('encoding model fits - passive');
        set(gca,'xticklabel',roiPlotNames);
                
        varargout = {T};
    case 'plot_encodingFits_active'
        % plots the normalized representational model fits     
        glm = 1;
        sn = 1:8;
        roi = [1:6];
        
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_repModelFits_active.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % find the null and noise ceiling models:
        numModels  = numel(unique(T.model));
        modelNull  = unique(T.model(strcmp(T.modelName,'null')));
        modelNCeil = unique(T.model(strcmp(T.modelName,'noiseCeiling')));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==modelNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==modelNCeil),ones(numModels,1));
        T.roiPlotNum = roiPlotNum(T.roi)';

        % nice plot
        sty = style.custom(modelClrs); 
        sty.general.linestyle = {'-','-','-','-','-.','-','-'};
        sty.general.markerfill = [{modelClrs{1:6}},{[1 1 1]}];
        
        plt.line(T.roiPlotNum,T.r_norm,'split',T.model,'plotfcn','median','style',sty,'subset',~ismember(T.model,[modelNull,modelNCeil]));
        ylabel(sprintf('normalized model fits\n(Pearson''s R)'));
        xlabel('Brodmann area');
        title('encoding model fits - passive');
        set(gca,'xticklabel',roiPlotNames);
                
        varargout = {T};
    
        
    case '0' % statistical tests for paper:
    case 'stats_encodingFits'
        % case to do statistical tests of model fits
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        
        % load data:
       % T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_Final.mat',glm)));
        T = load(fullfile(dataDir,sprintf('glm%d_repModelFits.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % find the null and noise ceiling models:
        numModels  = numel(unique(T.model));
        modelNull  = unique(T.model(strcmp(T.modelName,'null')));
        modelNCeil = unique(T.model(strcmp(T.modelName,'noise ceiling')));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==modelNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==modelNCeil),ones(numModels,1));
        T.roiPlotNum = roiPlotNum(T.roi)';
        
        % find the other model ids:
        modelLinear  = unique(T.model(strcmp(T.modelName,'linear')));
        model2f      = unique(T.model(strcmp(T.modelName,'2finger')));
        model3f      = unique(T.model(strcmp(T.modelName,'3finger')));
        model4f      = unique(T.model(strcmp(T.modelName,'4finger')));
        modelFlex    = unique(T.model(strcmp(T.modelName,'flexible linear')));
        model2f_noAdjacent   = unique(T.model(strcmp(T.modelName,'2finger no adjacent')));
        model2f_onlyAdjacent = unique(T.model(strcmp(T.modelName,'2finger no non-adjacent')));
        
        % do stats:
        % make pivottables:
        fprintf('MEAN model fits:\n');
        pivottable(T.model,T.roiPlotNum,T.r_norm,'mean','subset',~ismember(T.model,[modelNull,modelNCeil]));
        fprintf('SEM model fits:\n');
        pivottable(T.model,T.roiPlotNum,T.r_norm,'stderr','subset',~ismember(T.model,[modelNull,modelNCeil]));

        % LINEAR MODEL RESULTS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- LINEAR MODEL FITS ------------------\n');
        fprintf('\n----------------------------------------\n');
        % compare normalized null vs. linear model fits:
        fprintf('normalized null vs. linear model fits:\n');
        idx = ismember(T.model,[modelNull,modelLinear]);
        stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',idx);
        stats.eff.p
        fprintf('\n----------------------------------------\n');
        % compare normalized linear model fit across regions:
        fprintf('normalized linear model fit across regions:\n');
        stats = anovaMixed(T.r_norm,T.sn,'within',[T.roiPlotNum],{'roi'},'subset',T.model==modelLinear);
        stats.eff.p
        fprintf('\n----------------------------------------\n');
        % compare linear model fit in BA 3b vs. other SMc regions:
        % (two sided paired t-test, bonferroni corrected alpha=0.01)
        fprintf('linear model fit in BA 3b vs. other SMc regions:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==4 & T.model==modelLinear),T.r_norm(T.roiPlotNum==rr & T.model==modelLinear),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % compare noise ceiling vs. linear model fits:
        fprintf('noise ceiling vs. linear model fit (2-sided paired t-tests):\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelNCeil),T.r_norm(T.roiPlotNum==rr & T.model==modelLinear),2,'paired');
            fprintf('\n')
        end


        % 2FINGER MODEL RESULTS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- 2FINGER MODEL FITS -----------------\n');
        fprintf('\n----------------------------------------\n');
        % compare normalized linear and 2f int model fits:
        fprintf('normalized linear vs. 2f int model fits:\n');
        idx = ismember(T.model,[modelLinear,model2f]);
        stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',idx);
        stats.eff.p
        fprintf('\n----------------------------------------\n');
        % compare noise ceiling vs. 2f model fits:
        fprintf('noise ceiling vs. 2f model fit (2-sided paired t-tests):\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelNCeil),T.r_norm(T.roiPlotNum==rr & T.model==model2f),2,'paired');
            fprintf('\n')
        end


        % ADJACENT vs NON-ADJACENT FINGERS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- ADJACENT VS. NON-ADJACENT FITS -----\n');
        fprintf('\n----------------------------------------\n');
        % ttest all vs. neighbouring interactions:
        fprintf('two-sided paired ttest NO ADJACENT PAIRS vs. ALL PAIRS model fits:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model2f_noAdjacent),T.r_norm(T.roiPlotNum==rr & T.model==model2f),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % ttest all vs. non-neighbouring interactions:
        fprintf('two-sided paired ttest NO NON-ADJACENT PAIRS vs. ALL PAIRS model fits:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model2f_onlyAdjacent),T.r_norm(T.roiPlotNum==rr & T.model==model2f),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % compare neighbouring vs. non-neighbouring model fits:
        fprintf('NO ADJACENT vs. NO NON-ADJACENT model fits:\n');
        idx = ismember(T.model,[model2f_noAdjacent,model2f_onlyAdjacent]);
        stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',idx);
        stats.eff.p
        fprintf('\n----------------------------------------\n');
        % ttest neighbouring vs. non-neighbouring model fits:
        fprintf('two-sided paired ttest NO ADJACENT vs. NO NON-ADJACENT model fits:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model2f_noAdjacent),T.r_norm(T.roiPlotNum==rr & T.model==model2f_onlyAdjacent),2,'paired');
            fprintf('\n')
        end


        % 3 & 4FINGER MODEL RESULTS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- 3FINGER & 4FINGER MODEL FITS -------\n');
        fprintf('\n----------------------------------------\n');
        % compare noise ceiling vs. 3f model fits:
        fprintf('noise ceiling vs. 3f model fit (2-sided paired t-tests):\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelNCeil),T.r_norm(T.roiPlotNum==rr & T.model==model3f),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % compare noise ceiling vs. 4f model fits:
        fprintf('noise ceiling vs. 4f model fit (2-sided paired t-tests):\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelNCeil),T.r_norm(T.roiPlotNum==rr & T.model==model4f),2,'paired');
            fprintf('\n')
        end


        % FLEXIBLE MODEL RESULTS
        fprintf('\n----------------------------------------\n');
        fprintf('\n--- FLEXIBLE MODEL FITS ----------------\n');
        fprintf('\n----------------------------------------\n');
        % compare flexible vs. linear model fits in each region:
        % (two sided paired t-test, bonferroni corrected alpha=0.0083)
        fprintf('flexible vs. linear model fits in each region:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==modelFlex),T.r_norm(T.roiPlotNum==rr & T.model==modelLinear),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        % compare 2finger vs. flexible model fits in each region:
        % (two sided paired t-test, bonferroni corrected alpha=0.0083)
        fprintf('2finger vs. flexible model fits in each region:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model2f),T.r_norm(T.roiPlotNum==rr & T.model==modelFlex),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        D1=getrow(T,T.model==model2f);
        D2=getrow(T,T.model==modelFlex);
        D1.diff=D1.r_norm-D2.r_norm;
        fprintf('MEAN 2Finger-Flex fits:\n');
        pivottable([],D1.roiPlotNum,D1.diff,'mean')
        fprintf('SEM 2Finger-Flex fits:\n');
        pivottable([],D1.roiPlotNum,D1.diff,'stderr')
        fprintf('\n----------------------------------------\n');
        % compare 3finger vs. flexible model fits in each region:
        % (two sided paired t-test, bonferroni corrected alpha=0.0083)
        fprintf('3finger vs. flexible model fits in each region:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model3f),T.r_norm(T.roiPlotNum==rr & T.model==modelFlex),2,'paired');
            fprintf('\n')
        end
        fprintf('\n----------------------------------------\n');
        fprintf('4finger vs. flexible model fits in each region:\n');
        for rr=1:6
            fprintf('ROI %s\n',roiPlotNames{rr});
            ttest(T.r_norm(T.roiPlotNum==rr & T.model==model4f),T.r_norm(T.roiPlotNum==rr & T.model==modelFlex),2,'paired');
            fprintf('\n')
        end
            
        varargout = {T};
    
        
    case '0' % miscellanous cases (used for representational model analysis):    
    case 'region:loadPatterns'
        % Get betas for roi from subjects in PCM-friendly format.
        % Run means are NOT removed as this is a feature we want to retain.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi','betaType'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        % load betas
        betaType = 'betaW'; % multivariately prewhitened (or betaUW, raw_beta)
        B = load(fullfile(dataDir,sprintf('glm%d_roi%d_betas.mat',glm,roi)));
        B = getrow(B,B.roi==roi);
        % outputs
        Y = {};
        partVec = {};
        condVec = {};
        xyzCoord = {};
        for ii = 1:length(sn)
            % get subject data
            s = sn(ii);
            b = getrow(B,B.sn==s);
            bb = [];
            bb.run   = cell2mat(b.run);
            bb.chord = cell2mat(b.tt);
            eval(sprintf('bb.betas = cell2mat(b.%s);',betaType));
            bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive stimulation conditions (tt 32 is thumb press)
            % put subj data into pcm variables
            Y{ii}         = bb.betas;
            partVec{ii}   = bb.run;
            condVec{ii}   = bb.chord;
            Gcv(:,:,ii) = pcm_estGCrossval(Y{ii},partVec{ii},condVec{ii});
            Gcv(:,:,ii) = pcm_makePD(Gcv(:,:,ii)); % make subject G PD
            xyzCoord{ii}= b.xyzcoord{1};
        end
        varargout = {Y,partVec,condVec,Gcv,xyzCoord};     
    case 'region:makeG'
        % Estimate Region G (G is avg. semi-positive definite crossval G across
        % participants from roi):
        sn  = 2:11;
        glm = 4;
        roi = [1:9,12,15:21,29:30];
        D =[];
        for rr=roi
            G_pd = [];
            % load betas
            betaType = 'betaW'; % multivariately prewhitened (or betaUW, raw_beta)
            B = load(fullfile(dataDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            for ii = 1:length(sn)
                % get subject data
                s = sn(ii);
                b = getrow(B,B.sn==s & B.roi==rr);
                bb = [];
                bb.run   = cell2mat(b.run);
                bb.chord = cell2mat(b.tt);
                eval(sprintf('bb.betas = cell2mat(b.%s);',betaType));
                bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive stimulation conditions (tt 32 is thumb press)
                % calc subject G
                G_subj = pcm_estGCrossval(bb.betas,bb.run,bb.chord);
                G_pd(:,:,ii) = pcm_makePD(G_subj); % make subject G PD
            end
            G_pd = mean(G_pd,3);
            % add to output struct:
            d.roi=rr;
            d.g=rsa_vectorizeIPM(G_pd);
            D=addstruct(D,d);
        end   
       % save(fullfile(dataDir,sprintf('glm%d_regionG.mat',glm)),'-struct','D');   
       varargout = {D};
    case 'region:loadG'
        % Get Region G (G is avg. semi-positive definite crossval G across
        % participants from roi). Participants 2:11 are included in
        % estimate
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        
        D = load(fullfile(dataDir,sprintf('glm%d_regionG.mat',glm)));
        D = getrow(D,D.roi==roi);
        G = rsa_squareIPM(D.g);
        varargout={G}; 
    
    case 'region:loadPatterns_active'
        % Get betas for roi from subjects in PCM-friendly format.
        % Run means are NOT removed as this is a feature we want to retain.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        % load betas
        betaType = 'betaW'; % multivariately prewhitened (or betaUW, raw_beta)
        B = load(fullfile(dataDir,sprintf('cpd_glm%d_roi%d_betas.mat',glm,roi)));
        B = getrow(B,B.roi==roi);
        % outputs
        Y = {};
        partVec = {};
        condVec = {};
        for ii = 1:length(sn)
            % get subject data
            s = sn(ii);
            b = getrow(B,B.sn==s);
            bb = [];
            bb.sess  = cell2mat(b.sess);
            bb.run   = cell2mat(b.run) + (bb.sess-1)*8; % change run #s to reflect different sessions
            bb.chord = cell2mat(b.tt);
            eval(sprintf('bb.betas = cell2mat(b.%s);',betaType));
            bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive stimulation conditions (tt 32 is thumb press)
            % put subj data into pcm variables
            Y{ii}         = bb.betas;
            partVec{ii}   = bb.run;
            condVec{ii}   = bb.chord;
            Gcv(:,:,ii) = pcm_estGCrossval(Y{ii},partVec{ii},condVec{ii});
            Gcv(:,:,ii) = pcm_makePD(Gcv(:,:,ii)); % make subject G PD
        end
        varargout = {Y,partVec,condVec,Gcv};     
    case 'region:makeG_active'
        % Estimate Region G (G is avg. semi-positive definite crossval G across
        % participants from roi):
        sn  = 1:8;
        glm = 1;
        roi = 1:9;
        D =[];
        for rr=roi
            G_pd = [];
            % load betas
            betaType = 'betaW'; % multivariately prewhitened (or betaUW, raw_beta)
            B = load(fullfile(dataDir,sprintf('cpd_glm%d_roi%d_betas.mat',glm,rr)));
            for ii = 1:length(sn)
                % get subject data
                s = sn(ii);
                b = getrow(B,B.sn==s & B.roi==rr);
                bb = [];
                bb.run   = cell2mat(b.run);
                bb.chord = cell2mat(b.tt);
                eval(sprintf('bb.betas = cell2mat(b.%s);',betaType));
                bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive stimulation conditions (tt 32 is thumb press)
                % calc subject G
                G_subj = pcm_estGCrossval(bb.betas,bb.run,bb.chord);
                G_pd(:,:,ii) = pcm_makePD(G_subj); % make subject G PD
            end
            G_pd = mean(G_pd,3);
            % add to output struct:
            d.roi=rr;
            d.g=rsa_vectorizeIPM(G_pd);
            D=addstruct(D,d);
        end   
        save(fullfile(dataDir,sprintf('cpd_glm%d_regionG.mat',glm)),'-struct','D');   
    case 'region:loadG_active'
        % Get Region G (G is avg. semi-positive definite crossval G across
        % participants from roi). Participants 2:11 are included in
        % estimate
        glm = 1;
        roi = []; % only one roi supported
        vararginoptions(varargin,{'roi','glm'});
        
        D = load(fullfile(dataDir,sprintf('cpd_glm%d_regionG.mat',glm)));
        D = getrow(D,D.roi==roi);
        G = rsa_squareIPM(D.g);
        varargout={G};  
    
    case 'misc:chords'
        % returns indicator matrix for chords used in exp.
        chords = [eye(5);...                                                       % singles            5
             1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;...                        % doubles (thumb)    4
             0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;...                                   % doubles            3
             0 0 1 1 0; 0 0 1 0 1;...                                              % doubles            2
             0 0 0 1 1;...                                                         % doubles            1
             1 1 1 0 0; 1 1 0 1 0; 1 1 0 0 1; 1 0 1 1 0; 1 0 1 0 1; 1 0 0 1 1;...  % triples (thumb)    6
             0 1 1 1 0; 0 1 1 0 1; 0 1 0 1 1; 0 0 1 1 1;...                        % triples            4
             1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 0 1 1 1; 0 1 1 1 1;                % quadruples         5
             1 1 1 1 1];                                                           % all five           1
                                                                                   % total-------------31
        chord_strings = {'1','2','3','4','5',...
            '12','13','14','15','23','24','25','34','35','45',...
            '123','124','125','134','135','145','234','235','245','345',...
            '1234','1235','1245','1345','2345','12345'}; 
        varargout = {chords,chord_strings}; 
    case 'misc:chord_pairs'
        % returns indicator matrix for pairs of fingers use in each config:
        X=pp1_ana('misc:chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==2,:); % 2 finger pairs
        Xp = zeros(31,10);
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==2;
            Xp(pidx,ii) = 1;
        end
        varargout = {Xp};
    case 'misc:chord_triplets'
        % returns indicator matrix for sets of 3 fingers use in each config:
        X=pp1_ana('misc:chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==3,:); % 3 finger triplets
        Xp = zeros(31,10);
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==3;
            Xp(pidx,ii) = 1;
        end
        varargout = {Xp};
    case 'misc:chord_quads'
        % returns indicator matrix for set of 4 fingers use in each config:
        X=pp1_ana('misc:chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==4,:);
        Xp = zeros(31,5);
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==4;
            Xp(pidx,ii) = 1;
        end
        varargout = {Xp};
    case 'misc:chord_pairs_noAdjacent'
        % chord pairs, kicking out immediate neighours
        X=pp1_ana('misc:chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==2,:); % 2 finger pairs
        nIdx = sum(diff(pairs,1,2)==0,2); % which pairs are neighbours?
        pairs = pairs(~nIdx,:); % drop neighbour pairs
        Xp = zeros(31,size(pairs,1));
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==2;
            Xp(pidx,ii) = 1;
        end
        varargout = {Xp};
    case 'misc:chord_pairs_noNonAdjacent'
        % chord pairs, kicking out non-immediate neighour pairs
        X=pp1_ana('misc:chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==2,:); % 2 finger pairs
        nIdx = sum(diff(pairs,1,2)==0,2); % which pairs are neighbours?
        pairs = pairs(nIdx==1,:); % drop non-neighbour pairs
        Xp = zeros(31,size(pairs,1));
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==2;
            Xp(pidx,ii) = 1;
        end
        varargout = {Xp};
    
        
    case '0' % sanity check:
    case 'sim:testModelFitStability'
        % case to simulate stability of model fit normalization across different signal levels
        signal = [0.1:0.1:1];
        noise  = 1;
        numSim = 1000;
        roi    = 3; % what roi do we use for reference when simulating data?
        G = pp1_ana('region:loadG','roi',roi,'glm',4);
        D = [];
        for ii=signal
            fprintf('signal %1.2f\n',ii)
            % simulate data
            [Y,pV,cV] = pp1_ana('sim:Data',ii,noise,G,numSim);
            % fit encoding models:
            d = pp1_ana('model:fitModels',Y,pV,cV,G,'models',{'null','linear','flexible','2finger','noiseCeiling'});
            d.signal = ones(size(d.dataset)).*ii;
            d.noise  = ones(size(d.dataset)).*noise;
            D = addstruct(D,d);
        end
        
        D.r_norm = D.r_test - kron(D.r_test(D.model==1),ones(5,1));
        D.r_norm = D.r_norm./kron(D.r_norm(D.model==5),ones(5,1));
        
        
        sty = style.custom(plt.helper.get_shades(6,'hot'));
        plt.line(D.signal,D.r_norm,'split',D.model,'style',sty);
        xlabel('signal strength');
        ylabel('normalized model fit (Pearson''s r');
        
        
        keyboard
        varargout = {D};
    case 'sim:Data'
        % case to simulate multi-finger data from a covariance matrix
        signal = varargin{1};
        noise  = varargin{2};
        G      = varargin{3};
        numSim = varargin{4};
        numVox = 500;
        numRun = 10;
        numCond= 31;
        
        % define flexible PCM model:
        M{1}.type       = 'component';
        M{1}.numGparams = 1;
        M{1}.Gc         = G; 
        
        % simulate data:
        Z = kron(ones(numRun,1),eye(numCond));
        Y = pcm_makeDataset(M{1},1,'design',Z,'signal',signal,'noise',noise,...
            'numVox',numVox,'numSim',numSim);
        
        % make indicator vectors
        for ii=1:numSim
            cV{ii} = kron(ones(numRun,1),[1:numCond]');
            pV{ii} = kron([1:numRun]',ones(numCond,1));
        end

        varargout = {Y,pV,cV};

    case 'est_FlexibleParams' 
        % Ridge regression estimate of patterns
        % - use pcm with fixed model G to estimate signal and noise
        % parameters
        % - employ ridge regression with lambda = noise param / signal
        % param

        % inputs
        Y   = varargin{1}; % matrix of activity patterns [#conds*#runs x #vox]
        pV  = varargin{2}; % partition vector (assume cV and pV are same across subjs)
        cV  = varargin{3}; % condition vector (chord #s)
        G   = varargin{4}; % model prior
        
        % create single finger feature matrix
        Z0  = pp1_encoding('chords');
        Gm  = pinv(Z0)*G*pinv(Z0)';
        Zsf = kron(ones(numel(unique(pV)),1),Z0);
        
        % pcm non-linear model structure
        M{1}.type = 'nonlinear';
        M{1}.numGparams = 4;
        M{1}.theta0     = log([0.9 0.8 0.7 0.6])';
        M{1}.Ac         = pcm_diagonalize(Gm); 
        M{1}.modelpred  = @pcmSubj_modelpred_flexible;
               
        % fit model G to get noise and signal params:
        [~,theta_hat,~,tt] = pcm_fitModelIndivid({Y},M,pV,cV,'runEffect','none','verbose',0,'fitScale',1);
        
        % get scale parameters:
        scaleParam = theta_hat{1}(1:4)'; % scaling relative to single fingers
        
        varargout = {scaleParam};    
    
end
end % pp1_ana

function rss = modelLossRSS(theta,U,Utrain,modelName)
Ypred = pp1_ana('model:predictPatterns',U,modelName,theta); % predict patterns under perscribed model
Ypred = Ypred-mean(Ypred,1); % rmv voxel means
Utrain = Utrain-mean(Utrain,1);
rss   = sum(sum((Utrain-Ypred).^2)); % calculate L2 loss (RSS)
end

function [G,dGdtheta] = pcmSubj_modelpred_flexible(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model and returns the derivative in respect to the parameters
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
scaleParams = exp(theta(1:4));
% single finger features
A = Model.Ac(:,:,1); % A = pcm_diagonalize(G_singleFinger from group)
OM = A*A';
% activity scaling feature
chords     = pp1_encoding('chords');
M          = chords;
numFingers = sum(M,2);
for i = 1:4
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end         
G  = M*OM*M';  % Second moment matrix

for i = 1:4 % scale params
    dM                    = zeros(size(chords));
    dM(numFingers==i+1,:) = chords(numFingers==i+1,:);
    dGdtheta(:,:,i)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)       = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end
end

