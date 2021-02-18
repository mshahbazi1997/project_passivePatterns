function varargout = pp1_selectivity(what,varargin)
% case to test different models of single finger selectivity
switch what
    case 'chords'
        % returns matrix of finger chords
       varargout = {[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;...             % singles            5
             1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;...                        % doubles (thumb)    4
             0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;...                                   % doubles            3
             0 0 1 1 0; 0 0 1 0 1;...                                              % doubles            2
             0 0 0 1 1;...                                                         % doubles            1
             1 1 1 0 0; 1 1 0 1 0; 1 1 0 0 1; 1 0 1 1 0; 1 0 1 0 1; 1 0 0 1 1;...  % triples (thumb)    6
             0 1 1 1 0; 0 1 1 0 1; 0 1 0 1 1; 0 0 1 1 1;...                        % triples            4
             1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 0 1 1 1; 0 1 1 1 1;                % quadruples         5
             1 1 1 1 1]};                                                           % all five           1
    case 'usageG'
        % loads single finger nat stats G
        load(recipeFile);
        varargout = {Model(2).G_cent};
    case '0' % cases to calculate statistics on patterns    
    case 'SFT:calcFstat'
        % calculates F-statistic per voxel to determine if voxel is
        % significantly modulated by finger(s)
        Y        = varargin{1}; % N x P matrix of data. (N=numCond*numRun, P=numVox)
        condVec  = varargin{2}; % N x 1 vector of condition assignments
        runVec   = varargin{3}; % N x 1 vector of run assignments
        % housekeeping
        numVox   = size(Y,2);
        conds    = unique(condVec)';
        numCond  = numel(conds);
        runs     = unique(runVec)';
        numRun   = numel(runs);
        df1 = numCond-1;
        df2 = numCond*numRun - numCond - numRun;
        % remove run means 
        C0  = indicatorMatrix('identity',runVec);
        Y   = Y - C0*pinv(C0)*Y; 
        
        % compute mean and covariance matrices
        muK = zeros(numCond,numVox); % condition means
        SSR = zeros(1,numVox);       % ssr vector (common across conditions)
        n   = zeros(1,numCond);      % # observations per condition
        for ii=1:numCond
            c = conds(ii);
            idx  = find(condVec==c);
            n(ii)= numel(idx);
            muK(ii,:) = sum(Y(idx,:),1) ./ n(ii); % condition means
            res = bsxfun(@minus,Y(idx,:),muK(ii,:));
            SSR = SSR + sum(res.^2,1) ./ n(ii); % SSR (across conditions) scaled by number of observations
        end
        SSB = sum(muK.^2,1); 
        F   = (SSB./df1) ./ (SSR./df2);
        Fcrit = finv(0.95,df1,df2); % 95% cutoff for F-stat
        % calculate significance of F-stats:
        p=fcdf(F,df1,df2);
        eps=0.000001; 
        p(p>1-eps)=1-eps;
        p(p<eps)=eps;    
        varargout = {F,Fcrit,p};
    case 'SFT:calcFstatCV'
        
        % * * * REVIEW THIS! NOT SURE IF CORRECT! * * *
        
        error('please review ''calcFstatCV'' before use');
        % calculates F-statistic per voxel to determine if voxel is
        % significantly modulated by finger(s)
        Y        = varargin{1}; % N x P matrix of data. (N=numCond*numRun, P=numVox)
        condVec  = varargin{2}; % N x 1 vector of condition assignments
        runVec   = varargin{3}; % N x 1 vector of run assignments
        % housekeeping
        numVox   = size(Y,2);
        conds    = unique(condVec)';
        numCond  = numel(conds);
        runs     = unique(runVec)';
        numRun   = numel(runs);
        df1 = numCond-1;
        df2 = numCond*numRun - numCond - numRun +1; % account for leave-one-out cv
        % remove run means 
        C0  = indicatorMatrix('identity',runVec);
        Y   = Y - C0*pinv(C0)*Y; 
        
        % calculate CrossValidated F
        A = zeros(numCond,numVox,numRun); % zero-pad
        ApredCV = A;
        for ii = 1:numRun
            A(:,:,ii) = Y(runVec==ii,:);
        end
        muK = mean(A,3);
        res = bsxfun(@minus,A,muK);
        SSR = sum(sum(res.^2,1),3); % common residual covariance
        SSB = sum((muK.^2).*numRun,1);
        F   = (SSB./df1) ./ (SSR./df2);
        for ii = runs
            testRuns = runs~=ii;
            ApredCV(:,:,ii) = mean(A(:,:,testRuns),3);
        end
        % crossval f-stat
        %Apred = repmat(mean(A,3),1,1,numRun); % predicted voxel tunings
        TSS   = sum(sum((A).^2,3),1);         % total SS 
        %RSS   = sum(sum((A-Apred).^2,3),1);   % finger model (5 params)
        RSScv = sum(sum((A-ApredCV).^2,3),1); % unrestricted SSR
        SSBcv = TSS-RSScv;
        FstatCV = (SSBcv./df1) ./ (RSScv./df2);
        FcritCV = finv(0.95,df1,df2); % 95% cutoff for F-stat
        varargout = {FstatCV,FcritCV};
        %keyboard
    case 'SFT:estimateSFT'
        % calculate single finger tuning using normalzied distance approach
        X         = varargin{1}; % CxN matrix of data.
        numC      = size(X,1);
        maxX      = max(X,[],1);
        avgDistsN = (sum(maxX-X)./(numC-1)) ./ (maxX-min(X,[],1));
        varargout = {avgDistsN};
    case 'SFT:estimateVariances'
        % empirically estimates error variance in activity patterns across
        % runs. 
        % Estimate is more accurate when run means have been removed.
        % r = varSignal / (varSignal + varError)
        % varError = (varSignal / r) - varSignal
        % where 'r' is the signal correlation avg. across runs (correlation
        % across all voxels and conditions per run)
        
        % % NOTE: we integrate across conditions
        
        Y = varargin{1};    % patterns   [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        r = varargin{3};    % partitions [regressors x 1]
        
        nV = size(Y,2);         % # voxels
        nP = numel(unique(r));  % # partitions
        nC = numel(unique(c));  % # conditions
        Ya = zeros(nP,nV*nC);   % pre-allocate
        
        for pp = 1:nP
            y = Y(r==pp,:);     % vectorize patterns across conditions per run
            Ya(pp,:) = y(:)';
        end
        
        take = logical(tril(ones(nP),-1)); % lower-triangular index
        G    = cov(Ya');  % covariances between runs (each row = one run, has zero mean)
        R    = corr(Ya'); % correlations between runs
        
        var_S = sum(G(take))/sum(sum(take)); % signal variance (avg. across runs)
        r     = sum(R(take))/sum(sum(take)); % signal correlation (avg. across runs)
        var_E = var_S/r - var_S;             % error variance (avg. across runs)
        
        varargout = {var_E,var_S};
    case 'estErrVarMeans'
        % empirically estimates error variance in activity patterns across
        % runs. 
        % Estimate is more accurate without run mean removal.
        
        % When you remove run means, the error variance estimate is biased
        % to a lower estimate. 
        % This has important implications for whether you should remove run means or not.
        % If you want to show that some model is lower than the noise
        % ceiling, then you should submit betas WITHOUT run means removed.
        %       This is because run mean removal will result in lower error
        %       estimate, and so your model will have MORE signal (i.e.
        %       will be harder to fit).
        % And the opposite applies. If you want to show no difference,
        % estimate error variance using betas WITH run means removed.
        
        % % NOTE: patterns should NOT have run means removed, as this is needed for empirical noise estimates 
        
        Y = varargin{1};    % patterns   [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        r = varargin{3};    % partitions [regressors x 1]
        
        nV = size(Y,2);
        nP = numel(unique(r));
        nC = numel(unique(c));
        Yc = zeros(nP,nV*nC); % vectorized mean-removed patterns
        Ym = Yc;              % vectorized run-mean patterns
        Yo = Yc;              % vectorized original patterns
        
        C0 = indicatorMatrix('identity',r);
        Yr = C0*pinv(C0)*Y; % calc run-means
        for pp = 1:nP
            % vectorize condition patterns
            y = Y(r==pp,:)-Yr(r==pp,:); % remove run-means
            Yc(pp,:) = y(:)';
            % vectorize run-mean patterns
            y = Yr(r==pp,:);
            Ym(pp,:) = y(:)';
            % vectorize original patterns
            y = Y(r==pp,:);
            Yo(pp,:) = y(:)';
        end
        
        take = logical(tril(ones(nP),-1));
        Gc   = cov(Yc');
        Gm   = cov(Ym');
        R    = corr(Yo');
        
        denom = sum(sum(take));
        var_M = sum(Gm(take)) / denom;
        var_S = sum(Gc(take)) / denom;
        r     = sum(R(take)) / denom;
        
        var_E = (var_M+var_S)/r - var_M - var_S; 
        varargout = {var_E,var_S};
    case 'SFT:expectedValue_G'
        % calculates expected value (avg. sft across voxels) for voxels generated under specfified G
        
        % simulation params
        errvar  = varargin{1};
        sigvar  = varargin{2};
        G       = varargin{3};
        numVox  = varargin{4};
        numRun  = varargin{5};
        numSim  = varargin{6};

        % generate data
        D = pp1_selectivity('SFT:model_Gnoise',G,numVox,numRun,numSim,sigvar,errvar); % noise patterns
        % calc expected tuning on simulated datasets
        v = nan(1,numSim);
        sft = v;
        evar_est = v;
        svar_est = v; % preallocate
        for s = 1:numSim
            d      = getrow(D,D.sn==s);
            % remove run means and estimate simulated err and sig vars
            C0     = indicatorMatrix('identity',d.run); 
            d.Y    = d.Y -C0*pinv(C0)*d.Y; % remove run means
            [evar_est(s),svar_est(s)] = pp1_selectivity('SFT:estimateVariances',d.Y,d.cond,d.run);
            d = tapply(d,{'cond'},{'Y','mean'}); % avg. betas across runs for this simulated subject
            tmp_sft = pp1_selectivity('SFT:estimateSFT',d.Y);
            sft(s) = mean(tmp_sft); % mean tuning across voxels for this simulated dataset
        end
        varargout = {sft,mean(evar_est),mean(svar_est)};
    case 'SFT:expectedValue_Sparse'
        % calculates expected value (avg. sft across voxels) for voxels
        % with sparse tuning but corrupted by some degree of noise
        
        % simulation params
        errvar  = varargin{1};
        sigvar  = varargin{2};
        sparsity= varargin{3}; % 1-totally sparse, 2-two finger tuned, etc.
        numVox  = varargin{4};
        numRun  = varargin{5};
        numSim  = varargin{6};

        % generate data
        D = pp1_selectivity('SFT:model_Sparse',sparsity,numVox,numRun,numSim,sigvar,errvar); % sparse patterns with noise
        % calc expected tuning on simulated datasets
        v = nan(1,numSim);
        sft = v;
        evar_est = v;
        svar_est = v; % preallocate
        for s = 1:numSim
            d      = getrow(D,D.sn==s);
            % remove run means and estimate simulated err and sig vars
            C0     = indicatorMatrix('identity',d.run); 
            d.Y    = d.Y -C0*pinv(C0)*d.Y; % remove run means
            [evar_est(s),svar_est(s)] = pp1_selectivity('SFT:estimateVariances',d.Y,d.cond,d.run);
            d = tapply(d,{'cond'},{'Y','mean'}); % avg. betas across runs for this simulated subject
            tmp_sft = pp1_selectivity('SFT:estimateSFT',d.Y);
            sft(s) = mean(tmp_sft); % mean tuning across voxels for this simulated dataset
        end
        varargout = {sft,mean(evar_est),mean(svar_est)};
    case '0' % cases to simulate data under different generative distributions:
    case 'SFT:model_G'
        % true betas drawn from normal distribution
        muS     = 0; % signal mean (zero)- arbitrary
        G       = varargin{1};
        numVox  = varargin{2};
        numSim  = varargin{3};
        numCond = size(G,1);
        v = ones(numCond,1); % helper vector
        D = []; % output structure
        for s = 1:numSim
            U = mvnrnd_exact(G,numVox);
%             try
%                 U = mvnrnd(ones(numCond,1).*muS,G,numVox)'; % generate true patterns
%             catch % make G semi-positive definite
%                 A  = pcm_diagonalize(G); 
%                 Gd = A*A';
%                 U  = mvnrnd(ones(numCond,1).*muS,Gd,numVox)'; % generate true patterns
%             end
            d.Y    = U;
            d.sn   = v.*s;
            D = addstruct(D,d);
        end
        varargout = {D};
    case 'SFT:model_Gnoise'
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
    case 'SFT:model_Sparse'
        % true patterns are sparse (0's and 1's), with normally
        % distributed noise across runs
        sparsity = varargin{1};% sparsity level (1=totally sparse, 2-two conditions tuned, etc.)
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5}; % signal variance
        noise   = varargin{6}; % noise variance
        
        % get # conditions based on sparsity level
        chords = pp1_selectivity('chords');
        numDigits = sum(chords,2);
        chords = logical(chords(numDigits==sparsity,:));
        [numChord,numCond] = size(chords);
        numVoxPerChord = ceil(numVox/numChord); % voxels tuned per chord
        numVox = numVoxPerChord*numChord;
        signal = signal*(numCond/sparsity); % rescale the signal to account for the lack of signal variance aross non-responsive voxels
        % NOTE: 
        % The estimated signal will decay as sparsity decays.
        % This decay will not align with the decay that occurs for mvnrnd
        % patterns. 
        % This is because, as we decrease the sparsity level, we increase 
        % proportion of the overall signal that is accounted for as the run
        % means.
        
        % define signal generation
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        D = []; % output structure
        for s = 1:numSubj % per simulated dataset
            % draw uniform values for signal and noise (to be inverted
            % through any arbitrary function later)
            pNoise = unifrnd(0,1,numCond*numRun,numVox); 
            % Generate true sparse patterns
            U = kron(chords',ones(1,numVoxPerChord));
            % scale activity to match specified signal strength
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
    case 'SFT:model_Uniform'
        % true betas drawn from normal distribution, with added i.i.d. noise
        G       = varargin{1};
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5};
        noise   = varargin{6};
        
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        signalDist = @(x) (x*1);%unifinv(x,0,1);  % uniform inverse for signal generation
        
        numCond = size(G,1);
        D = []; % output structure
        for s = 1:numSubj
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
    case 'plotSimulatedHist'
        % make histogram plots of simulated data
        numCond = 5;
        svar = 1;
        evar = 0.1;
        G=pp1_selectivity('usageG');
        N = pp1_selectivity('SFT:model_Gnoise',eye(numCond),100*numCond,1,1,svar,evar);
        S = pp1_selectivity('SFT:model_Sparse',1,100*numCond,1,1,svar,evar);
        U = pp1_selectivity('SFT:model_Uniform',eye(numCond),100*numCond,1,1,svar,evar);
        
        % plot
        sty1 = style.custom({'blue'});
        sty2 = style.custom({'red'});
        sty3 = style.custom({'green'});
        subplot(2,3,1); plt.hist(N.Y(:),'style',sty1); title('mvnrnd'); xlabel('activity'); ylabel('count');
        subplot(2,3,2); plt.hist(S.Y(:),'style',sty2); title('sparse');
        subplot(2,3,3); plt.hist(U.Y(:),'style',sty3); title('uniform');
        
        subplot(2,3,4); plt.scatter(N.Y(1,:)',N.Y(2,:)','style',sty1,'regression','off'); title('mvnrnd'); xlabel('cond 1 act.'); ylabel('cond 2 act.');
        subplot(2,3,5); plt.scatter(S.Y(1,:)',S.Y(2,:)','style',sty2,'regression','off'); title('sparse');
        subplot(2,3,6); plt.scatter(U.Y(1,:)',U.Y(2,:)','style',sty3,'regression','off'); title('uniform');
    case '0' % cases to run test simulations on selectivity analysis:
    case 'testErrVarV1'
        % case to generate data at different levels of signal-to-noise,
        % plots results for evaluation
        numSim = 10; % 100 simulated datasets per model type
        evar   = [0.1:0.1:1];
        numCond = 5;
        G      = pp1_selectivity('usageG');%eye(numCond);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        svar = 1; % signal variance
        for ii = evar
            d_sparse = pp1_selectivity('SFT:model_Sparse',1,numVox,numRun,numSim,svar,ii);
            d_norm   = pp1_selectivity('SFT:model_Gnoise',G,numVox,numRun,numSim,svar,ii);
            %d_uni    = pp1_selectivity('SFT:model_Uniform',G,numVox,numRun,numSim,svar,ii);
            d = [];
            for jj = 1:numSim
                dn = getrow(d_norm,d_norm.sn==jj);  % normal
                ds = getrow(d_sparse,d_sparse.sn==jj); % sparse
                %du = getrow(d_uni,d_uni.sn==jj);    % uniform
                % indexing fields to output structure
                q.sn   = jj;
                q.evar = ii;
                q.svar = svar;

                % calc snr for normal data
                q.model       = 1;
                C0 = indicatorMatrix('identity',dn.run);
                dn.Y = dn.Y - C0*pinv(C0)*dn.Y;
                [var_e,var_s] = pp1_selectivity('SFT:estimateVariances',dn.Y,dn.cond,dn.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dnn           = tapply(dn,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_selectivity('SFT:estimateSFT',dnn.Y));
                d = addstruct(d,q);

                % calc snr for sparse data
                q.model       = 2;
                C0 = indicatorMatrix('identity',ds.run);
                ds.Y = ds.Y - C0*pinv(C0)*ds.Y;
                [var_e,var_s] = pp1_selectivity('SFT:estimateVariances',ds.Y,ds.cond,ds.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dss           = tapply(ds,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_selectivity('SFT:estimateSFT',dss.Y));
                d = addstruct(d,q);

            end
            D = addstruct(D,d);
        end
                
        sty = style.custom({'darkgray','red','orange'});
        
        % calculates snrs:
        D.snr_true = D.svar ./ D.evar;
        D.snr_est  = D.svar_Est ./ D.evar_Est;
        
        % apply log scale (slightly
        % tricky for log(0) so that's why the code below exists)
        D.evarLog = D.evar;
        D.evarLog(D.evar>0) = log(D.evar(D.evar>0));
        D.evarLog(D.evar==0)= -sum(diff(log(evar(evar>0))));

        % plot SIGNAL VARIANCE ESTIMATES (unchanged)
        subplot(1,4,1);
        plt.line(D.evar,D.svar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(signal)');
        xlabel('true var(error)');
        title('unchanged signal (run means rmv)');
        axis square
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        
        % plot ERROR VARIANCE ESTIAMTES (varies)
        subplot(1,4,2);
        plt.line(D.evar,D.evar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error)');
        xlabel('true var(error)');
        title('altered error (run means rmv)');
        axis square
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        plt.match('y');
        
        % plot SNRs:
        subplot(1,4,3);
        plt.scatter(D.snr_true,D.snr_est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('snr estimate');
        xlabel('snr true');
        title('true vs. estiated snr');
        axis square
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        %set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        
        % plot SFTs
        subplot(1,4,4);
        plt.line(log(D.snr_true),D.sftEst,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('selectivity index');
        xlabel('true snr (log scale)');
        title('selectivity (run means rmv)');
        axis square
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
                
        varargout = {D};
    case 'testSigVarV1'
        % case to generate data at different levels of signal-to-noise,
        % plots results for evaluation
        numSim = 100; % 100 simulated datasets per model type
        svar   = [0.1:0.1:1];
        numCond = 5;
        G      = pp1_selectivity('usageG');%eye(numCond);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        evar = 1; % error variance
        for ii = svar
            d_sparse = pp1_selectivity('SFT:model_Sparse',1,numVox,numRun,numSim,ii,evar);
            d_norm   = pp1_selectivity('SFT:model_Gnoise',G,numVox,numRun,numSim,ii,evar);
            %d_uni    = pp1_selectivity('SFT:model_Uniform',G,numVox,numRun,numSim,ii,evar);
            d = [];
            for jj = 1:numSim
                dn = getrow(d_norm,d_norm.sn==jj);  % normal
                ds = getrow(d_sparse,d_sparse.sn==jj); % sparse
                %du = getrow(d_uni,d_uni.sn==jj);    % uniform
                % indexing fields to output structure
                q.sn   = jj;
                q.evar = evar;
                q.svar = ii;

                % calc snr for normal data
                q.model       = 1;
                C0 = indicatorMatrix('identity',dn.run);
                dn.Y = dn.Y - C0*pinv(C0)*dn.Y;
                [var_e,var_s] = pp1_selectivity('SFT:estimateVariances',dn.Y,dn.cond,dn.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dnn           = tapply(dn,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_selectivity('SFT:estimateSFT',dnn.Y));
                d = addstruct(d,q);

                % calc snr for sparse data
                q.model       = 2;
                C0 = indicatorMatrix('identity',ds.run);
                ds.Y = ds.Y - C0*pinv(C0)*ds.Y;
                [var_e,var_s] = pp1_selectivity('SFT:estimateVariances',ds.Y,ds.cond,ds.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dss           = tapply(ds,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_selectivity('SFT:estimateSFT',dss.Y));
                d = addstruct(d,q);

            end
            D = addstruct(D,d);
        end
                
        sty = style.custom({'darkgray','red','orange'});
        % calculates snrs:
        D.snr_true = D.svar ./ D.evar;
        D.snr_est  = D.svar_Est ./ D.evar_Est;
        
        % apply log scale (slightly
        % tricky for log(0) so that's why the code below exists)
        D.svarLog = D.svar;
        D.svarLog(D.svar>0) = log(D.svar(D.svar>0));
        D.svarLog(D.svar==0)= -sum(diff(log(svar(svar>0))));
        D.svar_Est(D.svar_Est<0) = realmin;
        
        % plot SIGNAL VARIANCE ESTIMATES (varies)
        subplot(1,4,1);
        plt.line(D.svar,D.svar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(signal)');
        xlabel('true var(signal)');
        title('altered signal (run means rmv)');
        axis square
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);

        % plot ERROR VARIANCE ESTIMATES (unchanged)
        subplot(1,4,2);
        plt.line(D.svar,D.evar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error)');
        xlabel('true var(signal)');
        title('true vs. estimated snr');
        axis square
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);
        plt.match('y');
        
        % plot SNRs:
        subplot(1,4,3);
        plt.scatter(D.snr_true,D.snr_est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('snr estimate');
        xlabel('snr true');
        title('true vs. estimated snr');
        axis square
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        %set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        
        % plot SFTs
        subplot(1,4,4);
        plt.line(log(D.snr_true),D.sftEst,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('selectivity index');
        xlabel('true snr (log scale)');
        title('selectivity (run means rmv)');
        axis square
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);
                
        varargout = {D};
    case 'testErrVarV2_depreciated'
        % case to generate data at different levels of signal-to-noise,
        % plots results for evaluation
        numSim = 100; % 100 simulated datasets per model type
        evar   = [0,0.1,1,10:10:100];
        numCond = 5;
        G      = eye(numCond);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        svar = 1; % signal variance
        for ii = evar
            d_sparse = pp1_selectivity('SFT:model_Sparse_single',numCond,numVox,numRun,numSim,svar,ii);
            d_norm   = pp1_selectivity('SFT:model_Gnoise',G,numVox,numRun,numSim,svar,ii);
            d_uni    = pp1_selectivity('SFT:model_Uniform',G,numVox,numRun,numSim,svar,ii);
            d = [];
            for jj = 1:numSim
                dn = getrow(d_norm,d_norm.sn==jj);  % normal
                ds = getrow(d_sparse,d_sparse.sn==jj); % sparse
                du = getrow(d_uni,d_uni.sn==jj);    % uniform
                % indexing fields to output structure
                q.sn   = jj;
                q.evar = ii;
                q.svar = svar;

                % calc snr for normal data
                q.model       = 1;
                [var_e,var_s] = pp1_selectivity('estErrVarMeans',dn.Y,dn.cond,dn.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dnn           = tapply(dn,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_selectivity('estSingleFingerTuning',dnn.Y));
                d = addstruct(d,q);

                % calc snr for sparse data
                q.model       = 2;
                [var_e,var_s] = pp1_selectivity('estErrVarMeans',ds.Y,ds.cond,ds.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dss           = tapply(ds,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_selectivity('estSingleFingerTuning',dss.Y));
                d = addstruct(d,q);

%                 % calc snr for uniform data
%                 q.model       = 3;
%                 [var_e,var_s] = pp1_selectivity('estErrVarMeans',du.Y,du.cond,du.run);
%                 q.evar_Est    = var_e;
%                 q.svar_Est    = var_s;
%                 duu           = tapply(du,{'cond'},{'Y','mean'});
%                 q.sftEst      = mean(pp1_selectivity('estSingleFingerTuning',duu.Y));
                
                d = addstruct(d,q);
            end
            D = addstruct(D,d);
        end
                
        sty = style.custom({'darkgray','red','orange'});
        
        % plot estimate of the error (varies)
        subplot(1,3,2);
        plt.line(D.evar,D.evar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error)');
        xlabel('specified var(error)');
        title('altered error (v.2)');
        axis square
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);

        % plot estimate of the signal (unchanged)
        subplot(1,3,1);
        plt.line(D.evar,D.svar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(signal)');
        xlabel('specified var(error)');
        title('unchanged signal (v.2)');
        axis square
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        
        % plot estimated selectivities
        subplot(1,3,3);
        plt.line(D.evar,D.sftEst,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('sft');
        xlabel('specified var(error)');
        title('selectivity (v.2)');
        axis square
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
                
        varargout = {D};
    case 'testSigVarV2_depreciated'
        % case to generate data at different levels of signal-to-noise,
        % plots results for evaluation
        numSim = 20; % 100 simulated datasets per model type
        svar   = [0.1:0.1:1,2:10];
        sparsity= 1; % sparsity level (1=perfect sparse, 2=two conditions tuned, ..etc.)
        numCond = 5;
        G      = eye(numCond);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        evar = 1; % error variance
        for ii = svar
            %d_sparse = pp1_selectivity('SFT:model_Sparse_single',numCond,numVox,numRun,numSim,ii,evar); % this version does not allow chaning of sparsity
            d_sparse = pp1_selectivity('SFT:model_Sparse',sparsity,numVox,numRun,numSim,ii,evar); % allows changing sparsity levels
            d_norm   = pp1_selectivity('SFT:model_Gnoise',G,numVox,numRun,numSim,ii,evar);
            %d_uni    = pp1_selectivity('SFT:model_Uniform',G,numVox,numRun,numSim,ii,evar);
            d = [];
            for jj = 1:numSim
                dn = getrow(d_norm,d_norm.sn==jj);  % normal
                ds = getrow(d_sparse,d_sparse.sn==jj); % sparse
                %du = getrow(d_uni,d_uni.sn==jj);    % uniform
                % indexing fields to output structure
                q.sn   = jj;
                q.evar = evar;
                q.svar = ii;

                % calc snr for normal data
                q.model       = 1;
                [var_e,var_s] = pp1_selectivity('estErrVarMeans',dn.Y,dn.cond,dn.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dnn           = tapply(dn,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_selectivity('estSingleFingerTuning',dnn.Y));
                d = addstruct(d,q);

                % calc snr for sparse data
                q.model       = 2;
                [var_e,var_s] = pp1_selectivity('estErrVarMeans',ds.Y,ds.cond,ds.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dss           = tapply(ds,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_selectivity('estSingleFingerTuning',dss.Y));
                d = addstruct(d,q);

%                 % calc snr for uniform data
%                 q.model       = 3;
%                 [var_e,var_s] = pp1_selectivity('estErrVarMeans',du.Y,du.cond,du.run);
%                 q.evar_Est    = var_e;
%                 q.svar_Est    = var_s;
%                 duu           = tapply(du,{'cond'},{'Y','mean'});
%                 q.sftEst      = mean(pp1_selectivity('estSingleFingerTuning',duu.Y));
                
                d = addstruct(d,q);
            end
            D = addstruct(D,d);
        end
                
        sty = style.custom({'darkgray','red','orange'});
        
        % plot estimate of the signal (varies)
        subplot(1,3,1);
        plt.line(D.svar,D.svar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(signal)');
        xlabel('specified var(signal)');
        title('altered signal (v.2)');
        axis square
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);

        % plot estimate of the error (unchanged)
        subplot(1,3,2);
        plt.line(D.svar,D.evar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error)');
        xlabel('specified var(signal)');
        title('unchanged error (v.2)');
        axis square
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);
        
        % plot estimated selectivities
        subplot(1,3,3);
        plt.line(D.svar,D.sftEst,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('sft');
        xlabel('specified var(signal)');
        title('selectivity (v.2)');
        axis square
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);
                
        varargout = {D};
    case 'testFthres'
        % tests that restricting sft analysis to voxels with significant
        % responses does not introduce any bias
        numSim = [];
        svar = [];
        evar = [];
        vararginoptions(varargin,{'svar','evar','numSim'});
        numCond = 5;
        G       = eye(numCond);
        % use avg. single finger G from 3b
%         G = [ 0.11931    -0.014803    -0.035168    -0.037241     -0.02326
%     -0.014803      0.10113    -0.024087    -0.038672    -0.035275
%     -0.035168    -0.024087      0.08886    -0.015966    -0.024505
%     -0.037241    -0.038672    -0.015966     0.090358   -0.0042599
%      -0.02326    -0.035275    -0.024505   -0.0042599      0.10681];
        %G = diag([1,2,3,4,5]);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        models = [1,2]; % mvnrnd and sparse
        
        for mm=models % for each model
            for ee=evar % for each noise lvl
                for ii=1:numel(svar) % for each signal lvl
                    ss = svar(ii);
                    % simulate data
                    switch mm
                        case 1 % mvnrnd
                            simData = pp1_selectivity('SFT:model_Gnoise',G,numVox,numRun,numSim,ss,ee);
                        case 2
                            simData = pp1_selectivity('SFT:model_Sparse',1,numVox,numRun,numSim,ss,ee);
                    end
                    for jj = 1:numSim
                        % indexing fields to output structure
                        q.model= mm;
                        q.sn   = jj;
                        q.evar = ee;
                        q.svar = ss;
                        q.svarNum = ii;
                        % get data for this simulation
                        d = getrow(simData,simData.sn==jj);
                        % get f stat per voxel
                        [F,Fcrit,p] = pp1_selectivity('SFT:calcFstat',d.Y,d.cond,d.run);
                        % 1. do analysis without restricting to sig.
                        % voxels:
                        C0 = indicatorMatrix('identity',d.run);
                        d.Yall = d.Y - C0*pinv(C0)*d.Y; % remove run means
                        [q.evar_est,q.svar_est] = pp1_selectivity('SFT:estimateVariances',d.Yall,d.cond,d.run);
                        dall     = tapply(d,{'cond'},{'Yall','mean'}); % avg. across simulated runs
                        q.sftEst = mean(pp1_selectivity('SFT:estimateSFT',dall.Yall));
                        q.numVox = numVox;
                        q.F      = F;
                        q.Fcrit  = Fcrit;
                        q.fthres = 0;
                        D = addstruct(D,q);
                        % 2. do analysis, including only sig. voxels:
                        d.Yf = d.Y(:,F>=Fcrit); % using non-cv f stat as thresholder
                        d.Yf = d.Yf - C0*pinv(C0)*d.Yf; % remove run means
                        [q.evar_est,q.svar_est] = pp1_selectivity('SFT:estimateVariances',d.Yf,d.cond,d.run);
                        df     = tapply(d,{'cond'},{'Yf','mean'}); % avg. across simulated runs
                        q.sftEst = mean(pp1_selectivity('SFT:estimateSFT',df.Yf));
                        q.numVox = size(df.Yf,2);
                        q.F      = F;
                        q.Fcrit  = Fcrit;
                        q.fthres = 1;
                        D = addstruct(D,q);
                    end
                end
            end
        end
        
        % apply log scale (slightly
        % tricky for log(0) so that's why the code below exists)
        D.svarLog = D.svar;
        D.svarLog(D.svar>0) = log(D.svar(D.svar>0));
        D.svarLog(D.svar==0)= -sum(diff(log(svar(svar>0))));
        
        varargout = {D};
    case 'plotFtest'
        % Case to simulate data under different signal conditions and apply
        % F-test thresholding prior to calculating sparsity of tuning.
        % This case is to check that F-thresholding does not bias the
        % estimated sparsity.
        % At each signal level, data should show similar % voxels with
        % significant activity. (subplot 1)
        % But, with increasing signal, only the sparsely tuned data should
        % show increasing sparsity values. (subplot 2)
        numSim = 10;
        %svar = [0,0.001,0.01,0.1,0.5,1,10];
        svar = [0:0.1:1];
        evar = 1;
        D = pp1_selectivity('testFthres','svar',svar,'evar',evar,'numSim',numSim);
        D.se_true = D.svar./D.evar; % specified SE ratio
        D.se_hat  = D.svar_est./D.evar_est; % estimated SE ratio
        D = getrow(D,D.evar==1);

        % plot styling
        labels = string(svar);
        sty = style.custom({'darkgray','red','darkgray','red'});
        sty.general.linestyle  = {'-','-',':',':'};
        sty.general.markertype = {'o','o','^','^'};
        
        % plot se ratios:
        subplot(1,3,1);
        plt.line(D.se_true,D.se_hat,'split',[D.fthres D.model],'style',sty);
        set(gca,'xtick',[0:0.1:1],'xticklabel',[0:0.1:1],'ytick',[0:0.1:1],'xticklabelrotation',45);
        refline(1,0);
        xlabel('specified SE var (ratio)');
        ylabel('estimated SE var (ratio)');
        title('simulated signal strength');
        legend off
        
        % plot estimated sfts:
        subplot(1,3,2)
        plt.line(D.svarNum,D.sftEst,'split',[D.fthres D.model],'style',sty,'errorfcn','stderr');
        ylabel('selectivity index');
        xlabel('signal variance');
        title('estimated selectivity')
        set(gca,'xticklabel',labels,'xticklabelrotation',45);
        %legend off
        
        % plot % voxels retained after f-thresholding:
        warning off
        subplot(1,3,3);
        plt.bar(D.svarNum,(D.numVox/500).*100,'split',D.model,'style',sty,'subset',D.fthres==1);
        xlabel('signal variance');
        ylabel('% voxels significant');
        title('% significant voxels (omnibus F-test)');
        legend off
        setXTicksGroups(gca,numel(svar),2,labels);
        set(gca,'xticklabelrotation',45);
        warning on
        ylim([0 100]);
        
        varargout = {D};
    case '0' % do diff # of voxels influence selectivity (assuming signal var does not scale with more voxels)?
    case 'testNumVoxels'
        % tests that restricting sft analysis to voxels with significant
        % responses does not introduce any bias
        numSim = [];
        svar = [];
        evar = [];
        numVox = [];
        vararginoptions(varargin,{'svar','evar','numSim','numVox'});

        % use avg. single finger G from 3b
%         G = [ 0.11931    -0.014803    -0.035168    -0.037241     -0.02326
%     -0.014803      0.10113    -0.024087    -0.038672    -0.035275
%     -0.035168    -0.024087      0.08886    -0.015966    -0.024505
%     -0.037241    -0.038672    -0.015966     0.090358   -0.0042599
%      -0.02326    -0.035275    -0.024505   -0.0042599      0.10681];
        G = eye(5);

        numRun = 11;
        D      = [];
        models = [1,2]; % mvnrnd and sparse
        
        for mm=models % for each model
            for ee=evar % for each noise lvl
                for ii=1:numel(numVox) % for each signal lvl
                    pp = numVox(ii);
                    % simulate data
                    switch mm
                        case 1 % mvnrnd
                            simData = pp1_selectivity('SFT:model_Gnoise',G,numVox,numRun,numSim,svar,ee);
                        case 2
                            simData = pp1_selectivity('SFT:model_Sparse',1,numVox,numRun,numSim,svar,ee);
                    end
                    for jj = 1:numSim
                        % indexing fields to output structure
                        q.model= mm;
                        q.sn   = jj;
                        q.evar = ee;
                        q.svar = svar;
                        q.numVoxNum = ii;
                        q.numVox = pp;
                        % get data for this simulation
                        d = getrow(simData,simData.sn==jj);
                        % get f stat per voxel
                        [F,Fcrit] = pp1_selectivity('SFT:calcFstat',d.Y,d.cond,d.run);
                        % 1. do analysis without restricting to sig.
                        % voxels:
                        C0 = indicatorMatrix('identity',d.run);
                        d.Yall = d.Y - C0*pinv(C0)*d.Y; % remove run means
                        [q.evar_est,q.svar_est] = pp1_selectivity('SFT:estimateVariances',d.Yall,d.cond,d.run);
                        dall     = tapply(d,{'cond'},{'Yall','mean'}); % avg. across simulated runs
                        q.sftEst = mean(pp1_selectivity('SFT:estimateSFT',dall.Yall));
                        q.F      = F;
                        q.Fcrit  = Fcrit;
                        q.fthres = 0;
                        D = addstruct(D,q);
                        % 2. do analysis, including only sig. voxels:
                        d.Yf = d.Y(:,F>=Fcrit); % using non-cv f stat as thresholder
                        d.Yf = d.Yf - C0*pinv(C0)*d.Yf; % remove run means
                        [q.evar_est,q.svar_est] = pp1_selectivity('SFT:estimateVariances',d.Yf,d.cond,d.run);
                        df     = tapply(d,{'cond'},{'Yf','mean'}); % avg. across simulated runs
                        q.sftEst = mean(pp1_selectivity('SFT:estimateSFT',df.Yf));
                        q.numVox = size(df.Yf,2);
                        q.F      = F;
                        q.Fcrit  = Fcrit;
                        q.fthres = 1;
                        D = addstruct(D,q);
                    end
                end
            end
        end
        
        % apply log scale (slightly
        % tricky for log(0) so that's why the code below exists)
        D.svarLog = D.svar;
        D.svarLog(D.svar>0) = log(D.svar(D.svar>0));
        D.svarLog(D.svar==0)= -sum(diff(log(svar(svar>0))));
        
        varargout = {D};
    case 'plot_numVoxels'
        svar = 0.5;
        evar = 1;
        numVox = [100:100:1000];
        numSim = 10;
        D = pp1_selectivity('testFthres','svar',svar,'evar',evar,'numSim',numSim);
end
end

%% local functions
function xt = setXTicksGroups(ax,k1,k2,labels)
% assigns xtick to middle of X-group splits in a figure.
% ax = handle to current axis with graph
% k1 = number of main groups
% k2 = number of subgroups per group
xot = get(ax,'xtick');
xt = zeros(1,k1);
kidx = kron([1:k1]',ones(k2,1));
for k=1:k1
    xt(1,k)=mean(xot(kidx==k));
end
set(ax,'xtick',xt);
if ~isempty(labels)
    set(ax,'xticklabel',labels);
end
end
