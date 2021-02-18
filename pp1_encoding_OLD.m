function varargout=pp1_encoding(what,varargin)
%% details

% saarbuckle 2020

% encoding simulations for passivePatterns1 fmri experiment

% Directories for data and saving simulation results
%baseDir = '/Volumes/MotorControl/data/passivePatterns/passivePatterns1/fmri';
baseDir = '/Users/sarbuckle/DATA/passivePatterns1/fmri/'; % where we save model selection accuracies to
dataDir = fullfile(baseDir,'RegionOfInterest'); % where betas are saved (we use single finger betas from here for model simulations)
subj_name = {'pd01','pd02','s01','s02','s03','s04','s05','s06','s07','s08','s09'};
% roi details:
%   #   |   name
% ---------------
%   1   | ba 3a (hand area)
%   2   | ba 3b (hand area)
%   3   | ba 1 (hand area)
%   4   | ba 2 (hand area)
%   5   | rostral M1
%   6   | caudal M1
%   7   | SII
%   8   | ba 3a+3b+1+2 (al of s1 hand area)
%   9   | M1 hand area


% Models that we simulate data under:
% 'summation' - simply add constituent single finger patterns per chord
% 'summationNL' - take summated patterns, apply nonlinear scaling of activity with an exponent whose value >1

% We can also change how similar the models are- if the scaling param is
% close to 1, then the models are quite similar and thus it can be
% difficult to distinguish which dataset the true model comes from.

% CASES:
% 'chords' - returns chord indicator matrix
% 'getData' - returns activity patterns for subjects from a specified roi
% 'sim' - uses run-averaged single finger patterns from one subject to estimate multi-finger patterns + noise under one of two models.  
% 'doEncoding_noCV' - given simulated data, calculates model fits of summation and summation NL scaling model fits with simple crossvalidation. there is no crossval
% 'doEncoding_CV' - given simulated data, evaluates model fits using leave-one-run out crossval. 
% 'doPCM' - given simulated data, evaluate the fit of each model G (fixed model) in a crossvalidated individual fashion

% 'do' - wrapper case to simulate datasets, use three different model fitting approaches on the SAME simulated datasets, calculate model selection accuracies
% 'plot' - plot model selection accuracies

% 'compModelScaling' - compares the scaling of mean activity and mean variances of different models and the true data

% *factorize summation and summationNL
%% cases
switch what
    case 'do_DEPRECIATED'
        % wrapper to test different model fitting approaches and assess
        % model selection accuracies using different metrics
        
        % metrics:
        %   1 : r (enoding)
        %   2 : r2 (encoding)
        %   3 : r (pcm)
        %   4 : r2 (pcm)
        %   5 : likelihood (pcm)
        
        % sim params:
        snr = [0,0.01,0.1,1,10];
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        numSim = 50; % sims per subject per snr lvl under each model (2 models = 50 simulations * # snr lvls)
        C_inv = pinv(pcm_indicatorMatrix('identity',sum(pp1_encoding('chords'),2)));
        theta_baseline  = 0; % don't shift true baseline
        
        % define params for the linear and nonlinear models:
        modelNames = {'summation','summation_tanh'};
        modelThetas = {[theta_baseline], [theta_baseline]}; % {1}=linear, {2}=nonlinear
        
        % get data (load in once for all subjs to save time):
        [Y,partVec,condVec] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        
        % prep PCM model fitting struct:
        M = pp1_encoding('prep_pcmModels');
        
        % do simulations at different snrs for EACH subject:
        D = []; % overall output structure
        for ii=1:numel(sn)
            fprintf('s%02d snr: ',sn(ii));% verbose output to user
            % 1. simulate true Us under each model using this subject's
            % data:
            [U,G] = pp1_encoding('pcm_estTrueU',Y{ii},partVec{ii},condVec{ii});
            Usf = U(1:5,:);
            numVox = size(Usf,2);
            act_orig = mean(C_inv * U,2)';
            var_orig = (C_inv * diag(G))';
            
            % From now on, we use the simulated data!
            for jj=1:numel(snr)
                nn=snr(jj);
                fprintf('%2.2f...',nn);% verbose output to user
                % 2. simulate data under both models using Usf:            
                [U_hat_linear,cV,pV] = pp1_encoding('sim_MFPatterns',Usf,modelNames{1},modelThetas{1},'numSim',numSim,'snr',nn);
                [U_hat_nonlin,~,~]   = pp1_encoding('sim_MFPatterns',Usf,modelNames{2},modelThetas{2},'numSim',numSim,'snr',nn);
                Y_sim     = [U_hat_linear,U_hat_nonlin]; 
                trueModel = [ones(numSim,1); ones(numSim,1).*2];
                
                for is=1:numel(Y_sim) % per simulated dataset
                    % 3. split simualted data into even-odd partitions:
                    part = unique(partVec{ii});
                    for ip=1:2
                        partI{ip}=part(mod(part+ip,2)==0);
                    end
                    Dp = [];
                    for pp=1:numel(partI) 
                        dp=[];
                        dp.part = pp;
                        dp.sim  = is;
                        % Split into training and testing data: (make this
                        % clearer by assigning to appropriately named variables).
                        % Use training data to estimate true Us, then test models
                        % using estimated Us with left-out data:
                        trainIdx = ~ismember(partVec{ii},partI{pp}); 
                        testIdx = ismember(partVec{ii},partI{pp}); 
                        Y_train = Y_sim{is}(trainIdx,:);
                        Y_test  = Y_sim{is}(testIdx,:);

                        % 4. predict patterns from training data:
                        [Uest_train,G_train] = pp1_encoding('pcm_estTrueU',Y_train,pV(trainIdx),cV(trainIdx));
                        Usf_train  = Uest_train(1:5,:); % drop multi-finger patterns
                        Upred_linear  = pp1_encoding('pred_MFPatterns',Usf_train,modelNames{1},modelThetas{1});
                        Upred_tanh    = pp1_encoding('pred_MFPatterns',Usf_train,modelNames{2},modelThetas{2});
                        [Uest_test,G_test] = pp1_encoding('pcm_estTrueU',Y_test,pV(testIdx),cV(testIdx)); % estimates of Us from test data

                        % 5. do PCM likelihoods and encoding metrics per simulated dataset:
                        M{1}.Gc = Upred_linear*Upred_linear' ./ numVox; % summation G pred
                        M{2}.Gc = Upred_tanh*Upred_tanh' ./ numVox; % summation tanh G pred
                        T = pcm_fitModelIndivid({Y_test},M,pV(testIdx),cV(testIdx),'runEffect','random','verbose',0,'fitScale',1);
                        
                        % scale predicted patterns to test patterns:
                        scaleSum  = (Upred_linear(:)'*Uest_test(:))/(Upred_linear(:)'*Upred_linear(:));
                        scaleTanh = (Upred_tanh(:)'*Uest_test(:))/(Upred_tanh(:)'*Upred_tanh(:));
                        % remove mean patterns
                        Upred_linear = Upred_linear - mean(Upred_linear,1);
                        Upred_tanh   = Upred_tanh - mean(Upred_tanh,1);
                        Uest_test    = Uest_test - mean(Uest_test,1);
                        % calculate encoding metrics using centred patterns
                        [dp.r2(1,1), dp.r(1,1)] = pp1_encoding('calcModelFit',Upred_linear,Uest_test);
                        [dp.r2(1,2), dp.r(1,2)] = pp1_encoding('calcModelFit',Upred_tanh,Uest_test);
                        dp.like = T.likelihood;
                        
                        % save some info about avg. activity and avg.
                        % variances:
                        dp.act_orig  = act_orig;
                        dp.act_train = mean(C_inv * Uest_train,2)';
                        dp.act_test  = mean(C_inv * Uest_test,2)';
                        dp.act_linear= mean(C_inv * Upred_linear,2)';
                        dp.act_tanh  = mean(C_inv * Upred_tanh,2)';
                        dp.var_orig  = var_orig;
                        dp.var_train = (C_inv * diag(G_train))';
                        dp.var_test  = (C_inv * diag(G_test))';
                        dp.var_linear= (C_inv * diag(M{1}.Gc))';
                        dp.var_tanh  = (C_inv * diag(M{2}.Gc))';
                        
                        
                        Dp = addstruct(Dp,dp);
                    end
                    Dp = tapply(Dp,{'sim'},{'r2','mean'},{'r','mean'},{'like','mean'},...
                        {'act_orig','mean(x,1)'},{'act_train','mean(x,1)'},{'act_test','mean(x,1)'},{'act_linear','mean(x,1)'},{'act_tanh','mean(x,1)'},...
                        {'var_orig','mean(x,1)'},{'var_train','mean(x,1)'},{'var_test','mean(x,1)'},{'var_linear','mean(x,1)'},{'var_tanh','mean(x,1)'});
                    % 6. determine winning model (avg. metrics across CV
                    % splits)
                    evalCriterion = [Dp.r2; Dp.r; Dp.like];
                    winners       = pp1_encoding('chooseWinner',evalCriterion);
                    Dp.best_r2   = winners(1);
                    Dp.best_r    = winners(2);
                    Dp.best_like = winners(3);
                    Dp.r2_winner   = Dp.r2(Dp.best_r2);
                    Dp.r2_loser    = Dp.r2(Dp.r2~=Dp.r2_winner);
                    Dp.r_winner    = Dp.r(Dp.best_r);
                    Dp.r_loser     = Dp.r(Dp.r~=Dp.r_winner);
                    Dp.like_winner = Dp.like(Dp.best_like);
                    Dp.like_loser  = Dp.like(Dp.like~=Dp.like_winner);
                    if isempty(Dp.r2_loser)
                        Dp.r2_loser = Dp.r2_winner;
                    end
                    if isempty(Dp.r_loser)
                        Dp.r_loser = Dp.r_winner;
                    end
                    if isempty(Dp.like_loser)
                        Dp.like_loser = Dp.like_winner;
                    end
                    % indexing output fields:
                    Dp.sn   = sn(ii);
                    Dp.roi  = roi;
                    Dp.glm  = glm;
                    Dp.snr  = nn;
                    Dp.trueModel = trueModel(is);
                    
                    D = addstruct(D,Dp);
                end
            end
            fprintf('\n');
        end
    
        % evaluate metric accuracies:
        D.acc_r2   = D.best_r2==D.trueModel;
        D.acc_r    = D.best_r==D.trueModel;
        D.acc_like = D.best_like==D.trueModel;
        A = tapply(D,{'sn','roi','glm','snr'},{'acc_r2','sum'},{'acc_r','sum'},{'acc_like','sum'},...
            {'r2_winner','mean'},{'r_winner','mean'},{'like_winner','mean'},...
            {'r2_loser','mean'},{'r_loser','mean'},{'like_loser','mean'});
        A.acc_r2   = A.acc_r2 / (numSim*2);
        A.acc_r    = A.acc_r / (numSim*2);
        A.acc_like = A.acc_like / (numSim*2);
        A.logSnr   = log(A.snr);
        if ismember(0,snr)
            A.logSnr(A.snr==0) = log(snr(2))-log(snr(end));
        end

        varargout = {D,A};

    
       
    case 'chords'
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
    case 'getUsageG'
%         G=ef1_gloveData('calcDistances_chords');
%         G=mean(G.G,1);
        load(fullfile(dataDir,'chordUsageG.mat'));
        varargout = {usageG};
    case 'getData'
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
        B = load(fullfile(dataDir,sprintf('glm%d_roi%d_betas.mat',glm,roi)));
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
            bb.run   = cell2mat(b.run);
            bb.chord = cell2mat(b.tt);
            eval(sprintf('bb.betas = cell2mat(b.%s);',betaType));
            bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive stimulation conditions (tt 32 is thumb press)
            % put subj data into pcm variables
            Y{ii}         = bb.betas;
            partVec{ii}   = bb.run;
            condVec{ii}   = bb.chord;
            G_hat(:,:,ii) = pcm_estGCrossval(Y{ii},partVec{ii},condVec{ii});
        end
        varargout = {Y,partVec,condVec,G_hat};     
    case 'getData_restrictedMF'
        % Get betas for roi from subjects in PCM-friendly format.
        % IMPORTANTLY, we restrict the data to be only from voxels that have
        % significant responses to multi-finger chords (determined with an 
        % uncorrected non-crossvalidated omnibus F-test for each voxel).
        % Run means are NOT removed as this is a feature we want to retain.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        alpha = 0.05;
        % get data from all voxels
        [Yall,partVec,condVec]=pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % toss out non-sig (omnibus F-test) voxels for each subj:
        Y={};
        for ii=1:numel(Yall)
            mfIdx = condVec{ii}>5; % use multi-finger data only
            % get F-stats per voxel
            [vF,Fcrit,vP] = pp1_encoding('calcFstat',Yall{ii}(mfIdx,:),condVec{ii}(mfIdx,:),partVec{ii}(mfIdx,:));
            %[Z,F,p] = spmj_zFtest_imcalc(Yall{ii}(mfIdx,:),condVec{ii}(mfIdx,:)-5,partVec{ii}(mfIdx,:));
            % cull voxel data based on significant F-stat:
            sigIdx = vF>=Fcrit;
            Y{ii} = Yall{ii}(:,sigIdx);
        end
        varargout = {Y,partVec,condVec};    
    case 'getData_restrictedSF'
        % Get betas for roi from subjects in PCM-friendly format.
        % IMPORTANTLY, we restrict the data to be only from voxels that have
        % significant responses to single-finger stim (determined with an 
        % uncorrected non-crossvalidated omnibus F-test for each voxel).
        % Run means are NOT removed as this is a feature we want to retain.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        % get data from all voxels
        [Yall,partVec,condVec]=pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % toss out non-sig (omnibus F-test) voxels for each subj:
        Y={};
        for ii=1:numel(Yall)
            mfIdx = condVec{ii}<6; % use single-finger data only
            % get F-stats per voxel
            [vF,Fcrit,vP] = pp1_encoding('calcFstat',Yall{ii}(mfIdx,:),condVec{ii}(mfIdx,:),partVec{ii}(mfIdx,:));
            %[Z,F,p] = spmj_zFtest_imcalc(Yall{ii}(mfIdx,:),condVec{ii}(mfIdx,:)-5,partVec{ii}(mfIdx,:));
            % cull voxel data based on significant F-stat:
            sigIdx = vF>=Fcrit;
            Y{ii} = Yall{ii}(:,sigIdx);
        end
        varargout = {Y,partVec,condVec};    
    case 'getData_restricted'
        % Get betas for roi from subjects in PCM-friendly format.
        % IMPORTANTLY, we restrict the data to be only from voxels that have
        % significant responses to any condition (determined with an 
        % uncorrected non-crossvalidated omnibus F-test for each voxel).
        % Run means are NOT removed as this is a feature we want to retain.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        % get data from all voxels
        [Yall,partVec,condVec]=pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % toss out non-sig (omnibus F-test) voxels for each subj:
        Y={};
        for ii=1:numel(Yall)
            % get F-stats per voxel
            [vF,Fcrit,vP] = pp1_encoding('calcFstat',Yall{ii},condVec{ii},partVec{ii});
            % cull voxel data based on significant F-stat:
            sigIdx = vF>=Fcrit;
            Y{ii} = Yall{ii}(:,sigIdx);
        end
        varargout = {Y,partVec,condVec};    
    case 'getData_restrictedChord'
        % Get betas for roi from subjects in PCM-friendly format.
        % IMPORTANTLY, we restrict the data to be only from voxels that have
        % significant responses to multi-finger chords (determined with an 
        % uncorrected non-crossvalidated omnibus F-test for each voxel) &
        % are not significantly tuned to single finger chords. I.e., voxels
        % that are chord-specific only.
        % Run means are NOT removed as this is a feature we want to retain.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        alpha = 0.05;
        % get data from all voxels
        [Yall,partVec,condVec]=pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % toss out non-sig (omnibus F-test) voxels for each subj:
        Y={};
        for ii=1:numel(Yall)
            mfIdx = condVec{ii}>5; % use multi-finger data only
            sfIdx = condVec{ii}<6;
            % get F-stats per voxel
            [mf_vF,mf_Fcrit] = pp1_encoding('calcFstat',Yall{ii}(mfIdx,:),condVec{ii}(mfIdx,:),partVec{ii}(mfIdx,:));
            [sf_vF,sf_Fcrit] = pp1_encoding('calcFstat',Yall{ii}(sfIdx,:),condVec{ii}(sfIdx,:),partVec{ii}(sfIdx,:));
            % cull voxel data based on significant F-stat:
            sigIdx = mf_vF>=mf_Fcrit & sf_vF<sf_Fcrit;
            Y{ii} = Yall{ii}(:,sigIdx);
        end
        varargout = {Y,partVec,condVec};    
    
        
    case 'calcFstat'
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
    
    case 'do_new'
        Y = varargin{1}; 
        pV= varargin{2};
        cV= varargin{3};
        if numel(varargin)>3
            modelThetas = varargin{4}; % user passed EXACT model params
            bootstrapNull=0;
        end
        D = []; % overall output structure
        for ii=1:numel(Y)
            fprintf('%02d.',ii);% verbose output to user
            Y_pred_linear = [];
            Y_pred_tanh   = [];
            Y_pred_lnc    = [];
            G_pred_linear = [];
            G_pred_tanh   = [];
            G_pred_lnc    = [];
            G_cv_lnc      = [];
            U_overall = pp1_encoding('pcm_estTrueU',Y{ii},pV{ii},cV{ii});
            % split data into leave-one-out partitions:
            [dataTrain,dataTest] = pp1_encoding('partitionData',Y{ii},pV{ii},cV{ii});
            for pp=1:numel(dataTrain)
                % estimate true Us of training and test data:
                U_train = pp1_encoding('pcm_estTrueU',dataTrain{pp}.Y, dataTrain{pp}.pV, dataTrain{pp}.cV);
                Y_test = dataTest{pp}.Y;
                % predict data under each model using training data:
                [D_pred,M_subj] = pp1_encoding('modelPredict',U_train,Y_test,U_overall);
                Y_pred_linear = [Y_pred_linear; D_pred.U_pred{2}];
                Y_pred_tanh   = [Y_pred_tanh; D_pred.U_pred{3}];
                Y_pred_lnc    = [Y_pred_lnc; D_pred.U_pred{4}];
                G_pred_linear(:,:,pp) = M_subj{2}.Gc;
                G_pred_tanh(:,:,pp)   = M_subj{3}.Gc;
                G_pred_lnc(:,:,pp)    = M_subj{4}.Gc;
                G_cv_lnc(:,:,pp)      = pcm_estGCrossval(dataTrain{pp}.Y, dataTrain{pp}.pV, dataTrain{pp}.cV);
            end
            M_subj{2}.Gc = pcm_estGCrossval(Y_pred_linear,pV{ii},cV{ii});
            M_subj{3}.Gc = pcm_estGCrossval(Y_pred_tanh,pV{ii},cV{ii});
            M_subj{4}.Gc = pcm_estGCrossval(Y_pred_lnc,pV{ii},cV{ii});
            Y_subj.Y  = Y{ii};
            Y_subj.pV = pV{ii};
            Y_subj.cV = cV{ii};
            [D_pcm,T_subj] = pp1_encoding('modelEvaluate_pcm',M_subj,Y_subj);
            % add things to output structure
%             d=D_pred;
             v=ones(size(D_pcm.model));
%             
%             d.g_pred = D_pcm.g_pred;
%             d.like   = D_pcm.like;
            D_pcm.sn   = v.*ii;
            D=addstruct(D,D_pcm);
        end
        varargout = {D};
        
        
    case 'modelFittingWrapper'
        bootstrapNull=0;
        
        Y = varargin{1}; 
        pV= varargin{2};
        cV= varargin{3};
        if numel(varargin)>3
            modelThetas = varargin{4}; % user passed EXACT model params
            bootstrapNull=0;
        end
        D = []; % overall output structure
        for ii=1:numel(Y)
            fprintf('%02d.',ii);% verbose output to user
            U_overall = pp1_encoding('pcm_estTrueU',Y{ii},pV{ii},cV{ii});
            % split data into even-odd partitions:
            [dataTrain,dataTest] = pp1_encoding('partitionData',Y{ii},pV{ii},cV{ii});
            for pp=1:numel(dataTrain)
                % estimate true Us of training and test data:
                U_train = pp1_encoding('pcm_estTrueU',dataTrain{pp}.Y, dataTrain{pp}.pV, dataTrain{pp}.cV);
                U_test  = pp1_encoding('pcm_estTrueU',dataTest{pp}.Y, dataTest{pp}.pV, dataTest{pp}.cV);
                %U_test = dataTest{pp}.Y;
                % predict data under each model using training data:
                if ~exist('modelThetas','var')
                    [D_pred,M_subj] = pp1_encoding('modelPredict',U_train,U_test,U_overall);
                else
                    [D_pred,M_subj] = pp1_encoding('modelPredict',U_train,U_test,U_overall,modelThetas);
                end
                % do we want to add bootstrapped null summation models?
                if bootstrapNull
                    [D_pred_null,M_subj_null] = pp1_encoding('modelPredictNull',U_train,U_test);
                    % add models to other model structures:
                    D_pred_null.model = D_pred_null.model + max(D_pred.model);
                    D_pred = addstruct(D_pred,D_pred_null);
                    M_subj = [M_subj M_subj_null];
                end
                
                % evaluate model predictions against test data:
                D_fit1 = pp1_encoding('modelEvaluate_encoding',D_pred.U_pred,U_test);
                [D_fit2,T_pcm] = pp1_encoding('modelEvaluate_pcm',M_subj,dataTest{pp});
                
                if bootstrapNull
                    % need to avg. fits across different null models
                    D_fit1.model = D_pred.model;
                    D_fit2.model = D_pred.model;
                    D_fit1 = tapply(D_fit1,{'model'},{'r2','mean'},{'r','mean'});
                    D_fit2 = tapply(D_fit2,{'model'},{'g_pred','mean(x,1)'},{'like','mean'});
                    D_pred = tapply(D_pred,{'model'},{'modelThetas','mean(x,1)'},{'scale','mean(x,1)'},...
                                    {'act_pred','mean(x,1)'},{'var_pred','mean(x,1)'});
                end
                % add things to output structure
                d=D_pred;
                v=ones(numel(d.model),1);
                d.part = v.*pp;
                d.sn   = v.*ii;
                d.r    = D_fit1.r;
                d.r2   = D_fit1.r2;
                d.like = D_fit2.like;
                d.g_pred = D_fit2.g_pred;
                D=addstruct(D,d);
            end
        end
        % avg. across partitions per subject:
        D=tapply(D,{'sn','model'},{'modelThetas','nanmean(x,1)'},{'g_pred','nanmean(x,1)'},{'r2','mean'},{'r','mean'},{'like','mean'},...
            {'act_pred','mean(x,1)'},{'var_pred','mean(x,1)'});
        varargout = {D};
    case 'partitionData'
        % performs even-odd data partitioning.
        % returns data split into training and testing data structures:
        Y = varargin{1}; % data matrix of activity patterns from GLM [conds*runs x voxels]
        partVec = varargin{2}; % run assignment (this is what we use to split data) [conds*runs x 1]
        condVec = varargin{3}; % condition assignment [conds*runs x 1]
        
        % define splitting:
        splitType = 'leaveOneOut';%'leaveTwoOut';%'leaveTwoOut';
        
        % assign runs to each partition
        part = unique(partVec);
        numPart = numel(part);
        partI = {};
        switch splitType
            case 'evenOdd'
                for ip=1:2
                    partI{ip}=part(mod(part+ip,2)==0);
                end
            case 'leaveOneOut'
                for ip=1:numPart
                    partI{ip}=part(ip);
                end
            case 'leaveTwoOut'
                % leave two folds out, non-overlapping
                for ip=1:floor(numPart/2)
                    partI{ip}=part((ip-1)*2+[1:2]);
                end
            case 'leaveTwoOutOverlapping'
                % leave 2-folds out, with overlapping sets!!
                for ip=1:numPart
                    for ipp=ip+1:numPart
                        partI{end+1}=[part(ip) part(ipp)];
                    end
                end    
        end      
        % harvest data into partitions (Yes, this would not scale well for
        % bigger data sets b/c I copy the same data over and over, but I'm
        % not really concerned about that)
        for ip=1:numel(partI)
            trainIdx = ~ismember(partVec,partI{ip}); 
            testIdx  = ismember(partVec,partI{ip}); 
            % create training and test data structures:
            dataTrain{ip}.Y  = Y(trainIdx,:);
            dataTrain{ip}.pV = partVec(trainIdx,:);
            dataTrain{ip}.cV = condVec(trainIdx,:);
            
            dataTest{ip}.Y  = Y(testIdx,:);
            dataTest{ip}.pV = partVec(testIdx,:);
            dataTest{ip}.cV = condVec(testIdx,:);
        end
        % returns two data structures, each with two cells (dataTrain{1} is
        % used to train models that will be tested using dataTest{1})
        varargout = {dataTrain,dataTest};
    case 'modelPredict'
        % predicts patterns using U single finger Us from training
        % data.
        U_train = varargin{1};
        U_test  = varargin{2}; % used to scale mean of predicted patterns
        U_overall = varargin{3}; % used for upper noise ceiling
        if numel(varargin)>3
            trueModelThetas = varargin{4}; % user passed EXACT model params
        end
        % define which models we will predict data under:
        modelName   = {'null_scaling','summation','summation_tanh','lower_ceil','upper_ceil'};
        modelTheta0 = {[nan nan],[0 nan],[0 1],[nan nan],[nan nan],[nan nan]}; % theta(1) = baseline param, theta(2+)=model specific params
        if exist('trueModelThetas','var')
            % if true params passed through, assume we are testing model
            % simulated data...
            modelName   = {'summation','summation_tanh'};
            modelTheta0 = {[0 nan],[0 1]}; % theta(1) = baseline param, theta(2+)=model specific params
        end

        % housekeeping:
        C_inv = pinv(pcm_indicatorMatrix('identity',sum(pp1_encoding('chords'),2))); % to calculate mean activity or variance per # digits
        numVox = size(U_train,2);
        M_subj = {}; % pcm model structure
        D_pred = []; % output prediction structure
        G_preds= [];
        % Loop through models and do predictions:
        for mm=1:numel(modelName)
            d=[];
            % get appropriate model params and predict Us under
            % model.
            modelThetas = modelTheta0{mm};
            switch modelName{mm}
                case {'null','null_scaling','null_linear','usage'}
                    U_pred = nan(size(U_train));
                case 'lower_ceil'
                    U_pred = U_train;
                case 'upper_ceil'
                    U_pred = U_overall;
                otherwise
                    Usf_train = U_train(1:5,:);
                    if ~exist('trueModelThetas','var')
                        modelThetas = fminsearch(@(x) estModelThetas(x,Usf_train,U_train,modelName{mm}),modelTheta0{mm});
                    else
                        modelThetas = trueModelThetas{mm};
                    end
                    U_pred = pp1_encoding('pred_MFPatterns',Usf_train,modelName{mm},modelThetas);
            end
            % indexing for output struct:
            d.model       = mm;
            d.modelName   = {modelName{mm}};
            d.modelThetas = modelThetas;
            % scale mean of pred patterns to mean of test patterns
            scale = mean(U_test(:)) / mean(U_pred(:));%(U_pred(:)'*U_test(:)) / (U_pred(:)'*U_pred(:));
            U_pred = U_pred.*scale;
            d.scale = scale;
            G_pred = U_pred*U_pred' ./ numVox;
            % save some info about avg. activity and avg.
            % variances:
            d.act_pred = mean(C_inv * U_pred,2)';
            d.var_pred = (C_inv * diag(G_pred))';
            % create pcm model for null models definition:
            if strcmp(modelName{mm},'null_scaling')
                X = pp1_encoding('chords');
                %G_pred = X*eye(5)*X';
                G_pred = X*ones(5)*X';
            end
            if strcmp(modelName{mm},'null_linear')
                X = pp1_encoding('chords');
                G_pred = X*eye(5)*X';
            end
            if strcmp(modelName{mm},'null')
                G_pred = eye(31);
            end
            if strcmp(modelName{mm},'usage')
                G_pred = pp1_encoding('getUsageG');
            end
            M_subj{mm}      = pp1_encoding('prep_pcmModel');
            M_subj{mm}.name = modelName{mm};
            M_subj{mm}.Gc   = G_pred;
            G_preds(:,:,mm) = G_pred;
            
            % calculate encoding metrics:
            % add model predictions to output structure:
            d.U_pred = {U_pred};
            d.U_test = {U_test};
            d.g_pred = rsa_vectorizeIPM(G_pred);
            D_pred=addstruct(D_pred,d);
        end
        varargout = {D_pred,M_subj,G_preds};
    case 'modelPredictNull'
        % predicts patterns using U single finger Us from training
        % data. 
        % SPECIAL CASE TO DO LOTS OF NULL SUMMATION MODELS
        U_train = varargin{1};
        U_test  = varargin{2}; % used to scale mean of predicted patterns

        % housekeeping:
        C_inv = pinv(pcm_indicatorMatrix('identity',sum(pp1_encoding('chords'),2))); % to calculate mean activity or variance per # digits
        numVox = size(U_train,2);
        M_subj_null = {}; % pcm model structure
        D_pred_null = []; % output prediction structure
        % Loop through models and do predictions:
        OPT.MaxIter=10000;
        for mm=1:100
            d=[];
            % get appropriate model params and predict Us under
            % model.
            Usf_train = U_train(1:5,:);
            modelThetas = fminsearch(@(x) estModelThetas(x,Usf_train,U_train,'summation_null'),[0 nan],OPT);
            U_pred = pp1_encoding('pred_MFPatterns',Usf_train,'summation_null',modelThetas);
            d.modelThetas = modelThetas;
            % indexing for output struct:
            d.model       = 1;
            d.modelName   = {'summation_null'};
            % scale mean of pred patterns to mean of test patterns
            scale = mean(U_test(:)) / mean(U_pred(:));%(U_pred(:)'*U_test(:)) / (U_pred(:)'*U_pred(:));
            U_pred = U_pred.*scale;
            d.scale = scale;
            G_pred = U_pred*U_pred' ./ numVox;
            % save some info about avg. activity and avg.
            % variances:
            d.act_pred = mean(C_inv * U_pred,2)';
            d.var_pred = (C_inv * diag(G_pred))';
            % create pcm model for this model definition:
            M_subj_null{mm}      = pp1_encoding('prep_pcmModel');
            M_subj_null{mm}.name = 'summation_null';
            M_subj_null{mm}.Gc   = G_pred;
            
            % calculate encoding metrics:
            % add model predictions to output structure:
            d.U_pred = {U_pred};
            d.U_test = {U_test};
            d.g_pred = rsa_vectorizeIPM(G_pred);
            D_pred_null=addstruct(D_pred_null,d);
        end
        varargout = {D_pred_null,M_subj_null};
    case 'modelEvaluate_encoding'
        % evaluates predictions under models
        U_pred  = varargin{1}; % cell-array of model predictions (predictions for each model are in one cell)
        U_test  = varargin{2}; % matrix of data we are evaluating U_pred against
        D_fit = []; % output structure
        U_test_noMean = U_test-mean(U_test,1);
        for mm=1:numel(U_pred)
            d.model = mm;
            % calculate encoding metrics:
            U_pred_noMean = U_pred{mm}-mean(U_pred{mm},1); % remove means
            % R2:
            res  = U_test_noMean(:)-U_pred_noMean(:);
            rss  = sum(res.^2);
            tss1 = sum(U_test_noMean(:).*U_test_noMean(:)); % total sums of squares
            d.r2 = 1-(rss/tss1);
            % R:
            ssC  = sum(U_test_noMean(:).*U_pred_noMean(:)); % covariance of predicted with true patterns
            tss2 = sum(U_pred_noMean(:).*U_pred_noMean(:)); % predicted sums of squares
            d.r  = ssC/sqrt(tss1.*tss2);
            D_fit=addstruct(D_fit,d);
        end
        varargout = {D_fit};
    case 'modelEvaluate_pcm'
        % evaluates predictions under pcm models
        M_subj   = varargin{1}; % pcm model structure
        dataTest = varargin{2}; % matrix of data we are evaluating U_pred against
        % do PCM likelihoods for predicted models:
        [T,~,G_pcm] = pcm_fitModelIndivid({dataTest.Y},M_subj,dataTest.pV,dataTest.cV,...
            'runEffect','random','verbose',0,'maxIteration',10000,'fitScale',1);
        
        M_subj{4}.name='noiseceiling_cv'; M_subj{4}.type='noiseceiling';
        M_subj{5}.name='noiseceiling_overall'; M_subj{5}.type='noiseceiling';
%         [T_cv,DD,theta_hat,theta0] = pcm_fitModelIndividCrossval({dataTest.Y},M_subj,dataTest.pV,dataTest.cV,...
%             'runEffect','random','verbose',0,'MaxIteration',10000,'evaluation',{'R','R2','likelihood'});
        D_fit.model  = [1:numel(M_subj)]';
        D_fit.like   = T.likelihood';
        D_fit.g_pred = cell2mat(cellfun(@(x) rsa_vectorizeIPM(x),G_pcm,'uni',0)') .* T.scale';
%         D_fit.like_cv = T_cv.likelihood';
%         D_fit.r_cv    = T_cv.R';
%         D_fit.r2_cv   = T_cv.R2';
        varargout = {D_fit,T};
     
    case 'prep_pcmModel'
        M.type         = 'fixed';
        M.numGparams   = 0;
        M.name         = '';
        M.fitAlgorithm = 'NR';
        varargout = {M};
    case 'pcm_estTrueU'
        % use PCM to estimate the true Us for all chords using data from
        % all runs
        
        % inputs
        Y   = varargin{1}; % matrix of activity patterns [#conds*#runs x #vox]
        pV  = varargin{2}; % condition vector (chord #s)
        cV  = varargin{3}; % partition vector (assume cV and pV are same across subjs)
        % create PCM model structure M:
        M{1}.name = 'finger';
        M{1}.type = 'freedirect';
        M{1}.numGparams = 0;
        M{1}.theta0 = [];
        % fit model (crossvalidated second moment of data):
        [~,theta_hat,G_pred] = pcm_fitModelIndivid({Y},M,pV,cV,'runEffect','random','verbose',0,'fitScale',0);
        % reconstruct true patterns:
        [Z,B,X] = pcm_setUpFit({Y},pV,cV,'runEffect','random');
        U = pcm_estimateU(G_pred{1},theta_hat{1},Y,Z{1},X{1},'runEffect',B{1});
        G = U*U' ./ size(U,2);
        varargout = {U,G};              
    case 'sim_MFPatterns'  
        % case that accepts single finger patterns, generates multi finger
        % patterns under different models
        Usf    = varargin{1}; % single finger patterns {[numRuns*5 x numVox]}, assumes correctly ordered by condition [1;2;3;4;5;1;2;3;...]
        model  = varargin{2}; % 'summation' or 'summationNL'
        modelTheta  = varargin{3}; % theta(1) = baseline, theta(2:end) = model params
        numRun = 11;  % # runs we simulate data for
        snr    = inf; % signal-to-noise variance ratio (inf=noiseless, 0=patterns are pure i.i.d. noise)
        numSim = 1;   % # simulated datasets (use the same single finger patterns, but change the noise)
        vararginoptions(varargin(4:end),{'numRun','snr','numSim'});
        % housekeeping
        numConds  = 31; % # cconds we simulate data for
        numVox    = size(Usf,2);
        % Make true patterns according to perscribed model:
        X = kron(ones(numRun,1),eye(numConds)); % multi-run design
        if snr~=0 % simulated data has signal component
            U_true = pp1_encoding('pred_MFPatterns',Usf,model,modelTheta); % simulate data under this model
            noiseVar = var(U_true(:)) / snr; % use target snr to calc noise variance for this dataset
        elseif snr==0 % simulated data lacks signal component; patterns are purely noise
            U_true   = zeros(numConds,numVox);
            noiseVar = 1;
        end
        U_run = X*U_true; % pull true multi-finger Us through to each run
        % generate noise. Amount of noise added depends on snr variable.
        U_sim = {};
        N = size(U_run,1);
        noiseDist = @(x) norminv(x,0,1); % generate i.i.d. noise with zero mean
        for ii=1:numSim
            pNoise = unifrnd(0,1,N,numVox);            % add unique noise
            Noise  = noiseDist(pNoise)*sqrt(noiseVar); % scale noise component accordingly
            U_sim{ii} = U_run + Noise;
        end
        % update condition and partition vectors for output:
        cV = kron(ones(numRun,1),[1:numConds]');%find(numDigits>1));
        pV = kron([1:numRun]',ones(numConds,1));
        % calculate the true (noiseless) second moment for the model:
        G  = (U_true*U_true')./numVox; 
        
        varargout = {U_sim,cV,pV,U_true,G};
        % U_hat - model patterns with i.i.d. noise [numConds*numRun x numVox]
        % cV    - condition vector [numConds*numRun x 1]
        % pV    - partition vector [numConds*numRun x 1]
        % U     - true (noiseless) patterns [numConds x numVox]
        % G     - true (noiseless) second moment of multi finger patterns [numConds x numConds]
    case 'estimateVariances'
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
    case 'pred_MFPatterns'
        % factorization of encoding models:
        Ysf = varargin{1};
        model = varargin{2};
        theta = varargin{3}; % theta(1) = baseline, theta(2:end) = model params
        switch model
            case 'summation'
                % get multi-finger design:
                X = pp1_encoding('chords');
                Ymf_hat = X*(Ysf-theta(1)) + theta(1);
            case 'summation_null'
                % null linear summation model. Shuffle the single finger
                % pattern assignments and use shuffled assignments to
                % predict multi-finger patterns
                Ysf = Ysf(randperm(5),:); % shuffle fingers
                Ymf_hat = pp1_encoding('pred_MFPatterns',Ysf,'summation',theta(1)); % do linear summation
            case 'summationNL_OLD'
                % y = x^theta
                Ymf_hat = pp1_encoding('pred_MFPatterns',Ysf,'summation',theta(1)); % do summation
                % estimate theta
                %[theta,err,gEst] = fminsearch(@(x) fitPower(x,Ymf_hat(:),Ymf(:)),0.5);
                % squish patterns with theta param:
                sign_Uhat= sign(Ymf_hat);
                abs_Uhat = abs(Ymf_hat);
                Ymf_hat  = sign_Uhat.*(abs_Uhat.^theta(2));
            case 'summation_tanh'
                % y = tanh(x*theta)/theta
                % theta 2 is a scalar param to scale the patterns, apply
                % nonlinearity, then unscale
                Ymf_hat = pp1_encoding('pred_MFPatterns',Ysf,'summation',theta(1)); % do linear summation
                % nonlinearize patterns:
                Ymf_hat = (tanh(Ymf_hat.*theta(2)))./theta(2);
        end
        varargout = {Ymf_hat};
    case 'modelSelectionAccuracy'
        % case to calculate model selection accuracies
        % sloppy case, but works
        D = varargin{1}; % model fitted structure. output from 'modelFittingWrapper'
            % also has field called "trueModel" that defines which is the
            % true model for this simulated subject
        numDatasets = numel(unique(D.sn));
        A = [];
        for ss=1:numDatasets
            d = getrow(D,D.sn==ss);
            
            evalCriterion = [d.r2'; d.r'; d.like'];
            winners       = pp1_encoding('chooseWinner',evalCriterion); % this case can handle ties
            a.best_r2    = winners(1);
            a.best_r     = winners(2);
            a.best_like  = winners(3);
            a.r2_winner   = d.r2(a.best_r2);
            a.r2_loser    = d.r2(d.r2~=a.r2_winner);
            a.r_winner    = d.r(a.best_r);
            a.r_loser     = d.r(d.r~=a.r_winner);
            a.like_winner = d.like(a.best_like);
            a.like_loser  = d.like(d.like~=a.like_winner);
            % in case of tie...
            if isempty(a.r2_loser)
                a.r2_loser = a.r2_winner;
            end
            if isempty(a.r_loser)
                a.r_loser = a.r_winner;
            end
            if isempty(a.like_loser)
                a.like_loser = a.like_winner;
            end
            % indexing output fields:
            a.sn   = ss;
            a.trueModel = d.trueModel(1);
            A = addstruct(A,a);
        end
        
        % evaluate metric accuracies:
        A.acc_r2   = A.best_r2==A.trueModel;
        A.acc_r    = A.best_r==A.trueModel;
        A.acc_like = A.best_like==A.trueModel;
        A = tapply(A,{'trueModel'},{'acc_r2','sum'},{'acc_r','sum'},{'acc_like','sum'},...
            {'r2_winner','mean'},{'r_winner','mean'},{'like_winner','mean'},...
            {'r2_loser','mean'},{'r_loser','mean'},{'like_loser','mean'});
        A.acc_r2   = A.acc_r2 / (numDatasets/2); % half the datasets come from each model
        A.acc_r    = A.acc_r / (numDatasets/2);
        A.acc_like = A.acc_like / (numDatasets/2);
        
        varargout = {A};
    case 'chooseWinner'
        % case to choose winning model based on evaluation criterion.
        % if the best fitting models are tied, a random winner is choosen.
        % Can handle matrix input. Each col are different model fits, and
        % each row corresponds to one model evaluation metric.
        criterion = varargin{1};
        [c,m]  = size(criterion);
        winner = zeros(c,1);
        for ii=1:c % for each evaluation metric:
            % find best model:
            idx = criterion(ii,:)==max(criterion(ii,:)); 
            if sum(idx)>1
                % tiebreaker:
                ms = randperm(m);
                winner(ii,1) = ms(1);
            else
                % mark best model of this criterion:
                winner(ii,1) = find(idx);
            end
            
        end
        varargout = {winner};    

        
    case 'do_predSubj'
        % case to fit models to participant data
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        vararginoptions(varargin,{'roi'});
        
        % get data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        %[Y,pV,cV] = pp1_encoding('getData_restrictedChord','sn',sn,'roi',roi,'glm',glm);
        
        % do fitting:
        %D = pp1_encoding('modelFittingWrapper',Y,pV,cV);
        D = pp1_encoding('do_new',Y,pV,cV);
        
        % add appropriate indexing labels to each row of output data
        % structure
        v = ones(size(D.sn));
        D.sn  = pcm_indicatorMatrix('identity',D.sn)*sn';
        D.roi = v.*roi;
        D.glm = v.*glm;
        
        varargout = {D};
    case 'do_predSim'
        % case to fit models to simulated data
        snr = [0,0.01,0.1,1,10];
        glm = 4;
        numSim = 50; % sims per subject per snr lvl under each model (2 models = 50 simulations * # snr lvls)
%         
%         modelThetas = {[-0.0096875 nan],[-0.021841 11.07]}; sn=2; roi=2;
%         modelThetas = {[-0.00625 nan],[-0.015192 10.988]}; sn=3; roi=2;
%         modelThetas = {[-0.0053438 nan],[-0.013145 2.9121]}; sn=4; roi=2;
%         modelThetas = {[-0.018719 nan],[-0.057555 5.5776]}; sn=5; roi=2;
%         modelThetas = {[-0.0050938 nan],[-0.01609 7.1285]}; sn=6; roi=2;
%         modelThetas = {[-0.0185 nan],[-0.04224 2.6669]}; sn=7; roi=2;
%         modelThetas = {[-0.014335 nan],[-0.030332 10.004]}; sn=8; roi=2;
         modelThetas = {[-0.020719 nan],[-0.021714 1.5596]}; sn=9; roi=2;
%        modelThetas = {[-0.010281 nan],[-0.048676 13.554]}; sn=10; roi=2;
%         modelThetas = {[-0.014438 nan],[-0.035568 4.721]}; sn=11; roi=2;
        
        
        
        % get data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        
        U = pp1_encoding('pcm_estTrueU',Y{1},pV{1},cV{1});
        Usf = U(1:5,:);
        
        Ao=[];
        Do=[];
        for jj=1:numel(snr)
            nn=snr(jj);
            fprintf('%2.2f...',nn);% verbose output to user
            % simulate data under both models using Usf:            
            [U_hat_linear,cV,pV] = pp1_encoding('sim_MFPatterns',Usf,'summation',modelThetas{1},'numSim',numSim,'snr',nn);
            [U_hat_nonlin,~,~]   = pp1_encoding('sim_MFPatterns',Usf,'summation_tanh',modelThetas{2},'numSim',numSim,'snr',nn);
            Y_sim     = [U_hat_linear,U_hat_nonlin]; 
            for ii=1:numSim*2
                cV_sim{ii} = cV;
                pV_sim{ii} = pV;
            end
            trueModel = [ones(numSim,1); ones(numSim,1).*2];
            
            % check snr of simulated datasets:
            varE = [];
            varS = [];
            for ii=1:numel(Y_sim)
                [varE(ii,1),varS(ii,1)]=pp1_encoding('estimateVariances',Y_sim{ii},cV_sim{ii},pV_sim{ii});
            end
            snrEst = varS./varE;
            
            % do fitting:
            D = pp1_encoding('modelFittingWrapper',Y_sim,pV_sim,cV_sim,modelThetas);

            % calculate accuracies:
            D.snrEst    = pcm_indicatorMatrix('identity',D.sn) * snrEst;
            D.trueModel = pcm_indicatorMatrix('identity',D.sn) * trueModel;
            A = pp1_encoding('modelSelectionAccuracy',D);
            
            D.snr = ones(size(D.sn)).*snr(jj);
            A.snr = ones(size(A.acc_like)).*snr(jj);
            
            Do=addstruct(Do,D);
            Ao=addstruct(Ao,A);
            fprintf('\n')
        end
        Ao.logSnr = log(Ao.snr);
        Ao.logSnr(Ao.snr==0)=-6.9078;
        
        % make plotting-friendly accuracy structure:
        T=[];
        v=ones(3,1);
        for ii=1:numel(Ao.snr)
            t.trueModel=v.*Ao.trueModel(ii);
            t.acc = [Ao.acc_r2(ii); Ao.acc_r(ii); Ao.acc_like(ii)];
            t.metric = [1:3]';
            t.snr = v.*Ao.snr(ii);
            T=addstruct(T,t);
        end
        T.logSnr = log(T.snr);
        T.logSnr(T.snr==0)=-6.9078;
        
        % plot
        clrs = {[0 0 0.5625],[0 0.875 1],[1 0.8125 0]};
        CAT.markercolor=[0 0 0];
        CAT.markerfill=clrs;
        CAT.markersize=8;
        CAT.linecolor=clrs;
        lineplot(T.logSnr,T.acc,'plotfcn','nanmean','errorbars','none','split',T.metric,'CAT',CAT,'leg',{'r2','r','likelihood'});
        ylabel('model selection accuracy');
        xlabel('SNR (log scale)');
        set(gca,'xticklabel',unique(T.snr));
        
        
        varargout = {Ao,Do,T};
    
    case 'plot_pred'
        % case to plot activity and variance (per num digits) of original
        % data, and model predictions fitted to data
        D=varargin{1}; % ouptut structure from 'do_pred'
        nNull = 1; % which model # is the null model?
        nLC   = 4; % which model # is the lower noise ceiling?
        modelNames = {'nullScale','linear','tanh','lower NC','upper NC','usage'};
        % make plotting structure:
        Q=[];
        v=ones(5,1);
        for qq=1:size(D.sn,1)
            q.sn = v.*D.sn(qq);
            q.numD  = [1:5]';
            q.act = D.act_pred(qq,:)';
            q.var = D.var_pred(qq,:)';
            q.model = v.*D.model(qq);
%             for jj=1:3
%                 q.sn    = v.*D.sn(qq);
%                 q.numD  = [1:5]';
%                 if jj==1
%                     q.act = D.act_train(qq,:)';
%                     q.var = D.var_train(qq,:)';
%                     q.model = v.*1 + max(D.model);
%                 elseif jj==2
%                     q.act = D.act_test(qq,:)';
%                     q.var = D.var_test(qq,:)';
%                     q.model = v.*2 + max(D.model);
%                 elseif jj==3
%                     q.act = D.act_pred(qq,:)';
%                     q.var = D.var_pred(qq,:)';
%                     q.model = v.*D.model(qq);
%                 end
                Q=addstruct(Q,q);
%             end
        end
%         Q=tapply(Q,{'sn','numD','model'},{'act','mean'},{'var','mean'});
        
        % plot avgs model G:
        numModels = numel(unique(Q.model));
        for mm=1:numModels
            subplot(2,numModels,mm);
            G_pred = rsa_squareIPM(mean(D.g_pred(D.model==mm,:)));
            imagesc(G_pred);
            title([modelNames{mm} ' model']);
        end
        % plot avg. pred activity per # digits in each chord
        subplot(2,numModels,numModels+1);
        clrs = {[0 0 0],[0.71 0 0],[1 0.375 0],[1 1 0.0625],[0.7 0.7 0.7]};
        CAT.markercolor=[0 0 0];
        CAT.markerfill=clrs;
        CAT.markersize=8;
        CAT.linecolor=clrs;
        lineplot(Q.numD,Q.act,'split',Q.model,'plotfcn','nanmean',...
            'errorbars','plusminus','subset',Q.model,'CAT',CAT,'leg',modelNames);
        xlabel('# digits in chord');
        ylabel('avg. activity');
        title('avg. activity predictions');
        % plot avg. pred variance per # digits in each chord
        subplot(2,numModels,numModels+2);
        lineplot(Q.numD,Q.var,'split',Q.model,'plotfcn','nanmean',...
            'errorbars','plusminus','subset',Q.model,'CAT',CAT,'leg',modelNames);
        xlabel('# digits in chord');
        ylabel('avg. variance');
        title('avg. variance predictions');
        
        
        % plot model likelihoods:
        D.likeNorm = D.like-kron(D.like(D.model==nNull),ones(numModels,1));
        %D.likeNorm = D.likeNorm ./ kron(D.likeNorm(D.model==nLC+1),ones(numModels,1)); % normalize to upper noise ceil
        subplot(2,numModels,numModels*2);
        plt.bar(D.model,D.likeNorm);%,'plotfcn','median');%,'subset',D.model<nLC & D.model>nNull);
        %plt.box(D.model,D.likeNorm,'plotall',2);
        ylabel('log BF (vs. null model)');
        set(gca,'xticklabel',modelNames,'xticklabelrotation',45);
        % Draw a patch encompassing area between upper and lower noise ceilings
        drawline(mean(D.likeNorm(D.model==nLC)),'dir','horz','linestyle',':');
        drawline(mean(D.likeNorm(D.model==nLC+1)),'dir','horz','linestyle','-');
        varargout = {Q,D};
    
        

    case 'SEARCH_define'                                                    % STEP 4.1   :  Defines searchlights for 120 voxels in grey matter surface
        caretDir = '/Users/sarbuckle/DATA/passivePatterns1/fmri/surfaceCaret';
        glmDir = '/Users/sarbuckle/DATA/passivePatterns1/fmri/glm4';
        sn  = 4;
        vararginoptions(varargin,{'sn','glm'});
        
        mask = fullfile(glmDir,subj_name{sn},'mask.nii');
        mask = rsa.readMask(mask);

        LcaretDir = fullfile(caretDir,sprintf('x%s',subj_name{sn}),'LeftHem');
        RcaretDir = fullfile(caretDir,sprintf('x%s',subj_name{sn}),'RightHem');
        white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
        pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
        topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};
        S         = rsa_readSurf(white,pial,topo);

        L = rsa.defineSearchlight_surface(S,mask,'sphere',[15 160]);
        save(fullfile(glmDir,subj_name{sn},sprintf('%s_PCM_searchlight_160.mat',subj_name{sn})),'-struct','L');
    case 'SEARCH_run'                                                    % STEP 4.2   :  Runs LDC searchlight using defined searchlights (above)
        glmDir = '/Users/sarbuckle/DATA/passivePatterns1/fmri/glm4';
        % Requires java functionality unless running on SArbuckle's
        % computer.
        sn = 4;
        vararginoptions(varargin,{'sn','glm'});
        cwd = pwd;                                                      
        % go to subject's glm directory 
        cd(fullfile(glmDir,subj_name{sn}));
        % load their searchlight definitions and SPM file
        Searchlight = load(sprintf('%s_PCM_searchlight_160.mat',subj_name{sn}));
        load SPM;
        %SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{sn}));
        for ii=1:5
            outNames{ii} = fullfile(glmDir,subj_name{sn},sprintf('%s_PCM_glm4,%d',subj_name{sn},ii));
        end
        % define options:
        Opt.pV = kron([1:11]',ones(32,1));
        Opt.cV = kron(ones(11,1),[1:32]');
        Opt.idealBlock = 5e7;
        Opt.java = true;
        % run the searchlight
        rsa.runSearchlight(Searchlight,SPM.xY.P,outNames,@pcmModelFunc,'optionalParams',{SPM,Opt.pV,Opt.cV},'java',Opt.java,'idealBlock',Opt.idealBlock);
        cd(cwd);
    
    case 'modelFittingWrapper_searchlight'
        % case to call when doing searchlight PCM analysis
        % output is [1 x numModels] vector of model log likelihoods.
        % output(1) is null model likelihood
        % output(end-1) is lower-noise-ceiling
        % output(end) is upper-noise-ceiling
        Y = varargin{1}; % noise-normalized betas
        pV= varargin{2}; % partition vector
        cV= varargin{3}; % condition vector
        
        % estimate true Us using all data (for upper noise ceiling model)
        U_overall = pp1_encoding('pcm_estTrueU',Y,pV,cV);
        % split data into leave-two-out partitions:
        [dataTrain,dataTest] = pp1_encoding('partitionData',Y,pV,cV);
        numPart = numel(dataTrain);
        logLikelihood = nan(5,numPart);
        for pp=1:numPart
            % estimate true Us of training and test data:
            U_train = pp1_encoding('pcm_estTrueU',dataTrain{pp}.Y, dataTrain{pp}.pV, dataTrain{pp}.cV);
            U_test  = pp1_encoding('pcm_estTrueU',dataTest{pp}.Y, dataTest{pp}.pV, dataTest{pp}.cV);
            % predict data under each model using training data:
            [~,M_subj] = pp1_encoding('modelPredict',U_train,U_test,U_overall);
            d = pp1_encoding('modelEvaluate_pcm',M_subj,dataTest{pp});
            logLikelihood(pp,:) = d.like';
        end
        logLikelihood = sum(logLikelihood,1);
        varargout = {logLikelihood};



    case '0' % clean encoding cases below:
    case 'ENC:do'
        % case to fit models to participant data
        verbose=0;
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        vararginoptions(varargin,{'roi','verbose'});
        
        D=[]; % output
        chords = pp1_encoding('chords');
        splitType = 'leaveOneOut';
        patternModel = {'summation','summation_tanh','summation_flexible'}; % models where we actually estimate the patterns
        modelTheta0  = {[0],[0 1],[0 1 1/2 1/3 1/4 1/5]};
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % loop through subjects and fit each individually:
        for s=1:numel(sn)
            % split data into partitions (for estimation of patterns under each model)
            % assign runs to each partition
            part = unique(pV{s});
            numPart = numel(part);
            partI = {};
            switch splitType
                case 'evenOdd'
                    for ip=1:2
                        partI{ip}=part(mod(part+ip,2)==0);
                    end
                case 'leaveOneOut'
                    for ip=1:numPart
                        partI{ip}=part(ip);
                    end
                case 'leaveTwoOut'
                    % leave two folds out, non-overlapping
                    for ip=1:floor(numPart/2)
                        partI{ip}=part((ip-1)*2+[1:2]);
                    end
            end     
            % loop through partitions and estimate patterns
            Y_hat = nan([size(Y{s}), numel(patternModel)+1]);
            for ii=1:numel(partI)
                % estimate the true condition activity patterns using
                % training data:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii}); 
                U_est = pp1_encoding('ENC:estimateU',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx));
                % using estimated single-finger Us, predict data under each
                % model:
                for mm=1:numel(patternModel)
                    % estimate the model thetas:
                    fcn = @(x) estModelThetas(x,U_est(1:5,:),U_est,patternModel{mm});
                    [modelTheta,rss] = fminsearch(fcn, modelTheta0{mm}, optimset('MaxIter',10000));
                    Y_hat(testIdx,:,mm) = pp1_encoding('ENC:predictModelPatterns',U_est(1:5,:),patternModel{mm},modelTheta);
                    if verbose
                        Ytest = Y{s}(testIdx,:);
                        tss   = sum(sum((U_est-mean(U_est,1)).^2));
                        Upred = Y_hat(testIdx,:,mm);
                        fprintf('\nsn%d ii%d model%d | training R2: %1.3f | training R: %1.3f | test R: %1.3f | test R / relaiblity: %1.3f',...
                            sn(s),ii,mm,1-(rss/tss), corr(U_est(:),Upred(:)), corr(Ytest(:),Upred(:)), corr(Ytest(:),Upred(:)) / corr(U_est(:),Ytest(:)));
                    end
                end
                Y_hat(testIdx,:,end) = U_est; % patterns for lower bound of noise ceiling
            end
            
            % define PCM models:
            M = {};
            % null scaling model (all patterns identical, mean acitvity
            % scales)
            M{1}.type         = 'component'; % to allow G to scale across runs
            M{1}.numGparams   = 1;
            M{1}.name         = 'null_scaling';
            M{1}.fitAlgorithm = 'NR';
            M{1}.Gc           = chords*ones(5)*chords';
            % summation model
            M{2}.type         = 'component'; % to allow G to scale across runs
            M{2}.numGparams   = 1;
            M{2}.name         = 'summation';
            M{2}.fitAlgorithm = 'NR';
            M{2}.Gc           = pcm_estGCrossval(Y_hat(:,:,1),pV{s},cV{s});
            % summation tanh model
            M{3}.type         = 'component';
            M{3}.numGparams   = 1;
            M{3}.name         = 'summation_tanh';
            M{3}.fitAlgorithm = 'NR';
            M{3}.Gc           = pcm_estGCrossval(Y_hat(:,:,2),pV{s},cV{s});
            % summation tanh model
            M{4}.type         = 'component';
            M{4}.numGparams   = 1;
            M{4}.name         = 'noiseceiling_cv';
            M{4}.fitAlgorithm = 'NR';
            M{4}.Gc           = pcm_estGCrossval(Y_hat(:,:,4),pV{s},cV{s});
            % lower bound noise ceiling
            M{end+1}.type       = 'component';
            M{end}.numGparams   = 1;
            M{end}.name         = 'lower noise ceiling';
            M{end}.fitAlgorithm = 'NR';
            M{end}.Gc           = pcm_estGCrossval(Y_hat(:,:,4),pV{s},cV{s});
            % overall noise ceiling
            M{end+1}.type       = 'component'; % to allow G to scale across runs
            M{end}.numGparams   = 1;
            M{end}.name         = 'upper noise ceiling';
            M{end}.fitAlgorithm = 'NR';
            M{end}.Gc           = pcm_estGCrossval(Y{s},pV{s},cV{s});
            
            [T_cv,D_cv,theta_hat,theta0] = pcm_fitModelIndividCrossval(Y{s},M,pV{s},cV{s},...
                                        'runEffect','random','verbose',0,...
                                        'evaluation',{'R','Rpool','R2','likelihood'});
            % arrange data into output structure:
            v = ones(numel(M),1);
            d = [];
            d.model  = [1:numel(v)]';
            d.roi    = v.*roi;
            d.sn     = v.*sn(s);
            % get evaluation metrics
            for mm=1:numel(M)
                d.like_cv(mm,1)  = T_cv.likelihood(mm);
                d.r2_cv(mm,1)    = T_cv.R2(mm);
                d.r_cv(mm,1)     = T_cv.R(mm);
                d.rpool_cv(mm,1) = T_cv.Rpool(mm);
            end
            D=addstruct(D,d);

        end
        
        varargout = {D};
    case 'ENC:test'
        % case to calculate model selection accuracies using different
        % evaluation metrics
        sn = varargin{1};
        % get random subject data to use:
        roi = 2;
        glm = 4;
        switch sn
            case 2
                modelTheta = {[-0.0062188],[-0.009432 5.0292]}; % provide "real" model thetas
            case 3
                modelTheta = {[-0.0042557],[-0.010482 2.0918]};
            case 4
                modelTheta = {[-0.0042557],[-0.010482 2.0918]};
            case 5
                modelTheta = {[-0.0042557],[-0.010482 2.0918]};
            case 6
                modelTheta = {[-0.006733],[-0.01748 5.3082]};
            case 7
                modelTheta = {[-0.023446],[-0.039894 1.8946]};
            case 8
                modelTheta = {[-0.018398],[-0.020124 6.3584]};
            case 9
                modelTheta = {[-0.025324],[-0.025713 0.97957]};
            case 10
                modelTheta = {[-0.014778],[-0.047002 17.576]};
            case 11
                modelTheta = {[-0.019571],[-0.031856 4.6297]};
        end    
        
        A=[];
        D=[]; % output

        splitType = 'leaveOneOut';
        patternModel = {'summation','summation_tanh'}; % models where we actually estimate the patterns
        
        numModels = numel(patternModel);
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % estimate true Us for use in simulating patterns:
        U_est = pp1_encoding('ENC:estimateU',Y{1},pV{1},cV{1});
        
                
        snr = [0:0.1:1];
        for nn=snr
            % simulate data under both models using Usf: 
            numSim = 25; % sims per subject per snr lvl under each model (2 models = 50 simulations * # snr lvls)
            [U_hat_linear,cV,pV] = pp1_encoding('ENC:simulatePatterns',U_est(1:5,:),'summation',modelTheta{1},'numSim',numSim,'snr',nn);
            [U_hat_nonlin,~,~]   = pp1_encoding('ENC:simulatePatterns',U_est(1:5,:),'summation_tanh',modelTheta{2},'numSim',numSim,'snr',nn);
            Y_sim     = [U_hat_linear,U_hat_nonlin]; 
            trueModel = [ones(numSim,1); ones(numSim,1).*2];

            % measure snr of simulated datasets:
            varE = [];
            varS = [];
            for ii=1:numel(Y_sim)
                [varE(ii,1),varS(ii,1)]=pp1_encoding('estimateVariances',Y_sim{ii},cV,pV);
            end
            snrEst = varS./varE;
            
            % do model fitting (code copied directly from 'ENC:do':
            % split data into partitions (for estimation of patterns under each model)
            % assign runs to each partition
            part = unique(pV);
            numPart = numel(part);
            partI = {};
            switch splitType
                case 'evenOdd'
                    for ip=1:2
                        partI{ip}=part(mod(part+ip,2)==0);
                    end
                case 'leaveOneOut'
                    for ip=1:numPart
                        partI{ip}=part(ip);
                    end
                case 'leaveTwoOut'
                    % leave two folds out, non-overlapping
                    for ip=1:floor(numPart/2)
                        partI{ip}=part((ip-1)*2+[1:2]);
                    end
            end     
            for s=1:numel(Y_sim)
                % loop through partitions and estimate patterns
                Y_hat = nan([size(Y_sim{s}), numModels]);
                for ii=1:numel(partI)
                    % estimate the true condition activity patterns using
                    % training data:
                    trainIdx = ~ismember(pV,partI{ii}); 
                    testIdx  = ismember(pV,partI{ii}); 
                    U_est = pp1_encoding('ENC:estimateU',Y_sim{s}(trainIdx,:),pV(trainIdx),cV(trainIdx));
                    % using estimated single-finger Us, predict data under each
                    % model:
                    for mm=1:numModels
                        % estimate the patterns using the provided model thetas:
                        Y_hat(testIdx,:,mm) = pp1_encoding('ENC:predictModelPatterns',U_est(1:5,:),patternModel{mm},modelTheta{mm});
                    end
                end

                % define PCM models:
                M = {};
                % summation model
                M{1}.type         = 'component'; % to allow G to scale across runs
                M{1}.numGparams   = 1;
                M{1}.name         = 'summation';
                M{1}.fitAlgorithm = 'NR';
                M{1}.Gc           = pcm_estGCrossval(Y_hat(:,:,1),pV,cV);
                % summation tanh model
                M{2}.type         = 'component';
                M{2}.numGparams   = 1;
                M{2}.name         = 'summation_tanh';
                M{2}.fitAlgorithm = 'NR';
                M{2}.Gc           = pcm_estGCrossval(Y_hat(:,:,2),pV,cV);
                % do pcm fitting:
                T_cv = pcm_fitModelIndividCrossval(Y_sim{s},M,pV,cV,...
                          'runEffect','random','verbose',0,...
                          'evaluation',{'R','Rpool','R2','likelihood'});
                % determine winning models based on evaluation criterias:
                criterion = [T_cv.R; T_cv.Rpool; T_cv.R2; T_cv.likelihood];
                winners = pp1_encoding('ENC:chooseWinner',criterion);
                
                % save to output structure:
                numCriteria = 4;
                v = ones(numModels*numCriteria,1);
                d = [];
                d.snr = v.*nn;
                d.snrEst = v.*mean(snrEst);
                d.simNum = v.*s;
                d.modelTrue = v.*trueModel(s);
                d.modelFit  = kron(ones(numCriteria,1),[1:numModels]');
                d.modelWin  = kron(winners,ones(numModels,1));
                d.metricType = kron([1:numCriteria]',ones(numModels,1));
                d.metricValue= [T_cv.R, T_cv.Rpool, T_cv.R2, T_cv.likelihood]';
                D=addstruct(D,d);
            end
        end
        
        % calculate model selection accuracies:
        D.accurate = (D.modelTrue==D.modelWin) ./ numModels; % divide by numModels to ensure each fitted simulated dataset is only counted once
        A = tapply(D,{'snr','metricType'},{'metricValue','mean'},{'snrEst','mean'},{'accurate','sum'});
        A.accurate = A.accurate./(numSim*numModels);
        
        varargout = {A,D};
    case 'ENC:estimateU'
        % use PCM to estimate the true Us for all chords using data from
        % all runs
        
        % inputs
        Y   = varargin{1}; % matrix of activity patterns [#conds*#runs x #vox]
        pV  = varargin{2}; % condition vector (chord #s)
        cV  = varargin{3}; % partition vector (assume cV and pV are same across subjs)
        % create PCM model structure M:
        M{1}.name = 'finger';
        M{1}.type = 'freedirect';
        M{1}.numGparams = 0;
        M{1}.theta0 = [];
        % fit model (crossvalidated second moment of data):
        [~,theta_hat,G_pred] = pcm_fitModelIndivid({Y},M,pV,cV,'runEffect','random','verbose',0,'fitScale',0);
        % reconstruct true patterns:
        [Z,B,X] = pcm_setUpFit({Y},pV,cV,'runEffect','random');
        U = pcm_estimateU(G_pred{1},theta_hat{1},Y,Z{1},X{1},'runEffect',B{1});
        G = U*U' ./ size(U,2);
        varargout = {U,G};              
    case 'ENC:predictModelPatterns'
        % factorization of encoding models:
        Ysf = varargin{1};
        model = varargin{2};
        theta = varargin{3}; % theta(1) = baseline, theta(2:end) = model params
        switch model
            case 'summation'
                % get multi-finger design:
                X = pp1_encoding('chords');
                Ymf_hat = X*(Ysf-theta(1)) + theta(1);
            case 'summation_null'
                % null linear summation model. Shuffle the single finger
                % pattern assignments and use shuffled assignments to
                % predict multi-finger patterns
                Ysf = Ysf(randperm(5),:); % shuffle fingers
                Ymf_hat = pp1_encoding('pred_MFPatterns',Ysf,'summation',theta(1)); % do linear summation
            case 'summation_tanh'
                % y = tanh(x*theta)/theta
                % theta 2 is a scalar param to scale the patterns, apply
                % nonlinearity, then unscale
                Ymf_hat = pp1_encoding('pred_MFPatterns',Ysf,'summation',theta(1)); % do linear summation
                % nonlinearize patterns:
                Ymf_hat = (tanh(Ymf_hat.*theta(2)))./theta(2);
            case 'summation_flexible'
                % theta(1) = baseline param
                % theta(2:6) = finger scalar params, one per # fingers stimulated
                % get multi-finger design:
                X = pp1_encoding('chords');
                numD = sum(X,2);
                X = X.*theta(numD+1)'; % flexible scaling per # fingers stimulated
                Ymf_hat = X*(Ysf-theta(1)) + theta(1); 
        end
        varargout = {Ymf_hat};
    case 'ENC:simulatePatterns'  
        % case that accepts single finger patterns, generates multi finger
        % patterns under different models
        Usf    = varargin{1}; % single finger patterns {[numRuns*5 x numVox]}, assumes correctly ordered by condition [1;2;3;4;5;1;2;3;...]
        model  = varargin{2}; % 'summation' or 'summationNL'
        modelTheta  = varargin{3}; % theta(1) = baseline, theta(2:end) = model params
        numRun = 11;  % # runs we simulate data for
        snr    = inf; % signal-to-noise variance ratio (inf=noiseless, 0=patterns are pure i.i.d. noise)
        numSim = 1;   % # simulated datasets (use the same single finger patterns, but change the noise)
        vararginoptions(varargin(4:end),{'numRun','snr','numSim'});
        % housekeeping
        numConds  = 31; % # cconds we simulate data for
        numVox    = size(Usf,2);
        % Make true patterns according to perscribed model:
        X = kron(ones(numRun,1),eye(numConds)); % multi-run design
        if snr~=0 % simulated data has signal component
            U_true = pp1_encoding('ENC:predictModelPatterns',Usf,model,modelTheta); % simulate data under this model
            noiseVar = var(U_true(:)) / snr; % use target snr to calc noise variance for this dataset
        elseif snr==0 % simulated data lacks signal component; patterns are purely noise
            U_true   = zeros(numConds,numVox);
            noiseVar = 1;
        end
        U_run = X*U_true; % pull true multi-finger Us through to each run
        % generate noise. Amount of noise added depends on snr variable.
        U_sim = {};
        N = size(U_run,1);
        noiseDist = @(x) norminv(x,0,1); % generate i.i.d. noise with zero mean
        for ii=1:numSim
            pNoise = unifrnd(0,1,N,numVox);            % add unique noise
            Noise  = noiseDist(pNoise)*sqrt(noiseVar); % scale noise component accordingly
            U_sim{ii} = U_run + Noise;
        end
        % update condition and partition vectors for output:
        cV = kron(ones(numRun,1),[1:numConds]');%find(numDigits>1));
        pV = kron([1:numRun]',ones(numConds,1));
        % calculate the true (noiseless) second moment for the model:
        G  = (U_true*U_true')./numVox; 
        
        varargout = {U_sim,cV,pV,U_true,G};
        % U_hat - model patterns with i.i.d. noise [numConds*numRun x numVox]
        % cV    - condition vector [numConds*numRun x 1]
        % pV    - partition vector [numConds*numRun x 1]
        % U     - true (noiseless) patterns [numConds x numVox]
        % G     - true (noiseless) second moment of multi finger patterns [numConds x numConds]
    case 'ENC:estimateSNR'
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
    case 'ENC:chooseWinner'
        % case to choose winning model based on evaluation criterion.
        % if the best fitting models are tied, a random winner is choosen.
        % Can handle matrix input. Each col are different model fits, and
        % each row corresponds to one model evaluation metric.
        criterion = varargin{1};
        [c,m]  = size(criterion);
        winner = zeros(c,1);
        for ii=1:c % for each evaluation metric:
            % find best model:
            idx = criterion(ii,:)==max(criterion(ii,:)); 
            if sum(idx)>1
                % tiebreaker:
                ms = randperm(m);
                winner(ii,1) = ms(1);
            else
                % mark best model of this criterion:
                winner(ii,1) = find(idx);
            end
            
        end
        varargout = {winner};    
    
    end
end
% % % partial derivatives of tanh model thetas wrt L2 loss.
% % %   Ymf = tanh((W*(B-theta1)+theta1).*theta2) / theta2;
% % % To simplfy a bit, the linear summation portion of the model will be
% % % substitued with X:
% % %   X = (W*(B-theta1)+theta1);
% % 
% % X = (W*(B-theta1)+theta1); % do linear summation
% % Ypred = tanh(X.*theta2) / theta2; % apply nonlinearity (squishing)
% % % dtheta1 = (1-W) * sech(X.*theta2).^2; % partial derivative of baseline parameter
% % % dtheta2 = (X * sech(X.*theta2).^2 - tanh(X.*theta2)) ./ theta2; % partial derivative of scaling parameter
% % res = Y-Ypred;
% % L = sum(sum(res.^2)); % L2 loss (RSS)
% % dLdtheta1 = -2 * (1-W) * sech(X.*theta2).^2 * res;
% % dLdtheta2 = -2 * (tanh(X.*theta2)./theta2^2 - (X*sech(X.*theta2).^2)./theta2) * res;

function err = estModelThetas(theta,Usf_train,Utrain,modelName)
Upred = pp1_encoding('ENC:predictModelPatterns',Usf_train,modelName,theta); % predict patterns under perscribed model
err   = sum(sum((Utrain-Upred).^2)); % calculate L2 loss (RSS)
end

function logLikelihoods = pcmModelFunc(Y,SPM,pV,cV)
Y = rsa.spm.noiseNormalizeBeta(Y,SPM);  % Get noise normalised betas 
idx = cV<32; % keep only betas related to stimulation
Y = Y(idx,:); 
pV = pV(idx); % for consistency
cV = cV(idx);
logLikelihoods = pp1_encoding('modelFittingWrapper_searchlight',Y,pV,cV); % get model likelihoods
end