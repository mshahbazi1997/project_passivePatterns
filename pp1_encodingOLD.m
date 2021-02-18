function varargout=pp1_encodingOLD(what,varargin)
%% details

% saarbuckle 2020

% encoding simulations for passivePatterns1 fmri experiment

% Directories for data and saving simulation results
%baseDir = '/Volumes/MotorControl/data/passivePatterns/passivePatterns1/fmri';
baseDir = '/Users/sarbuckle/DATA/passivePatterns1/fmri/'; % where we save model selection accuracies to
dataDir = fullfile(baseDir,'RegionOfInterest'); % where betas are saved (we use single finger betas from here for model simulations)

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
        chords    = pp1_simulations('chords');
        numDigits = sum(chords,2);
        numConds  = 31;%sum(numDigits>1); % # of multi-fnger conditions we are simulating data for
        numVox    = size(Usf,2);
        % Make true patterns according to perscribed model:
        X = kron(ones(numRun,1),eye(numConds)); % multi-run design
        if snr~=0 % simulated data has signal component
            Umf = pp1_encoding('pred_MFPatterns',Usf,model,modelTheta);
            noiseVar = var(Umf(:)) / snr; % use target snr to calc noise variance for this dataset
        elseif snr==0 % simulated data lacks signal component; patterns are purely noise
            Umf   = zeros(numConds,numVox);
            noiseVar = 1;
        end
        U_run = X*Umf; % pull true multi-finger Us through to each run
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
        G  = (Umf*Umf')./numVox; 
        
        varargout = {U_sim,cV,pV,Umf,G};
        % U_hat - model patterns with i.i.d. noise [numConds*numRun x numVox]
        % cV    - condition vector [numConds*numRun x 1]
        % pV    - partition vector [numConds*numRun x 1]
        % U     - true (noiseless) patterns [numConds x numVox]
        % G     - true (noiseless) second moment of multi finger patterns [numConds x numConds]
    case 'pred_MFPatterns'
        % factorization of encoding models:
        Ysf = varargin{1};
        model = varargin{2};
        theta = varargin{3}; % theta(1) = baseline, theta(2:end) = model params
        % get multi-finger design:
        chords = pp1_encoding('chords');
        X = chords; 
        switch model
            case 'summation'
                Ymf_hat = X*(Ysf-theta(1)) + theta(1);
            case 'summationNL_OLD'
                % y = x^theta
                Ymf_hat = pp1_encoding('pred_MFPatterns',Ysf,'summation',theta(1)); % do summation
                % estimate theta
                %[theta,err,gEst] = fminsearch(@(x) fitPower(x,Ymf_hat(:),Ymf(:)),0.5);
                % squish patterns with theta param:
                sign_Uhat= sign(Ymf_hat);
                abs_Uhat = abs(Ymf_hat);
                Ymf_hat  = sign_Uhat.*(abs_Uhat.^theta(2));
            case 'summationNL_tanh'
                % y = tanh(x)*theta(2)
                % theta 2 is a scalar param to scale the patterns so they are no longer decimals, then squishes, then undo scaling
                Ymf_hat = pp1_encoding('pred_MFPatterns',Ysf,'summation',theta(1)); % do linear summation
                % nonlinearize patterns:
                Ymf_hat = tanh(Ymf_hat).*theta(2);
        end
        varargout = {Ymf_hat};
    case 'calcModelFit'
        % case to calculate r2 and r between predicted and true patterns.
        Ymf_pred = varargin{1};
        Ymf_true = varargin{2};
        % R2:
        res = Ymf_true(:)-Ymf_pred(:);
        rss = sum(res.^2);
        tss = sum(Ymf_true(:).*Ymf_true(:)); % total sums of squares
        r2  = 1-(rss/tss);
        % R:
        ssC = sum(Ymf_true(:).*Ymf_pred(:)); % covariance of predicted with true patterns
        pss = sum(Ymf_pred(:).*Ymf_pred(:)); % predicted sums of squares
        r   = ssC/sqrt(tss.*pss);
        varargout = {r2,r};
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
    
    case 'fit'
        % wrapper to test different model fitting approaches and assess
        % model selection accuracies using different metrics
        
        % metrics:
        %   1 : r2 (encoding CV)
        %   2 : r (encoding CV)
        %   3 : r2 (pcm individ cv)
        %   4 : r (pcm individ cv)
        %   5 : likelihood (pcm individ cv)
        
        % sim params:
        snr = [0,0.5,1];%[0,0.01,0.025,0.05,0.1,0.25,0.5,1];
        sn  = 9;%[2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        numSim          = 50; % sims per subject per snr lvl under each model (2 models = 50 simulations * # snr lvls)
        theta_baseline  = 0;    % true baseline
        theta_power     = 0.75; % squishing factor (exponent)
        theta_scale     = 100;   % scaling factor (applied before squishing, then divided out)
        
        % define params for the linear and nonlinear models:
        modelThetas = {[theta_baseline], [theta_baseline theta_power theta_scale]}; % {1}=linear, {2}=nonlinear
        
        % get data for this roi, glm, sn(s):
        [Ysubj,partVec,condVec,Gcv] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        
        
        % do simulations at different snrs for EACH subject:
        A = []; % overall output structure
        for ii=1:numel(sn)
            s=sn(ii);
            As = []; % output structure for this exp. subj's sims
            % estimate true single finger Us using PCM:
            [U,Gu] = pp1_encoding('pcm_estTrueU',Ysubj{ii},partVec{ii},condVec{ii});
            Usf    = U(1:5,:); clear U; % single finger pattern estimates
            % get true (noiseless) pcm model Gs:
            [~,~,~,~,G_lin] = pp1_encoding('sim_MFPatterns',Usf,'summation',modelThetas{1},'numSim',1,'snr',inf); % inf=no noise
            [~,~,~,~,G_nl]  = pp1_encoding('sim_MFPatterns',Usf,'summationNL',modelThetas{2},'numSim',1,'snr',inf);
            M = makePCM_M(eye(31),G_lin,G_nl);
            
            for jj=1:numel(snr)
                nn=snr(jj);
                fprintf('s%02d snr: %2.2f...',s,nn);% verbose output to user
                % simulate patterns under each model using the true single finger Us : 
                [U_hat_linear,cV,pV] = pp1_encoding('sim_MFPatterns',Usf,'summation',modelThetas{1},'numSim',numSim,'snr',nn); % make linear summation patterns
                U_hat_nonlin         = pp1_encoding('sim_MFPatterns',Usf,'summationNL',modelThetas{2},'numSim',numSim,'snr',nn); % make nonlinear scaled patterns
                Ysim      = [U_hat_linear,U_hat_nonlin]; % concatenate simulated datasets into one cell array
                trueModel = [ones(1,numSim), ones(1,numSim).*2]; % create vector denoting the true model of each dataset in Ymf (1=lin, 2=nl)
                
                % do encoding model fitting, and assess selection accuracy:
                % model selection using crossvalidated encoding (leave-one-out crossval scheme):
                fprintf('CV encoding...');
                [A1,D_encoding] = pp1_encoding('doEncoding',Ysim,cV,pV,trueModel,modelThetas);
                A1.metric = [1;2];
                A1.approach = [1;1];
                % model selection using PCM (individual leave-one-out crossval scheme):
                fprintf('PCM...');
                [A2,D_pcm] = pp1_encoding('doPCM',Ysim,cV,pV,trueModel,modelThetas);
                A2.metric = [3;4;5];
                A2.approach = [2;2;2];
                % collapse results into one struct, and only keep avg.
                % evaluation metric and the metric accuracy:
                An = addstruct(A1,A2);
                vA = ones(size(An.approach));
                An.snr = vA.*nn;
                An.snrGroup = vA.*jj;
                As = addstruct(As,An);
                fprintf('done\n');
            end
            vS = ones(size(As.snr));
            As.sn  = vS.*sn(ii);
            As.glm = vS.*glm;
            As.roi = vS.*roi;
            A = addstruct(A,As);
        end
        % determine if easy or hard model comparison (for filename):
        if truePower==1
            difficulty = 'identical';
        elseif truePower<=0.4 % arbitrary cutoff
            difficulty = 'easy';
        else 
            difficulty = 'hard';
        end
        outName = fullfile(baseDir,sprintf('pp1_modelAcc_roi%d_glm%d_%s.mat',roi,glm,difficulty));
        save(outName,'A','D_encoding','D_pcm');
    
        
    case 'do'
        % wrapper to test different model fitting approaches and assess
        % model selection accuracies using different metrics
        
        % metrics:
        %   1 : r2 (encoding CV)
        %   2 : r (encoding CV)
        %   3 : r2 (pcm individ cv)
        %   4 : r (pcm individ cv)
        %   5 : likelihood (pcm individ cv)
        
        % sim params:
        snr = [0,0.5,1];%[0,0.01,0.025,0.05,0.1,0.25,0.5,1];
        sn  = 9;%[2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        numSim = 50; % sims per subject per snr lvl under each model (2 models = 50 simulations * # snr lvls)
        theta_baseline  = 0; % don't shift true baseline
        theta_scale     = 1.01;
        
        % define params for the linear and nonlinear models:
        modelThetas = {[theta_baseline], [theta_baseline theta_scale]}; % {1}=linear, {2}=nonlinear
        
        % get data:
        [Y,partVec,condVec] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        
        %error('under revamp')
        % prep PCM model fitting struct:
        M{1}.type         = 'fixed';
        M{1}.numGparams   = 0;
        M{1}.name         = 'summation';
        M{1}.fitAlgorithm = 'NR';
        M{2}.type         = 'fixed';
        M{2}.numGparams   = 0;
        M{2}.name         = 'summationNL_tanh';
        M{2}.fitAlgorithm = 'NR';
        
        % do simulations at different snrs for EACH subject:
        A = []; % overall output structure
        for ii=1:numel(sn)
            % estimate true single finger Usf in crossval fashion using PCM:
            U = pp1_encoding('pcm_estTrueU',Y{ii},partVec{ii},condVec{ii});
            
            
            s=sn(ii);
            As = []; % output structure for this exp. subj's sims
            
            % get true (noiseless) model G for each model:
            [~,~,~,~,G_lin] = pp1_encoding('sim_MFPatterns',Usf{ii},'summation',modelThetas{1},'numSim',1,'snr',inf); % inf=no noise
            [~,~,~,~,G_nl]  = pp1_encoding('sim_MFPatterns',Usf{ii},'summationNL_tanh',modelThetas{2},'numSim',1,'snr',inf);
            M{1}.Gc = G_lin; % multi-finger G for PCM models
            M{2}.Gc = G_nl; 
            
            for jj=1:numel(snr)
                nn=snr(jj);
                fprintf('s%02d snr: %2.2f...',s,nn);% verbose output to user
                % make linear summation patterns:
                [U_hat_linear,cV,pV] = pp1_encoding('sim_MFPatterns',Usf{ii},'summation',modelThetas{1},'numSim',numSim,'snr',nn);
                % make nonlinear scaled patterns:
                [U_hat_nonlin,~,~] = pp1_encoding('sim_MFPatterns',Usf{ii},'summationNL',modelThetas{2},'numSim',numSim,'snr',nn);
                % concatenate simulated datasets into one cell array:
                Ymf = [U_hat_linear,U_hat_nonlin]; 
                % create vector denoting the true model of each dataset in Ymf (1=lin, 2=nl):
                trueModel = [ones(numSim,1); ones(numSim,1).*2];
                % model selection using crossvalidated encoding (leave-one-out crossval scheme):
                fprintf('CV encoding...');
                A2 = pp1_encoding('doEncoding',Ymf,Usf{ii},cV,pV,trueModel,modelThetas,'leaveOneOut');
                A2.metric = [1;2];
                A2.approach = [1;1];
                % model selection using PCM (individual leave-one-out crossval scheme):
                fprintf('PCM...');
                    % drop any single-finger patterns for this fitting 
                    % (in previous encoding mdoels, we needed these to 
                    % predict U_hats, but never included them in fit calculations)
                mfIdx = cV>5; % condition # 6-31 are multi-finger chords
                Y = cellfun(@(y) y(mfIdx,:),Y,'uni',0); 
                cV = cV(mfIdx); pV = pV(mfIdx);
                A3 = pp1_encoding('doPCM',Y,cV,pV,trueModel,M);
                A3.metric = [3;4;5];
                A3.approach = [2;2;2];
                % collapse results into one struct, and only keep avg.
                % evaluation metric and the metric accuracy:
                An = addstruct(A1,A2);
                An = addstruct(An,A3);
                vA = ones(size(An.approach));
                An.snr = vA.*nn;
                An.snrGroup = vA.*jj;
                As = addstruct(As,An);
                fprintf('done\n');
            end
            vS = ones(size(As.snr));
            As.sn  = vS.*sn(ii);
            As.glm = vS.*glm;
            As.roi = vS.*roi;
            A = addstruct(A,As);
        end
        % determine if easy or hard model comparison (for filename):
        if truePower==1
            difficulty = 'identical';
        elseif truePower<=0.4 % arbitrary cutoff
            difficulty = 'easy';
        else 
            difficulty = 'hard';
        end
        outName = fullfile(baseDir,sprintf('pp1_modelAcc_roi%d_glm%d_%s.mat',roi,glm,difficulty));
        save(outName,'-struct','A');
    case 'doPCM_OLD'
        % inputs
        Y  = varargin{1}; % cell array of datasets, [31 x numVox]
        cV = varargin{2}; % condition vector
        pV = varargin{3}; % partition vector (assume cV and pV are same across subjs)
        trueModel = varargin{4}; % vector denoting what is true model for each dataset in Y
        M = varargin{5}; % pcm model fitting structure

        [Tcv,DD,theta_cv] = pcm_fitModelIndividCrossval(Y,M,pV,cV,...
                'evaluation',{'likelihood','R2','R'},...
                'crossvalScheme','leaveOneOut',...
                'runEffect','random','verbose',0);
        % eval to decide best model:
        % assign best model (if tied, randomly choose a model from tied models):
        T = [];
        numDatasets = numel(Y);
        Tcv.trueModel = trueModel';
        Tcv.winner_r2 = zeros(numDatasets,1);
        Tcv.winner_r  = Tcv.winner_r2;
        Tcv.winner_like = Tcv.winner_r2;
        for ii=1:numDatasets
            % determine winning models:
            evalCriterion = [Tcv.R2(ii,:);Tcv.R(ii,:);Tcv.likelihood(ii,:)];
            winners = pp1_encoding('chooseWinner',evalCriterion);
            Tcv.winner_r2(ii,1)   = winners(1);
            Tcv.winner_r(ii,1)    = winners(2);
            Tcv.winner_like(ii,1) = winners(3);
        end

        % calculate model selection accuracy across datasets:
        r2_correct = Tcv.winner_r2==Tcv.trueModel; % how many times was true model the best fitting?
        r_correct  = Tcv.winner_r2==Tcv.trueModel;
        like_correct  = Tcv.winner_like==Tcv.trueModel;

        Tcv_r2 = getrow(Tcv,r2_correct);
        Tcv_r  = getrow(Tcv,r_correct);
        Tcv_like  = getrow(Tcv,like_correct);

        r2_trueIdx = logical(pcm_indicatorMatrix('identity',Tcv_r2.trueModel));
        r_trueIdx = logical(pcm_indicatorMatrix('identity',Tcv_r.trueModel));
        like_trueIdx = logical(pcm_indicatorMatrix('identity',Tcv_like.trueModel));

        A.metricName = {'R2';'R';'LIKE'};
        A.metric     = [1:3]';
        A.value      = [mean(Tcv_r2.R2(r2_trueIdx)); mean(Tcv_r.R(r_trueIdx)); mean(abs(Tcv_like.likelihood(:,1) - Tcv_like.likelihood(:,2)))];
        A.acc        = [sum(r2_correct); sum(r_correct); sum(like_correct)]./numDatasets;
        
        
        % crossval model fitting:
%         didIt = false;
%         try
%             [Tcv,DD,theta_cv] = pcm_fitModelIndividCrossval(Y,M,pV,cV,...
%                 'evaluation',{'likelihood','R2','R','Rpool'},...
%                 'crossvalScheme','leaveOneOut',...
%                 'runEffect','random','verbose',0);
%             didIt=true;
%         catch ME
%             % if there was an error, likely because matrix was singular.
%             % I'm not sure how to get around this other than skipping over
%             % this problem for now.
%             fprintf('skipping (error): %s...',ME.message);
%         end
%         % eval to decide best model:
%         % assign best model (if tied, randomly choose a model from tied models):
%         T = [];
%         if didIt
%             numDatasets = numel(Y);
%             Tcv.trueModel = trueModel;
%             Tcv.winner_r2 = zeros(numDatasets,1);
%             Tcv.winner_r  = Tcv.winner_r2;
%             Tcv.winner_like = Tcv.winner_r2;
%             for ii=1:numDatasets
%                 % determine winning models:
%                 evalCriterion = [Tcv.R2(ii,:);Tcv.R(ii,:);Tcv.likelihood(ii,:)];
%                 winners = pp1_encoding('chooseWinner',evalCriterion);
%                 Tcv.winner_r2(ii,1)   = winners(1);
%                 Tcv.winner_r(ii,1)    = winners(2);
%                 Tcv.winner_like(ii,1) = winners(3);
%             end
%             
%             % calculate model selection accuracy across datasets:
%             r2_correct = Tcv.winner_r2==Tcv.trueModel; % how many times was true model the best fitting?
%             r_correct  = Tcv.winner_r2==Tcv.trueModel;
%             like_correct  = Tcv.winner_like==Tcv.trueModel;
%             
%             Tcv_r2 = getrow(Tcv,r2_correct);
%             Tcv_r  = getrow(Tcv,r_correct);
%             Tcv_like  = getrow(Tcv,like_correct);
%             
%             r2_trueIdx = logical(pcm_indicatorMatrix('identity',Tcv_r2.trueModel));
%             r_trueIdx = logical(pcm_indicatorMatrix('identity',Tcv_r.trueModel));
%             like_trueIdx = logical(pcm_indicatorMatrix('identity',Tcv_like.trueModel));
%             
%             A.metricName = {'R2';'R';'LIKE'};
%             A.metric     = [1:3]';
%             A.value      = [mean(Tcv_r2.R2(r2_trueIdx)); mean(Tcv_r.R(r_trueIdx)); mean(abs(Tcv_like.likelihood(:,1) - Tcv_like.likelihood(:,2)))];
%             A.acc        = [sum(r2_correct); sum(r_correct); sum(like_correct)]./numDatasets;
%             
%         else
%             % if there was an error, skip this for this data batch.
%             % NOTE: this rarely (if never) happens, but I don't want to
%             % lose a full simulation set because there was an error on 1/1000 datasets.
%             A.metricName = {'R2';'R','LIKE'};
%             A.metric     = [1:3]';
%             A.value      = [nan; nan; nan];
%             A.acc        = [nan; nan; nan];
%             
%         end
        varargout = {A,Tcv,DD};
    case 'doEncoding_OLD'
        % inputs
        Y  = varargin{1}; % cell array of datasets, [31 x numVox]
        cV = varargin{2}; % condition vector
        pV = varargin{3}; % partition vector (assume cV and pV are same across subjs)
        trueModel = varargin{4}; % vector denoting what is true model for each dataset in Y
        modelTheta = varargin{5};
        crossvalScheme = varargin{6}; % 'none' or 'leaveOneOut'
        
        % housekeeping
        parts     = unique(pV);
        numPart   = numel(parts);
        model     = {'summation','summationNL'};
        numModel  = numel(model);
        numDatasets = numel(Y);
        chords    = pp1_encoding('chords');
        numDigits = sum(chords,2);
        sfIdx     = numDigits==1;
        mfIdx     = numDigits>1;
        D         = []; % output struct

        % define crossvalidation folds:
        partI   = {};
        switch crossvalScheme
            case 'none' % no crossvalidation
                partI={parts};
            case 'leaveOneOut'
                for p=1:numPart
                    partI{p}=parts(p);
                end
        end
        % fit 
        for ii=1:numDatasets % for each dataset
            d.sim = ii;
            d.trueModel = trueModel(ii);
            % loop through crossval folds and test models:
            numFolds = numel(partI);
            dp=[];
            for p=1:numFolds
                % training data is avg. single finger patterns across all
                % but one run
                % test data is multi-finger patterns from left-out run
                trainIdx = ~ismember(pV,partI{p}); 
                testIdx  = ismember(pV,partI{p}); 
                if strcmp(crossvalScheme,'none')
                    trainIdx = testIdx;
                end
                % make design matrices for training and test splits
                Xtrain = pcm_indicatorMatrix('identity',cV(trainIdx));
                Xtest  = pcm_indicatorMatrix('identity',cV(testIdx));
                % get data for training and test, avg. across conditions
                % within each split:
                Ytrain = pinv(Xtrain) * Y{ii}(trainIdx,:);
                Ytest  = pinv(Xtest) * Y{ii}(testIdx,:);
                % split into single finger and multi finger data
                Ysf_train = Ytrain(sfIdx,:); % single finger training patterns
                Ymf_test  = Ytest(mfIdx,:);  % multi finger patterns are test data
                % loop through models and predict mf chords accordingly:
                for mm=1:numModel % for each model
                    d.testModel(1,mm) = mm;
                    % predict multi finger chords according to model type:
                    Ymf_hat = pp1_encoding('pred_MFPatterns',Ysf_train,model{mm},modelTheta{mm});
                    % evaluate model fits:
                    [d.r2(1,mm), d.r(1,mm)] = pp1_encoding('calcModelFit',Ymf_hat,Ymf_test);
                end
                d.fold  = p;
                dp=addstruct(dp,d);
            end
            % avg. fits across CV folds:
            dp = tapply(dp,{'sim','trueModel'},{'r2','mean(x,1)'},{'r','mean(x,1)'},{'testModel','mean(x,1)'});
            % determine winning models:
            evalCriterion = [dp.r2; dp.r];
            winners = pp1_encoding('chooseWinner',evalCriterion);
            dp.winner_r2 = winners(1);
            dp.winner_r  = winners(2);
            D=addstruct(D,dp);
        end
        % calculate model selection accuracy across datasets:
        D.r2_correct = D.winner_r2==D.trueModel; % how many times was true model the best fitting?
        D.r_correct  = D.winner_r==D.trueModel;
        
        D_r2 = getrow(D,D.r2_correct);
        D_r  = getrow(D,D.r_correct);
        r2_trueIdx = logical(pcm_indicatorMatrix('identity',D_r2.trueModel));
        r_trueIdx = logical(pcm_indicatorMatrix('identity',D_r.trueModel));
        
        A.metricName = {'R2';'R'};
        A.metric     = [1:2]';
        A.value      = [mean(D_r2.r2(r2_trueIdx)); mean(D_r.r(r_trueIdx))];
        A.acc        = [sum(D.r2_correct); sum(D.r_correct)]./numDatasets;
        
        varargout = {A,D};
    case 'doEncoding'
        % case to do crossvalidated encoding models
        
        % inputs
        Y  = varargin{1}; % cell array of datasets, {[31*numRun x numVox]}
        cV = varargin{2}; % condition vector
        pV = varargin{3}; % partition vector (assume cV and pV are same across Y datasets)
        trueModel  = varargin{4}; % vector denoting what is true model for each dataset in Y
        modelTheta = varargin{5};
        
        % housekeeping
        model     = {'summation','summationNL'};
        numModel  = numel(model);
        numDatasets = numel(Y);
        D         = []; % output struct

        % define leave-one-out crossvalidation folds:
        parts   = unique(pV);
        numPart = numel(parts);
        partI   = {};
        for p=1:numPart
            partI{p}=parts(p);
        end
        % fit 
        for ii=1:numDatasets % for each dataset
            d.sim = ii;
            d.trueModel = trueModel(ii);
            % loop through crossval folds and test models:
            numFolds = numel(partI);
            dp=[];
            for p=1:numFolds
                % training data is avg. single finger patterns across all
                % but one run
                % test data is multi-finger patterns from left-out run
                trainIdx = ~ismember(pV,partI{p}); 
                testIdx  = ismember(pV,partI{p}); 
                % use training data to estimate true Us:
                Y_train  = Y{ii}(trainIdx,:);
                pV_train = pV(trainIdx);
                cV_train = cV(trainIdx);
                U_train  = pp1_encoding('pcm_estTrueU',Y_train,pV_train,cV_train);
                Usf_train  = U_train(1:5,:); clear U_train; % single finger pattern estimates
                % get patterns to test in this fold:
                Xtest  = pcm_indicatorMatrix('identity',cV(testIdx));
                Ytest  = pinv(Xtest) * Y{ii}(testIdx,:);
                % loop through models and predict mf chords accordingly:
                for mm=1:numModel % for each model
                    d.testModel(1,mm) = mm;
                    % predict multi finger chords according to model type:
                    Ypred = pp1_encoding('pred_MFPatterns',Usf_train,model{mm},modelTheta{mm});
                    % evaluate model fits:
                    [d.r2(1,mm), d.r(1,mm)] = pp1_encoding('calcModelFit',Ypred,Ytest);
                end
                d.fold  = p;
                dp=addstruct(dp,d);
            end
            % avg. fits across CV folds:
            dp = tapply(dp,{'sim','trueModel'},{'r2','mean(x,1)'},{'r','mean(x,1)'},{'testModel','mean(x,1)'});
            % determine winning models:
            evalCriterion = [dp.r2; dp.r];
            winners       = pp1_encoding('chooseWinner',evalCriterion);
            dp.winner_r2  = winners(1);
            dp.winner_r   = winners(2);
            D=addstruct(D,dp);
        end
        % calculate model selection accuracy across datasets:
        D.r2_correct = D.winner_r2==D.trueModel; % how many times was true model the best fitting?
        D.r_correct  = D.winner_r==D.trueModel;
        
        D_r2_good = getrow(D,D.r2_correct);
        D_r_good  = getrow(D,D.r_correct);
        D_r2_bad  = getrow(D,~D.r2_correct);
        D_r_bad   = getrow(D,~D.r_correct);
        r2_goodIdx = logical(pcm_indicatorMatrix('identity',D_r2_good.trueModel));
        r2_badIdx  = logical(pcm_indicatorMatrix('identity',D_r2_bad.trueModel));
        r_goodIdx  = logical(pcm_indicatorMatrix('identity',D_r_good.trueModel));
        r_badIdx   = logical(pcm_indicatorMatrix('identity',D_r_bad.trueModel));
        
        A.metricName = {'R2';'R'};
        A.metric     = [1:2]';
        A.value_correct = [mean(D_r2_good.r2(r2_goodIdx)); mean(D_r_good.r(r_goodIdx))];
        A.value_incorrect = [mean(D_r2_bad.r2(r2_badIdx)); mean(D_r_bad.r(r_badIdx))];
        A.acc        = [sum(D.r2_correct); sum(D.r_correct)]./numDatasets;
        
        varargout = {A,D};
    case 'doPCM'
        % inputs
        Y  = varargin{1}; % cell array of datasets, [31 x numVox]
        cV = varargin{2}; % condition vector
        pV = varargin{3}; % partition vector (assume cV and pV are same across subjs)
        trueModel = varargin{4}; % vector denoting what is true model for each dataset in Y
        modelTheta = varargin{5}; % pcm model fitting structure
        model = {'summation','summationNL'};
        
        numDatasets = numel(Y);
        numModel    = numel(model);
        
        % define leave-one-out crossvalidation folds:
        parts   = unique(pV);
        numPart = numel(parts);
        partI   = {};
        for pp=1:numPart
            partI{pp}=parts(pp);
        end
        D = []; % output struct
        for ii=1:numDatasets
            d.sim = ii;
            d.trueModel = trueModel(ii);
            % loop through crossval folds and test models:
            numFolds = numel(partI);
            dp=[];
            Gc = {};
            for pp=1:numFolds
                % training data is avg. single finger patterns across all
                % but one run
                % test data is multi-finger patterns from left-out run
                trainIdx = ~ismember(pV,partI{pp}); 
                testIdx  = ismember(pV,partI{pp}); 
                % use training data to estimate true Us:
                Y_train  = Y{ii}(trainIdx,:);
                pV_train = pV(trainIdx);
                cV_train = cV(trainIdx);
                U_train  = pp1_encoding('pcm_estTrueU',Y_train,pV_train,cV_train);
                Usf_train  = U_train(1:5,:); clear U_train; % single finger pattern estimates
                % get patterns to test in this fold:
                Xtest  = pcm_indicatorMatrix('identity',cV(testIdx));
                Ytest  = pinv(Xtest) * Y{ii}(testIdx,:);
                % loop through models and predict mf chords accordingly:
                for mm=1:numModel % for each model
                    % predict multi finger chords according to model type:
                    Ypred = pp1_encoding('pred_MFPatterns',Usf_train,model{mm},modelTheta{mm});
                    Gc{mm,pp} = Ypred * Ypred' ./ size(Ypred,2);
                end
                d.fold  = pp;
                dp=addstruct(dp,d);
            end
            % now build model structure to fit these models (we do this in
            % a funky way b/c we want to fit in crossval fashion:
            for mm=1:numModel 
                for pp=1:numFolds
                    idx = pp + numFolds * (mm-1);
                    % define PCM model for predicted patterns:
                    M{idx}.name = sprintf('model %d run %d',mm,pp);
                    M{idx}.type = 'fixed';
                    M{idx}.numGparams = 0;
                    M{idx}.theta0 = [];
                    M{idx}.fitAlgorithm = 'minimize';
                    M{idx}.Gc = Gc{mm,pp};
                end
            end
            % use left-out data to assess model Gs:
            % note: we allow a scalar here to account for different snr
            % values across runs
            [tcv_avg,tcv_fold,theta_hat] = pcm_fitModelIndividCrossval(Y(ii),M,pV,cV,'runEffect','random','verbose',0,'crossvalScheme','leaveOneOut');
            idx = logical(repmat(eye(numFolds),1,2));
            dp.like = reshape(tcv_fold.likelihood(idx),[numFolds,numModel]);
            dp.tss  = repmat(tcv_fold.TSS(:,1),1,2);
            dp.rss  = reshape(tcv_fold.RSS(idx),[numFolds,numModel]);
            dp.ssc  = reshape(tcv_fold.SSC(idx),[numFolds,numModel]);
            dp.ss1  = dp.tss; %reshape(tcv_fold.SS1(idx),[numFolds,numModel]);
            dp.ss2  = reshape(tcv_fold.SS2(idx),[numFolds,numModel]);
            % calc criterion per run
            dp.r2 = 1-(dp.rss./dp.tss);
            dp.r  = dp.ssc./sqrt(dp.ss1.*dp.ss2);
                
            % avg. fits across CV folds:
            dp = tapply(dp,{'sim','trueModel'},{'r2','mean(x,1)'},{'r','mean(x,1)'},{'like','sum'});
            % determine winning models:
            evalCriterion = [dp.r2; dp.r; dp.like];
            winners       = pp1_encoding('chooseWinner',evalCriterion);
            dp.winner_r2  = winners(1);
            dp.winner_r   = winners(2);
            dp.winner_like= winners(3);
            D=addstruct(D,dp);
        end

        % calculate model selection accuracy across datasets:
        D.r2_correct = D.winner_r2==D.trueModel; % how many times was true model the best fitting?
        D.r_correct  = D.winner_r2==D.trueModel;
        D.like_correct  = D.winner_like==D.trueModel;

        D_r2_good = getrow(D,D.r2_correct);
        D_r_good  = getrow(D,D.r_correct);
        D_r2_bad  = getrow(D,~D.r2_correct);
        D_r_bad   = getrow(D,~D.r_correct);
        D_like_good = getrow(D,D.like_correct);
        D_like_bad = getrow(D,~D.like_correct);
        
        r2_goodIdx = logical(pcm_indicatorMatrix('identity',D_r2_good.trueModel));
        r2_badIdx  = logical(pcm_indicatorMatrix('identity',D_r2_bad.trueModel));
        r_goodIdx  = logical(pcm_indicatorMatrix('identity',D_r_good.trueModel));
        r_badIdx   = logical(pcm_indicatorMatrix('identity',D_r_bad.trueModel));
        like_goodIdx = logical(pcm_indicatorMatrix('identity',D_like_good.trueModel));
        like_badIdx = logical(pcm_indicatorMatrix('identity',D_like_bad.trueModel));
        
        A.metricName = {'R2';'R';'Like'};
        A.metric     = [3:5]';
        A.value_correct = [mean(D_r2_good.r2(r2_goodIdx)); mean(D_r_good.r(r_goodIdx)); mean(D_like_good.like(like_goodIdx))];
        A.value_incorrect = [mean(D_r2_bad.r2(r2_badIdx)); mean(D_r_bad.r(r_badIdx)); mean(D_like_bad.like(like_badIdx))];
        A.acc        = [sum(D.r2_correct); sum(D.r_correct); sum(D.like_correct)]./numDatasets;
        
        varargout = {A,D};
    
        
    case 'plot'
        % plot accuracies for easy, hard, and same model comparisons
        roi = 2; % hardcoded
        glm = 4;
        % load simulation model selection accuracies
        Aeasy = load(fullfile(baseDir,sprintf('pp1_modelAcc_roi%d_glm%d_easy.mat',roi,glm))); % models are quite different from each other
        Ahard = load(fullfile(baseDir,sprintf('pp1_modelAcc_roi%d_glm%d_hard.mat',roi,glm))); % models are quite similar (small nl scaling effect)
        Asame = load(fullfile(baseDir,sprintf('pp1_modelAcc_roi%d_glm%d_identical.mat',roi,glm))); % models are identical (scaling param==1; as a comprehensive test)
        Aeasy.difficult = zeros(size(Aeasy.snr));
        Ahard.difficult = ones(size(Ahard.snr));
        Asame.difficult = ones(size(Asame.snr)).*2;
        % make plotting groups so scaling along x-axis is log10 (ugly but oh well):
        Aeasy.snrLog10 = log10(Aeasy.snr);
        Ahard.snrLog10 = log10(Ahard.snr);
        Aeasy.snrLog10(isinf(Aeasy.snrLog10)) = -2.3;
        Ahard.snrLog10(isinf(Ahard.snrLog10)) = -2.3;
        Asame.snrLog10 = nan(size(Asame.snr));
        A = addstruct(Aeasy,Ahard);
        A = addstruct(A,Asame);
        % metrics:
        % 1 : r2 (encoding CV)
        % 2 : r (encoding CV)
        % 3 : r2 (pcm)
        % 4 : r (pcm)
        % 5 : likelihood (pcm)
        
        % plot styling
        lgd = {'R2 (encod cv)','R (encod cv)',...
            'R2 (pcm)','R (pcm)','Like (pcm)'};
        %sty = style.custom({[0.8 0 0],[0 0.8 0],[0.8 0 0],[0 0.8 0],[0 0 0.8],[0.8 0 0],[0 0.8 0],[0 0 0.8],[0 0 0]});
        %sty.general.linestyle = {'-','-',':',':',':','--','--','--','--'};
        CAT.linestyle = {'-','-','-.','-.','-.'};
        CAT.linecolor = {[0.8 0 0],[0 0 0.8],[0.8 0 0],[0 0 0.8],[0.3 0.3 0.3]};
        CAT.markerfill = 'none';%CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        CAT.shadecolor = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        CAT.errorcap = 0;
        CAT.markersize = 8;
        CAT.markertype = {'s','^','s','^','o'};
        CAT.linewidth = {[1],[1],[1],[1],[1]};
        %PLOT
        titles= {'model selection accuracy (EASY)','model selection accuracy (HARD)','model selection accuracy (SAME)'};
        for ii=0:2
            snr = unique(A.snr(A.difficult==ii));
            subplot(2,2,ii+1);
            if ii<2
                lineplot(A.snrLog10,A.acc,'split',A.metric,'CAT',CAT,'plotfcn','nanmean','errorbars','plusminus','subset',A.difficult==ii);
                legend off
                ylim([0.4 1]);
            else
                lineplot(A.snr,A.acc,'split',A.metric,'CAT',CAT,'plotfcn','nanmean','errorbars','plusminus','subset',A.difficult==ii);
                leg = legend(lgd);
                leg.FontSize = 14;
                legend boxoff
                ylim([0 1]);
            end
            ylabel('proportion correct');
            xlabel('sig / noise');
            title(titles{ii+1});
            drawline(0.5,'dir','horz','linestyle','-');
            set(gca,'xticklabelrotation',45,'FontSize',12,'xticklabel',snr);
        end
       
        varargout = {A};
    case 'patchSNR'
        % overlay region specific snr patches on current figure axis:
        D = varargin{1}; % datastructure of region-specific snr values
        scale = varargin{2}; % scale (if not raw values) of the x-axis
        clrs = varargin{3}; % cell of colors per patch
        labels = varargin{4}; % region name(s) (cellstring)
        roi = unique(D.roi);
        for ii=1:numel(roi)
            d = getrow(D,D.roi==roi(ii));
            lowB = mean(d.snr)-stderr(d.snr);
            uppB = mean(d.snr)+stderr(d.snr);
            switch scale
                case 'raw'
                case 'log'
                    lowB = log(lowB);
                    uppB = log(uppB);
                case 'log10'
                    lowB = log10(lowB);
                    uppB = log10(uppB);
            end
            % create patch:
            ylims = ylim;
            X = [lowB uppB uppB lowB];
            Y = kron(ylims,ones(1,2));
            patch(X,Y,clrs{ii},'FaceAlpha',0.4,'EdgeAlpha',0);
            text(lowB,ylims(1)+0.025,labels{ii});
        end
        
    case 'compareModels'
        % compares mean activities and variances from actual patterns (X) to model
        % predicted patterns (Y)
        % sim params:
        snr = inf;
        sn  = 9; 
        roi = 2; % data from Ba 3b
        glm = 4;
        numSim = 1; % sims per subject per snr lvl under each model
        
        % model params:
        theta_baseline  = -0.025; % don't shift true baseline
        theta_power     = 0.5; % scaling factor (increases values <1, squishes values >1)
        theta_scale     = 50;
        % define params for the linear and nonlinear models:
        modelThetas = {[theta_baseline], [theta_baseline theta_power theta_scale]}; % {1}=linear, {2}=nonlinear
        
        scale = 0; % t/f scale mean activity and variances by scalar (estimated using model and true Gs)
        % easy = 0.25, hard = 0.75, same=1
        vararginoptions(varargin,{'truePower','snr','scale'});
        chords = pp1_encoding('chords');
        numDigits = sum(chords,2);
        
        %mfIdx     = numDigits>0; % which conditions to include when plotting?
        % get data:
        [Y,partVec,condVec,Gcv] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % estimate true single finger Usf using PCM:
        U = pp1_encoding('pcm_estTrueU',Y{1},partVec,condVec);
        Usf = cellfun(@(x) x(1:5,:), U, 'uni',0); % get single-finger chords
        
        % do simulations at different snrs for EACH subject:
        D = [];
        v = ones(31,1);
        G_hat_lin = nan(26,26,numel(sn));
        G_hat_nl  = nan(26,26,numel(sn));
        for ii=1:numel(sn)
            Usf_subj = Usf{ii}; % true single finger patterns for model sims
            Gemp = Usf_subj*Usf_subj' ./ size(Usf_subj,2); % empirical estimated G
            Gcv  = Gcv(:,:,ii);
            % simulate data under each model:
            [U_hat_lin,cV1,pV1] = pp1_encoding('sim_MFPatterns',Usf_subj,'summation',modelThetas{1},'numSim',numSim,'snr',snr); % inf=no noise
            [U_hat_nl,cV2,pV2]  = pp1_encoding('sim_MFPatterns',Usf_subj,'summationNL',modelThetas{2},'numSim',numSim,'snr',snr);
            % avg. betas across runs for simulated datasets:
            C0pinv = pinv(pcm_indicatorMatrix('identity',cV1));
            U_avg_lin = mean(cell2mat(cellfun(@(x) mean(C0pinv*x,2),U_hat_lin,'uni',0)),2); % avg. across runs and simulations
            C0pinv = pinv(pcm_indicatorMatrix('identity',cV2));
            U_avg_nl  = mean(cell2mat(cellfun(@(x) mean(C0pinv*x,2),U_hat_nl,'uni',0)),2);
%             U_avg_lin = U_avg_lin(mfIdx); % drop single finger conds
%             U_avg_nl  = U_avg_nl(mfIdx); 
            % estimate second moments from simulated noisy data & avg. across simulations (messy but works fast):
            g_lin_subj = rsa_squareIPM(mean(cell2mat(cellfun(@(x) rsa_vectorizeIPM(pcm_estGCrossval(x,pV1,cV1))', U_hat_lin, 'uni',0)),2)');
            g_nl_subj  = rsa_squareIPM(mean(cell2mat(cellfun(@(x) rsa_vectorizeIPM(pcm_estGCrossval(x,pV2,cV2))', U_hat_nl, 'uni',0)),2)');
             G_hat_lin(:,:,ii) = g_lin_subj; 
             G_hat_nl(:,:,ii)  = g_nl_subj;
            % scale the variances?:
            if scale
                % scale variances and mean activities (latter by
                % sqrt(variance scalar))
                g_lin_subj= G_hat_lin(:,:,ii);
                g_nl_subj = G_hat_nl(:,:,ii);
                lin_scale = pinv(g_lin_subj(:))*Gsubj(:);
                nl_scale  = pinv(g_nl_subj(:))*Gsubj(:);
                G_hat_lin(:,:,ii) = G_hat_lin(:,:,ii) .* lin_scale;
                G_hat_nl(:,:,ii)  = G_hat_nl(:,:,ii) .* nl_scale;
                U_avg_lin = U_avg_lin .* sqrt(lin_scale);
                U_avg_nl  = U_avg_nl .* sqrt(nl_scale);
            end
            % save to output struct
            d=[];
            d.sn        = v.*sn(ii);
            d.glm       = v.*glm;
            d.roi       = v.*roi;
            d.snr       = v.*snr;
            d.numDigits = numDigits;
            d.var       = diag(Gcv);
            d.act       = mean(U{ii},2);
            d.model     = v; % true=1
            D=addstruct(D,d);
            d.var       = [diag(Gemp(1:5,1:5)); diag(G_hat_lin(:,:,ii))];
            d.act       = [mean(Usf_subj,2); U_avg_lin];
            d.model     = v.*2; % linear=2
            D=addstruct(D,d);
            d.var       = [diag(Gemp(1:5,1:5)); diag(G_hat_nl(:,:,ii))];
            d.act       = [mean(Usf_subj,2); U_avg_nl];
            d.model     = v.*3; % nonlinear=3
            D=addstruct(D,d);
        end
        % avg across numDigits:
        D=tapply(D,{'sn','roi','glm','snr','model','numDigits'},{'act','mean'},{'var','mean'});
        % plot styling
        lgd = {'true','linear','nonlinear'};
        CAT.linecolor = {[0 0 0],[0.8 0 0],[0 0 0.8]};
        CAT.errorcolor = CAT.linecolor;
        CAT.errorcap = 0;
        CAT.markersize = 8;
        CAT.markertype = {'s','s','s'};
        CAT.markerfill = CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        % PLOT ACTIVITIES
        subplot(2,3,1);
        lineplot(D.numDigits,D.act,'split',D.model,'CAT',CAT,'errorbars','plusminus');
        drawline(theta_baseline,'linestyle',':','dir','horz');
        title(sprintf('MEAN ACTIVITY\nbaseline: %1.3f\npower: %1.3f',theta_baseline,theta_power));
        ylabel('avg beta');
        xlabel('# digits');
        legend off;
        % PLOT VARIANCES (diag G)
        subplot(2,3,2);
        lineplot(D.numDigits,D.var,'split',D.model,'CAT',CAT,'errorbars','plusminus');
        if scale
            title(sprintf('MEAN VARIANCE\nscaled variance'));
        else
            title(sprintf('MEAN VARIANCE\nunscaled variance'));
        end
        ylabel('avg variance');
        xlabel('# digits');
        leg = legend(lgd);
        leg.FontSize = 14;
        legend boxoff
        % PLOT MEAN vs. VARIANCE (diag G)
        subplot(2,3,3);
        xyplot(D.act,D.var,D.numDigits,'split',D.model,'CAT',CAT,'errorbars','plusminus');
        title('activity vs. variance');
        ylabel('mean variance');
        xlabel('mean activity');
        legend off;
        % PLOT MODEL Gs:
        for ii=1:3
            subplot(2,3,3+ii);
            switch ii
                case 1 % empirical G
                    imagesc(mean(Gcv,3));
                    title('crossval G (avg.)');
                    xlabel('chord #');
                    ylabel('chord #');
                case 2 % linear summation model
                    imagesc(blockdiag(nan(5),mean(G_hat_lin,3)));
                    title('linear summation G (avg.)');
                case 3 % summation nonlinear scaling
                    imagesc(blockdiag(nan(5),mean(G_hat_nl,3)));
                    title(sprintf('summation & nonlinear scaling G (avg.)'));
            end
            % draw #digit lines: 
            drawline(5.5,'dir','horz');  drawline(5.5,'dir','vert');
            drawline(15.5,'dir','horz'); drawline(15.5,'dir','vert');
            drawline(25.5,'dir','horz'); drawline(25.5,'dir','vert');
            drawline(30.5,'dir','horz'); drawline(30.5,'dir','vert');
        end
        
        varargout = {D};
end
end

function M = makePCM_M(varargin)
numG = numel(varargin);
for ii=1:numG
    M{ii}.type         = 'fixed';
    M{ii}.numGparams   = 0;
    M{ii}.name         = num2str(ii);
    M{ii}.fitAlgorithm = 'NR';
    M{ii}.Gc = varargin{ii};
end
end
