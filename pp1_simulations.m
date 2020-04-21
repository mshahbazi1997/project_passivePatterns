function varargout = pp1_simulations(what,varargin)
signal  = 0.3; % signal variance for pattern generation
verbose = 1; % display output for crossvalidated pcm fitting
recipeFile = '/Users/sarbuckle/MATLAB/pcm_toolbox/recipe_finger/data_recipe_finger7T.mat';
% models that you can simulate data from:
% 'muscle'
% 'usage'
% 'eye'
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
    case 'SIM_fingerPatterns'
        model = 'usage';
        vararginoptions(varargin,{'model'});
        % Simulates patterns where voxels are sparsely tuned to different
        % fingers.
        % Model params are defined in structure M.
        % Simulation stuff is defined in structure D.
        M.signal      = signal;
        M.noise       = 1;
        M.theta       = 0;          % theta parameters to construct model G- is empty because we provide a 'fixed' G
        M.numGparams  = 1;          % number of parameters to define G- we provide a 'fixed' G, so no params needed
        M.type        = 'fixed';    % type of model G- 'fixed' means we provide an unchanging second moment (G)
        % D: model simulation aux. param structure (these fields stay constant)
        D.numSim      = 1;          % number of simulations ("subjects")
        D.numVox      = 250;        % number of simulated voxels
        D.numPart     = 10;         % number of simulated runs/partitions
        D.N           = 1;          % number of simulated trials per partition
        
        switch model
            case 'sparse'
                M.Gc       = eye(5);
                fingers    = find(any(M.Gc));
                numConds   = size(fingers,2);
                numFingers = length(fingers);
                D.numVox   = floor(D.numVox/numFingers);
                % zero-pad output patterns matrix
                for s = 1:D.numSim
                    Y{s} = zeros(numConds*D.numPart,D.numVox*numFingers);
                end
                j = 1; % init. voxel index ticker
                for i = fingers
                    fidx      = [(j*D.numVox - D.numVox+1):((j*D.numVox))]; % indicies of output patterns where generated patterns are to be placed
                    M.Gc      = zeros(numConds);        % set zero tuning for all fingers
                    M.Gc(i,i) = 1;   % turn "on" tuning for specified finger
                    [yy,partVec,condVec] = pcm_generateData(M,M.theta,D,D.numSim,M.signal,M.noise);
                    %[yy{s},partVec,condVec] = makePatterns('G',M.Gc,'nPart',D.numPart,'nVox',D.numVox,'snr',M.signal/M.noise);
                    % harvest patterns for each subject
                    for s = 1:D.numSim 
                        Y{s}(:,fidx) = yy{s};
                    end
                    j = j+1; % update voxel index ticker
                end
            case {'usage','muscle','somat','eye'}
                % Load G from pcm toolbox finger recipe (hand use data)
                load(recipeFile);
                clear Y condVec partVec
                switch model 
                    case 'muscle'
                        %M.Gc = Model(1).G_cent;
                        M.Gc = Model(1).G;
                    case 'usage'
                        %M.Gc = Model(2).G_cent;
                        M.Gc = Model(2).G;
                    case 'somat'
                        error('no support for somat model.');
                        M.Gc = Model(3).G_cent;
                    case 'eye' % identity
                        M.Gc = eye(5);
                end
                % generate patterns for each subject
                for s = 1:D.numSim 
                    [Y,partVec,condVec] = pcm_generateData(M,M.theta,D,D.numSim,M.signal,M.noise);
                    %[Y{s},partVec,condVec] = makePatterns('G',M.Gc,'nPart',D.numPart,'nVox',D.numVox,'snr',M.signal/M.noise);
                end
        end
        % return patterns to user
        varargout= {Y,partVec,condVec};
    case 'SIM_chordPatterns_linearScale'
        model = 'usage';
        vararginoptions(varargin,{'model'});
        % Builds activity patterns for multi-finger chords as a function
        % of multiple single finger patterns
        [Y,partVec,~] = pp1_simulations('SIM_fingerPatterns','model',model);
        Y = Y{1};
        % housekeeping
        chords    = pp1_simulations('chords');
        numChords = size(chords,1);
        runs      = unique(partVec)';
        pV_chords = kron([1:length(runs)]',ones(numChords,1));
        cV_chords = kron(ones(length(runs),1),[1:numChords]');
        % make patterns
        Ychords   = [];
        for r = runs
            y  = Y(partVec==r,:);
            y  = y-10;
            yM = (y' * chords')'; % linearly scale
%            yM = yM-2;
            Ychords = [Ychords; yM]; 
        end
        % imagesc(Ychords(pV_chords==1,:))
        varargout = {Ychords,pV_chords,cV_chords};
    case 'SIM_chordPatterns_logScale'
        model = 'usage';
        vararginoptions(varargin,{'model'});
        scale = 1./[1:0.25:2];%1./log([2.72:2:11]);
        % Builds activity patterns for multi-finger chords as a function
        % of multiple single finger patterns
        [Y,partVec,~] = pp1_simulations('SIM_fingerPatterns','model',model);
        Y = Y{1};
        % housekeeping
        chords    = pp1_simulations('chords');
        numChords = size(chords,1);
        runs      = unique(partVec)';
        pV_chords = kron([1:length(runs)]',ones(numChords,1));
        cV_chords = kron(ones(length(runs),1),[1:numChords]');
        % activity scaling feature
        S          = chords;
        numFingers = sum(S,2);
        for i = 1:4
            S(numFingers==i+1,:) = S(numFingers==i+1,:) .* scale(i);
        end
        % make patterns
        Ychords   = [];
        for r = runs
            y  = Y(partVec==r,:);
            yM = (y' * S')'; % log scale
            Ychords = [Ychords; yM]; 
        end
        % imagesc(Ychords(pV_chords==1,:))
        varargout = {Ychords,pV_chords,cV_chords};
    case 'SIM_chordPatterns'
        % Builds activity patterns for multi-finger chords as a function
        % of multiple single finger patterns
        % Wrapper case to ensure avg. activity scales with positive slope
        % (if scaling is log or linear)
        model = 'usage';
        scale = 'linear';
        vararginoptions(varargin,{'model','scale'});
        switch scale
            case 'linear'
                fcn = @(x) pp1_simulations('SIM_chordPatterns_linearScale','model',x);
            case 'log'
                fcn = @(x) pp1_simulations('SIM_chordPatterns_logScale','model',x);
        end
        isBad = 1;
        while isBad
            [Y,partVec,condVec] = feval(fcn,model);
            % calc slope of linear regress line b/t activity and chord length
            [~,b] = pp1_simulations('plot_chordActivity',Y,partVec,0);
            if b>0 & sum(strcmp(scale,{'log','linear'})) % will still be positive if log scale
                isBad = 0;
            end
        end
        
        varargout = {Y,partVec,condVec};
        
    case 'pcm_sanityModelsFingers'
        runEffect = 'random';
        model     = 'usage';
        numSubj   = 5;
        vararginoptions(varargin,{'model','numSubj'});
        % Generates single finger patterns under muscle, usage, or somat
        % model, fits fixed pcm models to generated data. Sanity check.
        % get data
        for s = 1:numSubj
            [Y{s},partVec{s},condVec{s}] = pp1_simulations('SIM_fingerPatterns','model',model);
            Y{s} = Y{s}{1};
            % remove run mean
%             C0 = indicatorMatrix('identity',partVec{s});
%             Y{s} = Y{s}-C0*pinv(C0)*Y{s}; 
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3)
        % get models
        M = pp1_simulations('pcm_getSanityModelsFingers');
        % fit models
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',verbose);
        % plot model G correlations
        %pp1_simulations('plot_modelCorr',G_hat,G_pred);
        
        % plot fits
        pp1_simulations('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));  
        % plot G_hat
        subplot(2,numel(M),numel(M)+3);
        imagesc(G_hat);
        title('G hat (subj avg.)');
    case 'pcm_sanityModelsChords'
        runEffect = 'random';
        model     = 'usage';
        numSubj   = 5;
        scale     = 'linear';
        vararginoptions(varargin,{'model','numSubj'});
        % Generates single finger patterns under muscle, usage, or somat
        % model, fits fixed pcm models to generated data. Sanity check.
        % get data
        for s = 1:numSubj
            [Y{s},partVec{s},condVec{s}] = pp1_simulations('SIM_chordPatterns','model',model,'scale',scale);
            % remove run mean
%             C0 = indicatorMatrix('identity',partVec{s});
%             Y{s} = Y{s}-C0*pinv(C0)*Y{s}; 
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3);
        % get models
        M = pp1_simulations('pcm_getSanityModelsChords');
        % fit models
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',verbose);
        % plot fits
        pp1_simulations('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));
        % plot avg. activity per number of digits
        subplot(2,numel(M),numel(M)+2);
        pp1_simulations('plot_chordActivity',Y,partVec); 
        % plot G_hat
        subplot(2,numel(M),numel(M)+3);
        imagesc(G_hat);
        title('G hat (subj avg.)');
    case 'pcm_testModels'
        runEffect = 'random';
        model     = 'usage';
        scale     = 'linear';
        numSubj   = 5;
        vararginoptions(varargin,{'model','scale','numSubj'});
        % get data
        for s = 1:numSubj
            [Y{s},partVec{s},condVec{s}] = pp1_simulations('SIM_chordPatterns','model',model,'scale',scale);
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3);
        % get pcm models
        M{1} = pp1_simulations('pcm_null');
        M{2} = pp1_simulations('pcm_noScaleSummation');
        M{3} = pp1_simulations('pcm_linearSummation');
        M{4} = pp1_simulations('pcm_nonlinearSummation');
        M{5} = pp1_simulations('pcm_nonlinearSummationFingers');
        M{6} = pp1_simulations('pcm_freedirect');
        % set starting values for nonlinear scaling model
        scaleParams = log([1.1, 1.2, 1.3, 1.4])';
        M{4}.theta0 = [scaleParams];
        % set starting values for nonlinear finger scaling model (2nd moment for single fingers)
        Fx0 = pcm_free_startingval(G_hat(1:5,1:5));
        scaleParams = log([1.1, 1.2, 1.3, 1.4])';
        M{5}.theta0 = [Fx0;scaleParams];
        % fit models
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',0);
        % plot fits
        pp1_simulations('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));
        % plot avg. activity per number of digits
        subplot(2,numel(M),numel(M)+2);
        pp1_simulations('plot_chordActivity',Y,partVec);
        % plot G_hat
        subplot(2,numel(M),numel(M)+3);
        imagesc(G_hat);
        title('G hat (subj avg.)');
        
        %keyboard
    case 'pcm_testCompVsFeat'
        % sanity to check my assumptions about component vs. feature models
        runEffect = 'random';
        model     = 'eye';
        scale     = 'linear';
        numSubj   = 5;
        vararginoptions(varargin,{'model','scale','numSubj'});
        % get data
        for s = 1:numSubj
            [Y{s},partVec{s},condVec{s}] = pp1_simulations('SIM_chordPatterns','model',model,'scale',scale);
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3);
        % get pcm models
        M{1} = pp1_simulations('pcm_null');
        M{2} = pp1_simulations('pcm_linearSummation');
        M{2}.name = M{2}.type;
        M{4} = pp1_simulations('pcm_freedirect');
        M{3} = M{2};
        M{3}.Gc = M{3}.Ac*M{3}.Ac';
        M{3} = rmfield(M{3},{'Ac'});
        M{3}.type = 'component';
        M{3}.name = M{3}.type;
        % fit models
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',0);
        % plot fits
        pp1_simulations('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));
        % plot avg. activity per number of digits
        subplot(2,numel(M),numel(M)+2);
        pp1_simulations('plot_chordActivity',Y,partVec);
        % plot G_hat
        subplot(2,numel(M),numel(M)+3);
        imagesc(G_hat);
        title('G hat (subj avg.)');
    
    case 'pcm_testNC'
        % test noise ceiling models
        runEffect = 'random';
        model     = 'eye';
        scale     = 'linear';
        numSubj   = 5;
        vararginoptions(varargin,{'model','scale','numSubj'});
        % get data
        for s = 1:numSubj
            [Y{s},partVec{s},condVec{s}] = pp1_simulations('SIM_chordPatterns','model',model,'scale',scale);
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3);
        modelTimes = zeros(2);
        % get pcm model
        M{1} = pp1_simulations('pcm_null');
        M{2} = pp1_simulations('pcm_freedirect');
        % fit model
        tic;
        [Tfd,theta_hat_fd,G_pred_fd] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        modelTimes(1,1) = toc;
        tic;
        [Tcv_fd,theta_cv_fd]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat_fd,'fitScale',1,'verbose',0);
        modelTimes(1,2) = toc;
        % get pcm model
        M{1} = pp1_simulations('pcm_null');
        M{2} = pp1_simulations('pcm_freechol');
        % fit model
        tic;
        [Tfc,theta_hat_fc,G_pred_fc] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        modelTimes(2,1) = toc;
        tic;
        [Tcv_fc,theta_cv_fc]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat_fc,'fitScale',1,'verbose',0);
        modelTimes(2,2) = toc;
        
        keyboard
        
        figure('Color',[1 1 1]);
        subplot()
        
        % plot fits
        pp1_simulations('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));
        % plot avg. activity per number of digits
        subplot(2,numel(M),numel(M)+2);
        pp1_simulations('plot_chordActivity',Y,partVec);
        % plot G_hat
        subplot(2,numel(M),numel(M)+3);
        imagesc(G_hat);
        title('G hat (subj avg.)');
        
        %keyboard
        
        
    case 'pcm_null'
        % Model Null: chords patterns are fully independent
        M.type       = 'component';
        M.numGparams = 1;
        M.name       = 'null';
        M.Gc(:,:,1)  = eye(31);
        varargout = {M};
    case 'pcm_freechol'
        M.type       = 'freechol';  
        M.numCond    = 31;
        M.name       = 'noiseceiling';
        M            = pcm_prepFreeModel(M);
        varargout = {M};
    case 'pcm_freedirect'
        % Naive averaring model- noise ceiling method 1- totall free model
        M.type       = 'freedirect';  
        M.numGparams = 0;
        M.name       = 'noiseceiling';
        
        varargout = {M};
    case 'pcm_noScaleSummation'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns, but no increase in avg. activity w/ chord length
        chords = pp1_simulations('chords');
        chords = chords./(kron(sum(chords,2),ones(1,5)));
        M.type       = 'feature';
        M.numGparams = 1;
        M.name       = 'no scale summation';
        M.Ac(:,:,1)  = chords;
        varargout = {M};
    case 'pcm_linearSummation'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns
        M.type       = 'feature';
        M.numGparams = 1;
        M.name       = 'linear summation';
        M.Ac(:,:,1)  = pp1_simulations('chords');
        varargout = {M};
    case 'pcm_nonlinearSummation'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates ONLY scaling params- no finger params
        M.type       = 'nonlinear';
        M.name       = 'nonlinear summation';
        M.modelpred  = @pp1_modelpred_scale;
        M.numGparams = 4;
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_nonlinearSummationFingers'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates both finger params and scaling params.
        M.type       = 'nonlinear';
        M.name       = 'nonlinear fingers';
        M.modelpred  = @pp1_modelpred_singleFingerScale;
        M.numGparams = 18;
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_getSanityModelsFingers'
        load('/Users/sarbuckle/Documents/git_repos/pcm_toolbox/recipe_finger/data_recipe_finger7T.mat'); % Model struct
        clear Y condVec partVec
        % models
        M{1} = pp1_simulations('pcm_null');
        M{1}.Gc  = eye(5);
        % muscle
        M{2}.type       = 'component';
        M{2}.numGparams = 1;
        M{2}.name       = 'muscle';
        M{2}.Gc(:,:,1)  = Model(1).G;%Model(1).G_cent;
        % usage
        M{3}.type       = 'component';
        M{3}.numGparams = 1;
        M{3}.name       = 'usage';
        M{3}.Gc(:,:,1)  =  Model(2).G;%Model(2).G_cent;
        % somatotopic
        M{4}.type       = 'component';
        M{4}.numGparams = 1;
        M{4}.name       = 'somat';
        M{4}.Gc(:,:,1)  = Model(3).G_cent;
        % noiseceiling
        M{5} = pp1_simulations('pcm_freedirect');
        varargout = {M};
    case 'pcm_getSanityModelsChords'    
        load('/Users/sarbuckle/Documents/git_repos/pcm_toolbox/recipe_finger/data_recipe_finger7T.mat'); % Model struct
        clear Y condVec partVec
        chords = pp1_simulations('chords');
        % models
        M{1} = pp1_simulations('pcm_null');
        M{1}.Gc  = eye(31);
        % muscle
        M{2}.type       = 'fixed';
        M{2}.numGparams = 0;
        M{2}.name       = 'muscle';
        M{2}.Gc(:,:,1)  = chords*Model(1).G_cent*chords';
        % usage
        M{3}.type       = 'fixed';
        M{3}.numGparams = 0;
        M{3}.name       = 'usage';
        M{3}.Gc(:,:,1)  = chords*Model(2).G_cent*chords';
        % somatotopic
        M{4}.type       = 'fixed';
        M{4}.numGparams = 0;
        M{4}.name       = 'somat';
        M{4}.Gc(:,:,1)  = chords*Model(3).G_cent*chords';
        % noiseceiling
        M{5} = pp1_simulations('pcm_freedirect');
        varargout = {M};
    
    case 'plot_pcmFits'
        % assumes last model is noise ceiling
        M      = varargin{1};
        Tcv    = varargin{2};
        T      = varargin{3};
        G_pred = varargin{4};
        
        numModels = numel(M);
        figure('Color',[1 1 1]);
        % plot predicted 2nd moments
        for i = 1:numModels
            subplot(2,numModels,i);
            imagesc(G_pred{i});
            title(M{i}.name)
        end
        % plot scaled likelihoods
        subplot(2,numModels,numModels+1);
        T = pcm_plotModelLikelihood(Tcv,M,'upperceil',T.likelihood(:,end),'style','bar');
        set(gca,'xticklabelrotation',45);
        box off;
    case 'plot_modelCorr'
        % Correlate each subject's second moment with predicted model Gs.
        G_hat     = varargin{1}; % [KxKxS] second moment matrices (one per subj S)
        G_pred    = varargin{2}; % [KxKxM] predicted second moment matrices (one per model M)
        numSubj   = size(G_hat,3);
        numModel  = numel(G_pred);
        % arrange subject second moments into vectors
        Gs = [];
        for i = 1:size(G_hat,3)
            Gs(:,i) = rsa_vectorizeIPM(G_hat(:,:,i));
        end
        % correlate subject second moments to model predicted Gs
        R = [];
        for i = 1:numModel
            modelG  = rsa_vectorizeIPM(G_pred{i});
            r.corr  = corr(Gs',modelG');
            r.sn    = [1:numSubj]';
            r.model = ones(numSubj,1).*i;
            R = addstruct(R,r);
        end
        figure('Color',[1 1 1]);
        plt.bar(R.model,R.corr,'style',style.custom({[0 0 0],[.7 0 0],[0 0 .7],[.9 .6 0],[0 0 0]}),'split',R.model);
        ylim([-1 1]);
        drawline(0,'dir','horz');
        ylabel('pearson''s r');
        xlabel('model no.');
        title('Subject G correlation w/ G pred.')
    case 'plot_chordActivity'    
        % plots average activity per chord length across subjects.
        % Visualize if activity scales linearly/non-linearly
        Y      = varargin{1};
        pV     = varargin{2}; % partition vector
        plotIt = 1;
        if length(varargin)>2
            plotIt = varargin{3}; % plot figure?
        end
        % housekeeping
        chords  = pp1_simulations('chords');
        numD    = sum(chords,2);
        if iscell(Y)
            numSubj = numel(Y);
        else 
            numSubj = 1;
            Y  = {Y};
            pV = {pV};
        end
        D = [];
        b = [];
        % loop through "subjects", calc avg. activity per chord length
        for i = 1:numSubj
            numRun = length(unique(pV{i}));
            d      = []; 
            d.act  = mean(Y{i},2);
            d.numD = kron(ones(numRun,1),numD);
            d      = tapply(d,{'numD'},{'act','mean'});
            d.sn   = ones(5,1).*i;
            D      = addstruct(D,d);
            b      = [b, D.act\D.numD];
        end
        % plot
        if plotIt
            plt.line(D.numD,D.act);
            plt.labels('no. fingers','avg. beta weight (a.u.)');
        end
        varargout = {D,b};


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
        G       = varargin{1};
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5};
        noise   = varargin{6};
        
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation
        
        numCond = size(G,1);
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
    case 'SFT:model_SparseDEPRECIATED'
        % true betas drawn from sparse distribution (gamma), with normally
        % distributed noise across runs
        % Critically, neurons show sparse tuning (i.e. tuned to single
        % conditions)
        numCond = varargin{1};
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5};
        noise   = varargin{6};
        if numel(varargin)==7
            alpha = varargin{7};
        else
            alpha = 1;
        end
        %alpha = signal;
        % quick fix to control for lost signal variance after ensuring
        % sparse activity is positive (noise can be either pos/neg).
        %signal = signal + ((signal/numCond)/2);
        
        % define signal generation
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        %signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation, with abs to make positive sparse
        signalDist = @(x) gammaincinv(x,alpha);
        
        D = [];
        for s = 1:numSubj % per simulated dataset
            % draw uniform values for signal and noise (to be inverted
            % through any arbitrary function later)
            pSignal    = unifrnd(0,1,numCond,numVox/numCond); 
            pNoise     = unifrnd(0,1,numCond*numRun,numVox); 
            
            % Generate true sparse patterns
            U = [];
            for c = 1:numCond
                Gc      = zeros(numCond);
                Gc(c,c) = 1;
                u       = signalDist(pSignal); 
                E       = (u*u'); 
                Z       = abs(E^(-0.5) *u);   % Make random orthonormal vectors, BUT MAKE THEM POSITIVE! (so, mean is no longer zero for signal)
                A       = pcm_diagonalize(Gc); 
                if (size(A,2)>numVox)
                    error('not enough voxels to represent G'); 
                end
                trueU = A*Z(1:size(A,2),:)*sqrt(numVox);  % spread the signal out over all voxels
                trueU = bsxfun(@times,trueU,sqrt(signal));     % scale signal (note this scalar is a little bigger than what user specified to account for lost variance from abs() above)
                U     = [U, trueU];
            end
            % check: rank(cov(U'))==numCond
            X = kron(ones(numRun,1),eye(numCond)); % design matrix
            
            % Now add the random noise 
            Noise  = noiseDist(pNoise); 
            Noise  = bsxfun(@times,Noise,sqrt(noise)); % sacle noise
            d.Y    = X*U + Noise;                % pull through condition specific patterns and add i.i.d. noise
            
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
        
        % get conditions based on sparsity level
        chords = pp1_imana('chords');
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
    
    case 'SFT:model_Sparse_single'
        % true patterns are sparse (0's and 1's), with normally
        % distributed noise across runs
        numCond = varargin{1};
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5}; % signal variance
        noise   = varargin{6}; % noise variance
        
        % define signal generation
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        signal = signal*(numCond); % rescale the signal to account for the lack of signal variance aross non-responsive voxels
        D = [];
        for s = 1:numSubj % per simulated dataset
            % draw uniform values for signal and noise (to be inverted
            % through any arbitrary function later)
            pNoise     = unifrnd(0,1,numCond*numRun,numVox); 
            
            % Generate true sparse patterns
            U = [];
            Z = zeros(numCond,numVox/numCond);
            for c = 1:numCond
                u      = Z;
                u(c,:) = 1;
                trueU  = bsxfun(@times,u,sqrt(signal));     % scale signal (note this scalar is a little bigger than what user specified to account for lost variance from abs() above)
                U      = [U, trueU];
            end
            % check: rank(cov(U'))==numCond
            X = kron(ones(numRun,1),eye(numCond)); % design matrix
            
            % Now add the random noise 
            Noise  = noiseDist(pNoise); 
            Noise  = bsxfun(@times,Noise,sqrt(noise)); % sacle noise
            d.Y    = X*U + Noise;                % pull through condition specific patterns and add i.i.d. noise
            
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
        N = pp1_imana('SFT:model_Gnoise',eye(numCond),100*numCond,1,1,svar,evar);
        S = pp1_imana('SFT:model_Sparse',1,100*numCond,1,1,svar,evar);
        U = pp1_imana('SFT:model_Uniform',eye(numCond),100*numCond,1,1,svar,evar);
        
        % plot
        sty1 = style.custom({'darkgray'});
        sty2 = style.custom({'red'});
        sty3 = style.custom({'orange'});
        subplot(2,3,1); plt.hist(N.Y(:),'style',sty1); title('mvnrnd'); xlabel('activity'); ylabel('count');
        subplot(2,3,2); plt.hist(S.Y(:),'style',sty2); title('sparse');
        subplot(2,3,3); plt.hist(U.Y(:),'style',sty3); title('uniform');
        
        subplot(2,3,4); plt.scatter(N.Y(1,:)',N.Y(2,:)','style',sty1,'regression','off'); title('mvnrnd'); xlabel('cond 1 act.'); ylabel('cond 2 act.');
        subplot(2,3,5); plt.scatter(S.Y(1,:)',S.Y(2,:)','style',sty2,'regression','off'); title('sparse');
        subplot(2,3,6); plt.scatter(U.Y(1,:)',U.Y(2,:)','style',sty3,'regression','off'); title('uniform');
        
    case 'estErrVarOLD'
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
        
        % % NOTE: the signal variance estimate appears to underestimate the true signal for sparse data
        % this is probably due to how I enforce sparsity in these patterns: abs(data). I need to think about this? 
        
        Y = varargin{1};    % patterns   [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        
        take = tril(bsxfun(@eq,c,c'),-1);
        G    = cov(Y');
        R    = corr(Y');
        
        var_S = sum(G(take))/sum(sum(take));
        r     = sum(R(take))/sum(sum(take));
        
        var_E = var_S/r - var_S; % since all-possible pairs CV, can drop the last term in equation: (g/r - g)*(m-1), where m==# of partitions folded into training dataset
        varargout = {var_E,var_S};
    case 'estErrVar'
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
        
        % % NOTE: the signal variance estimate appears to underestimate the true signal for sparse data
        % this is probably due to how I enforce sparsity in these patterns: abs(data). I need to think about this? 
        
        Y = varargin{1};    % patterns   [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        r = varargin{3};    % partitions [regressors x 1]
        
        nV = size(Y,2);
        nP = numel(unique(r));
        nC = numel(unique(c));
        Ya = zeros(nP,nV*nC);
        
        for pp = 1:nP
            y = Y(r==pp,:);
            Ya(pp,:) = y(:)';
        end
        
        take = logical(tril(ones(nP),-1));
        G    = cov(Ya');
        R    = corr(Ya');
        
        var_S = sum(G(take))/sum(sum(take));
        r     = sum(R(take))/sum(sum(take));
        
        var_E = var_S/r - var_S; 
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
    
        
    case 'testErrVarV1'
        % case to generate data at different levels of signal-to-noise,
        % plots results for evaluation
        numSim = 100; % 100 simulated datasets per model type
        evar   = [0.01,0.1,1,10];
        numCond = 5;
        G      = eye(numCond);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        svar = 1; % signal variance
        for ii = evar
            d_sparse = pp1_imana('SFT:model_Sparse',1,numVox,numRun,numSim,svar,ii);
            d_norm   = pp1_imana('SFT:model_Gnoise',G,numVox,numRun,numSim,svar,ii);
            %d_uni    = pp1_imana('SFT:model_Uniform',G,numVox,numRun,numSim,svar,ii);
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
                [var_e,var_s] = pp1_imana('SFT:estimateVariances',dn.Y,dn.cond,dn.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dnn           = tapply(dn,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_imana('SFT:estimateSFT',dnn.Y));
                d = addstruct(d,q);

                % calc snr for sparse data
                q.model       = 2;
                C0 = indicatorMatrix('identity',ds.run);
                ds.Y = ds.Y - C0*pinv(C0)*ds.Y;
                [var_e,var_s] = pp1_imana('SFT:estimateVariances',ds.Y,ds.cond,ds.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dss           = tapply(ds,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_imana('SFT:estimateSFT',dss.Y));
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
        
        % plot ERROR VARIANCE ESTIAMTES (varies)
        subplot(1,4,2);
        plt.line(D.evarLog,log(D.evar_Est),'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error) - log scale');
        xlabel('specified var(error) - log scale');
        title('altered error (v.1)');
        axis square
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);

        % plot SIGNAL VARIANCE ESTIMATES (unchanged)
        subplot(1,4,1);
        plt.line(D.evarLog,D.svar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(signal)');
        xlabel('specified var(error) - log scale');
        title('unchanged signal (v.1)');
        axis square
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        
        % plot SNRs:
        subplot(1,4,3);
        plt.scatter(D.snr_true,D.snr_est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('snr estimate');
        xlabel('snr true');
        title('unchanged signal (v.1)');
        axis square
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        %set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        
        % plot SFTs
        subplot(1,4,4);
        plt.line(D.evarLog,D.sftEst,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('sft');
        xlabel('specified var(error) - log scale');
        title('selectivity (v.1)');
        axis square
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
                
        varargout = {D};
    case 'testSigVarV1'
        % case to generate data at different levels of signal-to-noise,
        % plots results for evaluation
        numSim = 100; % 100 simulated datasets per model type
        svar   = [0.01,0.1,1,10];
        numCond = 5;
        G      = eye(numCond); %diag([1,2,3,4,5])./mean(1:5);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        evar = 1; % error variance
        for ii = svar
            d_sparse = pp1_imana('SFT:model_Sparse',1,numVox,numRun,numSim,ii,evar);
            d_norm   = pp1_imana('SFT:model_Gnoise',G,numVox,numRun,numSim,ii,evar);
            %d_uni    = pp1_imana('SFT:model_Uniform',G,numVox,numRun,numSim,ii,evar);
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
                [var_e,var_s] = pp1_imana('SFT:estimateVariances',dn.Y,dn.cond,dn.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dnn           = tapply(dn,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_imana('SFT:estimateSFT',dnn.Y));
                d = addstruct(d,q);

                % calc snr for sparse data
                q.model       = 2;
                C0 = indicatorMatrix('identity',ds.run);
                ds.Y = ds.Y - C0*pinv(C0)*ds.Y;
                [var_e,var_s] = pp1_imana('SFT:estimateVariances',ds.Y,ds.cond,ds.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dss           = tapply(ds,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_imana('SFT:estimateSFT',dss.Y));
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
        plt.line(D.svarLog,log(D.svar_Est),'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(signal) - log scale');
        xlabel('specified var(error)');
        title('altered signal (v.1)');
        axis square
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);

        % plot ERROR VARIANCE ESTIMATES (unchanged)
        subplot(1,4,2);
        plt.line(D.svarLog,D.evar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error)');
        xlabel('specified var(signal) - log scale');
        title('unchanged error (v.1)');
        axis square
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);
        
        % plot SNRs:
        subplot(1,4,3);
        plt.scatter(D.snr_true,D.snr_est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('snr estimate');
        xlabel('snr true');
        title('unchanged signal (v.1)');
        axis square
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';
        %set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        
        % plot SFTs
        subplot(1,4,4);
        plt.line(D.svarLog,D.sftEst,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('sft');
        xlabel('specified var(signal) - log scale');
        title('selectivity (v.1)');
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
            d_sparse = pp1_imana('SFT:model_Sparse_single',numCond,numVox,numRun,numSim,svar,ii);
            d_norm   = pp1_imana('SFT:model_Gnoise',G,numVox,numRun,numSim,svar,ii);
            d_uni    = pp1_imana('SFT:model_Uniform',G,numVox,numRun,numSim,svar,ii);
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
                [var_e,var_s] = pp1_imana('estErrVarMeans',dn.Y,dn.cond,dn.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dnn           = tapply(dn,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_imana('estSingleFingerTuning',dnn.Y));
                d = addstruct(d,q);

                % calc snr for sparse data
                q.model       = 2;
                [var_e,var_s] = pp1_imana('estErrVarMeans',ds.Y,ds.cond,ds.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dss           = tapply(ds,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_imana('estSingleFingerTuning',dss.Y));
                d = addstruct(d,q);

%                 % calc snr for uniform data
%                 q.model       = 3;
%                 [var_e,var_s] = pp1_imana('estErrVarMeans',du.Y,du.cond,du.run);
%                 q.evar_Est    = var_e;
%                 q.svar_Est    = var_s;
%                 duu           = tapply(du,{'cond'},{'Y','mean'});
%                 q.sftEst      = mean(pp1_imana('estSingleFingerTuning',duu.Y));
                
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
            %d_sparse = pp1_imana('SFT:model_Sparse_single',numCond,numVox,numRun,numSim,ii,evar); % this version does not allow chaning of sparsity
            d_sparse = pp1_imana('SFT:model_Sparse',sparsity,numVox,numRun,numSim,ii,evar); % allows changing sparsity levels
            d_norm   = pp1_imana('SFT:model_Gnoise',G,numVox,numRun,numSim,ii,evar);
            %d_uni    = pp1_imana('SFT:model_Uniform',G,numVox,numRun,numSim,ii,evar);
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
                [var_e,var_s] = pp1_imana('estErrVarMeans',dn.Y,dn.cond,dn.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dnn           = tapply(dn,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_imana('estSingleFingerTuning',dnn.Y));
                d = addstruct(d,q);

                % calc snr for sparse data
                q.model       = 2;
                [var_e,var_s] = pp1_imana('estErrVarMeans',ds.Y,ds.cond,ds.run);
                q.evar_Est    = var_e;
                q.svar_Est    = var_s;
                dss           = tapply(ds,{'cond'},{'Y','mean'});
                q.sftEst      = mean(pp1_imana('estSingleFingerTuning',dss.Y));
                d = addstruct(d,q);

%                 % calc snr for uniform data
%                 q.model       = 3;
%                 [var_e,var_s] = pp1_imana('estErrVarMeans',du.Y,du.cond,du.run);
%                 q.evar_Est    = var_e;
%                 q.svar_Est    = var_s;
%                 duu           = tapply(du,{'cond'},{'Y','mean'});
%                 q.sftEst      = mean(pp1_imana('estSingleFingerTuning',duu.Y));
                
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
    case 'testSigVarSparse'
        % case to generate data at different levels of signal-to-noise,
        % plots results for evaluation
        numSim = 10; % 100 simulated datasets per model type
        svar   = [0.1:0.1:1,2];
        alpha  = 1;%[0.1:0.1:1,2];
        numCond = 5;
        G      = eye(numCond);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        evar = 1; % signal variance
        for ii = svar
            for a=1:numel(alpha)
                d_sparse = pp1_imana('SFT:model_Sparse_single',numCond,numVox,numRun,numSim,ii,evar,alpha(a));
                d = [];
                for jj = 1:numSim
                    ds = getrow(d_sparse,d_sparse.sn==jj);
                    % indexing fields to output structure
                    q.sn   = jj;
                    q.evar = evar;
                    q.svar = ii;

                    q.model       = a;
                    [var_e,var_s] = pp1_imana('estErrVarMeans',ds.Y,ds.cond,ds.run);
                    q.evar_Est    = var_e;
                    q.svar_Est    = var_s;
                    dss           = tapply(ds,{'cond'},{'Y','mean'});
                    q.sftEst      = mean(pp1_imana('estSingleFingerTuning',dss.Y));
                    d = addstruct(d,q);

                    d = addstruct(d,q);
                end
                D = addstruct(D,d);
            end
        end
                
        sty = style.custom(plt.helper.get_shades(numel(alpha),'cool','descend'));
        
        % plot estimate of the error (varies)
        subplot(1,3,2);
        plt.line(D.svar,D.evar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error)');
        xlabel('specified var(signal)');
        title('unchaged error (v.2)');
        axis square
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);

        % plot estimate of the signal (unchanged)
        subplot(1,3,1);
        plt.line(D.svar,D.svar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(signal)');
        xlabel('specified var(signal)');
        title('altered signal (v.2)');
        axis square
        set(gca,'xticklabel',string(svar),'xticklabelrotation',45);
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';

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
    case 'testErrVarSparse'
        % case to generate data at different levels of signal-to-noise,
        % plots results for evaluation
        numSim = 10; % 100 simulated datasets per model type
        alpha  = [0.1:0.1:1,2];
        evar   = [0,0.1,1,10:10:100];
        numCond = 5;
        G      = eye(numCond);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        svar = 1; % signal variance
        for ii = evar
            for a=1:numel(alpha)
                d_sparse = pp1_imana('SFT:model_Sparse_single',numCond,numVox,numRun,numSim,svar,ii,alpha(a));
                d = [];
                for jj = 1:numSim
                    ds = getrow(d_sparse,d_sparse.sn==jj);
                    % indexing fields to output structure
                    q.sn   = jj;
                    q.evar = ii;
                    q.svar = svar;

                    q.model       = a;
                    [var_e,var_s] = pp1_imana('estErrVarMeans',ds.Y,ds.cond,ds.run);
                    q.evar_Est    = var_e;
                    q.svar_Est    = var_s;
                    dss           = tapply(ds,{'cond'},{'Y','mean'});
                    q.sftEst      = mean(pp1_imana('estSingleFingerTuning',dss.Y));
                    d = addstruct(d,q);

                    d = addstruct(d,q);
                end
                D = addstruct(D,d);
            end
        end
                
        sty = style.custom(plt.helper.get_shades(numel(alpha),'cool','descend'));
        
        % plot estimate of the error (varies)
        subplot(1,3,2);
        plt.line(D.evar,D.evar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error)');
        xlabel('specified var(error)');
        title('altered error (v.2)');
        axis square
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';

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
    case 'testAlphaSig'
        % case to test indepedence of alpha and signal variance parameters
        numSim = 10; % 100 simulated datasets per model type
        alpha  = [0.1:0.1:1,2];
        evar   = [0,0.1,1,10:10:100];
        numCond = 5;
        G      = eye(numCond);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        svar = 1; % signal variance
        for ii = evar
            for a=1:numel(alpha)
                d_sparse = pp1_imana('SFT:model_SparseDEPRECIATED',numCond,numVox,numRun,numSim,alpha(a),ii,alpha(a));
                d = [];
                for jj = 1:numSim
                    ds = getrow(d_sparse,d_sparse.sn==jj);
                    % indexing fields to output structure
                    q.sn   = jj;
                    q.evar = ii;
                    q.svar = alpha(a);

                    q.model       = alpha(a);
                    [var_e,var_s] = pp1_imana('estErrVarMeans',ds.Y,ds.cond,ds.run);
                    q.evar_Est    = var_e;
                    q.svar_Est    = var_s;
                    dss           = tapply(ds,{'cond'},{'Y','mean'});
                    q.sftEst      = mean(pp1_imana('estSingleFingerTuning',dss.Y));
                    d = addstruct(d,q);

                    d = addstruct(d,q);
                end
                D = addstruct(D,d);
            end
        end
                
        sty = style.custom(plt.helper.get_shades(numel(alpha),'cool','descend'));
        
        % plot estimate of the error (varies)
        subplot(1,3,2);
        plt.line(D.evar,D.evar_Est,'split',D.model,'style',sty);
        plt.legend({'mvn','sparse'});
        ylabel('estimated var(error)');
        xlabel('specified var(error)');
        title('altered error (v.2)');
        axis square
        set(gca,'xticklabel',string(evar),'xticklabelrotation',45);
        % identity line
        h = refline(1,0);
        h.Color = 'k';
        h.LineWidth = 2;
        h.LineStyle = ':';

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
    case 'testFstatGauss'
        % tests that restricting sft analysis to voxels with significant
        % responses does not introduce any bias
        numSim = [];
        svar = [];
        evar = [];
        vararginoptions(varargin,{'svar','evar','numSim'});
        numCond = 5;
        G      = eye(numCond);
        %G = diag([1,2,3,4,5]);
        numVox = 100*numCond;
        numRun = 11;
        D      = [];
        
        for ee=evar
            for ii=1:numel(svar)
                ss = svar(ii);
                % simulate mvnrnd data (non-sparse)
                d_norm = pp1_imana('SFT:model_Gnoise',G,numVox,numRun,numSim,ss,ee);
                for jj = 1:numSim
                    % indexing fields to output structure
                    q.sparse = 0;
                    q.sn   = jj;
                    q.evar = ee;
                    q.svar = ss;
                    q.svarNum = ii;
                    % get data for this simulation
                    dn = getrow(d_norm,d_norm.sn==jj);
                    % get f stat per voxel
                    [F,Fcrit] = pp1_imana('SFT:calcFstat',dn.Y,dn.cond,dn.run);
                    %[Fcv,FcritCV] = pp1_imana('calcFstatCV',dn.Y,dn.cond,dn.run);
                    %[Z,Fj,p] = spmj_zFtest_imcalc(dn.Y,dn.cond,dn.run); % check correct F estiamtes
                    % restrict further analyses to significant voxels
                    dn.Y = dn.Y(:,F>=Fcrit); % using non-cv f stat as thresholder
                    q.numVoxOrig = numVox;
                    q.numVoxSig  = sum(F>=Fcrit);
                    if q.numVoxSig==0
                        q.evar_Est = nan;
                        q.svar_Est = nan;
                        q.sftEst   = nan;
                        D = addstruct(D,q);
                        continue
                    end
                    % remove run mean patterns of siginificant voxels only
                    C0 = indicatorMatrix('identity',dn.run);
                    dn.Y = dn.Y - C0*pinv(C0)*dn.Y;
                    % calc snr for significant voxels
                    [var_e,var_s] = pp1_imana('SFT:estimateVariances',dn.Y,dn.cond,dn.run);
                    q.evar_Est    = var_e;
                    q.svar_Est    = var_s;
                    % estimate sft of significant voxels
                    dnn           = tapply(dn,{'cond'},{'Y','mean'});
                    q.sftEst      = mean(pp1_imana('SFT:estimateSFT',dnn.Y));
                    
                    D = addstruct(D,q);
                end
            end
        end
        
        % apply log scale (slightly
        % tricky for log(0) so that's why the code below exists)
        D.svarLog = D.svar;
        D.svarLog(D.svar>0) = log(D.svar(D.svar>0));
        D.svarLog(D.svar==0)= -sum(diff(log(svar(svar>0))));
        
        varargout = {D};
    case 'testFstatSparse'
        % tests that restricting sft analysis to voxels with significant
        % responses does not introduce any bias
        numSim = [];
        svar = [];
        evar = [];
        vararginoptions(varargin,{'svar','evar','numSim'});
        sparsity = 1;
        numCond = 5;
        numVox = 100*numCond;
        numRun = 11;
        D      = [];

        for ee=evar
            for ii=1:numel(svar)
                ss = svar(ii);
                % simulate mvnrnd data (non-sparse)
                d_sparse = pp1_imana('SFT:model_Sparse',sparsity,numVox,numRun,numSim,ss,ee);
                for jj = 1:numSim
                    % indexing fields to output structure
                    q.sparse = 1;
                    q.sn   = jj;
                    q.evar = ee;
                    q.svar = ss;
                    q.svarNum = ii;
                    % get data for this simulation
                    ds = getrow(d_sparse,d_sparse.sn==jj);
                    % get f stat per voxel
                    [F,Fcrit] = pp1_imana('SFT:calcFstat',ds.Y,ds.cond,ds.run);
                    %[Fcv,FcritCV] = pp1_imana('calcFstatCV',dn.Y,dn.cond,dn.run);
                    %[Z,Fj,p] = spmj_zFtest_imcalc(dn.Y,dn.cond,dn.run); % check correct F estiamtes
                    % restrict further analyses to significant voxels
                    ds.Y = ds.Y(:,F>=Fcrit); % using non-cv f stat as thresholder
                    q.numVoxOrig = numVox;
                    q.numVoxSig  = sum(F>=Fcrit);
                    if q.numVoxSig==0
                        q.evar_Est = nan;
                        q.svar_Est = nan;
                        q.sftEst   = nan;
                        D = addstruct(D,q);
                        continue
                    end
                    % remove run mean patterns of siginificant voxels only
                    C0 = indicatorMatrix('identity',ds.run);
                    ds.Y = ds.Y - C0*pinv(C0)*ds.Y;
                    % calc snr for significant voxels
                    [var_e,var_s] = pp1_imana('SFT:estimateVariances',ds.Y,ds.cond,ds.run);
                    q.evar_Est    = var_e;
                    q.svar_Est    = var_s;
                    % estimate sft of significant voxels
                    dnn           = tapply(ds,{'cond'},{'Y','mean'});
                    q.sftEst      = mean(pp1_imana('SFT:estimateSFT',dnn.Y));
                    
                    D = addstruct(D,q);
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
        numSim = 20;
        svar = [0,0.001,0.01,0.1,1,10];
        evar = 1;
        D1 = pp1_simulations('testFstatGauss','svar',svar,'evar',evar,'numSim',numSim);  % simulate data that has gaussian tuning, & apply F-thresholding
        D2 = pp1_simulations('testFstatSparse','svar',svar,'evar',evar,'numSim',numSim); % simulate data with sparse tuning (with exact same signal & error variances) , & apply F-thresholding
        D = addstruct(D1,D2);
        D = getrow(D,D.evar==1);

        % plot
        labels = string(svar);
        sty = style.custom({'darkgray','red'});
        % plot % voxels retained in analyses
        subplot(1,2,1);
        plt.bar(D.svarLog,(D.numVoxSig./D.numVoxOrig).*100,'split',D.sparse,'style',sty);
        xlabel('signal-error variance ratio (logscale)');
        ylabel('% voxels');
        title('% significant voxels (omnibus F-test)');
        legend off
        setXTicksGroups(gca,numel(svar),2,labels);
        % plot estimated sfts
        subplot(1,2,2)
        plt.line(D.svarLog,D.sftEst,'split',D.sparse,'style',sty,'errorfcn','stderr');
        ylabel('selectivity');
        xlabel('signal / error variance (logscale)');
        title('estimated selectivity')
        set(gca,'xticklabel',labels);
        legend off
        varargout = {D};
    

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
