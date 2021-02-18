function varargin=pp1_pcmTest(what,varargin);
signal  = 0.3; % signal variance for pattern generation
verbose = 1; % display output for crossvalidated pcm fitting
recipeFile = '/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/git_repos/pcm_toolbox/recipe_finger/data_recipe_finger7T.mat';
% models that you can simulate data from:
% 'muscle'
% 'usage'
% 'eye'
switch what
    case 'chords'
        % returns indicator matrix for chords used in exp.
        % NOTE: this are in a different order compared to cpd1 and cpd2
        % codes.
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
        [Y,partVec,~] = pp1_selectivity('SIM_fingerPatterns','model',model);
        Y = Y{1};
        % housekeeping
        chords    = pp1_selectivity('chords');
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
        [Y,partVec,~] = pp1_selectivity('SIM_fingerPatterns','model',model);
        Y = Y{1};
        % housekeeping
        chords    = pp1_selectivity('chords');
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
                fcn = @(x) pp1_selectivity('SIM_chordPatterns_linearScale','model',x);
            case 'log'
                fcn = @(x) pp1_selectivity('SIM_chordPatterns_logScale','model',x);
        end
        isBad = 1;
        while isBad
            [Y,partVec,condVec] = feval(fcn,model);
            % calc slope of linear regress line b/t activity and chord length
            [~,b] = pp1_selectivity('plot_chordActivity',Y,partVec,0);
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
            [Y{s},partVec{s},condVec{s}] = pp1_selectivity('SIM_fingerPatterns','model',model);
            Y{s} = Y{s}{1};
            % remove run mean
%             C0 = indicatorMatrix('identity',partVec{s});
%             Y{s} = Y{s}-C0*pinv(C0)*Y{s}; 
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3)
        % get models
        M = pp1_selectivity('pcm_getSanityModelsFingers');
        % fit models
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',verbose);
        % plot model G correlations
        %pp1_selectivity('plot_modelCorr',G_hat,G_pred);
        
        % plot fits
        pp1_selectivity('plot_pcmFits',M,Tcv,T,G_pred);
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
            [Y{s},partVec{s},condVec{s}] = pp1_selectivity('SIM_chordPatterns','model',model,'scale',scale);
            % remove run mean
%             C0 = indicatorMatrix('identity',partVec{s});
%             Y{s} = Y{s}-C0*pinv(C0)*Y{s}; 
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3);
        % get models
        M = pp1_selectivity('pcm_getSanityModelsChords');
        % fit models
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',verbose);
        % plot fits
        pp1_selectivity('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));
        % plot avg. activity per number of digits
        subplot(2,numel(M),numel(M)+2);
        pp1_selectivity('plot_chordActivity',Y,partVec); 
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
            [Y{s},partVec{s},condVec{s}] = pp1_selectivity('SIM_chordPatterns','model',model,'scale',scale);
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3);
        % get pcm models
        M{1} = pp1_selectivity('pcm_null');
        M{2} = pp1_selectivity('pcm_noScaleSummation');
        M{3} = pp1_selectivity('pcm_linearSummation');
        M{4} = pp1_selectivity('pcm_nonlinearSummation');
        M{5} = pp1_selectivity('pcm_nonlinearSummationFingers');
        M{6} = pp1_selectivity('pcm_freedirect');
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
        pp1_selectivity('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));
        % plot avg. activity per number of digits
        subplot(2,numel(M),numel(M)+2);
        pp1_selectivity('plot_chordActivity',Y,partVec);
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
            [Y{s},partVec{s},condVec{s}] = pp1_selectivity('SIM_chordPatterns','model',model,'scale',scale);
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3);
        % get pcm models
        M{1} = pp1_selectivity('pcm_null');
        M{2} = pp1_selectivity('pcm_linearSummation');
        M{2}.name = M{2}.type;
        M{4} = pp1_selectivity('pcm_freedirect');
        M{3} = M{2};
        M{3}.Gc = M{3}.Ac*M{3}.Ac';
        M{3} = rmfield(M{3},{'Ac'});
        M{3}.type = 'component';
        M{3}.name = M{3}.type;
        % fit models
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',0);
        % plot fits
        pp1_selectivity('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));
        % plot avg. activity per number of digits
        subplot(2,numel(M),numel(M)+2);
        pp1_selectivity('plot_chordActivity',Y,partVec);
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
            [Y{s},partVec{s},condVec{s}] = pp1_selectivity('SIM_chordPatterns','model',model,'scale',scale);
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        G_hat = mean(G_hat,3);
        modelTimes = zeros(2);
        % get pcm model
        M{1} = pp1_selectivity('pcm_null');
        M{2} = pp1_selectivity('pcm_freedirect');
        % fit model
        tic;
        [Tfd,theta_hat_fd,G_pred_fd] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        modelTimes(1,1) = toc;
        tic;
        [Tcv_fd,theta_cv_fd]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat_fd,'fitScale',1,'verbose',0);
        modelTimes(1,2) = toc;
        % get pcm model
        M{1} = pp1_selectivity('pcm_null');
        M{2} = pp1_selectivity('pcm_freechol');
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
        pp1_selectivity('plot_pcmFits',M,Tcv,T,G_pred);
        title(sprintf('true model: %s',model));
        % plot avg. activity per number of digits
        subplot(2,numel(M),numel(M)+2);
        pp1_selectivity('plot_chordActivity',Y,partVec);
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
        chords = pp1_selectivity('chords');
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
        M.Ac(:,:,1)  = pp1_selectivity('chords');
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
        M{1} = pp1_selectivity('pcm_null');
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
        M{5} = pp1_selectivity('pcm_freedirect');
        varargout = {M};
    case 'pcm_getSanityModelsChords'    
        load('/Users/sarbuckle/Documents/git_repos/pcm_toolbox/recipe_finger/data_recipe_finger7T.mat'); % Model struct
        clear Y condVec partVec
        chords = pp1_selectivity('chords');
        % models
        M{1} = pp1_selectivity('pcm_null');
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
        M{5} = pp1_selectivity('pcm_freedirect');
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
        chords  = pp1_selectivity('chords');
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
