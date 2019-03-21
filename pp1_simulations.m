function varargout = pp1_simulations(what,varargin)

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
    case 'SIM_sparseVoxelTuning'
        % Simulates patterns where voxels are sparsely tuned to different
        % fingers.
        % Model params are defined in structure M.
        % Simulation stuff is defined in structure D.
        M.signal      = 0.1;
        M.noise       = 1;
        M.theta       = [];         % theta parameters to construct model G- is empty because we provide a 'fixed' G
        M.numGparams  = 0;          % number of parameters to define G- we provide a 'fixed' G, so no params needed
        M.type        = 'fixed';    % type of model G- 'fixed' means we provide an unchanging second moment (G)
        M.Gc          = eye(5);
        % D: model simulation aux. param structure (these fields stay constant)
        D.numSim      = 1;          % number of simulations ("subjects")
        D.numVox      = 250;        % number of simulated voxels
        D.numPart     = 10;         % number of simulated runs/partitions
        D.N           = 1;          % number of simulated trials per partition
        vararginoptions(varargin,{'M','D'});
        %signalDist = @(x) norminv(x,0,1); % default signal generation- signal mean is ~0 per condition per run
        
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
            % harvest patterns for each subject
            for s = 1:D.numSim 
                Y{s}(:,fidx) = yy{s};
            end
            j = j+1; % update voxel index ticker
        end
        % return patterns to user
        varargout= {Y,partVec,condVec};
    case 'SIM_sparseChordPatterns'
        % simulates activity patterns for multi-finger chords as a function
        % of multiple single finger patterns
        [Y,partVec,condVec] = pp1_simulations('SIM_sparseVoxelTuning');
        Y = Y{1};
        chords = pp1_simulations('chords');
        runs = unique(partVec)';
        Ychords = [];
        pV_chords = [];
        cV_chords = [];
        
        for r = runs
            y  = Y(partVec==r,:);
            for i = 1:31 
                yy = y(logical(chords(i,:)),:);
                yy = sum(yy,1);
                Ychords = [Ychords;yy];
                pV_chords = [pV_chords; r];
                cV_chords = [cV_chords; i];
            end
        end
        % imagesc(Ychords(pV_chords==1,:))
        varargout = {Ychords,pV_chords,cV_chords};
        
    case 'pcm_fitSparseChords'
        % get pcm models
        M{1} = pp1_simulations('pcm_null');
        M{2} = pp1_simulations('pcm_linearNoScale');
        M{3} = pp1_simulations('pcm_noiseceiling');
        % get data
        numSubj = 10;
        for s = 1:numSubj
            [Y{s},partVec{s},condVec{s}] = pp1_simulations('SIM_sparseChordPatterns');
        end
        
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect','fixed','fitScale',1);
        [Tcross,thetaCr]     = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect','fixed','groupFit',theta_hat,'fitScale',1);
       % [Tcross2,thetaCr2]   = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect','fixed','fitScale',1);
        
        keyboard
        pcm_plotModelLikelihood(Tcross,M,'upperceil',T.likelihood(:,3),'style','bar')
    case 'pcm_null'
        % Model Null: No chord patterns
        M.type       = 'feature';
        M.numGparams = 1;
        M.name       = 'null';
        M.Ac(:,:,1)  = zeros(31);
        varargout = {M};
    case 'pcm_noiseceiling'
        % Naive averaring model- noise ceiling method 1- totall free model
%         M.type       = 'freedirect';  
%         M.numGparams = 0;
%         M.name       = 'noiseceiling';
        
        M.type       = 'freechol';  
        M.numCond    = 31;
        M.name       = 'noiseceiling';
        M            = pcm_prepFreeModel(M);
        varargout = {M};
    case 'pcm_linearNoScale'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns
        M.type       = 'feature';
        M.numGparams = 1;
        M.name       = 'linear_noScale';
        M.Ac(:,:,1)  = pp1_simulations('chords');
        varargout = {M};
    case 'pcm_linearScale'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns
        M.type       = 'feature';
        M.numGparams = 1;
        M.name       = 'linear_Scale';
        M.Ac(:,:,1)  = pp1_simulations('chords');
        
        varargout = {M};
end