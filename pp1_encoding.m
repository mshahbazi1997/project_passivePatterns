function varargout=pp1_encoding(what,varargin)
%% details

% Directories for data
dataDir = '/Users/sarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_paper/pp1_encodingAnalyses';
%dataDir = '/Users/sarbuckle/DATA/passivePatterns1/fmri/RegionOfInterest'; % where betas are saved
%dataDir = '/Users/jdiedrichsen/Dropbox (Diedrichsenlab)/projects/encodingAnalyses'; % where betas are saved

subj_name = {'pd01','pd02','s01','s02','s03','s04','s05','s06','s07','s08','s09'}; 
% pd01 is missing fieldmaps, so analyses are on pd02-s09

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

%% cases
switch what   
    case 'plot_ldc'
        % plots avg. distance b/t specified conditions for specified rois
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiPlotNum = [3 4 5 6 1 2];
        conds = 1:5; % get paired distances between these conditions
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_distances.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        D = [];
        d = [];
        for ii = 1:size(T.sn,1)
            RDM = rsa_squareRDM(T.ldc(ii,:));
            RDM = RDM(conds,conds);
            d.rdm = rsa_vectorizeRDM(RDM);
            d.ldc = mean(d.rdm);
            d.roi = T.roi(ii);
            d.sn  = T.sn(ii);
            d.conds = conds;
            D = addstruct(D,d);
        end
    
        D.roiPlotNum = roiPlotNum(D.roi)';
        
        % plot single subject lines:
        CAT.markertype = 'none';
        CAT.linecolor  = {[0.75 0.75 0.75]};
        lineplot(D.roiPlotNum,D.ldc,'CAT',CAT,'split',D.sn,'errorfcn','');
        % plot group avg. line:
        hold on
        CAT.markertype = 'o';
        CAT.linewidth  = 2;
        CAT.errorwidth = 1.5;
        CAT.linecolor  = {[0 0 0]};
        CAT.markerfill = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        lineplot(D.roiPlotNum,D.ldc,'CAT',CAT,'errorfcn','stderr');
        drawline(0,'dir','horz','linestyle','-'); 
        xlabel('Brodmann area');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2'},'xticklabelrotation',0);
        ylabel('dissimilarity (a.u.)');
        title('distances between single finger patterns');
        hold off
        
        varargout = {D};
    case 'plot_selectivity'
        % plots normalized single finger selectivity values
        % Selectivity is from fthresholded voxels (responsive to any single
        % finger(s))
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiPlotNum = [3 4 5 6 1 2];
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_selectivity_fthres.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % compute normalized selectivity:
        Tdata   = getrow(T,T.isEV==0);
        Tgauss  = getrow(T,T.isEV==1);
        Tsparse = getrow(T,T.isEV==2);
        Tdata.sft_norm = Tdata.sft - Tgauss.sft;
        Tdata.sft_norm = Tdata.sft_norm./ (Tsparse.sft - Tgauss.sft);
        T=Tdata; clear Tgauss Tsparse Tdata
        T.roiPlotNum = roiPlotNum(T.roi)';
        
        % plot single subject lines:
        CAT.markertype = 'none';
        CAT.linecolor  = {[0.75 0.75 0.75]};
        lineplot(T.roiPlotNum,T.sft_norm,'CAT',CAT,'split',T.sn,'errorfcn','');
         % plot group avg. line:
        hold on
        CAT.markertype = 'o';
        CAT.linewidth  = 2;
        CAT.errorwidth = 1.5;
        CAT.linecolor  = {[0 0 0]};
        CAT.markerfill = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        lineplot(T.roiPlotNum,T.sft_norm,'CAT',CAT,'errorfcn','stderr');
        drawline(0,'dir','horz','linestyle','-'); 
        xlabel('Brodmann area');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2'},'xticklabelrotation',0);
        ylabel('normalized selectivity');
        title('single finger selectivity');
        hold off
        
        varargout = {T};
    case 'plot_encodingFits'
        % plots the normalized encoding model fits
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiPlotNum = [3 4 5 6 1 2];
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==1),ones(5,1)); % null is model 1
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==5),ones(5,1)); % lower noise ceiling is model 5
        T = getrow(T,T.model~=3); % drop the tanh model fits
        T.roiPlotNum = roiPlotNum(T.roi)';
        % plot group avg. line:
        CAT.markertype = 'o';
        CAT.linewidth  = 2;
        CAT.errorwidth = 1.5;
        CAT.linecolor  = {[0 0 0.9],[0.9 0 0]};
        CAT.markerfill = CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        lineplot(T.roiPlotNum,T.r_norm,'CAT',CAT,'errorfcn','stderr','split',T.model,'subset',T.model~=1 & T.model~=5);
        drawline(0,'dir','horz','linestyle','-');
        legend({'linear','flexible'});
        xlabel('Brodmann area');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2'},'xticklabelrotation',0);
        ylabel('normalized Pearson''s R');
        title('encoding model fits');
        hold off
        
        varargout = {T};
    case 'plot_pcmFits'
        % plots the normalized encoding model fits
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiPlotNum = [3 4 5 6 1 2];
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_pcmFits_passive.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized log likelihoods (0=null, 1=lower noise ceiling)
        T.like_norm = T.like - kron(T.like(T.model==1),ones(4,1)); % null is model 1
        T.like_norm = T.like_norm./kron(T.like_norm(T.model==4),ones(4,1)); % lower noise ceiling is model 4
        T.roiPlotNum = roiPlotNum(T.roi)';
        % plot group avg. line:
        CAT.markertype = 'o';
        CAT.linewidth  = 2;
        CAT.errorwidth = 1.5;
        CAT.linecolor  = {[0 0 0.9],[0.9 0 0]};
        CAT.markerfill = CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        lineplot(T.roiPlotNum,T.like_norm,'CAT',CAT,'errorfcn','stderr','split',T.model,'subset',T.model~=1 & T.model~=4);
        drawline(0,'dir','horz','linestyle','-');
        legend({'linear','flexible'});
        xlabel('Brodmann area');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2'},'xticklabelrotation',0);
        ylabel('pseudo R2 (normalized log like)');
        title('pcm model fits');
        hold off
        
        varargout = {T};
    
    case 'do_encoding_notPatterns'
        % case to fit models to participant data
        % encoding models are fit on the overall activity across chords
        % with the same # fingers (the marginal activities).
        % Here, we simply calculate fits for 5 values with 3 models:
        % linear, flexible, and lower noise ceiling.
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        estimateU = 'ridgeGgroup'; % 'ols' 'blupI' 'blupG'
        vararginoptions(varargin,{'roi','sn','estimateU'});
        numModels = 3;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % define Ytrain estimation:
        switch estimateU
            case 'ols'
                modelG=[];
                estFcn = @(x,y,z,g) pp1_encoding('estU_ols',x,y,z,g);
                estMethod = 1;
            case 'ridgeI' % ridge regression 
                modelG = eye(31);
                estFcn = @(x,y,z,g) pp1_encoding('estU_ridge',x,y,z,g);
                estMethod = 2;
            case 'ridgeGgroup' % group crossval G prior
                modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
                estFcn = @(x,y,z,g) pp1_encoding('estU_ridge',x,y,z,g);
                estMethod = 3;
        end
        % loop through subjects and fit each individually:
        D=[]; % output
        for s=1:numel(sn)
            fprintf('S%02d\n',s);
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
            G_pred = zeros(5,5,numModels);
            Y_avg = nan([numModels,5,numPart]); % avg. activity per # digits
            modelTheta = {};
            SS1 = zeros(numPart,numModels);
            SS2 = SS1; SSC = SS1; SS1t = SS1; SS2t = SS1; SSCt = SS1; 
            RSS = SS1; TSS=SS1; RSSt=SS1; TSSt=SS1;
            % loop through partitions and estimate patterns
            for ii=1:numel(partI)
                % estimate the true condition activity patterns using
                % training data:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii}); 
                Utrain   = estFcn(Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                Utrain   = mean(numD_inv*Utrain,2);
                Usf      = Utrain(1); % average single finger activity
                Ytest    = mean(numD_inv*Y{s}(testIdx,:),2);
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % linear summation model
                            theta0 = [0 1];
                            thetaFcn = @(x) modelLossRSS_activity(x,Usf,Utrain,'summation'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelActivity',Usf,'summation',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaScale = thetaEst(2);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 2 % flexible summation model
                            theta0 = [thetaBaseline log(2) log(3) log(4) log(5) thetaScale];
                            thetaFcn = @(x) modelLossRSS_activity(x,Usf,Utrain,'summation_flexible'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelActivity',Usf,'summation_flexible',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5)) thetaEst(6)];
                        case 3 % lower noise ceiling
                            Ypred = Utrain;
                            thetaEst = [nan nan nan nan nan nan];
                    end

                    modelTheta{mm}(ii,:) = thetaEst;
                    Y_avg_cent(mm,:,ii)  = Ypred'-mean(Ypred,1); % avg. activity per # digits
                    Y_avg(mm,:,ii)       = Ypred'; % avg. activity per # digits
                    % calculate metrics for R and R2 of prediction against TRAINING
                    % data:
                    [SS1t(ii,mm),SS2t(ii,mm),SSCt(ii,mm),RSSt(ii,mm),TSSt(ii,mm)] = pp1_encoding('evaluate',Ypred,Utrain); % corr b/t pred and train patterns
                    % calculate metrics for R and R2 of prediction against TEST data:
                    [SS1(ii,mm), SS2(ii,mm), SSC(ii,mm), RSS(ii,mm), TSS(ii,mm)] = pp1_encoding('evaluate',Ypred,Ytest); % corr b/t pred and test patterns
                end
            end
            d = [];
            % for each model, avg. model thetas across CV partitions:
            d.modelTheta = cell2mat(cellfun(@(x) mean(x,1),modelTheta,'uni',0)');
            % pearson R:
            d.r_train = [mean(SSCt./sqrt(SS1t.*SS2t))]';
            d.r_test  = [mean(SSC./sqrt(SS1.*SS2))]';
            % R2:
            d.r2_train = [1-sum(RSSt)./sum(TSSt)]';
            d.r2_test  = [1-sum(RSS)./sum(TSS)]';
            % arrange data into output structure:
            v = ones(numModels,1);
            d.avgAct= mean(Y_avg,3); % avg. activity per # digits
            d.avgAct_cent= mean(Y_avg_cent,3); % avg. activity per # digits

            d.model = [1:numModels]';
            d.roi   = v.*roi;
            d.sn    = v.*sn(s);
            d.estU  = v.*estMethod;
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};  
    case 'predictModelActivity'
        % factorization of encoding models:
        Ysf = varargin{1}; % training patterns (usually single finger U estimates from training data)
        model = varargin{2};
        theta = varargin{3}; % theta(1) = baseline, theta(2:end) = model params
        switch model
            case 'summation'
                % theta(1) = baseline
                % theta(2) = overall scaling param
                
                % get multi-finger design:
                X = [1:5]';
                Ymf_hat = theta(2) * (X*(Ysf-theta(1)) + theta(1));
            case 'summation_flexible'
                % theta(1) = baseline
                % theta(2:5) = finger combination param (per # fingers in
                % chords for 2:5 digits)
                % theta(6) = overall scaling param
                
                % flexible scaling of patterns per # finger stimulated
                % theta(1) = baseline param
                % theta(2:5) = finger scalar params, one per # fingers stimulated
                % get multi-finger design:
                
                X = [1 exp(theta(2:5))]'; % flexible scaling per # fingers stimulated (force positive values with exp)
                Ymf_hat = theta(6) * (X*(Ysf-theta(1)) + theta(1)); 
        end
        varargout = {Ymf_hat};
    case 'do_encoding'
        % case to fit models to participant data
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        estimateU = 'ridgeGgroup'; % 'ols' 'blupI' 'blupG'
        vararginoptions(varargin,{'roi','sn','estimateU'});
        numModels = 5;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % define Ytrain estimation:
        switch estimateU
            case 'ols'
                modelG=[];
                estFcn = @(x,y,z,g) pp1_encoding('estU_ols',x,y,z,g);
                estMethod = 1;
            case 'ridgeI' % ridge regression 
                modelG = eye(31);
                estFcn = @(x,y,z,g) pp1_encoding('estU_ridge',x,y,z,g);
                estMethod = 2;
            case 'ridgeGgroup' % group crossval G prior
                modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
                estFcn = @(x,y,z,g) pp1_encoding('estU_ridge',x,y,z,g);
                estMethod = 3;
        end
        % loop through subjects and fit each individually:
        D=[]; % output
        for s=1:numel(sn)
            fprintf('S%02d\n',s);
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
            Y_avg = nan([numModels,5,numPart]); % avg. activity per # digits
            modelTheta = {};
            SS1 = zeros(numPart,numModels);
            SS2 = SS1; SSC = SS1; SS1t = SS1; SS2t = SS1; SSCt = SS1; 
            RSS = SS1; TSS=SS1; RSSt=SS1; TSSt=SS1;
            % loop through partitions and estimate patterns
            for ii=1:numel(partI)
                % estimate the true condition activity patterns using
                % training data:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii}); 
                Utrain   = estFcn(Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                Usf      = Utrain(1:5,:); % single finger Us
                Ytest    = Y{s}(testIdx,:);
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaEst = [nan nan nan nan nan nan];
                        case 2 % linear summation model
                            theta0 = [0 1];
                            thetaFcn = @(x) modelLossRSS(x,Usf,Utrain,'summation'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Usf,'summation',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaScale = thetaEst(2);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % tanh summation model
                            theta0 = [thetaBaseline 1 thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Usf,Utrain,'summation_tanh'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Usf,'summation_tanh',thetaEst);
                            thetaEst = [thetaEst nan nan nan];
                        case 4 % flexible summation model
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5) thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Usf,Utrain,'summation_flexible'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Usf,'summation_flexible',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5)) thetaEst(6)];
                        case 5 % lower noise ceiling
                            Ypred = Utrain;
                            thetaEst = [nan nan nan nan nan nan];
                    end
                    % scale model predictions to left-out test data:
%                     scaling = (Ypred(:)'*Ytest(:))/(Ypred(:)'*Ypred(:));
%                     Ypred = Ypred.*scaling;
                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    Y_avg_cent(mm,:,ii)  = mean(numD_inv*(Ypred-mean(Ypred,1)),2)'; % avg. activity per # digits
                    Y_avg(mm,:,ii)       = mean(numD_inv*Ypred,2)'; % avg. activity per # digits
                    % calculate metrics for R and R2 of prediction against TRAINING
                    % data:
                    [SS1t(ii,mm),SS2t(ii,mm),SSCt(ii,mm),RSSt(ii,mm),TSSt(ii,mm)] = pp1_encoding('evaluate',Ypred,Utrain); % corr b/t pred and train patterns
                    % calculate metrics for R and R2 of prediction against TEST data:
                    [SS1(ii,mm), SS2(ii,mm), SSC(ii,mm), RSS(ii,mm), TSS(ii,mm)] = pp1_encoding('evaluate',Ypred,Ytest); % corr b/t pred and test patterns
                end
            end
            G_pred = G_pred./numel(partI); % avg. model predicted Gs across folds
            d = [];
            % for each model, avg. model thetas across CV partitions:
            d.modelTheta = cell2mat(cellfun(@(x) mean(x,1),modelTheta,'uni',0)');
            % pearson R:
            d.r_train = [mean(SSCt./sqrt(SS1t.*SS2t))]';
            d.r_test  = [mean(SSC./sqrt(SS1.*SS2))]';
            % R2:
            d.r2_train = [1-sum(RSSt)./sum(TSSt)]';
            d.r2_test  = [1-sum(RSS)./sum(TSS)]';
            % arrange data into output structure:
            v = ones(numModels,1);
            d.avgAct= mean(Y_avg,3); % avg. activity per # digits
            d.avgAct_cent= mean(Y_avg_cent,3); % avg. activity per # digits
            for mm=1:numModels
                d.gpred(mm,:) = rsa_vectorizeIPM(G_pred(:,:,mm));
            end
            d.model = [1:numModels]';
            d.roi   = v.*roi;
            d.sn    = v.*sn(s);
            d.estU  = v.*estMethod;
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    case 'predictModelPatterns'
        % factorization of encoding models:
        Ysf = varargin{1}; % training patterns (usually single finger U estimates from training data)
        model = varargin{2};
        theta = varargin{3}; % theta(1) = baseline, theta(2:end) = model params
        switch model
            case 'null'
                % Model predicts overall scaling of avg. activity with #
                % fingers. Scaling matches true mean scaling in training
                % data.
                % Ysf here are all 31 conditions from the training data
                % Set each condition pattern to the be mean pattern for all
                % chords with the same number of fingers.
                chords = pp1_encoding('chords'); 
                X = pcm_indicatorMatrix('identity',sum(chords,2)); % which patterns have the same # of fingers?
                Ymf_hat = X*pinv(X)*Ysf;
            case 'summation'
                % theta(1) = baseline
                % theta(2) = overall scaling param
                
                % get multi-finger design:
                X = pp1_encoding('chords');
                Ymf_hat = theta(2) * (X*(Ysf-theta(1)) + theta(1));
            case 'summation_tanh'
                % theta(1) = baseline
                % theta(2) = tanh scaling param
                % theta(3) = overall scaling param
                
                % y = tanh(x*theta2)/theta2
                % theta 2 is a scalar param to scale the patterns
%                 X = pp1_encoding('chords');
%                 Ymf_hat = X*(Ysf-theta(1)) + theta(1);
%                 Ymf_hat = tanh(Ymf_hat.*theta(2))./theta(2);
                X = pp1_encoding('chords');
                Ymf_hat = X*(Ysf-theta(1));
                Ymf_hat = theta(3) * (tanh(Ymf_hat.*theta(2))+theta(1));
            case 'summation_flexible'
                % theta(1) = baseline
                % theta(2:5) = finger combination param (per # fingers in
                % chords for 2:5 digits)
                % theta(6) = overall scaling param
                
                % flexible scaling of patterns per # finger stimulated
                % theta(1) = baseline param
                % theta(2:5) = finger scalar params, one per # fingers stimulated
                % get multi-finger design:
                
                X = pp1_encoding('chords');
                numD = sum(X,2);
                X = X.*[ones(1,5) exp(theta(numD(numD>1)))]'; % flexible scaling per # fingers stimulated (force positive values with exp)
                Ymf_hat = theta(6) * (X*(Ysf-theta(1)) + theta(1)); 
        end
        varargout = {Ymf_hat};
    case 'estU_ols' % estimate Us using ols
        % estimate condition Us by averaging patterns across runs
        
        % inputs
        Y   = varargin{1}; % matrix of activity patterns [#conds*#runs x #vox]
        pV  = varargin{2}; % partition vector (assume cV and pV are same across subjs)
        cV  = varargin{3}; % condition vector (chord #s)

        C0 = pcm_indicatorMatrix('identity',pV);
        X = pcm_indicatorMatrix('identity',cV);
        % Y = Y-C0*pinv(C0)*Y; % rmv run means (random effect)
        U = pinv(X)*Y; % avg. patterns per condition across runs
        G = U*U' ./ size(U,2);
        varargout = {U,G};                      
    case 'estU_ridge' % estimate patterns using ridge regression with a model prior
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
        % create PCM model structure M:
        
        % double-centre G:
        % H = eye(31)-ones(31)/31;
        % G=H*G*H; 
        % remove run mean from patterns:
        % X0 = pcm_indicatorMatrix('identity',pV);
        % Y = Y-X0*pinv(X0)*Y; % rmv run means (fixed effect)
        M{1}.type = 'component';
        M{1}.Gc=G;
        M{1}.numGparams = 1;
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,cV,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using ridge regression:
        Z = pcm_indicatorMatrix('identity',cV);
        U = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);
        
        % Do it the ridge regression route: 
        % reconstruct true patterns using ridge regression:
        % lambda = exp(theta_hat{1}(2))/exp(theta_hat{1}(1)); %lambda is noise/scale
        % U = ((Z'*Z + lambda*inv(G))^-1)*Z'*Y;

        varargout = {U};           
        
    case 'pcm_fitGroup'
        % fits pcm models to data from one region
        sn     = [2:11];
        glm    = 4;
        roi    = [];
        vararginoptions(varargin,{'sn','roi'});

        % load subject data:
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % load group G:
        Gcv = pp1_encoding('getRegionG','roi',roi,'glm',glm);
        Gsf = Gcv(1:5,1:5); % group avg. single-finger G
        
        % define group-level models:
        % 1. null model- fits overall scaling of avg. activity
%         M{1}.name    = 'null scaling';
%         M{1}.type       = 'component';
%         M{1}.numGparams = 1;
%         M{1}.Gc(:,:,1)  = pp1_encoding('pcm_Gnull');
        M{1}.name = 'null';
        M{1}.type       = 'nonlinear';
        M{1}.numGparams = 4;
        M{1}.theta0     = log([0.9 0.8 0.7 0.6])';
        M{1}.modelpred  = @pcmGroup_modelpred_flexibleNULL;
        
        % 2. linear model
        chords = pp1_encoding('chords');
        M{2}.name = 'linear';
        M{2}.type       = 'component';
        M{2}.numGparams = 1;
        M{2}.Gc(:,:,1)  = chords*Gsf*chords';
        
        % 3. flexible scaling linear model
        M{3}.name = 'flexible';
        M{3}.type       = 'nonlinear';
        M{3}.numGparams = 4;
        M{3}.theta0     = log([0.8 0.6 0.4 0.2])';
        M{3}.Ac         = pcm_diagonalize(Gsf); 
        M{3}.modelpred  = @pcmGroup_modelpred_flexible;
       
        % 4. noise ceiling model
        M{4}.name = 'noiseceiling';
        M{4}.type       = 'freedirect';
        M{4}.numGparams = 0;
        M{4}.theta0     = [];
                
        % choose proper way to deal with run effects
        runEffect = 'random'; % mean patterns not removed
        % fit all models
        [T,theta,G_pred]        = pcm_fitModelGroup(Y,M,pV,cV,'runEffect',runEffect,'fitScale',1,'isCheckDeriv',0,'verbose',1); % perform derivative checks (to ensure I've done correct deriv implementation for nonlinear models)
        [Tcv,theta_cv,Gcv_pred] = pcm_fitModelGroupCrossval(Y,M,pV,cV,'runEffect',runEffect,'groupFit',theta,'fitScale',1,'verbose',1,'isCheckDeriv',0);
        
        % normalize fits between null model and upper noise ceiling:
        Tcv.likelihood_norm = bsxfun(@minus,Tcv.likelihood,Tcv.likelihood(:,1)); % set null model likelihood==0
        Tcv.likelihood_norm = bsxfun(@rdivide,Tcv.likelihood_norm,T.likelihood(:,4)-Tcv.likelihood(:,1)); % set upper noise ceiling model likelihood==1
        
        % arrange into output structure:
        D=[];
        for ii=1:numel(sn)
            for mm=1:numel(M)
                d.roi = roi;
                d.sn = sn(ii);
                d.model = mm;
                
                d.noise = T.noise(ii,mm);
                d.noise_cv = Tcv.noise(ii,mm);
                d.scale = T.scale(ii,mm);
                d.scale_cv = Tcv.scale(ii,mm);
                d.run = T.run(ii,mm);
                d.run_cv = Tcv.run(ii,mm);
                if mm<4
                    d.theta_cv = {exp(theta_cv{mm}(:,ii)')};
                else
                    d.theta_cv = {[]};
                end
                d.gpred = rsa_vectorizeIPM(Gcv_pred{mm}(:,:,ii));
                
                d.like = T.likelihood(ii,mm);
                d.like_cv = Tcv.likelihood(ii,mm);
                d.like_norm = Tcv.likelihood_norm(ii,mm);
                
                D=addstruct(D,d);
            end
        end

        varargout = {D,Tcv,T,M,theta_cv,Gcv_pred};    
    case 'pcm_Gnull'
        sn     = [2:11];
        glm    = 4;
        roi    = 2;
        vararginoptions(varargin,{'sn','glm','roi'});
        % load subject data:
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % load group G:
        Gcv = pp1_encoding('getRegionG','roi',roi,'glm',glm);
        % estimate true Us:
        Gnull = zeros(31,31,numel(sn));
        for ii=1:numel(sn)
            Upred = pp1_encoding('estU_ridge',Y{ii},pV{ii},cV{ii},Gcv);
            Ypred = pp1_encoding('predictModelPatterns',Upred,'null',[]);
            Gnull(:,:,ii) = (Ypred*Ypred')./size(Ypred,2);
        end
        Gnull = mean(Gnull,3);
        varargout = {Gnull};
        
    case 'compare_activePassive'
%         roi = 2;
%         D1 = pp1_encoding('do_active','roi',roi);
%         D2 = pp1_encoding('do','roi',roi);
%         D1.dataset = ones(size(D1.sn));
%         D2.dataset = ones(size(D2.sn)).*2;
%         D=addstruct(D1,D2);
%         save(fullfile(dataDir,'active_passive_fits.mat'),'-struct','D');
        D=load(fullfile(dataDir,'active_passive_fits.mat'));
        
        titles={'3b active (n=8)','3b passive n=10)'};
        labels = {'null','linear','tanh','flexible','nceil'};
        for ii=1:2
            subplot(1,2,ii);
            barplot([D.roi D.model],D.r_test,'subset',D.dataset==ii);
            title(titles{ii});
            ylabel('Pearson''s R (on left-out test data)');
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
        end
        
        varargout = {D};
    case 'do_active'
        % case to fit models to participant data
        sn  = 1:8; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 1;
        estimateU = 'ridgeGgroup'; % 'ols' 'blupI' 'blupG'
        vararginoptions(varargin,{'roi','sn','estimateU'});
        numModels = 5;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData_active','sn',sn,'roi',roi,'glm',glm);
        % define Ytrain estimation:
        switch estimateU
            case 'ols'
                modelG=[];
                estFcn = @(x,y,z,g) pp1_encoding('estU_ols',x,y,z,g);
                estMethod = 1;
            case 'ridgeI' % ridge regression 
                modelG = eye(31);
                estFcn = @(x,y,z,g) pp1_encoding('estU_ridge',x,y,z,g);
                estMethod = 2;
            case 'ridgeGgroup' % group crossval G prior
                modelG = pp1_encoding('getRegionG_active','roi',roi,'glm',glm);
                estFcn = @(x,y,z,g) pp1_encoding('estU_ridge',x,y,z,g);
                estMethod = 3;
        end
        % loop through subjects and fit each individually:
        D=[]; % output
        for s=1:numel(sn)
            fprintf('S%02d\n',s);
            % split data into partitions (for estimation of patterns under each model)
            % assign runs to each partition
            part = unique(pV{s});
            numPart = numel(part);
            partI = {};
            for ip=1:numPart
                partI{ip}=part(ip);
            end  
            % allocate space for predicted patterns, evaluation metrics:
            Y_hat = nan([size(Y{s}), numModels]);
            Y_avg = nan([numModels,5,numPart]); % avg. activity per # digits
            modelTheta = {};
            SS1 = zeros(numPart,numModels);
            SS2 = SS1; SSC = SS1; SS1t = SS1; SS2t = SS1; SSCt = SS1; 
            RSS = SS1; TSS=SS1; RSSt=SS1; TSSt=SS1;
            % loop through partitions and estimate patterns
            for ii=1:numel(partI)
                % estimate the true condition activity patterns using
                % training data:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii}); 
                Utrain   = estFcn(Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                Usf      = Utrain(1:5,:); % single finger Us
                Ytest    = Y{s}(testIdx,:);
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaEst = [nan nan nan nan nan nan];
                        case 2 % linear summation model
                            theta0 = [0 1];
                            thetaFcn = @(x) modelLossRSS(x,Usf,Utrain,'summation'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Usf,'summation',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaScale = thetaEst(2);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % tanh summation model
                            theta0 = [thetaBaseline 1 thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Usf,Utrain,'summation_tanh'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Usf,'summation_tanh',thetaEst);
                            thetaEst = [thetaEst nan nan nan];
                        case 4 % flexible summation model
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5) thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Usf,Utrain,'summation_flexible'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Usf,'summation_flexible',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5)) thetaEst(6)];
                        case 5 % lower noise ceiling
                            Ypred = Utrain;
                            thetaEst = [nan nan nan nan nan nan];
                    end
                    modelTheta{mm}(ii,:) = thetaEst;
                    Y_avg(mm,:,ii)       = mean(numD_inv*(Ypred-mean(Ypred,1)),2)'; % avg. activity per # digits
%                     Y_avg(mm,:,ii)       = mean(numD_inv*Ypred,2)'; % avg. activity per # digits
                    Y_hat(testIdx,:,mm)  = Ypred;
                    % calculate metrics for R and R2 of prediction against TRAINING
                    % data:
                    [SS1t(ii,mm),SS2t(ii,mm),SSCt(ii,mm),RSSt(ii,mm),TSSt(ii,mm)] = pp1_encoding('evaluate',Ypred,Utrain); % corr b/t pred and train patterns
                    % calculate metrics for R and R2 of prediction against TEST data:
                    [SS1(ii,mm), SS2(ii,mm), SSC(ii,mm), RSS(ii,mm), TSS(ii,mm)] = pp1_encoding('evaluate',Ypred,Ytest); % corr b/t pred and test patterns
                end
            end
            d = [];
            % for each model, avg. model thetas across CV partitions:
            d.modelTheta = cell2mat(cellfun(@(x) mean(x,1),modelTheta,'uni',0)');
            % pearson R:
            d.r_train = [mean(SSCt./sqrt(SS1t.*SS2t))]';
            d.r_test  = [mean(SSC./sqrt(SS1.*SS2))]';
            % R2:
            d.r2_train = [1-sum(RSSt)./sum(TSSt)]';
            d.r2_test  = [1-sum(RSS)./sum(TSS)]';
            % arrange data into output structure:
            v = ones(numModels,1);
            d.avgAct= mean(Y_avg,3); % avg. activity per # digits
            d.model = [1:numModels]';
            d.roi   = v.*roi;
            d.sn    = v.*sn(s);
            d.estU  = v.*estMethod;
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    case 'pcm_fitGroup_active'  
    case 'getData_active'
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
    case 'makeRegionG_active'
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
    case 'getRegionG_active'
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
    
    case '0' % helper cases:
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
        betaType = 'betaW'; % multivariately prewhitened (or betaUW, raw_beta)
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
        end
        varargout = {Y,partVec,condVec,Gcv};     
    case 'makeRegionG'
        % Estimate Region G (G is avg. semi-positive definite crossval G across
        % participants from roi):
        sn  = 2:11;
        glm = 4;
        roi = [1:12,15,16,17,18,21,23,26];
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
        save(fullfile(dataDir,sprintf('glm%d_regionG.mat',glm)),'-struct','D');   
    case 'getRegionG'
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
    case 'evaluate'
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
        
    case '0' % cases to test/compare pattern estimation approaches:    
    case 'test_U_ridgeI' % compare ols and ridge regression with identity prior
        % create scatterplot of patterns estimated using ols and ridge regression with identity prior
        % Both estimation approaches should yeild identical values, scaled
        % by a constant offset
        sn=10;
        roi = 2; % data from Ba 3b
        glm = 4;
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        Utrain_ols = pp1_encoding('estU_ols',Y{1},pV{1},cV{1});
        Utrain_ridge = pp1_encoding('estU_ridge',Y{1},pV{1},cV{1},eye(31));
        scatter(Utrain_ols(:),Utrain_ridge(:));
        xlabel('ols');
        ylabel('ridge $\bf{I}$','Interpreter','latex');
        title('relationship between pattern estimates');
    case 'test_U_ridgeG' % compare ridge regression with Identity prior and with group G as prior
        % create scatterplot of patterns estimated using ridge regression
        % with different prior distributions.
        sn=10;
        roi = 2; % data from Ba 3b
        glm = 4;
        G = pp1_encoding('getRegionG','roi',roi,'glm',glm);
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        Utrain_G  = pp1_encoding('estU_ridge',Y{1},pV{1},cV{1},G);
        Utrain_I  = pp1_encoding('estU_ridge',Y{1},pV{1},cV{1},eye(31));
        scatter(Utrain_I(:),Utrain_G(:));
        xlabel('ridge $\bf{I}$','Interpreter','latex');
        ylabel('ridge $\bf{G}$','Interpreter','latex');
        title('relationship between pattern estimates');
    case 'test_noiseceiling' % compare fits of noise ceiling data across different pattern estimation approaches
        % case to get correlation fits of noise ceilings using different
        % methods to estimate Us on each crossval fold (ols, ridge, ridge
        % with model Gs):
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % loop through subjects and fit each individually:
        D=[]; % output
        v=ones(4,1); % helper array
        for s=1:numel(sn)
            % get CV folds
            part = unique(pV{s});
            numPart = numel(part);
            partI = {};
            for ip=1:numPart
                partI{ip}=part(ip);
            end  
            % loop through folds, and calculate noise ceiling fits for each
            % fold using 3 different estimation approaches for training Us
            for ii=1:numel(partI)
                % estimate the true condition activity patterns using
                % training data:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii});
                Ytest    = Y{s}(testIdx,:);
                for jj=1:4  % each of the 3 methods to esitmate training Us
                    switch jj
                        case 1 % estimate training patterns using OLS
                            Ypred = pp1_encoding('estU_ols',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx));
                        case 2 % estimate training patterns using ridge regression with Identity prior
                            modelG = eye(31);
                            Ypred = pp1_encoding('estU_ridge',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                        case 3 % " ridge regression with cross-val G estimated from subject's training data
                            modelG = pcm_makePD(pcm_estGCrossval(Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx)));
                            Ypred = pp1_encoding('estU_ridge',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                        case 4 % " ridge regression with group average cross-val G from this ROI
                            modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
                            Ypred = pp1_encoding('estU_ridge',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                    end
                    % evaluate estimated Us (Ypred) against Ytest patterns:
                    [SS1(ii,jj), SS2(ii,jj), SSC(ii,jj), RSS(ii,jj), TSS(ii,jj)] = pp1_encoding('evaluate',Ypred,Ytest);
                end
            end
            d.sn  = v.*sn(s);
            d.roi = v.*roi;
            d.glm = v.*glm;
            d.estMethod = [1:jj]';
            d.r   = [mean(SSC./sqrt(SS1.*SS2))]';
            d.r2  = [1-sum(RSS)./sum(TSS)]';
            D=addstruct(D,d);
        end
        % plot fits:
        barplot(D.estMethod, D.r);
        ylabel('Pearson''s R');
        xlabel('estimation of training data');
        set(gca,'xticklabel',{'OLS','ridge I','ridge Gtrain','ridge Ggroup'},'xticklabelrotation',45);
        varargout={D};
    
    case '0' % cases to test parameter estimation:
    case 'test_modelLinear' % test parameter estimation of linear model baseline param
        % simulate noiseless data with ground-truth baseline under linear model
        % attempt to estimate real baseline using fminsearch and minimize
        
        % get some single finger patterns:
        sn=10;
        roi = 2; % data from Ba 3b
        glm = 4;
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        Ysf = pp1_encoding('estU_ols',Y{1},pV{1},cV{1});
        Ysf = Ysf(1:5,:); % use only single finger patterns
        
        % generate multi-finger patterns under linear model (no noise)
        maxIter = 10000;
        thetaTrue = -0.05;
        Y = pp1_encoding('predictModelPatterns',Ysf,'summation',thetaTrue);
        
        % use fminsearch to find theta estimate:
        theta0 = 0;
        thetaFcn1 = @(x) modelLossRSS(x,Ysf,Y,'summation'); % minimize pattern RSS in parameter fitting
        thetaEst1 = fminsearch(thetaFcn1, theta0, optimset('MaxIter',maxIter));
        Ypred1 = pp1_encoding('predictModelPatterns',Ysf,'summation',thetaEst1);
        
        % plot mean activity of simulated patterns and estimated patterns:
        Cd = pinv(pcm_indicatorMatrix('identity',sum(pp1_encoding('chords'),2)));
        plot(mean(Cd*Y,2),'linewidth',1.25); % simulated patterns
        hold on;
        plot(mean(Cd*Ypred1,2),'linewidth',1.25,'linestyle','-.','color','r'); % estimated patterns
        xlabel('# digits in chord');
        ylabel('avg. activity');
        
        % use minimize to find theta estimate:
        thetaFcn2 = @(x) model_linearGradient(x,Ysf,Y); % minimize pattern RSS in parameter fitting
        [thetaEst2,~,fi] = minimize(theta0, thetaFcn2, maxIter);
        % % NOTE: my gradient calculation appears incorrect- minimize halts
        % after 2 iterations..
        
        keyboard
    case 'test_modelTanh' % test parameter estimation of tanh model baseline param
        % simulate noiseless data with specified params under tanh model
        % attempt to estimate params using fminsearch
        
        % get some single finger patterns:
        sn=10;
        roi = 2; % data from Ba 3b
        glm = 4;
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        Ysf = pp1_encoding('estU_ols',Y{1},pV{1},cV{1});
        Ysf = Ysf(1:5,:); % use only single finger patterns
        
        % generate multi-finger patterns under tanh model (no noise)
        maxIter = 10000;
        thetaTrue = [-0.25 10];
        Y = pp1_encoding('predictModelPatterns',Ysf,'summation_tanh',thetaTrue);
        
        % use fminsearch to find theta estimate:
        theta0 = [0 1];
        thetaFcn1 = @(x) modelLossRSS(x,Ysf,Y,'summation_tanh'); % minimize pattern RSS in parameter fitting
        thetaEst1 = fminsearch(thetaFcn1, theta0, optimset('MaxIter',maxIter));
        Ypred1 = pp1_encoding('predictModelPatterns',Ysf,'summation_tanh',thetaEst1);
        
        % plot mean activity of simulated patterns and estimated patterns:
        Cd = pinv(pcm_indicatorMatrix('identity',sum(pp1_encoding('chords'),2)));
        plot(mean(Cd*Y,2),'linewidth',1.25); % simulated patterns
        hold on;
        plot(mean(Cd*Ypred1,2),'linewidth',1.25,'linestyle','-.','color','r'); % estimated patterns
        xlabel('# digits in chord');
        ylabel('avg. activity');
        
        keyboard
        
        
        
    case '0' %-------------------------------------------------------------    
    case '0' % DEPRECIATED:
    case 'ENC:do_depreciated2'
        % case to fit models to participant data
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        vararginoptions(varargin,{'roi','verbose'});
        
        D=[]; % output
        chords = pp1_encoding('chords');
        patternModel = {'summation','summation_tanh'}; % models where we actually estimate the patterns
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
            for ip=1:numPart
                partI{ip}=part(ip);
            end  
            % loop through partitions and estimate patterns
            Y_hat = nan([size(Y{s}), numel(patternModel)+1]);
            modelTheta = {};
            SS1 = zeros(numPart,numel(patternModel)+1);
            SS2 = SS1; SSC = SS1; SS1t = SS1; SS2t = SS1; SSCt = SS1; 
            RSS = SS1; TSS=SS1; RSSt=SS1; TSSt=SS1;
            for ii=1:numel(partI)
                % estimate the true condition activity patterns using
                % training data:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii}); 
                Ytest = Y{s}(testIdx,:);
                Ytrain = pp1_encoding('ENC:estimateU',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx));
                % using estimated single-finger Us, predict data under each
                % model:
                for mm=1:numel(patternModel)
                    % estimate the model thetas:
                    fcn = @(x) modelLossRSS(x,Ytrain(1:5,:),Ytrain,patternModel{mm}); % minimize pattern RSS in parameter fitting
                    [modelTheta{mm}(ii,:),RSSt(ii,mm)] = fminsearch(fcn, modelTheta0{mm}, optimset('MaxIter',20000));
                    Ypred = pp1_encoding('predictModelPatterns',Ytrain(1:5,:),patternModel{mm},modelTheta{mm}(ii,:));
                    Y_hat(testIdx,:,mm) = Ypred;
                    % calc pearson R b/t pred and test patterns:
                    [SS1(ii,mm), SS2(ii,mm), SSC(ii,mm)] = pp1_encoding('do_corr',Ypred,Ytest); % corr b/t pred and test patterns
                    [SS1t(ii,mm),SS2t(ii,mm),SSCt(ii,mm)] = pp1_encoding('do_corr',Ypred,Ytrain); % corr b/t pred and train patterns
                    % calc R2:
                    RSS(ii,mm) = sum(sum(((Ytest-mean(Ytest,1))-(Ypred-mean(Ypred))).^2)); % rss b/t pred and test
                    TSS(ii,mm) = SS1(ii,mm);
                    TSSt(ii,mm) = SS1t(ii,mm);
                end
                mm=mm+1; % lower noise ceiling model:
                Y_hat(testIdx,:,end) = Ytrain; % patterns for lower bound of noise ceiling
                % calc simple correlation b/t pred and test patterns:
                Ypred = Ytrain;
                [SS1(ii,mm), SS2(ii,mm), SSC(ii,mm)] = pp1_encoding('do_corr',Ypred,Ytest); % corr b/t pred and test patterns
                [SS1t(ii,mm),SS2t(ii,mm),SSCt(ii,mm)] = pp1_encoding('do_corr',Ypred,Ytrain); % corr b/t pred and train patterns
                RSS(ii,mm) = sum(sum(((Ytest-mean(Ytest,1))-(Ypred-mean(Ypred))).^2)); % rss b/t pred and test
                TSS(ii,mm) = SS1(ii,mm);
                TSSt(ii,mm) = SS1t(ii,mm);
            end
            
            % define PCM models:
            M = {};
            % summation model
            M{1}.type         = 'component'; % to allow G to scale across runs
            M{1}.numGparams   = 1;
            M{1}.name         = 'summation';
            M{1}.fitAlgorithm = 'NR';
            M{1}.Gc           = pcm_estGCrossval(Y_hat(:,:,1),pV{s},cV{s});
            % summation tanh model
            M{2}.type         = 'component';
            M{2}.numGparams   = 1;
            M{2}.name         = 'tanh';
            M{2}.fitAlgorithm = 'NR';
            M{2}.Gc           = pcm_estGCrossval(Y_hat(:,:,2),pV{s},cV{s});
            
            % null scaling model (all patterns identical, mean acitvity
            % scales)
            M{3}.type         = 'component'; % to allow G to scale across runs
            M{3}.numGparams   = 1;
            M{3}.name         = 'null';
            M{3}.fitAlgorithm = 'NR';
            M{3}.Gc           = chords*ones(5)*chords';
            % lower bound noise ceiling
            M{end+1}.type       = 'component';
            M{end}.numGparams   = 1;
            M{end}.name         = 'lower noise ceiling';
            M{end}.fitAlgorithm = 'NR';
            M{end}.Gc           = pcm_estGCrossval(Y_hat(:,:,end),pV{s},cV{s});
            % overall noise ceiling
            M{end+1}.type       = 'component'; % to allow G to scale across runs
            M{end}.numGparams   = 1;
            M{end}.name         = 'upper noise ceiling';
            M{end}.fitAlgorithm = 'NR';
            M{end}.Gc           = pcm_estGCrossval(Y{s},pV{s},cV{s});
            
            [T_cv,D_cv,theta_hat] = pcm_fitModelIndividCrossval(Y{s},M,pV{s},cV{s},...
                                        'runEffect','random','verbose',0,...
                                        'evaluation',{'R','Rpool','R2','likelihood'});
            % do simple pearson with predicted patterns and test patterns:
            r_test = mean(SSC./sqrt(SS1.*SS2));
            r_test = [r_test(1:2),nan,r_test(3),nan]; % nan values for null and upper nc models
            r_train  = mean(SSCt./sqrt(SS1t.*SS2t));
            r_train  = [r_train(1:2),nan,r_train(3),nan];
            % do simple R2:
            r2_test = 1-sum(RSS)./sum(TSS);
            r2_test = [r2_test(1:2),nan,r2_test(3),nan];
            r2_train = 1-sum(RSSt)./sum(TSSt);
            r2_train = [r2_train(1:2),nan,r2_train(3),nan];
            % arrange data into output structure:
            v = ones(numel(M),1);
            d = [];
            d.model  = [1:numel(v)]';
            d.roi    = v.*roi;
            d.sn     = v.*sn(s);
            % get evaluation metrics
            for mm=1:numel(M)
                d.modelName(mm,1)= string(M{mm}.name);
                d.r_test(mm,1)   = r_test(mm);
                d.r_train(mm,1)  = r_train(mm);
                d.r2_test(mm,1)  = r2_test(mm);
                d.r2_train(mm,1) = r2_train(mm);
                d.like_pcm(mm,1) = T_cv.likelihood(mm);
                d.r2_pcm(mm,1)   = T_cv.R2(mm);
                d.r_pcm(mm,1)    = T_cv.R(mm);
                d.rpool_pcm(mm,1)= T_cv.Rpool(mm);
                d.g_pred(mm,:)   = rsa_vectorizeIPM(M{mm}.Gc);
                % save model parameters
                d.th_scale(mm,1) = exp(theta_hat{mm}(1));
                d.th_run(mm,1)   = exp(theta_hat{mm}(3));
                d.th_noise(mm,1) = exp(theta_hat{mm}(2));
                if mm<3 % models 1 and 2 have parameters fit for predicting the patterns
                    d.th_pattern(mm,1)= {mean(modelTheta{mm},1)};
                else
                    d.th_pattern(mm,1)= {[]};
                end
            end

            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};    
    case 'ENC:do_depreciated1'
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
                    fcn = @(x) modelLossRSS(x,U_est(1:5,:),U_est,patternModel{mm});
                    [modelTheta,rss] = fminsearch(fcn, modelTheta0{mm}, optimset('MaxIter',20000));
                    Y_hat(testIdx,:,mm) = pp1_encoding('predictModelPatterns',U_est(1:5,:),patternModel{mm},modelTheta);
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
            
            for mm=1:numPart
                trainIdx = ~ismember(pV{s},partI{ii}); 
                M{end+1}.type       = 'component';
                M{end}.numGparams   = 1;
                M{end}.name         = sprintf('lower noise ceiling %d',mm);
                M{end}.fitAlgorithm = 'NR';
                M{end}.Gc           = pcm_estGCrossval(Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx));
            end
            
            [T_cv,D_cv,theta_hat,theta0] = pcm_fitModelIndividCrossval(Y{s},M,pV{s},cV{s},...
                                        'runEffect','random','verbose',0,...
                                        'evaluation',{'R','Rpool','R2','likelihood'});
            % arrange data into output structure:
            cv_idx = logical([zeros(numPart,5),eye(numPart)]);
            v = ones(numel(M)-numPart+1,1);
            d = [];
            d.model  = [1:numel(v)]';
            d.roi    = v.*roi;
            d.sn     = v.*sn(s);
            % get evaluation metrics
            for mm=1:numel(M)-numPart+1
                d.like_cv(mm,1)  = T_cv.likelihood(mm);
                d.r2_cv(mm,1)    = T_cv.R2(mm);
                d.r_cv(mm,1)     = T_cv.R(mm);
                d.rpool_cv(mm,1) = T_cv.Rpool(mm);
                if mm>(numel(M)-numPart)
                    d.like_cv(mm,1)  = sum(D_cv.likelihood(cv_idx));
                    d.r2_cv(mm,1)    = 1 - sum(D_cv.RSS(cv_idx))/sum(D_cv.TSS(cv_idx));
                    d.r_cv(mm,1)     = sum(D_cv.SSC(cv_idx)) ./ sqrt(sum(D_cv.SS1(cv_idx))*sum(D_cv.SS2(cv_idx)));
                    d.rpool_cv(mm,1) = mean(D_cv.SSC(cv_idx)./sqrt(D_cv.SS1(cv_idx).*D_cv.SS2(cv_idx)));
                end
            end
            D=addstruct(D,d);

        end
        
        varargout = {D};
    case 'ENC:estimateU_depreciated'
        % use PCM to estimate the true Us for all chords using data from
        % all runs
        
        % inputs
        Y   = varargin{1}; % matrix of activity patterns [#conds*#runs x #vox]
        pV  = varargin{2}; % partition vector (assume cV and pV are same across subjs)
        cV  = varargin{3}; % condition vector (chord #s)
        % create PCM model structure M:
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
    case 'ENC:plotModelFits'
        % case to plot output from 'ENC:do' for a single roi
        % Currently, we plot the pearson model correlation R, normalized to the
        % correlation fit of the upper noise ceiling.
        % Assumes last model is the upper noise ceiling unless otherwise
        % indicated
        D=varargin{1};
        Lnc = 4; % model # of lower noise ceiling
        Unc = 5; % model # of upper noise ceiling
        vararginoptions(varargin{2:end},{'Lnc','Unc'});
        roi = unique(D.roi);
        if numel(roi)>1; error('more than 1 roi in structure'); end
        modelNames = unique(D.modelName,'stable');
        numModels  = numel(modelNames);
        % plot predicted second moments, and what is not captured by each
        % model:
        g_unc = mean(D.g_pred(D.model==Unc,:).*D.th_scale(D.model==Unc),1);
        for ii=1:numModels
            modelIdx = D.model==ii;
            subplot(2,numModels+1,ii);
            g_pred = mean(D.g_pred(modelIdx,:).*D.th_scale(modelIdx),1);
            imagesc(rsa_squareIPM(g_pred));
            title(sprintf('Gpred : %s',modelNames{ii}));
            
            subplot(2,numModels+1,ii+numModels+1);
            imagesc(rsa_squareIPM(g_unc - g_pred)); % positive plotted values reflect what is under-predicted by model (& vice versa)
            title('Gtrue - Gpred');
        end
        % plot corrected correlation fits:
        subplot(2,numModels+1,numModels+1);
        D.r_corrected = D.r_pcm ./ kron(D.r_pcm(D.model==Unc),ones(numModels,1));
        %barplot(D.model,D.r_corrected);
        plt.box(D.model,D.r_corrected);
        ylabel('corrected Pearson''s r');
        title('model fits');
        ylim([0 1]);
        set(gca,'xticklabel',modelNames,'xticklabelrotation',45);
        % plot raw correlation fits:
        subplot(2,numModels+1,2*(numModels+1));
        %barplot(D.model,D.r_cv);
        plt.box(D.model,D.r_pcm);
        ylabel('Pearson''s r');
        title('model fits');
        set(gca,'xticklabel',modelNames,'xticklabelrotation',45);
        varargout = {D};
    
    case 'ENC:test'
        % case to calculate model selection accuracies using different
        % evaluation metrics
        % MetricType | Name
        % 1 : pearson's R
        % 2 : pooled pearson's R
        % 3 : R2
        % 4 : likelihood
        sn = varargin{1};
        % get random subject data to use:
        roi = 2;
        glm = 4;
        % Get model params from this subject given roi==2
        % There params are the estimated params from the actual model fits,
        % avgered across crossvalidation folds.
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
        
                
        snr = [0,0.01,0.1:0.1:1];
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
                [varE(ii,1),varS(ii,1)]=pp1_encoding('ENC:estimateSNR',Y_sim{ii},cV,pV);
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
                        % estimate the patterns using the provided model
                        % thetas: (not testing how good our parameter
                        % selection is- we are testing how well can we
                        % distinguish the models, given the true params!)
                        Y_hat(testIdx,:,mm) = pp1_encoding('predictModelPatterns',U_est(1:5,:),patternModel{mm},modelTheta{mm});
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
        D.accuracy = (D.modelTrue==D.modelWin) ./ numModels; % divide by numModels to ensure each fitted simulated dataset is only counted once
        A = tapply(D,{'snr','metricType','modelTrue'},{'metricValue','nanmean'},{'snrEst','mean'},{'accuracy','sum'});
        A.accuracy = A.accuracy./numSim; % convert accuracy counts to percent accuracies
        % keep accuracies separate for each true model so that we can probe
        % if one model is more easily dismissed than the other.
        varargout = {A,D};
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
            U_true = pp1_encoding('predictModelPatterns',Usf,model,modelTheta); % simulate data under this model
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
    case 'test_nullmodel' % compare fits of noise ceiling data across different pattern estimation approaches
        % case to get correlation fits of noise ceilings using different
        % methods to estimate Us on each crossval fold (ols & ridge):
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % loop through subjects and fit each individually:
        D=[]; % output
        v=ones(2,1); % helper array
        for s=1:numel(sn)
            % get CV folds
            part = unique(pV{s});
            numPart = numel(part);
            partI = {};
            for ip=1:numPart
                partI{ip}=part(ip);
            end  
            % loop through folds, and calculate noise ceiling fits for each
            % fold using 3 different estimation approaches for training Us
            for ii=1:numel(partI)
                % estimate the true condition activity patterns using
                % training data:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii});
                Ytest    = Y{s}(testIdx,:);
                for jj=1:2  % each of the methods to esitmate training Us
                   switch jj
                        case 1 % estimate training patterns using OLS
                            Ytrain = pp1_encoding('estU_ols',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx));
                        case 2 % estimate training patterns using ridge regression with Identity prior
                            modelG = eye(31);
                            Ytrain = pp1_encoding('estU_ridge',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                    end
                    % use training Us to predict patterns:
                    for mm=1:2 % null, noise ceiling models   
                        switch mm
                            case 1 % null model
                                Ypred = pp1_encoding('predictModelPatterns',Ytrain,'null',[]);
                            case 2 % noise ceiling model
                                Ypred = Ytrain;
                        end
                        % evaluate estimated Us (Ypred) against Ytest patterns:
                        [SS1(ii,jj,mm), SS2(ii,jj,mm), SSC(ii,jj,mm), RSS(ii,jj,mm), TSS(ii,jj,mm)] = pp1_encoding('evaluate',Ypred,Ytest);
                    end
                end
            end
            for mm=1:2 % null, noise ceiling
                d.sn  = v.*sn(s);
                d.roi = v.*roi;
                d.glm = v.*glm;
                d.model = v.*mm;
                d.estMethod = [1:jj]';
                d.r   = [mean(SSC(:,:,mm)./sqrt(SS1(:,:,mm).*SS2(:,:,mm)))]';
                d.r2  = [1-sum(RSS(:,:,mm))./sum(TSS(:,:,mm))]';
                D=addstruct(D,d); 
            end
        end
        % plot fits:
        barplot([D.estMethod D.model], D.r);
        ylabel('Pearson''s R');
        set(gca,'xticklabel',{'OLS- null','OLS- nceil','ridge I- null','ridge I- nceil'},'xticklabelrotation',45);
        varargout={D};
    case 'test_models' % compare fits of noise ceiling data across different pattern estimation approaches
        % case to get correlation fits of noise ceilings using different
        % methods to estimate Us on each crossval fold (ols & ridge):
        sn  = 10;%[2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % loop through subjects and fit each individually:
        D=[]; % output
        v=ones(2,1); % helper array
        for s=1:numel(sn)
            % get CV folds
            part = unique(pV{s});
            numPart = numel(part);
            partI = {};
            for ip=1:numPart
                partI{ip}=part(ip);
            end  
            % loop through folds, and calculate noise ceiling fits for each
            % fold using 3 different estimation approaches for training Us
            for ii=1:numel(partI)
                % estimate the true condition activity patterns using
                % training data:
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii});
                Ytest    = Y{s}(testIdx,:);
                for jj=1:2  % each of the methods to esitmate training Us
                   switch jj
                        case 1 % estimate training patterns using OLS
                            Ytrain = pp1_encoding('estU_ols',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx));
                        case 2 % estimate training patterns using ridge regression with Identity prior
                            modelG = eye(31);
                            Ytrain = pp1_encoding('estU_ridge',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                    end
                    % use training Us to predict patterns:
                    for mm=1:4 % null, linear, noise ceiling models   
                        switch mm
                            case 1 % null model
                                Ypred = pp1_encoding('predictModelPatterns',Ytrain,'null',[]);
                                thetaEst = nan;
                            case 2 % linear summation
                                thetaFcn = @(x) modelLossRSS(x,Ytrain(1:5,:),Ytrain,'summation'); % minimize pattern RSS in parameter fitting
                                thetaEst = fminsearch(thetaFcn, [0], optimset('MaxIter',50000));
                                Ypred = pp1_encoding('predictModelPatterns',Ytrain(1:5,:),'summation',thetaEst);
                            case 3 % tanh summation
                                thetaFcn = @(x) modelLossRSS(x,Ytrain(1:5,:),Ytrain,'summation_tanh'); % minimize pattern RSS in parameter fitting
                                thetaEst = fminsearch(thetaFcn, [0 1], optimset('MaxIter',50000));
                                Ypred = pp1_encoding('predictModelPatterns',Ytrain(1:5,:),'summation',thetaEst);
                            case 4 % noise ceiling model
                                Ypred = Ytrain;
                                thetaEst = nan;
                        end
                        % evaluate estimated Us (Ypred) against Ytest patterns:
                        [SS1(ii,jj,mm), SS2(ii,jj,mm), SSC(ii,jj,mm), RSS(ii,jj,mm), TSS(ii,jj,mm)] = pp1_encoding('evaluate',Ypred,Ytest);
                    end
                    modelTheta(ii,jj,mm)=thetaEst;
                end
            end
            for mm=1:3 % null, linear, noise ceiling
                d.sn  = v.*sn(s);
                d.roi = v.*roi;
                d.glm = v.*glm;
                d.model = v.*mm;
                d.estMethod = [1:jj]';
                d.r   = [mean(SSC(:,:,mm)./sqrt(SS1(:,:,mm).*SS2(:,:,mm)))]';
                d.r2  = [1-sum(RSS(:,:,mm))./sum(TSS(:,:,mm))]';
                D=addstruct(D,d); 
            end
        end
        % plot fits:
        barplot([D.estMethod D.model], D.r);
        ylabel('Pearson''s R');
        set(gca,'xticklabel',{'OLS- null','OLS- nceil','ridge I- null','ridge I- nceil'},'xticklabelrotation',45);
        varargout={D};
    
    case 'test&plot'
        barClr = {[0.7 0.7 0.7],[0.3 0.3 0.3],[0 0 0],[1 1 1]};
        labels={'linear','','','tanh','','','nceil','','','null','',''};
        warning off
        roi = [1:4,9];
        ii = 1;
        T = [];
        figure('Color',[1 1 1]);
        for rr=roi
            D=[];
            D=pp1_encoding('ENC:do','roi',rr,'estimateU','ols');
            D2=pp1_encoding('ENC:do','roi',rr,'estimateU','blupI');
            D3=pp1_encoding('ENC:do','roi',rr,'estimateU','blupG');
            D=addstruct(D,D2);
            D=addstruct(D,D3);
            subplot(2,5,ii);
            barplot([D.model D.estU], D.r_test,'facecolor',barClr,'split',D.estU);
            ylabel('Pearson''s R');
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            subplot(2,5,ii+5);
            barplot([D.model D.estU], D.r2_test,'facecolor',barClr,'split',D.estU);
            ylabel('R2');
            set(gca,'xticklabel',labels,'xticklabelrotation',45);
            ii=ii+1;
            T=addstruct(T,D);
        end
        legend({'ols','blup I','blup G'});
        % plot model thetas:
        figure('Color',[1 1 1]);
        T.roiN = T.roi;
        T.roiN(T.roi==9) = 5;
        T1 = getrow(T,T.model==1);
        T2 = getrow(T,T.model==2);
        
        subplot(1,3,1);
        barplot([T1.roiN T1.estU],T1.modelTheta(:,1),'facecolor',barClr,'split',T1.estU);
        set(gca,'xticklabel',{'3a','','','3b','','','1','','','2','','','M1','',''});
        ylabel('theta 1 (baseline)');
        xlabel('region');
        title('linear model');
        
        subplot(1,3,2);
        barplot([T2.roiN T2.estU],T2.modelTheta(:,1),'facecolor',barClr,'split',T2.estU);
        set(gca,'xticklabel',{'3a','','','3b','','','1','','','2','','','M1','',''});
        ylabel('theta 1 (baseline)');
        xlabel('region');
        title('tanh model');
        
        subplot(1,3,3);
        barplot([T2.roiN T2.estU],T2.modelTheta(:,2),'facecolor',barClr,'split',T2.estU);
        set(gca,'xticklabel',{'3a','','','3b','','','1','','','2','','','M1','',''});
        ylabel('theta 2');
        xlabel('region');
        title('tanh model');
        legend({'ols','blup I','blup G'});
        
        % plot avg activity under each model:
        CAT.linecolor = {[0 0 0.7],[0.7 0 0],[0 0 0],[0 0 0.7]};
        CAT.patchcolor = {[0 0 0.7],[0.7 0 0],[0 0 0],[0 0 0.7]};
        CAT.markerfill = {[0 0 0.7],[0.7 0 0],[1 1 1],[0 0 0.7]};
        CAT.linestyle={'-','-','-','-.'};
        CAT.linewidth = 1.5;
        CAT.markersize = 8;
        CAT.markertype = 'o';
        
        figure('Color',[1 1 1]);
        t = getrow(T,T.roi==2);
        subplot(1,3,1);
        traceplot([1:5],t.avgAct,'split',t.model,'subset',t.estU==1,'CAT',CAT); % ols
        title('BA 3b ols');
        xlabel('# digits');
        ylabel('avg. activity');
        
        subplot(1,3,2);
        traceplot([1:5],t.avgAct,'split',t.model,'subset',t.estU==2,'CAT',CAT); % blup I
        title('BA 3b blup I');
        xlabel('# digits');
        ylabel('avg. activity');
        
        subplot(1,3,3);
        traceplot([1:5],t.avgAct,'split',t.model,'subset',t.estU==3,'CAT',CAT); % blup G
        title('BA 3b blup G');
        xlabel('# digits');
        ylabel('avg. activity');
        legend({'linear','tanh','nceil','null'});
        
        warning on
        varargout={T};
    case 'check_modelPatterns'
        % case to predict patterns under null and linear models, plot
        % patterns for visual inspection against true patterns:
        pltIdx1={[1,2],[4,5],[7,8]}; % subplot indicies
        pltIdx2={[3],[6],[9]};
        modelNames = {'null','linear','Y test'};
        sn=10;
        roi = 2; % data from Ba 3b
        glm = 4;
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        Uest      = pp1_encoding('estU_ols',Y{1},pV{1},cV{1});
        for mm=1:3 % null, linear, noise ceiling
            switch mm
                case 1 % null model
                    Ypred = pp1_encoding('predictModelPatterns',Uest,'null',[]);
                case 2 % linear model
                    thetaFcn = @(x) modelLossRSS(x,Uest(1:5,:),Uest,'summation'); % minimize pattern RSS in parameter fitting
                    thetaEst = fminsearch(thetaFcn, [0], optimset('MaxIter',50000));
                    Ypred = pp1_encoding('predictModelPatterns',Uest(1:5,:),'summation',thetaEst);
                case 3 % noise ceiling
                    Ypred = Uest;
            end
            % plot predicted patterns:
            subplot(3,3,pltIdx1{mm});
            imagesc(Ypred);
            title(modelNames{mm});
            ylabel('chords');
            xlabel('voxels');
            % scatterplot of predicted vs. actual patterns
            subplot(3,3,pltIdx2{mm});
            scatter(Uest(:),Ypred(:));
            title(modelNames{mm});
            ylabel('Y pred');
            xlabel('Y test');
        end
        keyboard
        
end
end

function rss = modelLossRSS(theta,Ysf_train,Ytrain,modelName)
Ypred = pp1_encoding('predictModelPatterns',Ysf_train,modelName,theta); % predict patterns under perscribed model
Ypred = Ypred-mean(Ypred,1); % rmv voxel means
Ytrain = Ytrain-mean(Ytrain,1);
rss   = sum(sum((Ytrain-Ypred).^2)); % calculate L2 loss (RSS)
end

function rss = modelLossRSS_activity(theta,Ysf_train,Ytrain,modelName)
Ypred = pp1_encoding('predictModelActivity',Ysf_train,modelName,theta); % predict patterns under perscribed model
Ypred = Ypred-mean(Ypred,1); % rmv voxel means
Ytrain = Ytrain-mean(Ytrain,1);
rss   = sum(sum((Ytrain-Ypred).^2)); % calculate L2 loss (RSS)
end


function [G,dGdtheta] = pcmGroup_modelpred_flexible(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model and returns the derivative in respect to the parameters
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
scaleParams = exp(theta(1:4));
% single finger features
A = Model.Ac(:,:,1); % A = pcm_diagonalize(G_singleFinger)
G = A*A'; 
[V,lam_G] = eig(full(G));
dS    = diag(lam_G);
idx   = dS>eps;
OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)';
% activity scaling feature
chords     = pp1_encoding('chords');
M          = chords;
numFingers = sum(M,2);
for i = 1:4
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end
%OM = A*A';          
G  = M*OM*M';  % Second moment matrix

for i = 1:4 % scale params
    dM                    = zeros(size(chords));
    dM(numFingers==i+1,:) = chords(numFingers==i+1,:);
    dGdtheta(:,:,i)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)       = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end
end

function [G,dGdtheta] = pcmGroup_modelpred_flexibleNULL(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model and returns the derivative in respect to the parameters
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
scaleParams = exp(theta(1:4));
% activity scaling feature
chords     = pp1_encoding('chords');
M          = chords;
numFingers = sum(M,2);
for i = 1:4
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end
OM = ones(5);
G  = M*OM*M';  % Second moment matrix

for i = 1:4 % scale params
    dM                    = zeros(size(chords));
    dM(numFingers==i+1,:) = chords(numFingers==i+1,:);
    dGdtheta(:,:,i)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)       = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end
end



function [rss,dLdtheta] = model_linearGradient(theta,Ysf,Y)
% linear model:
% Ypred = X*Y_sf - X*J.*q1 + Q.*q1
%   where X is chord indicator matrix [31x5],
%   Y_sf are single finger patterns [5xP], 
%   J is matrix of ones [5xP],
%   Q is matrix of ones [31xP],
%   q1 is baseline scalar parameter.

% dRSSdq1 = Q'*T1 - T0'*T1 - T1'*T0 + T1'*Q
%   where T0 = X*J,
%   T1 = Y - X*Y_sf - X*J.*q1 + Q.*q
%
% We can simplify T1:
%   T1 = Y - X*Y_sf - X*J.*q1 + Q.*q
%      = Y - Ypred
%      = residuals


Ypred = pp1_encoding('predictModelPatterns',Ysf,'summation',theta);
Q = ones(size(Ypred));
J = ones(size(Ysf));
X = pp1_encoding('chords');
T0 = X*J;
res = (Y-mean(Y,2))-(Ypred-mean(Ypred,2)); % residuals after mean pattern removal
dLdtheta = sum(sum( Q'*res - T0'*res - res'*T0 + res'*Q ));
rss = sum(sum(res.^2));
end







% use fminsearch for pattern parameter estimation. Calculating partial
% dervs of parameters wrt RSS (for gradient descent) is tricky here 
% and beyond me for the moment.
% % partial derivatives of tanh model thetas wrt L2 loss.
% %   Ymf = tanh((W*(B-theta1)+theta1).*theta2) / theta2;
% % To simplfy a bit, the linear summation portion of the model will be
% % substitued with X:
% %   X = (W*(B-theta1)+theta1);
% 
% X = (W*(B-theta1)+theta1); % do linear summation
% Ypred = tanh(X.*theta2) / theta2; % apply nonlinearity (squishing)
% res = Y-Ypred;
% L = sum(sum(res.^2)); % L2 loss (RSS)
% dLdtheta1 = -2 * (1-W) * sech(X.*theta2).^2 * res; % par deriv of baseline param wrt RSS- THIS IS WRONG! W is a matrix, not sure how to do this properly
% dLdtheta2 = -2 * (tanh(X.*theta2)./theta2^2 - (X*sech(X.*theta2).^2)./theta2) * res; % par deriv of tanh param wrt RSS
