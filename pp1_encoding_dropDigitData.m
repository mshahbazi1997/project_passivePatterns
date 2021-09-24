function varargout=pp1_encoding_dropDigitData(what,varargin)
%% details

% Encoding analyses done with dataset where we drop all data associated
% with one of the five digits.
digitToDrop = 1; 
keepDigits = true(1,5);
keepDigits(digitToDrop)=false;

% Directories for data
dataDir = '/Users/sarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_analysis/pp1_encodingAnalyses/data';
%dataDir = '/home/saarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_paper/pp1_encodingAnalyses/data';

subj_name = {'pd01','pd02','s01','s02','s03','s04','s05','s06','s07','s08','s09'}; 
% pd01 is missing fieldmaps, so analyses are on pd02-s0e9

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
    case 'plot_encodingFits'
        % plots the normalized encoding model fits
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiNames = {'3a','3b','1','2','4a','4p','S2'};
        roiPlotNum = [3 4 5 6 1 2 7];
        roiPlotNames = {'4a','4p','3a','3b','1','2','S2'};
        nCeil = 7;
        nNull = 1;
        numModels=7;
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_Final.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==nNull),ones(numModels,1)); % null is model 1
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1)); % lower noise ceiling is model 5
        T.roiPlotNum = roiPlotNum(T.roi)';
        
        
        % nice plot
        ch = plt.helper.get_shades(5,'gray');
        sty = style.custom([ch {[0.3 0 0.7]}]); 
        sty.general.linestyle = {'-','-','-','-','-','-.'};
        
        plt.line(T.roiPlotNum,T.r_norm,'split',T.model,'plotfcn','mean','style',sty,'subset',T.model>1 & T.model<8);
        ylabel(sprintf('normalized model fits\n(Pearson''s R)'));
        xlabel('Brodmann area');
        title('encoding model fits - passive');
        set(gca,'xticklabel',roiPlotNames);
        
%         % plot group avg. line:
%         CAT.linewidth  = 2;
%         CAT.errorwidth = 1.5;
%         CAT.markertype = 'none';
%         CAT.linecolor  = {[0 0 0]};
%         CAT.markerfill = CAT.linecolor;
%         CAT.markercolor = CAT.linecolor;
%         CAT.errorcolor = CAT.linecolor;
%         CAT.linestyle = {':','-'};
%      %   lineplot(T.roiPlotNum,T.r_test,'CAT',CAT,'errorfcn','','split',T.model,'subset',T.model~=2 & T.model~=3);
%         hold off
%         CAT.markertype = 'o';
%         CAT.linecolor  = {[0 0 0.9],[0.9 0 0],[0.9 0.6 0],[0.6 0 0.9]};
%         CAT.linestyle = {'-'};
%         CAT.markerfill = CAT.linecolor;
%         CAT.markercolor = CAT.linecolor;
%         CAT.errorcolor = CAT.linecolor;
%         lineplot(T.roiPlotNum,T.r_norm,'CAT',CAT,'plotfcn','median','errorfcn','stderr','split',T.model,'subset',T.model~=nNull & T.model~=nCeil);
%         set(gca,'xticklabel',{'4a','4p','3a','3b','1','2','S2'},'xticklabelrotation',0);
% %         sty = style.custom(CAT.linecolor);
% %         plt.dot(T.roiPlotNum,T.r_norm,'split',T.model,'subset',T.model~=nNull & T.model~=nCeil,'style',sty,'plotall',1);
% %         set(gca,'xticklabel',{'4a','','4p','','3a','','3b','','1','','2','','S2',''},'xticklabelrotation',0);
%        % drawline(0,'dir','horz','linestyle',':');
%        % drawline(1,'dir','horz','linestyle',':');
%         legend({'linear','flexible'});
%         xlabel('Brodmann area');
%         ylabel('normalized Pearson''s R');
%         title('encoding model fits');
%         hold off
        
        
%         % stats:
%         for rr=roi
%             d=getrow(T,T.roi==rr);
%             fprintf('\n-----ROI %s-----\n',roiNames{rr});
%             fprintf('--linear vs 0--\n');
%             ttest(d.r_norm(d.model==2),[],2,'onesample');
%             [p,~,s]=signtest(d.r_norm(d.model==2));
%             fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.r_norm(d.model==2)),s.sign,p);
%             
%             fprintf('--flex vs 0--\n');
%             ttest(d.r_norm(d.model==3),[],2,'onesample');
%             [p,~,s]=signtest(d.r_norm(d.model==3));
%             fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.r_norm(d.model==3)),s.sign,p);
%             
%             fprintf('--flex vs linear--\n');
%             ttest(d.r_norm(d.model==3),d.r_norm(d.model==2),2,'paired');
%             [p,~,s]=signtest(d.r_norm(d.model==3),d.r_norm(d.model==2),'tail','both');
%             fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\nmedian (m2)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.r_norm(d.model==3)),median(d.r_norm(d.model==2)),s.sign,p);
%             
%             fprintf('--ceil vs linear--\n');
%             ttest(d.r_norm(d.model==nCeil),d.r_norm(d.model==2),2,'paired');
%             [p,~,s]=signtest(d.r_norm(d.model==nCeil),d.r_norm(d.model==2),'tail','both');
%             fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\nmedian (m2)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.r_norm(d.model==nCeil)),median(d.r_norm(d.model==2)),s.sign,p);
%             
%             
%             fprintf('--ceil vs flex--\n');
%             ttest(d.r_norm(d.model==nCeil),d.r_norm(d.model==3),2,'paired');
%             [p,~,s]=signtest(d.r_norm(d.model==nCeil),d.r_norm(d.model==3),'tail','both');
%             fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\nmedian (m2)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.r_norm(d.model==nCeil)),median(d.r_norm(d.model==3)),s.sign,p);
%             
%             fprintf('\n')
%         end
        
        varargout = {T};
        
    case 'encoding_Passive'
        % case to fit models to participant data
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        vararginoptions(varargin,{'roi','sn'});
        numModels = 7;
        chords = pp1_encoding_dropDigitData('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding_dropDigitData('getData','sn',sn,'roi',roi,'glm',glm);
        modelG = pp1_encoding_dropDigitData('getRegionG','roi',roi,'glm',glm);
        
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
            G_pred = zeros(15,15,numModels);
            Y_avg = nan([numModels,4,numPart]); % avg. activity per # digits
            Y_avg_cent = Y_avg;
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
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_encoding_dropDigitData('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            modelName = 'null';
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding_dropDigitData('predictModelPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = nan(1,3);
                            
                        case 2 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding_dropDigitData('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            Ypred = pp1_encoding_dropDigitData('predictModelPatterns',Uf,'summation',thetaEst);
                            thetaEst = nan(1,3);
                        case 3 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf, thetaReg] = pp1_encoding_dropDigitData('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            Ypred = pp1_encoding_dropDigitData('predictModelPatterns',Uf,'2dInt',thetaEst);
                            thetaEst = nan(1,3);
                        case 4 % 3rd order polynomial (3-finger interactions)
                            modelName = '3finger';
                            [Uf, thetaReg] = pp1_encoding_dropDigitData('estU_tikhonov3F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            Ypred = pp1_encoding_dropDigitData('predictModelPatterns',Uf,'3dInt',thetaEst);
                            thetaEst = nan(1,3);
                        case 5 % 4th order polynomial (4-finger interactions, saturated model)
                            modelName = '4finger';
                            [Uf, thetaReg] = pp1_encoding_dropDigitData('estU_tikhonov4F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            Ypred = pp1_encoding_dropDigitData('predictModelPatterns',Uf,'4dInt',thetaEst);
                            thetaEst = nan(1,3);
                        
                        case 6 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding_dropDigitData('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [log(1/2) log(1/3) log(1/4)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding_dropDigitData('predictModelPatterns',Uf,'summation_flexible',thetaEst);
                            thetaEst = [exp(thetaEst)];
                        case 7 % lower noise ceiling
                            modelName = 'noise ceiling';
                            Ypred = Utrain;
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = nan(1,3);      
                    end
                    modNames{mm,1} = modelName;
                    % scale model predictions to left-out test data:
%                     scaling = (Ypred(:)'*Ytest(:))/(Ypred(:)'*Ypred(:));
%                     Ypred = Ypred.*scaling;
                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
                    Y_avg_cent(mm,:,ii)  = mean(numD_inv*(Ypred-mean(Ypred,1)),2)'; % avg. activity per # digits
                    Y_avg(mm,:,ii)       = mean(numD_inv*Ypred,2)'; % avg. activity per # digits
                    % calculate metrics for R and R2 of prediction against TRAINING
                    % data:
                    [SS1t(ii,mm),SS2t(ii,mm),SSCt(ii,mm),RSSt(ii,mm),TSSt(ii,mm)] = pp1_encoding_dropDigitData('evaluate',Ypred,Utrain); % corr b/t pred and train patterns
                    % calculate metrics for R and R2 of prediction against TEST data:
                    [SS1(ii,mm), SS2(ii,mm), SSC(ii,mm), RSS(ii,mm), TSS(ii,mm)] = pp1_encoding_dropDigitData('evaluate',Ypred,Ytest); % corr b/t pred and test patterns
                end
            end
            G_pred = G_pred./numel(partI); % avg. model predicted Gs across folds
            d = [];
            % for each model, avg. thetas across folds:
            d.modelName  = modNames;
            d.modelTheta = cell2mat(cellfun(@(x) mean(x,1),modelTheta,'uni',0)');
            d.regTheta   = cellfun(@(x) mean(x,1),regTheta,'uni',0)';
            d.regLambda  = nanmean(regLambda,2);
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
            d.droppedDigit = v.*digitToDrop;
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    
    case 'predictModelPatterns'
        % factorization of encoding models:
        Ysf = varargin{1}; % training patterns (usually single finger U estimates from training data)
        model = varargin{2};
        theta = varargin{3}; % thetas are for flexible scaling model
        
        switch model % NOTE: none of these models have an overall scaling param. Only a baseline param + any model specific params
            case 'null'
                % Model predicts overall scaling of avg. activity with #
                % fingers. Scaling matches true mean scaling in training
                % data.
                % Ysf here are all 31 conditions from the training data
                % Set each condition pattern to the be mean pattern for all
                % chords with the same number of fingers.
                chords = pp1_encoding_dropDigitData('chords'); 
                X = pcm_indicatorMatrix('identity',sum(chords,2)); % which patterns have the same # of fingers?
                Ymf_hat = X*pinv(X)*Ysf;
            case 'summation'
                % theta(1) = baseline
                % theta(2) = overall scaling param
                
                % get multi-finger design:
                X = pp1_encoding_dropDigitData('chords');
                Ymf_hat = X*Ysf;
            case 'summation_flexible'
                % theta(1) = baseline
                % theta(2:4) = finger combination param (per # fingers in
                % chords for 2:4 digits)
                X = pp1_encoding_dropDigitData('chords');
                numD = sum(X,2);
                X = X.*[ones(1,4) exp(theta(numD(numD>1)-1))]'; % flexible scaling per # fingers stimulated (force positive values with exp)
                Ymf_hat = X*Ysf; 
            case '2dInt' % linear summation and paired interactions
                % model that includes 2-finger interaction components
                
                X1 = pp1_encoding_dropDigitData('chords');
                X2 = pp1_encoding_dropDigitData('chord_pairs');
                X  = [X1 X2];
                Ymf_hat = X*Ysf; 
            case '3dInt' % linear summation and paired interactions
                % model that includes 2-finger interaction components
                
                X1 = pp1_encoding_dropDigitData('chords');
                X2 = pp1_encoding_dropDigitData('chord_pairs');
                X3 = pp1_encoding_dropDigitData('chord_triplets');
                X  = [X1 X2 X3];
                Ymf_hat = X*Ysf; 
            case '4dInt'
                X1 = pp1_encoding_dropDigitData('chords');
                X2 = pp1_encoding_dropDigitData('chord_pairs');
                X3 = pp1_encoding_dropDigitData('chord_triplets');
                X4 = pp1_encoding_dropDigitData('chord_quads');
                X  = [X1 X2 X3 X4];
                Ymf_hat = X*Ysf; 
        end
        varargout = {Ymf_hat};
    
    case 'estU_tikhonov' % estimate patterns using ridge regression with a model prior
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
        lambda = exp(theta_hat{1}(2))/exp(theta_hat{1}(1)); %lambda is noise/scale
        % U = ((Z'*Z + lambda*inv(G))^-1)*Z'*Y;
        
        varargout = {U,lambda,theta_hat{1}};           
    case 'estU_tikhonovSF' % estimate patterns using ridge regression with a model prior
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
        Z0 = pp1_encoding_dropDigitData('chords');
        Zsf = kron(ones(numel(unique(pV)),1),Z0);
        
        M{1}.type = 'component';
        %M{1}.Gc=G(1:5,1:5); % group average G
        M{1}.Gc = pinv(Z0)*G*pinv(Z0)';
        M{1}.numGparams = 1;
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,Zsf,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        Usf = pcm_estimateU(M{1},theta_hat{1},Y,Zsf,[]);

        lambda = exp(theta_hat{1}(2))/exp(theta_hat{1}(1)); %lambda is noise/scale
        
        varargout = {Usf,lambda,theta_hat{1}};    
    case 'estU_tikhonov2F' % estimate patterns using ridge regression with a model prior
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
        
        % create finger feature matrix
        c1 = pp1_encoding_dropDigitData('chords'); % single finger features
        c2 = pp1_encoding_dropDigitData('chord_pairs'); % finger pair features
        Z0 = [c1 c2];
        Z = kron(ones(numel(unique(pV)),1),Z0);      
                
        M{1}.type = 'component';
        M{1}.Gc   = pinv(Z0)*G*pinv(Z0)'; % group average G
        M{1}.numGparams = 1;
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,Z,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        U = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);

        lambda = exp(theta_hat{1}(2))/exp(theta_hat{1}(1)); %lambda is noise/scale
        
        varargout = {U,lambda,theta_hat{1}};    
    case 'estU_tikhonov3F' % estimate patterns using ridge regression with a model prior
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
        
        % create finger feature matrix
        c1 = pp1_encoding_dropDigitData('chords'); % single finger features
        c2 = pp1_encoding_dropDigitData('chord_pairs'); % finger pair features
        c3 = pp1_encoding_dropDigitData('chord_triplets'); % finger triplet features
        Z0 = [c1 c2 c3];
        Z = kron(ones(numel(unique(pV)),1),Z0);      
                
        M{1}.type = 'component';
        M{1}.Gc   = pinv(Z0)*G*pinv(Z0)'; % group average G
        M{1}.numGparams = 1;
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,Z,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        U = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);

        lambda = exp(theta_hat{1}(2))/exp(theta_hat{1}(1)); %lambda is noise/scale
        
        varargout = {U,lambda,theta_hat{1}};    
    case 'estU_tikhonov4F' % estimate patterns using ridge regression with a model prior
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
        
        % create finger feature matrix
        c1 = pp1_encoding_dropDigitData('chords'); % single finger features
        c2 = pp1_encoding_dropDigitData('chord_pairs'); % finger pair features
        c3 = pp1_encoding_dropDigitData('chord_triplets'); % finger triplet features
        c4 = pp1_encoding_dropDigitData('chord_quads'); % etc...
        Z0 = [c1 c2 c3 c4];
        Z = kron(ones(numel(unique(pV)),1),Z0);      
                
        M{1}.type = 'component';
        M{1}.Gc   = pinv(Z0)*G*pinv(Z0)'; % group average G
        M{1}.numGparams = 1;
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,Z,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        U = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);

        lambda = exp(theta_hat{1}(2))/exp(theta_hat{1}(1)); %lambda is noise/scale
        
        varargout = {U,lambda,theta_hat{1}};    
        
    case '0' % helper cases:
    case 'get_droppedIdx'
        % return column vector logical [31] that denotes which of the 31
        % chords has the dropped digit in it.
        c = pp1_encoding_dropDigitData('chordsAll');
        idx = c(:,digitToDrop)==1;
        varargout={idx};
    case 'chordsAll'
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
        hasDigit = chords(:,digitToDrop)==1;
        chords = chords(~hasDigit,logical(keepDigits));
        varargout = {chords}; 
    case 'chord_pairs'
        % returns indicator matrix for pairs of fingers use in each config:

        X=pp1_encoding_dropDigitData('chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==2,:); % 2 finger pairs
        Xp = zeros(size(X,1),sum(numD==2));
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==2;
            Xp(pidx,ii) = 1;
        end

        varargout = {Xp};
    case 'chord_triplets'
        % returns indicator matrix for sets of 3 fingers use in each config:
        X=pp1_encoding_dropDigitData('chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==3,:); % 3 finger triplets
        Xp = zeros(size(X,1),sum(numD==3));
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==3;
            Xp(pidx,ii) = 1;
        end
        
        varargout = {Xp};
    case 'chord_quads'
        % returns indicator matrix for set of 4 fingers use in each config:
        X=pp1_encoding_dropDigitData('chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==4,:); % 4 finger quads
        Xp = zeros(size(X,1),sum(numD==4));
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==4;
            Xp(pidx,ii) = 1;
        end
        
        varargout = {Xp};
       
    case 'getData'
        % Get betas for roi from subjects in PCM-friendly format.
        % Run means are NOT removed as this is a feature we want to retain.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        betaType = 'betaW'; % multivariately prewhitened (or betaUW, raw_beta)
        vararginoptions(varargin,{'sn','glm','roi','betaType'});
        dropIdx = pp1_encoding_dropDigitData('get_droppedIdx'); % which chords to drop?
        keepChords = [1:31];
        keepChords = keepChords(~dropIdx);
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
            bb = getrow(bb,ismember(bb.chord,keepChords)); % drop chords with data that included stimulation of the excluded digit
            % put subj data into pcm variables
            Y{ii}         = bb.betas;
            partVec{ii}   = bb.run;
            condVec{ii}   = bb.chord;
        end
        varargout = {Y,partVec,condVec};     
    case 'getRegionG'
        % Get Region G (G is avg. semi-positive definite crossval G across
        % participants from roi). Participants 2:11 are included in
        % estimate
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        
        dropIdx = pp1_encoding_dropDigitData('get_droppedIdx'); % which chords to drop?
        
        D = load(fullfile(dataDir,sprintf('glm%d_regionG.mat',glm)));
        D = getrow(D,D.roi==roi);
        G = rsa_squareIPM(D.g);
        G = G(1:31,1:31);
        G = G(~dropIdx,~dropIdx);
        
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
          
end
end

function rss = modelLossRSS(theta,Ysf_train,Ytrain,modelName)
Ypred = pp1_encoding_dropDigitData('predictModelPatterns',Ysf_train,modelName,theta); % predict patterns under perscribed model
Ypred = Ypred-mean(Ypred,1); % rmv voxel means
Ytrain = Ytrain-mean(Ytrain,1);
rss   = sum(sum((Ytrain-Ypred).^2)); % calculate L2 loss (RSS)
end
