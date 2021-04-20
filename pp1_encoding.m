function varargout=pp1_encoding(what,varargin)
%% details

% Directories for data
dataDir = '/Users/sarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_analysis/pp1_encodingAnalyses/data';
%dataDir = '/home/saarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_analysis/pp1_encodingAnalyses/data';
%dataDir = '/Users/sarbuckle/DATA/passivePatterns1/fmri/RegionOfInterest'; % where betas are saved
%dataDir = '/Users/jdiedrichsen/Dropbox (Diedrichsenlab)/projects/encodingAnalyses'; % where betas are saved

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
roiPlotClrs= {[69 114 180]./255, [116 173 209]./255, [254 224 144]./255, [253 174 97]./255, [244 109 67]./255, [215 48 39]./255};

%% cases
switch what       
    case 'plot_ldc'
        % plots avg. distance b/t specified conditions for specified rois
        glm = 4;
        sn = 2:11;
        roi = [1:7];
        roiPlotNum = [3 4 5 6 1 2 7];
        conds = 1:5; % get paired distances between these conditions
        % load data:
        D = pp1_encoding('getDists','sn',sn,'roi',roi,'conds',conds,'glm',glm);
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
        lineplot(D.roiPlotNum,D.ldc,'CAT',CAT,'errorfcn','stderr','plotfcn','median');
        drawline(0,'dir','horz','linestyle','-'); 
        xlabel('Brodmann area');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2','S2'},'xticklabelrotation',0);
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
        roiPlotNum = [3 4 5 6 1 2 7];
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
        lineplot(T.roiPlotNum,T.sft_norm,'CAT',CAT,'errorfcn','stderr','plotfcn','mean');
        drawline(0,'dir','horz','linestyle','-'); 
        xlabel('Brodmann area');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2','S2'},'xticklabelrotation',0);
        ylabel('normalized selectivity');
        title('single finger selectivity');
        hold off
        
        varargout = {T};
    case 'plot_tuningCurves'
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiPlotNum = [3 4 5 6 1 2];
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_sfTuningCurves.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        T.roiPlotNum = roiPlotNum(T.roi)';
        % plot styling
        CAT.linecolor   = {[0 0 1],[0 0.7 1],[0.5 0.5 0.5],[0 0 0],[1 0 0],[1 0.7 0]};
        CAT.linecolor = roiPlotClrs;
        CAT.errorcolor  = CAT.linecolor;
        CAT.markerfill  = CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        CAT.linestyle   = {'-','-','-','-.','-','-'};
        CAT.linewidth   = {1.25,1.25,1.25,2,1.25,1.25};
        CAT.markertype  = 'o';
        % plot:
        lineplot([T.digitMax T.digit],T.tuningNorm,'plotfcn','nanmean','CAT',CAT,'split',T.roiPlotNum,'leg',{'4a','4p','3a','3b','1','2'},'errorfcn','stderr','subset',T.digit~=T.digitMax)
        ylabel('normalized activity');
        set(gca,'xticklabel',{'1.2','1.3','1.4','1.5','2.1','2.3','2.4','2.5','3.1','3.2','3.4','3.5','4.1','4.2','4.3','4.5','5.1','5.2','5.3','5.4'});
        xlabel('tuned finger : other finger');
        ylim([0 1]);
        
        varargout = {T};    
    case 'plot_encodingFits'
        % plots the normalized encoding model fits
        
        doStats=1;
        
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiNames = {'3a','3b','1','2','4a','4p','S2'};
        roiPlotNum = [3 4 5 6 1 2 7];
        roiPlotNames = {'4a','4p','3a','3b','1','2','S2'};
        nCeil = 8;
        nNull = 1;
        numModels=8;
        
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_Final.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==nNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1));
        T.roiPlotNum = roiPlotNum(T.roi)';
        
        
        % nice plot
        ch = plt.helper.get_shades(5,'gray');
        sty = style.custom([{ch{1:4}} {[0.3 0 0.7]}]); 
        sty.general.linestyle = {'-','-','-','-','-.'};
        
        plt.line(T.roiPlotNum,T.r_norm,'split',T.model,'plotfcn','mean','style',sty,'subset',T.model>1 & T.model<8 & T.model~=6);
        ylabel(sprintf('normalized model fits\n(Pearson''s R)'));
        xlabel('Brodmann area');
        title('encoding model fits - passive');
        set(gca,'xticklabel',roiPlotNames);
        
        if doStats
            % make pivottables:
            fprintf('MEAN model fits:\n');
            pivottable(T.model,T.roiPlotNum,T.r_norm,'mean','subset',~ismember(T.model,[1,6,8]));
            fprintf('SEM model fits:\n');
            pivottable(T.model,T.roiPlotNum,T.r_norm,'stderr','subset',~ismember(T.model,[1,6,8]));
            fprintf('\n----------------------------------------\n');
            % compare normalized null vs. linear model fits:
            fprintf('normalized null vs. linear model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',T.model<3);
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare normalized linear model fit across regions:
            fprintf('normalized linear model fit across regions:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.roiPlotNum],{'roi'},'subset',T.model==2);
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare linear model fit in BA 3b vs. other SMc regions:
            % (two sided paired t-test, bonferroni corrected alpha=0.01)
            fprintf('linear model fit in BA 3b vs. other SMc regions:\n');
            for rr=1:6
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(T.r_norm(T.roiPlotNum==4 & T.model==2),T.r_norm(T.roiPlotNum==rr & T.model==2),2,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
            % compare normalized linear and 2f int model fits:
            fprintf('normalized linear vs. 2f int model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',T.model>1 & T.model<4);
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare predictive performance gain of 2f vs. linear model in SMc regions vs. BA 3b:
            % (two sided paired t-test, bonferroni corrected alpha=0.01)
            fprintf('predictive performance gain of 2f vs. linear model in SMc regions vs. BA 3b:\n');
            diff3b = T.r_norm(T.roiPlotNum==4 & T.model==3) - T.r_norm(T.roiPlotNum==4 & T.model==2);
            for rr=1:6
                diffFit = T.r_norm(T.roiPlotNum==rr & T.model==3) - T.r_norm(T.roiPlotNum==rr & T.model==2);
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(diffFit,diff3b,2,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
            % compare 2f vs. 3f model fits:
            fprintf('2f vs. 3f model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[3,4]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare 3f vs. 4f model fits:
            fprintf('3f vs. 4f model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[4,5]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare linear vs. flexible model fits:
            fprintf('linear vs. flexible model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[2,7]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare 2f int vs. flexible model fits:
            fprintf('2f int vs. flexible model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[3,7]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare flexible vs. 2f int model fits in each region:
            % (two sided paired t-test, bonferroni corrected alpha=0.01)
            fprintf('flexible vs. 2f int model fits in each region:\n');
            for rr=1:6
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(T.r_norm(T.roiPlotNum==rr & T.model==3),T.r_norm(T.roiPlotNum==rr & T.model==7),1,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
        end
        
        varargout = {T};
    
    case 'plot_encodingFits_neighbours'
        % plots the normalized encoding model fits from the neighbours
        % only/ no neighbours finger pair analyses
        doStats=1;
        
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
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_neighbours.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R (0=null, 1=noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==nNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1));
        T.roiPlotNum = roiPlotNum(T.roi)';
        
        % plot
        ch = plt.helper.get_shades(5,'gray');
        sty = style.custom({[ch{2}],[ch{3}],[1 0 0],[0 0 1]}); 
        sty.general.linestyle = {'-','-','-.','-.'};
        
        plt.line(T.roiPlotNum,T.r_norm,'split',T.model,'plotfcn','mean','style',sty,'subset',ismember(T.model,[2,3,4,5]));
        ylabel(sprintf('normalized model fits\n(Pearson''s R)'));
        xlabel('Brodmann area');
        title('encoding model fits - passive - neighbours');
        set(gca,'xticklabel',roiPlotNames);
        
        if doStats
            
            % make pivottables:
            fprintf('MEAN model fits [2f int, 2f no neighbours, 2f only neighbours]:\n');
            pivottable(T.model,T.roiPlotNum,T.r_norm,'mean','subset',ismember(T.model,[3,4,5]));
            fprintf('SEM model fits:\n');
            pivottable(T.model,T.roiPlotNum,T.r_norm,'stderr','subset',ismember(T.model,[3,4,5]));
            fprintf('\n----------------------------------------\n');
            
%             % evaluate neighbouring interactions
%             fprintf('linear vs. 2f only neighbours model fits:\n');
%             stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[2,5]));
%             stats.eff.p
%             fprintf('\n----------------------------------------\n');
%             % ttest linear vs. neighbouring interactions:
%             fprintf('two-sided paired ttest LINEAR vs. ONLY NEIGHBOURS model fits:\n');
%             for rr=1:6
%                 fprintf('ROI %s\n',roiPlotNames{rr});
%                 ttest(T.r_norm(T.roiPlotNum==rr & T.model==2),T.r_norm(T.roiPlotNum==rr & T.model==5),2,'paired');
%                 fprintf('\n')
%             end
%             fprintf('\n----------------------------------------\n');
%             
%             % evaluate non-neighbouring interactions
%             fprintf('linear vs. no neighbours model fits:\n');
%             stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[2,4]));
%             stats.eff.p
%             fprintf('\n----------------------------------------\n');
%             % ttest linear vs. non-neighbouring interactions:
%             fprintf('two-sided paired ttest LINEAR vs. NO NEIGHBOURS model fits:\n');
%             for rr=1:6
%                 fprintf('ROI %s\n',roiPlotNames{rr});
%                 ttest(T.r_norm(T.roiPlotNum==rr & T.model==2),T.r_norm(T.roiPlotNum==rr & T.model==4),2,'paired');
%                 fprintf('\n')
%             end
%             fprintf('\n----------------------------------------\n');
            
            % ttest all vs. neighbouring interactions:
            fprintf('two-sided paired ttest ONLY NEIGHBOURS vs. ALL model fits:\n');
            for rr=1:6
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(T.r_norm(T.roiPlotNum==rr & T.model==5),T.r_norm(T.roiPlotNum==rr & T.model==3),2,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
            % ttest all vs. non-neighbouring interactions:
            fprintf('two-sided paired ttest NO NEIGHBOURS vs. ALL model fits:\n');
            for rr=1:6
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(T.r_norm(T.roiPlotNum==rr & T.model==4),T.r_norm(T.roiPlotNum==rr & T.model==3),2,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
            
            % some ANOVAs:
            fprintf('compare ONLY NEIGHBOURS with All pairs model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[3,5]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            fprintf('compare NO NEIGHBOURS with All pairs model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[3,4]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            
            % compare neighbouring vs. non-neighbouring model fits:
            fprintf('NO vs. ONLY NEIGHBOURS model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[4,5]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % ttest neighbouring vs. non-neighbouring model fits:
            fprintf('two-sided paired ttest NO vs. ONLY NEIGHBOURs model fits:\n');
            for rr=1:6
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(T.r_norm(T.roiPlotNum==rr & T.model==4),T.r_norm(T.roiPlotNum==rr & T.model==5),2,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
        end
        
        varargout = {T};
    case 'plot_encodingFits_droppedDigit_dumb'
        % plots the normalized encoding model fits using models where we
        % dropped one of the five fingers:
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiNames = {'3a','3b','1','2','4a','4p','S2'};
        roiPlotNum = [3 4 5 6 1 2 7];
        roiPlotNames = {'4a','4p','3a','3b','1','2','S2'};
        nCeil = 7;
        nNull = 1;
        numModels=7; % with only 4 fingers, we drop one model because there is no 5-finger interaction model
        fingerName = {'Thumb','Index','Middle','Fourth','Little'};
        
        % plot style
        ch = plt.helper.get_shades(5,'gray');
        sty = style.custom([{[1 0 0]} ch{1:5}]); 
        sty = style.custom(ch);
        sty.general.linestyle = {'-','-','-','-','-','-'};
        
        % load undropped model data:
        T0 = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_Final.mat',glm)));
        T0 = getrow(T0,ismember(T0.sn,sn) & ismember(T0.roi,roi) & T0.model~=6);
        T0.model(T0.model==7)=6;
        T0.model(T0.model==8)=7;
        T0.droppedDigit=zeros(size(T0.sn));
        
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_droppedDigit.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        T=addstruct(T,T0);
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==nNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1));
        T.roiPlotNum = roiPlotNum(T.roi)';
        
        for ii=1:5
            % plot:
            subplot(1,5,ii);
            plt.line(T.roiPlotNum,T.r_norm,'split',T.model,'plotfcn','mean','errorfcn','stderr','style',sty,'subset',ismember(T.model,[2,3,6]) & T.droppedDigit==ii);
            ylabel(sprintf('normalized model fits\n(Pearson''s R)'));
            xlabel('Brodmann area');
            title(sprintf('dropped %s',fingerName{ii}));
            set(gca,'xticklabel',roiPlotNames);
            anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',T.model>1 & T.model<4 & T.droppedDigit==ii);
        end
        subplot(1,2,1);
        plt.line(T.roiPlotNum,T.r_norm,'split',T.droppedDigit,'style',sty,'subset',T.model==2,'plotfcn','mean','errorfcn','stderr');
        title('linear model');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2'});
        xlabel('Brodmann area');
        ylabel('model fit (median)');
        drawline(1,'dir','horz','linestyle',':','color','k');
        drawline(0,'dir','horz','linestyle','-','color','k');
        
        subplot(1,2,2);
        plt.line(T.roiPlotNum,T.r_norm,'split',T.droppedDigit,'style',sty,'subset',T.model==3,'plotfcn','mean','errorfcn','stderr');
        title('2-finger interaction model');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2'});
        xlabel('Brodmann area');
        ylabel('model fit (median)');
        drawline(1,'dir','horz','linestyle',':','color','k');
        drawline(0,'dir','horz','linestyle','-','color','k');
        plt.match('y');
        
        anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum T.droppedDigit],{'model','roi','droppedDigit'},'subset',T.model>1 & T.model<4)
        
        varargout = {T};
    
    case 'plot_selectivityVsModelFit'
        model = 2; % linear model fit
        vararginoptions(varargin,{'model'});
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiPlotNum = [3 4 5 6 1 2 7];
        regPlotNames = {'4a','4p','3a','3b','1','2'};
        
        % load selectivities:
        D = load(fullfile(dataDir,sprintf('glm%d_selectivity_fthres.mat',glm)));
        D = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        % compute normalized selectivity:
        Ddata   = getrow(D,D.isEV==0);
        Drand  = getrow(D,D.isEV==1);
        Dsparse = getrow(D,D.isEV==2);
        Ddata.sft_norm = Ddata.sft - Drand.sft;
        Ddata.sft_norm = Ddata.sft_norm./ (Dsparse.sft - Drand.sft);
        D=Ddata; clear Tgauss Tsparse Tdata
        
        % load model fits:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_Final.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==1),ones(8,1)); % null is model 1
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==8),ones(8,1)); % lower noise ceiling is model 5
        modelPlotName = T.modelName(find(T.model==model,1));
        % ensure rows match between structures:
        D=tapply(D,{'roi','sn'},{'sft_norm','mean'});
        T=tapply(T,{'roi','sn','model'},{'r_norm','mean'});
        % arrange into plot-friendly structure:
        D.modelFit = [];
        for mm=unique(T.model)'
            t=getrow(T,T.model==mm);
            D.modelFit = [D.modelFit,t.r_norm];
        end
        D.roiPlotNum = roiPlotNum(D.roi)';
        % plot:
        cla;
        sty = style.custom(plt.helper.get_shades(numel(roi)+1,'hot'));
        plt.xy(D.sft_norm,1-D.modelFit(:,model),D.roiPlotNum,'split',D.roiPlotNum,'style',sty) % linear model is second model
        
        ylabel(sprintf('%s model fit',modelPlotName{1}));
        xlabel('normalized selectivity');
        title('');
        %legend(regPlotNames);
        
        % correlate selectivities and fits per each region in each
        % participant:
        selectivity = pivottable(D.roi,D.sn,D.sft_norm,'mean'); % row is roi, col is subject
        modelFits   = pivottable(D.roi,D.sn,1-D.modelFit(:,model),'mean');
        r = diag(corr(selectivity,modelFits,'type','Spearman'));
        % compute mean and sem correlations using fisher z transform:
        rz = fisherz(r);
        meanR  = fisherinv(mean(rz));
        lowerB = fisherinv(mean(rz) - 1.96*stderr(rz));
        upperB = fisherinv(mean(rz) + 1.96*stderr(rz));
        fprintf('mean within-participant r=%1.3f [%1.3f - %1.3f]\n',meanR,lowerB,upperB);
        fprintf('ttest fisher-z Rs > 0:\n');
        ttest(rz,[],1,'onesample');
        
        varargout = {D,r};
    case 'plot_overlapVsInt'
        vararginoptions(varargin,{'model'});
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiPlotNum = [3 4 5 6 1 2 7];
        regPlotNames = {'4a','4p','3a','3b','1','2'};
        
        % load selectivities:
        D = load(fullfile(dataDir,sprintf('glm%d_selectivity_fthres.mat',glm)));
        D = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        % compute normalized selectivity:
        Ddata   = getrow(D,D.isEV==0);
        Drand  = getrow(D,D.isEV==1);
        Dsparse = getrow(D,D.isEV==2);
        Ddata.sft_norm = Ddata.sft - Drand.sft;
        Ddata.sft_norm = Ddata.sft_norm./ (Dsparse.sft - Drand.sft);
        D=Ddata; clear Tgauss Tsparse Tdata
        
        % load model fits:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_Final.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==1),ones(8,1)); % null is model 1
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==8),ones(8,1)); % lower noise ceiling is model 5
        % ensure rows match between structures:
        D=tapply(D,{'roi','sn'},{'sft_norm','mean'});
        T=tapply(T,{'roi','sn','model'},{'r_norm','mean'});
        % arrange into plot-friendly structure:
        linFit=getrow(T,T.model==2);
        noiseceilingFit=getrow(T,T.model==8);
        D.intMagnitude = noiseceilingFit.r_norm - linFit.r_norm; %1-fit
        
        D.roiPlotNum = roiPlotNum(D.roi)';
        % plot:
        cla;
        %sty = style.custom(plt.helper.get_shades(numel(roi)+1,'hot'));
        sty = style.custom(roiPlotClrs);
        plt.xy(D.sft_norm,D.intMagnitude,D.roiPlotNum,'split',D.roiPlotNum,'style',sty) % linear model is second model
        
        ylabel('1 - fit_{linear}');
        xlabel('normalized selectivity');
        title('');
        %legend(regPlotNames);
        
        % correlate selectivities and fits per each region in each
        % participant:
        selectivity = pivottable(D.roi,D.sn,D.sft_norm,'mean'); % row is roi, col is subject
        modelFits   = pivottable(D.roi,D.sn,D.intMagnitude,'mean');
        r = diag(corr(selectivity,modelFits,'type','Pearson'));
        % compute mean and sem correlations using fisher z transform:
        rz = fisherz(r);
        meanR  = fisherinv(mean(rz));
        lowerB = fisherinv(mean(rz) - 1.96*stderr(rz));
        upperB = fisherinv(mean(rz) + 1.96*stderr(rz));
        fprintf('mean within-participant r=%1.3f [%1.3f - %1.3f]\n',meanR,lowerB,upperB);
        fprintf('ttest fisher-z Rs > 0:\n');
        ttest(rz,[],1,'onesample');
        
        varargout = {D,r};
    
    case 'corr_selectivityVsFeatures'
        % correlate selectivity profiles with linear model fit, predictive
        % gain of flexible model, and predictive gain of 2f interaction
        % model
        vararginoptions(varargin,{'model'});
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        
        % Get Data:
        % load selectivities:
        D = load(fullfile(dataDir,sprintf('glm%d_selectivity_fthres.mat',glm)));
        D = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        % compute normalized selectivity:
        Ddata   = getrow(D,D.isEV==0);
        Drand  = getrow(D,D.isEV==1);
        Dsparse = getrow(D,D.isEV==2);
        Ddata.sft_norm = Ddata.sft - Drand.sft;
        Ddata.sft_norm = Ddata.sft_norm./ (Dsparse.sft - Drand.sft);
        D=Ddata; clear Tgauss Tsparse Tdata
        
        % load model fits:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_Final.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==1),ones(8,1)); % null is model 1
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==8),ones(8,1)); % lower noise ceiling is model 5
        % ensure rows match between structures:
        D=tapply(D,{'roi','sn'},{'sft_norm','mean'});
        T=tapply(T,{'roi','sn','model'},{'r_norm','mean'});
        % arrange into plot-friendly structure:
        D.modelFit = [];
        for mm=unique(T.model)'
            t=getrow(T,T.model==mm);
            D.modelFit = [D.modelFit,t.r_norm];
        end
        selectivity = pivottable(D.roi,D.sn,D.sft_norm,'mean'); % row is roi, col is subject
        clear T t
        % calculate correlations under each model:
        T = []; t=[];
        v = ones(numel(sn),1);
        for mm=1:5
            switch mm
                case 1 % linear model
                    modelProfile = pivottable(D.roi,D.sn,D.modelFit(:,2),'mean'); % linear model fit
                case 2 % flexible gain
                    modelProfile = pivottable(D.roi,D.sn,D.modelFit(:,7)-D.modelFit(:,2),'mean'); % [flexible - linear model fit]
                case 3 % interaction gain 
                    modelProfile = pivottable(D.roi,D.sn,D.modelFit(:,3)-D.modelFit(:,7),'mean'); % [2fint - flexible model fit]
                case 4 % lower noise ceiling (selectivities)
                    modelProfile = [];
                    for ii=1:numel(sn)
                        idx = sn~=sn(ii);
                        modelProfile(:,ii) = mean(selectivity(:,idx),2);
                    end
                case 5 % upper noise ceiling (selectivities)
                    modelProfile = mean(selectivity,2)*ones(1,numel(sn));
            end
            t.r = diag(corr(selectivity,modelProfile,'type','Pearson'));
            t.model = v.*mm;
            t.rois  = v*roi;
            t.sn = sn';
            t.glm   = v.*glm;
            T=addstruct(T,t);
        end
        
        varargout = {T};
        
    case 'plot_encodingFits_active'
        % plots the normalized encoding model fits
        
        doStats=0;
        
        glm = 1;
        sn = 1:8;
        roi = [1:6];
        roiNames = {'3a','3b','1','2','4a','4p','S2'};
        roiPlotNum = [3 4 5 6 1 2 7];
        roiPlotNames = {'4a','4p','3a','3b','1','2','S2'};
        nCeil = 8;
        nNull = 1;
        numModels=8;
        
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_active_Final.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R (0=null, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==nNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1));
        T.roiPlotNum = roiPlotNum(T.roi)';
        
        
        % nice plot
        
        ch = plt.helper.get_shades(5,'gray');
        sty = style.custom([ch {[0.3 0 0.7]}]); 
        sty.general.linestyle = {'-','-','-','-','-','-.'};
        
        plt.line(T.roiPlotNum,T.r_norm,'split',T.model,'plotfcn','mean','style',sty,'subset',T.model>1 & T.model<8);
        drawline(0,'dir','horz');
        drawline(1,'dir','horz');
        ylabel(sprintf('normalized model fits\n(Pearson''s R)'));
        xlabel('Brodmann area');
        title('encoding model fits - active');
        set(gca,'xticklabel',roiPlotNames);
        
        if doStats
            % make pivottables:
            fprintf('MEAN model fits:\n');
            pivottable(T.model,T.roiPlotNum,T.r_norm,'mean','subset',~ismember(T.model,[1,6,8]));
            fprintf('SEM model fits:\n');
            pivottable(T.model,T.roiPlotNum,T.r_norm,'stderr','subset',~ismember(T.model,[1,6,8]));
            fprintf('\n----------------------------------------\n');
            % compare unnormalized null vs. linear model fits:
            fprintf('unnormalized null vs. linear model fits:\n');
            stats = anovaMixed(T.r_test,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',T.model<3);
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare normalized linear model fit across regions:
            fprintf('normalized linear model fit across regions:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.roiPlotNum],{'roi'},'subset',T.model==2);
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare linear model fit in BA 3b vs. other SMc regions:
            % (two sided paired t-test, bonferroni corrected alpha=0.01)
            fprintf('linear model fit in BA 3b vs. other SMc regions:\n');
            for rr=1:6
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(T.r_norm(T.roiPlotNum==4 & T.model==2),T.r_norm(T.roiPlotNum==rr & T.model==2),2,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
            % compare normalized linear and 2f int model fits:
            fprintf('normalized linear vs. 2f int model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',T.model>1 & T.model<4);
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare predictive performance gain of 2f vs. linear model in SMc regions vs. BA 3b:
            % (two sided paired t-test, bonferroni corrected alpha=0.01)
            fprintf('predictive performance gain of 2f vs. linear model in SMc regions vs. BA 3b:\n');
            diff3b = T.r_norm(T.roiPlotNum==4 & T.model==3) - T.r_norm(T.roiPlotNum==4 & T.model==2);
            for rr=1:6
                diffFit = T.r_norm(T.roiPlotNum==rr & T.model==3) - T.r_norm(T.roiPlotNum==rr & T.model==2);
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(diffFit,diff3b,2,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
            % compare 2f vs. 3f model fits:
            fprintf('2f vs. 3f model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[3,4]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare 3f vs. 4f model fits:
            fprintf('3f vs. 4f model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[4,5]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare linear vs. flexible model fits:
            fprintf('linear vs. flexible model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[2,7]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare linear vs. flexible model fits:
            fprintf('linear vs. flexible model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[2,7]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare 2f int vs. flexible model fits:
            fprintf('2f int vs. flexible model fits:\n');
            stats = anovaMixed(T.r_norm,T.sn,'within',[T.model T.roiPlotNum],{'model','roi'},'subset',ismember(T.model,[3,7]));
            stats.eff.p
            fprintf('\n----------------------------------------\n');
            % compare flexible vs. 2f int model fits in each region:
            % (two sided paired t-test, bonferroni corrected alpha=0.01)
            fprintf('flexible vs. 2f int model fits in each region:\n');
            for rr=1:6
                fprintf('ROI %s\n',roiPlotNames{rr});
                ttest(T.r_norm(T.roiPlotNum==rr & T.model==3),T.r_norm(T.roiPlotNum==rr & T.model==7),1,'paired');
                fprintf('\n')
            end
            fprintf('\n----------------------------------------\n');
        end
        
        varargout = {T};
    
    case 'plot_pcmFits'
        % plots the normalized encoding model fits
        numModels = 4; % how many models?
        nNull = 1; % which is null?
        nCeil = 4; % which is noise ceiling?
        
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiNames = {'3a','3b','1','2','4a','4p','S2'};
        roiPlotNum = [3 4 5 6 1 2 7];
        roiPlotNames = {'4a','4p','3a','3b','1','2','S2'};
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_pcmFits_passive_Final.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized log likelihoods (0=null, 1=lower noise ceiling)
        T.like_norm_null = T.like_cv - kron(T.like_cv(T.model==nNull),ones(numModels,1)); % null is model 1
        T.like_norm = T.like_norm_null./kron(T.like_norm_null(T.model==nCeil),ones(numModels,1)); % lower noise ceiling is model 4
        T.roiPlotNum = roiPlotNum(T.roi)';
        T.roiPlotNames = roiPlotNames(T.roi)';
        %T=getrow(T,T.model~=4); % drop super flex model
        % plot group avg. line:
        CAT.markertype = 'o';
        CAT.linewidth  = 2;
        CAT.errorwidth = 1.5;
        CAT.linecolor  = {[0 0 0.9],[0.9 0 0],[0.9 0 0.9]};
        CAT.markerfill = CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        lineplot(T.roiPlotNum,T.like_norm,'CAT',CAT,'errorfcn','stderr','split',T.model,'subset',T.model~=nNull & T.model~=nCeil);
%        sty = style.custom({[0 0 0.9],[0.9 0 0],[0.5 0 0.9]});
%        plt.dot(T.roiPlotNum,T.like_norm,'split',T.model,'subset',T.model~=nNull & T.model~=nCeil,'style',sty);
        
        drawline(0,'dir','horz','linestyle','-');
        drawline(1,'dir','horz','linestyle','-');
        legend({'linear','flexible','super flexible'});
        xlabel('Brodmann area');
        set(gca,'xticklabel',{'4a','','4p','','3a','','3b','','1','','2','','S2',''},'xticklabelrotation',0);
        ylabel('pseudo R2 (normalized log like)');
        title('pcm model fits');
        hold off
        
        
        % stats:
        for rr=roi
            d=getrow(T,T.roi==rr);
            fprintf('\n-----ROI %s-----\n',roiNames{rr});
            fprintf('--linear vs 0--\n');
            ttest(d.like_norm(d.model==2),[],2,'onesample');
            [p,~,s]=signtest(d.like_norm(d.model==2));
            fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.like_norm(d.model==2)),s.sign,p);
            
            fprintf('--flex vs 0--\n');
            ttest(d.like_norm(d.model==3),[],2,'onesample');
            [p,~,s]=signtest(d.like_norm(d.model==3));
            fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.like_norm(d.model==3)),s.sign,p);
            
            fprintf('--flex vs linear--\n');
            ttest(d.like_norm(d.model==3),d.like_norm(d.model==2),2,'paired');
            [p,~,s]=signtest(d.like_norm(d.model==3),d.like_norm(d.model==2),'tail','both');
            fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\nmedian (m2)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.like_norm(d.model==3)),median(d.like_norm(d.model==2)),s.sign,p);
            
            fprintf('--ceil vs linear--\n');
            ttest(d.like_norm(d.model==4),d.like_norm(d.model==2),2,'paired');
            [p,~,s]=signtest(d.like_norm(d.model==4),d.like_norm(d.model==2),'tail','both');
            fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\nmedian (m2)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.like_norm(d.model==4)),median(d.like_norm(d.model==2)),s.sign,p);
            
            
            fprintf('--ceil vs flex--\n');
            ttest(d.like_norm(d.model==4),d.like_norm(d.model==3),2,'paired');
            [p,~,s]=signtest(d.like_norm(d.model==4),d.like_norm(d.model==3),'tail','both');
            fprintf('Wilcoxon signed-rank test:\nmedian (m1)=%1.3f\nmedian (m2)=%1.3f\n W = %1.3f\tp = %1.8f\n\n',median(d.like_norm(d.model==4)),median(d.like_norm(d.model==3)),s.sign,p);
            
            fprintf('\n')
        end
        
        varargout = {T};
    case 'plot_pcmFits_active'
        % plots the normalized encoding model fits
        numModels = 5; % how many models?
        nNull = 1; % which is null?
        nCeil = 5; % which is noise ceiling?
        
        sn = 1:8;
        roi = [1:6];
        roiPlotNum = [3 4 5 6 1 2];
        % load data:
        T = load(fullfile(dataDir,'glm1_pcmFits_active_Final.mat'));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized log likelihoods (0=null, 1=lower noise ceiling)
        T.like_norm = T.like_cv - kron(T.like_cv(T.model==nNull),ones(numModels,1)); % null is model 1
        T.like_norm = T.like_norm./kron(T.like_norm(T.model==nCeil),ones(numModels,1)); % lower noise ceiling is model 4
        T.roiPlotNum = roiPlotNum(T.roi)';
        % plot group avg. line:
        CAT.markertype = 'o';
        CAT.linewidth  = 2;
        CAT.errorwidth = 1.5;
        CAT.linecolor  = {[0 0 0.9],[0.9 0 0],[0.9 0 0.9]};
        CAT.markerfill = CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        lineplot(T.roiPlotNum,T.like_norm,'CAT',CAT,'errorfcn','stderr','split',T.model,'subset',T.model~=nNull & T.model~=nCeil);
%         sty = style.custom({[0 0 0.9],[0.9 0 0],[0.5 0 0.9]});
%         plt.dot(T.roiPlotNum,T.like_norm,'split',T.model,'subset',T.model~=nNull & T.model~=nCeil,'style',sty);
        
        drawline(0,'dir','horz','linestyle','-');
        legend({'linear','flexible','superFlex'});
        xlabel('Brodmann area');
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2'},'xticklabelrotation',0);
        ylabel('pseudo R2 (normalized log like)');
        title('pcm model fits - ACTIVE');
        hold off
        
        
        % stats:
        for rr=roi
            d=getrow(T,T.roi==rr);
            fprintf('\nROI %d\n',rr);
            fprintf('flex vs linear\n');
            ttest(d.like_norm(d.model==3),d.like_norm(d.model==2),2,'paired');
            fprintf('ceil vs flex\n');
            ttest(d.like_norm(d.model==5),d.like_norm(d.model==3),2,'paired');
            fprintf('superflex vs flex\n');
            ttest(d.like_norm(d.model==4),d.like_norm(d.model==3),2,'paired');
            fprintf('ceil vs superflex\n');
            ttest(d.like_norm(d.model==5),d.like_norm(d.model==3),2,'paired');
        end
        
        varargout = {T};
    case 'plot_globalNonlinearity'
        % plot fit of flexible model relative to the linear model fit
        glm = 4;
        sn = 2:11;
        roi = [1:6];
        roiNames = {'3a','3b','1','2','4a','4p','S2'};
        roiPlotNum = [3 4 5 6 1 2 7];
        nCeil    = 8;
        nLinear  = 2; % which model is linear summation?
        nFlex    = 7; % which is flexible scaling?
        numModels=8;
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_encodingFits_passive_Final.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % calculate normalized pearson's R relative to linear model (0=linear, 1=lower noise ceiling)
        T.r_norm = T.r_test - kron(T.r_test(T.model==nLinear),ones(numModels,1)); % null is model 1
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1)); % lower noise ceiling is model 5
        T = getrow(T,T.model==nFlex);
        T.roiPlotNum = roiPlotNum(T.roi)';
        % plot
        CAT.linewidth  = 2;
        CAT.errorwidth = 1.5;
        CAT.markertype = 'o';
        CAT.linecolor  = {[0 0 0],[0.9 0 0]};
        CAT.linestyle = {'-.','-.'};
        CAT.markerfill = CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        CAT.patchcolor = {[0 0 0],[0.9 0 0]};
        % plot
        lineplot(T.roiPlotNum,T.r_norm,'CAT',CAT,'errorfcn','stderr','plotfcn','median');
       % plt.dot(T.roiPlotNum,T.r_norm);
        %hold on;
        
        set(gca,'xticklabel',{'4a','4p','3a','3b','1','2','S2'},'xticklabelrotation',0);
        legend({'R_{flex}'});
        xlabel('Brodmann area');
        ylabel('normalized Pearson''s R (relative to linear fit)');        
        ylim([0 1]);
        varargout = {T};
    
    case 'plot_simStability'
        % plot clrs:
        clrs = {[0.0416671 0 0],[0.45833 0 0],[0.875 0 0],[1 0.29167 0],[0.0416671 0 0]};
        CAT.linecolor = clrs;
        CAT.markerfill = CAT.linecolor;
        CAT.markercolor = CAT.linecolor;
        CAT.errorcolor = CAT.linecolor;
        CAT.facecolor = CAT.linecolor;
        CAT.edgecolor = {'k'};
        
        S = load(fullfile(dataDir,'modelFitStability_simulations.mat'));
        subplot(2,1,1); % un-normalized fits
        barplot(S.signal,S.r_test,'CAT',CAT,'split',S.model,'plotfcn','mean','errorfcn','std');
        ylabel(sprintf('cross-validated model fits\n(Pearson''s R)'));
        xlabel('signal strength');
        title('simulated data from BA1');
        subplot(2,1,2); % normalized fits
        lineplot(S.signal,S.r_norm,'CAT',CAT,'split',S.model,'plotfcn','mean','errorfcn','std');
        ylabel('normalized model fits');
        xlabel('signal strength');
        varargout = {S};
        
    case 'encoding_Passive'
        % case to fit models to participant data
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        vararginoptions(varargin,{'roi','sn'});
        numModels = 8;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
        
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
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            modelName = 'null';
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];
                            
                        case 2 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                             
                            theta0 = [0];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_ns',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 4 % 3rd order polynomial (3-finger interactions)
                            modelName = '3finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov3F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'3dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'3dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 5 % 4th order polynomial (4-finger interactions)
                            modelName = '4finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov4F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'4dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'4dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 6 % 5th order polynomial (5-finger interactions, saturated model)
                            modelName = '5finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov5F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'5dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'5dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        
                        case 7 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible_ns',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5))];
                        case 8 % lower noise ceiling
                            modelName = 'noise ceiling';
                            Ypred = Utrain;
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];      
                    end
                    modNames{mm,1} = modelName;

                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
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
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    case 'encoding_Passive_neighbours'
        % case to fit models to participant data
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        vararginoptions(varargin,{'roi','sn'});
        numModels = 7;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
        
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
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            modelName = 'null';
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];
                            
                        case 2 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                             
                            theta0 = [0];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_ns',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 4 % 2-finger pairs, no neighbours
                            modelName = '2finger no neighbours';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov2F_noNeighbours',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns_nn'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns_nn',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 5 % 2-finger pairs, only neighbours
                            modelName = '2finger only neighbours';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov2F_onlyNeighbours',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns_on'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns_on',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
             
                        case 6 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible_ns',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5))];
                        case 7 % lower noise ceiling
                            modelName = 'noise ceiling';
                            Ypred = Utrain;
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];      
                    end
                    modNames{mm,1} = modelName;

                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
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
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    case 'encoding_Passive_force'
        % case to fit models to participant data with models using
        % stimulation forces
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        vararginoptions(varargin,{'roi','sn'});
        numModels = 7;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
        F = pp1_imana('BEHA_forceMatrix');
        % loop through subjects and fit each individually:
        D=[]; % output
        for s=1:numel(sn)
            fprintf('S%02d\n',s);
            % get force matrix for subj:
            Zs = F.forceN(F.sn==sn(s),:);
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
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                
                Fs0 = tapply(F,{'chordNum'},{'forceN','mean(x,1)'},'subset',F.sn==sn(s) & ~ismember(F.run,partI{ii}));
                Zs0 = Fs0.forceN;
                
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            modelName = 'null';
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];
                            
                        case 2 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                             
                            theta0 = [0];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_ns',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible_ns',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5))];
                        case 4 % linear force model
                            modelName = 'linear force';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSFforce',Y{s}(trainIdx,:),pV{s}(trainIdx),Zs(trainIdx,:),Zs0,modelG);
                             
                            theta0 = [0];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_ns',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 5 % flexible force model
                            modelName = 'flexible force';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSFforce',Y{s}(trainIdx,:),pV{s}(trainIdx),Zs(trainIdx,:),Zs0,modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible_ns',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5))];
                        
                        case 6 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 7 % lower noise ceiling
                            modelName = 'noise ceiling';
                            Ypred = Utrain;
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];      
                    end
                    modNames{mm,1} = modelName;

                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
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
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    case 'encoding_Passive_tessels'
        % case to fit models to participant data
        Y = varargin{1};
        pV = varargin{2};
        cV = varargin{3};
        tesselG = varargin{4};
        tesselNum = varargin{5};

        
        numModels = 8;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        
        % loop through tessels and fit each individually:
        D=[]; % output
        for s=1:numel(Y)
            qq = tesselNum(s);
            fprintf('\ntessel # %d',qq);
            if isempty(Y{s})
                fprintf('\n********** NO DATA **********');
                continue % few participants don't have coverage in all tessels (rare)
            end
            % define model G:
            modelG = rsa_squareIPM(tesselG.g(tesselG.roi==qq,:));
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
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                % predict patterns under each model:
                fprintf('\nModel:'); 
                for mm=1:numModels
                    fprintf('%d.',mm); 
                    switch mm
                        case 1 % null summation model
                            modelName = 'null';
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];
                            
                        case 2 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                             
                            theta0 = [0];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_ns',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 4 % 3rd order polynomial (3-finger interactions)
                            modelName = '3finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov3F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'3dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'3dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 5 % 4th order polynomial (4-finger interactions)
                            modelName = '4finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov4F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'4dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'4dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 6 % 5th order polynomial (5-finger interactions, saturated model)
                            modelName = '5finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov5F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'5dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'5dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        
                        case 7 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible_ns',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5))];
                        case 8 % lower noise ceiling
                            modelName = 'noise ceiling';
                            Ypred = Utrain;
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];      
                    end
                    modNames{mm,1} = modelName;

                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
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
            d.tessel= v.*qq;
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    case 'encoding_Active'
        % case to fit models to participant data from active chord dataset
        sn  = [1:8]; 
        roi = 2; % data from Ba 3b
        glm = 1;
        vararginoptions(varargin,{'roi','sn'});
        numModels = 8;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data
        [Y,pV,cV] = pp1_encoding('getData_active','sn',sn,'roi',roi,'glm',glm);
        modelG = pp1_encoding('getRegionG_active','roi',roi,'glm',glm);
        
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
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            modelName = 'null';
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];
                            
                        case 2 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                             
                            theta0 = [0];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_ns',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 4 % 3rd order polynomial (3-finger interactions)
                            modelName = '3finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov3F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'3dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'3dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 5 % 4th order polynomial (4-finger interactions)
                            modelName = '4finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov4F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'4dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'4dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 6 % 5th order polynomial (5-finger interactions, saturated model)
                            modelName = '5finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov5F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'5dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'5dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        
                        case 7 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible_ns',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5))];
                        case 8 % lower noise ceiling
                            modelName = 'noise ceiling';
                            Ypred = Utrain;
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];      
                    end
                    modNames{mm,1} = modelName;

                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
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
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
       
    case 'predictModelPatterns'
        % factorization of encoding models:
        Ysf = varargin{1}; % training patterns (usually single finger U estimates from training data)
        model = varargin{2};
        theta = varargin{3}; % theta(1) = baseline, theta(2:end) = model params
        if numel(varargin)==4
            subjNum=varargin{4};
        end
        
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
            case '2dInt' % linear summation and paired interactions
                % model that includes 2-finger interaction components
                
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs');
                X  = [X1 X2];
                Ymf_hat = theta(2) * (X*(Ysf-theta(1)) + theta(1)); 
            case '3dInt' % linear summation and paired interactions
                % model that includes 2-finger interaction components
                
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs');
                X3 = pp1_encoding('chord_triplets');
                X  = [X1 X2 X3];
                Ymf_hat = theta(2) * (X*(Ysf-theta(1)) + theta(1)); 
            case '4dInt'
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs');
                X3 = pp1_encoding('chord_triplets');
                X4 = pp1_encoding('chord_quads');
                X  = [X1 X2 X3 X4];
                Ymf_hat = theta(2) * (X*(Ysf-theta(1)) + theta(1)); 
            case '5dInt'
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs');
                X3 = pp1_encoding('chord_triplets');
                X4 = pp1_encoding('chord_quads');
                X  = [X1 X2 X3 X4 [zeros(30,1);1]];
                Ymf_hat = theta(2) * (X*(Ysf-theta(1)) + theta(1)); 
            case '0' % models below do not have a scalar parameter, only baseline parameter
            case 'summation_ns'
                % theta(1) = baseline
                
                % get multi-finger design:
                X = pp1_encoding('chords');
                Ymf_hat = X*(Ysf-theta(1)) + theta(1);
            case 'summation_flexible_ns'
                % theta(1) = baseline
                % theta(2:5) = finger combination param (per # fingers in
                % chords for 2:5 digits)
                
                X = pp1_encoding('chords');
                numD = sum(X,2);
                X = X.*[ones(1,5) exp(theta(numD(numD>1)))]'; % flexible scaling per # fingers stimulated (force positive values with exp)
                Ymf_hat = X*(Ysf-theta(1)) + theta(1); 
            case '2dInt_ns' % linear summation and paired interactions
                % model that includes 2-finger interaction components
                
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs');
                X  = [X1 X2];
                Ymf_hat = X*(Ysf-theta(1)) + theta(1); 
            case '3dInt_ns' % linear summation and paired interactions
                % model that includes 2-finger interaction components
                
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs');
                X3 = pp1_encoding('chord_triplets');
                X  = [X1 X2 X3];
                Ymf_hat = X*(Ysf-theta(1)) + theta(1); 
            case '4dInt_ns'
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs');
                X3 = pp1_encoding('chord_triplets');
                X4 = pp1_encoding('chord_quads');
                X  = [X1 X2 X3 X4];
                Ymf_hat = X*(Ysf-theta(1)) + theta(1); 
            case '5dInt_ns'
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs');
                X3 = pp1_encoding('chord_triplets');
                X4 = pp1_encoding('chord_quads');
                X  = [X1 X2 X3 X4 [zeros(30,1);1]];
                Ymf_hat =X*(Ysf-theta(1)) + theta(1); 
            case '2dInt_ns_nn' % linear summation and paired interactions, no immediately neighbouring pairs
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs_noNeighbours');
                X  = [X1 X2];
                Ymf_hat = X*(Ysf-theta(1)) + theta(1); 
            case '2dInt_ns_on' % linear summation and paired interactions, only immediately neighbouring pairs
                X1 = pp1_encoding('chords');
                X2 = pp1_encoding('chord_pairs_onlyNeighbours');
                X  = [X1 X2];
                Ymf_hat = X*(Ysf-theta(1)) + theta(1); 
                
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
        Z0 = pp1_encoding('chords');
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
        c1 = pp1_encoding('chords'); % single finger features
        c2 = pp1_encoding('chord_pairs'); % finger pair features
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
        c1 = pp1_encoding('chords'); % single finger features
        c2 = pp1_encoding('chord_pairs'); % finger pair features
        c3 = pp1_encoding('chord_triplets'); % finger triplet features
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
        c1 = pp1_encoding('chords'); % single finger features
        c2 = pp1_encoding('chord_pairs'); % finger pair features
        c3 = pp1_encoding('chord_triplets'); % finger triplet features
        c4 = pp1_encoding('chord_quads'); % etc...
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
    case 'estU_tikhonov5F' 
        % Ridge regression estimate of patterns
        % - use pcm with fixed model G to estimate signal and noise
        % parameters
        % - employ ridge regression with lambda = noise param / signal
        % param
        % Explicitly model single fingers and all possible interactions

        % inputs
        Y   = varargin{1}; % matrix of activity patterns [#conds*#runs x #vox]
        pV  = varargin{2}; % partition vector (assume cV and pV are same across subjs)
        cV  = varargin{3}; % condition vector (chord #s)
        G   = varargin{4}; % model prior [31x31]
        
        % create finger feature matrix
        c1 = pp1_encoding('chords'); % single finger features
        c2 = pp1_encoding('chord_pairs'); % finger pair features
        c3 = pp1_encoding('chord_triplets'); % finger triplet features
        c4 = pp1_encoding('chord_quads'); % etc...
        Z0 = [c1 c2 c3 c4 [zeros(30,1);1]];
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
    
    case 'estU_tikhonov2F_noNeighbours' 
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
        c1 = pp1_encoding('chords'); % single finger features
        c2 = pp1_encoding('chord_pairs_noNeighbours'); % finger pair features (for non-neighbouring chords)
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
    case 'estU_tikhonov2F_onlyNeighbours' 
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
        c1 = pp1_encoding('chords'); % single finger features
        c2 = pp1_encoding('chord_pairs_onlyNeighbours'); % finger pair features (for non-neighbouring chords)
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
        
    case 'estU_tikhonovSFforce' 
        % Ridge regression estimate of patterns
        % - use pcm with fixed model G to estimate signal and noise
        % parameters
        % - employ ridge regression with lambda = noise param / signal
        % param

        % inputs
        Y   = varargin{1}; % matrix of activity patterns [#conds*#runs x #vox]
        pV  = varargin{2}; % partition vector (assume cV and pV are same across subjs)
        Z   = varargin{3}; % finger force matrix (all runs)
        Z0  = varargin{4}; % finger force matrix (avg. across runs)
        G   = varargin{5}; % model prior
        
        M{1}.type = 'component';
        %M{1}.Gc=G(1:5,1:5); % group average G
        M{1}.Gc = pinv(Z0)*G*pinv(Z0)';
        M{1}.numGparams = 1;
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,Z,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        Usf = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);

        lambda = exp(theta_hat{1}(2))/exp(theta_hat{1}(1)); %lambda is noise/scale
        
        varargout = {Usf,lambda,theta_hat{1}};    
    
        
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
        M{1}.name = 'null';
        M{1}.type       = 'nonlinear';
        M{1}.numGparams = 4;
        M{1}.theta0     = log([0.9 0.8 0.7 0.6])';
        M{1}.modelpred  = @pcmGroup_modelpred_flexibleNULL;
        
%         M{1}.name = 'null_direct';
%         M{1}.type       = 'fixed';
%         M{1}.numGparams = 0;
%         M{1}.theta0     = [];
%         M{1}.Gc  = pp1_encoding('pcm_null',Y,pV,cV,Gcv);
        
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
        
        % 4. SUPER flexible scaling linear model
%         M{4}.name = 'superDuperFlex';
%         M{4}.type       = 'nonlinear';
%         M{4}.numGparams = 75;
%         M{4}.theta0     = [ones(20,1)*log(0.8); ones(30,1)*log(0.6); ones(20,1)*log(0.4); ones(5,1)*log(0.2);];
%         M{4}.Ac         = pcm_diagonalize(Gsf); 
%         M{4}.modelpred  = @pcmGroup_modelpred_superDuperFlex;
       
        % 5. noise ceiling model
        M{4}.name = 'noiseceiling_direct';
        M{4}.type       = 'freedirect';
        M{4}.numGparams = 0;
        M{4}.theta0     = [];
                
        % choose proper way to deal with run effects
        runEffect = 'random'; % mean patterns not removed
        % fit all models
        [T,theta,G_pred]        = pcm_fitModelGroup(Y,M,pV,cV,'runEffect',runEffect,'fitScale',1,'isCheckDeriv',0,'verbose',1); % perform derivative checks (to ensure I've done correct deriv implementation for nonlinear models)
        [Tcv,theta_cv,Gcv_pred] = pcm_fitModelGroupCrossval(Y,M,pV,cV,'runEffect',runEffect,'groupFit',theta,'fitScale',1,'verbose',1,'isCheckDeriv',0);
        
        % arrange into output structure:
        D=[];
        numModels = numel(M);
        for ii=1:numel(sn)
            for mm=1:numModels
                d.roi = roi;
                d.sn = sn(ii);
                d.model = mm;
                d.modelName = {M{mm}.name};
                
                d.noise = T.noise(ii,mm);
                d.noise_cv = Tcv.noise(ii,mm);
                d.scale = T.scale(ii,mm);
                d.scale_cv = Tcv.scale(ii,mm);
                d.run = T.run(ii,mm);
                d.run_cv = Tcv.run(ii,mm);
                if mm<numModels %&& mm>1
                    d.theta_cv = {exp(theta_cv{mm}(:,ii)')};
                else
                    d.theta_cv = {[]};
                end
                d.gpred = rsa_vectorizeIPM(Gcv_pred{mm}(:,:,ii));
                
                d.like = T.likelihood(ii,mm);
                d.like_cv = Tcv.likelihood(ii,mm);
%                 d.like_norm = Tcv.likelihood_norm(ii,mm);
                
                D=addstruct(D,d);
            end
        end

        varargout = {D,Tcv,T,M,theta_cv,Gcv_pred};    
    case 'pcm_null'
        Y=varargin{1};
        pV=varargin{2};
        cV=varargin{3};
        modelG=varargin{4};
        Z = pcm_indicatorMatrix('identity',sum(pp1_encoding('chords'),2));
        Gnull=[];
        for ii=1:numel(Y)
            U=pp1_encoding('estU_tikhonov',Y{ii},pV{ii},cV{ii},modelG);
            Ua=pinv(Z)*U;
            Up=Z*Ua;
            Gnull(:,:,ii)=(Up*Up')./size(Up,2);
        end
        Gnull=mean(Gnull,3);
        varargout = {Gnull};    
    case 'pcm_fitGroup_active'
        % fits pcm models to data from one region
        sn     = 1:8;
        glm    = 1;
        roi    = [];
        vararginoptions(varargin,{'sn','roi'});

        % load subject data:
        [Y,pV,cV] = pp1_encoding('getData_active','sn',sn,'roi',roi,'glm',glm);
        % load group G:
        Gcv = pp1_encoding('getRegionG_active','roi',roi,'glm',glm);
        Gsf = Gcv(1:5,1:5); % group avg. single-finger G
        
        % define group-level models:
        % 1. null model- fits overall scaling of avg. activity
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
        
        % 4. SUPER flexible scaling linear model
        M{4}.name = 'superFlex';
        M{4}.type       = 'nonlinear';
        M{4}.numGparams = 20;
        M{4}.theta0     = kron(ones(1,5),log([0.8 0.6 0.4 0.2]))';
        M{4}.Ac         = pcm_diagonalize(Gsf); 
        M{4}.modelpred  = @pcmGroup_modelpred_superFlex;
       
        % 5. noise ceiling model
        M{5}.name = 'noiseceiling';
        M{5}.type       = 'freedirect';
        M{5}.numGparams = 0;
        M{5}.theta0     = [];
                
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
        numModels = numel(M);
        for ii=1:numel(sn)
            for mm=1:numModels
                d.roi = roi;
                d.sn = sn(ii);
                d.model = mm;
                
                d.noise = T.noise(ii,mm);
                d.noise_cv = Tcv.noise(ii,mm);
                d.scale = T.scale(ii,mm);
                d.scale_cv = Tcv.scale(ii,mm);
                d.run = T.run(ii,mm);
                d.run_cv = Tcv.run(ii,mm);
                if mm<numModels
                    d.theta_cv = {exp(theta_cv{mm}(:,ii)')};
                else
                    d.theta_cv = {[]};
                end
                d.gpred = rsa_vectorizeIPM(Gcv_pred{mm}(:,:,ii));
                
                d.like = T.likelihood(ii,mm);
                d.like_cv = Tcv.likelihood(ii,mm);
%                 d.like_norm = Tcv.likelihood_norm(ii,mm);
                
                D=addstruct(D,d);
            end
        end

        varargout = {D,Tcv,T,M,theta_cv,Gcv_pred};    
    
        
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
    case 'chord_pairs'
        % returns indicator matrix for pairs of fingers use in each config:
        X=pp1_encoding('chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==2,:); % 2 finger pairs
        Xp = zeros(31,10);
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==2;
            Xp(pidx,ii) = 1;
        end
        varargout = {Xp};
    case 'chord_triplets'
        % returns indicator matrix for sets of 3 fingers use in each config:
        X=pp1_encoding('chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==3,:); % 3 finger triplets
        Xp = zeros(31,10);
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==3;
            Xp(pidx,ii) = 1;
        end
        varargout = {Xp};
    case 'chord_quads'
        % returns indicator matrix for set of 4 fingers use in each config:
        X=pp1_encoding('chords');
        numD = sum(X,2);
        X(X==0)=nan;
        pairs = X(numD==4,:);
        Xp = zeros(31,5);
        for ii=1:size(pairs,1) % for each pair, find where in chords it is used
            pidx = sum(X==pairs(ii,:),2)==4;
            Xp(pidx,ii) = 1;
        end
        varargout = {Xp};
    
    case 'chord_pairs_noNeighbours'
        % chord pairs, kicking out immediate neighours
        X=pp1_encoding('chords');
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
    case 'chord_pairs_onlyNeighbours'
        % chord pairs, kicking out non-immediate neighour pairs
        X=pp1_encoding('chords');
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
        roi = [1:9,12,15:21];
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
    case 'makeRegionG_subj'
        % Estimate Region G (G is avg. semi-positive definite crossval G across
        % participants from roi):
        sn  = 2:11;
        glm = 4;
        roi = [1:9];
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
                G_subj = pcm_makePD(G_subj); % make subject G PD
                d.g = rsa_vectorizeIPM(G_subj);
                % calc PCs from subj G:
                [~,L] = eig((G_subj+G_subj')./2);
                l = sort(diag(L),1,'descend')';   % Sort the eigenvalues from most-to-least important
                % get variance explained by each dimension:
                d.varExp = l./sum(l);
                d.varCum = cumsum(d.varExp);
                d.sn = s;
                d.roi = rr;
                d.glm = glm;
                D=addstruct(D,d);
            end
        end   
        save(fullfile(dataDir,sprintf('glm%d_regionG_subj.mat',glm)),'-struct','D');   
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
    
    case 'calcDists'
        % Estimate pattern dissimilarities
        sn  = 2:11;
        glm = 4;
        roi = [1:7,15:21];

        D =[];
        for rr=roi
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
                bb.numDigits = cell2mat(b.numDigits);
                eval(sprintf('bb.beta = cell2mat(b.%s);',betaType));
                bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive stimulation conditions (tt 32 is thumb press)
                % remove run means
                C0 = indicatorMatrix('identity',bb.run);
                bb.beta_noRunMean = bb.beta - C0*pinv(C0)*bb.beta;
                % remove run means per number of fingers (for cosine
                % distances)
                Cc = indicatorMatrix('identity',bb.numDigits+[bb.run-1].*5);
                bb.beta_noChordMean = bb.beta - Cc*pinv(Cc)*bb.beta;
                bb.beta_noChordMean(bb.chord==31,:) = bb.beta_noRunMean(bb.chord==31,:);
                % calc dissimilarities
                d.ldc = rsa.distanceLDC(bb.beta,bb.run,bb.chord);
                d.ldc_noRunMean = rsa.distanceLDC(bb.beta_noRunMean,bb.run,bb.chord);
                d.ldc_noChordMean = rsa.distanceLDC(bb.beta_noChordMean,bb.run,bb.chord);
                
                G_subj = pcm_estGCrossval(bb.beta_noChordMean,bb.run,bb.chord);
                G_subj = pcm_makePD(G_subj);
                d.cos_noChordMean = corr_crossval(G_subj,'reg','abs');
                
                d.psc = mean(b.psc{1}(1:31,:),2)'; % avg. psc across voxels per chord
                
                % indexing fields:
                d.sn = b.sn;
                d.roi = b.roi;
                d.region = b.regType;
                d.hemi = b.regSide;
                D=addstruct(D,d);
            end
        end   
      %  save(fullfile(dataDir,sprintf('glm%d_dissimilarities.mat',glm)),'-struct','D');   
        varargout = {D};
    case 'getDists'
        glm=[];
        sn =[];
        roi=[];
        conds=[]; % get paired distances between these conditions
        vararginoptions(varargin,{'sn','glm','roi','conds'});
        % load data:
        T = load(fullfile(dataDir,sprintf('glm%d_distances.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        D = [];
        d = [];
        for ii = 1:size(T.sn,1)
            RDM_crossnobis = rsa_squareRDM(T.ldc(ii,:)); 
            RDM_crossnobis = RDM_crossnobis(conds,conds);
            
            d.rdm_ldc = rsa_vectorizeRDM(RDM_crossnobis);
            d.ldc = mean(d.rdm_ldc);
            
            if sum(ismember(conds,[1:5]))==5 && length(conds)==5
                c=indicatorMatrix('allpairs',1:5);
                c(c~=0)=1;
                d.ldc_avgPerFinger = d.rdm_ldc*c./4;
            end
%             d.rdm_corrdist = rsa_vectorizeRDM(RDM_corrdist);
%             d.corrdist = mean(d.rdm_corrdist);
            d.psc = mean(T.psc(ii,conds));
            d.pscCond = T.psc(ii,conds);
            d.roi = T.roi(ii);
            d.reg = T.regType(ii);
            d.hem = T.regSide(ii);
            d.sn  = T.sn(ii);
            d.conds = conds;
            D = addstruct(D,d);
        end
        varargout = {D};
         
    case 'getUsageG'
        % usage second moments:
        G = [      1      0.79727      0.78972      0.78523      0.77084
      0.79727            1       0.9299      0.87708      0.83762
      0.78972       0.9299            1      0.93939      0.88286
      0.78523      0.87708      0.93939            1      0.95249
      0.77084      0.83762      0.88286      0.95249            1];
        G_cent=[ 0.22777    -0.034723    -0.062276    -0.069227    -0.061543
    -0.034723      0.10824     0.018142    -0.037138    -0.054524
    -0.062276     0.018142     0.068245    0.0051692     -0.02928
    -0.069227    -0.037138    0.0051692     0.063315      0.03788
    -0.061543    -0.054524     -0.02928      0.03788      0.10747];
        varargout = {G,G_cent};
    case 'getChordUsageG'
        D=load(fullfile(dataDir,'chordUsageG.mat'));
        G=rsa_squareIPM(mean(D.G));
        varargout = {G};
    
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
        

    case '0' % cases to test/compare encoding approaches:    
    case 'test_noiseceiling' % compare fits of noise ceiling data across different pattern estimation approaches
        % case to get correlation fits of noise ceilings using different
        % methods to estimate single finger Us on each crossval fold (ols, ridge, tikhonov):
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        % loop through subjects and fit each individually:
        D=[]; % output
        v=ones(4,1); % helper array
        for s=1:numel(sn)
            % drop multi-finger chords:
            cidx = cV{s}<6;
            Y{s} = Y{s}(cidx,:);
            cV{s}= cV{s}(cidx);
            pV{s}= pV{s}(cidx);
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
                for jj=1:5 % each of the methods to esitmate training Us
                    switch jj
                        case 1 % estimate training patterns using OLS
                            Ypred = pp1_encoding('estU_ols',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx));
                        case 2 % ridge regression 
                            modelG = eye(5);
                            Ypred = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                        case 3 % tikhonov regression w/ subject Gcv
                            modelG = pcm_makePD(pcm_estGCrossval(Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx)));
                            modelG = modelG(1:5,1:5);
                            Ypred = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                        case 4 % tikhonov regression w/ group Gcv 
                            modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
                            modelG = modelG(1:5,1:5);
                            Ypred = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                        case 5 % tikhonov regression w/ usage G
                            T=load('/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/git_repos/pcm_toolbox/recipe_finger/data_recipe_finger7T.mat');
                            modelG = T.Model(2).G_cent;
                            Ypred = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
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
        set(gca,'xticklabel',{'OLS','ridge','Gsubj','Ggroup','Gusage'},'xticklabelrotation',45);
        varargout={D};
    case 'test_noiseceiling_Gpriors'
        % how does the correlation between single finger patterns from
        % training and test data change when the regularization prior that
        % is used to estimate SF patterns changes? 
        
        % How do we estimate 2f interactions if the prior distribution used
        % in regularized regression will lead to reduced single finger pattern
        % fits to test data??
        
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
                sfidx_train = cV{s}<6 & trainIdx==1; % which are only single-finger conds?
                sf2f_idx    = cV{s}<16 & trainIdx==1;
                sfidx_test  = cV{s}<6 & testIdx==1; % which are only single-finger conds?
                Ytest    = Y{s}(sfidx_test,:);
                for jj=1:3 % each of the methods to esitmate training Us
                    switch jj
                        case 1 % tikhonov regression w/ group Gcv 
                            modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
                            Ypred = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            Ypred = Ypred(1:5,:); % take only single finger patterns
                        case 2 % tikhonov regression w/ group Gcv for single fingers
                            modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
                            modelG = modelG(1:5,1:5);
                            Ypred = pp1_encoding('estU_tikhonov',Y{s}(sfidx_train,:),pV{s}(sfidx_train),cV{s}(sfidx_train),modelG);
                        case 3
                            modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
                            Ypred = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
%                         case 3 % tikhonov regression w/ group Gcv for single and 2finger configs
%                             modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
%                             modelG = modelG(1:15,1:15);
%                             Ypred = pp1_encoding('estU_tikhonov',Y{s}(sf2f_idx,:),pV{s}(sf2f_idx),cV{s}(sf2f_idx),modelG);
%                             Ypred = Ypred(1:5,:); % take only single finger patterns
%                         case 4 % tikhonov regression w/ group Gcv for single and eye for 2finger configs
%                             modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
%                             modelG = [modelG(1:5,1:5) zeros(5,10); zeros(10,5) eye(10).*diag(modelG(6:15,6:15))];
%                             Ypred = pp1_encoding('estU_tikhonov',Y{s}(sf2f_idx,:),pV{s}(sf2f_idx),cV{s}(sf2f_idx),modelG);
%                             Ypred = Ypred(1:5,:);
%                     
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
        set(gca,'xticklabel',{'whole G','sf G','sf+2f G','sf+eye2f G'},'xticklabelrotation',45);
        varargout={D};
    case 'test_scaleParam'
        % case to test if model correlation fits are altered by excluding
        % the model scaling parameter
        % case to fit models to participant data
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        vararginoptions(varargin,{'roi','sn'});
        numModels = 12;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        % load subject data (load in once for all subjs to save time):
        [Y,pV,cV] = pp1_encoding('getData','sn',sn,'roi',roi,'glm',glm);
        modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
        
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
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                % predict patterns under each model:
                for mm=1:numModels
                    fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            modelName = 'null';
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan nan];
                            
                        case 2 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                             
                            theta0 = [0 1];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaScale = thetaEst(2);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf, thetaReg] = pp1_encoding('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 4 % 3rd order polynomial (3-finger interactions)
                            modelName = '3finger';
                            [Uf, thetaReg] = pp1_encoding('estU_tikhonov3F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'3dInt'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'3dInt',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 5 % 4th order polynomial (4-finger interactions)
                            modelName = '4finger';
                            [Uf, thetaReg] = pp1_encoding('estU_tikhonov4F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'4dInt'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'4dInt',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 6 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5) thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5)) thetaEst(6)];
                            thetaFlex = thetaEst(2:5);
                        
                        case 7 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                             
                            theta0 = [0 1];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_ns',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaScale = thetaEst(2);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 8 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf, thetaReg] = pp1_encoding('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 9 % 3rd order polynomial (3-finger interactions)
                            modelName = '3finger';
                            [Uf, thetaReg] = pp1_encoding('estU_tikhonov3F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'3dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'3dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 10 % 4th order polynomial (4-finger interactions)
                            modelName = '4finger';
                            [Uf, thetaReg] = pp1_encoding('estU_tikhonov4F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline thetaScale];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'4dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'4dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 11 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible_ns',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5)) nan];
                            thetaFlex = thetaEst(2:5);    
                        
                        case 12 % lower noise ceiling
                            modelName = 'noise ceiling';
                            Ypred = Utrain;
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan nan];      
                    end
                    modNames{mm,1} = modelName;
                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
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
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
    
    case 'test_stability'
        % case to simulate stability of normalization across different signal levels
        signal = [0.1:0.1:1];
        noise  = 1;
        numSim = 1000;
        roi    = 3; % what roi do we use for reference when simulating data?
        G = pp1_encoding('getRegionG','roi',roi,'glm',4);
        D = [];
        for ii=signal
            fprintf('signal %1.2f\n',ii)
            % simulate data
            [Y,pV,cV] = pp1_encoding('simulate_true',ii,noise,G,numSim);
            % fit encoding models:
            d = pp1_encoding('test_encoding',Y,pV,cV,G);
            d.signal = ones(size(d.simNum)).*ii;
            d.noise  = ones(size(d.simNum)).*noise;
            D = addstruct(D,d);
        end
        
        D.r_norm = D.r_test - kron(D.r_test(D.model==1),ones(5,1));
        D.r_norm = D.r_norm./kron(D.r_norm(D.model==5),ones(5,1));
        
        
        sty = style.custom(plt.helper.get_shades(6,'hot'));
        plt.line(D.signal,D.r_norm,'split',D.model,'style',sty);
        
        save('/Users/sarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_paper/pp1_encodingAnalyses/data/modelFitStability_simulations.mat','-struct','D')
        %keyboard
        varargout = {D};
    case 'simulate_flex'
        % case to simulate multi-finger data under the flexible summation model
        signal = varargin{1};
        noise  = varargin{2};
        G      = varargin{3};
        numSim = varargin{4};
        numVox = 100;
        numRun = 10;
        numCond= 31;
        
        % define flexible PCM model:
        M{1}.name = 'flexible';
        M{1}.type       = 'nonlinear';
        M{1}.numGparams = 4;
        M{1}.theta0     = log([0.8 0.6 0.4 0.2])';
        M{1}.Ac         = pcm_diagonalize(G(1:5,1:5)); 
        M{1}.modelpred  = @pcmGroup_modelpred_flexible;
        
        % load subject data:
%         [Y,pV,cV] = pp1_encoding('getData','sn',2:11,'roi',roi,'glm',4);
%         [~,theta]        = pcm_fitModelGroup(Y,M,pV,cV,'runEffect','random','fitScale',1,'isCheckDeriv',0,'verbose',1); % perform derivative checks (to ensure I've done correct deriv implementation for nonlinear models)
%         theta = theta{1}(1:4);
        theta = [-0.29761; -0.52097; -0.72444; -0.92675];
        % simulate data:
        Z = kron(ones(numRun,1),eye(numCond));
        Y = pcm_makeDataset(M{1},theta,'design',Z,'signal',signal,'noise',noise,...
            'numVox',numVox,'numSim',numSim);
        % get indicator vectors
        for ii=1:numSim
            cV{ii} = kron(ones(numRun,1),[1:numCond]');
            pV{ii} = kron([1:numRun]',ones(numCond,1));
        end
        theta = theta';
        varargout = {Y,pV,cV,theta};
    case 'simulate_true'
        % case to simulate multi-finger data under the flexible summation model
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
        
        % load subject data:
%         [Y,pV,cV] = pp1_encoding('getData','sn',2:11,'roi',4,'glm',4);
%         [~,theta] = pcm_fitModelGroup(Y,M,pV,cV,'runEffect','random','fitScale',1,'isCheckDeriv',0,'verbose',1); % perform derivative checks (to ensure I've done correct deriv implementation for nonlinear models)
        % simulate data:
        Z = kron(ones(numRun,1),eye(numCond));
        Y = pcm_makeDataset(M{1},1,'design',Z,'signal',signal,'noise',noise,...
            'numVox',numVox,'numSim',numSim);
        % get indicator vectors
        for ii=1:numSim
            cV{ii} = kron(ones(numRun,1),[1:numCond]');
            pV{ii} = kron([1:numRun]',ones(numCond,1));
        end

        varargout = {Y,pV,cV};
    case 'test_encoding'
        Y = varargin{1};
        pV = varargin{2};
        cV = varargin{3};
        modelG = varargin{4};
        
        numModels = 5;
        chords = pp1_encoding('chords');
        numD_inv = pinv(indicatorMatrix('identity',sum(chords,2)));
        
        % loop through subjects and fit each individually:
        D=[]; % output
        for s=1:numel(Y)
            fprintf('simNum%02d\n',s);
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
                Ytest    = Y{s}(testIdx,:);
                [Utrain,lambda0, thetaReg0] = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                % predict patterns under each model:
                for mm=1:numModels
                    %fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % null summation model
                            modelName = 'null';
                            % model scaling of mean activity, independent
                            % of finger
                            Ypred = pp1_encoding('predictModelPatterns',Utrain,'null',[]);
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan]; 
                        case 2 % linear summation model
                            modelName = 'linear';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                             
                            theta0 = [0];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_ns',thetaEst);
                            thetaBaseline = thetaEst(1);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 3 % 2nd order polynomial (2-finger interactions)
                            modelName = '2finger';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonov2F',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                                                        
                            theta0 = [thetaBaseline];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'2dInt_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'2dInt_ns',thetaEst);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 4 % flexible summation model
                            modelName = 'flexible';
                            [Uf,lambdaReg,thetaReg] = pp1_encoding('estU_tikhonovSF',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                            
                            theta0 = [thetaBaseline log(1/2) log(1/3) log(1/4) log(1/5)];
                            thetaFcn = @(x) modelLossRSS(x,Uf,Utrain,'summation_flexible_ns'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelPatterns',Uf,'summation_flexible_ns',thetaEst);
                            thetaEst = [thetaEst(1) exp(thetaEst(2:5))];
                        case 5 % lower noise ceiling
                            modelName = 'noise ceiling';
                            Ypred = Utrain;
                            thetaReg = thetaReg0;
                            lambdaReg = lambda0;
                            thetaEst = [nan nan nan nan nan];      
                    end
                    modNames{mm,1} = modelName;

                    G = Ypred*Ypred'./ numVox;
                    G_pred(:,:,mm) = G_pred(:,:,mm) + G;
                    modelTheta{mm}(ii,:) = thetaEst;
                    regTheta{mm}(ii,:)   = thetaReg;
                    regLambda(mm,ii)     = lambdaReg;
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
            d.model  = [1:numModels]';
            d.simNum = v.*s;
            
            D=addstruct(D,d);

        end % for each subject
        
        varargout = {D};
     
        
        
        
        
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
                            Ytrain = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
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
                            Ytrain = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
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
    case 'estU_tikhonov2F_FEATURE' % estimate patterns using ridge regression with a model prior
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
        c1 = pp1_encoding('chords');
        c2 = pp1_encoding('chord_pairs'); 
        Z = kron(ones(numel(unique(pV)),1),[c1 c2]);      
        
        
        M{1}.type = 'feature';
        % define single finger features using second moment of single
        % fingers:
        Gsf = G(1:5,1:5);
        Asf = pcm_diagonalize(Gsf);
        M{1}.Ac(:,:,1) = [Asf, zeros(5,10);zeros(10,15)];
        % features for two-finger interactions:
        for ii=1:10
            M{1}.Ac(:,:,ii+1) = zeros(15); %[zeros(5,15);zeros(10,5) eye(10)];
            M{1}.Ac(ii+5,ii+5,ii+1) = 1;
        end
      %  M{1}.Ac(:,:,2) = [zeros(5,15);zeros(10,5) eye(10)];
        M{1}.numGparams = size(M{1}.Ac,3);
        M{1}.theta0 = ones(M{1}.numGparams,1);
        M{1}.fitAlgorithm = 'NR';
        % fit model G to get noise and signal params:
        [T,theta_hat,Gpred,info] = pcm_fitModelIndivid({Y},M,pV,Z,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        U = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);
        lambda = nan;
        varargout = {U,lambda,theta_hat{1}};    
    case 'estU_tikhonov3F_FEATURE' % estimate patterns using ridge regression with a model prior
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
        c1 = pp1_encoding('chords');
        c2 = pp1_encoding('chord_pairs'); 
        c3 = pp1_encoding('chord_triplets'); 
        Z = kron(ones(numel(unique(pV)),1),[c1 c2 c3]);      
        
        
        M{1}.type = 'feature';
        % define single finger features using second moment of single
        % fingers:
        Gsf = G(1:5,1:5);
        Asf = pcm_diagonalize(Gsf);
        M{1}.Ac(:,:,1) = [Asf, zeros(5,20); zeros(20,25)];
        for ii=1:20
            M{1}.Ac(:,:,ii+1) = zeros(25);
            M{1}.Ac(ii+5,ii+5,ii+1) = 1;
        end

        M{1}.numGparams = size(M{1}.Ac,3);
        M{1}.theta0 = ones(M{1}.numGparams,1);
        M{1}.fitAlgorithm = 'NR';
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,Z,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        U = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);
        lambda = nan;
        varargout = {U,lambda,theta_hat{1}};    
    case 'estU_tikhonov4F_FEATURE' % estimate patterns using ridge regression with a model prior
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
        c1 = pp1_encoding('chords');
        c2 = pp1_encoding('chord_pairs'); 
        c3 = pp1_encoding('chord_triplets'); 
        c4 = pp1_encoding('chord_quads'); 
        Z = kron(ones(numel(unique(pV)),1),[c1 c2 c3 c4]);      
        
        
        M{1}.type = 'feature';
        % define single finger features using second moment of single
        % fingers:
        Gsf = G(1:5,1:5);
        Asf = pcm_diagonalize(Gsf);
        M{1}.Ac(:,:,1) = [Asf, zeros(5,25); zeros(25,30)];
        for ii=1:25
            M{1}.Ac(:,:,ii+1) = zeros(30);
            M{1}.Ac(ii+5,ii+5,ii+1) = 1;
        end
        
        M{1}.numGparams = size(M{1}.Ac,3);
        M{1}.theta0 = ones(M{1}.numGparams,1);
        M{1}.fitAlgorithm = 'NR';
        % fit model G to get noise and signal params:
        [~,theta_hat] = pcm_fitModelIndivid({Y},M,pV,Z,'runEffect','none','verbose',0,'fitScale',0);
        % reconstruct true patterns using regularized regression:
        U = pcm_estimateU(M{1},theta_hat{1},Y,Z,[]);
        lambda = nan;
        varargout = {U,lambda,theta_hat{1}};    
        
    case 'testRegularization'
        sn = 2:11;
        glm = 4;
        roi = 2;
        lambdaGrid = [0.01, 0.1:0.1:1, 2:10, 20:10:100];
        % get data
        [Y,pV,cV]=pp1_encoding('getData','roi',roi,'glm',glm,'sn',sn);
        modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);

        % do PCM lambda estimation:
        D = [];
        for s=1:numel(Y)
            % get CV folds
            part = unique(pV{s});
            numPart = numel(part);
            partI = {};
            for ip=1:numPart
                partI{ip}=part(ip);
            end

            % do regularization with gridded lambda values
            SS1=[]; SS2=[]; SSC=[]; RSS=[]; TSS=[];
            for jj=1:numel(lambdaGrid) 
                for ii=1:numel(partI)
                    trainIdx = ~ismember(pV{s},partI{ii}); 
                    testIdx  = ismember(pV{s},partI{ii}); 
                    Ytest    = Y{s}(testIdx,:);
                    Z = pcm_indicatorMatrix('identity',cV{s}(trainIdx)); % feature matrix
                    Ypred = ((Z'*Z + lambdaGrid(jj)*inv(modelG))^-1)*Z'*Y{s}(trainIdx,:); % regularized feature estimation
                    % calculate fit:
                    [SS1(ii,jj), SS2(ii,jj), SSC(ii,jj), RSS(ii,jj), TSS(ii,jj)] = pp1_encoding('evaluate',Ypred,Ytest);
                end
            end
            % avg. fits across cv folds for each strength of regularization:
            v = ones(jj,1);
            d.sn   = v.*sn(s);
            d.glm  = v.*glm;
            d.roi  = v.*roi;
            d.type = v.*1;
            d.typeName = repmat({'direct'},jj,1);
            d.lambda   = lambdaGrid';
            d.r_test   = [mean(SSC./sqrt(SS1.*SS2))]'; % pearson R
            d.r2_test  = [1-sum(RSS)./sum(TSS)]';      % R2
            D=addstruct(D,d);


            % do regularization with PCM estimate of lambda value
            lambdaPCM = [];
            SS1=[]; SS2=[]; SSC=[]; RSS=[]; TSS=[];
            for ii=1:numel(partI)
                trainIdx = ~ismember(pV{s},partI{ii}); 
                testIdx  = ismember(pV{s},partI{ii}); 
                Ytest    = Y{s}(testIdx,:);
                [Ypred,lambdaPCM(ii)] = pp1_encoding('estU_tikhonov',Y{s}(trainIdx,:),pV{s}(trainIdx),cV{s}(trainIdx),modelG);
                % calculate fit:
                [SS1(ii,1), SS2(ii,1), SSC(ii,1), RSS(ii,1), TSS(ii,1)] = pp1_encoding('evaluate',Ypred,Ytest);
            end
            % avg. fits and lambda across cv folds:
            d=[];
            d.sn   = sn(s);
            d.glm  = glm;
            d.roi  = roi;
            d.type = 2; % PCM
            d.typeName = {'pcm'};
            d.lambda   = mean(lambdaPCM);
            d.r_test   = mean(SSC./sqrt(SS1.*SS2)); % pearson R
            d.r2_test  = 1-sum(RSS)./sum(TSS);      % R2
            D=addstruct(D,d);
        end
        
        varargout = {D};
        
    case 'simulate_sf'
        % simulate chord data from model that only has single finger
        % features
        signal = 0.3;
        noise  = 1;
        % get G
%         G = pp1_encoding('getRegionG','roi',2,'glm',4);
%         G = G./mean(diag(G(1:5,1:5)));
        G = pp1_encoding('getUsageG');
        G = G(1:5,1:5); % single finger G
        % single finger feature matrix
        Z = pp1_encoding('chords');
        % model struct:
        M.type = 'fixed';
        M.Gc=Z*G*Z';
        M.numGparams = 0;
        % simulate:
        [Ysim,pV,cV] = pcm_generateData(M,[],'signal',signal,'noise',noise,...
            'numPart',11,'numVox',500,'numSim',100);
        
        % check true vs simulated G:
%         Gcv = pcm_estGCrossval(Ysim{1},pV,cV);
%         subplot(1,2,1); imagesc(M.Gc);
%         subplot(1,2,2); imagesc(Gcv);
        
        varargout = {Ysim,pV,cV};
    case 'simulate_sf2int'
        % simulate chord data from model that has single finger & 2finger
        % interaction features
        signal = 0.3;
        noise  = 1;
        % get G
%         G = pp1_encoding('getRegionG','roi',2,'glm',4);
%         G = G./mean(diag(G(1:5,1:5)));
        G = pp1_encoding('getUsageG');
        c2_var = mean(diag(G(1:5,1:5)))*0.5;
        G  = [G(1:5,1:5) zeros(5,10); zeros(10,5) eye(10).*c2_var];
        % feature matrix
        c1 = pp1_encoding('chords');
        c2 = pp1_encoding('chord_pairs');
        Z = [c1 c2];
        % model struct:
        M.type = 'fixed';
        M.Gc=Z*G*Z';
        M.numGparams = 0;
        % simulate:
        [Ysim,pV,cV] = pcm_generateData(M,[],'signal',signal,'noise',noise,...
            'numPart',11,'numVox',500,'numSim',100);
        
        % check true vs simulated G:
%         Gcv = pcm_estGCrossval(Ysim{1},pV,cV);
%         subplot(1,2,1); imagesc(M.Gc);
%         subplot(1,2,2); imagesc(Gcv);
        
        varargout = {Ysim,pV,cV};
    case 'test_intModel'
        % case to test that 2 finger interaction model works
        G = pp1_encoding('getUsageG');
        [Y1,pV1,cV1] = pp1_encoding('simulate_sf'); % data under single finger features
        [Y2,pV2,cV2] = pp1_encoding('simulate_sf2int'); % data under single finger +interaction features
        Y = [Y1 Y2];
        pV = [{pV1} {pV2}];
        cV = [{cV1} {cV2}];
        model = [ones(numel(Y1),1),ones(numel(Y2),1).*2];
        D=[];
        for ii=1:numel(Y) % per data type (true model)
            d=[];
            d.modelTrue = model(ii);
            [U,theta_hat] = pp1_encoding('estU_tikhonov2F',Y{ii},pV{d.modelTrue},cV{d.modelTrue},G);
            
            d.U = {U};
            d.theta = theta_hat';
            d.g = rsa_vectorizeIPM(U*U'./size(U,2));
            D=addstruct(D,d);
        end
        D.ratio = D.theta(:,1)./D.theta(:,2); % 1 finger theta / 2 finger theta
        
        varargout = {D};
    case 'plot_intEst_2f'
        % calculate 2nd order finger interactions and plot the resulting 
        % pattern estimates per participant for visual inspection
        [Y,pV,cV]=pp1_encoding('getData','sn',2:11,'roi',6,'glm',4);
        G = pp1_encoding('getRegionG','roi',6,'glm',4);
        for ii=1:numel(Y)
            subplot(4,3,ii);
            U{ii}=pp1_encoding('estU_tikhonov2F',Y{ii},pV{ii},cV{ii},G);
            imagesc(U{ii});
            title(subj_name{ii+1});
        end
        varargout = {U};
    case 'plot_intEst_3f'
        % calculate 2nd & 3rd order finger interactions and plot the resulting 
        % pattern estimates per participant for visual inspection
%         [Y,pV,cV]=pp1_encoding('getData','sn',2:11,'roi',4,'glm',4);
%         G = pp1_encoding('getRegionG','roi',4,'glm',4);
        [Y,pV,cV]=pp1_encoding('getData_active','sn',1:8,'roi',2,'glm',1);
        G = pp1_encoding('getRegionG_active','roi',2,'glm',1);
        for ii=1:numel(Y)
            subplot(4,3,ii);
            U{ii}=pp1_encoding('estU_tikhonov3F',Y{ii},pV{ii},cV{ii},G);
            imagesc(U{ii});
            title(subj_name{ii+1});
        end
        varargout = {U};
          
    case 'do_encoding_activity'
        % case to fit models to participant data
        % encoding models are fit on the overall activity across chords
        % with the same # fingers (the marginal activities).
        % Here, we simply calculate fits for 5 values with 3 models:
        % linear, flexible, and lower noise ceiling.
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        roi = 2; % data from Ba 3b
        glm = 4;
        estimateU = 'Ggroup'; 
        vararginoptions(varargin,{'roi','sn','estimateU'});
        numModels = 3; % linear and noise ceiling (flexible nonlinear model)
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
            case 'Identity' % ridge regression 
                modelG = eye(31);
                estFcn = @(x,y,z,g) pp1_encoding('estU_tikhonov',x,y,z,g);
                estMethod = 2;
            case 'Ggroup' % tikhonov regularization with group crossval G distribution
                modelG = pp1_encoding('getRegionG','roi',roi,'glm',glm);
                estFcn = @(x,y,z,g) pp1_encoding('estU_tikhonov',x,y,z,g);
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
%                     fprintf('Model: %d\n',mm); 
                    switch mm
                        case 1 % linear summation model
%                             X = [ones(5,1) [1:5]'];
%                             b = pinv(X)*Utrain;
%                             Ypred = X*b;
                            theta0 = [0 1];
                            thetaFcn = @(x) modelLossRSS_activity(x,Utrain,'summation'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo] = fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelActivity','summation',thetaEst);
                            thetaInt = thetaEst(1);
                            thetaEst = [thetaEst nan nan nan nan];
                        case 2 % flexible summation model
%                             X = [ones(5,1) eye(5)];
%                             b = pinv(X)*Utrain;
%                             Ypred = X*b;
                            theta0 = real(log(Utrain))';
                            thetaFcn = @(x) modelLossRSS_activity(x,Utrain,'summation_flexible'); % minimize pattern RSS in parameter fitting
                            [thetaEst,feval,ef,fitInfo]= fminsearch(thetaFcn, theta0, optimset('MaxIter',50000));
                            Ypred = pp1_encoding('predictModelActivity','summation_flexible',thetaEst);
                            thetaEst = [exp(thetaEst(1:5)) nan];
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
        model = varargin{1};
        theta = varargin{2}; % theta(1) = baseline, theta(2:end) = model params
        theta = theta';
        switch model
            case 'summation'
                % theta(1) = baseline intercept
                % theta(2) = slope
                X = [ones(5,1) [1:5]']; % design matrix
                Ypred = X*theta;
            case 'summation_flexible'               
                % flexible scaling of patterns per # finger stimulated
                % theta(1) = intercept
                % theta(2:5) = finger scalar params, one per # fingers stimulated
                % get multi-finger design:
                
                X = [eye(5)]; % flexible scaling per # fingers stimulated (force positive values with exp)
                Ypred = X*[exp(theta(1:5))]; 
        end
        varargout = {Ypred};
        
end
end

function rss = modelLossRSS(theta,Ysf_train,Ytrain,modelName)
Ypred = pp1_encoding('predictModelPatterns',Ysf_train,modelName,theta); % predict patterns under perscribed model
Ypred = Ypred-mean(Ypred,1); % rmv voxel means
Ytrain = Ytrain-mean(Ytrain,1);
rss   = sum(sum((Ytrain-Ypred).^2)); % calculate L2 loss (RSS)
end

function rss = modelLossRSS_activity(theta,Ytrain,modelName)
Ypred = pp1_encoding('predictModelActivity',modelName,theta); % predict activity under perscribed model
Ypred = Ypred-mean(Ypred,1); % rmv means
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
OM = A*A';
% G = A*A'; 
% [V,lam_G] = eig(full(G));
% dS    = diag(lam_G);
% idx   = dS>eps;
% OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)';
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

function [G,dGdtheta] = pcmGroup_modelpred_superFlex(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model and returns the derivative in respect to the parameters
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
digitParams = exp(theta(1:20));
% single finger features
A = Model.Ac(:,:,1); % A = pcm_diagonalize(G_singleFinger)
OM = A*A';
% G = A*A'; 
% [V,lam_G] = eig(full(G));
% dS    = diag(lam_G);
% idx   = dS>eps;
% OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)';
% activity scaling feature
chords     = pp1_encoding('chords');
X0         = chords;
numFingers = sum(X0,2);
theta_didx = {[1:4],[5:8],[9:12],[13:16],[17:20]};
M = [];
for dd=1:5
    x = X0(:,dd); % this digit's indicator vector across conds
    theta_digit = [1; digitParams(theta_didx{dd})];
    x(x>0) = theta_digit(numFingers(x>0));
    M = [M,x];
end
     
G  = M*OM*M';  % Second moment matrix

t=1;
for i = 1:5 % per finger
    for j = 1:4 % per number of digits in chord (2:5)- single finger conds are not scaled
       dM                    = zeros(size(chords));
       dM(numFingers==j+1,i) = chords(numFingers==j+1,i);
       dGdtheta(:,:,t)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
       dGdtheta(:,:,t)       = dGdtheta(:,:,t)*digitParams(t); % scaled derivative  
       t=t+1;
    end
end
end

function [G,dGdtheta] = pcmGroup_modelpred_superDuperFlex(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model and returns the derivative in respect to the parameters
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
weightParams = exp(theta(1:75));
% single finger features
A = Model.Ac(:,:,1); % A = pcm_diagonalize(G_singleFinger)
OM = A*A';
% G = A*A'; 
% [V,lam_G] = eig(full(G));
% dS    = diag(lam_G);
% idx   = dS>eps;
% OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)';

% fully flexible weighting of fingers in each chord
chords     = pp1_encoding('chords');
M = chords';
idx = find(M);
idx = idx(6:end); % don't change single finger weights
M(idx) = weightParams;
M = M';
G = M*OM*M';  % Second moment matrix

% parameter index matrix (for derivative calculation)
Widx = [zeros(5) chords(6:end,:)'];
Widx(idx) = 1:75;
Widx = Widx';
for i = 1:75 % per weight param
    dM                    = zeros(size(chords));
    dM(Widx==i)           = 1;
    dGdtheta(:,:,i)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)       = dGdtheta(:,:,i)*weightParams(i); % scaled derivative 
end
end

function [G,dGdtheta] = pcmGroup_modelpred_normDistance(theta,Model)
% Predicts G-matrix under the normalization model that applies
% normalization based on single-finger cortical distances.

% Conceptually, local inhibition will occur moreso when to cortical
% patterns (representations) overlap more strongly than when the patterns
% overlap less.

% To estimate the degree of overlap, we use group-averaged single finger
% dissimilarities.

% The normalization applied is a constant factor whose scale is dependent on the cortical
% dissimiarity. Here we estimate one normalization factor

% Harvest appropriate params
scaleParam = exp(theta(1));
% single finger features
A = Model.Ac(:,:,1); % A = pcm_diagonalize(G_singleFinger)
% G = A*A'; 
% [V,lam_G] = eig(full(G));
% dS    = diag(lam_G);
% idx   = dS>eps;
% OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)'; % make semi-pos definite
OM = A*A'; 

% normalization features
chords = pp1_encoding('chords');
numFingers = sum(chords,2);
M      = chords;
distMat= rsa_squareRDM(Model.distances);
maxDist= max(Model.distances)+0.001;
for ii=6:31 
    idx = find(M(ii,:)>0); % which fingers are stimulated?
    for jj = idx % find avg. distance between this stimulated finger and all other stimulated fingers
        otherFingers = idx(idx~=jj);
        d_avg = mean(distMat(jj,otherFingers));
        M(ii,jj)= d_avg;
    end
end
M(6:31,:) = M(6:31,:).*scaleParam;% rescale normalization by the average distance between simultaneously stimulated fingers
G = M*OM*M';  % Second moment matrix

dM                    = zeros(size(chords));
dM(numFingers>1,:) = chords(numFingers>1,:);
dGdtheta(:,:,1)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
dGdtheta(:,:,1)       = dGdtheta(:,:,1)*scaleParam(1); % scaled derivative 
% 
% for i = 1:4 % scale param
%     dM                    = zeros(size(chords));
%     dM(numFingers==i+1,:) = chords(numFingers==i+1,:);
%     dGdtheta(:,:,i)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
%     dGdtheta(:,:,i)       = dGdtheta(:,:,i)*scaleParam(i); % scaled derivative 
% end
end

function [G,dGdtheta] = pcmGroup_modelpred_flexNormDistance(theta,Model)
% Predicts G-matrix under the normalization model that applies
% normalization based on single-finger cortical distances.

% Conceptually, local inhibition will occur moreso when to cortical
% patterns (representations) overlap more strongly than when the patterns
% overlap less.

% To estimate the degree of overlap, we use group-averaged single finger
% dissimilarities.

% The normalization applied is a constant factor whose scale is dependent on the cortical
% dissimiarity. Here we estimate one normalization factor

% Harvest appropriate params
scaleParams = exp(theta(1:4));
% single finger features
A = Model.Ac(:,:,1); % A = pcm_diagonalize(G_singleFinger)
% G = A*A'; 
% [V,lam_G] = eig(full(G));
% dS    = diag(lam_G);
% idx   = dS>eps;
% OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)'; % make semi-pos definite
OM = A*A'; 

% activity scaling feature
chords     = pp1_encoding('chords');
M          = chords;
numFingers = sum(M,2);
for i = 1:4
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end

% distance feature
M2 = chords;
distMat= rsa_squareRDM(Model.distances);
for ii=6:31 
    idx = find(M2(ii,:)>0); % which fingers are stimulated?
    for jj = idx % find avg. distance between this stimulated finger and all other stimulated fingers
        otherFingers = idx(idx~=jj);
        d_avg = mean(distMat(jj,otherFingers));
        M2(ii,jj)= d_avg;
    end
end
M=M.*M2;% rescale normalization by the average distance between simultaneously stimulated fingers
G = M*OM*M';  % Second moment matrix

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


% % F=pp1_imana('getForceData');
% % F=getrow(F,F.stimulated==1);
% % 
% % T=load('/Users/sarbuckle/Dropbox (Diedrichsenlab)/passivePatterns_paper/pp1_encodingAnalyses/glm4_pcmFits_passive_2.mat');
% % 
% % idx = find(pp1_encoding('chords')');
% % idx = idx(6:end); % don't change single finger weights
% % 
% % % model thetas:
% % W(idx) = cell2mat(T.theta_cv(T.model==4 & T.sn==p.sn));
% % f = pivottable(F.chordNum,F.finger,F.peakF_raw,'nanmean','subset',F.stimulated==1 & F.sn==p.sn);
% % s_idx = ~isnan(f);
% % 
% % P=[];
% % for rr=unique(T.roi)'
% %     for ii=1:10
% %         p.sn=ii+1;
% %         p.roi=rr;
% %         W = zeros(5,31);
% %         % model thetas:
% %         W(idx) = cell2mat(T.theta_cv(T.model==4 & T.sn==p.sn & T.roi==rr)); % place in matrix so we know EXACT arrangement
% %         w = W(:,6:end)'; % drop single finger (all zeros)
% %         % finger forces per chord:
% %         f = pivottable(F.chordNum,F.finger,F.peakF_raw,'nanmean','subset',F.stimulated==1 & F.sn==p.sn);
% %         f = f./diag(f(1:5,1:5))'; % normalize forces relative to single finger stimulation
% %         f = f(6:end,:); % drop single fingers
% %         s_idx = ~isnan(f);
% %         p.forces = f(s_idx)';
% %         p.thetas = w(s_idx)';
% %         p.corr   = corr(p.forces',p.thetas');
% %         P=addstruct(P,p);
% %     end
% % end