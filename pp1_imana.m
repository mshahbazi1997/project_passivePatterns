function varargout = pp1_imana(what,varargin)
% function    varargout = pp1_imana(what,varargin)
%  
% passive patterns 1 fmri analyses
%
% SArbuckle, Motor Control Group, 2019, UWO
% saarbuckle@gmail.com

%% ------------------------- Directories ----------------------------------
filePrefix = 'pp1_fmri'; % prefix for output files from experiment script
cwd        = cd; % get current directory when called, and return at end of script
% paths to project directories
codeDir         = '/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_passivePatterns';
baseDir         = '/Users/sarbuckle/DATA/passivePatterns1/fmri';   % base directory for analysis
behavDir        = [baseDir '/data'];             % behavioural data directory
dicomDir        = [baseDir '/data_dicom'];       % imgs hot off the scanner
imagingDirRaw   = [baseDir '/imaging_data_raw']; % raw nifti backups 
fieldmapDir     = [baseDir '/fieldmaps/'];       
imagingDir      = [baseDir '/imaging_data'];               
anatomicalDir   = [baseDir '/anatomicals'];  
cerebAnatDir    = [baseDir '/c_anatomicals'];  
scAnatDir       = [baseDir '/sc_anatomicals']; 
freesurferDir   = [baseDir '/surfaceFreesurfer'];  
wbDir           = [baseDir '/surfaceWB'];
atlasDir        = '/Users/sarbuckle/DATA/Atlas_templates';
caretDir        = [baseDir '/surfaceCaret'];     
gpCaretDir      = [caretDir '/fsaverage_sym'];
regDir          = [baseDir '/RegionOfInterest/'];   
pcmDir          = [baseDir '/PCM_models'];
glmDir          = {[baseDir '/glm1'],[baseDir '/glm2'],[baseDir '/glm3']};              

% set default plotting style
style.file(fullfile(codeDir,'pp1_style.m'));
style.use('default');
%% ------------------------- Exp info -------------------------------------
TR_length  = 1.5;  % seconds
numDummys  = 2;  % dummy images at the start of each run (these are discarded)
numTRs     = {[410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410],...
              [410,410,410,410,410,410]};  % total # of images per run (including dummies) per subj
run{1}     = {[1:6],...
             [1:5],...
             [1:5],...
             [1:6],...
             [1:5],...
             [1:5],...
             [1:6],...
             [1:6]};
run{2}     = {[7:11],...
             [6:11],...
             [6:11],...
             [7:11],...
             [6:11],...
             [6:11],...
             [],...
             []};
%% ------------------------- ROI info -------------------------------------
hem        = {'lh','rh'};                                                   % left & right hemi folder names/prefixes
regname    = {'ba3A','ba3B','ba1','ba2','rM1','cM1','S1','M1','SPLa','SPLp','TH'};        % Cortical ROIs, 5 = S1, 6 = M1; Thalamus = 12                                             % roi names, independent of hemisphere    
regSide    = [ones(size(regname)),...                                       % Hemisphere of the roi
                ones(size(regname)).*2];                                    % [1 = left hemi (contra), 2 = right hemi (ipsi)]
regType    = [1:length(regname),...                                         % roi # (within hemisphere)
                1:length(regname)];
numregions = max(regType);                                                  % total number of regions 

% cerebellum rois (don't over interpret these feature labels)
% Region 1 - Left Hand Presses
% Region 2 - Right Hand Presses
% Region 3 - Saccades
% Region 4 - Action Observation
% Region 5 - Divided Attention
% Region 6 - Active Maintenance
% Region 7 - Narrative
% Region 8 - Semantic Knowledge
% Region 9 - Verbal Fluency
% Region 10- Autobiographical Recall
%% ------------------------- FS info --------------------------------------
atlasA    = 'x';                      % freesurfer filename prefix
atlasname = 'fsaverage_sym';          % freesurfer average atlas
hemName   = {'LeftHem','RightHem'};   % freesurfer hemisphere folder names    
%% ------------------------- Subj info ------------------------------------
% The variables in this section must be updated for every new subject.
%       DiconName  :  first portion of the raw dicom filename
%       NiiRawName :  first protion of the nitfi filename (get after 'PREP_4d_nifti')
%       fscanNum   :  series # for corresponding functional runs. Enter in run order
%       anatNum    :  series # for anatomical scans (~208 or so imgs/series)
%       loc_AC     :  location of the anterior commissure. For some reason,
%                      files from this dataset were not recentred prior to
%                      surface reconstruction (even though 'PREP_centre_AC'
%                      was run for each subject). Thus, AC coords are not
%                      [0 0 0] in this dataset.
%
% The values of loc_AC should be acquired manually prior to the preprocessing
%   Open .nii file with MRIcron and manually find AC and read the xyz coordinate values
%           (note: there values are not [0 0 0] in the MNI coordinate)
%   In MRIcron, the values you want are presented at the top in the program
%   bar.
subj_name  = {'pd01','pd02','s01','s02','s03','s04','s05','s06'};
numSess    = [2,2,2,2,2,2,1,1];
anatNum{1} = {[],[40],[],[],[27:31],[25:29],[29:33],[31:35]};
anatNum{2} = {[],[],[29,32],[20:24],[],[],[],[]};
loc_AC     = {[-78  -122 -126],...
              [-80 -128 -133],...
              [-80 -122 -123],...
              [-82 -116 -132],...
              [-81 -120 -135],...
              [-84 -112 -122],...
              [-82 -122 -133],...
              [-81 -138 -137]};
          
DicomName{1}  = {'2019_03_01_PD01.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_12_pd02_sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_26_PP1_S01.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_15_S02_sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_22_S03_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_24_S04_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_30_S05_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_05_02_S06_Sess1.MR.Diedrichsen_PassivePatterns'};
DicomName{2}  = {'2019_04_22_PD01_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_13_pd02_sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_28_PP1_S01_SESS2.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_16_S02_sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_24_S03_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_30_S04_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '',...
                 ''};
             
NiiRawName{1} = {'2019_03_01_PD01',...
                 '2019_03_12_pd02_sess1',...
                 '2019_03_26_PP1_S01',...
                 '2019_04_15_S02_sess1',...
                 '2019_04_22_S03_Sess1',...
                 '2019_04_24_S04_Sess1',...
                 '2019_04_30_S05_Sess1',...
                 '2019_05_02_S06_Sess1'};
NiiRawName{2} = {'2019_04_22_PD01_Sess2',...
                 '2019_03_13_pd02_sess2',...
                 '2019_03_28_PP1_S01_SESS2',...
                 '2019_04_16_S02_sess2',...
                 '2019_04_24_S03_Sess2',...
                 '2019_04_30_S04_Sess2',...
                 '',...
                 ''};
             
fscanNum{1}   = {[21,24,27,30,33,36],...
                 [24,27,30,33,36],...
                 [14,17,20,23,26],...
                 [13,15,17,19,21,23],...
                 [12,14,16,18,20],...
                 [12,18,20,22,24],...
                 [12,14,16,18,20,22],...
                 [14,16,18,20,22,24]};   
fscanNum{2}   = {[12,14,16,18,20],...
                 [12,15,18,21,24,27],...
                 [12,14,16,18,20,22],...
                 [29,31,35,37,39],...
                 [12,14,16,18,20,22],...
                 [14,16,18,20,22,24],...
                 [],...
                 []};  
             
fieldNum{1}   = {[],...
                 [43,44],...
                 [31:34],...
                 [26:29],...
                 [23:26],...
                 [34:37],...
                 [25:28],...
                 [27:30]};                                               
fieldNum{2}   = {[23:26],...
                 [29,30],...
                 [25:28],...
                 [42:45],...
                 [25:28],...
                 [27:30],...
                 [],...
                 []};
             
dataPrefix    = {'r','u','u','u','u','u','u','u'};             

%% ------------------------- ANALYSES -------------------------------------
switch(what)
    case 'LIST_subjs'               
        D = dload(fullfile(baseDir,'subj_info.txt'));
        
        if nargout==0
            fprintf('\nSN\torigSN\tID\t\tAge\tGender\tHandedness\tNum_FMRI_Sessions');
            fprintf('\n--\t------\t--\t\t---\t------\t----------\t-----------------\n');
            for s = unique(D.sn)'
                S = getrow(D,ismember(D.sn,s));
                fprintf('%02d\t%s\t%s\t\t%d\t%d\t%s\t\t%d',S.sn,S.origSN{1},S.ID{1},S.age,S.gender,S.handedness{1},S.fmri_sessions);
                fprintf('\n');
            end
            fprintf('\n');
        else
            D = rmfield(D,{'ID'});
            varargout = {D};
        end
    case 'getSubjs'
        % returns vector of subject numbers for which we have fully
        % analyzed fmri data available.
        D = pp1_imana('LIST_subjs');
        sn = D.sn(D.fmri_sessions==2)';
        varargout = {sn};
    case 'LIST_rois'
        fprintf('\nROI#\tName\tHemipshere');
        fprintf('\n----\t----\t----------\n');
        for r = 1:length(regSide)
            fprintf('%d\t%s\t%s\n',r,regname{r-(regSide(r)-1)*length(regname)},hemName{regSide(r)});
        end
    case 'CHECK_startTimes'         
        % Check that the startimes match TR times
        % example: sc1_sc2_imana('CHECK:run_times',2,1)
        sn = 1;
        vararginoptions(varargin,{'sn'});

        D = dload(fullfile(behavDir,sprintf('%s_%s.dat',filePrefix,subj_name{sn})));
        figure('Color',[1 1 1])
        for r = unique(D.BN)'
            d = getrow(D,D.BN==r);
            subplot(2,6,r);
            plot(d.startTime,d.startTimeMeas,'k','LineWidth',1.5);
            hold on
            plot(d.startTime,d.startTime,'LineStyle',':','Color','r','LineWidth',2);
            hold off
            title(sprintf('run %02d',r));
            xlabel('expected start time (ms)');
            ylabel('measured start time (ms)');
            axis equal
            grid on
            box off
        end
    case 'CHECK_movement'           
        vararginoptions(varargin,{'sn'});
        glm = 1;
        load(fullfile(glmDir{glm},subj_name{sn},'SPM.mat'));
        spm_rwls_resstats(SPM)          
    
    case '0' % ------------ MISC: some aux. things ------------------------
    case 'MISC_scatterplotMDS'                                              % Called by ROI_MDS_overall to plot scaled representational structures.
        Y     = varargin{1};
        split = [];
        label = [];
        vararginoptions(varargin(2:end),{'split','label'});
        
        % style setup
        if numel(unique(split))==5
            color = {[0 0 0] [0.5 0 0] [0.9 0 0] [1 0.6 0] [0 0 0]};
        elseif numel(unique(split))==2
            color = {[0 0 0] [0 0 0.8]};
        end
        CAT.markercolor = color;
        CAT.markerfill  = color;
        CAT.markertype  = 'o';
        CAT.markersize  = 7;
        % do plot
        scatterplot3(Y(1:31,1),Y(1:31,2),Y(1:31,3),'split',split,'label',label,'CAT',CAT);
        % link the single digit chords for clarification
        indx=[1:5 1]';
        line(Y(indx,1),Y(indx,2),Y(indx,3),'color',color{1});

        if (size(Y,1)==32) % thumb response
            hold on;
            plot3(Y(32,1),Y(32,2),Y(32,3),'o','MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]);
            hold off;
        end
         % rest crosshairs
        hold on;
        plot3(0,0,0,'+','MarkerFaceColor',[0.75, 0, 0.75],'MarkerEdgeColor',[0.75, 0, 0.75],'MarkerSize',8);
        hold off;
        axis equal;
        xlabel('eig 1');
        ylabel('eig 2');
        zlabel('eig 3');
        
        %__________________________________________________________________
          
    case '0' % ------------ BEHA: behavioural force data cases. ----------- % These functions are from fdf2/3 so need editing for your paradigm
    case 'BEHA_analyzeTrials'
        % Process force traces for fmri sessions.
        % - loop through subjects
        % - per subject per block, loop through trials and harvest force traces
        % - save output structure per subject separately
        sn = 1;
        vararginoptions(varargin,{'sn'});
        for s = sn
            dataFileIn  = fullfile(behavDir,sprintf('%s_%s.dat',filePrefix,subj_name{sn}));
            dataFileOut = fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}));
            T           = []; % output structure for subject
            D           = dload(dataFileIn);
            runs        = unique(D.BN);
            for r = runs'
                trials  = D.TN(D.BN==r);
                try
                    % deal with missing mov file from s04 run 1
                    MOV = movload(fullfile(behavDir,sprintf('%s_%s_%02d.mov',filePrefix,subj_name{sn},r)));
                catch
                    continue
                end
                for i = trials' 
                    d = getrow(D,D.TN==i & D.BN==r);
                    t = pp1_fmri_trial(MOV{1,i},d,1);
                    if ismember(r,run{1}{sn}) && ~isempty(t)
                        t.sess = ones(5,1).*1;
                        t.sn   = ones(5,1).*sn;
                        T = addstruct(T,t);
                    elseif ismember(r,run{2}{sn}) && ~isempty(t)
                        t.sess = ones(5,1).*2;  
                        t.sn   = ones(5,1).*sn;
                        T = addstruct(T,t);
                    end
                end
            end
            %save(dataFileOut,'T');
        end
    case 'SUBJ_plotStimForce'
        % plots max stimulation force per finger 
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % load data
        load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T
        D = tapply(T,{'run','numDigits','finger','stimulated'},...
            {'peakF_raw','nanmean'},{'peakF_filt','nanmean'},{'time_stimOnset','nanmean'},{'forceStim','mean'});
        % plot
        style.use('5fingers');
        
        subplot(1,3,1);
        plt.line(D.numDigits,D.peakF_filt./D.forceStim,'split',D.finger,'errorfcn','std','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 1.25]);
        plt.labels('# fingers stimulated','% target force',sprintf('stimulated\ntrials'));
        drawline(1,'dir','horz','linestyle',':','linewidth',1.5);
        legend off
        
        subplot(1,3,2);
        plt.line(D.numDigits,D.peakF_raw./D.forceStim,'split',D.finger,'errorfcn','std','subset',D.stimulated==0);
        plt.set('xlim',[0.5 4.5],'ylim',[0 1.25]);
        plt.labels('# fingers stimulated','% target force',sprintf('non-stimulated\ntrials'));
        drawline(1,'dir','horz','linestyle',':','linewidth',1.5);
        plt.legend('east',{'thumb','index','middle','ring','little'});
        legend off
        
        subplot(1,3,3);
        plt.line(D.numDigits,D.time_stimOnset,'split',D.finger,'errorfcn','std','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 1000]);
        plt.labels('# fingers stimulated','stimulus onset (ms)',sprintf('stimulus onset\ntime'));
        plt.legend('east',{'thumb','index','middle','ring','little'});
        
        varargout = {D};
    case 'SUBJ_getDprime'
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % load data
        % arrange data into plotting structure
        load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T
        % calculate percent incorrect trials
        T.numTrials = ones(size(T.sn));
        P = getrow(T,T.falseResp==0);
        P = tapply(P,{'sn','numDigits','chordNum','sess'},{'isError','sum'},{'numTrials','sum'}); % T has 5 rows per trial (one per finger)
        P.isError   = P.isError./5; % correct for number of rows per trial
        P.numTrials = P.numTrials./5;
        P.perErr    = P.isError./P.numTrials;
        falseAlarm  = sum(P.isError)/sum(P.numTrials);
        % calculate percent correct for misleading trials (mismatch
        % chord)
        F = getrow(T,T.falseResp==1);
        F = tapply(F,{'run','trial','sn','numDigits','chordNum','sess','falseResp'},{'isError','sum'});
        F.isError = F.isError./5;
        hitRate   = sum(F.isError==0)/size(F.sn,1);
        % calculate D prime
        D.sn     = sn;
        D.dprime = norminv(hitRate) - norminv(falseAlarm);
        D.perErr = 1-hitRate;
        
        varargout = {D,P,F,T};
    case 'SUBJ_getTaskPerformance'
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % load data
        % arrange data into plotting structure
        load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T
        % calculate percent incorrect trials
        T.numTrials = ones(size(T.sn));
        P = tapply(T,{'sn','numDigits','chordNum','sess'},{'isError','sum'},{'numTrials','sum'}); % T has 5 rows per trial (one per finger)
        P.isError   = P.isError./5; % correct for number of rows per trial
        P.numTrials = P.numTrials./5;
        P.perErr    = P.isError./P.numTrials;
        
        varargout = {P,T};
    case 'SUBJ_getStimTime'
        % plots max stimulation peak times per finger
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % load data
        % arrange data into plotting structure
        load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T

        D = [];
        v = ones(2,1);
        for i = 1:length(T.run)
            d.sn          = v.*sn;
            d.run         = v.*T.run(i);
            d.trial       = v.*T.trial(i);
            d.finger      = v.*T.finger(i);
            d.chordNum    = v.*T.chordNum(i);
            d.numDigits   = v.*T.numDigits(i);
            d.stimulated  = v.*T.stimulated(i);
            d.peakF_stims = T.peakF_stims(i,:)';
            d.peakF_times = T.peakF_times(i,:)';
            d.stimNum     = [1,2]'; % faster than [1;2]
            D = addstruct(D,d);
        end
        
        varargout = {D};
    case 'SUBJ_plotStimTime'
        % plots max stimulation peak times per finger
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % get data
        D = pp1_imana('SUBJ_getStimTime','sn',sn);
        
        % plot
        style.use('5fingers');
        
        subplot(1,2,1);
        plt.hist(D.peakF_times,'split',D.finger,'subset',D.stimulated & D.stimNum==1);
        plt.labels('time (ms)','count',sprintf('peakStim 1'));
        xlim([0 4000]);
        legend off
        
        subplot(1,2,2);
        plt.hist(D.peakF_times,'split',D.finger,'subset',D.stimulated & D.stimNum==2);
        plt.labels('time (ms)','count',sprintf('peakStim 2'));
        xlim([0 4000]);
        plt.legend('east',{'thumb','index','middle','ring','little'});
        legend off
        
        varargout = {D};
    case 'GROUP_plotStimTime'
        % plots max stimulation peak times per finger
        sn = 1:length(subj_name);
        vararginoptions(varargin,{'sn'});
        % get data
        D = [];
        for s = sn
            d = pp1_imana('SUBJ_getStimTime','sn',s);
            D = addstruct(D,d);
        end
        % avg. across trials per subj
        D = tapply(D,{'sn','finger','stimulated','stimNum'},{'peakF_times','nanmean'});
        % plot
        style.use('5fingers');
        plt.bar(D.stimNum,D.peakF_times,'errorfcn','stderr','split',D.finger,'subset',D.stimulated==1);
        plt.set('xticklabel',{'thumb','index','middle','fourth','little','thumb','index','middle','fourth','little'},'xticklabelrotation',45);
        plt.labels('stim 1   |   stim 2','time (ms)');
        ylim([0 4000]);
        legend off
        
        varargout = {D};
    case 'GROUP_plotStimForce'
        % plots max stimulation force per finger 
        sn = 1:6;%length(subj_name);
        vararginoptions(varargin,{'sn'});
        % get data
        D = [];
        for s = sn
            load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s}))); % T
            D = addstruct(D,T);
        end
        D = tapply(D,{'sn','numDigits','finger','stimulated'},...
            {'peakF_raw','nanmean'},{'peakF_filt','nanmean'},{'time_stimOnset','nanmean'},{'forceStim','mean'});
        
        % plot
        style.use('5fingers');
        
        subplot(1,3,1);
        plt.line(D.numDigits,(D.peakF_filt./D.forceStim).*100,'split',D.finger,'errorfcn','stderr','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 125]);
        plt.labels('# fingers stimulated','% target force',sprintf('stimulated\ntrials'));
        drawline(100,'dir','horz','linestyle',':','linewidth',1.5);
        legend off
        
        subplot(1,3,2);
        plt.line(D.numDigits,(D.peakF_raw./D.forceStim).*100,'split',D.finger,'errorfcn','stderr','subset',D.stimulated==0);
        plt.set('xlim',[0.5 4.5],'ylim',[0 125]);
        plt.labels('# fingers stimulated','% target force',sprintf('non-stimulated\ntrials'));
        drawline(100,'dir','horz','linestyle',':','linewidth',1.5);
        legend off
        
        subplot(1,3,3);
        plt.line(D.numDigits,D.time_stimOnset,'split',D.finger,'errorfcn','stderr','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 1000]);
        plt.labels('# fingers stimulated','stimulus onset (ms)',sprintf('stimulus onset\ntime'));
        plt.legend('east',{'thumb','index','middle','ring','little'});
        
        varargout = {D};
    case 'GROUP_taskPerformance'
        % plots the percent correct trial responses (this is a little
        % cheeky as most trials don't require a response so it is easy to
        % get a high success rate)
        sn = 1:6;
        vararginoptions(varargin,{'sn'});
        P = []; % percent correct
        D = []; % d prime
        for s = sn
            p = pp1_imana('SUBJ_getTaskPerformance','sn',s);
            d = pp1_imana('SUBJ_getDprime','sn',s);
            P = addstruct(P,p);
            D = addstruct(D,d);
        end
        P = tapply(P,{'sn','numDigits'},{'perErr','mean'}); % avg. across sessions and chords
        subplot(1,5,[1:3]); % plot percent correct
        sty = style.custom({[0 0 0] [0.5 0 0] [0.9 0 0] [1 0.6 0] [0.4 0.4 0.4]});
        plt.dot(P.numDigits,1-P.perErr,'split',P.numDigits,'style',sty);
        ylabel('% correct (all trial types)');
        xlabel('number of fingers stimulated');
        drawline(0.5,'dir','horz','linestyle',':');
        ylim([0.25 1]);
        
        subplot(1,5,4); % plot % correct for mismatch trials
        sty = style.custom({[0.6 0 0.6]});
        plt.box([],1-D.perErr,'style',sty);
        ylabel('% correct (mistmatch trials)');
        set(gca,'xtick',[]);
        drawline(0.5,'dir','horz','linestyle',':');
        ylim([0.25 1]);

        subplot(1,5,5); % plot d prime per subject (sort of meaningless without comparison group but whatever)
        sty = style.custom({'blue'});
        plt.dot([],D.dprime,'style',sty);
        ylabel('d''');
        title(sprintf('signal\ndetection'));
        set(gca,'xtick',[]);
        drawline(0,'dir','horz','linestyle',':');
        varargout = {P,D};

    case '0' % ------------ PREP: preprocessing. Expand for more info. ----
        % The PREP cases are preprocessing cases. 
        % You should run these in the following order:
        %       'WRAPPER_dicom_import' :  does step below for all
        %                                  series_types.
        %       'PREP_dicom_import'*   :  call with 'series_type','functional', 
        %                                  and again with
        %                                  'series_type','anatomical'.
        %       'WRAPPER_preprocess1'  :  Runs steps 1.2 - 1.7 (see below).
        %       'PREP_coreg'*          :  Registers meanepi to anatomical img. (step 1.8)
        %       'WRAPPER_preprocess2'* :  Runs steps 1.9 - 1.11 (see below).
        %
        %   * requires user input/checks after running BEFORE next steps.
        %       See corresponding cases for more info about required
        %       user input.
        %
        % When calling any case, you can submit only one subj#.
        % Use 'WRAPPER_...' to loop through multiple subjects.
        %
        % Notes on prefixes appended to data in the preprocessing stages:
        %   - standard realignment ('_realign'; sans fieldmap correction) 
        %      appends prefix 'r' to all filenames
        %   - fieldmap correction ('_fieldmap_make' & '_fieldmap_RealignUnwarp') appends
        %      prefix 'u' to all filenames (fieldmap corr. is suggested!)
        %   - 'meanimage_bias_correction' appends prefix 'bb' + one of the
        %      prefixs mentioned above to filenames
        %
        %   ** Be sure you are aware what prefixes are appropriate. To
        %   ensure smooth sailing, either submit prefixes to 'PREP' stages
        %   with option 'prefix', or adjust local function
        %   'getCorrectPrefix' accordingly (this function assumes there is
        %   only one appropriate prefix per subject).
        %
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    case 'WRAPPER_dicomImport'                                             
        % Converts dicom to nifti files w/ spm_dicom_convert.
        % Comment out series_type not required.
        sn   = [];
        sess = [];
        vararginoptions(varargin,{'sn','sess'});
        for s = sn
            for i = sess
                % Did we acquire this subject's anatomical during this scan
                % session? If not, don't process.
                if ~isempty(anatNum{i}{s})
                    fprintf('importing ANATOMICAL runs- %s\n',subj_name{sn});
                    pp1_imana('PREP_dicomImport','sn',s,'seriesType','anatomical','sess',i);
                end
                if ~isempty(fscanNum{i}{s})
                    fprintf('importing FUNCTIONAL runs- %s\n',subj_name{sn});
                    pp1_imana('PREP_dicomImport','sn',s,'seriesType','functional','sess',i);
                    fprintf('concatinating 3d niftis to 4d per run\n');
                    pp1_imana('PREP_make4dNifti','sn',s,'sess',i);
                end
                if ~isempty(fieldNum{i}{s})
                    fprintf('importing FIELDMAP runs- %s\n',subj_name{sn});
                    pp1_imana('PREP_dicomImport','sn',s,'seriesType','fieldmap','sess',i);
                end
            end
        end
    case 'PREP_dicomImport'                                                
        % converts dicom to nifti files w/ spm_dicom_convert
        seriesType = 'functional';
        vararginoptions(varargin,{'sn','seriesType','sess'});
        cwd = pwd;
        switch seriesType
            case 'functional'
                seriesNum = fscanNum;
            case 'anatomical'
                seriesNum = anatNum;  
            case 'fieldmap'
                seriesNum = fieldNum;
        end

        subjSessDir = fullfile(dicomDir,sprintf('%s_sess%02d',subj_name{sn},sess));
        cd(subjSessDir);
        for i = 1:length(seriesNum{sess}{sn}) % per run
            % Get DICOM FILE NAMES
            r      = seriesNum{sess}{sn}(i);
            folder = fullfile(subjSessDir,sprintf('%4.4d',r));
            cd(folder)
            DIR    = dir(sprintf('%s.%4.4d.*.dcm',DicomName{sess}{sn},r));  
            Names  = vertcat(DIR.name);
            % Convert the dicom files
            if (~isempty(Names))
                % Load dicom headers
                HDR = spm_dicom_headers(Names,1);  
                % Make a directory for series{r} for this subject.
                % The nifti files will be saved here.
                dirname = fullfile(subjSessDir,sprintf('series%2.2d',r));
                dircheck(dirname);
                cd(dirname);
                % Convert the data to nifti
                spm_dicom_convert(HDR,'all','flat','nii');                  
                cd ..
            else
                error('no files in directory. did you forget to update the dicomName for this session?')
            end
            fprintf('Series %d done \n\n',seriesNum{sess}{sn}(i))
        end
        % Display verbose messages to user. 
        % Lazy and won't include none-verbose version here.
        switch seriesType
            case 'functional'
                fprintf('Subject %02d functional runs imported. Copy the unique .nii name for subj files and place into ''Subject Things''.\n',sn)
            case 'anatomical'
                fprintf('Anatomical runs have been imported for subject %d.\n',sn); 
                fprintf('Please locate the T1 weighted anatomical img. Copy it to the anatomical folder.\n')
                fprintf('Rename this file to ''s%02d_anatomical_raw.nii'' in the anatomical folder.\n',sn); 
            case 'fieldmap'
                fprintf('Subject %02d fieldmaps imported.\n',sn)
        end
        cd(cwd); 
    case 'PREP_make4dNifti'                                               
        vararginoptions(varargin,{'sn','sess'});
        subjSessDir = fullfile(dicomDir,sprintf('%s_sess%02d',subj_name{sn},sess));
        % For each functional run
        for i = 1:length(fscanNum{sess}{sn})                                      
            outfilename = fullfile(imagingDirRaw,subj_name{sn},sprintf('sess_%02d_raw',sess),sprintf('%s_run_%02d.nii',subj_name{sn},run{sess}{sn}(i)));
            % Create a 4d nifti of all functional imgs in this run.
            % Don't include the first few dummy scans in this 4d nifti.
            for j = 1:(numTRs{sn}(i)-numDummys)                                        
                P{j}=fullfile(subjSessDir,sprintf('series%2.2d',fscanNum{sess}{sn}(i)),...
                    sprintf('f%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{sess}{sn},fscanNum{sess}{sn}(i),j+numDummys,j+numDummys));
            end
            dircheck(fullfile(imagingDirRaw,subj_name{sn},sprintf('sess_%02d_raw',sess)));
            spm_file_merge(char(P),outfilename);
            fprintf('Sess %02d Run %02d done\n',sess,run{sess}{sn}(i));
        end
    
    case 'PREP_moveRaw4d'
        % move unrealigned 4d niftis to a processing directory (realignVDM)
        % Moves image data from imaging_dicom_raw into a "working dir":
        % imaging_dicom.
        postfx = 'realignVDM';
        vararginoptions(varargin,{'sn','postfx'});
        % get filenames to copy and where to copy
        for q = 1:numSess(sn)
            sourceDir = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},sprintf('sess_%02d_raw',q));
            % check target output directory exists
            destDir = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},sprintf('sess_%02d_%s',q,postfx));
            dircheck(destDir);
            for r = (run{q}{sn})
                % copy 4d niftis
                source = fullfile(sourceDir, sprintf('%s_run_%02d.nii',subj_name{sn},r));
                dest   = fullfile(destDir, sprintf('%s_run_%02d.nii',subj_name{sn},r));
                display(dest)
                copyfile(source,dest);
            end % run
        end % session
        fprintf('done.\n')
    case 'PREP_realignEst'
        % This ESTIMATES realignment only (no reslicing).
        vararginoptions(varargin,{'sn'});
        data = {};
        for q = 1:numSess(sn)
            sourceDir = fullfile(imagingDirRaw,subj_name{sn},sprintf('sess_%02d_realignVDM',q));
            for r = run{q}{sn}
                for i = 1:(numTRs{sn}(r)-numDummys);
                    data{r}{i,1} = fullfile(sourceDir,sprintf('%s_run_%2.2d.nii,%d',subj_name{sn},r,i));
                end
            end
        end
        % set spm batch options
        J.data = data;
        J.eoptions.quality      = 0.9;                                                                            
        J.eoptions.sep          = 2;                                                                                  
        J.eoptions.fwhm         = 5;                                                                                 
        J.eoptions.rtm          = 1;  % two pass realignment (register to first image in run, then all images to mean of the images after first realignment)                                                                       
        J.eoptions.interp       = 2;                                                                               
        J.eoptions.wrap         = [0 0 1];                                                                           
        J.eoptions.weight       = '';          
        % do
        matlabbatch{1}.spm.spatial.realign.estimate = J;
        spm_jobman('run',matlabbatch);
        fprintf('Subj %d realignment estimated\n',sn); 
    case 'PREP_calcVDM'
        vararginoptions(varargin,{'sn'});
        
        % To calculate total readout time (tert) for SPM:
        % tert = 1000/('Bandwith Per Pixel Phase Encode (dicom tag 0019,1028)' * GRAPPA accel. factor)
        % See: https://groups.google.com/forum/#!topic/7t_trt/fDR35kgU4mU
        
        % helpful file for many of these options: tbx_cfg_fieldmap.m
        
        % Set options for batch job
        spm_dir= fileparts(which('spm'));
        %J.defaults.defaultsfile = {'/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_passivePatterns/pm_pp1.m'};
        J.defaults.defaultsval.et              = [4.08 5.1];    % [shortest, longest echotimes]                                                       
        J.defaults.defaultsval.etd             = J.defaults.defaultsval.et(2) - J.defaults.defaultsval.et(1);
        J.defaults.defaultsval.tert            = 1000/(22.904*2);
        J.defaults.defaultsval.maskbrain       = 1;             % mask brain w/ f magnitude img before processing                                                  
        J.defaults.defaultsval.blipdir         = -1;                                                                                                                                     
        J.defaults.defaultsval.epifm           = 0;             % are field maps based on epi acquistion (1) or no (0)?                                                                    
        J.defaults.defaultsval.ajm             = 0;             % apply Jacobian modulation to unwarped epi img -> suggested always NO                                    
        J.defaults.defaultsval.uflags.method   = 'Mark3D';      % default phase-unwarping method                                                     
        J.defaults.defaultsval.uflags.fwhm     = 10;                                                                
        J.defaults.defaultsval.uflags.pad      = 0;                                                                  
        J.defaults.defaultsval.uflags.ws       = 1;   
        J.defaults.defaultsval.mflags.template = {fullfile(spm_dir,'canonical','avg152T1.nii')};
        J.defaults.defaultsval.mflags.fwhm     = 5;                                                                 
        J.defaults.defaultsval.mflags.nerode   = 2;                                                               
        J.defaults.defaultsval.mflags.ndilate  = 4;                                                               
        J.defaults.defaultsval.mflags.thresh   = 0.5;                                                             
        J.defaults.defaultsval.mflags.reg      = 0.02;        
        J.matchvdm                             = 1;              % coregister field map data to the selected EPI for each run                                                                                                                                                                
        J.writeunwarped                        = 1;                                                                                    
        J.anat                                 = [];                  
        J.matchanat                            = 0; 
        % loop through sessions
        for q = 1:numSess(sn)
            % calculate VDM for each run per session of this subject
            rawdataDir = fullfile(imagingDirRaw,subj_name{sn},sprintf('sess_%02d_realignVDM',q));
            J.session  = struct;
            J.sessname = {};
            for i = 1:numel(run{q}{sn})
                r = run{q}{sn}(i);
                % VDM file will be created for each run (based on the first image in each run)
                J.session(i).epi = {fullfile(rawdataDir,sprintf('%s_run_%02d.nii,1',subj_name{sn},r))};
                J.sessname(i)    = {sprintf('%s_run_%02d',subj_name{sn},r)};
            end
            J.data.presubphasemag.phase     = {fullfile(fieldmapDir,subj_name{sn},sprintf('sess_%02d',q),[subj_name{sn},'_phase.nii,1'])};   
            J.data.presubphasemag.magnitude = {fullfile(fieldmapDir,subj_name{sn},sprintf('sess_%02d',q),[subj_name{sn},'_magnitude.nii,1'])};
            % run for each session seperately 
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj = J;
            spm_jobman('run',matlabbatch);
        end
    case 'PREP_applyVDM'
        % apply voxel displacement map to each session.
        % Each image in sessions should be realigned to the first image in
        % the timeseries (realign:estimate), but don't need to be resliced.
        vararginoptions(varargin,{'sn'});
        % set spm job options for applyVDM
        J.roptions.pedir   = 3; % apply in Z phase-encoding direction
        J.roptions.which   = [2 1]; % 
        J.roptions.rinterp = 4;
        J.roptions.wrap    = [0 0 0]; % no unwrapping in any specified direction (no major aliasing present)
        J.roptions.mask    = 1;
        J.roptions.prefix  = 'u';
        % loop through sessions
        for q = 1:numSess(sn)
            % calculate VDM for each run per session of this subject
            rawdataDir = fullfile(imagingDirRaw,subj_name{sn},sprintf('sess_%02d_realignVDM',q));
            % apply per each run
            for i = 1:numel(run{q}{sn})
                r = run{q}{sn}(i);
                scans = {};
                for j = 1:numTRs{sn}(i)-numDummys
                    scans{j,1} = fullfile(rawdataDir,sprintf('%s_run_%02d.nii,%d',subj_name{sn},r,j));                               
                end
                % apply each run's VDM file to scans from this run
                J.data(r).vdmfile = {fullfile(fieldmapDir,subj_name{sn},sprintf('sess_%02d',q),['vdm5_sc',subj_name{sn},'_phase_session',num2str(i),'.nii,1'])};
                J.data(r).scans   = scans;
            end
        end
        
        matlabbatch{1}.spm.tools.fieldmap.applyvdm = J;
        spm_jobman('run',matlabbatch);
    
    case 'PREP_realignUnwarp'
        % fieldmap_RealignUnwarp_multiSess
        vararginoptions(varargin,{'sn'});
        % Estimation options
        J.eoptions.quality    = 0.9;        %    0.9                   - quality vs. speed trade-off. range: 0 (least) to 1 (best quality)
        J.eoptions.sep        = 2;          %    4                     - separation (mm) b/t points sampled in ref. image (smaller is more accurate, but slower)                                                                                  
        J.eoptions.fwhm       = 5;          %    5                     - FWHM of smoothing kernel (mm) applied to imgs prior to estimating REALIGNMENT params                                                                            
        J.eoptions.rtm        = 0;          %    0 (one pass)          - what img all other imgs are registered to. One-pass registers to first img (good for fMRI)                                                                       
        J.eoptions.einterp    = 3;          %    2 (degree b-spline)   - interpolation method (imgs are sampled when estimating optimum transformation). Higher degree = better interp, slower). range:0-7                                                                      
        J.eoptions.ewrap      = [0 0 0];    %    [0 0 0]               - wrap-around in the [x y z] direction during the estimation (ewrap).                                                                     
        J.eoptions.weight     = {''};       %    ''                    - weighting img to weight voxels of reference img                                                                  
        % Unwarp estimation options
        J.uweoptions.basfcn   = [12 12];    %    [12 12 *empty*]       - number of basis functions used in each dimension                                                                 
        J.uweoptions.regorder = 1;          %    1                     - order of the regularisation- unwarp solution minimizes variance while maxing smoothness of estimated field. Reg param determines balance between the two. range: 0-3                                                                  
        J.uweoptions.lambda   = 100000;     %    100000 (medium)       - regularisation factor (not sure...). 3 options: 'a little','medium','a lot'.                                                                    
        J.uweoptions.jm       = 0;          %    0 (no)                - include Jacobian intensity modulation when estimating fields. Recommend NOT to use.                                                          
        J.uweoptions.fot      = [4 5];      %    [4 5] (pitch & roll)  - first-order effects: model movements effects. Movement types referred to by number. 1:3 = x,y,z translation; 4:6 = x,y,z rotation.                                                                      
        J.uweoptions.sot      = [1];        %    []                    - second-order (derivative) effects: entered similar as above.                                                                          
        J.uweoptions.uwfwhm   = 4;          %    4                     - FWHM of smoothing kernel (mm) applied to imgs prior to estimating DEFORMATION fields                                                                                                           
        J.uweoptions.rem      = 1;          %    1                     - re-estimate movement params: 'YES' means movement-params will be re-estimated at each unwarping iteration                                                                       
        J.uweoptions.noi      = 5;          %    5                     - max number of unwarp iterations                                                           
        J.uweoptions.expround = 'Average';  %    'Average'             - taylor expansion point: point in position space to perform expansion around. 'First','Last','Average'. (avg should yeild best variance reduction). If fieldmap acquired before timeseries, consider 'First'...
        % Unwarp reslicing options
        J.uwroptions.uwwhich  = [2 1];      %    [2 1]                 - Reslice imgs: first val specifies to reslice all imgs, second val creates mean of resliced imgs
        J.uwroptions.rinterp  = 4;          %    4                     - reslicing interpolation. see estimation options for description (above)                                                                   
        J.uwroptions.wrap     = [0 0 0];    %    []                    - wrap-around in the [x y z] direction during the reslicing (wrap)                                                                                
        J.uwroptions.mask     = 1;          %    1 (mask imgs)         - masks voxels that exist outside original imgs (reference to subject motion in certain timeseries imgs)                                                
        J.uwroptions.prefix   = 'u';        %    'u'                   - filename prefix added to names of smoothed imgs
        
        sCount = 1;
        for q = 1:numSess(sn)
            rawdataDir = fullfile(imagingDirRaw,subj_name{sn},sprintf('sess_%02d_unwarp',q));
            for j = 1:numel(run{q}{sn})
                r = run{q}{sn}(j);
                scans = {};
                for i = 1:numTRs{sn}(j)-numDummys
                    scans{i,1} = fullfile(rawdataDir,sprintf('%s_run_%02d.nii,%d',subj_name{sn},r,i));                               
                end
                J.data(sCount).scans  = scans;
                J.data(sCount).pmscan = {fullfile(fieldmapDir,subj_name{sn},sprintf('sess_%02d',q),['vdm5_sc',subj_name{sn},'_phase_session',num2str(j),'.nii,1'])};
                sCount = sCount+1;
            end 
        end;

        matlabbatch{1}.spm.spatial.realignunwarp = J;
        spm_jobman('run',matlabbatch);   
    case 'PREP_realignEstReslice'                                                    
        % SPM realigns first volume in each run to first volume of first
        % run, and then registers each image in that run to the first
        % volume of that run. Hence also why it's often better to run
        % anatomical before functional scans.

        % SPM does this with 4x4 affine transformation matrix in nifti
        % header (see function 'coords'). These matrices convert from voxel
        % space to world space(mm). If the first image has an affine
        % transformation matrix M1, and image two has one (M2), the mapping
        % from 1 to 2 is: M2/M1 (map image 1 to world space-mm - and then
        % mm to voxel space of image 2).

        % Registration determines the 6 parameters that determine the rigid
        % body transformation for each image (described above). Reslice
        % conducts these transformations; resampling each image according
        % to the transformation parameters. This is for functional only!
        
        % This ESTIMATES realignment and preforms RESLICING.
        
        % Appends prefix 'r' to realigned imgs.
        vararginoptions(varargin,{'sn'});
        data = {};
        for q = 1:numSess(sn)
            sourceDir = fullfile(imagingDirRaw,subj_name{sn},sprintf('sess_%02d_realigned',q));
            for r = (run{q}{sn});
                for i = 1:(numTRs{sn}(r)-numDummys);
                    data{r}{i,1} = fullfile(sourceDir,sprintf('%s_run_%2.2d.nii,%d',subj_name{sn},r,i));
                end
            end
        end
        spmj_realign(data);
        fprintf('Subj %d realigned\n',sn);
    
    case 'PREP_moveData'                                                   
        % Moves image data from imaging_dicom_raw into a "working dir":
        % imaging_dicom.
        postfx = 'realignVDM';
        vararginoptions(varargin,{'sn','postfx'});
        prefix = dataPrefix{sn};
        % check target output directory exists
        destDir = fullfile(baseDir, 'imaging_data',subj_name{sn});
        dircheck(destDir);
        % get filenames to copy and where to copy
        for q = 1:numel(run)
            sourceDir = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},sprintf('sess_%02d_%s',q,postfx));
            for r = (run{q}{sn})
                % copy 4d niftis
                source = fullfile(sourceDir, sprintf('%s%s_run_%02d.nii',prefix,subj_name{sn},r));
                dest   = fullfile(destDir, sprintf('%s%s_run_%02d.nii',prefix,subj_name{sn},r));
                copyfile(source,dest);
                disp(dest);
                % copy realignment txt files
                source = fullfile(sourceDir, sprintf('rp_%s_run_%02d.txt',subj_name{sn},r));
                dest   = fullfile(destDir, sprintf('rp_%s_run_%02d.txt',subj_name{sn},r));
                copyfile(source,dest);
            end % run
        end % session
        % copy mean epi img (from sess_01 raw data dir)
        if strcmp(prefix,'u') % fieldmap vs. affine alignment have different naming conventions
            source = fullfile(imagingDirRaw,subj_name{sn},['sess_01_' postfx], sprintf('mean%s%s_run_01.nii',prefix,subj_name{sn}));
        elseif strcmp(prefix,'r')
            source = fullfile(imagingDirRaw,subj_name{sn},['sess_01_' postfx], sprintf('mean%s_run_01.nii',subj_name{sn}));
        end
        dest   = fullfile(destDir, sprintf('%smeanepi_%s.nii',prefix,subj_name{sn}));
        copyfile(source,dest);
        fprintf('Moved niftis to working dir.\n')

    case 'PREP_resliceLPI'                                                 
        vararginoptions(varargin,{'sn'});
        % (1) Reslice anatomical image to set it within LPI co-ordinate frames
        source  = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical_raw','.nii']);
        dest    = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii']);
        spmj_reslice_LPI(source,'name', dest);
        % (2) In the resliced image, set translation to zero
        V               = spm_vol(dest);
        dat             = spm_read_vols(V);
        V.mat(1:3,4)    = [0 0 0];
        spm_write_vol(V,dat);
        disp 'Done'
    case 'PREP_centreAC'                                                   
        % Set origin of anatomical to anterior commissure (must provide
        % coordinates in section (4)).
        vararginoptions(varargin,{'sn'});
        img    = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii']);
        V               = spm_vol(img);
        dat             = spm_read_vols(V);
        V.mat(1:3,4)    = loc_AC{sn};
        spm_write_vol(V,dat);
        disp 'Done'
    case 'PREP_segmentation'                                               
        vararginoptions(varargin,{'sn'});

        SPMhome=fileparts(which('spm.m'));

        J.channel.vols     = {fullfile(anatomicalDir,subj_name{sn},[subj_name{sn},'_anatomical.nii,1'])};
        J.channel.biasreg  = 0.001;
        J.channel.biasfwhm = 60;
        J.channel.write    = [0 0];
        
        % (1) gray matter
        J.tissue(1).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,1')};   % tissue probability map
        J.tissue(1).ngaus  = 1;                                     % # gaussians (?)- leave alone
        J.tissue(1).native = [1 1];                                 % writes segmented images at resolution of native (1) and dartel (2)
        J.tissue(1).warped = [0 0];                                 % don't write warped images
        % (2) white matter
        J.tissue(2).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,2')};
        J.tissue(2).ngaus  = 1;
        J.tissue(2).native = [1 1];
        J.tissue(2).warped = [0 0];
        % (3) CSF
        J.tissue(3).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,3')};
        J.tissue(3).ngaus  = 2;
        J.tissue(3).native = [1 0];
        J.tissue(3).warped = [0 0];
        % (4) skull
        J.tissue(4).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,4')};
        J.tissue(4).ngaus  = 3;
        J.tissue(4).native = [1 0];
        J.tissue(4).warped = [0 0];
        % (5) soft tissue outside brain
        J.tissue(5).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,5')};
        J.tissue(5).ngaus  = 4;
        J.tissue(5).native = [1 0];
        J.tissue(5).warped = [0 0];
        % (6) other stuff, largely outside brain
        J.tissue(6).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,6')};
        J.tissue(6).ngaus  = 2;
        J.tissue(6).native = [0 0];
        J.tissue(6).warped = [0 0];
        
        J.warp.mrf         = 1;
        J.warp.cleanup     = 1;
        J.warp.reg         = [0 0.001 0.5 0.05 0.2];
        J.warp.affreg      = 'mni';
        J.warp.fwhm        = 0;
        J.warp.samp        = 3;
        J.warp.write       = [0 0];
        
        matlabbatch{1}.spm.spatial.preproc=J;
        spm_jobman('run',matlabbatch);
        fprintf('Check segmentation results for %s\n', subj_name{sn})
    
    case 'PREP_coreg'                                                      
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi
        %   image
        vararginoptions(varargin,{'sn'});
        prefix = dataPrefix{sn};
        cd(fullfile(anatomicalDir,subj_name{sn}));
        coregtool;
        keyboard();
        
        % (2) Automatically co-register functional and anatomical images
        
        J.ref    = {fullfile(anatomicalDir,subj_name{sn},[ subj_name{sn}, '_anatomical','.nii'])};
        J.source = {fullfile(imagingDir,subj_name{sn},['r' char(prefix) 'meanepi_' subj_name{sn} '.nii'])}; 
        J.other  = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
        % (3) Manually check again
        coregtool;
        keyboard();
        
        % NOTE:
        % Overwrites meanepi, unless you update in step one, which saves it
        % as rmeanepi.
        % Each time you click "update" in coregtool, it saves current
        % alignment by appending the prefix 'r' to the current file
        % So if you continually update rmeanepi, you'll end up with a file
        % called r...rrrmeanepi.
        
        %__________________________________________________________________
    case 'PREP_makeSameAlign'                                             
        % apply rigid-body transform from mean epi coregistration to all
        % functional images
        vararginoptions(varargin,{'sn'});
        prefix  = dataPrefix{sn};
        subjDir = fullfile(imagingDir,subj_name{sn});
        % Select coregistered epi image to be image we align to
        P{1} = fullfile(subjDir,sprintf('r%smeanepi_%s.nii',prefix,subj_name{sn}));
        % Select images to be realigned to subject's coregistered mean epi
        N = findFiles(subjDir,[prefix subj_name{sn} '_run_']); % find run files for subject
        Q = {};
        for r = 1:numel(N)
            for i = 1:(numTRs{sn}(r)-numDummys)
                % algin each img within 4d nifti to coregistered epi
                Q{end+1} = fullfile(subjDir,[N{r},',',num2str(i)]);
            end
        end
        spmj_makesamealign_nifti(char(P),char(Q));
        fprintf('Done. Run spmj_checksamealign to check alignment.\n')
        spmj_checksamealign
    case 'PREP_makeMaskImage'                                             
        vararginoptions(varargin,{'sn'});
        prefix = dataPrefix{sn};
        cd(fullfile(imagingDir,subj_name{sn}));

        nam{1}  = fullfile(imagingDir,subj_name{sn}, ['r' prefix 'meanepi_' subj_name{sn} '.nii']);
        nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
        nam{3}  = fullfile(anatomicalDir, subj_name{sn}, ['c2' subj_name{sn}, '_anatomical.nii']);
        nam{4}  = fullfile(anatomicalDir, subj_name{sn}, ['c3' subj_name{sn}, '_anatomical.nii']);
        spm_imcalc(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2');
        
        nam{1}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
        nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c2' subj_name{sn}, '_anatomical.nii']);
        nam{3}  = fullfile(anatomicalDir, subj_name{sn}, ['c3' subj_name{sn}, '_anatomical.nii']);
        spm_imcalc(nam, 'rmask_noskull.nii', '(i1+i2+i3)>0.2');

        nam     = {};
        nam{1}  = fullfile(imagingDir,subj_name{sn}, ['r' prefix 'meanepi_' subj_name{sn} '.nii']);
        nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
        spm_imcalc(nam, 'rmask_gray.nii', 'i1>1 & i2>0.4');

    case '0' % ------------ SURF: Freesurfer funcs. Expand for more info. -
        % The SURF cases are the surface reconstruction functions. Surface
        % reconstruction is achieved via freesurfer.
        % All functions can be called with ('WRAPPER_SURF','sn',[Subj#s]).
        % You can view reconstructed surfaces with Caret software.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'WRAPPER_WB'                                                 
        vararginoptions(varargin,{'sn'});
        % You can call this case to do all the freesurfer & workbench processing.
        % 'sn' can be an array of subjects.
        for s = sn 
            pp1_imana('SURF_freesurfer','sn',s);
            pp1_imana('SURF_WBresample','sn',s);
        end
    case 'WRAPPER_CARET'                                                 
        vararginoptions(varargin,{'sn'});
        % You can call this case to do all the freesurfer & caret processing.
        % 'sn' can be an array of subjects.
        for s = sn 
            pp1_imana('SURF_freesurfer','sn',s);
            pp1_imana('SURF_xhemireg','sn',s);
            pp1_imana('SURF_map_ico','sn',s);
            pp1_imana('SURF_make_caret','sn',s);
        end
    case 'SURF_freesurfer'   % run reconall   
        vararginoptions(varargin,{'sn'});
        freesurfer_reconall(freesurferDir,subj_name{sn},fullfile(anatomicalDir,subj_name{sn},[subj_name{sn} '_anatomical.nii']));
    
    case 'SURF_WBresample'   % Reslice indiv surfaces into fs_lr standard mesh
        % This reslices from the individual surfaces into the the fs_lr
        % standard mesh - This replaces calls to freesurfer_registerXhem,
        % freesurfer_mapicosahedron_xhem, & caret_importfreesurfer. It
        % requires connectome wb to be installed, added to the bash_profile
        % (on terminal), and updated on the startup.m file
        atlasDir = fullfile(atlasDir, 'standard_mesh');
        vararginoptions(varargin, {'sn'});
        fprintf('reslicing %s...',subj_name{sn});
        surf_resliceFS2WB(subj_name{sn}, freesurferDir, wbDir,'resolution','32k'); 
        fprintf('done\n');
        
    case 'SURF_xhemireg'     % Cross-Register surfaces left / right hem  
        vararginoptions(varargin,{'sn'});
        freesurfer_registerXhem({subj_name{sn}},freesurferDir,'hemisphere',[1 2]); % For debug... [1 2] orig
    case 'SURF_map_ico'      % Align subj surface to atlas surface 
        vararginoptions(varargin,{'sn'});
        freesurfer_mapicosahedron_xhem(subj_name{sn},freesurferDir,'smoothing',1,'hemisphere',[1:2]);
    case 'SURF_make_caret'   % convert freesurfer to caret  
        vararginoptions(varargin,{'sn'});
        caret_importfreesurfer(['x' subj_name{sn}],freesurferDir,caretDir);
            
    case '0' % ------------ CM: Cerebellum analyses. ----------------------
        % The CM cases are the cerebellum analysis cases.
        % Outputs for segmentations are saved in the c_anatomicals folder.
        % Uses the SUIT toolbox.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'SUIT_isoSeg'                  % isolate & segment cerebellum
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % make copy of anatomical in cerebellar anatomical directory
        subjCAnatDir = fullfile(cerebAnatDir,subj_name{sn});
        dircheck(subjCAnatDir);
        source = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn} '_anatomical.nii']);
        dest   = fullfile(subjCAnatDir,[subj_name{sn} '_anatomical.nii']);
        copyfile(source,dest);
        % segment cerebellum
        suit_isolate_seg({dest});
        % detele copied anatomical
        delete(dest);
        fprintf('Done. Check cerebellum segmentation for %s. Make appropriate corrections if needed.\n',subj_name{sn});
    case 'SUIT_normalize&Reslice'       % calc deformation fields and reslice subject cerebellum into atlas space
        % Normalize segmented cerebellum (estimate non-lienar deformation
        % flowfield & reslice segemnted cereb into SUIT template space).
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % estimate flow-fields
        job1.subjND(1).gray      = {fullfile(cerebAnatDir,subj_name{sn},[subj_name{sn} '_anatomical_seg1.nii'])};
        job1.subjND(1).white     = {fullfile(cerebAnatDir,subj_name{sn},[subj_name{sn} '_anatomical_seg2.nii'])};
        job1.subjND(1).isolation = {fullfile(cerebAnatDir,subj_name{sn},['c_' subj_name{sn} '_anatomical_pcereb_corr.nii'])}; % corrected version
        suit_normalize_dartel(job1);
        % reslice cerebellum into SUIT atlas space
        job2.subj(1).affineTr  = {fullfile(cerebAnatDir,subj_name{sn},['Affine_' subj_name{sn} '_anatomical_seg1.mat'])};
        job2.subj(1).flowfield = {fullfile(cerebAnatDir,subj_name{sn},['u_a_' subj_name{sn} '_anatomical_seg1.nii'])};
        job2.subj(1).resample  = {fullfile(cerebAnatDir,subj_name{sn},['c_' subj_name{sn} '_anatomical.nii'])};
        job2.subj(1).mask      = {fullfile(cerebAnatDir,subj_name{sn},['c_' subj_name{sn} '_anatomical_pcereb_corr.nii'])}; % corrected version
        suit_reslice_dartel(job2);
        fprintf('Done. Check resliced cerebellum for %s.\n',subj_name{sn});
    case 'SUIT_mapAtlasRois'            % use inverse flow field to map atlas rois to single subject
        % Use inverse subject flow fields to warp cerebellar atlas into
        % subject space.
        sn = 1;
        vararginoptions(varargin,{'sn'});
        atlasName = 'MDTB_10_1mm.nii';
        % copy parcellation atlas to subject directory
        global defaults
        source = fullfile(defaults.tbx.dir{1},'suit','atlas',atlasName);
        dest   = fullfile(cerebAnatDir,subj_name{sn},atlasName);
        copyfile(source,dest);
        % do inverse warping
        job.Affine    = {fullfile(cerebAnatDir,subj_name{sn},['Affine_' subj_name{sn} '_anatomical_seg1.mat'])};
        job.flowfield = {fullfile(cerebAnatDir,subj_name{sn},['u_a_' subj_name{sn} '_anatomical_seg1.nii'])};
        job.resample  = {dest};
        job.ref       = {fullfile(cerebAnatDir,subj_name{sn},['c_' subj_name{sn} '_anatomical.nii'])}; % corrected version
        suit_reslice_dartel_inv(job);
        delete(dest);
        % rename output image to more friendly name
        tmp    = strsplit(atlasName,'.');
        source = fullfile(cerebAnatDir,subj_name{sn},['iw_' tmp{1} '_u_a_' subj_name{sn} '_anatomical_seg1.nii']);
        dest   = fullfile(cerebAnatDir,subj_name{sn},[subj_name{sn} '_MDTB_10.nii']);
        copyfile(source,dest);
        delete(source);
        fprintf('Done. Check warped atlas rois for %s.\n',subj_name{sn});
    case 'depreciated_CM_getGrayBetas'  % resample patterns from all gray matter voxels (irrespective of region)
        sn  = 1;
        glm = 2;
        % load in beta vol data
        HCPDir=fullfile(atlasDir,'HCP/');
        cd(HCPDir);
        HCPFiles=dir('*HCP_*');
        
        % get V and volIndx
        load(fullfile(studyDir{2},'encoding','glm4','cereb_avrgDataStruct.mat'));
        
        % get HCP contrasts
        for i=1:length(HCPFiles),
            VA{i}=spm_vol(fullfile(HCPDir,HCPFiles(i).name));
        end
        
        % now sample the contrasts into the same space
        [i,j,k]=ind2sub(V.dim,volIndx);
        [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
        [i1,j1,k1]=spmj_affine_transform(x,y,z,inv(VA{1}.mat));
        for i=1:length(VA),
            map(i,:) = spm_sample_vol(VA{i},i1,j1,k1,0);
            colNames{i,1}=HCPFiles(i).name(5:end-9);
        end
        
        % normalise data
        X_C=bsxfun(@minus,map,nanmean(map));
        
        C{1}.dim=V.dim;
        C{1}.mat=V.mat;
        
        % make volume
        Vi=zeros(size(map,1), [V.dim(1)*V.dim(2)*V.dim(3)]);
        Vi(:,volIndx)=map;
        
        % map vol2surf
        for i=1:size(map,1),
            data=reshape(Vi(i,:,:,:),[C{1}.dim]);
            C{i}.dat=data;
        end
        S=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',colNames);  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
        
        caret_save(fullfile(HCPDir,'HCPContrasts.metric'),S);
        varargout={X_C,colNames};    
    case 'depreciated_CM_logicalGrayMatterMask'
        % define t/f gray matter mask for cerebellum
        sn = 1;
        vararginoptions(varargin,{'sn'});
        G                       = spm_vol(fullfile(cerebAnatDir,subj_name{sn},[subj_name{sn} '_anatomical_seg1.nii']));
        G.data                  = spm_read_vols(G);
        G.data(G.data>0.25)     = 1; % mask out gray matter voxels
        G.data(G.data<1)        = 0;
        G.fname                 = fullfile(cerebAnatDir,subj_name{sn},'gray_mask.nii');
        G                       = rmfield(G,{'pinfo'});
        spm_write_vol(G,G.data);
        fprintf('Done %s cerebellar gray matter mask.\n',subj_name{sn});
    case 'depreciated_CM_defineFunctionalMask'
        % make functional mask for volumetric searchlight of cerebellar
        % data.
        glm = 2;
        sn  = 1;
        vararginoptions(varargin,{'sn','glm'});
        
        % load mask files to combine (functional and anatomical gray)
        F      = spm_vol(fullfile(glmDir{glm},subj_name{sn},'mask.nii'));
        F.data = spm_read_vols(F);
        G      = spm_vol(fullfile(cerebAnatDir,subj_name{sn},'gray_mask.nii'));
        G.data = spm_read_vols(G);
        Vin(1) = F;
        Vin(2) = G;
        Vo       = Vin(1); 
        Vo.fname = 'func_grayMask.nii'; 
        % find gray matter voxels of cerebellum for which we have
        % functional imaging data (make mask img in functional coords)
        spm_imcalc(Vin,Vo,'i1.*i2');
        fprintf('Done.\n')
    case 'CM_defineSearchVol'           % define volumetric searchlight in cerebellum 
        % Defines searchlights for 200 voxels in gray matter of cerebellum
        glm = 2;
        sn  = 1;
        numvox = 200;
        fig = 1;
        vararginoptions(varargin,{'sn','glm','numvox','fig'});
        
        % load masks
        M           = spm_vol(fullfile(glmDir{glm},subj_name{sn},'mask.nii')); % load functional mask structure that defines the voxels available for the searchlight and the voxel space.
        M.data      = spm_read_vols(M);
        ROI         = spm_vol(fullfile(cerebAnatDir,subj_name{sn},[subj_name{sn} '_anatomical_seg1.nii'])); % will define a volumetric searchlight over the volume of the anatomical grey matter
        ROI.data    = spm_read_vols(ROI);
        % define volumetric searchlight
        L = rsa.defineSearchlight_volume(ROI,M,'sphere',[50 numvox]);
        save(fullfile(cerebAnatDir,subj_name{sn},sprintf('%s_searchVol_%d.mat',subj_name{sn},numvox)),'-struct','L'); 
        
        if fig==1
            % visualize searchlight to check
            figure('Color',[1 1 1]);
            centers = surfing_inds2subs(M.dim,L.voxel);
            [coords(:,1),coords(:,2),coords(:,3)] = spmj_affine_transform(centers(:,1),centers(:,2),centers(:,3),M.mat);
            plot3(coords(:,1),coords(:,2),coords(:,3),'sk'); axis equal;
            hold on;clear coords
            centers = surfing_inds2subs(M.dim,L.LI{1000});
            centers = double(centers);
            [coords(:,1),coords(:,2),coords(:,3)] = spmj_affine_transform(centers(:,1),centers(:,2),centers(:,3),M.mat);
            plot3(coords(:,1),coords(:,2),coords(:,3),'.r'); axis equal;
            clear coords
            centers = surfing_inds2subs(M.dim,L.LI{end});
            centers = double(centers);
            [coords(:,1),coords(:,2),coords(:,3)] = spmj_affine_transform(centers(:,1),centers(:,2),centers(:,3),M.mat);
            plot3(coords(:,1),coords(:,2),coords(:,3),'.b'); axis equal;
            clear coords
        end
        
        fprintf('Done.\n')
    case 'CM_runVolSearchlight'         % run volumetric cerebellum searchlight analysis
        % Requires java functionality unless running on SArbuckle's
        % computer.
        glm = 2;
        sn  = 1;
        vararginoptions(varargin,{'sn','glm'});
        numvox = 200;
        
        block = 5e7; % set block size for searchlight function (larger = faster)
        cwd   = pwd; % copy current directory (to return to later)
        
        % copy spm directory to cerebellum analysis folder
        cd(fullfile(glmDir{glm},subj_name{sn}));
        D = load('SPM_info.mat');
        % load subject's searchlight definitions and SPM file
        L = load(fullfile(cerebAnatDir,subj_name{sn},sprintf('%s_searchVol_%d.mat',subj_name{sn},numvox)));
        % output filename
        name = sprintf('%s_cVolSearch_glm%d',subj_name{sn},glm);
        % run the searchlight
        rsa.runSearchlightLDC(L,'conditionVec',D.chord,'partition',D.run,'analysisName',name,'idealBlock',block,'java',0);

        cd(cwd);     
    
    case 'SUIT_mapVolLDC'                 % write volumetric file for results from cerebellar searchlight 
        % maps avg. paired ldc to a nifti file
        sn  = 1;
        glm = 2;
        con = {'1finger','2finger','3finger','4finger'};
        name = 'LDC';
        vararginoptions(varargin,{'sn','glm','con','name'});
        % Use 'con' option to define different contrasts.
        %   'avg_dist'   :  Average LDC nii for all 5 conds
        %   'avg_pattern':  Submits contrast matrix to MISC_SEARCH_calculate_contrast
        for s = sn
            % Load subject surface searchlight results (1 vol per paired conds)
            LDC_file            = fullfile(glmDir{glm},subj_name{s},sprintf('%s_cVolSearch_glm%d_%s.nii',subj_name{sn},glm,name)); % searchlight nifti
            [subjDir,fname,ext] = fileparts(LDC_file);
            cd(subjDir);
            vol  = spm_vol([fname ext]);
            vdat = spm_read_vols(vol); % is searchlight data
            
            % For each of the predefined contrast types (see above)...
            for c = 1:length(con)
                switch con{c}
                    case 'avgDist' % just average across all paired distances
                        gidx    = 1:465;
                    case '1finger' 
                        gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',1); 
                    case '2finger' 
                        gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',2);
                    case '3finger'
                        gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',3);
                    case '4finger' 
                        gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',4);
                    case 'numDigits'
                        gidx    = pp1_imana('GET_idxAcrossNumDigits','glm',glm);
                    otherwise
                        error('unknown contrast');
                end
                % avg. distances according to contrast selected
                Y.LDC   = vdat(:,:,:,gidx);
                Y.LDC   = ssqrt(Y.LDC);
                Y.LDC   = nanmean(Y.LDC,4); 
                % prep output file
                Y.dim   = vol(1).dim;
                Y.dt    = vol(1).dt;
                Y.mat   = vol(1).mat;    
                % save output
                Y.fname   = sprintf('%s_cVolSearch_glm%d_%s_%s.nii',subj_name{sn},glm,con{c},name);
                Y.descrip = sprintf('exp: ''pp1_fmri'' \nglm: %d \ncontrast: ''%s''',glm,con{c});

                spm_write_vol(Y,Y.LDC);
                fprintf('Done %s \n',Y.fname)
                clear Y
            end
            clear vol vdat LDC Y
        end
    case 'SUIT_vol2SUIT'
        % take labelled volume and reslice to SUIT space
        % Resliced file output named: 'wd<filename>'
        sn  = 1;
        glm = 2;
        con = 'avgDist';
        vararginoptions(varargin,{'sn','con','glm'});
        % get subject LDC volumetric file
        job2.subj(1).affineTr  = {fullfile(cerebAnatDir,subj_name{sn},['Affine_' subj_name{sn} '_anatomical_seg1.mat'])};
        job2.subj(1).flowfield = {fullfile(cerebAnatDir,subj_name{sn},['u_a_' subj_name{sn} '_anatomical_seg1.nii'])};
        job2.subj(1).resample  = {fullfile(glmDir{glm},subj_name{sn},sprintf('%s_cVolSearch_glm%d_%s_LDC.nii',subj_name{sn},glm,con))};
        job2.subj(1).mask      = {fullfile(cerebAnatDir,subj_name{sn},['c_' subj_name{sn} '_anatomical_pcereb_corr.nii'])}; % corrected version
        suit_reslice_dartel(job2);
        fprintf('Done. Check resliced volume for %s.\n',subj_name{sn});
    case 'SUIT_plotflatmap'
        % this function takes any labelled volume (already in SUIT space)
        % and plots to the surface
        sn  = 1;
        glm = 2;
        con = 'avgDist';
        vararginoptions(varargin,{'sn','con','glm'});
        % get subject LDC volumetric file
        volFile = fullfile(glmDir{glm},subj_name{sn},sprintf('wd%s_cVolSearch_glm%d_%s_LDC.nii',subj_name{sn},glm,con));
        map = colormap('parula');
        data    = suit_map2surf(volFile,'space','SUIT');
        suit_plotflatmap(data,'threshold',0.01,'cscale',[0.01 0.05],'cmap',map);
        
%         volFile   = fullfile(glmDir{glm},subj_name{sn},sprintf('wd%s_cVolSearch_glm%d_%s_LDC.nii',subj_name{sn},glm,con));
%         Vo        = spm_vol(volFile);
%         Vi        = spm_read_vols(Vo);
%         Vv{1}.dat = Vi;
%         Vv{1}.dim = Vo.dim;
%         Vv{1}.mat = Vo.mat;
%         % map file to SUIT surface and save
%         outFile = fullfile(caretDir,subj_name{sn},sprintf('%s_SUIT_glm%d_%s_LDC.nii',subj_name{sn},glm,con));
%         M       = caret_suit_map2surf(Vv,'space','SUIT');
%         caret_save(outFile,M,'ignore_zeros',1);
        
    case '0' % ------------ SC: Subcortical analyses. ---------------------
        % The SC cases are the subcortical (mostly thalamus) analyses.
        % Outputs for segmentations are saved in the sc_anatomicals folder.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'notused_ANAT_native2mni'
        % run SPM jobman to align anatomical to MNI space 
        vararginoptions(varargin,{'sn'});
        
        spm('defaults','fmri');
        spm_jobman('initcfg');
        for s=sn
            cd(fullfile(anatomicalDir,subj_name{s}));
            anat                = fullfile(anatomicalDir,subj_name{s},sprintf('%s_anatomical.nii,1',subj_name{s}));
            J.subj.vol          = {anat};
            J.subj.resample     = {anat};
            J.eoptions.biasreg  = 0.0001;
            J.eoptions.biasfwhm = 60;
            J.eoptions.tpm      = {'/Users/sarbuckle/Documents/MotorControl/matlab/spm12/tpm/TPM.nii'};
            J.eoptions.affreg   = 'mni';
            J.eoptions.reg      = [0 0.001 0.5 0.05 0.2];
            J.eoptions.fwhm     = 0;
            J.eoptions.samp     = 3;
            J.woptions.bb       = [-78 -112 -70; 78 76 85];
            J.woptions.vox      = [1 1 1]; % 1mm isotropic atlas
            J.woptions.interp   = 4;
            
            matlabbatch{1}.spm.spatial.normalise.estwrite = J;
            spm_jobman('run',matlabbatch);
            fprintf('MNI alignment for %s done\n',subj_name{s});    
        end
        % creates MNI-aligned anatomical ws01_anatomical.nii and
        % transformation matrix y_s01_anatomical.nii -> used for coreg      
    case 'SC_segment'                  
        % uses FSL to do subcortical segmentation
        % option3 - space: native / FNIRT or SPM space
        % need to run MATLAB through terminal
        % command: /Applications/MATLAB_R2015b.app/bin/matlab
        sn = 1;
        vararginoptions(varargin,{'sn'});
        space = 'native'; % native or mni
        switch(space)
            case 'native' % subject space         
                prefix = '';
            case 'mni' % mni aligned space
                prefix = 'w';
        end
        % run FSL routine
        for s = sn
            IN     = fullfile(anatomicalDir, subj_name{s}, [prefix subj_name{s}, '_anatomical.nii']);
            outDir = fullfile(scAnatDir, [prefix subj_name{s}]);
            OUT    = fullfile(outDir, [prefix subj_name{s},'_SC.nii']);
            dircheck(outDir);
            % calc with FSL
            cmd = sprintf('run_first_all -i %s -o %s', IN, OUT);
            fprintf('%s\n',cmd);
            %[status,result] = call_fsl(cmd);
            [status,result] = system(cmd);
            if status; error(result); end
            % extract .nii.gz to .nii
            fname = fullfile(outDir,[prefix subj_name{s},'_SC_all_fast_firstseg.nii']);
            cmd   = sprintf('mri_convert %s.gz %s', fname, fname);
            [status,result] = system(cmd);
            if status; error(result); end
            fprintf('Done %s.\n',subj_name{s});
        end
    case 'SC_makeRegionVols'
        % Uses output from 'SC_segment' to map specified subcortical
        % structures to subject's anatomical (in subject space).
        % Saves each roi as vol .nii logical.
        sn = 1;
        vararginoptions(varargin,{'sn'});
        space = 'native';
        switch(space)
            case 'native' % subject space         
                prefix = '';
            case 'mni' % mni aligned space
                prefix = 'w';
        end
        % 10 L Thalamus
        % 11 L Caudate
        % 12 L Putamen
        % 13 L Pallidum
        % 49 R Thalamus
        % 50 R Caudate
        % 51 R Putamen
        % 52 R Pallidum
        roi = [10:13; 49:52];
        
        for s = sn
            IN  = fullfile(scAnatDir, [prefix subj_name{s}], [prefix subj_name{s},'_SC_all_fast_firstseg.nii']);
            outDir = fullfile(scAnatDir, [prefix subj_name{s}]);
            for h = 1:2
                for r = roi(h,:)
                    OUT = fullfile(outDir,sprintf('%s_%d_%d.nii',subj_name{s},h,r));
                    spm_imcalc(IN,OUT,['i1==' num2str(r)]);
                end
            end
        end
        
    case '0' % ------------ GLM: SPM GLM fitting. Expand for more info. ---
        % The GLM cases fit general linear models to subject data with 
        % SPM functionality.
        %
        % All functions can be called with ('GLM_processAll','sn',[Subj#s]).
        %
        % You can view reconstructed surfaces with Caret software.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'WRAPPER_GLM'                                                  
        vararginoptions(varargin,{'sn','glm','sess'});
        % You can call this case to do all the GLM estimation and contrasts.
        for s = sn
            for g = glm
                pp1_imana('GLM_make','sn',s,'glm',g);
                pp1_imana('GLM_estimate','sn',s,'glm',g);
                if g==1
                    pp1_imana('GLM_contrastglm1','sn',s);
                elseif g>1
                    if g==2
                        pp1_imana('GLM_contrastglm2','sn',s);
                    elseif g==3
                        pp1_imana('GLM_contrastglm3','sn',s);
                    end
                    pp1_imana('PSC_calcChord','sn',s,'glm',g);
                end
            end
        end
    case 'GLM_make'                                                         % STEP 3.1   :  Make the SPM.mat and SPM_info.mat files (prep the GLM)
        % makes the GLM file for each subject, and a corresponding aux.
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        % housekeeping
        prefix		 = dataPrefix{sn};
        T			 = [];
        % Gather appropriate GLM presets.
        [~,d_hrf]      = spm_hrf(TR_length); % default hrf params
        subj_hrfParams = {[3.6873 11.893 d_hrf(3) d_hrf(4) 0.23299 d_hrf(6) d_hrf(7)],... % subject-specific hrf params (estimated from glm1)
                          [5.5018 15.717 d_hrf(3) d_hrf(4) 6.2265 d_hrf(6) d_hrf(7)],...
                          [4.986  15.452 d_hrf(3) d_hrf(4) 6.1327 d_hrf(6) d_hrf(7)],...
                          [5.0097 16.379 d_hrf(3) d_hrf(4) 5.6836 d_hrf(6) d_hrf(7)],...
                          [4.9406 13.158 d_hrf(3) d_hrf(4) 2.6733 d_hrf(6) d_hrf(7)],...
                          [4.1708 12.603 d_hrf(3) d_hrf(4) 4.6819 d_hrf(6) d_hrf(7)]};
        % Define number of regressors in glm
        switch glm 
            case 1
                % model all conditions together and use to optimize hrf fit
                hrf_params = d_hrf; % use defaults
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 1; % to optimize hrf fits
            case 2
                % model all chords, don't exclude error trials
                % define the 7 parameters of the HRF
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 31; % 31 chords
            case 3
                % model all chords (31) and a regressor for thumb movements
                % (reg number 32)
                % define the 7 parameters of the HRF
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 32; % 31 chords
        end
        
        % Load subject's .dat file (has info on each trial)
        D = dload(fullfile(behavDir,sprintf('pp1_fmri_%s.dat',subj_name{sn})));
        D.onset = D.startTimeMeas./1000 - TR_length*numDummys + D.cueTime./1000; % correct timing for removed dummy scans
        % Do some subject structure fields.
        J.dir 			 = {fullfile(glmDir{glm}, subj_name{sn})};
        J.timing.units   = 'secs';                                          % timing unit that all timing in model will be
        J.timing.RT 	 = TR_length;                                       % TR (in seconds, as per 'J.timing.units')
        J.timing.fmri_t  = 16;
        J.timing.fmri_t0 = 1;
        % Make glm model job structure 
        subjDir = fullfile(imagingDir,subj_name{sn});
        N       = findFiles(subjDir,[prefix subj_name{sn} '_run']); % find run files for subject (organized in order from run 1 to X)
        runs = unique(D.BN)';
        for r = runs
            % get image names for run
            numImgs = (numTRs{sn}(r)-numDummys);
            Q = {};
            for i = 1:numImgs
                % algin each img within 4d nifti to coregistered epi
                Q{i} = fullfile(subjDir,[N{r},',',num2str(i)]);
            end
            J.sess(r).scans = Q; % images to model for this run
            % Loop through conditions.
            R = getrow(D,D.BN==r);
            for c = 1:numConds
                switch glm
                    case 1
                        idx	= logical(R.chordNum>0);     
                        J.sess(r).cond(c).name     = 'stim_on';
                        J.sess(r).cond(c).onset    = R.onset(idx);    
                        J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
                    case 2
                        % Model chords regressors
                        % include all trials regardless of judgement
                        idx	= find(R.chordNum==c); % find indx of all trials in run of that condition 
                        J.sess(r).cond(c).name     = sprintf('chord_%d',R.chordNum(idx(1)));  % make condition name (for user readability)
                        J.sess(r).cond(c).onset    = R.onset(idx);    
                        J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
                        S.isError = 0;
                    case 3
                        S.isError = 0;
                        if c <32
                            % Model chords regressors
                            % include all trials regardless of judgement
                            idx	= find(R.chordNum==c); % find indx of all trials in run of that condition 
                            J.sess(r).cond(c).name     = sprintf('chord_%d',R.chordNum(idx(1)));  % make condition name (for user readability)
                            J.sess(r).cond(c).onset    = R.onset(idx);    
                            J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
                        elseif c==32
                            % Model thumb presses
                            idx = find(R.RT>0);
                            J.sess(r).cond(c).name     = 'thumb_move';  % make condition name (for user readability)
                            J.sess(r).cond(c).onset    = R.onset(idx) + R.RT(idx)/1000;    
                            J.sess(r).cond(c).duration = 1;
                        end
                end
                J.sess(r).cond(c).tmod = 0;
                J.sess(r).cond(c).orth = 0;
                J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                % Do some subject info for fields in SPM_info.mat.
                S.sn  = sn;
                S.run = r;
                if (c==1 && glm==1)
                    S.chord 		= 0;                          
                    S.numDigits     = 0; 
                    S.targetForce   = R.numStim(idx(1));
                    S.numStim       = R.numStim(idx(1));
                elseif c<32
                    S.chord 		= R.chordNum(find(R.chordNum==c,1));                       
                    S.numDigits     = R.numDigits(find(R.chordNum==c,1)); 
                    S.targetForce   = R.targetForceStim(find(R.chordNum==c,1));
                    S.numStim       = R.numStim(find(R.chordNum==c,1));
                elseif c
                    S.chord 		= 32;                          
                    S.numDigits     = 0; 
                    S.targetForce   = 0;
                    S.numStim       = 0;
                end  
                S.tt	  = c;
                S.regtype = 'avgTask';
                T		  = addstruct(T,S);
            end
            % Add any additional aux. regressors here.
            J.sess(r).multi 	= {''};
            J.sess(r).regress 	= struct('name', {}, 'val', {});
            J.sess(r).multi_reg = {''};                                
            % Define high pass filter cutoff (in seconds): see glm cases.
            J.sess(r).hpf 		= hrf_cutoff;
        end
        J.fact 			   = struct('name', {}, 'levels', {});
        J.bases.hrf.derivs = [0 0];
        J.bases.hrf.params = hrf_params;
        J.volt 			   = 1;
        J.global 		   = 'None';
        J.mask 	           = {fullfile(baseDir, 'imaging_data',subj_name{sn}, 'rmask_noskull.nii')};
        J.mthresh 		   = 0.05;
        J.cvi_mask 		   = {fullfile(baseDir, 'imaging_data',subj_name{sn},'rmask_gray.nii')};
        J.cvi 			   = cvi_type;
        % Save the GLM file for this subject.
        spm_rwls_run_fmri_spec(J);
        % Save the aux. information file (SPM_info.mat).
        % This file contains user-friendly information about the glm
        % model, regressor types, condition names, etc.
        save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
    case 'GLM_estimate'                                                     % STEP 3.2   :  Run the GLM according to model defined by SPM.mat
        % Estimate the GLM from the appropriate SPM.mat file. 
        % Make GLM files with case 'GLM_make'.
        vararginoptions(varargin,{'sn','glm'});
        % Load files
        load(fullfile(glmDir{glm},subj_name{sn},'SPM.mat'));
        SPM.swd = fullfile(glmDir{glm},subj_name{sn});
        % Run the GLM.
        spm_rwls_spm(SPM);
        % for checking -returns img of head movements and corrected sd vals
        % spm_rwls_resstats(SPM)
    case 'GLM_contrastglm1'
        % Make t-stat contrasts for single task regressor (localizer)
        % enter sn, glm #
        % 1: task vs. rests
        vararginoptions(varargin,{'sn'});
        glm = 1;
        cwd = pwd;
        % Go to subject's directory
        cd(fullfile(glmDir{glm}, subj_name{sn}));
        load SPM;
        SPM = rmfield(SPM,'xCon');
        T   = load('SPM_info.mat');

        con                = zeros(1,size(SPM.xX.X,2));
        con(:,T.tt>0)      = 1;
        con                = con/sum(con);
        SPM.xCon(1)        = spm_FcUtil('Set','task', 'T', 'c',con',SPM.xX.xKXs);

        %____do the constrast
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save('SPM.mat','SPM');
        cd(cwd);
    case 'GLM_contrastglm2'                                                 % STEP 3.3   :  Make t-stat contrasts for specified GLM estimates.
        % Make t-stat contrasts for specified GLM estimates.
        % enter sn, glm #
        % models chords, no error trials modeled in glm2 so no contrast
        vararginoptions(varargin,{'sn'});
        glm = 2;
        cwd = pwd;
        % Go to subject's directory
        cd(fullfile(glmDir{glm}, subj_name{sn}));
        load SPM;
        SPM = rmfield(SPM,'xCon');
        T   = load('SPM_info.mat');

        %_____t contrast for chords
        for d = 1:31
            con                = zeros(1,size(SPM.xX.X,2));
            con(:,T.tt==d)     = 1;
            con                = con/sum(con);
            SPM.xCon(d)        = spm_FcUtil('Set',sprintf('chord_%d',d), 'T', 'c',con',SPM.xX.xKXs);
        end;

        %_____t contrast overall chords
        con                    = zeros(1,size(SPM.xX.X,2));
        con(:,T.tt>0)          = 1;
        con                    = con/sum(con);
        SPM.xCon(32)            = spm_FcUtil('Set',sprintf('overall'), 'T', 'c',con',SPM.xX.xKXs);

        %____do the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save('SPM.mat','SPM');
        cd(cwd);
    case 'GLM_contrastglm3'                                                 % STEP 3.3   :  Make t-stat contrasts for specified GLM estimates.
        % Make t-stat contrasts for specified GLM estimates.
        % enter sn, glm #
        % models each chord and also error trials
        vararginoptions(varargin,{'sn'});
        glm = 3;
        % Go to subject's directory
        cd(fullfile(glmDir{glm}, subj_name{sn}));
        load SPM;
        SPM = rmfield(SPM,'xCon');
        T   = load('SPM_info.mat');
        %_____t contrast for chords & thumb movement
        for d = 1:32
            con               = zeros(1,size(SPM.xX.X,2));
            con(:,T.chord==d) = 1;
            con               = con/sum(con);
            if d<32 % chords
                SPM.xCon(d) = spm_FcUtil('Set',sprintf('chord_%d',d), 'T', 'c',con',SPM.xX.xKXs);
            else % thumb movement
                SPM.xCon(d) = spm_FcUtil('Set','thumb_response', 'T', 'c',con',SPM.xX.xKXs);
            end
        end

        %_____t contrast overall chords
        con                    = zeros(1,size(SPM.xX.X,2));
        con(:,T.chord>0 & T.chord<32) = 1;
        con                    = con/sum(con);
        SPM.xCon(33)           = spm_FcUtil('Set',sprintf('overall'), 'T', 'c',con',SPM.xX.xKXs);

        %____do the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save('SPM.mat','SPM');
    case 'PSC_calcChord'                                                    % calculates % signal change for chords
        % calculate psc for all digits vs. rest - based on betas from glm 1    
        sn  = 1;
        glm = 2;
        vararginoptions(varargin,{'sn','glm'});
        if glm==2
            numImgs = 31;
        elseif glm==3
            numImgs = 32;
        end
        % assumes first 31 contrasts are for the chords
        for s=sn
            cd(fullfile(glmDir{glm}, subj_name{s}));
            load SPM;
            T = load('SPM_info.mat');
            X = (SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
            h = median(max(X));               % Height of response is defined as median of max regressor height (across conds and runs) for this subject
            P = {};                           % Filenames of input images
            numB = length(SPM.xX.iB);         % Partitions - runs
            for p = SPM.xX.iB
                P{end+1} = sprintf('beta_%4.4d.nii',p);       % get the intercepts (for each run) and use them to calculate the baseline (mean images) * max height of design matrix regressor
            end
            for con = 1:numImgs   % 31 chords + other regressors 
                P{numB+1} = sprintf('con_%04d.nii',con);
                outname   = sprintf('psc_%02d.nii',con); % ,subj_name{s}
                formula   = '100.*%f.*i%1.0f./((';
                % construct formula dynamically incase numRuns changes
                for i = 1:numB
                    if i~=numB
                        fadd = sprintf('i%1.0f+',i);
                    else
                        fadd = sprintf('i%1.0f)/',i);
                    end
                    formula = [formula fadd];
                end
                formula = [formula num2str(numB) ')'];
                formula = sprintf(formula,h,numB+1);

                spm_imcalc_ui(P,outname,formula,{0,[],spm_type(16),[]});        % Calculate percent signal change
            end
            fprintf('Subject %d: %3.3f\n',s,h);
        end
        
    case '0' % ------------ SEARCH: searchlight analyses. Expand for more info.  REQUIRES FURTHER EDITING in group_cSPM (editing for comments on what is happening)!!!
        % The SEARCH cases are used to conduct surface-searchlight analyses 
        % using the rsa toolbox from JDiedrichsen and NEjaz (among others).
        %
        % All functions can be called with ('SEARCH_processAll','sn',[Subj#s]).
        %
        % You can view reconstructed surfaces with Caret software.
        %
        % The contrast metrics calculated from the full condition
        % searchlight are used to estimate boundaries for the rois.
        % The values used to estimate boundaries are the avg. paired
        % distances.
        %
        % See blurbs in each SEARCH case to understand what they do.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
%     case 'WRAPPER_SEARCH'                                               
%         glm = 2;
%         vararginoptions(varargin,{'sn','glm'});
%         % You can call this case to do all the searchlight analyses.
%         % 'sn' can be an array of subjects because each processing case
%         % contained within loops through the subject array submitted to the
%         % case.
%         for s = sn 
%             for g = glm
%                 pp1_imana('SEARCH_define','sn',s,'glm',g);
%                 pp1_imana('SEARCH_run_LDC','sn',s,'glm',g);
%                 pp1_imana('SEARCH_mapSurfLDC','sn',s,'glm',glm);
%             end
%         end
%         %fmri_imana('SEARCH_group_make');   % require group data
%         %fmri_imana('SEARCH_group_cSPM');
% %     case 'SEARCH_define'                                                    % STEP 4.1   :  Defines searchlights for 120 voxels in grey matter surface
%         glm = 2;
%         vararginoptions(varargin,{'sn','glm'});
%         
%         mask       = fullfile(glmDir{glm},subj_name{sn},'mask.nii');
%         Vmask      = spm_vol(mask);
%         Vmask.data = spm_read_vols(Vmask);
% 
%         LcaretDir = fullfile(caretDir,sprintf('x%s',subj_name{sn}),'LeftHem');
%         RcaretDir = fullfile(caretDir,sprintf('x%s',subj_name{sn}),'RightHem');
%         white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
%         pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
%         topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};
%         S         = rsa_readSurf(white,pial,topo);
% 
%         L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[15 120]);
%         save(fullfile(anatomicalDir,subj_name{sn},sprintf('%s_searchlight_120.mat',subj_name{sn})),'-struct','L');
%     case 'SEARCH_run_LDC'                                                   % STEP 4.2   :  Runs LDC searchlight using defined searchlights (above)
%         % Requires java functionality unless running on SArbuckle's
%         % computer.
%         glm  = 2;
%         vararginoptions(varargin,{'sn','glm'});
%         
%         if glm==2
%             numConds = 31; % no errors modeled
%         elseif glm==3
%             numConds = 32; % error modeled
%         elseif glm==4
%             numConds = 32; % foot movement modeled
%         end
%         
%         block = 5e7;
%         cwd   = pwd;                                                        % copy current directory (to return to later)
%         % go to subject's glm directory 
%         cd(fullfile(glmDir{glm},subj_name{sn}));
%         % load their searchlight definitions and SPM file
%         L = load(fullfile(anatomicalDir,subj_name{sn},sprintf('%s_searchlight_120.mat',subj_name{sn})));
%         load SPM;
%         SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{sn}));
%         name = sprintf('%s_glm%d',subj_name{sn},glm);
%         % make index vectors
%         conditionVec  = kron(ones(numel(SPM.Sess),1),[1:numConds]');
%         partitionVec  = kron([1:numel(SPM.Sess)]',ones(numConds,1));
%         % run the searchlight
%         rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partitionVec,'analysisName',name,'idealBlock',block);
%         cd(cwd);
    case 'SEARCH_mapSurfLDC'                                                % STEP 4.3   :  Averaged LDC values for specified contrasts
        % Calls 'MISC_SEARCH_calculate_contrast'
        sn  = 1;
        glm = 1;
        con = {'avg'};
        vararginoptions(varargin,{'sn','glm','con'});
        % Use 'con' option to define different contrasts.
        %   'avg'    :  Average LDC nii for all 20 conds,
        %                invariant of speed (avg of 10 pairwise distances)

        % Load subject surface searchlight results (1 vol per paired conds)
        LDC_file            = fullfile(glmDir{glm},subj_name{sn},sprintf('%s_glm%d_LDC.nii',subj_name{sn},glm)); % searchlight nifti
        [subjDir,fname,ext] = fileparts(LDC_file);
        cd(subjDir);
        vol  = spm_vol([fname ext]);
        vdat = spm_read_vols(vol); % is searchlight data
        % For each of the predefined contrast types (see above)...
        for c = 1:length(con)
            switch con{c}
                case 'avg' % just average across all paired distances
                    gidx    = 1:465;
                case '1digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',1); 
                case '2digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',2);
                case '3digitChords'
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',3);
                case '4digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',4);
                case 'numDigits'
                    gidx    = pp1_imana('GET_idxAcrossNumDigits','glm',glm);
            end
            % avg. distances according to contrast selected
            Y.LDC   = vdat(:,:,:,gidx);
            Y.LDC   = ssqrt(Y.LDC);
            Y.LDC   = nanmean(Y.LDC,4); 
            % prep output file
            Y.dim   = vol(1).dim;
            Y.dt    = vol(1).dt;
            Y.mat   = vol(1).mat;    
            % save output
            Y.fname   = sprintf('%s_glm%d_%sLDC.nii',subj_name{sn},glm,con{c});
            Y.descrip = sprintf('exp: ''pp1'' \nglm: ''FAST'' \ncontrast: ''%s''',con{c});

            spm_write_vol(Y,Y.LDC);
            fprintf('Done %s\n',Y.fname);

            clear Y
        end
    case 'SEARCH_map_nii'                                                   % STEP 4.4   :  Map searchlight results (.nii) onto surface (.metric)
        % map volume images to metric file and save them in individual surface folder
        sn  = 1;
        glm = 1;
        con = {'avg','1digitChords','2digitChords','3digitChords','4digitChords','numDigits'}; % does all con imgs as default
        vararginoptions(varargin,{'sn','con','glm'});
        % 'con' option defines each contrast.
        %   'avg'    :  Average LDC nii for all 20 conds

        hemisphere = 1:2;
        for s = sn
            for c = 1:length(con)
                ctype = con{c};
                for h=hemisphere
                    caretSDir = fullfile(caretDir,[atlasA,subj_name{s}],hemName{h});
                    white     = caret_load(fullfile(caretSDir,[hem{h} '.WHITE.coord']));
                    pial      = caret_load(fullfile(caretSDir,[hem{h} '.PIAL.coord']));
                    images    = fullfile(glmDir{glm},subj_name{s},sprintf('%s_glm%d_%sLDC.nii',subj_name{s},glm,ctype));
                    outfile   = sprintf('%s_%sfunc_%d.metric',subj_name{s},ctype,glm);
                    M         = caret_vol2surf_own(white.data,pial.data,images,'ignore_zeros',1);
                    caret_save(fullfile(caretSDir,outfile),M);
                    fprintf('Done subj %d con %s hemi %d \n',s,ctype,h)
                end
            end
        end
    case 'SEARCH_group_make'                                                % STEP 4.5   :  Make group metric files by condensing subjec contrast metric files
        % Calculate group metric files from the searchlight results. 
        % Takes the 3 contrast results ('avg','speed', & 'digit') across
        % subjects and makes a group level metric file that contains each
        % subject's data for that contrast type.
        
        % Some presets
        INname     = {'avgfunc_3','speedfunc_3','digitfunc_3'};
        OUTname    = {'group_avg_3','group_speed_3','group_digit_3'};
        inputcol   = [1 1 1];
        replaceNaN = [0 0 0];     
        % Loop over hemispheres.
        for h = 1:2
            % Go to the directory where the group surface atlas resides
            surfaceGroupDir = [caretDir filesep atlasname filesep hemName{h}];
            cd(surfaceGroupDir);
            % Loop over each input metric file in 'INname' and make a group metric file
            for j = 1:length(INname); 
                % Loop over subjects...
                for i = 1:length(subj_name); 
                    % ...and define the names of their metric files
                    infilenames{j}{i} = [caretDir filesep atlasA subj_name{i} filesep hemName{h} filesep subj_name{i} '_' INname{j} '.metric'];
                end;
                % Name the output filename for this group metric file in average surface folder
                outfilenames{j} = [surfaceGroupDir filesep hem{h} '.' OUTname{j} '.metric'];
                % Finally, make the group metric file for this metric type/contrast
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
                % Verbose display to user
                fprintf('hem: %i  image: %i \n', h,j);
            end;
        end;
    case 'SEARCH_group_cSPM'                                                % STEP 4.6   :  Generate a statistical surface map (onesample_t test) from smoothed group metric files. Also avgs. distances across subjs.
        % Calculate group stats files from the group metric files. 
        % Takes the 3 group metric files ('avg','speed', & 'digit') and 
        % calculates group level stats (one sample t test that the mean 
        % effect is bigger than zero). 
        % 
        % Although we calculate the t-score, corresponding p-value for that
        % t-score, and finally the z-score 
        % subject's data for that contrast type.
        % ******* finish blurb here and comments + styling code below.
        
        sn = 1:9;
        SPMname={'group_avg_3','group_speed_3','group_digit_3'};
        
        sqrtTransform=[1,1,1]; % Should you take ssqrt before submitting? 
                                % Yes, b/c used rsa.distanceLDC to
                                % calculate distances. This function
                                % returns squared cv mahalanobis distance.
        SummaryName = '.summary.metric';
        hemi = [1 2];
        
        for h=hemi
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
            %----get the full directory name of the metric files and the NONsmoothed metric files that we create below
            for i=1:length(SPMname);
                %sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' SPMname{i} '.metric']; % smoothed
                sfilenames{i}=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '.metric']; % no smoothing
            end;
            %----loop over the metric files and calculate the cSPM of each with the non-smoothed metrics
            for i=1:length(SPMname);
                Data=caret_load(sfilenames{i});
                if sqrtTransform(i)
                    Data.data=ssqrt(Data.data);
                end;
                cSPM=caret_getcSPM('onesample_t','data',Data.data(:,sn),'maskthreshold',0.5); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                caret_savecSPM([surfaceGroupDir filesep hem{h} '.' SPMname{i} '_stats.metric'],cSPM);
                save([surfaceGroupDir  filesep   'cSPM_' SPMname{i} '.mat'],'cSPM');
                data(:,i)=cSPM.con(1).con; % mean
                data(:,i+length(SPMname))=cSPM.con(1).Z; % T
                column_name{i}=['mean_' SPMname{i}];
                column_name{i+length(SPMname)}=['T_' SPMname{i}];
            end;
            C = caret_struct('metric','data',data,'column_name',column_name);
            caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C);
        end;
        fprintf('Done \n')
    case 'SEARCH_vol2Surf'                                                  % STEP 4.4   :  Map searchlight results (.nii) onto surface (.metric)
        % map volume images to metric file and save them in individual surface folder
        sn  = 1;
        glm = 2;
        %con = [1:31];
        hemisphere = 1:2;
        for s = sn
            for h=hemisphere
                outfile   = sprintf('%s_%s_glm%d.metric',subj_name{s},hemName{h},glm);
                caretSDir = fullfile(caretDir,[atlasA,subj_name{s}],hemName{h});
                white     = caret_load(fullfile(caretSDir,[hem{h} '.WHITE.coord']));
                pial      = caret_load(fullfile(caretSDir,[hem{h} '.PIAL.coord']));
                images    = {};
                % image names for searchlight analyses (cols 1-4)
%                 for c = 1:length(con)
%                     ctype = con{c};
%                     images{c} = fullfile(glmDir{glm},subj_name{s},sprintf('s0%d_glm%d_%s_LDC.nii',s,glm,ctype));
%                 end
                % image names flex > ext and ext > flex contrasts (cols 5 & 6)
                for c = 1:31
                    images{c} = fullfile(glmDir{glm},subj_name{s},sprintf('spmT_00%02d.nii',c));
                end
                % map volumes to surface
                M = caret_vol2surf_own(white.data,pial.data,images,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,outfile),M);
                fprintf('Done subj %02d\n',s)
            end
        end
        
    case '0' % ------------ ROI: roi analyses. Expand for more info. ------
        % The ROI cases are used to:
        %       - map ROIs to each subject
        %       - harvest timeseries from each roi for each condition
        %       - harvest activity patterns (i.e. beta weights for each voxel 
        %          in roi of that subject)
        %       - conduct statistical analyses on activity patterns and distances
        %       - assess pattern consistencies (for each subject for a glm)
        %       - assess reliability of distance estimates across subejcts
        %
        % There is no 'processAll' case here. However, the following cases
        % must be called to utilize other roi cases:
        %       'ROI_makePaint'   :  Creates roi paint files (see case)
        %       'ROI_define'      :  Maps rois to surface of each subject-
        %                             requires paint files from above case.
        %       'ROI_timeseries'  :  Only if you wish to plot timeseries
        %       'ROI_getBetas'    :  Harvest patterns from roi
        %       'ROI_stats'       :  Estimate distances, etc. 
        %                               This is the big kahuna as it is
        %                               often loaded by future cases.
        %
        % You can view roi maps by loading paint files and subject surfaces
        % in Caret (software).
        % 
        % Most functionality is achieved with rsa toolbox by JDiedrichsen
        % and NEjaz (among others).
        %
        % See blurbs in each SEARCH case to understand what they do.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    case 'ROI_makePaintFile' % defines rostral & caudal M1, ba 3a, 3b, 1, and 2, all cut to hand areas.                                                   
    % Creates ROI boundaries on the group template/atlas (fsaverage_sym) for rostral and caudal M1.
    % ROIs are defined using:
    %       - probabilistic atlas
    %       - boundaries for 4 major lobes
    %       - flatmap coordinates for X Y coordinate restriction 
    
    for h=1:2 % loop over hemispheres
        
        % load previously defined M1 S1 rois, and probabilistic atlas (both
        % found in caret group directory).
        groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
        cd(groupDir);
        R  = caret_load('ROI.paint');                  
		M  = caret_load([hem{h} '.propatlas.metric']); 
		
        % assign roi labels to verticies
        Ba1 = M.data(:,1); 
        Ba2 = M.data(:,2); 
        B3a = M.data(:,3); 
        B3b = M.data(:,4); 
        rM1 = M.data(:,7); % BA4a
        cM1 = M.data(:,8); % BA4p
        
        % cut these rois to fit within bounds of previously defined M1 & S1
        % hand rois.
        S1  = R.data==1;
        M1  = R.data==2; 
        Ba1(~S1) = 0;
        Ba2(~S1) = 0;
        B3a(~S1) = 0;
        B3b(~S1) = 0;
        rM1(~M1) = 0;
        cM1(~M1) = 0;
        
        % vertex is assigned to roi based on max probabilitiy (above a
        % threshold)
        [Prop,ROI]     = max([B3a B3b Ba1 Ba2 rM1 cM1],[],2); 
        ROI(Prop<0.25) = 0;
        
        % save paint files of the rois
        areas = {{'ba3A'},{'ba3B'},{'ba1'},{'ba2'},{'rM1'},{'cM1'}};
        names = {'ba3A','ba3B','ba1','ba2','rM1','cM1'};
        color = [0 153 153;...      % ba3A
                 0 89 171;...       % ba3B
                 89 122 189;...     % ba1
                 140 186 235;...    % ba2
                 230 0 0;...        % rM1
                 255 153 0];        % cM1
            
        Paint = caret_struct('paint','data',ROI,'paintnames',names,'column_name',{'ROI'});
        caret_save(['ROI_pp1.paint'],Paint);
        caret_combinePaint('ROI_pp1.paint','ROI_pp1.paint','ROI_pp1.areacolor',...
            'areas',areas,'names',names,'colors',color);
        
    end    
    case 'ROI_makeBApaint'                                                  % Make paint file for brodmann area ROIs cut to hand area (saves as ROI_3.paint)
    % Modified from df1_imana. (in current state this is ugly..)
    % Creates ROI boundaries on the group template/atlas (fsaverage).
    % ROIs are defined using:
    %       - probabilistic atlas
    %       - boundaries for 4 major lobes
    %       - group surface searchlight results
    %       - flatmap coordinates for X Y coordinate restriction 
    
    for h=1:2 % loop over hemispheres
        % - - - - - - - - - - - - Load some files - - - - - - - - - - - - -
        groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
        cd(groupDir);
        %P    = caret_load([hem{h} '.lobes.paint']);                  % paint file of 4 major lobes (1 = col 2 of paintnames)
        %C    = caret_load([hem{h} '.FLAT.coord']);                   % Caret flatmap (for X Y coordinate restrictions)
		M = caret_load([hem{h} '.propatlas.metric']);             % probabilistic atlas
		Defined_rois = caret_load('ROI.paint');                     % ROI.paint file from probabalistic atlas- file has nice, standard ROIs for motor, sensory, aux. rois, and some visual rois
        % - - - - - - - - - - Assign roi labels to vertices - - - - - - - -
        % get data for new rois (from probability atlas)
        % renumber rois so they are arranged in anatomical order (3a,b,1,2)    
        Ba1 = M.data(:,1); %logical(M.data(:,1)).*3;
        Ba2 = M.data(:,2); %logical(M.data(:,2)).*4;
        B3a = M.data(:,3); %logical(M.data(:,3)).*1;
        B3b = M.data(:,4); %logical(M.data(:,4)).*2;
        
        % coordinate is ROI w/ for which it has greatest associated probability
        [Prop3,ROI3] = max([B3a B3b Ba1 Ba2],[],2);
        
        % Define ROIS with:
        %...cytoarchitectonic prob (>0.2) - - - - - - - - - - - - - - - - -
            ROI3(Prop3<0.2) = 0;
        %...and cut to hand-knob size (mask with S1 roi)- - - - - - - - - -
            ROI3(ROI3==1 & Defined_rois.data~=1)   = 0;  % Ba3a :  must be within S1 roi already
            ROI3(ROI3==2 & Defined_rois.data~=1)   = 0;  % Ba3b :  must be within S1 roi already
            ROI3(ROI3==3 & Defined_rois.data~=1)   = 0;  % Ba1  :  must be within S1 roi already
            ROI3(ROI3==4 & Defined_rois.data~=1)   = 0;  % Ba2  :  must be within S1 roi already
            

        % Save paint file for Ba1:3b 
        clear areas colors names Paint
        areas{1}{1} = {'sB3a'}; 
        areas{2}{1} = {'sB3b'}; 
        areas{3}{1} = {'sBa1'}; 
        areas{4}{1} = {'sBa2'};
        
        names = {'sB3a','sB3b','sBa1','sBa2'};
        
        colors = [0 100 200; 255 255 0; 150 50 255; 0 153 153]; 
        
        % save paint structure (roi map)
        C = caret_struct('paint','data',ROI3,'paintnames',names,'column_name',{'ROI'});
        caret_save('ROI_3.paint',C);
        % save areacolor (for roi map)
        M = [];
        M.encoding             = {'BINARY'};
        M.column_name          = C.paintnames;
        M.column_color_mapping = repmat([-5 5],4,1);
        M.paintnames           = C.paintnames;
        M.data                 = colors;
        caret_save('ROI_3.areacolor',M);
    end;        
    fprintf('Done.\n');
    
    case 'ROI_define'                                                       % Define rois: BA rois, M1/S1 cut to hand area, BA rois cut to hand area
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 

        I = pp1_imana('LIST_subjs');
        I = getrow(I,I.fmri_sessions==2);
        %vararginoptions(varargin,{'sn'});

        linedef = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
        tmpMaskFile = {};  % for volumetric sub-cortical regions (these tmp masks are deleted after mapping to subject)
        for s = unique(I.sn)'
            
            caretSubjDir = fullfile(caretDir,[atlasA I.origSN{s}]);
            mask         = fullfile(glmDir{1},I.origSN{s},'mask.nii,1');  
            
            for h = 1:2 
                
                D1 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI_pp1.paint'])); % ba3A, ba3B, ba1, ba2, rM1, cM1
                D2 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI.paint']));     % premade ROIs from fsaverage_sym: M1 and S1
                
                for r = 1:numregions
                    if r<7; D = D1; rr = r;         % pp1 rois
                    elseif (r==7 || r==7+numregions); D = D2; rr = 1;    % S1
                    elseif (r==8 || r==8+numregions); D = D2; rr = 2;    % M1
                    elseif (r==9 || r==9+numregions); D = D2; rr = 7;    % SPLa
                    elseif (r==10 || r==10+numregions); D = D2; rr = 8;  % SPLp
                    end
                    idx = r+(h-1)*numregions;
                    % make R region structure for participant
                    R{idx}.name  = [I.origSN{s} '_' regname{r} '_' hem{h}];
                    if r<11 % cortical surface regions
                        R{idx}.type     = 'surf_nodes';
                        R{idx}.location = find(D.data(:,1)==rr);
                        R{idx}.white    = fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                        R{idx}.pial     = fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                        R{idx}.topo     = fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                        R{idx}.linedef  = linedef;
                        R{idx}.flat     = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                        R{idx}.image    = mask; % functional mask
                    elseif r==11 % thalamus volumetric region
                        if h==1; rr = 10; elseif h==2; rr = 49; end % 10 is left thalamus, 49 is right thalamus
                        R{idx}.type     = 'roi_image';
                        R{idx}.file     = fullfile(scAnatDir,I.origSN{s},sprintf('%s_%d_%d.nii',I.origSN{s},h,rr));
                        R{idx}.value    = 1;
                        R{idx}.image    = fullfile(glmDir{1},I.origSN{s},sprintf('mask_tmp%d%d.nii',h,rr)); % functional mask
                        % make temporary functional mask of thalamus
                        % This is done b/c region_calcregions assumes mask
                        % is same dimension as roi image. Thus,
                        % region_calcregions will basically redo the mask,
                        % but that's okay.
                        F      = spm_vol(R{idx}.file);
                        F.data = spm_read_vols(F);
                        M      = spm_vol(mask);
                        M.data = spm_read_vols(M);
                        V(1)   = F;
                        V(2)   = M;
                        spm_imcalc(V,R{idx}.image,'i1.*i2');
                        tmpMaskFile{end+1} = R{idx}.image;
                        clear V
                    end
                end    
            end
            %R = region_calcregions(R,'exclude',[2,6; 5,6],'exclude_thres',0.75);
            exculdeRoi = [1,5; 1,6; 2,5; 2,6; 3,5; 3,6; 7,8];
            R = region_calcregions(R,'exclude',[exculdeRoi; exculdeRoi+numregions],'exclude_thres',0.75);
            cd(regDir);
            save(['regions_' I.origSN{s} '.mat'],'R');
            fprintf('\n %s done\n',I.origSN{s})
            clear R
        end
        for i = 1:numel(tmpMaskFile) 
            delete(tmpMaskFile{i});
        end
    
    case 'ROI_getTimeseries'                                                % (optional) :  Harvest ROI timeseries for specified region.
        % Use this and 'ROI_plot_timeseries' to ensure good GLM fits with
        % measured BOLD in rois.
        
        % Defaults
        sn  = 1;
        glm = 1;
        roi = 1:11;
        vararginoptions(varargin,{'sn','glm','roi'});
        pre  = 4;                                                                 % how many TRs before trial onset (2.8 secs)
        post = 20;                                                                % how many TRs after trial onset (11.2 secs)
        % (2) Load SPM and region.mat files, extract timeseries, save file
        T=[];

        for s=sn
                cd(fullfile(glmDir{glm},subj_name{s}));                                  % cd to subject's GLM dir
                load SPM;
                load(fullfile(regDir,['regions_' subj_name{s} '.mat']));  % load subject's region_define info- variable loaded is R
                for reg=roi
                    [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R{reg});          % get SPM info for voxels contained in specified region
                    D = spmj_get_ons_struct(SPM);                                     % get trial onsets in TRs- because model was in secs, spmj converts onsets to TR #s by dividing time/TR length (320 trials to 4872 TRs)
                    for r=1:size(y_raw,2)
                        for i=1:size(D.block,1);                                    % extract the timeseries of each trial from y_adj, y_hat, & y_res
                            D.y_adj(i,:) = cut(y_adj(:,r),pre,round(D.ons(i)) -1,post,'padding','nan')';
                            D.y_hat(i,:) = cut(y_hat(:,r),pre,round(D.ons(i)) -1,post,'padding','nan')';
                            D.y_res(i,:) = cut(y_res(:,r),pre,round(D.ons(i)) -1,post,'padding','nan')';
                            D.y_raw(i,:) = cut(y_raw(:,r),pre,round(D.ons(i)) -1,post,'padding','nan')';
                        end;
                        D.roi = ones(size(D.event,1),1)*reg;
                        D.sn  = ones(size(D.event,1),1)*s;
                        T=addstruct(T,D);
                    end;
                    fprintf('Done %s (region # %d) for glm %d \n',reg_title{reg},reg,glm);
                end
        end
        save(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)),'-struct','T');
        
        %__________________________________________________________________
    case 'ROI_fitHRF'                                                       % (optional) :  Fit the hrf to the data of a region.
        % p    - parameters of the response function (two Gamma functions)
        %                                                           defaults
        %                                                          (seconds)  
        %        p(1) - delay of response (relative to onset)          6
        %        p(2) - delay of undershoot (relative to onset)       16
        %        p(3) - dispersion of response                         1
        %        p(4) - dispersion of undershoot                       1
        %        p(5) - ratio of response to undershoot                6
        %        p(6) - onset (seconds)                                0
        %        p(7) - length of kernel (seconds)                    32
        %
        sn  = 1;
        glm = 1;
        P0  = [6  16  1   1    6]';               % Default parameters for the SPM hrf
        LB  = [0 0 0 0 0]';%[0  9   0.2 0.2  3  -2 0.2]';     
        UB  = [7  16  10   10    10]'; 
        duration = 1;
        onsetshift = 0;
        fit = [1,2,5]'; % hrf parameter to be fitted
        roi = 1;
        eCriteria = 0.975;
        numIter = 10;
        vararginoptions(varargin,{'sn','glm','roi','fit','LB','UB','P0','duration','eCriteria'});
        
        hemi = 1;
        pre     = 4;
        post    = 20;
        
        T = []; Ts = [];
        for s=sn
            warning off
            % display to user which subject we are fitting
            subj = subj_name{s};
            fprintf('%s\n',subj);
            % load appropriate subject data
            cd(fullfile(glmDir{glm},subj_name{s}));
            load('SPM.mat');
            load(fullfile(regDir,sprintf('regions_%s',subj_name{sn})));
            % default onset and duration
            for r = 1:length(SPM.nscan)
                for u=1:length(SPM.Sess(r).U)
                    SPM.Sess(r).U(u).dur = ones(size(SPM.Sess(r).U(u).dur))*duration; % 1
                    SPM.Sess(r).U(u).ons = SPM.Sess(r).U(u).ons+onsetshift; % return to TR at announce trial
                end
                SPM.Sess(r).U=spm_get_ons(SPM,r);
            end
            % loop through rois for subject
            for r = roi
                fprintf('%s\n',R{r}.name);
                % Get data
                Y = region_getts(SPM,R{r});
                % Check Error before
                Ypre        = spm_filter(SPM.xX.K,SPM.xX.W*Y);
                Yres        = spm_sp('r',SPM.xX.xKXs,Ypre);
                Epre        = sum(sum(Yres.^2))/numel(Yres(:));
                % Fit a common hrf
                e = 1;
                P0_ = P0;
                iter = 1;
                fprintf('Iteration...');
                while ((e >= eCriteria)&&(iter<numIter))
                    fprintf('%d.',iter);
                    % fit hrf
                    [P,SPM,Yhat,Yres] = spmj_fit_hrf(SPM,Y,...
                        'fit',fit,'LB',LB,'UB',UB,'P0',P0_);
                    % update initial value
                    P0_(fit) = P0(fit)+0.1*rand*(UB(fit)-LB(fit));
                    iter = iter+1;
                    % Check Error after
                    Epost       = sum(sum(Yres.^2))/numel(Yres(:));
                    e           = Epost/Epre;
                end
                fprintf('Epost/Epre: %d\n',e);
                % Parameter values
                T_.fit = fit';
                T_.P = P(fit)';
                T_.R2 = 1-(trace(Yres*Yres')/(trace(Yres*Yres')+trace(Yhat*Yhat')));
                T_.Eratio = Epost/Epre;
                T_.regSide = hemi;
                T_.regType = r;
                T_.SN = s;
                T_.hemname = {hemName{hemi}};
                T_.regname = {regname{r}};
                T = addstruct(T,T_);
                % Timeseries
                D = spmj_get_ons_struct(SPM);
                y_hat = mean(Yhat,2);
                y_res = mean(Yres,2);
                for i=1:size(D.block,1)
                    D.y_hat(i,:)=cut(y_hat,pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_res(i,:)=cut(y_res,pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_adj(i,:)=D.y_hat(i,:)+D.y_res(i,:);
                end
                D.regType   = ones(size(D.event,1),1)*r;
                D.regSide   = ones(size(D.event,1),1)*hemi;
                D.sn        = ones(size(D.event,1),1)*s;
                Ts          = addstruct(Ts,D);
            end
            warning on
        end
        
        figure('name',[what,'_timeseries'])
        title(hemName{hemi});
        for reg = roi
            subset = (Ts.regType==reg & Ts.regSide==hemi);
            traceplot([-pre:post],Ts.y_adj,'errorfcn','stderr',...
                'split',[Ts.regType],'leg',regname(reg),'subset',subset); % ,
            hold on;
            traceplot([-pre:post],Ts.y_hat,'linecolor',[1 0 0],...
                'split',[Ts.regType], 'linewidth',3,'subset',subset); % ,
            %hold off;
        end
        xlabel('Seconds');
        ylabel('activation');
        xlim([-pre post]);
        drawline(0);

        keyboard
        % save data
        %         save(fullfile(regDir,sprintf('ROI_timeseries_fit_%d_%d.mat',glm_num{glm},fithrf)),'T','Ts');
        %
        %         varargout = {T,Ts};
    case 'ROI_plotTimeseriesAvg'                                            % (optional) :  Plots timeseries for specified region (avg. across digits, separated by pressing speed)
        glm = 1;
        sn  = 1;
        roi = 5;
        vararginoptions(varargin,{'sn','glm','roi'});
        
        figure('Color',[1 1 1]);
        D=load(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)));
        i = 1;
        for s=sn
            T=getrow(D,D.sn==s & D.roi==roi); % only plot chord regressors
            
                % plot timeseries
                plt.trace([-4:20],T.y_adj,'errorfcn','stderr');
                hold on;
                traceplot([-4:20],T.y_hat,'linewidth',2,'linecolor',[1 0 0]);
                hold off;
                xlabel('TR');
                ylabel('adjusted activation');
                xlim([-4 11]);
                title(sprintf('subj: %d  .', s));
                legend off
                i = i+1;
                
            
                drawline(0,'dir','vert');
                drawline(10,'dir','vert');
                drawline(15,'dir','vert');
                drawline(0,'dir','horz');

            
        end
        
        
        %__________________________________________________________________
    case 'ROI_plotTimeseries'                                               % (optional) :  Plots timeseries for specified region (avg. across digits, separated by pressing speed)
        glm = 2;
        sn  = 1;
        roi = 5;
        vararginoptions(varargin,{'sn','glm','roi'});
        
        numDigits = stimulationChords;
        numDigits = sum(numDigits,2)';
        
        shades = plt.helper.get_shades(4,'cool','increase',0);
        shades{5} = [1 0 0];
        shades{6} = [0 0.7 0];
        
        figure('Color',[1 1 1]);
        D=load(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)));
        i = 1;
        for s=sn
            T=getrow(D,D.sn==s & D.roi==roi); % only plot chord regressors
            for d = 1:5 % plot separately for numDigits
                subplot(length(sn),5,i)
                % plot timeseries
                sty = style.custom({'black'});
                sty.general.markersize = 1;
                plt.trace([-4:20],T.y_adj,'errorfcn','stderr','subset',ismember(T.event,find(numDigits==d)),'style',sty);
                hold on;
                traceplot([-4:20],T.y_hat,'linewidth',2,'subset',ismember(T.event,find(numDigits==d)),'linecolor',shades{roi});
                hold off;
                xlabel('TR');
                ylabel('adjusted activation');
                xlim([-4 10]);
                title(sprintf('subj: %d numD: %d  .', s,d));
                legend off
                i = i+1;
                
            end
%             if glm>2
%                 subplot(length(sn),5,i)
%                 plt.trace([-4:20],T.y_adj,'errorfcn','stderr','subset',T.event==32);
%                 hold on;
%                 traceplot([-4:20],T.y_hat,'linestyle',':','linewidth',2,'subset',T.event==32);
%                 hold off;
%                 xlabel('TR');
%                 ylabel('adjusted activation');
%                 xlim([-4 10]);
%                 title('foot reg  .');
%                 legend off
%                 i = i+1;
%             end
            plt.match('y');
            for j = 1:i-1
                subplot(length(sn),5,j);
                drawline(0,'dir','vert','linestyle',':');
                drawline(2.7,'dir','vert','linestyle',':');
                drawline(0,'dir','horz');
            end
            
        end
        
        
        %__________________________________________________________________
    case 'ROI_plotTimeseriesCond'                                           % (optional) :  Plots timeseries for specified region (avg. across digits, separated by pressing speed)
        glm = 1;
        sn  = 1;
        roi = 5;
        conds = 1:5;
        vararginoptions(varargin,{'sn','glm','roi','conds'});
        
        figure('Color',[1 1 1]);
        D=load(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)));
        i = 1;
        for s=sn
            T=getrow(D,D.sn==s & D.roi==roi); % only plot chord regressors
            for d = conds % plot separately for each cond
                subplot(2,5,i)
                % plot timeseries
                plt.trace([-4:20],T.y_adj,'errorfcn','stderr','subset',T.event==d);
                hold on;
                traceplot([-4:20],T.y_hat,'linestyle',':','linewidth',2,'subset',T.event==d);
                hold off;
                xlabel('TR');
                ylabel('adjusted activation');
                xlim([-4 11]);
                title(sprintf('cond: %d  .', s,d));
                legend off
                i = i+1;
                
            end
            plt.match('y');
            for j = i-length(conds):i-1
                subplot(2,5,j);
                drawline(0,'dir','vert');
                drawline(10,'dir','vert');
                drawline(0,'dir','horz');
            end
        end
        
        %__________________________________________________________________
        
    case 'ROI_getBetas'                                                     % STEP 5.3   :  Harvest activity patterns from specified rois
        
        I = pp1_imana('LIST_subjs');
        I = getrow(I,I.fmri_sessions==2);
        
        glm = 3;
        roi = [1:10,12:21];
        append = 0; % just add betas to currently existing datastructure
        vararginoptions(varargin,{'sn','glm','roi','append'});
        
        excludeROIs = [];%[12,24];
        T=[];
        if append
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        end
        % harvest
        for s=unique(I.sn)' % for each subj
            fprintf('\nSubject: %d\n',s) % output to user
            
            % load files
            load(fullfile(glmDir{glm}, I.origSN{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            D = load(fullfile(glmDir{glm}, I.origSN{s},'SPM_info.mat'));
            load(fullfile(regDir,sprintf('regions_%s.mat',subj_name{s})));          % load subject's region parcellation & depth structure (R)
            
            % add percent signal change imgs for subject
            Q = {}; 
            if glm==2
                numImgs = 31; 
            elseif glm==3
                numImgs = 32; 
            end
            for q = 1:numImgs
                Q{q} = (fullfile(glmDir{glm}, I.origSN{s}, sprintf('psc_%02d.nii',q))); 
            end
            Q = spm_vol(char(Q));
            
            % TR img info
            V = SPM.xY.VY; 

            % remove run means from patterns
            C0         = indicatorMatrix('identity',D.run); 
            ofInterest = 1:size(C0,1); % indicies for regressors of interest
            
            for r = roi % for each region
                % get raw data/psc for voxels in region
                Y = region_getdata(V,R{r});  % Data Y is N x P (P is in order of transpose of R{r}.depth)
                P = region_getdata(Q,R{r});
                % estimate region betas
                if sum(ismember(r,excludeROIs))==0
                    [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','runwise');
                else
                    [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','overall');
                end
                % toss stuff into output structure
                S.sn                 = s;
                S.roi                = r;
                S.tt                 = {D.tt};
                S.run                = {D.run};
                S.numDigits          = {D.numDigits};
                % remove nuisance regressor betas
                betaUW               = bsxfun(@rdivide,beta,sqrt(resMS));
                betaUW               = betaUW(ofInterest,:);
                betaW                = betaW(ofInterest,:);
                raw_beta             = beta(ofInterest,:);
                % add data to output structure
                S.betaW_noRunMean    = {betaW-C0*pinv(C0)*betaW};
                S.betaUW_noRunMean   = {betaUW-C0*pinv(C0)*betaUW};
                S.betaW              = {betaW};        
                S.betaUW             = {betaUW};  
                S.raw_beta           = {raw_beta};
                S.psc                = {P};
                S.resMS              = {resMS};
                S.xyzcoord           = {R{r}.data'}; % excl already applied
                S.depth = {[NaN]};     try S.depth = {R{r}.depth(R{r}.excl==0,:)'}; catch end
                S.flatcoord = {[NaN]}; try S.flatcoord = {R{r}.flatcoord(R{r}.excl==0,1:2)'}; catch end
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
        end
        % save T
        save(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)),'-struct','T'); 
        fprintf('\n')
    case 'ROI_nnClassify'
        % do nearest neighbour classifications. Also estimate chance level
        numIters = 1000; % number of chance shuffles
        betas    = varargin{1};
        condVec  = varargin{2};
        runVec   = varargin{3};
        
        % reorient betas, condVec, & runVec to column vectors (function
        % requires these dimensions)
        runs    = unique(runVec);
        conds   = unique(condVec);
        numRuns = numel(runs);
        numConds= numel(conds); 
        betas   = betas';
        condVec = condVec';
        runVec  = runVec';
        % estimate chance level
        chanceAcc = nan(numIters,1);
        for ii = 1:numIters
            % shuffle condition labels within each run
            shuffledCondVec = sample_wor(conds',numConds,numRuns);
            shuffledCondVec = shuffledCondVec(:)';
            % do classification with shuffled labels
            chanceAcc(ii) = crossval_classify(@classify_NearestNeighbour,betas,shuffledCondVec,runVec);
        end
        
        % estimate actual classification accuracy (across conditions)
        estAcc    = crossval_classify(@classify_NearestNeighbour,betas,condVec,runVec);
        pval      = 1-(find(estAcc<=sort(chanceAcc),1)/numIters); % proportion of chance accuracies greater than estimated accuracy
        chanceAcc = mean(chanceAcc);
        if isempty(pval)
            pval = realmin;
        end
        
        varargout = {estAcc,chanceAcc,pval};
        
    case 'ROI_stats'                                                        % STEP 5.4   :  Calculate stats/distances on activity patterns
        glm = 3;
        I = pp1_imana('LIST_subjs');
        I = getrow(I,I.fmri_sessions==2);
        vararginoptions(varargin,{'glm'});
        % housekeeping
        numDigits = stimulationChords;
        numDigits = sum(numDigits,2)';
        if glm==2 % 32 chords
            numConds = 31;
        elseif glm==3 % 31 chords + thumb response regressor
            numConds = 32;
        end
        % output structures
        To = [];
        Td = [];
        % get data
        T   = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        roi = unique(T.roi)';
        % do stats
        for s = unique(I.sn)' % for each subject
            D  = load(fullfile(glmDir{glm}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
            C0 = indicatorMatrix('identity',D.run);
            fprintf('\nSubject: %d\n',s)
            num_run = length(unique(D.run));
            ofInterest = 1:(numConds*num_run); % indicies for regressors of interest
            
            for r = roi % for each region
                S = getrow(T,(T.sn==s & T.roi==r)); % subject's region data
                betaW       = S.betaW{1}; 
                betaW_nmean = betaW(ofInterest,:)-C0*pinv(C0)*betaW(ofInterest,:); % run mean subtraction  
                % % Toverall structure stats
                % crossval second moment matrix
                [G,Sig]      = pcm_estGCrossval(betaW_nmean(ofInterest,:),D.run,D.tt);
                So.sig       = rsa_vectorizeIPM(Sig);
                So.G         = rsa_vectorizeIPM(G);
                So.G_wmean   = rsa_vectorizeIPM(pcm_estGCrossval(betaW(ofInterest,:),D.run,D.tt));
                % calculate empirical G (not cv)
                tmp.betaW    = betaW_nmean;
                tmp.tt       = D.tt;
                tmp          = tapply(tmp,{'tt'},{'betaW','mean'});
                So.G_emp     = rsa_vectorizeIPM(cov(tmp.betaW')); 
                % squared dissimilarities
                So.ldc_wmean = rsa.distanceLDC(betaW,D.run,D.tt);        % rdm crossvalidated, on patterns without run mean patterns removed
                So.ldc       = rsa.distanceLDC(betaW_nmean,D.run,D.tt);  % rdm crossvalidated, patterns with run means removed
                % do condition classification for all chords (NOTE: this pools data across sessions)
                H = D;
                H.beta = betaW;
                H  = getrow(H,H.chord<32);
                [So.nnEstAcc,So.nnChanceAcc,So.nnPval] = pp1_imana('ROI_nnClassify',H.beta,H.chord,H.run);
                % do condition classification for single fingers 
                H  = getrow(H,H.chord<6);
                [So.nnEstAccSF,So.nnChanceAccSF,So.nnPvalSF] = pp1_imana('ROI_nnClassify',H.beta,H.chord,H.run);
                
                % PSC
                So.psc = nanmean(S.psc{1},2)';
                if glm==2 % only regressors for chords
                    So.psc_chord = [1:31]; 
                    So.psc_numD  = numDigits;
                elseif glm==3 % regressors for chords and thumb response
                    So.psc_chord = [1:32]; 
                    So.psc_numD  = [numDigits,99];
                end
                % Calculate avg. betas for each condition
                Q                       = [];
                Q.raw_beta              = S.raw_beta{1};
                Q.tt                    = D.tt;
                Q                       = tapply(Q,{'tt'},{'raw_beta','mean'});
                So.avg_betas            = mean(Q.raw_beta,2)';
                So.avg_tt               = Q.tt';
                % indexing fields
                So.sn       = s;
                So.roi      = r;
                So.numVox   = size(betaW,2);
                So.regSide  = regSide(r);
                So.regType  = regType(r);
                To          = addstruct(To,So);
                % calc avg. chord patterns for each number of digits 
                d           = [];
                d.betaW     = betaW_nmean;
                d.numDigits = D.numDigits;
                d.run       = D.run;
                d.roi       = ones(size(d.run)).*r;
                d.chord     = D.chord;
                d           = getrow(d,d.chord<32);
                d5          = getrow(d,d.numDigits==5);
                d           = getrow(d,d.numDigits<5 & d.numDigits>0);
                d           = tapply(d,{'numDigits','run','roi'},{'betaW','mean'});
                d           = addstruct(d,d5);
                d           = rmfield(d,{'chord'});
                % calc distance between avg patterns for 1 finger up to
                % 5 finger chords:
                td.ldc = rsa.distanceLDC(d.betaW,d.run,d.numDigits)';
                td.distPair  = [1:10]';
                td.digitDiff = [1;2;3;4;1;2;3;1;2;1];
                td.roi       = ones(10,1).*r;
                td.sn        = ones(10,1).*s;
                Td           = addstruct(Td,td);
                fprintf('%d.',r)
            end % each region
        end % each subject

        % % save
        save(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)),'-struct','To');
        save(fullfile(regDir,sprintf('glm%d_reg_TnumDigits.mat',glm)),'-struct','Td');
        fprintf('done.\n')
   
    case 'ROI_pattConsist'                                                
        % Crossvalidated Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % This stat is useful for determining which GLM model yeilds least
        % variable beta estimates. That GLM should be the one you use for
        % future analysis cases.

        glm = 3;
        sn  = 1:6;
        roi = 1:8;
        removeMean = 1;
        conds = 1:31;
        vararginoptions(varargin,{'sn','glm','roi','removeMean'});
        
        % % Calculate pattern consistency for each roi, each subj.
        % Do so separately per session per subject.
        R = []; % output structure
        for g = glm
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',g))); % loads in struct 'T'
            for r = roi
                for s = sn
                    % get subj data
                    S      = getrow(T,(T.sn==s & T.roi==r));
                    b      = [];
                    b.beta = S.raw_beta{1};
                    b.run  = S.run{1};
                    b.tt   = S.tt{1};
                    b.sess = zeros(size(b.run));
                    b.sess(ismember(b.run,run{1}{s})) = 1;
                    b.sess(ismember(b.run,run{2}{s})) = 2;
                    b      = getrow(b,ismember(b.tt,conds)); % restrict to specific conditions
                    %for ii = 1:2 % per session
                    bs = b;    
                    %bs = getrow(b,b.sess==ii);
                        % calculate the pattern consistency
                        rs.r2              = rsa_patternConsistency(bs.beta,bs.run,bs.tt,'removeMean',removeMean);
                        [rs.r2_cv,rs.r_cv] = rsa_patternConsistency_crossval(bs.beta,bs.run,bs.tt,'removeMean',removeMean);
                        rs.sn              = s;
                        rs.roi             = r;
                        rs.glm             = g;
                        rs.numConds        = numel(conds);
                        rs.passive         = 1;
                        %rs.sess            = ii;
                        rs.removeMean      = removeMean;
                        R = addstruct(R,rs);
                    %end
                end
            end
        end
        %pivottable(R.glm,R.sn,R.r2,'mean','numformat','%0.4f');
        %pivottable(R.glm,R.sn,R.r2_cv,'mean','numformat','%0.4f');
        varargout = {R};
        % output arranged such that each row is an roi, each col is subj
    
    case 'ROI_MDS_overall'                                                  % (optional) :  Plots the scaled representational structure. 
        % enter region, glm #, sn (if desired)
        cplot = 'one';
        glm   = 3;
        roi   = 5; % default primary sensory cortex   
        sn    = 1;
        clrCode = 0; % if >0, color all chords with digit X red, and all other chords black.
        vararginoptions(varargin,{'roi','glm','cplot','clrCode','sn'});
        % cplot = 'all' to plot all 4 MDS figures (i.e. no contrast and 3 contrasts)- default
        % cplot = 'one'  to plot only no contrast MDS figure        

        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        T = getrow(T,T.roi==roi & ismember(T.sn,sn));
        IPM = mean(T.G_wmean,1); 
        
        chords    = stimulationChords;
        numDigits = sum(chords,2);

        if clrCode>0
            split = zeros(size(numDigits));
            split(chords(:,clrCode)>0) = 2;
            split(split==0) = 1;
        else
            split = numDigits;
        end
        
        % use to create rotation matrix for left hemi M1
        % r = 2
        % IPM2 = mean(T.IPM(T.region==r,:)); 
        % Y{2} = rsa_classicalMDS(IPM2,'mode','IPM');
        % [D,Z,Transform]=procrustes(Y{1},Y{2},'Scaling',false);
        % Y{2}=Y{2}*Transform.T;
        
        switch cplot
            case 'all' % do and plot 
                if glm>2
                    % account for thumb response regressor
                    CnumDigits = indicatorMatrix('identity',[numDigits;0]);
                    CnumDigits = bsxfun(@minus,CnumDigits,mean(CnumDigits,2));
                    Call   = eye(32)-ones(32)/32;
                else
                    CnumDigits = indicatorMatrix('identity',[numDigits]);
                    CnumDigits = bsxfun(@minus,CnumDigits,mean(CnumDigits,2));
                    Cdigit     = indicatorMatrix('identity',[1:5]);
                    Cdigit     = bsxfun(@minus,Cdigit,mean(Cdigit,2));
                    Call       = eye(31)-ones(31)/31;
                end

                G = rsa_squareIPM(IPM);
                G = G([1:5],[1:5]);
                IPMdigit = rsa_vectorizeIPM(G);
                Y{1} = rsa_classicalMDS(IPMdigit,'mode','IPM','contrast',Cdigit);
                
                Y{2} = rsa_classicalMDS(IPM,'mode','IPM');
                Y{3} = rsa_classicalMDS(IPM,'mode','IPM','contrast',Call);
                Y{4} = rsa_classicalMDS(IPM,'mode','IPM','contrast',CnumDigits);
                
                figure('Color',[1 1 1]);
                subplot(1,4,1);
                pp1_imana('MISC_scatterplotMDS',Y{1}(1:end,1:3),'split',split,'label',[1:5]','fig',gca);
                title('contrast: digit');
                subplot(1,4,2);
                pp1_imana('MISC_scatterplotMDS',Y{2}(1:end,1:3),'split',split,'label',[1:31]','fig',gca);
                title('no contrast');
                subplot(1,4,3);
                pp1_imana('MISC_scatterplotMDS',Y{3}(1:end,1:3),'split',split,'label',[1:31]','fig',gca);
                title('contrast: diff b/t all conds');
                subplot(1,4,4);
                pp1_imana('MISC_scatterplotMDS',Y{4}(1:end,1:3),'split',split,'label',[1:31]','fig',gca);
                title('contrast: num digits');
                
            case 'one' % only do and plot no contrast MDS
                Y{1} = rsa_classicalMDS(IPM,'mode','IPM');
                pp1_imana('MISC_scatterplotMDS',Y{1}(1:end,1:3),'split',split,'label',[1:31]');
        end
%         keyboard
    
    case 'ROI_patternReliability'                                           % plot w/in subj, w/in speed rdm reliability (Across two partitions), compare with across-speed correlations. Insights into RSA stability    
        % Splits data for each session into two partitions (even and odd runs).
        % Calculates correlation coefficients between each condition pair 
        % between all partitions.
        % Default setup includes subtraction of each run's mean
        % activity pattern (across conditions).
        glm           = 2;
        roi           = 5; % default roi
        sn            = 1;
        mean_subtract = 1; % subtract run means
        betaType      = 'raw'; % use raw, not normalized betas
        conds         = 1:31; % restrict to specific conditions?
        split         = 'oddEven'; % split runs into odd and even partitions
        % Correlate patterns across even-odd run splits within subjects.
        vararginoptions(varargin,{'roi','glm','sn','mean_subtract','betaType','conds','split'});
        
        % Select function to harvest approrpaite betas
        switch betaType
            case 'raw'
                betaFcn = 't.raw_beta{1}(1:length(D.tt),:);';
            case 'uni'
                betaFcn = 't.betaUW{1}(1:length(D.tt),:);';%'bsxfun(@rdivide,t.beta{1}1:length(D.tt),ssqrt(t.resMS{1}));';
            case 'mlt'
                betaFcn = 't.betaW{1}(1:length(D.tt),:);';
        end
        % Load subject's betas in all rois
        T   = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));   
        % Prep some variables (output structure, indexing variables, counter)
        Q = [];
        % loop across and correlate
        for s = sn % for each subject
            % cluster run numbers into partition splits
            switch split
                case 'oddEven' % split odd & even runs
                    runs       = unique(T.run{1});
                    oddRuns    = logical(mod(runs,2));
                    partitions = [{runs(oddRuns)},{runs(~oddRuns)}];
                case 'sess'    % split across session
                    partitions = [run{1}(sn),run{2}(sn)];
                case 'sess_cv_odd' % compare odd splits across sessions
                    sess1      = [cell2mat(run{1}(sn))];
                    sess2      = [cell2mat(run{2}(sn))];
                    runs       = [sess1 sess2];
                    sess       = [ones(1,length(sess1)) ones(1,length(sess2)).*2];
                    oddRuns    = logical(mod(runs,2));
                    runs       = runs(oddRuns);
                    sess       = sess(oddRuns);
                    partitions = [{runs(sess==1)},{runs(sess==2)}];
                case 'sess_cv_even' % compare even splits across sessions
                    sess1      = [cell2mat(run{1}(sn))];
                    sess2      = [cell2mat(run{2}(sn))];
                    runs       = [sess1 sess2];
                    sess       = [ones(1,length(sess1)) ones(1,length(sess2)).*2];
                    oddRuns    = logical(mod(runs,2));
                    runs       = runs(~oddRuns);
                    sess       = sess(~oddRuns);
                    partitions = [{runs(sess==1)},{runs(sess==2)}];
                case 'sess_cv_evenOdd' % compare even splits sess 1 to odd splits sess 2
                    sess1      = [cell2mat(run{1}(sn))];
                    sess2      = [cell2mat(run{2}(sn))];
                    runs       = [sess1 sess2];
                    sess       = [ones(1,length(sess1)) ones(1,length(sess2)).*2];
                    oddRuns    = logical(mod(runs,2));
                    partitions = [{runs(sess==1 & ~oddRuns)},{runs(sess==2 & oddRuns)}];
                case 'sess_cv_oddEven' % compare odd splits sess 1 to Even splits sess 2
                    sess1      = [cell2mat(run{1}(sn))];
                    sess2      = [cell2mat(run{2}(sn))];
                    runs       = [sess1 sess2];
                    sess       = [ones(1,length(sess1)) ones(1,length(sess2)).*2];
                    oddRuns    = logical(mod(runs,2));
                    partitions = [{runs(sess==1 & oddRuns)},{runs(sess==2 & ~oddRuns)}];
                case 'sess1'   % odd-even split within session 1
                    runs       = cell2mat(run{1}(sn));
                    oddRuns    = logical(mod(runs,2));
                    partitions = [{runs(oddRuns)},{runs(~oddRuns)}];
                case 'sess2'   % odd-even split within session 2
                    runs       = cell2mat(run{2}(sn));
                    oddRuns    = logical(mod(runs,2));
                    partitions = [{runs(oddRuns)},{runs(~oddRuns)}];
            end
            D = load(fullfile(glmDir{glm}, subj_name{s}, 'SPM_info.mat')); % load subject's trial structure
            for r = roi % for each roi
                t     = getrow(T,T.sn==s & T.roi==r); 
                betas = eval(betaFcn); % get specified patterns
                % remove run means?
                if mean_subtract
                    C0  = indicatorMatrix('identity',D.run);
                    betas = betas - C0*pinv(C0)* betas; % run mean subtraction  
                end
                % create partition indexes
                prepBetas = {};
                for i = 1:2%numel(partitions)
                    % select conditions for specified force levels
                    partitionIdx = logical(ismember(D.run,partitions{i}))';
                    condIdx{i}   = D.tt(partitionIdx);
                    prepBetas{i} = betas(partitionIdx,:);
                end
                % correlate patterns across partitions
                for c1 = conds % for each condition
                    % condition mean activity pattern for this run partition
                    oddCon   = (condIdx{1})==c1;
                    oddBetas = mean(prepBetas{1}(oddCon,:),1); 
                    for c2 = conds % for each condition
                        % condition mean activity pattern for the other partition
                        evenCon   = condIdx{2}==c2;
                        evenBetas = mean(prepBetas{2}(evenCon,:),1); 
                        % correlate condition patterns across partitions
                        tmp = corrcoef(evenBetas,oddBetas);
                        % harvest fields and allocate to output structure
                        q.corr  = tmp(1,2);
                        q.sn    = s;
                        q.roi   = r;
                        q.type  = ~isequal(c1,c2) + 1; % within condition correlations are type 1, between-condition are type 2
                        q.betas = betaType;
                        Q = addstruct(Q,q);
                    end
                end
            end % each region
        end
        varargout = {Q};  
    case 'ROI_compareReliability'
        % compares within-session pattern reliabilities to across-session
        % reliabilities. 
        % Reliabilities are odd-even splits.
        % corr(A_s1,B_s1) = r_s1
        % corr(A_s2,B_s2) = r_s2
        % corr(A_s1,A_s2) = r_A
        % corr(B_s1,B_s2) = r_B
        % ssqrt(r_s1 * r_s2) vs. ssqrt(r_A * r_B)
        sn  = 1:6;
        glm = 3;
        roi = 5;
        vararginoptions(varargin,{'roi','glm','sn'});
        
        conds = 1:31;
        mean_sub = 1;
        betaType = 'raw';
        
        R = [];
        for s = sn
            r_s1 = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'glm',glm,'split','sess1','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            r_s2 = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'glm',glm,'split','sess2','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            r_A  = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'glm',glm,'split','sess_cv_odd','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            r_B  = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'glm',glm,'split','sess_cv_even','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            %r_C  = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'glm',glm,'split','sess_cv_oddEven','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            %r_D  = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'glm',glm,'split','sess_cv_evenOdd','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            % keep only within condition comparisons:
            r_s1 = getrow(r_s1,r_s1.type==1);
            r_s2 = getrow(r_s2,r_s2.type==1);
            r_A  = getrow(r_A,r_A.type==1);
            r_B  = getrow(r_B,r_B.type==1);
            % add to output structure
            within  = ssqrt(mean(r_s1.corr) * mean(r_s2.corr));
            between = ssqrt(mean(r_A.corr) * mean(r_B.corr));
            r.corr  = [within; between];
            r.type  = [1;2];
            r.sn    = [s;s];
            r.roi   = [roi;roi];
            r.glm   = [glm;glm];
            R = addstruct(R,r);
        end
        % plt.dot(R.sn,R.corr,'split',R.type);
        %save(fullfile(regDir,sprintf('patternReliability_roi%d_glm%d.mat',roi,glm)),'-struct','R');
        varargout = {R};
    
    case 'ROI_getChordPSC'
        % makes a nearly-universal plotting structure from ROI_stats output
        glm = 3;
        vararginoptions(varargin,{'glm'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        D = [];
        d = [];
        v = ones(size(T.psc,2),1);
        for i = 1:size(T.sn,1)
            d.psc       = T.psc(i,:)';
            d.chord     = T.psc_chord(i,:)';
            d.numDigits = T.psc_numD(i,:)';
            d.roi       = v.*T.roi(i);
            d.sn        = v.*T.sn(i);
            D = addstruct(D,d);
            d = [];
        end
        varargout = {D};
    case 'ROI_getSingleFingerRatio'
        glm = 3;
        vararginoptions(varargin,{'glm'});
        % calculate variance-to-covarance of single finger second moment
        % per subject per roi:
        if glm==2
            Gtmp = zeros(31);
        elseif glm==3
            Gtmp = zeros(32); % account for thumb response condition
        end
        Gtmp(1:5,1:5) = eye(5);
        varIdx   = logical(rsa_vectorizeIPM(Gtmp));
        Gtmp(1:5,1:5) = ones(5)-eye(5);
        covarIdx = logical(rsa_vectorizeIPM(Gtmp));
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        D.sn  = T.sn;
        D.roi = T.roi;
        var     = mean(T.G(:,varIdx),2);
        covar   = mean(T.G(:,covarIdx),2);
        D.ratio = var./covar;
        varargout = {D};
    
    case 'calcFstat'
        % calculates F-statistic per voxel to determine if voxel is
        % significantly modulated by finger(s)
        X        = varargin{1}; % C*RxN matrix of data.
        runVec   = varargin{2};
        condVec  = varargin{3};
        % housekeeping
        numVox   = size(X,2);
        numCond  = length(unique(condVec));
        numRun   = length(unique(runVec));
        % calc F-stat per voxel
        C0  = indicatorMatrix('identity',runVec);
        X0  = X - C0*pinv(C0)*X; % remove run means
        A   = zeros(numCond,numVox,numRun);
        ApredCV = A;
        for i = 1:numRun
            A(:,:,i) = X0(runVec==i,:);
        end
        runs = 1:numRun;
        for i = 1:numRun
            testRuns = runs~=i;
            ApredCV(:,:,i) = mean(A(:,:,testRuns),3);
        end

        % non-crossval f-stat
        Apred = repmat(mean(A,3),1,1,numRun); % predicted voxel tunings
        TSS   = sum(sum((A).^2,3),1);         % total SS (null model- intercept)
        RSS   = sum(sum((A-Apred).^2,3),1);   % finger model (5 params)
        RSScv = sum(sum((A-ApredCV).^2,3),1); % unrestricted SSR
        dfN   = numCond - 1;                  % numerator DF is num params (5 fingers) - 1 (intercept)    
        dfD   = numCond*numRun - numCond;
        Fstat = ((TSS-RSS)./dfN) ./ (RSS./dfD);
        % crossval f-stat
        FstatCV = ((TSS-RSScv)./dfN) ./ (RSS./dfD);
        Fcrit   = finv(0.95,dfN,dfD);
        varargout = {Fstat,FstatCV,Fcrit};
        %keyboard
    
    case 'ROI_getSingleFingerG'
        glm = 3;
        roi = 1:8;
        I = pp1_imana('LIST_subjs');
        I = getrow(I,I.fmri_sessions==2);
        sn = unique(I.sn);
        vararginoptions(varargin,{'glm','roi','sn'});
        % helper vector
        if glm==3
            numCond = 32;
        elseif glm==2
            numCond = 31;
        end
        take = zeros(numCond);
        take(1:5,1:5) = ones(5);
        take1 = logical(rsa_vectorizeIPM(take));
        take = zeros(numCond);
        take(1:31,1:31) = ones(31);
        take2 = logical(rsa_vectorizeIPM(take));
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % cut data to output structure
        D = []; 
        D.glm        = ones(size(T.sn)).*glm;
        D.sn         = T.sn;
        D.roi        = T.roi;
        D.Gcv        = T.G(:,take1);
        D.Gemp       = T.G_emp(:,take1);
        D.Gcv_multi  = T.G(:,take2);
        D.Gemp_multi = T.G_emp(:,take2);
        D.numVox     = T.numVox;
        varargout = {D};       
    
    case 'ROI_getSingleFingerTuningSess'
        % do sft analysis within sessions per subj, then avg, results
        % across sessions
        glm = 3;
        roi = 1:4;
        I = pp1_imana('LIST_subjs');
        I = getrow(I,I.fmri_sessions==2);
        sn = unique(I.sn);
        vararginoptions(varargin,{'glm','roi','sn'});
        conds = 1:5;
        numCond = numel(conds);
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        D = []; % output structure
        v = ones(2,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b           = [];
            b.beta      = T.raw_beta{ii};
            b.tt        = T.tt{ii};
            b.run       = T.run{ii};
            % loop through and do separately per session
            for ss = 1:2
                bs  = getrow(b,ismember(b.tt,conds) & ismember(b.run,run{ss}{T.sn(ii)}));
                bs  = tapply(bs,{'tt'},{'beta','mean'});
                Gtt = cov(bs.beta');
                % rescale betaNoMax by setting the lowest value to 0,
                % calculating distance from maxB to each value.
                sftBeta         = pp1_imana('estSingleFingerTuning',bs.beta);
                sftBeta         = mean(sftBeta);
                [sftEV,sftDist] = pp1_imana('SFT:calcExpectedValue',Gtt,size(bs.beta,2)); % expected value of the null
                % do prob test for each subject on their expected sft distribution
                p = sum(sftBeta<=sftDist)/length(sftDist);
                if isempty(p)
                    p = realmin;
                end
                fprintf('s%02d roi%02d : p=%1.5f\n',T.sn(ii),T.roi(ii),p);
                % add to output structure
                d.passive   = v;
                d.sn        = v.*T.sn(ii);
                d.sess      = v.*ss;
                d.roi       = v.*T.roi(ii);
                d.glm       = v.*glm;
                d.sft       = [sftBeta;sftEV];
                d.sftProb   = v.*p;
                d.sftSD     = v.*std(sftDist);
                d.isEV      = [0;1];
                D = addstruct(D,d);
            end
        end
        save(fullfile(regDir,sprintf('sft_glm%d_sess',glm)),'-struct','D');
        varargout = {D};    
    case 'ROI_getSingleFingerTuning'
        glm = 3;
        roi = 1:4;
        I = pp1_imana('LIST_subjs');
        I = getrow(I,I.fmri_sessions==2);
        sn = unique(I.sn);
        vararginoptions(varargin,{'glm','roi','sn'});
        conds = 1:5;
        numCond = numel(conds);
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        D = []; % output structure
        v = ones(2,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b           = [];
            b.beta      = T.raw_beta{ii};
            b.tt        = T.tt{ii};
            b.run       = T.run{ii};
            b           = getrow(b,ismember(b.tt,conds));
            [fstat,fstatCV,fcrit] = pp1_imana('calcFstat',b.beta,b.run,b.tt);
            b           = tapply(b,{'tt'},{'beta','mean'});
            Gtt         = cov(b.beta');
            % calc tuning
            sftBeta         = pp1_imana('estSingleFingerTuning',b.beta);
            sftSigF         = mean(sftBeta(fstat>fcrit));       % voxels significantly tuned to fingers
            sftSigFcv       = mean(sftBeta(fstatCV>fcrit));     % voxels significantly (crossvalidated) tuned to fingers 
            sftBeta         = mean(sftBeta);                    % all voxels
            [sftEV,sftDist] = pp1_imana('SFT:calcExpectedValue',Gtt,size(b.beta,2)); % expected value of the null
            % do prob test for each subject on their expected sft distribution
            p = sum(sftBeta<=sftDist)/length(sftDist);
            if isempty(p)
                p = realmin;
            end
            fprintf('s%02d roi%02d : p=%1.5f\n',T.sn(ii),T.roi(ii),p);
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.glm       = v.*glm;
            d.sft       = [sftBeta;sftEV];
            d.sftF      = [sftSigF;sftEV];
            d.sftFcv    = [sftSigFcv;sftEV];
            d.sftProb   = v.*p;
            d.sftSD     = v.*std(sftDist);
            d.isEV      = [0;1];
            D = addstruct(D,d);
        end
        save(fullfile(regDir,sprintf('sft_glm%d_fcrit',glm)),'-struct','D');
        varargout = {D};
    case 'estSingleFingerTuning'
        % calculate single finger tuning using normalzied distance approach
        X         = varargin{1}; % CxN matrix of data.
        numC      = size(X,1);
        maxX      = max(X,[],1);
        avgDistsN = (sum(maxX-X)./(numC-1)) ./ (maxX-min(X,[],1));
        varargout = {avgDistsN};
    case 'SFT:calcExpectedValue'
        % calculates expected value (avg. sft across voxels) for voxels generated under specfified G
        
        % simulation params
        numSim  = 100;
        modelG  = varargin{1};
        numVox  = varargin{2};
        % generate data
        D = pp1_imana('SFT:model_G',modelG,numVox,numSim); % noiseless patterns
        % calc expected tuning on simulated datasets
        sft    = zeros(1,numSim);
        sftF   = sft;
        sftFcv = sft;
        runVec = ones(size(modelG,1),1);
        condVec = [1:5]';
        for s = 1:numSim
            d         = getrow(D,D.sn==s);
            tmp       = pp1_imana('estSingleFingerTuning',d.data);
%             [fstat,fstatCV,fcrit] = pp1_imana('calcFstat',d.data,runVec,condVec);
            sft(s)    = mean(tmp); % mean tuning across voxels for this simulated dataset
%             sftF(s)   = mean(tmp(fstat>fcrit));
%             sftFcv(s) = mean(tmp(fstatCV>fcrit));
        end
        sftAll = sft;
%         sftF   = mean(sftF);
%         sftFcv = mean(sftFcv);
        sft    = mean(sft); % mean tuning across datasets for this subject's simulations
        varargout = {sft,sftAll};
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
            d.data = U;
            d.sn   = v.*s;
            D = addstruct(D,d);
        end
        varargout = {D};
 
    case 'ROI_rdmStability'
        glm = 3;
        roi = 1:8;
        I = pp1_imana('LIST_subjs');
        I = getrow(I,I.fmri_sessions==2);
        sn = unique(I.sn);
        vararginoptions(varargin,{'glm','roi','sn'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        % housekeeping
        T = getrow(T,ismember(T.sn,sn));
        D = []; % output structure
        take = logical(tril(ones(numel(sn)),-1));
        for r = roi
            t = getrow(T,T.roi==r);
            R = corr(t.ldc');
            d = [];
            d.numSN = numel(sn); 
            d.roi   = r;
            d.corrs = R(take)';
            d.corr  = mean(R(take));
            % calc confidence bounds
            rz        = fisherz(d.corrs)'; 
            d.is_mean = fisherinv(mean(rz));
            d.is_LB   = fisherinv(d.is_mean - 1.96*stderr(rz));
            d.is_UB   = fisherinv(d.is_mean + 1.96*stderr(rz));
            D = addstruct(D,d);
        end
        varargout = {D};
        
    case '0' % ------------ HARVEST: condense/make new datastructures. ----
        % The STATS cases usually harvest some data from ROI stats structures.
        % These cases usually correspond to FIG cases of the same name.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'GET_idxPerNumDigits'
        glm = 2;
        numDigits = 1:5;
        vararginoptions(varargin,{'glm','numDigits'});

        % condition numbers
        if glm==2
            numConds = 31;
        elseif glm==3
            numConds = 32;           
        end
        
        numD = stimulationChords;
        numD = sum(numD,2)';
        
        if length(numD)~=numConds
            % set non-digit contrasts to zero
            numD(end+1:end+(numConds-length(numD))) = 0;
        end
        
        C = indicatorMatrix('allpairs',1:numConds);
        P = [];
        for ii = 1:size(C,1)
            if(sum(any(C(ii,numD~=0 & ismember(numD,numDigits)),1))==2)
                P = [P;ii];
            end
        end
        varargout = {P};
    case 'GET_LDCperNumDigits'    
        sn  = 1;
        glm = 4;
        roi = 5;
        vararginoptions(varargin,{'sn','glm','roi'});
        
        % % Gather distances between chords split by number of digits
        % stimulated (1,2,3,or 4).
        % - find column indicies of distance in full rdm for distance pairs
        % of X number of digits
        % - pull those distances from the rdm
        % - loop through until complete
        
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        D = [];
        d = [];
        for i = 1:size(T.sn,1)
            for dd = 1:4 % number of digits in chord
                P = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',dd); % returns row vector of indicies
                v = ones(size(P,1),1);
                d.ldc       = T.ldc(i,P')';
                d.distpair  = [1:length(P)]';
                d.roi       = v.*T.roi(i);
                d.numDigits = v.*dd;
                d.sn        = v.*T.sn(i);
                D = addstruct(D,d);
                d = [];
            end
        end
    
        D = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        varargout = {D};
    
    case 'GET_idxAcrossNumDigits'
        glm = 2;
        vararginoptions(varargin,{'glm'});
        
        % pulls out rdm column indicies for distances between chords with
        % different numbers of fingers stimulated.
        
        % condition numbers
        if glm==2
            numConds = 31;
        elseif glm==3
            numConds = 32;
        elseif glm==4
            numConds = 32;            
        end
        
        numD = stimulationChords;
        numD = sum(numD,2)';
        
        if length(numD)~=numConds
            % set non-digit contrasts to zero
            numD(end+1:end+(numConds-length(numD))) = 0;
        end
        
        C = indicatorMatrix('allpairs',1:numConds);
        P = [];
        for ii = 1:size(C,1)
            for dd=1:5
                if(sum(any(C(ii,numD~=0 & numD==dd),1))==2) % b/t chords of same number of digits
                    P = [P;ii];
                elseif(sum(any(C(ii,numD==0),1))==1) && (sum(any(C(ii,numD==dd),1))==1) % b/t error/foot and any condtion
                    P = [P;ii];
                end
            end
        end
        
        % now exclude indicies for within numdigits and comparisons between
        % chord and error/foot regressor
        Q = 1:(numConds*((numConds-1)/2));
        Q(P)=0;
        Q = unique(Q(Q>0));
        varargout = {Q};
    case 'GET_LDCacrossNumDigits'    
        sn  = 1;
        glm = 4;
        roi = 5;
        vararginoptions(varargin,{'sn','glm','roi'});
        
        % % Gather distances between chords split by number of digits
        % stimulated (1,2,3,or 4).
        % - find column indicies of distance in full rdm for distance pairs
        % of X number of digits
        % - pull those distances from the rdm
        % - loop through until complete
        
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        D = [];
        d = [];
        for i = 1:size(T.sn,1);
            P = pp1_imana('GET_idxAcrossNumDigits','glm',glm)'; % returns column vector of indicies, so transpose
            v = ones(size(P,1),1);
            d.ldc       = T.ldc(i,P')';
            d.distpair  = [1:length(P)]';
            d.roi       = v.*T.roi(i);
            d.sn        = v.*T.sn(i);
            D = addstruct(D,d);
            d = [];
        end
    
        D = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        varargout = {D};
    
 
    case '0' % ------------ STATS: statistical analyses. ------------------
        % The STATS cases usually harvest some data from ROI stats structures.
        % These cases usually correspond to FIG cases of the same name.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    case 'STATS_activityLinearScaling'
        sn  = 1:6;
        glm = 3;
        roi = 1:4;
        vararginoptions(varargin,{'sn','roi','glm'});
        % get psc data
        T  = pp1_imana('ROI_getChordPSC','glm',glm);
        T  = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        T  = getrow(T,T.numDigits<6);
        T = tapply(T,{'sn','roi','numDigits'},{'psc','mean'});
        % do linear fits
        numSubj = length(sn);
        S = [];
        for r = roi
           d = getrow(T,T.roi==r);
           y = d.psc;
           x = d.numDigits;
           % do linear regression in matrix form
           c = indicatorMatrix('identity',d.sn);
           Y = c.*repmat(y,1,numSubj);
           % do linear finger scaling fits
           X = c.*repmat(x,1,numSubj); % block diagonalize predictors
           X = [X c.*ones(size(c,1),numSubj)];   % add intercepts (one per subject)
           B = inv(X'*X)*X'*Y;                   % do linear regression
           % calculate fits
           Yp  = X*B;   % predicted data
           res = Y-Yp;  % residuals 
           s = [];
           s.rss = diag(res'*res);   % residual sums of squares
           s.tss = diag(Y'*Y);   % total sums of squares
           s.r2  = ones(numSubj,1) - (s.rss./s.tss); % r2 coefficient
           s.b_finger = diag(B([1:numSubj],[1:numSubj])); % snag finger-scaling betas
           s.b_int    = diag(B([numSubj+1:end],[1:numSubj])); % snag intercept betas
           s.sn  = d.sn(1:numSubj,1);
           s.roi = ones(numSubj,1).*r;
           S = addstruct(S,s);
        end
        varargout = {S};
    case 'STATS_activityLogScaling'
        sn  = 1:6;
        glm = 3;
        roi = 1:4;
        vararginoptions(varargin,{'sn','roi','glm'});
        % get psc data
        T  = pp1_imana('ROI_getChordPSC','glm',glm);
        T  = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        T  = getrow(T,T.numDigits<6);
        T = tapply(T,{'sn','roi','numDigits'},{'psc','mean'});
        % do log-linear fits
        numSubj = length(sn);
        S = [];
        for r = roi
           d = getrow(T,T.roi==r);
           y = d.psc;
           x = log(d.numDigits);
           % do linear regression in matrix form
           c = indicatorMatrix('identity',d.sn);
           Y = c.*repmat(y,1,numSubj);
           % do linear finger scaling fits
           X = c.*repmat(x,1,numSubj); % block diagonalize predictors
           X = [X c.*ones(size(c,1),numSubj)];   % add intercepts (one per subject)
           B = inv(X'*X)*X'*Y;                   % do linear regression
           % calculate fits
           Yp  = X*B;   % predicted data
           res = Y-Yp;  % residuals 
           s = [];
           s.rss = diag(res'*res);   % residual sums of squares
           s.tss = diag(Y'*Y);   % total sums of squares
           s.r2  = ones(numSubj,1) - (s.rss./s.tss); % r2 coefficient
           s.b_finger = diag(B([1:numSubj],[1:numSubj])); % snag finger-scaling betas
           s.b_int    = diag(B([numSubj+1:end],[1:numSubj])); % snag intercept betas
           s.sn  = d.sn(1:numSubj,1);
           s.roi = ones(numSubj,1).*r;
           S = addstruct(S,s);
        end
        varargout = {S};
    case 'dev_STATS_activityLinearScaling'
        sn  = 1:6;
        glm = 3;
        roi = 1:4;
        vararginoptions(varargin,{'sn','roi','glm'});
        % get psc data
        T  = pp1_imana('ROI_getChordPSC','glm',glm);
        T  = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        T  = getrow(T,T.numDigits<6);
        T = tapply(T,{'sn','roi','numDigits'},{'psc','mean'});
        % do linear fits
        numSubj = length(sn);
        S = [];
        for r = roi
           d = getrow(T,T.roi==r);
           y = d.psc;
           x = d.numDigits;
           % do linear regression in matrix form
           c = indicatorMatrix('identity',d.sn);
           Y = c.*repmat(y,1,numSubj);
           % first fit mean activity as predictor
           meanPSC = kron(mean(Y),ones(1,5))';
           X   = c.*repmat(meanPSC,1,numSubj);  % block diagonalize predictors
           B   = inv(X'*X)*X'*Y;                % do linear regression
           Yp  = X*B;                           % predicted data
           res = Y-Yp;                          % residuals 
           int_rss = diag(res'*res);
           % now do linear finger scaling fits
           X = c.*repmat(x,1,numSubj); % block diagonalize predictors
           %X = [X c.*ones(size(c,1),numSubj)];   % add intercepts (one per subject)
           B = inv(X'*X)*X'*Y;                   % do linear regression
           % calculate fits
           Yp  = X*B;   % predicted data
           res = Y-Yp;  % residuals 
           s = [];
           s.rss = diag(res'*res);   % residual sums of squares
           s.tss = diag(Y'*Y);   % total sums of squares
           s.int_rss = int_rss;
           s.r2  = ones(numSubj,1) - (s.rss./int_rss); % r2 coefficient
           %s.r2  = ones(numSubj,1) - (s.rss./s.tss); % r2 coefficient
           s.b_finger = diag(B([1:numSubj],[1:numSubj])); % snag finger-scaling betas
           %s.b_int    = diag(B([numSubj+1:end],[1:numSubj])); % snag intercept betas
           s.sn  = d.sn(1:numSubj,1);
           s.roi = ones(numSubj,1).*r;
           S = addstruct(S,s);
        end
        varargout = {S};
    case 'dev_STATS_activityLogScaling'
        sn  = 1:6;
        glm = 3;
        roi = 1:4;
        vararginoptions(varargin,{'sn','roi','glm'});
        % get psc data
        T  = pp1_imana('ROI_getChordPSC','glm',glm);
        T  = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        T  = getrow(T,T.numDigits<6);
        T = tapply(T,{'sn','roi','numDigits'},{'psc','mean'});
        % do log-linear fits
        numSubj = length(sn);
        S = [];
        for r = roi
           d = getrow(T,T.roi==r);
           y = d.psc;
           x = log(d.numDigits);
           % do linear regression in matrix form
           c = indicatorMatrix('identity',d.sn);
           Y = c.*repmat(y,1,numSubj);
           % first fit mean activity as predictor
           meanPSC = kron(mean(Y),ones(1,5))';
           X   = c.*repmat(meanPSC,1,numSubj);  % block diagonalize predictors
           B   = inv(X'*X)*X'*Y;                % do linear regression
           Yp  = X*B;                           % predicted data
           res = Y-Yp;                          % residuals 
           int_rss = diag(res'*res);
           % now do linear finger scaling fits
           X = c.*repmat(x,1,numSubj); % block diagonalize predictors
           %X = [X c.*ones(size(c,1),numSubj)];   % add intercepts (one per subject)
           B = inv(X'*X)*X'*Y;                   % do linear regression
           % calculate fits
           Yp  = X*B;   % predicted data
           res = Y-Yp;  % residuals 
           s = [];
           s.rss = diag(res'*res);   % residual sums of squares
           s.tss = diag(Y'*Y);   % total sums of squares
           s.int_rss = int_rss;
           s.r2  = ones(numSubj,1) - (s.rss./int_rss); % r2 coefficient
           %s.r2  = ones(numSubj,1) - (s.rss./s.tss); % r2 coefficient
           s.b_finger = diag(B([1:numSubj],[1:numSubj])); % snag finger-scaling betas
           %s.b_int    = diag(B([numSubj+1:end],[1:numSubj])); % snag intercept betas
           s.sn  = d.sn(1:numSubj,1);
           s.roi = ones(numSubj,1).*r;
           S = addstruct(S,s);
        end
        varargout = {S};
        
    case '0' % ------------ FIG: figures. ---------------------------------
        % The FIG cases usually harvest some data from ROI stats structures.
        % Sometimes they also do stats (since the data is harvested
        % accordingly).
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'FIG_avgBetas'
        sn  = 1;
        glm = 4;
        roi = 5;
        fig = [];
        vararginoptions(varargin,{'sn','glm','roi','fig'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
       
        D = [];
        d = [];
        v = ones(size(T.avg_tt,2),1);
        for i = 1:size(T.sn,1);
            d.avg_beta = T.avg_betas(i,:)';
            d.tt       = T.avg_tt(i,:)';
            d.roi   = v.*T.roi(i);
            d.sn    = v.*T.sn(i);
            D = addstruct(D,d);
            d = [];
        end

        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn));
        % plot
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.tt,D.avg_beta,'split',D.roi);
        plt.labels('roi','avg raw beta')
        drawline(0,'dir','horz');
        
        varargout = {D};
    case 'FIG_RDMline_LDC'
        sn           = 1:6;
        glm          = 3;
        roi          = 1:4;
        numDigits    = 1;
        fig          = [];
        split        = 'roi';
        vararginoptions(varargin,{'sn','glm','roi','numDigits','fig','split'});
        
        % % LDC lineplot for specified roi(s).
        % - load distances form ROI_stats
        % - determine indicies of distance pairs from overall RDM
        % - harvest specified distances from indicies
        % - arrange into output structure & plot
        
        D = pp1_imana('GET_LDCperNumDigits','sn',sn,'glm',glm,'roi',roi);
        D = getrow(D,D.numDigits==numDigits);
        % plot
        sty = style.custom(plt.helper.get_shades(length(eval(split)),'gray','increase',0));
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.distpair,ssqrt(D.ldc),'split',eval(sprintf('D.%s',split)),'style',sty);
        plt.labels('distance pair','ssqrt(ldc^2)')
        drawline(0,'dir','horz');
        
        varargout = {D};

    case 'PLOT_pattReliability'
        % plots pattern reliability w/in session and across session for
        % specified conditions (default = all)
        sn    = 1:5;
        glm   = 3;
        roi   = 5;
        vararginoptions(varargin,{'sn','roi','glm','conds'});
        
        % R = pp1_imana('ROI_compareReliability','sn',sn,'glm',glm,'roi',roi);
        R = load(fullfile(regDir,sprintf('patternReliability_roi%d_glm%d.mat',roi,glm)));
        R = getrow(R,ismember(R.sn,sn));
        plt.box(R.type,R.corr,'split',R.type);
        plt.labels('','Pearson''s r','raw pattern reliability');
        plt.set('xticklabel',{'within session','between session'},'xticklabelrotation',45);
        legend off
        ylims = ylim;
        ylim([0 ylims(2)]);
        
        varargout = {R};
    case 'PLOT_PSCperNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn  = 1:6;
        glm = 3;
        roi = [1:4]; % [3a, 3b, 1, 2]
        vararginoptions(varargin,{'sn','roi','glm'});
        
        D  = pp1_imana('ROI_getChordPSC','glm',glm);
        D  = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        D  = getrow(D,D.numDigits<6);
        Dr = tapply(D,{'sn','roi','numDigits'},{'psc','mean'});
        % plot
        style.use('numDigits');
        plt.box([Dr.roi Dr.numDigits],Dr.psc,'split',Dr.numDigits);
        xtick = {};
        for r = roi
            xtick = {xtick{:} '','',regname{r},'',''};
        end
        plt.set('xticklabel',xtick,'xticklabelrotation',45);
        plt.legend('east',{'1 digit','2 digits','3 digits','4 digits','5 digits'});
        plt.labels('region','percent signal change');
        drawline(0,'dir','horz');
        
        varargout = {Dr,D};
    case 'PLOT_singleFingerRatio'
        sn  = 1:6;
        glm = 3;
        roi = [2:4];
        vararginoptions(varargin,{'sn','roi','glm'});
        % get ratios
        D = pp1_imana('ROI_getSingleFingerRatio','glm',glm);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn));
        % plot
        sty = style.custom(plt.helper.get_shades(length(sn),'gray','descend'));
        plt.line(D.roi,D.ratio,'style',sty);
        xlabels = {};
        for i = 1:length(roi)
            xlabels{end+1} = reg_title{roi(i)};
        end
        plt.set(gca,'xticklabel',xlabels,'xticklabelrotation',45);
        plt.labels('','variance/covariance ratio','single-finger specificity');
        drawline(0,'dir','horz');
        varargout = {D};
    case 'PLOT_activityScaling'    
        sn  = 1:6;
        roi = 1:4;
        glm = 3;
        fig = [];
        vararginoptions(varargin,{'sn','roi','glm','fig'});
        % get fits
        S1 = pp1_imana('STATS_activityLinearScaling','sn',sn,'roi',roi,'glm',glm);
        S1.type = ones(size(S1.sn));
        S2 = pp1_imana('STATS_activityLogScaling','sn',sn,'roi',roi,'glm',glm);
        S2.type = ones(size(S2.sn)).*2;
        S = addstruct(S1,S2);
        % plot
        if isempty(fig)
            figure('Color',[1 1 1]);
        else
            fig;
        end
        plt.bar(S.roi,S.r2,'split',S.type);
        plt.labels('','r2 coefficient','psc scaling');
        xtick = {};
        for r = roi
            xtick = {xtick{:} '',reg_title{r}};
        end
        plt.set('xticklabel',xtick,'xticklabelrotation',45);
        varargout = {S};
        
    case 'PLOT_RDMsquare_avg'
        sn  = 1:6;
        glm = 3;
        roi = 1;
        chords = 1:5;
        vararginoptions(varargin,{'sn','roi','glm','chords'});
        % load data
        T = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat']));
        T = getrow(T,ismember(T.sn,sn) & T.roi==roi);
        % get rdm
        rdmVec = ssqrt(T.ldc);
        rdm = [];
        for i = 1:size(rdmVec,1)
            tmp = rdmVec(i,:);
            tmp = rsa_squareRDM(tmp);
            tmp = tmp(chords,chords);
            rdm(i,:) = rsa_vectorizeRDM(tmp);
        end
        % quickly calculate inter-subject rdm correlations
        r = corr(rdm');
        interSubjCorr = 1-squareform(1-r)';
        % plot
        rdm = rsa_squareRDM(mean(rdm));
        patchimg(rdm);
        title(reg_title{roi});
        set(gca,'xtick',[0.5:length(chords)],'xticklabel',chords);
        set(gca,'ytick',[0.5:length(chords)],'yticklabel',chords);
        ylim([0 length(chords)]);
        xlim([0 length(chords)]);
        
        varargout = {rdm,interSubjCorr};
        
    case 'FIG_LDCperNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn  = 1:6;
        glm = 3;
        roi = 1:4; % [3a, 3b, 1, 2, M1]
        fig = [];
        vararginoptions(varargin,{'sn','roi','glm','fig'});
        
        D = pp1_imana('GET_LDCperNumDigits','sn',sn,'glm',glm,'roi',roi);
        D.ldc = ssqrt(D.ldc);
        D = tapply(D,{'sn','roi','numDigits'},{'ldc','mean'});
        % plot
        style.use('numDigits');
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.roi,D.ldc,'split',D.numDigits);
        plt.labels('roi','ldc^2 between chords with same # digits')
        drawline(0,'dir','horz');
        varargout = {D};
    case 'FIG_LDCacrossNumDigits'
        % plot distance between chords with differing number of digits.
        sn    = 1;
        glm   = 4;
        roi = [1:4,6]; % [3a, 3b, 1, 2, M1]
        fig = [];
        vararginoptions(varargin,{'sn','roi','glm','fig'});
        
        D = pp1_imana('GET_LDCacrossNumDigits','sn',sn,'glm',glm,'roi',roi);
        D = tapply(D,{'sn','roi'},{'ldc','mean'});
        % plot
        style.use('default');
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.roi,D.ldc);
        plt.labels('roi','ldc^2')
        drawline(0,'dir','horz');
        
        varargout = {D};
    case 'FIG_LDCacrossAvgNumDigits'
        % plot distance between chords with differing number of digits.
        % distances are calculated between avg. chord patterns (done in
        % ROI_stats);
        sn    = 1;
        glm   = 4;
        roi = [1:4,6]; % [3a, 3b, 1, 2, M1]
        fig = [];
        vararginoptions(varargin,{'sn','roi','glm','fig'});
        
        D = load(fullfile(regDir,sprintf('glm%d_reg_TnumDigits.mat',glm)));
        D = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        % plot
        style.use('default');
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.roi,D.ldc,'split',D.digitDiff);
        plt.labels('roi','ldc^2')
        drawline(0,'dir','horz');
        
        varargout = {D};
    case 'FIG_plotLDCS1'
        sn = 1;
        glm = 4;
        roi = [1,2,3,4];
        figure('Color',[1 1 1]);
        subplot(1,2,1); pp1_imana('FIG_LDCperNumDigits','sn',sn,'glm',glm,'roi',roi,'fig',gca);
        plt.set('xticklabel',{'ba 3a','ba 3b','ba 1','ba 2'});
        plt.labels('roi','ldc^2','chords with same # digits  .');
        subplot(1,2,2); pp1_imana('FIG_LDCacrossAvgNumDigits','sn',sn,'glm',glm,'roi',roi,'fig',gca);
        plt.set('xticklabel',{'ba 3a','ba 3b','ba 1','ba 2'});
        plt.labels('roi','ldc^2','chords with diff # digits  .');
        plt.match('y');
        %subplot(1,3,3); pp1_imana('FIG_LDCperNumDigits','sn',sn,'glm',glm,'roi',roi,'fig',gca);

    case 'PLOT_sft'
        % plot single-finger preference analysis results
        roi = [1:4];
        glm = 3;
        vararginoptions(varargin,{'roi','glm'});
        % load data
        T = load(fullfile(regDir,sprintf('sft_glm%d_fcrit.mat',glm)));
        %T = tapply(T,{'passive','sn','roi','glm','isEV'},{'sft','mean'},{'sftProb','mean'},{'sftSD','mean'}); % avg. across sessions (if appropriate)
        T = getrow(T,ismember(T.roi,roi));
        % get roi labels
        labels1 = {};
        labels2 = {};
        for r=roi
            labels1{end+1} = '';
            labels1{end+1} = regname{r};
            labels2{end+1} = regname{r};
        end
        % plot styles
        sty1 = style.custom({'blue','darkgray'}); 
        sty2 = style.custom({'black'});           
        sty1.general.markersize = 4;
        sty2.general.markersize = 4;
        
        % plot sft and mean (of simulated subject samples) expected value
        subplot(1,3,1);
        plt.box(T.roi,T.sft,'split',T.isEV,'style',sty1);
        ylabel('single-finger tuning index');
        set(gca,'xticklabel',labels1,'xticklabelrotation',45);
        title('PASSIVE');
        % plot difference between calculated and simulated subject avg. sft
        subplot(1,3,2);
        plt.box(T.roi(T.isEV==0),T.sft(T.isEV==0)-T.sft(T.isEV==1),'style',sty2);
        drawline(0,'dir','horz','linestyle',':');
        set(gca,'xticklabel',labels2,'xticklabelrotation',45);
        ylabel('calculated - expected sft');
        title('PASSIVE (blue - gray bars)');
        % plot p-value of each subject's sft according to their own null
        % distribution (expected sft distribution)
%         subplot(1,3,3);
%         plt.dot(T.roi,T.sftProb,'style',sty2,'subset',T.isEV==0);
%         drawline(0,'dir','horz','linestyle',':');
%         set(gca,'xticklabel',labels2,'xticklabelrotation',45);
%         ylabel('p(H0)');
%         ylim([0 1]);
%         drawline(0.05,'dir','horz','linestyle',':');
%         title('subject-specific null prob');

        varargout = {T};
    case 'PLOT_nnClassification'
        % plots nearest neighbour classification accuracies for regions
        % (accuracies for all chords and single fingers in separate
        % subplots)
        glm = 3;
        roi = [1:4];
        I   = pp1_imana('LIST_subjs');
        I   = getrow(I,I.fmri_sessions==2);
        sn  = unique(I.sn)';
        vararginoptions(varargin,{'glm','sn','roi'});
        
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        
        % arrange into plotting structure
        D = [];
        v = ones(4,1);
        for ii = 1:length(T.sn)
            d.roi    = v.*T.roi(ii);
            d.sn     = v.*T.sn(ii);
            d.numVox = v.*T.numVox(ii);
            d.acc    = [T.nnEstAcc(ii);T.nnChanceAcc(ii);T.nnEstAccSF(ii);T.nnChanceAccSF(ii)];
            d.pval   = [T.nnPval(ii);nan;T.nnPvalSF(ii);nan];
            d.isEV   = [0;1;0;1];
            d.sf     = [0;0;1;1]; % all chords or single finger classification?
            D = addstruct(D,d);
        end
        
        sty = style.custom({'blue','black'});
        % plot all condition accuracy
        subplot(1,2,1);
        plt.box(D.roi,D.acc.*100,'split',D.isEV,'style',sty,'subset',D.sf==0);
        title('all conditions');
        ylabel('classification accuracy (%)');
        xlabel('region');
        %drawline((1/31)*100,'dir','horz','linestyle',':');
        % plot single finger accuracy
        subplot(1,2,2);
        plt.box(D.roi,D.acc.*100,'split',D.isEV,'style',sty,'subset',D.sf==1);
        title('single fingers');
        ylabel('classification accuracy (%)');
        xlabel('region');
        %drawline(20,'dir','horz','linestyle',':');

        varargout = {D};

    case '0' % ------------ PCM: pcm analyses. ----------------------------
    case 'usageG'
        U = load('/Users/sarbuckle/DATA/passivePatterns1/fmri/PCM_models/usageG.mat');
        varargout = {U.G,U.G_cent};
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
        M.name       = 'fc noiseceiling';
        M            = pcm_prepFreeModel(M);
        varargout = {M};
    case 'pcm_freedirect'
        % Naive averaring model- noise ceiling method 1- totall free model
        M.type       = 'freedirect';  
        M.numGparams = 0;
        M.name       = 'fd noiseceiling';
        varargout = {M};
    case 'pcm_noScale'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns, but no increase in avg. activity w/ chord length
        chords = pp1_simulations('chords');
        chords = chords./(kron(sum(chords,2),ones(1,5)));
        M.type       = 'feature';
        M.numGparams = 1;
        M.name       = 'no scale I fingers';
        M.Ac(:,:,1)  = chords;
        varargout = {M};
    case 'pcm_linearScale'
        % Model Linear: chord patterns are linear summations of single
        % fingers (independent)
        M.type       = 'feature';
        M.numGparams = 1;
        M.name       = 'L scale I fingers';
        M.Ac(:,:,1)  = pp1_simulations('chords');
        varargout = {M};
    case 'pcm_linearScale_FixedG'
        % Model Linear: chord patterns are linear summations of single
        % fingers (independent)
        M.type       = 'component';
        M.numGparams = 0;
        M.theta0     = [];
        M.name       = 'L scale single finger G';
%         M.Ac(:,:,1)  = pp1_simulations('chords');
        varargout = {M};
    case 'pcm_linearScale_nonlinearFingers'
        % Model nonLinear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates both finger params. Scaling params are fixed according
        % to the number of fingers in the chord.
        M.type       = 'nonlinear';
        M.name       = 'L scale NL fingers';
        M.modelpred  = @pp1_modelpred_linearScale_nonlinearFinger;
        M.numGparams = 14;
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_nonlinearScale'
        % Model nonLinear: chord patterns are linear summations of single
        % fingers (independent) scaled by number of fingers pressed.
        % Estimates ONLY scaling params- no finger params
        M.type       = 'nonlinear';
        M.name       = 'NL scale I fingers';
        M.modelpred  = @pp1_modelpred_nonlinearScale_fixedFinger;
        M.numGparams = 4;
        M.Gc(:,:,1)  = eye(5);
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_nonlinearScale_linearUsage'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates both scaling params with usage finger model.
        M.type       = 'nonlinear';
        M.name       = 'NL scale Usage fingers';
        M.modelpred  = @pp1_modelpred_nonlinearScale_fixedFinger;
        M.numGparams = 4;
        %M.theta0     = [];
        M.Gc(:,:,1)  = pp1_imana('usageG');
        varargout = {M};
    case 'pcm_nonlinearScale_linearH1'
        % nonlinear scaling with G as previous roi's single finger G
        M.type       = 'nonlinear';
        M.name       = 'NL scale H1 fingers';
        M.modelpred  = @pp1_modelpred_nonlinearScale_linearFinger;
        M.numGparams = 4;
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_nonlinearScale_linearH2'
        % nonlinear scaling with G as next roi's single finger G
        M.type       = 'nonlinear';
        M.name       = 'NL scale H2 fingers';
        M.modelpred  = @pp1_modelpred_nonlinearScale_linearFinger;
        M.numGparams = 4;
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_nonlinearScale_nonlinearFingers'
        % Model nonLinear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates both finger params and scaling params.
        M.type       = 'nonlinear';
        M.name       = 'NL scale NL fingers';
        M.modelpred  = @pp1_modelpred_nonlinearScale_nonlinearFinger;
        M.numGparams = 18;
        %M.theta0     = [];
        varargout = {M};
        
    case 'PCM_getData'
        % Get betas for roi from subjects in PCM-friendly format.
        % Betas do not have run means removed.
        sn  = 1;
        glm = 2;
        roi = 5; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        % load betas
        betaType = 'betaW'; % or betaUW, raw_beta
        B = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        B = getrow(B,B.roi==roi);
        % outputs
        Y = {};
        partVec = {};
        condVec = {};
        for i = 1:length(sn)
            % get subject data
            s = sn(i);
            b = getrow(B,B.sn==s);
            bb = [];
            bb.run   = cell2mat(b.run);
            bb.chord = cell2mat(b.tt);
            eval(sprintf('bb.betas = cell2mat(b.%s);',betaType));
            bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive patterns only
            % put subj data into pcm variables
            Y{i}         = bb.betas;
            partVec{i}   = bb.run;
            condVec{i}   = bb.chord;
            G_hat(:,:,i) = pcm_estGCrossval(Y{i},partVec{i},condVec{i});
        end
        varargout = {Y,partVec,condVec,G_hat};
    case 'PCM_defineModels'
        % case to define models fitted with PCM
        G_hat = varargin{1}; % 31x31 group avg. second moment matrix
        
        % define model params and starting values
        Asf         = pcm_diagonalize(G_hat(1:5,1:5));
        chords      = pp1_simulations('chords');
        Fx0         = pcm_free_startingval(G_hat(1:5,1:5));
        scaleParams = log([0.8 0.6 0.4 0.2])';
        
        % 1. null model
        M{1} = pp1_imana('pcm_null');
        % - - - - -
        % 2. no overall activity scaling, but patterns are linear finger
        % summations
%         M{end+1} = pp1_imana('pcm_noScale'); % this is dumb model, but eh.
        % - - - - -
        % 3. linear scaling, identity single-finger G
%         M{end+1} = pp1_imana('pcm_linearScale');
        % - - - - -
        % 4. linear scaling of avg. single-finger G
        M{end+1}  = pp1_imana('pcm_linearScale_FixedG');
        M{end}.Gc = chords*(Asf*Asf')*chords';
        % - - - - -
        % 5. linear scaling of flexible single-finger G
%         M{end+1}      = pp1_imana('pcm_linearScale_nonlinearFingers');
%         M{end}.theta0 = Fx0;
        % - - - - -
        % 6. non-linear scaling of identity single-finger G
%         M{end+1}      = pp1_imana('pcm_nonlinearScale');
%         M{end}.theta0 = scaleParams;
        % - - - - -
        % 8. non-linear scaling of avg. single-finger G
        M{end+1}      = pp1_imana('pcm_nonlinearScale');
        M{end}.theta0 = scaleParams;
        M{end}.Gc     = G_hat(1:5,1:5);
        M{end}.name   = 'NL scale single finger G';
        % - - - - -
        % 7. non-linear scaling of usage single-finger G
        M{end+1}      = pp1_imana('pcm_nonlinearScale_linearUsage');
        M{end}.theta0 = scaleParams;
        % - - - - -
        % 9. non-linear scaling of flexible single-finger G
%         M{end+1}      = pp1_imana('pcm_nonlinearScale_nonlinearFingers');
%         M{end}.theta0 = [Fx0;scaleParams];
        % - - - - -
        % X. non-linear scaling of single-finger G + additive background
        % pattern
        M{end+1}         = pp1_imana('pcm_nonlinearScale');
        M{end}.modelpred = @pp1_modelpred_nonlinearScale_fixedFinger_backgroundPattern;
        M{end}.theta0    = [scaleParams;flipud(scaleParams)];
        M{end}.numGparams= 8;
        M{end}.Gc        = G_hat(1:5,1:5);
        M{end}.name      = 'NL scale sf G add patt';

        % - - - - -
        % 10 & 11. noiseceilings
        M{end+1} = pp1_imana('pcm_freedirect'); 
%         M{end+1} = pp1_imana('pcm_freechol'); 
        
        varargout = {M};
    case 'PCM_fitModels_oneROI'
        % fits pcm models to data from one region
        sn     = 3;
        glm    = 2;
        roi    = 2;
        plotit = 0; % plot fits
        saveit = 0; % save fits
        runEffect = 'random';
        vararginoptions(varargin,{'sn','glm','roi','plotit','saveit'});
        
        % get data
        [Y,partVec,condVec,G_hat] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm);
        G_hat = mean(G_hat,3);
        
        % get models
        M = pp1_imana('PCM_defineModels',G_hat);

        % fit models
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',1);
        % plot fits?
        if plotit
            figure('Color',[1 1 1]);
            % plot fits
            pp1_simulations('plot_pcmFits',M,Tcv,T,G_pred);
            % plot avg. activity per number of digits
            subplot(2,numel(M),numel(M)+2);
            pp1_simulations('plot_chordActivity',Y,partVec);
            % plot G_hat
            subplot(2,numel(M),numel(M)+3);
            imagesc(G_hat);
            title('G hat (subj avg.)');
        end
        % save fits?
        if saveit
           outfile = fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,roi));
           save(outfile,'M','T','theta_hat','G_pred','Tcv','theta_cv'); 
        end
        
        %keyboard
        varargout = {Tcv,T,M,theta_cv,G_pred,Y,partVec};
    case 'PCM_fit'
        % Does pcm fitting for multiple rois
        sn  = 1:6;
        roi = [1:8];
        glm = 3;
        vararginoptions(varargin,{'sn','roi','glm'});
        for r = roi
            fprintf('\nroi: %d...',r);
            pp1_imana('PCM_fitModels_oneROI','sn',sn,'roi',r,'glm',glm,'plotit',0,'saveit',1);
        end
    case 'PCM_getFits'
        % Gets models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        glm = 3;
        roi = [];
        nNull = 1; % which model is null?
        nCeil = 4; % which model is noise ceiling ?
        vararginoptions(varargin,{'glm','roi','nNull','nCeil'});
        D   = []; % output structure
        for r = roi
            Tcv = [];
            % if exists, load pcm fits for region (otherwise, skip region)
            fprintf('\nroi %d...',r);
            try
                load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,r)));
            catch
                fprintf('no file.');
                continue
            end
            % scale likelihoods
            Tcv.likelihood_norm = bsxfun(@minus,Tcv.likelihood,Tcv.likelihood(:,nNull));
            % arrange into plotting structure
            numSubjs   = size(Tcv.SN,1);
            numModels  = numel(M);
            nameModels = {};
            Q = [];
            for m = 1:numModels
                % get model names
                nameModels{end+1,1} = M{m}.name;
                % get thetas
                if isempty(strfind(M{m}.name,'noiseceiling'))
                    q.thetaCV = num2cell(theta_cv{m}',2);
                else
                    q.thetaCV = num2cell(nan(numSubjs,1),2);
                end
                q.model = ones(numSubjs,1).*m;
                q.sn    = [1:numSubjs]';
                Q = addstruct(Q,q);
            end
            v = ones(numModels,1);
            for j = 1:numSubjs
                d.sn  = v.*Tcv.SN(j);
                d.roi = v.*r;
                d.model = [1:numModels]';
                d.modelName  = nameModels;
                d.likeNormCV = Tcv.likelihood_norm(j,:)';
                d.likeCV     = Tcv.likelihood(j,:)';
                d.likeNorm   = v.*nan;
                d.likeNorm(nCeil) = T.likelihood(j,nCeil) - Tcv.likelihood(j,nNull); % upper noise ceiling
                d.logLikeDiff= Tcv.likelihood(j,:)'- Tcv.likelihood(j,nCeil);
                d.thetaCV    = Q.thetaCV(Q.sn==j);
                D = addstruct(D,d);
            end
            fprintf('done.');
        end
        fprintf('\n');
        varargout = {D};
    case 'PCM_plotPseudoR2'
        % plots pseudo R2 value (0 = null, 1 = upper noise ceiling)
        glm = 3;
        roi = [1:4];
        sn  = 1:6;
        nNull = 1;
        nCeil = 4;
        modelsToPlot = [nNull,2,3,nCeil];
        vararginoptions(varargin,{'glm','roi','sn'});
        % get pcm fits
        D = pp1_imana('PCM_getFits','glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn) & ismember(D.model,modelsToPlot));
        % calc pseudo R2
        D.pseudoR2 = D.likeNormCV./kron(D.likeNorm(D.model==nCeil),ones(numel(unique(D.model)),1));
        % plot pcm fits
        numModels  = length(unique(D.model));
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        figure('Color',[1 1 1]);
        for i = 1:numPlots
            r = roi(i);
            subplot(1,numPlots,i);
            sty = style.custom(plt.helper.get_shades(numModels-2,'gray','descend'));
            plt.box(D.model,D.pseudoR2,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty);
%             plt.box(D.model,D.likeNormCV,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
            plt.labels('','pseudo model R2',regname{r});
%             plt.labels('','normalized log-likelihood',regname{r});
            % plot noise ceilings
            drawline(mean(D.pseudoR2(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-.'); % lower noise ceiling
            drawline(1,'dir','horz','linestyle','-');
%             drawline(mean(D.likeNormCV(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-.'); % lower noise ceiling
%             drawline(mean(D.likeNorm(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-'); % lower noise ceiling
            legend off
            ylims = ylim;
            ylim([0 ylims(2)]);
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotFitsBox'  
        % loads fit results per roi and plots them.
        glm = 3;
        roi = [1:4];
        sn  = 1:6;
        nNull = 1;
        nCeil = 10;
        modelsToPlot = [nNull,4,8,9,nCeil];
        vararginoptions(varargin,{'glm','roi','sn'});
        % get pcm fits
        D = pp1_imana('PCM_getFits','glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn) & ismember(D.model,modelsToPlot));
        % plot pcm fits
        numModels  = length(unique(D.model));
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        figure('Color',[1 1 1]);
        for i = 1:numPlots
            r = roi(i);
            subplot(1,numPlots,i);
            sty = style.custom(plt.helper.get_shades(numModels-2,'gray','descend'));
            plt.box(D.model,D.likeNormCV,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
            plt.labels('','relative log likelihood',regname{r});
            % plot noise ceilings
            drawline(mean(D.likeNormCV(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-.');
            drawline(mean(D.likeNorm(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-');
            legend off
            ylims = ylim;
            %ylim([0 ylims(2)]);
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotFitsBar'  
        % loads fit results per roi and plots them.
        glm = 3;
        roi = [1:4];
        sn  = 1:6;
        nNull = 1;
        nCeil = 10;
        modelsToPlot = [nNull,4,8,9,nCeil]; % first model is plotted as the null
        vararginoptions(varargin,{'glm','roi','sn'});
        numPlots = numel(roi);
        figure('Color',[1 1 1]);
        for i = 1:numPlots
            r = roi(i);
            load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,r)));
            T   = rmfield(T,{'reg'}); 
            T   = getrow(T,ismember(T.SN,sn));
            Tcv = rmfield(Tcv,{'reg'});
            Tcv = getrow(Tcv,ismember(Tcv.SN,sn));
            
            % % for plotting specified models
            Tp.SN           = T.SN;
            Tp.noise        = T.noise(:,modelsToPlot);
            Tp.scale        = T.scale(:,modelsToPlot);
            Tp.run          = T.run(:,modelsToPlot);
            Tp.likelihood   = T.likelihood(:,modelsToPlot);
                
            Tpcv.SN         = Tcv.SN;
            Tpcv.noise      = Tcv.noise(:,modelsToPlot);
            Tpcv.scale      = Tcv.scale(:,modelsToPlot);
            Tpcv.run        = Tcv.run(:,modelsToPlot);
            Tpcv.likelihood = Tcv.likelihood(:,modelsToPlot);
            
            Mp = {};
            for m = modelsToPlot
                Mp{end+1} = M{m};
            end
            % %
            Tpcv.likelihood_norm = bsxfun(@minus,Tpcv.likelihood,Tpcv.likelihood(:,1));
            % plot fits (errorbars are stderr)
            subplot(1,numPlots,i);
            pcm_plotModelLikelihood(Tpcv,Mp,'upperceil',Tp.likelihood(:,end),'style','bar','Nceil',length(modelsToPlot));
            set(gca,'xticklabelrotation',45);
            ylabel('relative log-likelihood')
            title(regname{r});
            box off;
        end
        plt.match('y');
        varargout = {Tpcv};
    case 'PCM_plotGpred'  
        % loads fit results per roi and plots them.
        glm = 3;
        roi = 2;
        nNull = 1;
        nCeil = 10;
        modelsToPlot = [nNull,4,8,9,nCeil];
        vararginoptions(varargin,{'glm','roi','sn'});
        if length(roi)>1
           error('can only call case with one roi'); 
        end
        % load fits
        load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,roi)));
        figure('Color',[1 1 1]);
        numPlots = numel(modelsToPlot);
        j = 1;
        for i = modelsToPlot
            % plot fits
            subplot(2,ceil(numPlots/2),j);
            imagesc(G_pred{i});
            title(sprintf('%s : %d params',M{i}.name,M{i}.numGparams));
            axis square
            j = j+1;
        end 
        varargout = {G_pred};    
    case 'PCM_plotThetas'
        % loads fit results per roi and plots them.
        glm   = 3;
        roi   = [1:4];
        sn    = 1:5;
        model = 8;
        thetas= 1:4;
        vararginoptions(varargin,{'glm','roi','sn'});
        % load plotting-friendly data structure
        D = pp1_imana('PCM_getFits','glm',glm,'roi',roi);
        D = getrow(D,D.model==model & ismember(D.roi,roi) & ismember(D.sn,sn));
        D.thetaCV = cell2mat(D.thetaCV); 
        % plot
        sty = style.custom(plt.helper.get_shades(length(roi),'parula','descend'));
        plt.trace([thetas],exp(D.thetaCV(:,thetas)),'split',D.roi,'style',sty);
        plt.labels('parameter no.','theta_cv value','single finger theta parameters');
    case 'PCM_singleFingerRatio'
        % loads fit results per roi and plots them.
        glm   = 3;
        roi   = 2:4;
        sn    = 1:5;
        model = 7; % most flexible model
        vararginoptions(varargin,{'glm','roi','sn','model'});
        varIdx = logical(rsa_vectorizeIPM(eye(5)));
        % load plotting-friendly data structure
        D = pp1_imana('PCM_getFits','glm',glm);
        D = getrow(D,D.model==model & ismember(D.roi,roi) & ismember(D.sn,sn));
        D.thetaCV = cell2mat(D.thetaCV); 
        D.thetaCV = [ones(size(D.sn)) D.thetaCV(:,1:14)]; % thumb variance is fixed to 1 in model fitting
        % calculate variance-to-covariance ratios    
        var   = mean(D.thetaCV(:,varIdx),2);
        covar = mean(D.thetaCV(:,~varIdx),2);
        D.ratio = var./covar;
        % plot
        %sty = style.custom({'black'});
        sty = style.custom(plt.helper.get_shades(length(sn),'gray','descend'));
        %plt.line(D.roi,D.ratio,'style',sty,'split',D.sn);
        plt.box(D.roi,D.ratio,'style',sty);
        xlabels = {};
        for i = 1:length(roi)
            xlabels{end+1} = reg_title{roi(i)};
        end
        plt.set(gca,'xticklabel',xlabels,'xticklabelrotation',45);
        plt.labels('','variance/covariance theta ratio','single-finger specificity');

    case 'PCM_addPatterns'
        sn  = 1:4;
        roi = 5;
        glm = 3;
        vararginoptions(varargin,{'glm','roi','sn'});
        if length(sn)>1
            error('can only call for one subect');
        elseif length(roi)>1
            error('can only call for one roi');
        elseif length(glm)>1
            error('can only call for one glm');
        end
        % get data
        [Y,partVec,condVec] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm);
        t.y      = Y{1};
        t.part   = partVec{1};
        t.chords = condVec{1};
        t = getrow(t,ismember(t.chords,1:31)); % restrict to chord regressors
        t = tapply(t,{'chords'},{'y','mean'}); % avg. patterns across runs
        t.sn        = ones(31,1).*sn;
        t.glm       = ones(31,1).*glm;
        t.roi       = ones(31,1).*roi;
        chords      = pp1_simulations('chords');
        t.numDigits = sum(chords,2);
        % make patterns
        y = t.y(ismember(t.chords,1:5),:);
        t.yA = (y' * chords')'; % linearly scale
        varargout = {t};
    case 'PCM_plotPatterns'
        sn  = 1:4;
        roi = 1:4;
        glm = 3;
        vararginoptions(varargin,{'glm','roi','sn'});
        T = [];
        for r = roi
            for s = sn
                t    = pp1_imana('PCM_addPatterns','sn',s,'roi',r,'glm',glm);
                % avg. true and additive patterns across voxels
                t.y  = mean(t.y,2);
                t.yA = mean(t.yA,2);
                T    = addstruct(T,t);
            end
        end
        % avg. each chord across subjects
        Ta = tapply(T,{'chords','numDigits','glm','roi'},{'y','mean'},{'yA','mean'});
        % plot
        style.use('numDigits');
        plt.dot(Ta.roi,Ta.yA,'split',Ta.numDigits);
        xtick = {};
        for r = roi
            xtick = {xtick{:} '','',reg_title{r},'',''};
        end
        plt.set('xticklabel',xtick,'xticklabelrotation',45);
        plt.legend('east',{'1 digit','2 digits','3 digits','4 digits','5 digits'});
        plt.labels('region','avg. beta');
        drawline(0,'dir','horz');
        
        varargout = {T,Ta};
    
    case 'PCM_plotActPassG'
        % plots active and passive Gs, and their associated difference.
        roi = 2;
        vararginoptions(varargin,{'roi'});
        glm_passive = 3;
        % load data
        load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm_passive,roi)));
        G_act = cpd2_imana('PCM_plotGpred','roi',roi); close gcf
        G_diff = G_act{end} - G_pred{end};
        chords = pp1_simulations('chords');
        % plot
        figure('Color',[1 1 1]);
        subplot(1,4,1);
        imagesc(chords);
        title('configurations');
        axis equal
        xlim([0.5 5.5]);
        subplot(1,4,2);
        imagesc(G_pred{end});
        title('avg. passive Gcv');
        subplot(1,4,3);
        imagesc(G_act{end});
        title('avg. active Gcv');
        subplot(1,4,4);
        imagesc(G_diff);
        title('difference (a-p)');
        for ii = 2:4
            subplot(1,4,ii);
            axis square
            drawline(5.5,'dir','vert');
            drawline(5.5,'dir','horz');
            drawline(15.5,'dir','vert');
            drawline(15.5,'dir','horz');
            drawline(25.5,'dir','vert');
            drawline(25.5,'dir','horz');
            drawline(30.5,'dir','vert');
            drawline(30.5,'dir','horz');
        end
        varargout = {G_diff};
        
    case '0' % ------------ fingerpics: project patterns onto surface & take pictures.
        % You will absolutely need to edit these cases. 
        % These versions of the cases are frome ef1_imana
        % (extensionflexion).
        % 1 values of coordinates in xlims and ylims correspond to 1/10mm
        % distances.
        % So a range of 20 = 2mm distance on surface projection.
    case 'fingerpics'                                                       % Makes jpegs of finger activity patterns on cortical surface M1/S1
        sn  = 1;
        glm = 3;
        vararginoptions(varargin,{'sn','glm'});
        
        for s=sn;
            for g=glm;
                %pp1_imana('surf_mapFingers','sn',s,'glm',g)
                pp1_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[1,5]);
                pp1_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[2]);
                %pp1_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[3]);
                %pp1_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[4]);
            end
        end
    case 'surf_mapFingers'                                                % Map locations of finger patterns- run after glm estimation
        % map volume images to metric file and save them in individual surface folder
        sn  = 1;
        glm = 3;
        vararginoptions(varargin,{'sn','glm'});
        hemisphere = 1;   % left hemi

        for c = 1:31 % contrast #
            fileList{c} = sprintf('spmT_00%02d.nii',c); % see case 'contrast' for contrast number index
        end;
        for s = sn
            for h = hemisphere
                caretSDir = fullfile(caretDir,[atlasA,subj_name{s}],hemName{h});
                %specname  = fullfile(caretSDir,[atlasA,subj_name{s} '.' hem{h}   '.spec']);
                white     = fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial      = fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                topo      = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                
                C1 = caret_load(white);
                C2 = caret_load(pial);
                
                images = {};
                for f=1:length(fileList)
                    images{f} = fullfile(glmDir{glm},subj_name{s},fileList{f});
                end
                
                M  = caret_vol2surf_own(C1.data,C2.data,images,'topo',topo,'exclude_thres',0.75,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,sprintf('%s_glm%d_hemi%d_finger.metric',subj_name{s},glm,h)),M);
            end
        end;   
    case 'surf_fingerpatterns'             % Make finger pattern jpegs
        %close all;
        sn  = 1;
        glm = 3;
        numDigits = [1,5];
        vararginoptions(varargin,{'sn','glm','numDigits'});
        
        h = 1; % left hemi
        groupDir = [gpCaretDir filesep hemName{h} ];
        cd(groupDir);
        switch(h)
            case 1 % left hemi
                coord  = 'lh.FLAT.coord';
                topo   = 'lh.CUT.topo';
                %data  = 'lh.surface_shape';  
%                 xlims=[-4 15]; % may need to adjust locations for pics
%                 ylims=[-10 9];
                xlims=[-4 18]; % may need to adjust locations for pics
                ylims=[-9 20];
            case 2 % right hemi
                coord  = 'rh.FLAT.coord';
                topo   = 'rh.CUT.topo';
                %data  = 'rh.surface_shape';
                xlims  = [-10 20];
                ylims  = [-15 30];
        end;
        
        % load Central Sulcus border line (to plot as dashed line in pics)
        border = fullfile(caretDir,'fsaverage_sym',hemName{h},['BA_borders.border']);
        B      = caret_load(border);
        % set path to caret surface patterns
        data = fullfile(caretDir,['x' subj_name{sn}],hemName{h},sprintf('%s_glm%d_hemi%d_finger.metric',subj_name{sn},glm,h));

        % plot topographic image of surface reconstruction (w/out patterns)
        figure('Color',[1 1 1]); % make figure
        sshape = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.surface_shape']);
        M      = caret_plotflatmap('col',2,'data',sshape,'border',B.Border,...
                 'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims,'bordersize',10);
        colormap('bone');
        %close gcf;
        
        % plot patterns for single finger stimulation
        figure('Color',[1 1 1]); % make figure 
        digits    = sum(stimulationChords,2);
        digitCols = find(ismember(digits,numDigits));
        j = 1;
        for i = digitCols'
            subplot(1,length(digitCols),j);
            [M,d]   = caret_plotflatmap('M',M,'col',i,'data',data,...
                        'border',B.Border,'bordersize',10,'topo',topo,'coord',coord,'bordercolor',{'k.'});
            colormap('parula');
            j = j+1;
        end;
        
        mm = 3; % force colour scaling on patterns
        % loop through both figures and all subplots to: force colour
        % scaling and label conditions
        for i = 1:length(digitCols)
            subplot(1,length(digitCols),i);
            caxis([-mm/2 mm]);   % scale color across plots
            set(gca,'XTick',[]); % remove X and Y axis ticks
            set(gca,'YTick',[]);
            %axis equal;
            box on
            ax = get(gca);
            ax.XAxis.LineWidth = 4;
            ax.YAxis.LineWidth = 4;
        end % for each subplot
    
    case 'surf_map_patterns'                                                % Map locations of finger patterns- run after glm estimation
        % map volume images to metric file and save them in individual surface folder
        sn  = 1;
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        hemisphere = 1;   % left hemi

        for c = 1:5 % contrast #
            fileList{c} = sprintf('spmT_00%02d.nii',c); % see case 'contrast' for contrast number index
        end;
        for s = sn
            for h = hemisphere
                caretSDir = fullfile(caretDir,[atlasA,subj_name{s}],hemName{h});
                %specname  = fullfile(caretSDir,[atlasA,subj_name{s} '.' hem{h}   '.spec']);
                white     = fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial      = fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                topo      = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                
                C1 = caret_load(white);
                C2 = caret_load(pial);
                
                for f=1:length(fileList)
                    images{f} = fullfile(glmDir{glm},subj_name{s},fileList{f});
                end
                
                M  = caret_vol2surf_own(C1.data,C2.data,images,'topo',topo,'exclude_thres',0.75,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,sprintf('s%02d_glm%d_hemi%d_finger.metric',s,glm,h)),M);
            end
        end;
    case 'surf_condpatterns'                                                % Make finger pattern jpegs
        close all;
        sn  = 1;
        glm = 1;
        xlims = [-12 7]; % optimized for left-hemi (contralateral to movement)
        ylims = [-5 14];
        label = 1;
        forces = [1.5,2];
        vararginoptions(varargin,{'sn','glm','xlims','ylims','forces'});

        h = 1; % left hemi
        groupDir = [gpCaretDir filesep hemName{h} ];
        cd(groupDir);
        switch(h)
            case 1 % left hemi
                coord  = 'lh.FLAT.coord';
                topo   = 'lh.CUT.topo';
                %data  = 'lh.surface_shape';  
%                 xlims=[-12 7]; % may need to adjust locations for pics
%                  ylims=[-5 14];
                xlims = [-17 2];
                ylims = [-2 17];
            case 2 % right hemi
                coord  = 'rh.FLAT.coord';
                topo   = 'rh.CUT.topo';
                %data  = 'rh.surface_shape';
                xlims  = [-10 20];
                ylims  = [-15 30];
        end;
        
        % load Central Sulcus border line (to plot as dashed line in pics)
        border = fullfile(caretDir,'fsaverage_sym',hemName{h},['CS.border']);
        B      = caret_load(border);
        % set path to caret surface patterns
        data = fullfile(caretDir,['x' subj_name{sn}],hemName{h},sprintf('s%02d_glm%d_hemi%d_finger.metric',sn,glm,h));

        % plot topographic image of surface reconstruction (w/out patterns)
        figure('Color',[1 1 1]); % make figure
        sshape = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.surface_shape']);
         M      = caret_plotflatmap('col',2,'data',sshape,'border',B.Border,...
                 'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims,'bordersize',10);
        colormap('bone');
        %close gcf;
        
        % plot patterns for flexion
        flex_forces     = forces(ismember(forces,cond_forces{1}));     
        if ~isempty(flex_forces)
            flex_con_forces = [ones(5,1).*1.5; ones(5,1).*2; ones(5,1).*2.5];
            force_col_num   = 1:15;
            force_flex_num   = force_col_num(ismember(flex_con_forces,forces)); % column numbers for forces to plot
            f{1} = figure('Color',[1 1 1]); % make figure 
            for i = 1:(5*length(flex_forces))
                subplot(length(flex_forces),5,i);
                [M,d]   = caret_plotflatmap('M',M,'col',force_flex_num(i),'data',data,...
                            'border',B.Border,'bordersize',10,'topo',topo,'coord',coord);
                colormap('jet');
            end;
        end
        % plot patterns for extension
        ext_forces      = forces(ismember(forces,cond_forces{2}));
        if ~isempty(ext_forces)
            ext_con_forces  = [ones(5,1); ones(5,1).*1.5; ones(5,1).*2];
            force_col_num   = 16:30;
            force_ext_num   = force_col_num(ismember(ext_con_forces,forces)); % column numbers for forces to plot
            f{2} = figure('Color',[1 1 1]); % make figure w/ black background
            for i = 1:(5*length(ext_forces))
                subplot(length(ext_forces),5,i);
                [M,d]   = caret_plotflatmap('M',M,'col',force_ext_num(i),'data',data,...
                            'border',B.Border,'bordersize',10,'topo',topo,'coord',coord);
                colormap('jet');
            end;
        end
        
        % set coordinates for condition labels (as text) on figures
        xText = xlims(1):xlims(2);
        xText = xText(floor(length(xText)/2));
        yText = ylims(2)-2;
        forces = {flex_forces,ext_forces};
        dir_label = {'Fx','Ex'};
        dig_label = {repmat([1:5]',length(forces{1}),1),repmat([1:5]',length(forces{2}),1)};
        
        mm = 10; % force colour scaling on patterns
        % loop through both figures and all subplots to: force colour
        % scaling and label conditions
        for j = 1:2
            if ~isempty(forces{j})
                % bring figure to front
                set(0,'currentfigure',f{j});
                % get subplot titles
                force_titles = kron(forces{j}',ones(5,1));
                for i=1:length(force_titles)
                    subplot(length(forces{j}),5,i);
                    caxis([-mm/2 mm]);   % scale color across plots
                    set(gca,'XTick',[]); % remove X and Y axis ticks
                    set(gca,'YTick',[]);
                    %axis equal;
                    box on
                    ax = get(gca);
                    ax.XAxis.LineWidth = 4;
                    ax.YAxis.LineWidth = 4;
                    if label
                        text(xText,yText,sprintf('\b %s d%d %.1fN',dir_label{j},dig_label{j}(i),force_titles(i)),...
                            'interpreter','tex','Color',[0 0 0],'HorizontalAlignment','center','FontSize',11);
                    end
                end % for each subplot
            end
        end; % for each figure

        keyboard
        % save figures
        movement_dir = {'flex','ext'};
        for j = 1:2
            set(0,'currentfigure',f{j})
            saveas(gcf, [movement_dir{j},'_',subj_name{sn},'_',hemName{h},'_',sprintf('%d',mm)], 'jpg')
        end
    case 'surf_meanpatterns'
        close all;
        sn  = 1;
        glm = 1;
        xlims = [-12 7]; % optimized for left-hemi (contralateral to movement)
        ylims = [-5 14];
        label = 1;
        forces = [1.5,2];
        vararginoptions(varargin,{'sn','glm','xlims','ylims','forces'});
        
        h = 1; % left hemi
        groupDir = [gpCaretDir filesep hemName{h} ];
        cd(groupDir);
        switch(h)
            case 1 % left hemi
                coord  = 'lh.FLAT.coord';
                topo   = 'lh.CUT.topo';
                %data  = 'lh.surface_shape';  
                xlims=[-12 7]; % may need to adjust locations for pics
                ylims=[-5 14];
            case 2 % right hemi
                coord  = 'rh.FLAT.coord';
                topo   = 'rh.CUT.topo';
                %data  = 'rh.surface_shape';
                xlims  = [-10 20];
                ylims  = [-15 30];
        end;
        
        % load Central Sulcus border line (to plot as dashed line in pics)
        border = fullfile(caretDir,'fsaverage_sym',hemName{h},['CS.border']);
        B      = caret_load(border);
        % set path to caret surface patterns
        data = fullfile(caretDir,['x' subj_name{sn}],hemName{h},sprintf('s%02d_glm%d_hemi%d_finger.metric',sn,glm,h));

        % plot topographic image of surface reconstruction (w/out patterns)
        figure('Color',[1 1 1]); % make figure
        sshape = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.surface_shape']);
         M      = caret_plotflatmap('col',2,'data',sshape,'border',B.Border,...
                 'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims,'bordersize',15);
        colormap('bone');
        close gcf;
        
        % harvest contrast numbers to plot 
        % forces for each contrast are found in script variable 'cond_forces'
        contrast_nums   = {[31,32,33],[34,35,36]}; % avg flex patterns (at each force), avg ext patterns (at each force)
        subplot_nums    = {[2 3 4],[5 6 7]};
        subplot_nums    = {subplot_nums{1}(ismember(cond_forces{1},forces)),...
                           subplot_nums{2}(ismember(cond_forces{2},forces))};
        contrast_labels = {{'Avg Fx 1.5N','Avg Fx 2.0N','Avg Fx 2.5N'},...
                           {'Avg Ex 1.0N','Avg Ex 1.5N','Avg Ex 2.0N'}};
        contrast_labels = {{contrast_labels{1}{ismember(cond_forces{1},forces)}},...
                           {contrast_labels{2}{ismember(cond_forces{2},forces)}}};  
        flex_forces     = forces(ismember(forces,cond_forces{1})); 
        ext_forces      = forces(ismember(forces,cond_forces{2}));
        forces_to_plot  = {flex_forces,ext_forces};
        cols_to_plot    = {contrast_nums{1}(ismember(cond_forces{1},forces)),...
                           contrast_nums{2}(ismember(cond_forces{2},forces))};
        
        % plot avg movement patterns (for specified forces)
        mm       = 7; % force colour scaling on patterns
        num_rows = size(forces_to_plot{1},1) + size(forces_to_plot{2},1);
        f{1}     = figure('Color',[1 1 1]); % make figure 
        counter  = 1;
        xText = xlims(1):xlims(2);
        xText = xText(floor(length(xText)/2));
        yText = ylims(2)-2;
        for j = 1:2 % for flexion/extension
            for i = 1:length(forces_to_plot{j})
                subplot(num_rows,4,subplot_nums{j}(i));
                [M,d]   = caret_plotflatmap('M',M,'col',cols_to_plot{j}(i),'data',data,'cscale',[-mm/2 mm],...
                            'border',B.Border,'bordersize',10,'topo',topo,'coord',coord);
                colormap('jet');
                % styling
                set(gca,'XTick',[]); % remove X and Y axis ticks
                set(gca,'YTick',[]);
                %axis equal;
                box on
                ax = get(gca);
                ax.XAxis.LineWidth = 4;
                ax.YAxis.LineWidth = 4;
                % text
                if label
                    text(xText,yText,contrast_labels{j}{i},...
                        'interpreter','tex','Color',[0 0 0],'HorizontalAlignment','center','FontSize',11);
                end
                counter = counter + 1;
            end;
        end
        
        keyboard
        % save figure
        %saveas(gcf, [movement_dir{j},'_',subj_name{sn},'_',hemName{h},'_',sprintf('%d',mm)], 'jpg')
             
    case 'figure_voxdepth'                  % Plot voxel depth for each subject
        
        %roi = [4];
        roi=[1:4 8:11];    % S1-1,2,3   M1-4
        sn = [1:7];
        vararginoptions(varargin,{'roi','sn'});
        
        figure;
        f=1; % subplot counter
        for j = 1:length(sn)
            s=sn(j); % get subj number (useful if specifying non-incremental subjs)
            load(fullfile(regDir,sprintf('s0%d_regions.mat',s)));
            for r=roi
                d = R{r}.depth; % voxel depths
                % % plotting stuff
                subplot(length(sn),length(roi),f);
                title(sprintf('Subj %d roi %d',s,r))
                hold on
                xlim([-1.5,2.5]);
                ylim([-0.15,0.15]);
                set(gca,'YTick',[]);
                drawline(0,'dir','horz');
                plot(d,0,'kx','markersize',10);
                drawline([0,1],'dir','vert','color',[1 0 0],'lim',[-0.05,0.05]);
                drawline([0.5],'dir','vert','color',[0 0 1],'lim',[-0.05,0.05]);
                % add text for numvox at certain depth thresholds
                text([-0.5,0.5,1.5,0.25,0.75],[0.07,0.07,0.07,-0.07,-0.07],{sprintf('%d',sum(d<0)),... % above pial(0)
                    sprintf('%d',sum((d>=0)&(d<=1))),...  % between pial and white
                    sprintf('%d',sum(d>1)),...            % below white
                    sprintf('%d',sum((d>=0)&(d<=0.4999))),...  % superficial layer
                    sprintf('%d',sum((d>=0.5001)&(d<=1)))},... % deep layer
                    'interpreter','latex'); %latex makes nicer text
                % add text for surface types
                text([0,1],[-0.1,-0.1],{'pial','white'},'interpreter','latex');
                %--------------------------------%
                f=f+1; % increase subplot counter
                hold off
                %fprintf('vox between q2-q3: %d\n',sum((d>=0.25)&(d<=0.75)))
            end
        end
        
        %set(gcf,'PaperPosition',[2 2 3*length(sn) 3*length(roi)+0.5]);
        wysiwyg;
    case 'figure_cortical_depth'            % Plot histogram of cortical depth in roi for each subject
        % draws red line on last plot to indicate average depth acorss
        % subjects
        sn = [1:7];
        numbins = 20; % number of bins to concatinate data for histplot
        hemi = 2; % left hemi default (2 is right hemisphere)
        roi = 11; % set to 0 for no specified roi
        vararginoptions(varargin,{'sn','numbins','hemi','roi'})
        
        f1=figure;
        spindx = 1; % subplot number
        for r=roi % make enough figures for difference from mean roi plots
            f2{r+1} = figure; title(sprintf('ROI %d',r));
        end
        
        for s=sn;
            % load pial and white-gray surfaces
            a = caret_load(fullfile(caretDir,['x',subj_name{s}],hemName{hemi},'rh.PIAL.coord'));
            b = caret_load(fullfile(caretDir,['x',subj_name{s}],hemName{hemi},'rh.WHITE.coord'));
            % Both structures contain .data subfield, which has [x,y,z]
            % coords in mm of each vertex location.
            % Calc distance between surfaces w/ d = sqrt((x2 - x1)^2)
            d = sqrt(sum((a.data-b.data).^2,2));
            % load region data- want R.location arrays
            load(fullfile(regDir,[subj_name{s},'_regions.mat']));
            for r=roi;
                if r==0 % roi = 0 if we take all rois together
                    dindx = [1:length(d)];
                    dat = d(dindx);
                else
                    dindx = R{r}.location;
                    dat = d(dindx);
                end
                figure(f1);
                subplot(numel(sn),numel(roi),spindx);
                [N,~]= histplot(dat,'numcat',numbins);
                title(sprintf('Subj %d ROI %d Cortical Thickness (mm)',s,r));
                ax{s,r+1} = gca;
                d_max(1,s) = max(dat);
                d_min(1,s) = min(dat);
                y_lim(1,s) = max(N);
                spindx = spindx+1;
                
                g_dat{r+1}(:,s) = dat;
                
                all_depths(s,r) = mean(dat);
            end
        end
        subplot(numel(sn),numel(roi),spindx-1); hold on;
        drawline(mean(all_depths(:,roi)),'dir','vert','color',[1 0 0]);
        hold off;     
    
    case 'WRAPPER_searchlight'                                               
        glm = 3;
        vararginoptions(varargin,{'sn','glm'});
        % You can call this case to run searchlight analyses.
        % 'sn' can be an array of subjects.
        for s = sn 
            pp1_imana('SEARCH_define','sn',s,'glm',glm);
            pp1_imana('SEARCH_ldcRun','sn',s,'glm',glm);
            %pp1_imana('SEARCH_ldcContrasts','sn',s,'glm',glm);
        end
    case 'SEARCH_define'                                                    % STEP 4.1   :  Defines searchlights for 120 voxels in grey matter surface
        glm = 3;
        vararginoptions(varargin,{'sn','glm'});
        
        mask       = fullfile(glmDir{glm},subj_name{sn},'mask.nii');
        Vmask       = rsa.readMask(mask);

        LcaretDir = fullfile(caretDir,sprintf('xs0%d',sn),'LeftHem');
        RcaretDir = fullfile(caretDir,sprintf('xs0%d',sn),'RightHem');
        white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
        pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
        topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};
        S         = rsa_readSurf(white,pial,topo);

        L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[15 120]);
        save(fullfile(anatomicalDir,subj_name{sn},sprintf('s%d_searchlight_120.mat',sn)),'-struct','L');
    case 'SEARCH_ldcRun'                                                    % STEP 4.2   :  Runs LDC searchlight using defined searchlights (above)
        % Requires java functionality unless running on SArbuckle's
        % computer.
        glm = 3;
        vararginoptions(varargin,{'sn','glm'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        % go to subject's glm directory 
        spmDir = fullfile(glmDir{glm},subj_name{sn});
        cd(spmDir);
        % load their searchlight definitions and SPM file
        L = load(fullfile(anatomicalDir,subj_name{sn},sprintf('s%d_searchlight_120.mat',sn)));
        load SPM;
        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{sn}));
        % make index vectors
        D = load('SPM_info.mat');
        conditionVec  = D.tt;
        partition     = D.run;
        name = sprintf('s0%d_glm%d',sn,glm);
        % run the searchlight
        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,...
                                'analysisName',name,'idealBlock',block,...
                                'spmDir',spmDir);
        cd(cwd);
    case 'SEARCH_ldcContrasts'                                        
        % Calls 'MISC_SEARCH_calculate_contrast'
        sn  = 1;
        glm = 3;
        con = {'passive','singleFinger','multiFinger','thumbPress'};
        vararginoptions(varargin,{'sn','glm','con'});
        % Use 'con' option to define different contrasts.
        %   'avg'    :  Average LDC nii for all 20 conds,
        %                invariant of speed (avg of 10 pairwise distances)

        % Load subject surface searchlight results (1 vol per paired conds)
        LDC_file            = fullfile(glmDir{glm},subj_name{sn},sprintf('%s_glm%d_LDC.nii',subj_name{sn},glm)); % searchlight nifti
        [subjDir,fname,ext] = fileparts(LDC_file);
        cd(subjDir);
        vol  = spm_vol([fname ext]);
        vdat = spm_read_vols(vol); % is searchlight data
        % For each of the predefined contrast types (see above)...
        for c = 1:length(con)
            switch con{c}
                case 'passive' % just average across all paired distances
                    gidx    = 1:size(vdat,4);
                case '1digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',1); 
                case '2digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',2);
                case '3digitChords'
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',3);
                case '4digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',4);
                case 'numDigits'
                    gidx    = pp1_imana('GET_idxAcrossNumDigits','glm',glm);
            end
            % avg. distances according to contrast selected
            Y.LDC   = vdat(:,:,:,gidx);
            Y.LDC   = ssqrt(Y.LDC);
            Y.LDC   = nanmean(Y.LDC,4); 
            % prep output file
            Y.dim   = vol(1).dim;
            Y.dt    = vol(1).dt;
            Y.mat   = vol(1).mat;    
            % save output
            Y.fname   = sprintf('%s_glm%d_%sLDC.nii',subj_name{sn},glm,con{c});
            Y.descrip = sprintf('exp: ''pp1'' \nglm: ''FAST'' \ncontrast: ''%s''',con{c});

            spm_write_vol(Y,Y.LDC);
            fprintf('Done %s\n',Y.fname);

            clear Y
        end
    
    

    otherwise
        fprintf('%s: no such case.\n',what)

end

%% ------------------------- Local funcs ----------------------------------
function dircheck(dir)
% Checks existance of specified directory. Makes it if it does not exist.
% SA 01/2016
if ~exist(dir,'dir');
    %warning('%s didn''t exist, so this directory was created.\n',dir);
    mkdir(dir);
end
end

function Q = findFiles(folder,pattern)
   % find files in 'folder' with names that match 'pattern' (to a degree)
   L     = dir(folder);
   Q = {};
   for i = 1:numel(L); 
       try
           if strfind(L(i).name,pattern)
               Q{end+1,1} = L(i).name;
           end
       catch
       end
   end
   
end

function AllCombos = stimulationChords
AllCombos = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;...             % singles            5
             1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;...                        % doubles (thumb)    4
             0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;...                                   % doubles            3
             0 0 1 1 0; 0 0 1 0 1;...                                              % doubles            2
             0 0 0 1 1;...                                                         % doubles            1
             1 1 1 0 0; 1 1 0 1 0; 1 1 0 0 1; 1 0 1 1 0; 1 0 1 0 1; 1 0 0 1 1;...  % triples (thumb)    6
             0 1 1 1 0; 0 1 1 0 1; 0 1 0 1 1; 0 0 1 1 1;...                        % triples            4
             1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 0 1 1 1; 0 1 1 1 1;                % quadruples         5
             1 1 1 1 1];                                                           % all five           1
end
%% ------------------------- Old code -------------------------------------
% old GLM cases (may be useful in future)
% case 3
%     % model error trials as separate regressor
%     if c<numConds
%         idx	= find(R.chordNum==c & R.isError==0); % find indx of all trials in run of that condition 
%         if ~isempty(idx)
%             % Correct start time for numDummys removed & announce time, & convert to seconds
%             J.sess(r).cond(c).onset    = [R.startTimeMeas(idx)/1000 - J.timing.RT*numDummys + R.cueTime(idx)/1000];    
%             J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
%             J.sess(r).cond(c).tmod     = 0;
%             J.sess(r).cond(c).orth     = 0;
%             J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
%             S.incl = 1;
%         else
%             % All trials of this cond type were errors in
%             % this block. Zero-pad regressor column.
%             % Correct start time for numDummys removed & announce time, & convert to seconds
%             J.sess(r).cond(c).onset    = 0;    
%             J.sess(r).cond(c).duration = 0;                          
%             J.sess(r).cond(c).tmod     = 0;
%             J.sess(r).cond(c).orth     = 0;
%             J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
%             S.incl = 0;
%         end
%         S.isError = 0;
%     else
%         % these are error trials to exclude
%         idx	= find(R.isError==1);
%         % Correct start time for numDummys removed & announce time, & convert to seconds
%         J.sess(r).cond(c).onset    = [R.startTimeMeas(idx)/1000 - J.timing.RT*numDummys + R.cueTime(idx)/1000];    
%         J.sess(r).cond(c).duration = R.stimTime(idx)/1000                           % durations of task we are modeling (not length of entire trial)
%         J.sess(r).cond(c).tmod     = 0;
%         J.sess(r).cond(c).orth     = 0;
%         J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
%         J.sess(r).cond(c).name = 'error_trials'; 
%         S.isError = 1;
%         S.incl    = 1;
%     end
%     J.sess(r).cond(c).name = sprintf('chord_%d',R.chordNum(find(R.chordNum==c,1)));  % make condition name (for user readability)
% case 4
%     % model error trials as separate regressor
%     if c<numConds
%         idx	= find(R.chordNum==c); % find indx of all trials in run of that condition 
%         % Correct start time for numDummys removed & announce time, & convert to seconds
%         J.sess(r).cond(c).onset    = [R.startTimeMeas(idx)/1000 - J.timing.RT*numDummys + R.cueTime(idx)/1000];    
%         J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
%         J.sess(r).cond(c).tmod     = 0;
%         J.sess(r).cond(c).orth     = 0;
%         J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
% 
%         S.incl = 1;
%         S.isFoot  = 0;
%     else
%         % these are trials with foot movement
%         idx	= find(R.isError==1 & R.RT>0);
%         if ~isempty(idx)
%             % Correct start time for numDummys removed & announce time, & convert to seconds
%             J.sess(r).cond(c).onset    = [R.startTimeMeas(idx)/1000 - J.timing.RT*numDummys + R.cueTime(idx)/1000 + R.RT(idx)/1000];    
%             J.sess(r).cond(c).duration = 1;                           % duration of foot movement- let's say 1 second
%             J.sess(r).cond(c).tmod     = 0;
%             J.sess(r).cond(c).orth     = 0;
%             J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
%             J.sess(r).cond(c).name = 'footmovement_trials'; 
%             S.incl   = 1;
%         else
%             % No foot movement this block
%             J.sess(r).cond(c).onset    = 0;    
%             J.sess(r).cond(c).duration = 0;                          
%             J.sess(r).cond(c).tmod     = 0;
%             J.sess(r).cond(c).orth     = 0;
%             J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
%             S.incl = 0;
%         end
%         S.isFoot = 1;
%     end
%     J.sess(r).cond(c).name = sprintf('chord_%d',R.chordNum(find(R.chordNum==c,1)));  % make condition name (for user readability)
% case 5
%     error('not yet implemented')
%     % model error trials as separate regressor
%     if c<32
%         idx	= find(R.chordNum==c & R.isError==0); % find indx of all trials in run of that condition 
%         if ~isempty(idx)
%             % Correct start time for numDummys removed & announce time, & convert to seconds
%             J.sess(r).cond(c).onset    = [R.startTimeMeas(idx)/1000 - J.timing.RT*numDummys + R.cueTime(idx)/1000];    
%             J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
%             J.sess(r).cond(c).tmod     = 0;
%             J.sess(r).cond(c).orth     = 0;
%             J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
%             S.incl = 1;
%         else
%             % All trials of this cond type were errors in
%             % this block. Zero-pad regressor column.
%             % Correct start time for numDummys removed & announce time, & convert to seconds
%             J.sess(r).cond(c).onset    = 0;    
%             J.sess(r).cond(c).duration = 0;                          
%             J.sess(r).cond(c).tmod     = 0;
%             J.sess(r).cond(c).orth     = 0;
%             J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
%             S.incl = 0;
%         end
%         S.isError = 0;
%     else
%         % these are error trials to exclude
%         idx	= find(R.isError==1);
%         % Correct start time for numDummys removed & announce time, & convert to seconds
%         J.sess(r).cond(c).onset    = [R.startTimeMeas(idx)/1000 - J.timing.RT*numDummys + R.cueTime(idx)/1000];    
%         J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
%         J.sess(r).cond(c).tmod     = 0;
%         J.sess(r).cond(c).orth     = 0;
%         J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
%         J.sess(r).cond(c).name = 'error_trials'; 
%         S.isError = 1;
%         S.incl    = 1;
%     end
% J.sess(r).cond(c).name = sprintf('chord_%d',R.chordNum(find(R.chordNum==c,1)));  % make condition name (for user readability)
%                 


% take = logical(tril(ones(5),-1));
% c = rsa.util.pairMatrix(5);
% D = [];
% for i = 1:size(T.sn,1)
%     g = T.G_wmean(i,:);
%     g = rsa_squareIPM(g);
%     g = g(1:5,1:5);
%     d.dist = zeros(2,1);
%     d.dist(1) = mean(diag(g)) - mean(g(take));
%     d.dist(2) = mean(sum((c*g).*c,2));
%     d.type    = [1;2];
%     d.sn   = ones(2,1).*T.sn(i);
%     d.roi  = ones(2,1).*T.roi(i);
%     D = addstruct(D,d);
% end

cd(cwd); % return user to their directory when function was called
end