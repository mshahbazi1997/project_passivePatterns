function varargout = pp1_imana(what,varargin)
% function    varargout = fmri_imana(what,varargin)
%  
% Operates with vararginoptions, eg.:
%       fmri_imana(analysis_step,'sn',1,'roi',2)
%
% This is scaffolding for imaging analysis code.
%
% WRAPPER cases submit batch processing stages. WRAPPER cases can accept an
% array of subject numbers with 'sn' option (this is often not the case
% when calling individual preprocessing cases).
%
% IMPORTANT NOTES:
%   With every new subject, you must update variables in 'Subject Things'.
%
%      'dircheck' is a local function that creates a new directories so it
%      can save appropriate datafile structures. 
%
%   Of course, for full functionality you will need many functions from
%      various Diedrichsen Lab toolboxes (sourced through github and the
%      lab website).

% SArbuckle, Motor Control Group, 2019, UWO
% saarbuckle@gmail.com


% -------------------------------------------------------------------------
cwd = cd; % get current directory when called, and return at end of script
% ------------------------- Directories -----------------------------------
codeDir         ='/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_passivePatterns';
baseDir         ='/Users/sarbuckle/DATA/passivePatterns1/fmri';   % base directory for analysis
% directories for "raw" data (i.e. data during/from preprocessing)
behavDir        =[baseDir '/data'];          % behavioural data directory
dicomDir        =[baseDir '/data_dicom'];    % imgs hot off the scanner
imagingDirRaw   =[baseDir '/imaging_data_raw']; 
% working data directories (this data has been preprocessed and is ready to go)
fieldmapDir     =[baseDir '/fieldmaps/'];    % this explicit path is not required here, but is helpful to know it is required by realign_unwarp
imagingDir      =[baseDir '/imaging_data'];               
anatomicalDir   =[baseDir '/anatomicals'];  
cerebAnatDir    =[baseDir '/c_anatomicals'];  
scAnatDir       =[baseDir '/sc_anatomicals']; 
freesurferDir   =[baseDir '/surfaceFreesurfer'];          
caretDir        =[baseDir '/surfaceCaret'];     
gpCaretDir      =[caretDir '/fsaverage_sym'];
regDir          =[baseDir '/RegionOfInterest/'];   % where most of the important analysis data structures will be saved       
pcmDir          =[baseDir '/PCM_models'];
% update glmDir when adding new glms
glmDir          ={[baseDir '/glm1'],[baseDir '/glm2'],[baseDir '/glm3']};              
    
filePrefix = 'pp1_fmri';
% check that these directories exist (may be best to comment out once run
% the first time)
% dircheck(behavDir);
% dircheck(dicomDir);
% dircheck(imagingDirRaw);
% dircheck(phaseDirRaw);

% dircheck(fieldmapDir);
% dircheck(phaseDir);
% dircheck(imagingDir);

% dircheck(anatomicalDir);
% dircheck(freesurferDir);
% dircheck(caretDir);
% dircheck(gpCaretDir);   
% dircheck(regDir);
% dircheck(glmDir{1});

% set default plotting style
style.file(fullfile(codeDir,'pp1_style.m'));
style.use('default');

% ------------------------- Experiment Info -------------------------------
TR_length  = 1.5;  % seconds
numDummys  = 2;  % dummy images at the start of each run (these are discarded)
numTRs     = {[410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410]};  % total # of images per run (including dummies) per subj
run{1}     = {[1:6],...
             [1:5],...
             [1:5],...
             [1:6],...
             [1:5],...
             [1:5]};
run{2}     = {[7:11],...
             [6:11],...
             [6:11],...
             [7:11],...
             [6:11],...
             []};
% ------------------------- ROI things ------------------------------------
hem        = {'lh','rh'};                                                   % left & right hemi folder names/prefixes
regname    = {'sB3a','sB3b','sB1','sB2','S1','M1','B1','B2','B3','B4','B6','TH'}; % Cortical ROIs, 5 = S1, 6 = M1; Thalamus = 12                                             % roi names, independent of hemisphere    
regSide    = [ones(size(regname)),...                                       % Hemisphere of the roi
                ones(size(regname)).*2];                                    % [1 = left hemi (contra), 2 = right hemi (ipsi)]
regType    = [1:length(regname),...                                         % roi # (within hemisphere)
                1:length(regname)];
numregions = max(regType);                                                  % total number of regions 
% title of regions ordered according to numerical call id (eg. 2 = Lh M1)
reg_title  = {'left sB3a','left sB3b','left sB1','left sB2','left S1','left M1','left B1','left B2','left B3','left B4','left B6','left TH',...
            'right sB3a','right sB3b','right sB1','right sB2','right S1','right M1','right B1','right B2','right B3','rh B4','right B6','right TH'}; 

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

% ------------------------- Freesurfer things -----------------------------         
atlasA    = 'x';                                                            % freesurfer filename prefix
atlasname = 'fsaverage_sym';                                                % freesurfer average atlas
hemName   = {'LeftHem','RightHem'};                                         % freesurfer hemisphere folder names    

% ------------------------- Subject things --------------------------------
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
%   Step 1: get .nii file of anatomical data by running "spmj_tar2nii(TarFileName,NiiFileName)"
%   Step 2: open .nii file with MRIcron and manually find AC and read the xyz coordinate values
%           (note: there values are not [0 0 0] in the MNI coordinate)
%   Step 3: set those values into loc_AC (subtract from zero)
subj_name  = {'pd01','pd02','s01','s02','s03','s04'};
numSess    = [2,2,2,2,2,1];
anatNum{1} = {[],[40],[],[],[27:31],[25:29]};
anatNum{2} = {[],[],[29,32],[20:24],[],[]};
loc_AC     = {[-78  -122 -126],...
              [-80 -128 -133],...
              [-80 -122 -123],...
              [-82 -116 -132],...
              [-81 -120 -135],...
              [-84 -112 -122]};
          
DicomName{1}  = {'2019_03_01_PD01.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_12_pd02_sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_26_PP1_S01.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_15_S02_sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_22_S03_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_24_S04_Sess1.MR.Diedrichsen_PassivePatterns'};
DicomName{2}  = {'2019_04_22_PD01_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_13_pd02_sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_28_PP1_S01_SESS2.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_16_S02_sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_24_S03_Sess2.MR.Diedrichsen_PassivePatterns',...
                 ''};
             
NiiRawName{1} = {'2019_03_01_PD01',...
                 '2019_03_12_pd02_sess1',...
                 '2019_03_26_PP1_S01',...
                 '2019_04_15_S02_sess1',...
                 '2019_04_22_S03_Sess1',...
                 '2019_04_24_S04_Sess1'};
NiiRawName{2} = {'2019_04_22_PD01_Sess2',...
                 '2019_03_13_pd02_sess2',...
                 '2019_03_28_PP1_S01_SESS2',...
                 '2019_04_16_S02_sess2',...
                 '2019_04_24_S03_Sess2',...
                 ''};
             
fscanNum{1}   = {[21,24,27,30,33,36],...
                 [24,27,30,33,36],...
                 [14,17,20,23,26],...
                 [13,15,17,19,21,23],...
                 [12,14,16,18,20],...
                 [12,18,20,22,24]};   
fscanNum{2}   = {[12,14,16,18,20],...
                 [12,15,18,21,24,27],...
                 [12,14,16,18,20,22],...
                 [29,31,35,37,39],...
                 [12,14,16,18,20,22],...
                 []};  
             
fieldNum{1}   = {[],...
                 [43,44],...
                 [31:34],...
                 [26:29],...
                 [23:26],...
                 [36,37]};                                               
fieldNum{2}   = {[23:26],...
                 [29,30],...
                 [25:28],...
                 [42:45],...
                 [25:28],...
                 []};
             
dataPrefix    = {'r','u','u','u','u','u'};
             
% ------------------------- Analysis Cases --------------------------------
switch(what)
    case 'LIST_subjs'               
        D = dload(fullfile(baseDir,'subj_info.txt'));
        
        if nargout==0
            fprintf('\nSN\tID\t\tAge\tGender\tHandedness\tNum_FMRI_Sessions');
            fprintf('\n--\t--\t\t---\t------\t----------\t-----------------\n');
            for s = unique(D.SN)'
                S = getrow(D,ismember(D.SN,s));
                fprintf('%s\t%s\t\t%d\t%d\t%s\t\t%d',S.SN{1},S.ID{1},S.age,S.gender,S.handedness{1},S.fmri_sessions);
                fprintf('\n');
            end;
            fprintf('\n');
        else
            varargout = {D};
        end;
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
        fig   = [];
        vararginoptions(varargin(2:end),{'split','label','fig'});
        
        if isempty(fig); figure('Color',[1 1 1]); else; fig; end
        
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
                MOV     = movload(fullfile(behavDir,sprintf('%s_%s_%02d.mov',filePrefix,subj_name{sn},r)));
                for i = trials' 
                    d = getrow(D,D.TN==i & D.BN==r);
                    t = pp1_fmri_trial(MOV{1,i},d,0);
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
            save(dataFileOut,'T');
        end
    case 'SUBJ_plotStimForce'
        % plots max stimulation force per finger 
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % load data
        load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T
        D = tapply(T,{'run','numDigits','finger','stimulated'},{'peakF_raw','nanmean'},{'peakF_filt','nanmean'},{'time_stimOnset','nanmean'});
        % plot
        style.use('5fingers');
        
        subplot(1,3,1);
        plt.line(D.numDigits,D.peakF_filt,'split',D.finger,'errorfcn','std','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 4]);
        plt.labels('# fingers stimulated','force (N)',sprintf('stimulated\ntrials'));
        drawline(unique(T.forceStim),'dir','horz','linestyle',':','linewidth',1.5);
        legend off
        
        subplot(1,3,2);
        plt.line(D.numDigits,D.peakF_raw,'split',D.finger,'errorfcn','std','subset',D.stimulated==0);
        plt.set('xlim',[0.5 4.5],'ylim',[0 4]);
        plt.labels('# fingers stimulated','force (N)',sprintf('non-stimulated\ntrials'));
        drawline(unique(T.forceStim),'dir','horz','linestyle',':','linewidth',1.5);
        plt.legend('east',{'thumb','index','middle','ring','little'});
        legend off
        
        subplot(1,3,3);
        plt.line(D.numDigits,D.time_stimOnset,'split',D.finger,'errorfcn','std','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 1000]);
        plt.labels('# fingers stimulated','stimulus onset (ms)',sprintf('stimulus onset\ntime'));
        plt.legend('east',{'thumb','index','middle','ring','little'});
        
        varargout = {D};
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
        sn = 1:length(subj_name);
        vararginoptions(varargin,{'sn'});
        % get data
        D = [];
        for s = sn
            load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s}))); % T
            D = addstruct(D,T);
        end
        D = tapply(D,{'sn','numDigits','finger','stimulated'},{'peakF_raw','nanmean'},{'peakF_filt','nanmean'},{'time_stimOnset','nanmean'});
        
        % plot
        style.use('5fingers');
        
        subplot(1,3,1);
        plt.line(D.numDigits,D.peakF_filt,'split',D.finger,'errorfcn','stderr','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 4]);
        plt.labels('# fingers stimulated','force (N)',sprintf('stimulated\ntrials'));
        drawline(unique(T.forceStim),'dir','horz','linestyle',':','linewidth',1.5);
        legend off
        
        subplot(1,3,2);
        plt.line(D.numDigits,D.peakF_raw,'split',D.finger,'errorfcn','stderr','subset',D.stimulated==0);
        plt.set('xlim',[0.5 4.5],'ylim',[0 4]);
        plt.labels('# fingers stimulated','force (N)',sprintf('non-stimulated\ntrials'));
        drawline(unique(T.forceStim),'dir','horz','linestyle',':','linewidth',1.5);
        plt.legend('east',{'thumb','index','middle','ring','little'});
        legend off
        
        subplot(1,3,3);
        plt.line(D.numDigits,D.time_stimOnset,'split',D.finger,'errorfcn','stderr','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 1000]);
        plt.labels('# fingers stimulated','stimulus onset (ms)',sprintf('stimulus onset\ntime'));
        plt.legend('east',{'thumb','index','middle','ring','little'});
        
        varargout = {D};
    
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
                    pp1_imana('PREP_make_4dNifti','sn',s);
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
            end
            display(sprintf('Series %d done \n',seriesNum{sess}{sn}(i)))
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
    case 'WRAPPER_SURF'                                                 
        vararginoptions(varargin,{'sn'});
        % You can call this case to do all the freesurfer processing.
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
                          [4.9406 13.158 d_hrf(3) d_hrf(4) 2.6733 d_hrf(6) d_hrf(7)]};
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
    case 'WRAPPER_SEARCH'                                               
        glm = 2;
        vararginoptions(varargin,{'sn','glm'});
        % You can call this case to do all the searchlight analyses.
        % 'sn' can be an array of subjects because each processing case
        % contained within loops through the subject array submitted to the
        % case.
        for s = sn 
            for g = glm
                pp1_imana('SEARCH_define','sn',s,'glm',g);
                pp1_imana('SEARCH_run_LDC','sn',s,'glm',g);
                pp1_imana('SEARCH_mapSurfLDC','sn',s,'glm',glm);
            end
        end
        %fmri_imana('SEARCH_group_make');   % require group data
        %fmri_imana('SEARCH_group_cSPM');
    case 'SEARCH_define'                                                    % STEP 4.1   :  Defines searchlights for 120 voxels in grey matter surface
        glm = 2;
        vararginoptions(varargin,{'sn','glm'});
        
        mask       = fullfile(glmDir{glm},subj_name{sn},'mask.nii');
        Vmask      = spm_vol(mask);
        Vmask.data = spm_read_vols(Vmask);

        LcaretDir = fullfile(caretDir,sprintf('x%s',subj_name{sn}),'LeftHem');
        RcaretDir = fullfile(caretDir,sprintf('x%s',subj_name{sn}),'RightHem');
        white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
        pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
        topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};
        S         = rsa_readSurf(white,pial,topo);

        L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[15 120]);
        save(fullfile(anatomicalDir,subj_name{sn},sprintf('%s_searchlight_120.mat',subj_name{sn})),'-struct','L');
    case 'SEARCH_run_LDC'                                                   % STEP 4.2   :  Runs LDC searchlight using defined searchlights (above)
        % Requires java functionality unless running on SArbuckle's
        % computer.
        glm  = 2;
        vararginoptions(varargin,{'sn','glm'});
        
        if glm==2
            numConds = 31; % no errors modeled
        elseif glm==3
            numConds = 32; % error modeled
        elseif glm==4
            numConds = 32; % foot movement modeled
        end
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        % go to subject's glm directory 
        cd(fullfile(glmDir{glm},subj_name{sn}));
        % load their searchlight definitions and SPM file
        L = load(fullfile(anatomicalDir,subj_name{sn},sprintf('%s_searchlight_120.mat',subj_name{sn})));
        load SPM;
        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{sn}));
        name = sprintf('%s_glm%d',subj_name{sn},glm);
        % make index vectors
        conditionVec  = kron(ones(numel(SPM.Sess),1),[1:numConds]');
        partitionVec  = kron([1:numel(SPM.Sess)]',ones(numConds,1));
        % run the searchlight
        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partitionVec,'analysisName',name,'idealBlock',block);
        cd(cwd);
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

        sn  = [1,3:7];
        vararginoptions(varargin,{'sn'});

        linedef = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
        tmpMaskFile = {};  % for volumetric sub-cortical regions (these tmp masks are deleted after mapping to subject)
        for s=sn
            caretSubjDir = fullfile(caretDir,[atlasA subj_name{s}]);
            mask         = fullfile(glmDir{1},subj_name{s},'mask.nii,1');  
            for h = 1:2 
                D1 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI_3.paint']));   % sBa1 sBa2 sB3a sB3b
                D2 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI.paint']));     % premade ROIs from fsaverage_sym: M1 and S1
                D3 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROIsm.paint']));   % premade ROIs from fsaverage_sym: brodmann areas
                for r = 1:numregions
                    if r<5; D = D1; rr = r;         % premade ROIs from fsaverage_sym: brodmann areas
                    elseif r<7 D = D2; rr = r-4;    % premade ROIs from fsaverage_sym: S1 and M1
                    else D = D3; rr = r-6; end      % sBa1 sBa2 sB3a sB3b: brodmann areas cut with S1 boundary
                    idx = r+(h-1)*numregions;
                    % make R region structure for participant
                    R{idx}.name  = [subj_name{s} '_' regname{r} '_' hem{h}];
                    if r<12 % cortical surface regions
                        R{idx}.type     = 'surf_nodes';
                        R{idx}.location = find(D.data(:,1)==rr);
                        R{idx}.white    = fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                        R{idx}.pial     = fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                        R{idx}.topo     = fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                        R{idx}.linedef  = linedef;
                        R{idx}.flat     = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                        R{idx}.image    = mask; % functional mask
                    elseif r==12 % thalamus volumetric region
                        if h==1; rr = 10; elseif h==2; rr = 49; end % 10 is left thalamus, 49 is right thalamus
                        R{idx}.type     = 'roi_image';
                        R{idx}.file     = fullfile(scAnatDir,subj_name{s},sprintf('%s_%d_%d.nii',subj_name{s},h,rr));
                        R{idx}.value    = 1;
                        R{idx}.image    = fullfile(glmDir{1},subj_name{s},sprintf('mask_tmp%d%d.nii',h,rr)); % functional mask
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
            R = region_calcregions(R,'exclude',[2,6; 5,6],'exclude_thres',0.75);
            cd(regDir);
            save(['regions_' subj_name{s} '.mat'],'R');
            fprintf('\n %s done\n',subj_name{s})
            clear R
        end
        for i = 1:numel(tmpMaskFile) 
            delete(tmpMaskFile{i});
        end
    case 'ROI_makePaint'                                                    % STEP 5.1   :  Make paint file for ROIs (saves as ROI_2.paint)
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
        P    = caret_load([hem{h} '.lobes.paint']);                  % paint file of 4 major lobes (1 = col 2 of paintnames)
        % A  = caret_load('ROI.paint');                             % ROI.paint file from ??- file has nice, standard ROIs for motor, sensory, aux. rois, and some visual rois: Naveed maybe made this?
        C    = caret_load([hem{h} '.FLAT.coord']);                   % Caret flatmap (for X Y coordinate restrictions)
		M    = caret_load([hem{h} '.propatlas.metric']);             % probabilistic atlas
		% Load searchlight metrics
		Avg  = caret_load([hem{h} '.group_avg_3_stats.metric']);     % group surface searchlight results (FAST glm)- data(:,10) is group avg LDC
		Dgt  = caret_load([hem{h} '.group_digit_3_stats.metric']);   % group surface searchlight results, contrasted for distance b/t fingers, independent of pressing speed
		Spd  = caret_load([hem{h} '.group_speed_3_stats.metric']);   % group surface searchlight results, contrasted for distance b/t pressing speeds, independent of fingers
        % Make index for voxels that are speed-specific or digit-specific from Dgt and Spd metrics (loaded above)
		Didx = logical(Dgt.data(:,10)>Spd.data(:,10));				 % index for voxels that demonstrate greater sensitivity (distance) between fingers than speeds
		Sidx = logical(Spd.data(:,10)>Dgt.data(:,10));		         % index for voxels that demonstrate greater sensitivity (distance) between speeds than fingers
		
        % - - - - - - - - - - Assign roi labels to vertices - - - - - - - -
        % get data for new rois (from probability atlas)
        M1  = sum(M.data(:,[7,8]),2); % M1 is BA4a + 4p
        S1  = sum(M.data(:,[1:4]),2); % S1 is BA1 + 2 + 3a + 3b
        MT  = M.data(:,10);         
        SMA = M.data(:,9);           
        V1  = M.data(:,12);
        V2  = M.data(:,13);
        Ba1 = M.data(:,1);
        Ba2 = M.data(:,2);
        B3a = M.data(:,3);
        B3b = M.data(:,4);
        ALL = ones(size(M.data,1),1);       % takes all vertices available (full cortex roi)- is for UWO BrainHack dataset, not analyses
        
        % coordinate is ROI w/ for which it has greatest associated probability
        [Prop,ROI]   = max([S1 M1 SMA MT V1 V2],[],2); 
                    % ROI-> 1...2..3..4..5...6  
        [Prop3,ROI3] = max([Ba1 Ba2 B3a B3b],[],2);
                    % ROI->  7...8...9..10

        % Define ROIS with:
        %...cytoarchitectonic prob (>0.2) - - - - - - - - - - - - - - - - -
            ROI(Prop<0.2)   = 0;
            ROI3(Prop3<0.2) = 0;
        %...boundaries of 4 major lobes (F.P.O.T.)- - - - - - - - - - - - -
            ROI(ROI==1 & P.data~=2 & Prop<0.25)    = 0;  % sS1  :  parietal and higher prob
            ROI(ROI==2 & P.data~=1 & Prop<0.25)    = 0;  % sM1  :  frontal and higher prob
            ROI(ROI==3 & P.data~=1)                = 0;  % sSMA :  frontal
            ROI(ROI==4 & P.data~=2)                = 0;  % sMT  :  frontal
            ROI(ROI==5 & P.data~=3)                = 0;  % sV1  :  occipital
            ROI(ROI==6 & P.data~=3)                = 0;  % sV2  :  occipital
            ROI3(ROI3==1 & P.data~=2)              = 0;  % Ba1  :  parietal
            ROI3(ROI3==2 & P.data~=2)              = 0;  % Ba2  :  parietal
            ROI3(ROI3==3 & P.data~=2)              = 0;  % B3a  :  parietal
            ROI3(ROI3==4 & P.data~=2)              = 0;  % B3b  :  parietal

        %...average distances of all condition pairs (for SMA, MT, V1, and V2)
        % (and the flatmap coords for left hemi rois)
        switch h   % hemispheres
            case 1 % left  
                ROI(ROI==1 & Avg.data(:,10)<0.15)    = 0;  % sS1- these dissimiarlity cutoffs lead to ROIs that focus on hand knob M1 S1
                ROI(ROI==2 & Avg.data(:,10)<0.15)    = 0;  % sM1
                ROI(ROI==3 & Avg.data(:,10)<0.1)     = 0;  % sSMA
                ROI(ROI==4 & Avg.data(:,10)<0.1)     = 0;  % sMT
                ROI(ROI==5 & Avg.data(:,10)<0.1)     = 0;  % sV1
                ROI(ROI==6 & Avg.data(:,10)<0.1)     = 0;  % sV2
                
                %...flatmap X Y coordinates - - - - - - - - - - - - - - - -
                % SMAcoords
                    % Y coords
                    ROI(ROI==3 & C.data(:,2)>42)   = 0;
                    ROI(ROI==3 & C.data(:,2)<25)   = 0;
                    % X coords
                    ROI(ROI==3 & C.data(:,1)<-42)  = 0;
                    ROI(ROI==3 & C.data(:,1)>-21)  = 0;   
                % M1 (little blob near temporal lobe that still persists)
                    ROI(ROI==2 & C.data(:,2)<-1.5) = 0;
				
				% make speed or finger-sensitive rois within already defined sS1 and sM1.
                ROI4 = ROI;			% for making digit and speed- specific voxel selections within S1 and M1
                ROI4(ROI4==3 | ROI4==4 | ROI4==5 | ROI4==6) = 0;
				ROI4(ROI==1 & Didx) 			   = 1;	   % sS1 :  finger-sensitive
				ROI4(ROI==1 & Sidx)                = 2;	   % sS1 :  speed-sensitive
				ROI4(ROI==2 & Didx)                = 3;	   % sM1 :  |finger-sensitive
				ROI4(ROI==2 & Sidx)			       = 4;    % sM1 :  speed-sensitive
            case 2 % right 
                ROI(ROI==1 & Avg.data(:,10)<0.08)    = 0;  % sS1
                ROI(ROI==2 & Avg.data(:,10)<0.07)    = 0;  % sM1
                %ROI(ROI==3 & Avg.data(:,10)<0.001)=0;     % sSMA- most voxels don't survive distance criteria in right hemi
                ROI(ROI==4 & Avg.data(:,10)<0.05)    = 0;  % sMT
                ROI(ROI==5 & Avg.data(:,10)<0.05)    = 0;  % sV1
                ROI(ROI==6 & Avg.data(:,10)<0.05)    = 0;  % sV2
        end;
        
        % - - - - - - - - - - - Save ROIpaint files - - - - - - - - - - - -  
        % Save Paint files for first six ROIs - - - - - - - - - - - - - - -
        areas{1}{1}  = {'sS1'}; 
        areas{2}{1}  = {'sM1'}; 
        areas{3}{1}  = {'sSMA'}; 
        areas{4}{1}  = {'sMT'};
        areas{5}{1}  = {'sV1'}; 
        areas{6}{1}  = {'sV2'}; 
        names        = {'sS1','sM1','sSMA','sMT','sV1','sV2'};
        
        colors=[255 0 0;...  % sS1
            0 255 0;...      % sM1
            0 0 255;...      % sSMA
            160 0 160;...    % sMT
            75 150 0;...     % sV1
            200 0 100];      % sV2
            
        Paint = caret_struct('paint','data',ROI,'paintnames',names,'column_name',{'ROI'});
        caret_save(['ROI_2.paint'],Paint);
        caret_combinePaint('ROI_2.paint','ROI_2.paint','ROI_2.areacolor',...
            'areas',areas,'names',names,'colors',colors);
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Save paint file for Ba1:3b (cannot index same vertex for many
        % rois in same file..or at least not sure how to do so
        % currently)-SA
        clear areas colors names Paint
        areas{1}{1} = {'sBa1'}; 
        areas{2}{1} = {'sBa2'}; 
        areas{3}{1} = {'sB3a'}; 
        areas{4}{1} = {'sB3b'};
        
        names = {'sBa1','sBa2','sB3a','sB3b'};
        
        colors = [0 100 200;... % Ba1
            255 255 0;...       % Ba2
            150 50 255;...      % Ba3
            0 153 153];         % Ba3b
        
        Paint = caret_struct('paint','data',ROI3,'paintnames',names,'column_name',{'ROI'});
        caret_save(['ROI_3.paint'],Paint);
        caret_combinePaint('ROI_3.paint','ROI_3.paint','ROI_3.areacolor',...
            'areas',areas,'names',names,'colors',colors);
        
		% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Save finger and speed-sensitive M1 and S1 rois
        clear areas colors names Paint
        areas{1}{1} = {'sS1_digit'}; 
        areas{2}{1} = {'sS1_speed'}; 
        areas{3}{1} = {'sM1_digit'}; 
        areas{4}{1} = {'sM1_speed'};
        
        names = {'sS1_digit','sS1_speed','sM1_digit','sM1_speed'};
        
        colors = [0 100 200;... % sS1_digit
            255 255 0;...       % sS1_speed
            150 50 255;...      % sM1_digit
            0 153 153];         % sM1_speed
        
        Paint = caret_struct('paint','data',ROI4,'paintnames',names,'column_name',{'ROI'});
        caret_save(['ROI_4.paint'],Paint);
        caret_combinePaint('ROI_4.paint','ROI_4.paint','ROI_4.areacolor',...
            'areas',areas,'names',names,'colors',colors);
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Save full cortex roi (used to harvest full brain dataset for
        % brainhack event).
        clear areas colors names Paint
        areas{1}{1} = {'Cortex'};
        names       = {'Cortex'};
        colors      = [255 0 0];
        Paint       = caret_struct('paint','data',ALL,'paintnames',names,'column_name',{'ROI'});
        caret_save(['ROI_all.paint'],Paint);
        caret_combinePaint('ROI_all.paint','ROI_all.paint','ROI_all.areacolor',...
            'areas',areas,'names',names,'colors',colors);
        
    end;    
    
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
        glm = 2;
        sn  = 1;
        roi = [1:8,13:20];
        append = 0; % just add betas to currently existing datastructure
        vararginoptions(varargin,{'sn','glm','roi','append'});
        excludeROIs = [];%[12,24];
        T=[];
        if append
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        end
        % harvest
        for s=sn % for each subj
            fprintf('\nSubject: %d\n',s) % output to user
            
            % load files
            load(fullfile(glmDir{glm}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            D = load(fullfile(glmDir{glm}, subj_name{s},'SPM_info.mat'));
            load(fullfile(regDir,sprintf('regions_%s.mat',subj_name{s})));          % load subject's region parcellation & depth structure (R)
            
            % add percent signal change imgs for subject
            Q = {}; 
            if glm==2
                numImgs = 31; 
            elseif glm==3
                numImgs = 32; 
            end
            for q = 1:numImgs
                Q{q} = (fullfile(glmDir{glm}, subj_name{s}, sprintf('psc_%02d.nii',q))); 
            end
            Q = spm_vol(char(Q));
            
            % TR img info
            V = SPM.xY.VY; 

            % remove run means from patterns
            C0   = indicatorMatrix('identity',D.run); 
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
    case 'ROI_stats'                                                        % STEP 5.4   :  Calculate stats/distances on activity patterns
        glm = 2;
        sn  = 1:5;
        vararginoptions(varargin,{'sn','glm'});
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
        for s = sn % for each subject
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
                % squared dissimilarities
                So.ldc_wmean = rsa.distanceLDC(betaW,D.run,D.tt);        % rdm crossvalidated, on patterns without run mean patterns removed
                So.ldc       = rsa.distanceLDC(betaW_nmean,D.run,D.tt);  % rdm crossvalidated, patterns with run means removed
                % PSC
                So.psc       = nanmean(S.psc{1},2)';
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
   
    case 'ROI_pattConsist'                                                  % (optional) :  Calculates pattern consistencies for each subject in roi across glms.
        % Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % This stat is useful for determining which GLM model yeilds least
        % variable beta estimates. That GLM should be the one you use for
        % future analysis cases.

        glm = 4;
        sn  = 1;
        roi = 5; % default contra S1 (lh)
        removeMean = 1;
        betaType = 3; % 1: raw, 2: uni-whitened, 3: multi-whitened
        vararginoptions(varargin,{'sn','glm','roi','removeMean','betaType'});
        
        % consistencies only calculate for chord regressors, regardless of
        % any other foot/error regressors included in model.
        
        R = []; % output structure

        for g = glm
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',g))); % loads in struct 'T'
            % number of conditions per glm
            if g==2
                %conds = [1:31]';
                conds = [1:5,zeros(1,26)]';
            elseif g==3
                %conds = [1:31,0]';
                conds = [1:5,zeros(1,26),0]';
            elseif g==4
                %conds = [1:31,0]';
                conds = [1:5,zeros(1,26),0]';
            end
            
            for r = roi

                for s = sn
                    
                    % % Calculate pattern consistency (avg. across
                    % conditions) for each roi, each subj, each glm:
                    % - get corresponding data from beta structure
                    % - index patterns for conditions and runs
                    % - exclude any foot/error/intercept regressors
                    % - calculate non-CV pattern consistency
                    
                    S = getrow(T,(T.sn==s & T.roi==r));
                    runs = unique(S.run{1});
                    switch betaType 
                        case 1 % raw betas
                            betas = S.raw_beta{1};
                        case 2 % univariately prewhitened betas
                            betas = S.betaUW{1};
                        case 3 % multivariately prewhitened betas
                            betas = S.betaW{1};
                    end
                    % make vectors for pattern consistency func
                    condVec = kron(ones(numel(runs),1),conds);
                    partVec = kron(runs',ones(length(conds),1));
                    partVec(condVec==0) = 0;
                    % calculate the pattern consistency
                    rs.r2  = rsa_patternConsistency(betas,partVec,condVec,'removeMean',removeMean);
                    rs.sn  = s;
                    rs.roi = r;
                    rs.glm = g;
                    rs.betaType = betaType;
                    rs.removeMean = removeMean;
                    R = addstruct(R,rs);
                end
            end
        end
        pivottable(R.glm,R.sn,R.r2,'mean','numformat','%0.4f');
        varargout = {R};
        % output arranged such that each row is an roi, each col is subj
        
        %_______________
    case 'ROI_pattConsistCV'                                                % (optional) :  Calculates pattern consistencies for each subject in roi across glms.
        % Crossvalidated Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % This stat is useful for determining which GLM model yeilds least
        % variable beta estimates. That GLM should be the one you use for
        % future analysis cases.

        glm = 4;
        sn  = 1;
        roi = 5; % default contra S1 (lh)
        removeMean = 1;
        betaType = 3; % 1: raw, 2: uni-whitened, 3: multi-whitened
        vararginoptions(varargin,{'sn','glm','roi','removeMean','betaType'});
        
        % consistencies only calculate for chord regressors, regardless of
        % any other foot/error regressors included in model.
        
        R = []; % output structure

        for g = glm
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',g))); % loads in struct 'T'
            % number of conditions per glm
            if g==2
                %conds = [1:31]';
                conds = [1:5,zeros(1,26)]';
            elseif g==3
                %conds = [1:31,0]';
                conds = [1:5,zeros(1,26),0]';
            elseif g==4
                %conds = [1:31,0]';
                conds = [1:5,zeros(1,26),0]';
            end
            
            for r = roi

                for s = sn
                    
                    % % Calculate pattern consistency (avg. across
                    % conditions) for each roi, each subj, each glm:
                    % - get corresponding data from beta structure
                    % - index patterns for conditions and runs
                    % - exclude any foot/error/intercept regressors
                    % - calculate non-CV pattern consistency
                    
                    S = getrow(T,(T.sn==s & T.roi==r));
                    runs = unique(S.run{1});
                    switch betaType 
                        case 1 % raw betas
                            betas = S.raw_beta{1};
                        case 2 % univariately prewhitened betas
                            betas = S.betaUW{1};
                        case 3 % multivariately prewhitened betas
                            betas = S.betaW{1};
                    end
                    % make vectors for pattern consistency func
                    condVec = kron(ones(numel(runs),1),conds);
                    partVec = kron(runs',ones(length(conds),1));
                    partVec(condVec==0) = 0;
                    % calculate the pattern consistency
                    [rs.r2,rs.r] = rsa_patternConsistency_crossval(betas,partVec,condVec,'removeMean',removeMean);
                    rs.sn  = s;
                    rs.roi = r;
                    rs.glm = g;
                    rs.betaType = betaType;
                    rs.removeMean = removeMean;
                    R = addstruct(R,rs);
                end
            end
        end
        pivottable(R.glm,R.sn,R.r2,'mean','numformat','%0.4f');
        pivottable(R.glm,R.sn,R.r,'mean','numformat','%0.4f');
        varargout = {R};
        % output arranged such that each row is an roi, each col is subj
        
        %_______________
    
    case 'ROI_dist_stability'                                               % (optional) :  Correlate distances of conditions in roi across subjects. 
        glm = 1;
        roi = 2; % default primary motor cortex
        % Correlate distances of conditions in roi across subjs.
        % Regression line is not forced through origin so ignores distance 
        % scaling across subjects.
        vararginoptions(varargin,{'roi','glm'});
        
        D   = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));      
        D   = getrow(D,D.region==roi);
        Cs  = corr(D.RDM');
        
        varargout = {Cs};
    case 'ROI_MDS_overall'                                                  % (optional) :  Plots the scaled representational structure. 
        % enter region, glm #, sn (if desired)
        cplot = 'one';
        glm   = 2;
        roi   = 5; % default primary sensory cortex   
        sn    = 1;
        clrCode = 0; % if >0, color all chords with digit X red, and all other chords black.
        fig = [];
        vararginoptions(varargin,{'roi','glm','cplot','clrCode','sn','fig'});
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
                pp1_imana('MISC_scatterplotMDS',Y{1}(1:end,1:3),'split',split,'label',[1:31]','fig',fig);
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
        sn  = 2;
        glm = 2;
        roi = 5;
        vararginoptions(varargin,{'roi','glm','sn'});
        
        conds = 1:31;
        mean_sub = 1;
        betaType = 'raw';
        
        R = [];
        for s = sn
            r_s1 = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'split','sess1','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            r_s2 = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'split','sess2','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            r_A  = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'split','sess_cv_odd','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
            r_B  = pp1_imana('ROI_patternReliability','sn',s,'roi',roi,'split','sess_cv_even','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
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
        varargout = {R};
        
        
    case '0' % ------------ HARVEST: condense/make new datastructures. ----
        % The STATS cases usually harvest some data from ROI stats structures.
        % These cases usually correspond to FIG cases of the same name.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'GET_idxPerNumDigits'
        glm = 2;
        numDigits = 1:5;
        vararginoptions(varargin,{'glm','numDigits'});
        
%         if numel(numDigits)>1
%             error('NOTE: numDigits is exclusive! Cannot provide range.')
%         end
        
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
        for i = 1:size(T.sn,1);
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
        sn           = 1;
        glm          = 4;
        roi          = 5;
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
    
    case 'old_FIG_pattReliabilitySingleFinger'
        % plot split-half reliabilities for patterns evoked by single finger
        % chords.
        sn  = 1;
        glm = 4;
        roi = 5;
        fig = [];
        betaType = 'raw';
        vararginoptions(varargin,{'sn','glm','roi','fig','betaType'});
        
        % get split-half pattern reliabilities
        D = pp1_imana('ROI_splithalfPattReliability','sn',sn,'glm',glm,'roi',roi,'betaType',betaType,'conds',1:5);
        % plot
        if isempty(fig); figure('Color',[1 1 1]); else; fig; end
        style.use('5shades');
        plt.box(D.type,D.corr,'split',D.roi,'plotall',2);
        %plt.set('xticklabel',{'within conds','across conds'});
        plt.labels('','pearson''s r');
        drawline(0,'dir','horz');
        
        varargout = {D};
    case 'old_FIG_pattReliabilityAllChords'
        % plot split-half reliabilities for patterns evoked by single finger
        % chords.
        sn  = 1;
        glm = 4;
        roi = 5;
        fig = [];
        betaType = 'raw';
        vararginoptions(varargin,{'sn','glm','roi','fig','betaType'});
        
        % get split-half pattern reliabilities
        D = pp1_imana('ROI_splithalfPattReliability','sn',sn,'glm',glm,'roi',roi,'betaType','raw','conds',1:31);
        % plot
        if isempty(fig); figure('Color',[1 1 1]); else; fig; end
        style.use('5shades');
        plt.box(D.type,D.corr,'split',D.roi,'plotall',2);
        %plt.set('xticklabel',{'within conds','across conds'});
        plt.labels('','pearson''s r');
        drawline(0,'dir','horz');
        
        varargout = {D};
    
    case 'PLOT_pattReliability'
        % plots pattern reliability w/in session and across session for
        % sispecified conditions (default = all)
        sn    = 2;
        glm   = 2;
        roi   = 5;
        conds = 1:31; 
        vararginoptions(varargin,{'sn','roi','glm','conds'});
        
        % sess 1 odd-even reliability
        R1 = pp1_imana('ROI_patternReliability','sn',sn,'glm',glm,'roi',roi,'conds',conds,'split','sess1');
        R1 = tapply(R1,{'roi','type'},{'corr','mean'});
        R1.split = [1;1];
        % sess 2 odd-even reliability
        R2 = pp1_imana('ROI_patternReliability','sn',sn,'glm',glm,'roi',roi,'conds',conds,'split','sess2');
        R2 = tapply(R2,{'roi','type'},{'corr','mean'});
        R2.split = [2;2];
        % across-session reliability
        R3 = pp1_imana('ROI_patternReliability','sn',sn,'glm',glm,'roi',roi,'conds',conds,'split','sess');
        R3 = tapply(R3,{'roi','type'},{'corr','mean'});
        R3.split = [3;3];
        % plot
        R = addstruct(R1,R2);
        R = addstruct(R,R3);
        plt.bar([R.split,R.type],R.corr,'style',style.custom({'black','red'}));
        plt.labels('split','pearson''s r');
        
    case 'FIG_PSCperNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn    = 1;
        glm   = 2;
        roi = [1:4]; % [3a, 3b, 1, 2]
        fig = [];
        vararginoptions(varargin,{'sn','roi','glm','fig'});
        
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
        %D  = getrow(D,D.numDigits<6);
        Dr = tapply(D,{'sn','roi','numDigits'},{'psc','mean'});
        Dr = getrow(Dr,ismember(Dr.roi,roi) & ismember(Dr.sn,sn));
        % plot
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        style.use('numDigits');
        plt.box([Dr.roi Dr.numDigits],Dr.psc,'split',Dr.numDigits);
        xtick = {};
        for r = roi
            xtick = {xtick{:} '','',reg_title{r},'',''};
        end
        plt.set('xticklabel',xtick,'xticklabelrotation',45);
        plt.legend('east',{'1 digit','2 digits','3 digits','4 digits','5 digits'});
        plt.labels('region','percent signal change');
        drawline(0,'dir','horz');
        
        varargout = {Dr,D};
        
    case 'FIG_LDCperNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn    = 1;
        glm   = 4;
        roi = [1:4,6]; % [3a, 3b, 1, 2, M1]
        fig = [];
        vararginoptions(varargin,{'sn','roi','glm','fig'});
        
        D = pp1_imana('GET_LDCperNumDigits','sn',sn,'glm',glm,'roi',roi);
        D.ldc = ssqrt(D.ldc);
        D = tapply(D,{'sn','roi','numDigits'},{'ldc','mean'});
        % plot
        %style.use('5shades')
        sty = style.custom(plt.helper.get_shades(length(roi),'jet','decrease',10));
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.roi,D.ldc,'split',D.numDigits,'style',sty);
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
        
    case '0' % ------------ PCM: pcm analyses. ----------------------------
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
        M.name       = 'eye no scale';
        M.Ac(:,:,1)  = chords;
        varargout = {M};
    case 'pcm_linearScale'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns
        M.type       = 'feature';
        M.numGparams = 1;
        M.name       = 'eye L scale';
        M.Ac(:,:,1)  = pp1_simulations('chords');
        varargout = {M};
    case 'pcm_nonlinearScale'
        % Model nonLinear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates ONLY scaling params- no finger params
        M.type       = 'nonlinear';
        M.name       = 'eye nL scale';
        M.modelpred  = @pp1_modelpred_scale;
        M.numGparams = 4;
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_linearScaleNonlinearFingers'
        % Model nonLinear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates both finger params. Scaling params are fixed according
        % to the number of fingers in the chord.
        M.type       = 'nonlinear';
        M.name       = 'nL fingers L scale';
        M.modelpred  = @pp1_modelpred_singleFingerLinear;
        M.numGparams = 14;
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_nonlinearScaleNonlinearFingers'
        % Model nonLinear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates both finger params and scaling params.
        M.type       = 'nonlinear';
        M.name       = 'nL fingers nL scale';
        M.modelpred  = @pp1_modelpred_singleFingerScale;
        M.numGparams = 18;
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_nonlinearScaleLinearUsage'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates both scaling params with usage finger model.
        M.type       = 'nonlinear';
        M.name       = 'usage nL scale';
        M.modelpred  = @pp1_modelpred_usageScale;
        M.numGparams = 4;
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
    case 'PCM_fitModels_oneROI'
        % fits pcm models to data from one region
        sn     = 3;
        glm    = 2;
        roi    = 5;
        plotit = 0; % plot fits
        saveit = 0; % save fits
        runEffect = 'random';
        vararginoptions(varargin,{'sn','glm','roi','plotit','saveit'});
        % get data
        [Y,partVec,condVec,G_hat] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm);
        G_hat = mean(G_hat,3);
        
        % get pcm models
        M{1} = pp1_imana('pcm_null');
        M{end+1} = pp1_imana('pcm_noScale');                        % 2
        M{end+1} = pp1_imana('pcm_linearScale');                    % 3
        M{end+1} = pp1_imana('pcm_nonlinearScale');                 % 4
        M{end+1} = pp1_imana('pcm_nonlinearScaleLinearUsage');      % 5
        M{end+1} = pp1_imana('pcm_linearScaleNonlinearFingers');    % 6
        M{end+1} = pp1_imana('pcm_nonlinearScaleNonlinearFingers'); % 7
        M{end+1} = pp1_imana('pcm_freedirect'); 
        %M{end+1} = pp1_imana('pcm_freechol'); 
        % set starting values for nonlinear scaling model
        scaleParams = log([1.1, 1.2, 1.3, 1.4])';
        M{4}.theta0 = scaleParams;
        M{5}.theta0 = scaleParams;
        % set starting values for nonlinear finger scaling model (2nd moment for single fingers)
        Fx0 = pcm_free_startingval(G_hat(1:5,1:5));
        scaleParams = log([1.1, 1.2, 1.3, 1.4])';
        M{6}.theta0 = Fx0;
        M{7}.theta0 = [Fx0;scaleParams];
        
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
        sn  = 1:5;
        roi = [1:4,6,13:16,18];
        glm = 2;
        vararginoptions(varargin,{'sn','roi','glm'});
        for r = roi
            fprintf('roi: %d\n',r);
            pp1_imana('PCM_fitModels_oneROI','sn',sn,'roi',r,'glm',glm,'plotit',0,'saveit',1);
        end
    case 'PCM_getFits'
        % Gets models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        glm = 3;
        vararginoptions(varargin,{'glm'});
        D   = []; % output structure
        for r = 1:length(reg_title)
            % if exists, load pcm fits for region (otherwise, skip region)
            try
                load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,r)));
            catch
                continue
            end
            % scale likelihoods
            Tcv.likelihood_norm = bsxfun(@minus,Tcv.likelihood,Tcv.likelihood(:,1));
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
                d.likeNorm   = [nan(numModels-1,1); T.likelihood(j,end) - Tcv.likelihood(j,1)]; % upper noise ceiling
                d.thetaCV    = Q.thetaCV(Q.sn==j);
                D = addstruct(D,d);
            end
        end
        varargout = {D};
    case 'PCM_plotFitsBox'  
        % loads fit results per roi and plots them.
        glm = 3;
        roi = [1:4,6];
        sn  = 1:5;
        vararginoptions(varargin,{'glm','roi','sn'});
        % get pcm fits
        D = pp1_imana('PCM_getFits','glm',glm);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn));
        % plot pcm fits
        numModels  = length(unique(D.model));
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        figure('Color',[1 1 1]);
        for i = 1:numPlots
            r = roi(i);
            subplot(1,numPlots,i);
            sty = style.custom(plt.helper.get_shades(numModels-2,'gray','descend'));
            plt.box(D.model,D.likeNormCV,'subset',D.model>1 & D.model<numModels & D.roi==r,'split',D.model,'style',sty);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
            plt.labels('','relative log likelihood',reg_title{r});
            % plot noise ceilings
            drawline(mean(D.likeNormCV(D.model==numModels & D.roi==r)),'dir','horz','linestyle','-.');
            drawline(mean(D.likeNorm(D.model==numModels & D.roi==r)),'dir','horz','linestyle','-');
            legend off
            ylims = ylim;
            ylim([0 ylims(2)]);
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotFitsBar'  
        % loads fit results per roi and plots them.
        glm = 3;
        roi = [1:4,6];
        sn  = 1:5;
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
            Tcv.likelihood_norm = bsxfun(@minus,Tcv.likelihood,Tcv.likelihood(:,1));
            % plot fits (errorbars are stderr)
            subplot(1,numPlots,i);
            pcm_plotModelLikelihood(Tcv,M,'upperceil',T.likelihood(:,end),'style','bar','Nceil',numel(M));
            set(gca,'xticklabelrotation',45);
            ylabel('relative log-likelihood')
            title(reg_title{r});
            box off;
        end
        plt.match('y');
    case 'PCM_plotGpred'  
        % loads fit results per roi and plots them.
        glm = 2;
        roi = 4;
        sn  = 2:4;
        vararginoptions(varargin,{'glm','roi','sn'});
        if length(roi)>1
           error('can only call case with one roi'); 
        end
        % load fits
        load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,roi)));
        figure('Color',[1 1 1]);
        numPlots = numel(G_pred);
        for i = 1:numPlots
            % plot fits
            subplot(1,numPlots,i);
            imagesc(G_pred{i});
            title(sprintf('%s : %d params',M{i}.name,M{i}.numGparams));
        end    
    case 'PCM_plotThetas'
        % loads fit results per roi and plots them.
        glm   = 3;
        roi   = [1:4];
        sn    = 1:5;
        model = 6;
        vararginoptions(varargin,{'glm','roi','sn'});
        % load plotting-friendly data structure
        D = load(fullfile(pcmDir,sprintf('pcmModelFits_glm%d',glm)));
        D = getrow(D,D.model==model & ismember(D.roi,roi));
        D.theta_cv = cell2mat(D.theta_cv); 
        % plot
        sty = style.custom(plt.helper.get_shades(length(roi),'jet','descend'));
        plt.trace([1:14],D.theta_cv,'split',D.roi,'style',sty);

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
        
    case '0' % ------------ fingerpics: project patterns onto surface & take pictures.
        % You will absolutely need to edit these cases. 
        % These versions of the cases are frome ef1_imana
        % (extensionflexion).
        % 1 values of coordinates in xlims and ylims correspond to 1/10mm
        % distances.
        % So a range of 20 = 2mm distance on surface projection.
    case 'fingerpics'                                                       % Makes jpegs of finger activity patterns on cortical surface M1/S1
        sn  = 1;
        glm = 2;
        vararginoptions(varargin,{'sn','glm'});
        
        for s=sn;
            for g=glm;
                pp1_imana('surf_mapFingers','sn',s,'glm',g)
                pp1_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[1,5]);
                pp1_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[2]);
                pp1_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[3]);
                pp1_imana('surf_fingerpatterns','sn',s,'glm',g,'numDigits',[4]);
            end
        end
    case 'surf_mapFingers'                                                % Map locations of finger patterns- run after glm estimation
        % map volume images to metric file and save them in individual surface folder
        sn  = 1;
        glm = 2;
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
        glm = 4;
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
                xlims=[-5 18]; % may need to adjust locations for pics
                ylims=[-8 20];
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
                        'border',B.Border,'bordersize',10,'topo',topo,'coord',coord);
            colormap('jet');
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
         
    otherwise
        fprintf('%s: no such case.\n',what)

end


% % Local Functions
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

cd(cwd); % return user to their directory when function was called
end