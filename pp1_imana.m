function varargout = pp1_imana(what,varargin)
% function    varargout = pp1_imana(what,varargin)
%  
% passive patterns 1 fmri analyses
% 7T scanner with penumatic stimulation box. June 2018-ish to Dec. 2019
%
% SArbuckle 2019, UWO
% saarbuckle@gmail.com

%% ------------------------- Directories ----------------------------------
filePrefix = 'pp1_fmri'; % prefix for output files from experiment script
cwd        = cd; % get current directory when called, and return at end of script
% paths to project directories
codeDir         = '/Users/mahdiyarshahbazi/Documents/GitHub/project_passivePatterns';
baseDir         = '/Volumes/Diedrichsen_data$/data/passivePatterns_NoiseNormalisation';   % base directory for analysis
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
atlasDir        = '/Users/mahdiyarshahbazi/Documents/data/Atlas_templates/FS_LR_164';
caretDir        = [baseDir '/surfaceCaret'];     
gpCaretDir      = [caretDir '/fsaverage_sym'];
regDir          = [baseDir '/RegionOfInterest/'];   
pcmDir          = [baseDir '/PCM_models'];
glmDir          = {[baseDir '/glm1'],[baseDir '/glm2'],[baseDir '/glm3'],[baseDir '/glm4'],[baseDir '/glm5']};   
% glm4 = main analysis
% glm5 = modeled out all stimulation trials that preceeded a thumb press

% set default plotting style
style.file(fullfile(codeDir,'pp1_style.m'));
style.use('default');
%% ------------------------- ROI info -------------------------------------
hem        = {'lh','rh'};                                                   % left & right hemi folder names/prefixes
hemLetter  = {'L','R'};
hemName    = {'CortexLeft','CortexRight'};
%hemNum     = [1 2];
regname    = {'ba3A','ba3B','ba1','ba2','rM1','cM1','SII','S1','M1','SPLa','SPLp','VPL','MGN','LGN'};        % Cortical ROIs, 5 = S1, 6 = M1; Thalamus = 12                                             % roi names, independent of hemisphere    
cortical   = repmat([1 1 1 1 1 1 1 1 1 1 1 0 0 0],1,2);
thalamic   = repmat([0 0 0 0 0 0 0 0 0 0 0 1 1 1],1,2);
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
%% ------------------------- Exp info -----edit this for new subjs---------
trialDur  = 6; % duration (s) of trial from start to end
TR_length  = 1.5;  % seconds
numDummys  = 2;  % dummy images at the start of each run (these are discarded)
numTRs     = {[410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410],...
              [410,410,410,410,410,410,410,410,410,410,410]};  % total # of images per run (including dummies) per subj
run{1}     = {[1:6],...
             [1:5],...
             [1:5],...
             [1:6],...
             [1:5],...
             [1:5],...
             [1:6],...
             [1:6],...
             [1:5],...
             [1:4],...
             [1:2]};
run{2}     = {[7:11],...
             [6:11],...
             [6:11],...
             [7:11],...
             [6:11],...
             [6:11],...
             [7:11],...
             [7:11],...
             [6:11],...
             [5:11],...
             [3:5]};
% some participants needed to readjust, so this is what these extra
% "sessions" are from (before readjustment, we acquired fieldmaps)
run{3}     = {[],[],[],[],[],[],[],[],[],[],[6:8]};
run{4}     = {[],[],[],[],[],[],[],[],[],[],[9:11]};
%% ------------------------- Subj info ----edit this for new subjs---------
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
subj_name  = {'pd01','pd02','s01','s02','s03','s04','s05','s06','s07','s08','s09'};
numSess    = [2,2,2,2,2,2,2,2,2,2,4];
anatNum{1} = {[],[40],[],[],[27:31],[25:29],[29:33],[31:35],[26],[35],[]};
anatNum{2} = {[],[],[29,32],[20:24],[],[],[],[],[],[8],[43]};
anatNum{3} = {[],[],[],[],[],[],[],[],[],[],[]};
anatNum{4} = {[],[],[],[],[],[],[],[],[],[],[]};
loc_AC     = {[-78  -122 -126],...
              [-80 -128 -133],...
              [-80 -122 -123],...
              [-82 -116 -132],...
              [-81 -120 -135],...
              [-84 -112 -122],...
              [-82 -122 -133],...
              [-81 -138 -137],...
              [-77 -116 -140],...
              [-77 -124 -128],...
              [-89 -130 -125]};
          
DicomName{1}  = {'2019_03_01_PD01.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_12_pd02_sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_26_PP1_S01.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_15_S02_sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_22_S03_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_24_S04_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_30_S05_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_05_02_S06_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_11_27_S07_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_11_29_S08_Sess1.MR.Diedrichsen_PassivePatterns',...
                 '2019_12_03_S09_Sess1.MR.Diedrichsen_PassivePatterns'};
DicomName{2}  = {'2019_04_22_PD01_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_13_pd02_sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_03_28_PP1_S01_SESS2.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_16_S02_sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_24_S03_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_04_30_S04_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_10_31_S05.MR.Diedrichsen_PassivePatterns',...
                 '2019_11_13_S06_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_11_28_S07_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_12_05_S08_Sess2.MR.Diedrichsen_PassivePatterns',...
                 '2019_12_03_S09_Sess1.MR.Diedrichsen_PassivePatterns'};
% Most participants had two separate sessions, but some participants needed
% to readjust in the scanner during these sessions. Therefore, before
% adjustment, we  we acquired fieldmaps, and so it's easiest to treat these
% as separate sessions (when applying fieldmap correction). It's still two
% total sessions, but we split the runs within these sessions into
% "sub-sessions".
DicomName{3}  = {'','','','','','','','','','','2019_12_04_S09_Sess3.MR.Diedrichsen_PassivePatterns'};
DicomName{4}  = {'','','','','','','','','','','2019_12_04_S09_Sess3.MR.Diedrichsen_PassivePatterns'};

NiiRawName{1} = {'2019_03_01_PD01',...
                 '2019_03_12_pd02_sess1',...
                 '2019_03_26_PP1_S01',...
                 '2019_04_15_S02_sess1',...
                 '2019_04_22_S03_Sess1',...
                 '2019_04_24_S04_Sess1',...
                 '2019_04_30_S05_Sess1',...
                 '2019_05_02_S06_Sess1',...
                 '2019_11_27_S07_Sess1',...
                 '2019_11_29_S08_Sess1',...
                 '2019_12_03_S09_Sess1'};
NiiRawName{2} = {'2019_04_22_PD01_Sess2',...
                 '2019_03_13_pd02_sess2',...
                 '2019_03_28_PP1_S01_SESS2',...
                 '2019_04_16_S02_sess2',...
                 '2019_04_24_S03_Sess2',...
                 '2019_04_30_S04_Sess2',...
                 '2019_10_31_S05',...
                 '2019_11_13_S06_Sess2',...
                 '2019_11_28_S07_Sess2',...
                 '2019_12_05_S08_Sess2',...
                 '2019_12_03_S09_Sess1'};
NiiRawName{3} = {'','','','','','','','','','','2019_12_04_S09_Sess3'};
NiiRawName{4} = {'','','','','','','','','','','2019_12_04_S09_Sess3'};
% series numbers for epis         
fscanNum{1}   = {[21,24,27,30,33,36],...
                 [24,27,30,33,36],...
                 [14,17,20,23,26],...
                 [13,15,17,19,21,23],...
                 [12,14,16,18,20],...
                 [12,18,20,22,24],...
                 [12,14,16,18,20,22],...
                 [14,16,18,20,22,24],...
                 [3,6,9,12,15],...
                 [9,12,18,31],...
                 [6,9]};   
fscanNum{2}   = {[12,14,16,18,20],...
                 [12,15,18,21,24,27],...
                 [12,14,16,18,20,22],...
                 [29,31,35,37,39],...
                 [12,14,16,18,20,22],...
                 [14,16,18,20,22,24],...
                 [3,5,8,15,21],...
                 [10,13,25,19,22],...
                 [3,6,9,12,15,18],...
                 [11,14,17,20,23,26,29],...
                 [23,26,29]};  
fscanNum{3}   = {[],[],[],[],[],[],[],[],[],[],[3,6,9]};
fscanNum{4}   = {[],[],[],[],[],[],[],[],[],[],[15,18,21]};
% series numbers for fieldmaps
fieldNum{1}   = {[],...
                 [43,44],...
                 [31:34],...
                 [26:29],...
                 [23:26],...
                 [34:37],...
                 [25:28],...
                 [27:30],...
                 [22,23],...
                 [28:29],...
                 [19,20]};                                               
fieldNum{2}   = {[23:26],...
                 [29,30],...
                 [25:28],...
                 [42:45],...
                 [25:28],...
                 [27:30],...
                 [25,26],...
                 [32,33],...
                 [25,26,31:34],...
                 [36,37],...
                 [36,37]};
fieldNum{3}   =  {[],[],[],[],[],[],[],[],[],[],[11,12]};
fieldNum{4}   =  {[],[],[],[],[],[],[],[],[],[],[28,29]};
% prefixes of preprocessed images. r-realigned only, u-fieldmap corrected   
dataPrefix    = {'r','u','u','u','u','u','u','u','u','u','u'};             

%% ------------------------- ANALYSES -------------------------------------
switch(what)
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
    case 'LIST_subjs'               
        D = dload(fullfile(baseDir,'subj_info.txt'));
        
        if nargout==0
            fprintf('\nSN\torigSN\tID\t\tAge\tGender\tHandedness\t#fMRI Sessions\tFieldmaps\tChecked ROIs?');
            fprintf('\n--\t------\t--\t\t---\t------\t----------\t--------------\t---------\t-------------\n');
            for s = unique(D.sn)'
                S = getrow(D,ismember(D.sn,s));
                fprintf('%02d\t%s\t%s\t\t%d\t%d\t%s\t\t%d\t\t%d\t\t%d',S.sn,S.origSN{1},S.ID{1},S.age,S.gender,S.handedness{1},S.fmri_sessions,S.fmap,S.checkedROIs);
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
        I = pp1_imana('LIST_subjs');
        %sn = I.sn(I.fmri_sessions>1)';
        sn = I.sn(I.fmri_sessions>1 & I.fmap==1)'; % pd01 does not have fieldmaps- their data is much noisier
        varargout = {sn};                                                                        
    case 'LIST_rois'
        fprintf('\nROI#\tName\tHemipshere\tCortical\tThalamic');
        fprintf('\n----\t----\t----------\t--------\t-------\n');
        for r = 1:length(regSide)
            fprintf('%d\t%s\t%s\t\t%d\t\t%d\n',r,regname{r-(regSide(r)-1)*length(regname)},hem{regSide(r)},cortical(r-(regSide(r)-1)*length(regname)),thalamic(r-(regSide(r)-1)*length(regname)));
        end
    case 'MISC_check_time'                                                  % Check alignment of scanner and recorded time (sanity check): enter sn
        sn = 1;
        vararginoptions(varargin,{'sn'});

        D = dload(fullfile(behavDir,sprintf('%s_%s.dat',filePrefix,subj_name{sn})));
        % plot alignment of TR time and trial onset time
        runNum = unique(D.BN);
        numRun = numel(runNum);
        for ii=1:numRun
            r = runNum(ii);
            d = getrow(D,D.BN==r);
            subplot(2,numRun,ii);
            plot(d.startTimeMeas/1000,(d.mStartTR-1)*(TR_length) + d.mStartTRTime/1000,'LineWidth',1.5,'Color','k')
            title(sprintf('run %02d',r));
            xlabel('trial start time (s)');
            ylabel('tr start time (s)');
            grid on
            axis equal
        end
        % plot difference of TR time and trial onset time
        subplot(2,numRun,[numRun+1 : numRun*2]); 
        plot(D.startTimeMeas/1000 - ((D.mStartTR-1)*TR_length + D.mStartTRTime/1000),'LineWidth',1.5,'Color','k')
        ylabel('trial onset - tr time (s)');
        xlabel('trial number');
        title('Difference of Trial onset time and TR time')
        xlim([0 numRun*62]);
        hold on
        % draw line marking each new run
        for r = 2:numRun
            drawline((r-1)*62,'dir','vert','linestyle',':');
        end
        hold off    
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
        glm = 3;
        load(fullfile(glmDir{glm},subj_name{sn},'SPM.mat'));
        spm_rwls_resstats(SPM)          
    case 'CHECK_design'
        vararginoptions(varargin,{'sn'});
        glm = 3;
        load(fullfile(glmDir{glm},subj_name{sn},'SPM.mat'));
        imagesc(SPM.xX.X);
        colormap(bone);
        varargout = {SPM};
    case 'numVoxels'
        sn  = 1:8;
        glm = 3;
        roi = [1:4,7,8];
        vararginoptions(varargin,{'sn','glm','roi','fig'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        D = [];
        d = [];
        for i = 1:size(T.sn,1)
            d.numVox = T.numVox(i,:)';
            d.roi    = T.roi(i);
            d.sn     = T.sn(i);
            D = addstruct(D,d);
            d = [];
        end
        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn));
        
        varargout = {D};
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
        scatterplot3(Y(1:end,1),Y(1:end,2),Y(1:end,3),'split',split,'label',label,'CAT',CAT);
        % link the single digit chords for clarification
        indx=[1:5 1]';
        line(Y(indx,1),Y(indx,2),Y(indx,3),'color',color{1});

         % rest crosshairs
        hold on;
        plot3(0,0,0,'+','MarkerFaceColor',[0.75, 0, 0.75],'MarkerEdgeColor',[0.75, 0, 0.75],'MarkerSize',8);
        hold off;
        axis equal;
        xlabel('pc 1');
        ylabel('pc 2');
        zlabel('pc 3');
        
        %__________________________________________________________________  
    case '0' % ------------ BEHA: behavioural force data cases. ----------- % These functions are from fdf2/3 so need editing for your paradigm
    case 'getForceTraces'
        sn = pp1_imana('getSubjs');
        fwhm = 10; % (ms) fwhm of smoothing kernel applied to force traces
        sigma = fwhm / (2 * sqrt(2*log(2)));
        T=[]; % output structure across subjects
        for s=sn
            fprintf('%s\n',subj_name{s});
            dataFileIn  = fullfile(behavDir,sprintf('%s_%s.dat',filePrefix,subj_name{s}));
            
            D    = dload(dataFileIn);
            runs = unique(D.BN);
            for r = runs'
                trials  = D.TN(D.BN==r);
                try
                    % deal with missing mov file from sn6 (subjname=s04) run 1
                    MOV = movload(fullfile(behavDir,sprintf('%s_%s_%02d.mov',filePrefix,subj_name{s},r)));
                catch
                    fprintf('skipping %s_%s_%02d.mov - cannot load\n',filePrefix,subj_name{s},r);
                    continue
                end
                for ii = trials' 
                    d = getrow(D,D.TN==ii & D.BN==r);
                    mov = MOV{1,ii};
                    % 1. Extract data:
                    state = mov(:,1);        % trial state (4 = stimulation, 6 = response)
                    Force = mov(:,5:9);      % 5 = thumb, 9 = little finger
                    if length(unique(state))==1
                        continue
                    end
                    didx = logical([d.d1,d.d2,d.d3,d.d4,d.d5]);
                    % state timing events  
                    Ia1 = find(state==4,1,'first'); % stimulation start
                    Ia2 = find(state==4,1,'last');  % stimulation end
                    sidx = Ia1:Ia2;
                    % smooth force trace
                    f = smooth_kernel(Force,sigma); % Smoothing with Gaussian kernel
                    f = mean(f(sidx,didx),2);
                    t.f = cut(f,0,1,800,'padding','nan')';
                    % make output structure
                    t.trial    = ii;
                    t.run      = r;
                    t.chordNum = d.chordNum;
                    t.chord    = didx;
                    t.sn       = s;
                    if ismember(r,run{1}{s}) && ~isempty(t)
                        t.sess = 1;
                    elseif ismember(r,run{2}{s}) && ~isempty(t)
                        t.sess = 2;
                    end
                    T=addstruct(T,t);
                end % trials
            end % runs
        end % subjs
        dataFileOut = fullfile(behavDir,sprintf('%s_forceTraces.mat',filePrefix));
        save(dataFileOut,'-struct','T');
        varargout = {T};
        
    case 'BEHA_simpleAna'
        % simple analysis of force data:
        % - load .dat files
        % - use .dat files to calculate stimulation force
        % - use .dat files to calculate behavioural performance in each run
        % - save performance structure per participant
        sn = pp1_imana('getSubjs');
        T=[]; % output structure across subjects
        for s = sn
            fprintf('%s\n',subj_name{s});
            dataFileIn  = fullfile(behavDir,sprintf('%s_%s.dat',filePrefix,subj_name{s}));
            
            D    = dload(dataFileIn);
            runs = unique(D.BN);
            for r = runs'
                d = getrow(D,D.BN==r); % trials from this run
                % indexing fields:
                numT = numel(d.TN);
                t.sn          = ones(numT,1).*s;
                t.subjName    = repmat({subj_name{s}},numT,1);
                t.run         = d.BN;
                t.trialNum    = d.TN;
                t.chordNum    = d.chordNum;
                t.chordStim   = [d.d1 d.d2 d.d3 d.d4 d.d5];
                t.chordScreen = [d.c1 d.c2 d.c3 d.c4 d.c5];
                t.numDigits   = sum(t.chordStim,2);
                
                % analyze forces:
                t.forceN         = [d.peakF_d1, d.peakF_d2, d.peakF_d3, d.peakF_d4, d.peakF_d5]; % force per finger
                t.forceN_stim    = sum(t.forceN.*t.chordStim,2) ./ t.numDigits; % avg. force across stimulated fingers
                t.forceN_notStim = sum(t.forceN.*~t.chordStim,2) ./ (5-t.numDigits); % avg. force across non-stimulated fingers
                t.forcePer         = t.forceN./d.targetForceStim; % force expressed as % of target stimulation force
                t.forcePer_stim    = t.forceN_stim ./ d.targetForceStim; % "
                t.forcePer_notStim = t.forceN_notStim ./ d.targetForceStim; % "
                
                % analyze behaviour:
                t.trialType = d.falseResponse + 1; % 1=match, 2=mis-match
                t.correct   = ~d.isError;
                t.rxntime   = d.RT;
                
                % calculate points based on behaviour:
                t.points = zeros(size(t.correct));
                idx1 = d.RT==0 & t.trialType==1; % correct rejection (not mismatch, didn't press)
                idx2 = d.RT>0 & t.trialType==1; % false alarms (not mismatch, but pressed)
                idx3 = d.RT>0 & t.trialType==2; % hits (was mismatch, pressed)
                idx4 = d.RT==0 & t.trialType==2; % misses (was mismatch, didn't press)
                t.points(idx1) = 1; % correct rejection
                t.points(idx2) = -1; % false alarms
                t.points(idx3) = 10; % hits
                t.points(idx4) = -10; % misses
                
                % categorize behavioural performance into response types:
                v = zeros(size(t.points));
                t.resp_CR   = v;
                t.resp_FA   = v;
                t.resp_hit  = v;
                t.resp_miss = v;
                t.resp_CR(idx1)   = 1;
                t.resp_FA(idx2)   = 1;
                t.resp_hit(idx3)  = 1;
                t.resp_miss(idx4) = 1;
                t.resp_type = v;
                t.resp_type(idx1) = 1;
                t.resp_type(idx2) = 2;
                t.resp_type(idx3) = 3;
                t.resp_type(idx4) = 4;
                
                
                if ismember(r,run{1}{s})
                    t.sess = ones(numT,1).*1;
                elseif ismember(r,run{2}{s})
                    t.sess = ones(numT,1).*2;  
                end
                
                
                T=addstruct(T,t);
            end    
        end
        dataFileOut = fullfile(behavDir,sprintf('%s_simpleAna.mat',filePrefix));
        save(dataFileOut,'-struct','T');
        varargout = {T};
    case 'BEHA_simple_errorRates'
        % calculate error rates
        T = load(fullfile(behavDir,sprintf('%s_simpleAna.mat',filePrefix)));
        T.numTrials = ones(size(T.sn));
        T.mismatch = T.trialType==2;
        D=tapply(T,{'sn','trialType'},{'resp_CR','sum'},{'resp_FA','sum'},...
            {'resp_hit','sum'},{'resp_miss','sum'},{'numTrials','sum'});
        % calc error rates:
        D.prop_FA = D.resp_FA./D.numTrials;
        D.prop_CR = D.resp_CR./D.numTrials;
        D.prop_hit = D.resp_hit./D.numTrials;
        D.prop_miss = D.resp_miss./D.numTrials;
        
        % display error rates:
        prop_FA = D.prop_FA(D.trialType==1); % not mismatch
        prop_CR = D.prop_CR(D.trialType==1); % not mismatch
        prop_hit = D.prop_hit(D.trialType==2); % is mismatch
        prop_miss= D.prop_miss(D.trialType==2); % is mismatch
        prop_error = (D.resp_miss(D.trialType==2) + D.resp_FA(D.trialType==1)) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        
        tt_FA = D.resp_FA(D.trialType==1) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        tt_CR = D.resp_CR(D.trialType==1) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        tt_hit = D.resp_hit(D.trialType==2) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        tt_miss = D.resp_miss(D.trialType==2) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        
        tt_thumb = (D.resp_FA(D.trialType==1) + D.resp_hit(D.trialType==2)) ./ (D.numTrials(D.trialType==1) + D.numTrials(D.trialType==2));
        
        fprintf('mean FALSE ALARM:       %1.2f +- %1.2f (of not mismatch trials) [%1.2f +- %1.2f of total trials]\n',mean(prop_FA)*100,stderr(prop_FA)*100, mean(tt_FA)*100,stderr(tt_FA)*100);
        fprintf('mean CORRECT REJECTION: %1.2f +- %1.2f (of not mismatch trials) [%1.2f +- %1.2f of total trials]\n',mean(prop_CR)*100,stderr(prop_FA)*100, mean(tt_CR)*100,stderr(tt_CR)*100);
        fprintf('mean HIT RATE:          %1.2f +- %1.2f (of     mismatch trials) [%1.2f +- %1.2f of total trials]\n',mean(prop_hit)*100,stderr(prop_hit)*100, mean(tt_hit)*100,stderr(tt_hit)*100);
        fprintf('mean MISS RATE:         %1.2f +- %1.2f (of     mismatch trials) [%1.2f +- %1.2f of total trials]\n',mean(prop_miss)*100,stderr(prop_miss)*100, mean(tt_miss)*100,stderr(tt_miss)*100);
        
        fprintf('mean ERROR rate (FA & misses): %1.2f +- %1.2f [of total trials]\n',mean(prop_error)*100,stderr(prop_error)*100);
        
        fprintf('mean THUMB PRESS rate: %1.2f +- %1.2f [of total trials]\n',mean(tt_thumb)*100,stderr(tt_thumb)*100);
        
        varargout = {D};
    case 'BEHA_check_thumbPress'
        % counts how many instances we have a chord that is NOT followed by
        % thumb press in each run per participant:
        
        T = load(fullfile(behavDir,sprintf('%s_simpleAna.mat',filePrefix)));
        T.numTrials = ones(size(T.sn));
        T = tapply(T,{'sn','chordNum','run'},{'resp_miss','sum'},{'resp_CR','sum'});
        T.noThumb = T.resp_miss + T.resp_CR;
%         T.yesThumb = T.resp_hit + T.resp_FA;
        x=pivottable(T.chordNum,[T.sn],T.noThumb==0,'sum');
        varargout = {x,T};
    case 'BEHA_count_thumbPress'
        % counts how many instances we have a chord that is NOT followed by
        % thumb press in each run per participant:
        
        T = load(fullfile(behavDir,sprintf('%s_simpleAna.mat',filePrefix)));
        T.numTrials = ones(size(T.sn));
%        T = tapply(T,{'sn','chordNum','run'},{'resp_miss','sum'},{'resp_CR','sum'});
%        T.noThumb = T.resp_miss + T.resp_CR;
        T.yesThumb = T.resp_hit + T.resp_FA;
        x=pivottable(T.chordNum,T.sn,T.yesThumb,'sum');
        varargout = {x,T};
    
        
    case 'BEHA_simple_forceAna'
        % two analyses:
        % 1. calculate avg stimulated force (across fingers, across trials)
        % 2. finger x #fingers anova of stimulated forces:
        T = load(fullfile(behavDir,sprintf('%s_simpleAna.mat',filePrefix)));
        
        % loop through trials, and arrange forces into anova-friendly
        % structure:
        D=[];
        v=ones(5,1);
        for ii=1:size(T.sn,1)
            d.digit      = [1:5]';
            d.forceN     = T.forceN(ii,:)';
            d.forcePer   = T.forcePer(ii,:)';
            d.stimulated = T.chordStim(ii,:)';
            d.numDigits  = v.*T.numDigits(ii);
            % indexing fields:
            d.sn = v.*T.sn(ii);
            d.run = v.*T.run(ii);
            d.chordNum = v.*T.chordNum(ii);
            d.trialNum = v.*T.trialNum(ii);
            D=addstruct(D,d);
        end
        
        % 1. avg. stimulated force (across fingers, across trials):
        t = tapply(D,{'sn'},{'forceN','mean'},{'forcePer','mean'},'subset',D.stimulated==1);
        fprintf('mean force (across fingers, across combos, across subjs): %1.2fN +- %1.2fN\n',mean(t.forceN),stderr(t.forceN));
        fprintf('mean percent force: %1.2fN +- %1.2fN\n\n',mean(t.forcePer),stderr(t.forcePer));
        % 2. anova
        anovaMixed(D.forcePer,D.sn,'within',[D.digit D.numDigits],{'digit','numD'},'subset',D.stimulated==1);
%         for s=unique(D.sn)'
%             d=getrow(D,D.sn==s & D.stimulated==1);
%             MANOVA2([d.digit d.numDigits],d.forceN);
%         end

%         % plot:
%         d = tapply(D,{'sn','digit','numDigits'},{'forceN','mean'},{'forcePer','mean'},'subset',D.stimulated==1);
%         sty = style.custom(plt.helper.get_shades(5,'jet'));
%         plt.line(d.numDigits,d.forceN,'split',d.digit,'style',sty);
        
        varargout = {D};
    case 'corr_forceVsBold'
        % correlate summed finger forces with PSC per participant per
        % region:
        B = pp1_encoding('getDists','conds',1:31,'sn',2:11,'glm',4,'roi',1:6);
        F = pp1_imana('BEHA_simple_forceAna');
        F = tapply(F,{'sn','chordNum','digit'},{'forceN','mean'},'subset',F.stimulated==1);
        F = tapply(F,{'sn','chordNum'},{'forceN','sum'});
        D =[];
        for rr=1:6
            f = pivottable(F.sn,F.chordNum,F.forceN,'mean');
            b = B.pscCond(B.roi==rr,:);
            r = corr(f',b');
            d.sn   = [2:11]';
            d.corr = diag(r);
            d.roi  = ones(10,1).*rr;
            D=addstruct(D,d);
        end
        varargout = {D};
    case 'BEHA_forceMatrix'
        T = load(fullfile(behavDir,sprintf('%s_simpleAna.mat',filePrefix)));
        T = tapply(T,{'sn','run','chordNum'},{'forceN','mean'});
        varargout = {T};
        
    case 'getForceData'
        sn = pp1_imana('getSubjs');
        vararginoptions(varargin,{'sn'});
        % get data
        D = [];
        for s = sn
            load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s}))); % T
            D = addstruct(D,T);
        end
        varargout = {D};
    case 'BEHA_analyzeTrials'
        % Process force traces for fmri sessions.
        % - loop through subjects
        % - per subject per block, loop through trials and harvest force traces
        % - save output structure per subject separately
        sn = pp1_imana('getSubjs');
        vararginoptions(varargin,{'sn'});
        for s = sn
            fprintf('%s\n',subj_name{s});
            dataFileIn  = fullfile(behavDir,sprintf('%s_%s.dat',filePrefix,subj_name{s}));
            dataFileOut = fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s}));
            T           = []; % output structure for subject
            D           = dload(dataFileIn);
            runs        = unique(D.BN);
            for r = runs'
                trials  = D.TN(D.BN==r);
                try
                    % deal with missing mov file from sn6 run 1
                    MOV = movload(fullfile(behavDir,sprintf('%s_%s_%02d.mov',filePrefix,subj_name{s},r)));
                catch
                    fprintf('skipping %s_%s_%02d.mov - cannot load\n',filePrefix,subj_name{s},r);
                    continue
                end
                for ii = trials' 
                    d = getrow(D,D.TN==ii & D.BN==r);
                    t = pp1_fmri_trial(MOV{1,ii},d,0); % don't display figure per trial
                    if ismember(r,run{1}{s}) && ~isempty(t)
                        t.sess = ones(5,1).*1;
                        t.sn   = ones(5,1).*s;
                        T = addstruct(T,t);
                    elseif ismember(r,run{2}{s}) && ~isempty(t)
                        t.sess = ones(5,1).*2;  
                        t.sn   = ones(5,1).*s;
                        T = addstruct(T,t);
                    end
                end
            end
            save(dataFileOut,'-struct','T');
        end
    case 'update'
        sn = pp1_imana('getSubjs');
        vararginoptions(varargin,{'sn'});
        for s = [3:11]
            fprintf('%s\n',subj_name{s});
            dataFileIn  = fullfile(behavDir,sprintf('%s_%s.dat',filePrefix,subj_name{s}));
            dataFileOut = fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s}));
            T           = []; % output structure for subject
            D           = dload(dataFileIn);
            runs        = unique(D.BN);
            for r = runs'
                trials  = D.TN(D.BN==r);
                try
                    % deal with missing mov file from s04 run 1
                    MOV = movload(fullfile(behavDir,sprintf('%s_%s_%02d.mov',filePrefix,subj_name{s},r)));
                catch
                    fprintf('skipping %s_%s_%02d.mov - cannot load\n',filePrefix,subj_name{s},r);
                    continue
                end
                for ii = trials' 
                    d = getrow(D,D.TN==ii & D.BN==r);
                    t = pp1_fmri_trial(MOV{1,ii},d,0); % don't display figure per trial
                    if ismember(r,run{1}{s}) && ~isempty(t)
                        t.sess = ones(5,1).*1;
                        t.sn   = ones(5,1).*s;
                        T = addstruct(T,t);
                    elseif ismember(r,run{2}{s}) && ~isempty(t)
                        t.sess = ones(5,1).*2;  
                        t.sn   = ones(5,1).*s;
                        T = addstruct(T,t);
                    end
                end
            end
            
            D = load(dataFileOut);

            d = pivottable(D.numDigits,D.finger,D.peakF_raw,'nanmean','subset',D.stimulated==1);
            f = pivottable(T.numDigits,T.finger,T.peakF_raw,'nanmean','subset',T.stimulated==1);

            t = d-f;
            for ii=1:5
                for jj=1:5
                    idx = T.numDigits==ii & T.finger==jj & T.stimulated==1;
                    T.peakF_raw(idx) = T.peakF_raw(idx) + t(ii,jj);
                end
            end

            T=rmfield(T,{'peakF_filt'});

            
            
            save(dataFileOut,'-struct','T');
        end
    case 'BEHA_loadForce'
        % loads force traces
        sn = pp1_imana('getSubjs');
        % get data
        D = [];
        for s = sn
            T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s})));
            D = addstruct(D,T);
        end
        % avg. across trials and runs:
        % data for when finger is stimulated, when it is not stimulated,
        % for each appropriate number of fingers in chord: 45 values per
        % subj
        D = tapply(D,{'sn','numDigits','finger','stimulated'},...
            {'peakF_raw','nanmean'},{'peakF_filt','nanmean'},{'time_stimOnset','nanmean'},{'forceStim','mean'});
        varargout = {D};
    case 'BEHA_getErrorRate'
        % case to load and return error rate data per subject:
        % loads force traces
        sn = pp1_imana('getSubjs');
        % get data
        D = [];
        for s=sn
            T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s}))); % T
            % calculate percent incorrect trials
            T.numTrials = ones(size(T.sn));
            d = tapply(T,{'sn','numDigits','chordNum','sess'},{'isError','sum'},{'numTrials','sum'}); % T has 5 rows per trial (one per finger)
            d.isError   = d.isError./5; % correct for number of rows per trial (each trial has 5 entries- one per finger. so counts are 5x bigger)
            d.numTrials = d.numTrials./5;
            d.perErr    = d.isError./d.numTrials;
            D=addstruct(D,d);
        end
        varargout = {D};
    case 'BEHA_getDprime'
        % case to calculate Dprime scores per participant:
        sn = pp1_imana('getSubjs');
        % get data
        D = [];
        P = [];
        for s=sn
            T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s}))); % T
            % calculate percent incorrect trials
            T.numTrials = ones(size(T.sn));
            p = getrow(T,T.falseResp==0);
            p = tapply(p,{'sn','numDigits','chordNum','sess'},{'isError','sum'},{'numTrials','sum'}); % T has 5 rows per trial (one per finger)
            p.isError   = p.isError./5; % correct for number of rows per trial
            p.numTrials = p.numTrials./5;
            p.perErr    = p.isError./p.numTrials;
            falseAlarm  = sum(p.isError)/sum(p.numTrials);
            % calculate percent correct for misleading trials (mismatch
            % chord)
            f = getrow(T,T.falseResp==1);
            f = tapply(f,{'run','trial','sn','numDigits','chordNum','sess','falseResp'},{'isError','sum'});
            f.isError = f.isError./5;
            hitRate   = sum(f.isError==0)/size(f.sn,1);
            if hitRate==1 % if perfect hit rate, remove small constant to ensure norminv returns a real value (assumes no probability associated with 100% performance)
                hitRate = hitRate - (1e-16);
            end
            % calculate D prime
            d.sn     = s;
            d.dprime = norminv(hitRate) - norminv(falseAlarm);
            d.hitRate    = hitRate;
            d.falseAlarm = falseAlarm;
            D=addstruct(D,d);
        end
        varargout = {D};
    case 'PLOT_performance'
        % wrapper to plot force on stimulated vs. not stimulated trials,
        % and task performance (behaviour):

        % NOTE: all cases load data for all subjects EXCEPT sn 1 (first
        % pilot subj, who lacks fieldmaps for one session so data is noisy)
        
        % plotting styles:
        sty1 = style.custom({'black','lightgray'});
        
        % 1. plot stimlation force info:
        subplot(1,3,1);
        F = pp1_imana('BEHA_loadForce');
        F = tapply(F,{'sn','numDigits','stimulated'},{'peakF_raw','nanmean'},{'forceStim','mean'}); % avg. across fingers for simpler plot
        plt.line(F.numDigits,(F.peakF_raw./F.forceStim).*100,'split',~F.stimulated,'errorfcn','stderr','style',sty1);
        ylim([0 125]);
        drawline(100,'dir','horz','linestyle',':');
        xlabel('# digits in chord');
        ylabel('% of target force');
        title('stimulation accuracy');
        
        % 2. plot error rate per numdigits (across match and mismatch trials):
        subplot(1,3,2);
        P = pp1_imana('BEHA_getErrorRate');
        P = tapply(P,{'sn','numDigits'},{'isError','sum'},{'numTrials','sum'}); % integrate error rate across chords with same number of digits
        P.perError = (P.isError./P.numTrials).*100; 
        plt.bar(P.numDigits,P.perError,'style',sty1);
        ylim([0 100]);
        xlabel('# digits in chord');
        ylabel(sprintf('%% correct\n(summed over chords)'));
        title('error rate');
        
        % 3. plot error rate for mismatch trials only:
        subplot(1,3,3);
        D = pp1_imana('BEHA_getDprime');
        T.sn = D.sn;
        T.perError = [D.falseAlarm; 1-D.hitRate].*100; % how many repsonded yes when answer was no & responded no when answer was yes
        T.mismatch = [zeros(10,1);ones(10,1)];
        plt.bar(T.mismatch,100-T.perError);
        ylim([0 100]);
        drawline(50,'dir','horz','linestyle',':');
        ylabel('% correct');
        set(gca,'xticklabel',{'match','mismatch'});
        title('judgement accuracy');
        
        % 4. dprime:
%         subplot(2,2,4);
%         plt.bar([],D.dprime);
%         ylabel('d-prime');
%         title('sensitivity');
%         ylims = ylim;
%         ylim([0 ylims(2)]);
    case 'SUBJ_plotStimForce'
        % plots max stimulation force per finger 
        sn = [];
        vararginoptions(varargin,{'sn'});
        % load data
        T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn})));
        D = tapply(T,{'run','numDigits','finger','stimulated'},...
            {'peakF_raw','nanmean'},{'time_stimOnset','nanmean'},{'forceStim','mean'});
        % plot
        style.use('5fingers');
        
        subplot(1,3,1);
        plt.line(D.numDigits,D.peakF_raw,'split',D.finger,'errorfcn','std','subset',D.stimulated==1);
        plt.set('xlim',[0.5 5.5],'ylim',[0 4]);
        plt.labels('# fingers stimulated','stimulated force',sprintf('stimulated\ntrials'));
        drawline(unique(D.forceStim),'dir','horz','linestyle',':','linewidth',1.5);
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
        T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T
        % calculate percent incorrect trials for mismatch trials
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
        D.hitRate    = hitRate;
        D.falseAlarm = falseAlarm;
        if hitRate==1 % if perfect hit rate, remove small constant to ensure norminv returns a real value (assumes no probability associated with 100% performance)
            hitRate = hitRate - (1e-16);
        end
        % calculate D prime
        D.sn     = sn;
        D.dprime = norminv(hitRate) - norminv(falseAlarm);
        D.perErr = 1-hitRate;
        
        varargout = {D,P};
    case 'SUBJ_getDprimePerRun'
        sn = 1;
        vararginoptions(varargin,{'sn'});
        % load data
        % arrange data into plotting structure
        T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T
        % calculate percent incorrect trials for mismatch trials
        T.numTrials = ones(size(T.sn));
        P = getrow(T,T.falseResp==0);
        P = tapply(P,{'sn','numDigits','chordNum','sess','run'},{'isError','sum'},{'numTrials','sum'}); % T has 5 rows per trial (one per finger)
        P.isError   = P.isError./5; % correct for number of rows per trial
        P.numTrials = P.numTrials./5;
        P.perErr    = P.isError./P.numTrials;
        falseAlarm  = sum(P.isError)/sum(P.numTrials);
        % calculate percent correct for misleading trials (mismatch
        % chord)
        F = getrow(T,T.falseResp==1);
        F = tapply(F,{'run','trial','sn','numDigits','chordNum','sess','run','falseResp'},{'isError','sum'});
        F.isError = F.isError./5;
        hitRate   = sum(F.isError==0)/size(F.sn,1);
        D.hitRate    = hitRate;
        D.falseAlarm = falseAlarm;
        if hitRate==1 % if perfect hit rate, remove small constant to ensure norminv returns a real value (assumes no probability associated with 100% performance)
            hitRate = hitRate - (1e-16);
        end
        % calculate D prime
        D.sn     = sn;
        D.dprime = norminv(hitRate) - norminv(falseAlarm);
        D.perErr = 1-hitRate;
        
        varargout = {D,P};
    case 'SUBJ_getTaskPerformance'
        sn = 7;
        vararginoptions(varargin,{'sn'});
        % load data
        % arrange data into plotting structure
        T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T
        % calculate percent incorrect trials
        T.numTrials = ones(size(T.sn));
        P = tapply(T,{'sn','numDigits','chordNum','falseResp','sess'},{'isError','sum'},{'numTrials','sum'}); % T has 5 rows per trial (one per finger)
        P.isError   = P.isError./5; % correct for number of rows per trial
        P.numTrials = P.numTrials./5;
        P.perErr    = P.isError./P.numTrials;
        
        varargout = {P,T};
    case 'SUBJ_getTaskPerformancePerRun'
        sn = 7;
        vararginoptions(varargin,{'sn'});
        % load data
        % arrange data into plotting structure
        T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T
        % calculate percent incorrect trials
        T.numTrials = ones(size(T.sn));
        P = tapply(T,{'sn','falseResp','sess','run'},{'isError','sum'},{'numTrials','sum'}); % T has 5 rows per trial (one per finger)
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
        T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{sn}))); % T

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
        sn = pp1_imana('getSubjs');
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
        sn = pp1_imana('getSubjs');
        vararginoptions(varargin,{'sn'});
        % get data
        D = [];
        for s = sn
            T = load(fullfile(behavDir,sprintf('%s_%s_ana.mat',filePrefix,subj_name{s}))); % T
            D = addstruct(D,T);
        end
        D = tapply(D,{'sn','numDigits','finger','stimulated'},...
            {'peakF_raw','nanmean'},{'time_stimOnset','nanmean'},{'forceStim','mean'});
        
        % plot
        style.use('5fingers');
        
        subplot(1,3,1);
        plt.line(D.numDigits,(D.peakF_raw./D.forceStim).*100,'split',D.finger,'errorfcn','stderr','subset',D.stimulated==1);
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
        % plots behavioural performance during scanning session
        sn = pp1_imana('getSubjs');
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
%         subplot(1,4,[1:3]); % plot percent correct
%         sty = style.custom({[0 0 0] [0.5 0 0] [0.9 0 0] [1 0.6 0] [0.4 0.4 0.4]});
%         plt.bar(P.numDigits,1-P.perErr,'split',P.numDigits,'style',sty);
%         ylabel('% correct (all trial types)');
%         xlabel('number of fingers stimulated');
%         drawline(0.5,'dir','horz','linestyle',':');
%         ylim([0.25 1]);
        
        subplot(1,2,1); % plot % correct for mismatch trials
        sty = style.custom({[0.6 0 0.6]});
        plt.bar([],D.hitRate.*100,'style',sty,'plotall',1);
        ylabel('% trials');
        title(sprintf('correct\nrejections'));
        set(gca,'xtick',[]);
        drawline(50,'dir','horz','linestyle',':');
        ylim([0 100]);
        
        subplot(1,2,2); % plot % false rejections on non-mismatch trials
        sty = style.custom({'red'});
        plt.bar([],D.falseAlarm.*100,'style',sty,'plotall',1);
        ylabel('% trials');
        title(sprintf('false\nalarms'));
        set(gca,'xtick',[]);
        drawline(50,'dir','horz','linestyle',':');
        ylim([0 100]);
        
% 
%         subplot(1,5,5); % plot d prime per subject (sort of meaningless without comparison group but whatever)
%         sty = style.custom({'blue'});
%         plt.dot([],D.dprime,'style',sty);
%         ylabel('d''');
%         title(sprintf('signal\ndetection'));
%         set(gca,'xtick',[]);
%         drawline(0,'dir','horz','linestyle',':');
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
    case 'WRAPPER1_dicomImport'                                             
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
                    fprintf('importing FUNCTIONAL runs- %s sess %d\n',subj_name{sn},sess);
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
    case 'WRAPPER2_anatomical'
        sn = [];
        vararginoptions(varargin,{'sn'});
        if numel(sn)>1; error('one subj at a time.'); end
        pp1_imana('PREP_resliceLPI','sn',sn);
        pp1_imana('PREP_centreAC','sn',sn);
        pp1_imana('PREP_segmentation','sn',sn);
    case 'WRAPPER3_preprocessFunc'
        sn = [];
        vararginoptions(varargin,{'sn'});
        if numel(sn)>1; error('one subj at a time.'); end
        pp1_imana('PREP_moveRaw4d','sn',sn);
        pp1_imana('PREP_realignEst','sn',sn);
        pp1_imana('PREP_calcVDM','sn',sn);
        pp1_imana('PREP_applyVDM','sn',sn);
        pp1_imana('PREP_moveData','sn',sn);
        pp1_imana('PREP_coreg','sn',sn);
    case 'WRAPPER4_coregAlignFunc'
        sn = [];
        vararginoptions(varargin,{'sn'});
        if numel(sn)>1; error('one subj at a time.'); end
        pp1_imana('PREP_coreg','sn',sn);
        pp1_imana('PREP_makeSameAlign','sn',sn);
        pp1_imana('PREP_makeMaskImage','sn',sn);
    
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
        % https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#Chap:FieldMap
        % https://osf.io/k6rm5/wiki/1.1_Field_map_correction/
        % https://www.fil.ion.ucl.ac.uk/spm/data/fieldmap/
        
        % Set options for batch job
        spm_dir= fileparts(which('spm'));
        %J.defaults.defaultsfile = {'/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_passivePatterns/pm_pp1.m'};
        J.defaults.defaultsval.et              = [4.08 5.1];    % [shortest, longest echotimes]                                                       
        J.defaults.defaultsval.etd             = J.defaults.defaultsval.et(2) - J.defaults.defaultsval.et(1);
        J.defaults.defaultsval.tert            = 1000/(22.904*2); % this is GRAPPA accel. factor (2), not multi-band accel. factor (3)
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
        
%         nam{1}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
%         nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c2' subj_name{sn}, '_anatomical.nii']);
%         nam{3}  = fullfile(anatomicalDir, subj_name{sn}, ['c3' subj_name{sn}, '_anatomical.nii']);
%         spm_imcalc(nam, 'rmask_noskull.nii', '(i1+i2+i3)>0.2');

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
    case 'SURF_freesurfer'   
        vararginoptions(varargin,{'sn'});
        freesurfer_reconall(freesurferDir,subj_name{sn},fullfile(anatomicalDir,subj_name{sn},[subj_name{sn} '_anatomical.nii']));
    case 'SURF_WBresample'   
        % This reslices from the individual surfaces into the the fs_lr
        % standard mesh - This replaces calls to freesurfer_registerXhem,
        % freesurfer_mapicosahedron_xhem, & caret_importfreesurfer. It
        % requires connectome wb to be installed, added to the bash_profile
        % (on terminal), and updated on the startup.m file
        vararginoptions(varargin, {'sn'});
        fprintf('reslicing %s...',subj_name{sn});
        surf_resliceFS2WB(subj_name{sn}, freesurferDir, wbDir,'resolution','32k'); 
        fprintf('done\n');  
    case 'SURF_segThalamus'
        % segment the thalamus nuclei
        I = pp1_imana('LIST_subjs');
        sn = [];
        vararginoptions(varargin,{'sn'});
        % do segmentation
        % % ***SEGMENTATION MUST BE DONE VIA TERMINAL
%         cmd = ['segmentThalamicNuclei.sh ' subj_name{sn} ' ' freesurferDir];
%         fprintf('%s\n',cmd);
%         [status,result] = system(cmd);
%         if status; error(result); end
        % convert segmentation to nii and move segmentation file to
        % anatomicals directory
        inName  = fullfile(freesurferDir,I.origSN{sn},'mri','ThalamicNuclei.v10.T1.mgz');
        outName = fullfile(anatomicalDir,I.origSN{sn},[I.origSN{sn} '_thalamicNuclei_raw.nii']);
        cmd = sprintf('mri_convert %s %s', inName, outName);
        fprintf('%s\n',cmd);
        [status,result] = system(cmd);
        if status; error(result); end
        % reslice the thalamic masks into same space as epis (this will
        % also mask them to be within bounds of epi FOV)
        mask = fullfile(imagingDir,I.origSN{sn},['rmask_gray.nii']);
        T = spm_vol(outName);
        outName = fullfile(anatomicalDir,I.origSN{sn},[I.origSN{sn} '_thalamicNuclei.nii']);
        M = spm_vol(mask);
        spmj_reslice_vol(T,M.dim,M.mat,outName);
        % finally, mask the thalamus rois by the FOV of the epis so we
        % don't get voxels without signal in the definitions
%         mask = fullfile(imagingDir,I.origSN{sn},['rmask_noskull.nii']);
%         outName2 = fullfile(anatomicalDir,I.origSN{sn},[I.origSN{sn} '_thalamicNuclei_masked.nii']);
%         cmd = sprintf('fslmaths %s -mas %s %s', outName, mask, outName2);
%         [status,result] = system(cmd);
%         if status; error(result); end
        fprintf('done.\n');   
    case 'SURF_check_alignment'        % check the alignment for the identified sensorimotor cortex for controls/patients
        regNum      = [1:6];
        % 0. Load sulcal geometries from the chording data to be
        % used as a reference for both controls and patients
        s = spm_vol(fullfile(atlasDir,'/fs_LR.164k.LR.sulc.dscalar.nii'));
        s = s.private.dat(:);
        ref(:,1) = s(1:163842);
        ref(:,2) = s(163843:end);
        for h=1:2
            % get regions of interest from atlas label file
            s           = gifti(fullfile(atlasDir,['ROI_pp1.164K.',hemLetter{h},'.label.gii']));
            roi(:,h)    = ismember(s.cdata,regNum);
        end
        
        % 1. Get list of subjects
        D = pp1_imana('LIST_subjs');
        S = [];
        for ii=1:length(D.sn)
            sn = D.origSN{ii};
            disp(sn);
            
            for h=1:2
                % load sulcal file for subject
                s           = gifti(fullfile(wbDir,sn,[sn '.' hemLetter{h} '.sulc.164k.shape.gii']));
                subjSulc    = -double(s.cdata); % signs seem flipped?
                mse         = mean((ref(roi(:,h))-subjSulc(roi(:,h))).^2);  % difference between ref sulc and subj sulc in roi

                % making/saving alldat structure
                Dh              = getrow(D,strcmp(D.origSN,sn));
                Dh.hem          = h;
                Dh.mse          = mse;                
                S               = addstruct(S,Dh);
            end
        end
        % Plot error in registration across hemispheres for patients
        meanC   = mean(S.mse);
        stdC    = std(S.mse);
        
        % plot controls as dashed lines
        lineplot(S.sn,S.mse,'split',S.hem,'style_thickline',...
                 'leg',hemName,'leglocation','northwest','markersize',6);
        hold on;
        % plot confidence bounds as solid lines    
        drawline(meanC+2*stdC,'dir','horz','color','k');
        drawline(meanC-2*stdC,'dir','horz','color','k');
       % set(gca,'FontSize',14);
        set(gca,'xtick',unique(S.sn));
        xlabel('subject #');
        ylabel('mse');
        title('deviation of regions between subject and reference surf');
        hold off
        
        varargout = {S};
%         % List bad participants/save all participants in file
%         S.bad = S.mse>=(meanC+2*stdC);
%         Sp = getrow(S,S.bad==1);
%         fprintf('\nBad participants\n');
%         fprintf('----------------\n');
%         for i=1:length(Sp.SN)
%             fprintf('SN: %s, Hem: %d\n',Sp.ID{i},Sp.hem(i));
%         end
%         dsave(fullfile(baseDir,'subj_list_check_alignment.txt'),S);
    
    case '0' % ------------ GLM: SPM GLM fitting. Expand for more info. ---
        % The GLM cases fit general linear models to subject data with 
        % SPM functionality.
        %
        % All functions can be called with ('GLM_processAll','sn',[Subj#s]).
        %
        % You can view reconstructed surfaces with Caret software.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'GLM_hrfParams'
        % case to store hrf params per participant
        [~,d_hrf]      = spm_hrf(TR_length); % default hrf params
        subj_hrfParams = {[3.6873 11.893 d_hrf(3) d_hrf(4) 0.23299],... % subject-specific hrf params (estimated from glm1)
                          [5.5018 15.717 d_hrf(3) d_hrf(4) 6.226],...
                          [4.986  15.452 d_hrf(3) d_hrf(4) 6.1327],...
                          [5.0097 16.379 d_hrf(3) d_hrf(4) 5.6836],...
                          [4.9406 13.158 d_hrf(3) d_hrf(4) 2.6733],...
                          [4.1708 12.603 d_hrf(3) d_hrf(4) 4.6819],...
                          [4.5518 14.833 d_hrf(3) d_hrf(4) 5.2753],...
                          [5.329  14.959 d_hrf(3) d_hrf(4) 5.0944],...
                          [4.1893 13.015 d_hrf(3) d_hrf(4) 1.2336],...
                          [4.9436 15.137 d_hrf(3) d_hrf(4) 0.56789],...
                          [5.8241 16.161 d_hrf(3) d_hrf(4) 4.6441]};
        varargout = {subj_hrfParams};
    case 'GLM_hrfParamsGLM4'
        % case to store hrf params per participant
        [~,d_hrf]      = spm_hrf(TR_length); % default hrf params
        subj_hrfParams = {[1.5 12 d_hrf(3) d_hrf(4) 6 0],... % subject-specific hrf params (estimated from glm1)
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0],...
                          [1.5 12 d_hrf(3) d_hrf(4) 6 0]};
        varargout = {subj_hrfParams};
    case 'WRAPPER_GLM'                                                  
        vararginoptions(varargin,{'sn','glm'});
        % You can call this case to do all the GLM estimation and contrasts.
        for s = sn
            for g = glm
                if s>1
                    pp1_imana('GLM_make','sn',s,'glm',g);
                    pp1_imana('GLM_estimate','sn',s,'glm',g);
                    if g==1
                        pp1_imana('GLM_contrastglm1','sn',s);
                    elseif g>1
                        if g==2
                            pp1_imana('GLM_contrastglm2','sn',s);
                        elseif g==3
                            pp1_imana('GLM_contrastglm3','sn',s);
                        elseif g==4
                            pp1_imana('GLM_contrastglm4','sn',s);
                        elseif g==5
                            pp1_imana('GLM_contrastglm5','sn',s);
                        end
                
                    end
                end
                pp1_imana('PSC_calcChord','sn',s,'glm',g);
                mySendmail(sprintf('done glm%d s%02d',g,s));
            end
        end
    case 'GLM_make'                                                         
        % makes the GLM file for each subject, and a corresponding aux.
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        % housekeeping
        prefix		 = dataPrefix{sn};
        T			 = [];
                      
        % Define number of regressors in glm
        switch glm 
            case 1
                % model all conditions together and use to optimize hrf fit
                [~,d_hrf]  = spm_hrf(TR_length); % default hrf params
                hrf_params = d_hrf; % default hrf params
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 1; % to optimize hrf fits
            case 2
                % model all chords, don't exclude error trials
                % define the 7 parameters of the HRF
                subj_hrfParams = pp1_imana('GLM_hrfParams');
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 31; % 31 chords
            case 3
                % model all chords (31) and a regressor for thumb movements
                % (reg number 32)
                % define the 7 parameters of the HRF
                subj_hrfParams = pp1_imana('GLM_hrfParams');
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 32; % 31 chords + thumb responses
            case 4 % USE THIS ONE!
                % model all chords (31) and a regressor for thumb movements
                % (reg number 32)
                % define the 7 parameters of the HRF
                subj_hrfParams = pp1_imana('GLM_hrfParamsGLM4');
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 32; % 31 chords + thumb responses
            case 5
                % model all chords (31), regressor for the thumb presses
                % (reg#33), and regressor for the stimulations that
                % preceeded thumb presses (reg#32)
                subj_hrfParams = pp1_imana('GLM_hrfParamsGLM4');
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 33;
        end
        
        % Load subject's .dat file (has info on each trial)
        D = dload(fullfile(behavDir,sprintf('pp1_fmri_%s.dat',subj_name{sn})));
        % calculate onset times (in seconds), correcting for removed dummy scans
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
            % determine session #
            for kk = 1:numel(run)
                if sum(ismember(run{kk}{sn},r))
                    sess = kk;
                    break
                end
            end
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
                    case {3,4}
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
                    case 5
                        S.isError = 0;
                        if c <32
                            % Model chords regressors
                            % include all trials regardless of judgement
                            idx	= find(R.chordNum==c & R.RT==0); % find indx of all trials in run of that condition and did not preceed thumb press
                            J.sess(r).cond(c).name     = sprintf('chord_%d',R.chordNum(idx(1)));  % make condition name (for user readability)
                            J.sess(r).cond(c).onset    = R.onset(idx);    
                            J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
                        elseif c==32
                            % Model thumb presses
                            idx = find(R.RT>0);
                            J.sess(r).cond(c).name     = 'thumb_move';  % make condition name (for user readability)
                            J.sess(r).cond(c).onset    = R.onset(idx) + R.RT(idx)/1000;    
                            J.sess(r).cond(c).duration = 1;
                        elseif c==33
                            % Model stimulation that preceeded thumb press
                            idx	= find(R.RT>0);
                            J.sess(r).cond(c).name     = 'bad_stim';
                            J.sess(r).cond(c).onset    = R.onset(idx);    
                            J.sess(r).cond(c).duration = R.stimTime(idx)/1000;
                            S.isError=1;
                        end
                end
                J.sess(r).cond(c).tmod = 0;
                J.sess(r).cond(c).orth = 0;
                J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                % Do some subject info for fields in SPM_info.mat.
                S.sn   = sn;
                S.run  = r;
                S.sess = sess;
                if (c==1 && glm==1)
                    S.chord 		= 0;                          
                    S.numDigits     = 0; 
                    S.targetForce   = R.targetForceStim(idx(1));
                    S.numStim       = R.numStim(idx(1));
                    S.regtype       = 'avgTask';
                elseif c<32
                    S.chord 		= R.chordNum(find(R.chordNum==c,1));                       
                    S.numDigits     = R.numDigits(find(R.chordNum==c,1)); 
                    S.targetForce   = R.targetForceStim(find(R.chordNum==c,1));
                    S.numStim       = R.numStim(find(R.chordNum==c,1));
                    S.regtype       = sprintf('chord%02d',R.chordNum(idx(1)));
                elseif c==32 % thumb press
                    S.chord 		= 32;                          
                    S.numDigits     = 0; 
                    S.targetForce   = 0;
                    S.numStim       = 0;
                    S.regtype       = 'respDg1';
                elseif c==33 % stimulations that preceeded thumb press
                    S.chord 		= 33;                          
                    S.numDigits     = 0; 
                    S.targetForce   = 0;
                    S.numStim       = 0;
                    S.regtype       = 'badStim';
                end  
                S.tt = c;
                T	 = addstruct(T,S);
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
    case 'GLM_estimate'                                                     
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
    case 'GLM_contrastglm2'                                                 
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
    case 'GLM_contrastglm3'                                                 
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
    case 'GLM_contrastglm4'                                                 
        % Make t-stat contrasts for specified GLM estimates.
        % enter sn, glm #
        % models each chord and also error trials
        vararginoptions(varargin,{'sn'});
        glm = 4;
        % Go to subject's directory
        I = pp1_imana('LIST_subjs');
        subjName = I.origSN{I.sn==sn};
        fprintf('%s : ',subjName);
        % Go to subject's directory
        subjDir = fullfile(glmDir{glm},subjName);
        cd(subjDir);
        load SPM;
        SPM = rmfield(SPM,'xCon');
        T   = load('SPM_info.mat');
        % erase any current con_, spmT_, and psc_ files to avoid confusion
        delete(fullfile(subjDir,'con_*.nii'));
        delete(fullfile(subjDir,'spmT_*.nii'));
        delete(fullfile(subjDir,'psc_*.nii'));
        
        %_____t contrast for chords & thumb movement (1:32)
        for dd = 1:32
            con               = zeros(1,size(SPM.xX.X,2));
            con(:,T.chord==dd) = 1;
            con               = con/sum(con);
            if dd<32 % chords
                SPM.xCon(dd) = spm_FcUtil('Set',sprintf('chord_%d',dd), 'T', 'c',con',SPM.xX.xKXs);
            elseif dd==32 % thumb movement
                SPM.xCon(dd) = spm_FcUtil('Set','thumb_response', 'T', 'c',con',SPM.xX.xKXs);
            end
        end

        %_____t contrast overall chords vs. rest (33)
        con                    = zeros(1,size(SPM.xX.X,2));
        con(:,T.chord>0 & T.chord<32) = 1;
        con                    = con/sum(con);
        SPM.xCon(33)           = spm_FcUtil('Set','overall', 'T', 'c',con',SPM.xX.xKXs);
        
        %_____t contrast for chords split by # digits (34:38)
        for dd=1:5
            con                    = zeros(1,size(SPM.xX.X,2));
            con(:,T.numDigits==dd) = 1;
            con                    = con/sum(con);
            SPM.xCon(end+1)        = spm_FcUtil('Set',sprintf('numDigits_%d',dd), 'T', 'c',con',SPM.xX.xKXs);
        end
        
        %____do the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save('SPM.mat','SPM');
    case 'GLM_contrastglm5'                                                 
        % Make t-stat contrasts for specified GLM estimates.
        % enter sn, glm #
        % models each chord and also error trials
        vararginoptions(varargin,{'sn'});
        glm = 4;
        % Go to subject's directory
        I = pp1_imana('LIST_subjs');
        subjName = I.origSN{I.sn==sn};
        fprintf('%s : ',subjName);
        % Go to subject's directory
        subjDir = fullfile(glmDir{glm},subjName);
        cd(subjDir);
        load SPM;
        SPM = rmfield(SPM,'xCon');
        T   = load('SPM_info.mat');
        % erase any current con_, spmT_, and psc_ files to avoid confusion
        delete(fullfile(subjDir,'con_*.nii'));
        delete(fullfile(subjDir,'spmT_*.nii'));
        delete(fullfile(subjDir,'psc_*.nii'));
        
        %_____t contrast for chords & thumb movement (1:32)
        for dd = 1:32
            con               = zeros(1,size(SPM.xX.X,2));
            con(:,T.chord==dd) = 1;
            con               = con/sum(con);
            if dd<32 % chords
                SPM.xCon(dd) = spm_FcUtil('Set',sprintf('chord_%d',dd), 'T', 'c',con',SPM.xX.xKXs);
            elseif dd==32 % thumb movement
                SPM.xCon(dd) = spm_FcUtil('Set','thumb_response', 'T', 'c',con',SPM.xX.xKXs);
            end
        end

        %_____t contrast overall chords vs. rest (33)
        con                    = zeros(1,size(SPM.xX.X,2));
        con(:,T.chord>0 & T.chord<32) = 1;
        con                    = con/sum(con);
        SPM.xCon(33)           = spm_FcUtil('Set','overall', 'T', 'c',con',SPM.xX.xKXs);
        
        %_____t contrast for chords split by # digits (34:38)
        for dd=1:5
            con                    = zeros(1,size(SPM.xX.X,2));
            con(:,T.numDigits==dd) = 1;
            con                    = con/sum(con);
            SPM.xCon(end+1)        = spm_FcUtil('Set',sprintf('numDigits_%d',dd), 'T', 'c',con',SPM.xX.xKXs);
        end
        
        %____do the constrasts
        SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save('SPM.mat','SPM');
    case 'PSC_calcChord'                                                    
        % calculate psc for all digits vs. rest - based on betas from glm 1 
        % we want psc = [max predicted / baseline activity] * 100
        % To do that appropraitely, here we do:
        % 100* median(max(regressor height of each condition)) * contrast
        % of interest / mean(intercept beta images = baselines)
        I = pp1_imana('LIST_subjs');
        sn  = [];
        glm = [];
        vararginoptions(varargin,{'sn','glm'});
        if glm==2
            numImgs = 31;
        elseif glm==3
            numImgs = 32;
        elseif glm==4
            numImgs = 38; % 31 chords, thumb response (32), overall avg (33), & per num digits (34:38)
        end
        % assumes first 31 contrasts are for the chords
        for s=sn
            subjName = I.origSN{I.sn==s};
            cd(fullfile(glmDir{glm}, subjName));
            load SPM;
            T = load('SPM_info.mat');
            X = (SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
            h = median(max(X));               % Height of response is defined as median of max regressor height (across conds and runs) for this subject
            P = {};                           % Filenames of input images
            numB = length(SPM.xX.iB);         % Partitions - runs
            for p = SPM.xX.iB
                P{end+1} = sprintf('beta_%4.4d.nii',p);       % get the intercepts (for each run) and use them to calculate the baseline (mean images) * max height of design matrix regressor
            end
            for con = 1:numImgs   
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
            fprintf('%s: %3.3f\n',subjName,h);
        end
    case 'PSC_calcMF'
        % make multi-finger psc contrast using already made psc contrasts:
        sn  = [];
        glm = [];
        vararginoptions(varargin,{'sn','glm'});
        % Use 'con' option to define different contrasts.
        I = pp1_imana('LIST_subjs');
        for s=sn
            subjName = I.origSN{I.sn==s};
            subjGLM  = fullfile(glmDir{glm}, subjName);
            % Load subject contrasts:
            fileIn = fullfile(subjGLM, 'psc_06.nii,1'); % load file for chord 6 
            vol = spm_vol(fileIn);
            vdat = spm_read_vols(vol);
            for ii=7:31 % chord numbers
                fileIn = fullfile(subjGLM, sprintf('psc_%02d.nii,1',ii)); 
                vol = spm_vol(fileIn);
                vdat_chord = spm_read_vols(vol);
                vdat = vdat + vdat_chord;
            end
            Y.data = vdat./26; % average psc across multi-finger chords per voxel
            % prep output file
            Y.dim   = vol(1).dim;
            Y.dt    = vol(1).dt;
            Y.mat   = vol(1).mat;    
            % save output
            Y.fname   = 'psc_mf.nii';
            Y.descrip = sprintf('exp: ''pp1'' \nglm: ''FAST'' \ncontrast: multi-finger chords');
            cd(subjGLM);
            spm_write_vol(Y,Y.data);
            fprintf('Done %s %s\n',subjName,Y.fname);
        end
        
    case '0' % ------------ Cases to make surface PSC maps. ---------------
    case 'WB:mapPSC'
        % mapping PSC maps to surface
        % wrapper for cases 'vol2surf_indiv' & 'vol2surf_stats'
        res = '164k';
        glm = 4;
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        map = 'PSC';
        
        % map chords split by number of digits per participant
        pp1_imana('WB:vol2surf_subj','sn',sn,'glm',glm,'res',res,'map',map);
        
        % now do group map
        groupDir = fullfile(wbDir, ['group_' res]);
        groupFiles = cell(1,numel(sn)); % pad filename cell array (rows=hemis, cols=subjs)
        for ii=1:2 % loop through hemis
            outFile = {};
            for jj=1:numel(sn) % get input files from subjects
                subjName = I.origSN{I.sn==sn(jj)};
                groupFiles{1,jj} = fullfile(wbDir, subjName, sprintf('%s.%s.glm%d.%s.%s.func.gii',subjName,hemLetter{ii},glm,map,res));
            end
            outFile{1} = fullfile(groupDir, sprintf('%s.glm%d.%s.avg.%s.func.gii', hemLetter{ii}, glm, map, res));
            outFile{2} = fullfile(groupDir, sprintf('%s.glm%d.%s.mf_config.%s.func.gii', hemLetter{ii}, glm, map, res));
            for dd=1:5
                outFile{end+1} = fullfile(groupDir, sprintf('%s.glm%d.%s.%dd_config.%s.func.gii', hemLetter{ii}, glm, map, dd, res));
            end
            surf_groupGiftis(groupFiles, 'outcolnames',I.origSN', 'outfilenames',outFile);
        end

        % do summary stats:
        maps = {'PSC.avg','PSC.mf_config','PSC.1d_config','PSC.2d_config','PSC.3d_config','PSC.4d_config','PSC.5d_config'};
        pp1_imana('WB:vol2surf_stats','res',res,'glm',glm,'maps',maps,'summaryMap',map);
    case 'WB:mapLDC'
        % mapping LDC (average distance) maps to surface
        res = '164k';
        glm = 4;
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        map = 'LDC';
        
        % map chords split by number of digits per participant
         pp1_imana('WB:vol2surf_subj','sn',sn,'glm',glm,'res',res,'map',map);
        
        % now do group map
        groupDir = fullfile(wbDir, ['group_' res]);
%         dircheck(groupDir);
        groupFiles = cell(1,numel(sn)); % pad filename cell array (rows=hemis, cols=subjs)
        for ii=1:2 % loop through hemis
            outFile = {};
            for jj=1:numel(sn) % get input files from subjects
                subjName = I.origSN{I.sn==sn(jj)};
                groupFiles{1,jj} = fullfile(wbDir, subjName, sprintf('%s.%s.glm%d.%s.%s.func.gii',subjName,hemLetter{ii},glm,map,res));
            end
            outFile{1} = fullfile(groupDir, sprintf('%s.glm%d.%s.avg.%s.func.gii', hemLetter{ii}, glm, map, res));
            outFile{2} = fullfile(groupDir, sprintf('%s.glm%d.%s.mf_config.%s.func.gii', hemLetter{ii}, glm, map, res));
            for dd=1:4
                outFile{end+1} = fullfile(groupDir, sprintf('%s.glm%d.%s.%dd_config.%s.func.gii', hemLetter{ii}, glm, map, dd, res));
            end
            surf_groupGiftis(groupFiles, 'outcolnames',I.origSN', 'outfilenames',outFile);
        end

        % do summary stats:
        maps = {'LDC.avg','LDC.mf_config','LDC.1d_config','LDC.2d_config','LDC.3d_config','LDC.4d_config'};
        pp1_imana('WB:vol2surf_stats','res',res,'glm',glm,'maps',maps,'summaryMap',map);
        
    case 'WB:mapChords'
        % mapping chord t-maps to surface
        % wrapper for cases 'vol2surf_indiv' & 'vol2surf_stats'
        res = '164k';
        glm = 4;
        sn  = pp1_imana('getSubjs');
        map = 't';
        
        % map chords split by number of digits per participant
        pp1_imana('WB:vol2surf_subj','sn',sn,'glm',glm,'res',res,'map',map);
    
    case 'WB:vol2surf_subj'                                               
        I = pp1_imana('LIST_subjs');
        % map indiv vol contrasts (.nii) onto surface (.gifti)
        sn    = [];       
        glm   = [];              
        hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'
        h     = [1 2];
        map   = []; % 'psc'; 'avg';
        res   = []; % resolution
        vararginoptions(varargin,{'sn', 'glm', 'h', 'map', 'res'});
        
        % load topo files (for restricting voxels across central sulcus)
        if strcmp(res,'32k')
            topo{1} = gifti(sprintf('/Users/sarbuckle/DATA/Atlas_templates/FS_LR_32/fs_LR.%s.%s.flat.surf.gii',res,hemLetter{1}));
            topo{2} = gifti(sprintf('/Users/sarbuckle/DATA/Atlas_templates/FS_LR_32/fs_LR.%s.%s.flat.surf.gii',res,hemLetter{2}));
        elseif strcmp(res,'164k')
            topo{1} = gifti(sprintf('/Users/sarbuckle/DATA/Atlas_templates/FS_LR_164/fs_LR.%s.%s.flat.surf.gii',res,hemLetter{1}));
            topo{2} = gifti(sprintf('/Users/sarbuckle/DATA/Atlas_templates/FS_LR_164/fs_LR.%s.%s.flat.surf.gii',res,hemLetter{2}));
        end
        
        for s=sn
            subjName = I.origSN{I.sn==s};
            
            for ii=h
                surfDir     = fullfile(wbDir, subjName);
                white       = fullfile(surfDir, sprintf('%s.%s.white.%s.surf.gii',subjName, hemLetter{ii}, res));
                pial        = fullfile(surfDir, sprintf('%s.%s.pial.%s.surf.gii',subjName, hemLetter{ii}, res));
                C1          = gifti(white);
                C2          = gifti(pial);
                subjGLM     = fullfile(glmDir{glm}, subjName);
                fnames = {};
                switch map
                    case 'LDC' % searchlight dissimilarity maps
                        con_name = {'avg','mf_config','1d_config','2d_config','3d_config','4d_config'};% 
                        fnames{1}   = fullfile(subjGLM,sprintf('%s_glm%d_LDC_avg.nii',subjName,glm)); % avg. paired crossnobis
                        fnames{2}   = fullfile(subjGLM,sprintf('%s_glm%d_LDC_mf_config.nii',subjName,glm)); % avg. crossnobis between multi-finger chords
                        for cc=1:4
                            fnames{end+1} = fullfile(subjGLM,sprintf('%s_glm%d_LDC_%dd_config.nii',subjName,glm,cc));
                        end
                        outfile = fullfile(surfDir, sprintf('%s.%s.glm%d.LDC.%s.func.gii', subjName, hemLetter{ii}, glm, res));
                    case 'PSC' % maps percent-signal-change maps onto surfaces
                        con_name = {'avg','mf_config','1d_config','2d_config','3d_config','4d_config','5d_config'};
                        fnames{1} = fullfile(subjGLM, 'psc_33.nii,1'); % psc33==avg. psc across stimulation conditions
                        fnames{2} = fullfile(subjGLM, 'psc_mf.nii,1'); % avg psc across multi-finger stimulation conditions
                        for dd=1:5
                            fnames{end+1} = fullfile(subjGLM, sprintf('psc_%d.nii,1',33+dd)); %psc_34==numDigits1, etc.
                        end
                        outfile = fullfile(surfDir, sprintf('%s.%s.glm%d.PSC.%s.func.gii', subjName, hemLetter{ii}, glm, res));
                    case 't'   % t-values maps (univariate GLM)
                        load(fullfile(subjGLM,'SPM.mat'));
                        fnames      = cell(1,32);
                        con_name    = cell(1,32);
                        for jj=1:numel(fnames)
                            con_name{jj} = SPM.xCon(jj).name;
                            fnames{jj}   = fullfile(subjGLM, sprintf('spmT_%04d.nii', jj));
                        end
                        outfile = fullfile(surfDir, sprintf('%s.%s.glm%d.TMAPS.%s.func.gii', subjName, hemLetter{ii}, glm, res));
                end

                G = surf_vol2surf(double(C1.vertices), double(C2.vertices), fnames, 'column_names',con_name, 'anatomicalStruct',hname{ii},...
                    'exclude_thres',0.75,'faces',double(topo{ii}.faces),'ignore_zeros',1);
                
                save(G, outfile);
                fprintf('%s\n',outfile);
            end
        end
    case 'WB:vol2surf_stats'                                                
        % do stats on mapped group surface contrasts (.gifti)
        glm   = [];              
        hem   = {'L', 'R'};     % hemisphere: 1=LH 2=RH
        h     = [1 2];
        hname = {'CortexLeft', 'CortexRight'};
        maps  = {''};
        res   = [];
        summaryMap = [];
        %sm    = 0; % smoothing kernel in mm (optional)
        vararginoptions(varargin,{'sn', 'glm', 'h', 'maps','group','res','summaryMap'});

        groupDir = fullfile(wbDir, ['group_' res]);
        numMaps  = numel(maps);
        % Loop over the metric files and calculate the cSPM of each
        for ii=h
            for m = 1:numMaps
                fprintf('stats %s : %s',h,maps{m});
                groupfiles{m}   = fullfile(groupDir, sprintf('%s.glm%d.%s.%s.func.gii', hem{ii},glm,maps{m},res));
                metric          = gifti(groupfiles{m});
                cSPM            = surf_getcSPM('onesample_t', 'data',metric.cdata, 'maskthreshold',0.7);
                C.data(:,m)     = cSPM.con.con; % mean
                C.c_name{m}     = ['mean_' maps{m}];
                C.data(:,m+numMaps) = cSPM.con.Z; % t
                C.c_name{m+numMaps} = ['t_' maps{m}];
                fprintf('...done.\n');
            end
            % Save output
            O = surf_makeFuncGifti(C.data, 'columnNames',C.c_name, 'anatomicalStruct',hname{ii});
            summaryfile = fullfile(groupDir, sprintf('summary.%s.%s.glm%d.%s.func.gii', summaryMap,hem{ii},glm,res));
            save(O, summaryfile);
        end   
    
    case 'surf:profile_plot'
        
        clr = [0.2 0.14 0.53];
        linestyle = '-';
        % activation profile for stimulation
        map = [];
        location = 'middle'; % bottom, top
        vararginoptions(varargin,{'map','location'});
        saveBorder = 1;
        
        groupDir = fullfile(wbDir, ['group_164k']);
        switch map
            case 'LDC.avg'
                % dissimilarities avg. across all stimulation
                % configurations
                dataFile = fullfile(groupDir,'L.glm4.LDC.avg.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'average dissimilarity';
                clr = [0.8 0.27 0.5];
            case 'LDC.1d'
                % dissimilarities b/t single fingers
                dataFile = fullfile(groupDir,'L.glm4.LDC.1d_config.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'single finger dissimilarity';
                clr = [0.2 0.14 0.53];
            case 'LDC.2d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.LDC.2d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'finger pair dissimilarity';
                clr = [0.19048 0.80952 1];
            case 'LDC.3d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.LDC.3d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'finger triplet dissimilarity';
                clr = [0.38095 0.61905 1];
            case 'LDC.4d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.LDC.4d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'finger quad dissimilarity'; 
                clr = [0.57143 0.42857 1];
            case 'LDC.mf'
                % dissimilarities b/t single fingers
                dataFile = fullfile(groupDir,'L.glm4.LDC.mf_config.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'multi finger dissimilarity';
                clr = [0.5 0.5 0.8];
            case 'LDC.avg.ipsi'
                % signal change across all stimulation configurations from
                % ipsilateral side
                dataFile = fullfile(groupDir,'R.glm4.LDC.avg.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'average finger dissimilarity';
                clr = [0.8 0.27 0.5];
                linestyle = ':';
            case 'LDC.1d.ipsi'
                % signal change across all stimulation configurations from
                % ipsilateral side
                dataFile = fullfile(groupDir,'R.glm4.LDC.1d_config.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'single finger dissimilarity';
                clr = [0 0 0];%[0.8 0.27 0.5];
                linestyle = ':';
            case 'LDC.mf.ipsi'
                % signal change across all stimulation configurations from
                % ipsilateral side
                dataFile = fullfile(groupDir,'R.glm4.LDC.mf_config.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'multi finger dissimilarity';
                clr = [0.8 0.27 0.5];
                linestyle = ':';
            case 'PSC.avg'
                % signal change across all stimulation configurations
                dataFile = fullfile(groupDir,'L.glm4.PSC.avg.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'average activity';
            case 'PSC.mf'
                % signal change across all multi-finger configurations
                dataFile = fullfile(groupDir,'L.glm4.PSC.mf_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'average activity';
            case 'PSC.1d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.1d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
                clr = [0.2 0.14 0.53];
            case 'PSC.2d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.2d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
                clr = [0.19048 0.80952 1];
                %clr = [0.402889876768966,0.000755861497321285,0.656615372714286];
            case 'PSC.3d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.3d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
                clr = [0.38095 0.61905 1];
                %clr = [0.670513481428572,0.143621099597743,0.581658503922099];
            case 'PSC.4d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.4d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity'; 
                clr = [0.57143 0.42857 1];
                %clr = [0.858834543714286,0.359293758510337,0.407883583231319];
            case 'PSC.5d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.5d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
                clr = [0.7619 0.2381 1];
                %clr = [0.976835837428571,0.598145892487094,0.243802292730916];
            case 'PSC.avg.ipsi'
                % signal change across all stimulation configurations from
                % ipsilateral side
                dataFile = fullfile(groupDir,'R.glm4.PSC.avg.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'average activity';
                linestyle = ':';
                clr = [0 0 0];
            case 'PSC.1d.ipsi'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'R.glm4.PSC.1d_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
                linestyle = ':';
            case 'PSC.mf.ipsi'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'R.glm4.PSC.mf_config.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'multi finger activity';
                linestyle = ':';
        end

        switch location
            case 'top'
                from = [8 95];
                to = [55 92];
                width = 10; 
            case 'middle'
%                 from = [3 78];
%                 to = [46 72];
%                 width = 10; 
                from = [4 80];
                to = [46 68];
                width = 26;
            case 'bottom'
                from = [0 58];
                to = [44 54];
                width = 10; 
            case 'OP1'
                from = [24 -6];
                to = [51 -6];
                width = 15;
        end
        
        % sample surface data
        [Y,P] = surf_cross_section(fullfile(atlasDir,'fs_LR.164k.L.flat.surf.gii'),dataFile,'from',from,'to',to,'width',width,'fcn',@(x) mean(x));
        Yb = surf_cross_section(fullfile(atlasDir,'fs_LR.164k.L.flat.surf.gii'),fullfile(atlasDir,'ROI_pp1.164k.L.label.gii'),'from',from,'to',to,'width',width,'fcn',@(x) mode(x));
        % Y=data
        % Yb=approximation of region boundaries
        
        % plot
        traceplot(1:100,Y','errorfcn','stderr','linecolor',clr,'linewidth',2,'linestyle',linestyle,...
            'patchcolor',clr,'transp',0.2);
        drawline(0,'dir','horz');
        % draw lines at region boundaries:
        drawline(find(diff(Yb))','dir','vert');
        ylabel(yaxisLabel);
        title(titleLabel);
        
        % save cross section file to check what vertices are mapped
        if saveBorder
            P = surf_makeLabelGifti(P);
            save(P,fullfile(groupDir,sprintf('L.glm4.profile.%s.label.gii',location)));
        end
        
        %set(gca,'xtick',[6,21,38,55,71,88],'xticklabel',{'4a','4p','3a','3b','1','2'})
        
        varargout = {Y,Yb};
    case 'surf:profile_strip_OLD'
        rois = [5,6,1,2,3,4];
        
        % activation profile for stimulation
        map = 'LDC.avg';
        vararginoptions(varargin,{'map','location'});
        saveBorder = 1;
        groupDir = fullfile(wbDir, ['group_164k']);
        regLabel = gifti(fullfile(atlasDir,'ROI_pp1.164k.L.label.gii')); % load label file to mask mask for profile plot
        switch map
            case 'LDC.avg'
                % dissimilarities avg. across all stimulation
                % configurations
                dataFile = fullfile(groupDir,'L.glm4.LDC.avg.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'average dissimilarity';
            case 'LDC.1d'
                % dissimilarities b/t single fingers
                dataFile = fullfile(groupDir,'L.glm4.LDC.1d_config.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'single finger dissimilarity';
            case 'PSC.avg'
                % signal change across all stimulation configurations
                dataFile = fullfile(groupDir,'L.glm4.PSC.avg.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'average activity';
            case 'PSC.1d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.1digit.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
            case 'PSC.2d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.2digit.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
            case 'PSC.3d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.3digit.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
            case 'PSC.4d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.4digit.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
            case 'PSC.5d'
                % signal change for single fingers
                dataFile = fullfile(groupDir,'L.glm4.PSC.5digit.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
            case 'PSC.avg.ipsi'
                % signal change across all stimulation configurations from
                % ipsilateral side
                dataFile = fullfile(groupDir,'R.glm4.PSC.avg.164k.func.gii');
                yaxisLabel = '% signal change';
                titleLabel = 'average activity';
            case 'LDC.avg.ipsi'
                % dissimilarities avg. across all stimulation
                % configurations
                dataFile = fullfile(groupDir,'R.glm4.LDC.avg.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'average dissimilarity';
            case 'LDC.1d.ipsi'
                % dissimilarities b/t single fingers
                dataFile = fullfile(groupDir,'R.glm4.LDC.1d_config.164k.func.gii');
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'single finger dissimilarity';
            
        end
        
        % plotting colours
        clrs = {[0 0 0.8] [0.3 0.6 0.9] [0.6 0.6 0.6] [0.5 0 0] [0.9 0 0] [1 0.6 0]};
        % get data for regions
        from = [31 104];
        to = [19 45];
        width = 25;        
        D=[];
        for ii=1:numel(rois)
            mask  = regLabel.cdata==rois(ii);
            % sample surface data
            [Yr,P] = surf_cross_section(fullfile(atlasDir,'fs_LR.164k.L.flat.surf.gii'),dataFile,'from',from,'to',to,'width',width,'fcn',@(x) mean(x),'mask',mask);
            % save cross section file to check what vertices are mapped
            if saveBorder
                P = surf_makeLabelGifti(P);
                save(P,fullfile(groupDir,sprintf('L.glm4.profile.%d.label.gii',rois(ii))));
            end
            d=[];
            d.sn = [2:11]';
            d.roi = ones(10,1).*rois(ii);
            d.glm = ones(10,1).*4;
            d.data = Yr';
            D=addstruct(D,d);
            % plot:
            subplot(1,numel(rois),ii);
            traceplot(1:100,d.data,'linecolor',clrs{ii},'patchcolor',clrs{ii},'errorfcn','stderr','linewidth',2,'transp',0.2);
            drawline(0,'dir','horz');
            ylabel(yaxisLabel);
            % flip figure to be vertical
            view(-90,90);
            set(gca, 'ydir', 'reverse');
            title(regname{rois(ii)});
            if ii>1
                a=get(gca);
                a.XAxis.Visible = 'off';
            end
        end    
        plt.match('y')
        varargout={D};
    case 'surf:profile_strip'
        rois = [5,6,1,2,3,4];
        
        % activation profile for stimulation
        map = 'LDC.avg';
        vararginoptions(varargin,{'map','location'});
        saveBorder = 1;
        groupDir = fullfile(wbDir, ['group_164k']);
        regLabel = gifti(fullfile(atlasDir,'ROI_pp1.164k.L.label.gii')); % load label file to mask mask for profile plot
        switch map
            case 'LDC.avg'
                % dissimilarities avg. across all stimulation
                % configurations
                dataFile{1} = fullfile(groupDir,'L.glm4.LDC.avg.164k.func.gii'); % contra
                dataFile{2} = fullfile(groupDir,'R.glm4.LDC.avg.164k.func.gii'); % ipsi
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'average dissimilarity';
            case 'LDC.1d'
                % dissimilarities b/t single fingers
                dataFile{1} = fullfile(groupDir,'L.glm4.LDC.1d_config.164k.func.gii'); % contra
                dataFile{2} = fullfile(groupDir,'R.glm4.LDC.1d_config.164k.func.gii'); % ipsi
                yaxisLabel = 'dissimilarity (a.u.)';
                titleLabel = 'single finger dissimilarity';
            case 'PSC.avg'
                % signal change across all stimulation configurations
                dataFile{1} = fullfile(groupDir,'L.glm4.PSC.avg.164k.func.gii'); % contra
                dataFile{2} = fullfile(groupDir,'R.glm4.PSC.avg.164k.func.gii'); % ipsi
                yaxisLabel = '% signal change';
                titleLabel = 'average activity';
            case 'PSC.1d'
                % signal change for single fingers
                dataFile{1} = fullfile(groupDir,'L.glm4.PSC.1digit.164k.func.gii'); % contra
                dataFile{2} = fullfile(groupDir,'R.glm4.PSC.1digit.164k.func.gii'); % ipsi
                yaxisLabel = '% signal change';
                titleLabel = 'single finger activity';
        end
        
        % plotting colours
        clrs = {[0 0 0.8] [0.3 0.6 0.9] [0.6 0.6 0.6] [0.5 0 0] [0.9 0 0] [1 0.6 0]};
        % get data for regions
        from = [31 104];
        to = [19 45];
        width = 25;        
        D=[];
        for ii=1:numel(rois)
            Dh = [];
            mask  = regLabel.cdata==rois(ii);
            for hh=1:2 % contra, ipsi hemispheres
                % sample surface data
                Yr = surf_cross_section(fullfile(atlasDir,'fs_LR.164k.L.flat.surf.gii'),dataFile{hh},'from',from,'to',to,'width',width,'fcn',@(x) mean(x),'mask',mask);
                
                d=[];
                d.sn = [2:11]';
                d.roi = ones(10,1).*rois(ii);
                d.hemi = ones(10,1).*hh; % 1=contra, 2=ipsi
                d.glm = ones(10,1).*4;
                d.data = Yr';
                Dh=addstruct(Dh,d);
            end
            
            % plot style:
            CAT.linestyle = {'-',':'};
            CAT.linecolor = clrs{ii};
            CAT.patchcolor = clrs{ii};
            CAT.linewidth = 2;
            
            % plot:
            subplot(1,numel(rois),ii);
            traceplot(1:100,Dh.data,'split',Dh.hemi,'errorfcn','stderr','CAT',CAT,'transp',0.2);
            drawline(0,'dir','horz');
            ylabel(yaxisLabel);
            
            % flip plot to be vertical
            view(-90,90);
            set(gca, 'ydir', 'reverse');
            title(regname{rois(ii)});
            if ii>1
                a=get(gca);
                a.XAxis.Visible = 'off';
            end
            
            D=addstruct(D,Dh);
        end    
        plt.match('y')
        varargout={D};
        
    case 'surf:patterns'
        I = pp1_imana('LIST_subjs');
        % create finger pattern pictures from workbench files
        sn  = 4;
        glm = 4;
        conds = 1:5;
        map = 0;
        vararginoptions(varargin,{'sn','glm','conds','map'});
        
        subjName = I.origSN{I.sn==sn};
        warning off
        h = 1; % left hemi
        groupDir = [gpCaretDir filesep hemName{h} ];
        cd(groupDir);
        switch(h)
            case 1 % left hemi
                F = gifti(fullfile(atlasDir,sprintf('fs_LR.164k.%s.flat.surf.gii',hemLetter{h})));
%                 xlims=[-4 44];
%                 ylims=[55 95];
                xlims=[0 50];
                ylims=[50 100];
        end
        
        % Get the X-Y coordinates for all tiles 
        faces = double(F.faces);
        verts = double(F.vertices);
        for ii=1:3
            X(ii,:) = verts(faces(:,ii),1);
            Y(ii,:) = verts(faces(:,ii),2);
        end

        % Find all tiles that have a single vertex (or more) in the image 
        M.k=find(any(X>xlims(1) & X<xlims(2),1) & any(Y>ylims(1) & Y<ylims(2),1)); 
        M.X=X(:,M.k); 
        M.Y=Y(:,M.k); 

        M.xlims=xlims; 
        M.ylims=ylims; 
        
        T.data = faces;

        % get region borders
        Bg = gifti(fullfile(atlasDir,'ROI_pp1_border.func.gii'));
        borderVert = logical(sum(double(Bg.cdata),2));
        B(1).data(:,1) = verts(borderVert,1); % border x-coord;
        B(1).data(:,2) = verts(borderVert,2); % border y-coord;
        
        % plot patterns:
        mm = 2.25; % force colour scaling on patterns
        G = gifti(fullfile(wbDir, subjName, sprintf('%s.%s.glm%d.TMAPS.164k.func.gii', subjName, hemLetter{h}, glm)));
        for ii=1:numel(conds)
            subplot(1,numel(conds),ii);
            caret_plotflatmap('data',double(G.cdata(:,conds(ii))),'M',M,'T',T,'xlims',xlims,'ylims',ylims,'border',B,'bordercolor',{'w.'});
            title(sprintf('chord %d',conds(ii)))
            colormap('parula');
            caxis([-mm/2 mm]);   % scale color across plots
        end

        % plot schematic?
        if map
            % Get sulcual depths and patterns to plot:
            ref = spm_vol(fullfile(atlasDir,'/fs_LR.164k.LR.sulc.dscalar.nii'));
            ref = ref.private.dat(:);
            sulc(:,1) = ref(1:163842); % left hem
            sulc(:,2) = ref(163843:end); % right hem
            % plot
            figure('Color',[1 1 1]);
            caret_plotflatmap('data',sulc(:,h),'M',M,'T',T,'xlims',xlims,'ylims',ylims,'border',B,'bordercolor',{'r.'});
            colormap('bone');
        end
        warning on
        
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
    case 'WRAPPER_searchlight'                                               
        % TO DO for glm4: sn 3,5,6,7,9,10,11
        
        glm = 4;
        sn = [];
        vararginoptions(varargin,{'sn','glm'});
        % You can call this case to run searchlight analyses.
        % 'sn' can be an array of subjects.
        %addpath('/Users/sarbuckle/MATLAB/email/')
        for s = sn 
            pp1_imana('SEARCH_define','sn',s,'glm',glm);
            %mySendmail(sprintf('done searchlight definition %s',subj_name{sn}));
            pp1_imana('SEARCH_ldcRun','sn',s,'glm',glm);
            %mySendmail(sprintf('done searchlight %s',subj_name{sn}));
            %pp1_imana('SEARCH_ldcContrasts','sn',s,'glm',glm);
        end
    case 'SEARCH_define'                                                    % STEP 4.1   :  Defines searchlights for 160 voxels in grey matter surface
        glm = [];
        sn = [];
        vararginoptions(varargin,{'sn','glm'});
        
        mask       = fullfile(glmDir{glm},subj_name{sn},'mask.nii');
        Vmask      = rsa.readMask(mask);

        LcaretDir = fullfile(caretDir,sprintf('x%s',subj_name{sn}),'LeftHem');
        RcaretDir = fullfile(caretDir,sprintf('x%s',subj_name{sn}),'RightHem');
        white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
        pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
        topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};
        S         = rsa_readSurf(white,pial,topo);

        L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[15 160]);
        save(fullfile(anatomicalDir,subj_name{sn},sprintf('%s_searchlight_160.mat',subj_name{sn})),'-struct','L');
    case 'SEARCH_ldcRun'                                                    % STEP 4.2   :  Runs LDC searchlight using defined searchlights (above)
        % Requires java functionality unless running on SArbuckle's
        % computer.
        glm = [];
        sn = [];
        vararginoptions(varargin,{'sn','glm'});
        
        block = 5e7;
        cwd   = pwd;                                                        % copy current directory (to return to later)
        % go to subject's glm directory 
        spmDir = fullfile(glmDir{glm},subj_name{sn});
        cd(spmDir);
        % load their searchlight definitions and SPM file
        L = load(fullfile(anatomicalDir,subj_name{sn},sprintf('%s_searchlight_160.mat',subj_name{sn})));
        load SPM;
        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{sn}));
        % make index vectors
        D = load('SPM_info.mat');
        conditionVec  = D.tt;
        partition     = D.run;
        name = sprintf('%s_glm%d',subj_name{sn},glm);
        % run the searchlight
        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partition,...
                                'analysisName',name,'idealBlock',block,...
                                'spmDir',spmDir);
        cd(cwd);
    case 'SEARCH_ldcContrasts'                                        
        % Calls 'MISC_SEARCH_calculate_contrast'
        sn  = [];
        glm = [];
        con = {'avg'};
        vararginoptions(varargin,{'sn','glm','con'});
        % Use 'con' option to define different contrasts.
        %   'avg'    :  Average LDC nii for all conds

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
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',1:5); 
                case '1d_config' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',1); 
                case '2d_config' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',2);
                case '3d_config'
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',3);
                case '4d_config' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',4);
                case 'mf_config' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',2:5);
                case 'numDigits'
                    gidx    = pp1_imana('GET_idxAcrossNumDigits','glm',glm);
            end
            % avg. distances according to contrast selected
            Y.LDC   = vdat(:,:,:,gidx);
            %Y.LDC   = ssqrt(Y.LDC);
            Y.LDC   = nanmean(Y.LDC,4); 
            % prep output file
            Y.dim   = vol(1).dim;
            Y.dt    = vol(1).dt;
            Y.mat   = vol(1).mat;    
            % save output
            Y.fname   = sprintf('%s_glm%d_LDC_%s.nii',subj_name{sn},glm,con{c});
            Y.descrip = sprintf('exp: ''pp1'' \nglm: ''FAST'' \ncontrast: ''%s''',con{c});

            spm_write_vol(Y,Y.LDC);
            fprintf('Done %s\n',Y.fname);

            clear Y
        end   
            
        
    case 'SEARCH_group_make'
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
    case 'SEARCH_group_cSPM'                                                
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
    case 'SEARCH_vol2Surf'                                                  
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
    case 'ROI_makeLabelGii'
        % takes the pp1 label gifti & adds SII roi that is made from SPM
        % anatomy toolbox maps (prob atlases) which are already aligned to the WB atlases.
        
        atlasDir = '/Users/sarbuckle/DATA/Atlas_templates/FS_LR_164';
        hem      = {'L','R'}; % both hemispheres
        hemName  = {'CortexLeft','CortexRight'};
        spmAnatDir = [fileparts(which('Anatomy.m')) filesep 'PMaps']; % where the anatomy niftis for spm toolbox are saved (we use the operculum images to make sII roi)
        
        % (1) convert caret paint files to .label.gii:
        inFile = '/Users/sarbuckle/DATA/passivePatterns1/fmri/surfaceCaret/fsaverage_sym/LeftHem/ROI_pp1.paint'; % only left hemi since we mirror left for right hemi
        surf_resamplePaint('/Users/sarbuckle/DATA/Atlas_templates',inFile,'ROI_pp1');
        
        % (2) harvest the spm anatomy toolbox images (toss out niftis 
        % of amygdala, cerebellum, hippocampus, thalamus, & basil
        % forebrain rois):
        D=dir(spmAnatDir);
        D=D(4:end); % get nifti names only
        exclude = {'Amygdala','Bforebrain','Cerebellum','Hippocampus','wThalamus'};
        idx = zeros(size(D));
        for ii=1:numel(exclude)
            idx = idx + cell2mat(arrayfun(@(x) startsWith(x,exclude{ii}),{D.name}','uni',0));
        end
        D = D(~logical(idx)); % files to keep
        % gather filenames into appropriate format:
        fnames = {};
        col_names = {};
        for ii=1:numel(D)
            fnames{end+1}    = [fullfile(D(ii).folder,D(ii).name) ',1'];
            col_names{end+1} = D(ii).name;
        end
        
        % loop through hemis
        for hh=1:2 % left and right hemi
            % (3) find vertices for SII roi on this hemi
            % load atlas surfaces
            C1 = gifti(fullfile(atlasDir,['fs_LR.164k.' hem{hh} '.white.surf.gii'])); % white surf
            C2 = gifti(fullfile(atlasDir,['fs_LR.164k.' hem{hh} '.pial.surf.gii']));  % pial surf
            % map roi probabilities to functional gifti:
            G = surf_vol2surf(C1.vertices, C2.vertices, fnames,...
                'column_names',col_names, 'anatomicalStruct',hemName{hh});
            % as per Caret roi definition, since these images from spm
            % anatomy toolbox are probabilistic atlases, we ensure the
            % vertices are assigned to roi with highest probability
            [Prob,ROI] = max(G.cdata,[],2); 
            % also ensure cytoarchitectonic prob (>0.2)
            ROI(Prob<0.25) = 0;
            % now restrict verticies to be only those that are part of
            % Operculum 1-4 rois:
            roi = find(cell2mat(arrayfun(@(x) startsWith(x,'Operculum'),{D.name}','uni',0))==1); % which roi# are the operculum subdivisions?
            s2_vert = ismember(ROI,roi);
            
            % (3) load current rois and include them in new roi file:
            inFile = fullfile(atlasDir,['ROI_pp1.164k.' hem{hh} '.label.gii']);
            C = gifti(inFile);
            labels = C.labels.name;
            labels{end+1} = 'SII';
            rgba = [0 0 0 0;...
                255 41 255 255;...
                41 41 255 255;...
                41 255 41 255;...
                0 204 255 255;...
                255 255 41 255;...
                204 153 0 255;...
                255 95 41 255;...
                230 41 95 255]; 
            rgba=rgba./255;

            numRoi = max(unique(C.cdata));
            newIdx = C.cdata==0 & s2_vert==1;
            C.cdata(newIdx) = int32(numRoi+1);

            Cn = surf_makeLabelGifti(C.cdata,...
                    'columnNames',labels,...
                    'labelNames',labels,...
                    'labelRGBA',rgba,...
                    'anatomicalStruct',hemName{hh});

            save(Cn,inFile);
            fprintf('done %s\n',hem{hh});
        end
    case 'ROI_makeLabelGii_OP1'
        % makes paint file for OP1 roi
        
        atlasDir = '/Users/sarbuckle/DATA/Atlas_templates/FS_LR_164';
        hem      = {'L','R'}; % both hemispheres
        hemName  = {'CortexLeft','CortexRight'};
        spmAnatDir = [fileparts(which('Anatomy.m')) filesep 'PMaps']; % where the anatomy niftis for spm toolbox are saved (we use the operculum images to make sII roi)
        
        % (1) convert caret paint files to .label.gii:
        inFile = '/Users/sarbuckle/DATA/passivePatterns1/fmri/surfaceCaret/fsaverage_sym/LeftHem/ROI_pp1.paint'; % only left hemi since we mirror left for right hemi
        %surf_resamplePaint('/Users/sarbuckle/DATA/Atlas_templates',inFile,'ROI_pp1');
        
        % (2) harvest the spm anatomy toolbox images (toss out niftis 
        % of amygdala, cerebellum, hippocampus, thalamus, & basil
        % forebrain rois):
        D=dir(spmAnatDir);
        D=D(4:end); % get nifti names only
        exclude = {'Amygdala','Bforebrain','Cerebellum','Hippocampus','wThalamus'};
        idx = zeros(size(D));
        for ii=1:numel(exclude)
            idx = idx + cell2mat(arrayfun(@(x) startsWith(x,exclude{ii}),{D.name}','uni',0));
        end
        D = D(~logical(idx)); % files to keep
        % gather filenames into appropriate format:
        fnames = {};
        col_names = {};
        for ii=1:numel(D)
            fnames{end+1}    = [fullfile(D(ii).folder,D(ii).name) ',1'];
            col_names{end+1} = D(ii).name;
        end
        
        % loop through hemis
        for hh=1:2 % left and right hemi
            % (3) find vertices for SII roi on this hemi
            % load atlas surfaces
            C1 = gifti(fullfile(atlasDir,['fs_LR.164k.' hem{hh} '.white.surf.gii'])); % white surf
            C2 = gifti(fullfile(atlasDir,['fs_LR.164k.' hem{hh} '.pial.surf.gii']));  % pial surf
            % map roi probabilities to functional gifti:
            G = surf_vol2surf(double(C1.vertices), double(C2.vertices), fnames,...
                'column_names',col_names, 'anatomicalStruct',hemName{hh});
            % as per Caret roi definition, since these images from spm
            % anatomy toolbox are probabilistic atlases, we ensure the
            % vertices are assigned to roi with highest probability
            [Prob,ROI] = max(G.cdata,[],2); 
            % also ensure cytoarchitectonic prob (>0.2)
            ROI(Prob<0.25) = 0;
            % now restrict verticies to be only those that are part of
            % Operculum 1-4 rois:
            roi = find(cell2mat(arrayfun(@(x) startsWith(x,'Operculum_OP1'),{D.name}','uni',0))==1); % which roi# are the operculum subdivisions?
            s2_vert = ismember(ROI,roi);
            
            % (3) load current rois and include them in new roi file:
            inFile = fullfile(atlasDir,['ROI_pp1.164k.' hem{hh} '.label.gii']);
            C = gifti(inFile);
            labels = C.labels.name;
            labels{end} = 'OP1';
            rgba = [0 0 0 0;...
                69 114 180 255;...
                116 173 209 255;...
                254 224 144 255;...
                253 174 97 255;...
                244 109 67 255;...
                215 48 39 255;...
                255 95 41 255]; 
            rgba=rgba./255;
            
            % change current OP1+OP2+OP3+OP4 region into just OP1
            numRoi = max(unique(C.cdata));
            C.cdata(C.cdata==numRoi) = 0;
            newIdx = C.cdata==0 & s2_vert==1;
            C.cdata(newIdx) = int32(numRoi);

            Cn = surf_makeLabelGifti(C.cdata,...
                    'columnNames',labels,...
                    'labelNames',labels,...
                    'labelRGBA',rgba,...
                    'anatomicalStruct',hemName{hh});

            outFile = fullfile(atlasDir,['ROI_pp1_OP1.164k.' hem{hh} '.label.gii']);
            save(Cn,outFile);
            fprintf('done %s\n',hem{hh});
        end
    case 'ROI_define'                                                       % Define rois: BA rois, M1/S1 cut to hand area, BA rois cut to hand area
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 
        I = pp1_imana('LIST_subjs');
        sn = pp1_imana('getSubjs');
        vararginoptions(varargin,{'sn'});
        
        thalamicMarker = [8133,8115,8109,8233,8215,8209]; % values of voxels assigned to each of the 3 thalamic regions (vpl,mgn,lgn) split by hemisphere
        corticalMarker = [1,2,3,4,5,6,7,1,2,7,8]; % values of surface nodes assigned to the cortical rois
        corticalFile   = {'D1','D1','D1','D1','D1','D1','D1','D2','D2','D2','D2'}; % different files for coritcal rois
        
        for s = sn % for each subject
            R = {};
            j = 1; % overall region ticker
            t = 1; % thalamic ticker
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            mask      = fullfile(glmDir{1},subjName,'mask.nii,1');  % load mask file now 
            for h = 1:2 % per hemisphere
                wbFilePrefix = fullfile(wbDir,subjName,[subjName '.' hemLetter{h}]);
                regFile1 = fullfile(atlasDir,['ROI_pp1.164k.' hemLetter{h} '.label.gii']);
                regFile2 = fullfile(atlasDir,['ROI.164k.' hemLetter{h} '.label.gii']);
                D1       = gifti(regFile1);
                D2       = gifti(regFile2);
                for r = 1:numregions % make regions
                    % make R region structure for participant
                    R{j}.name     = [subjName '_' regname{r} '_' hem{h}];
                    R{j}.regNum   = j;
                    R{j}.regType  = r;
                    R{j}.cortical = cortical(r);
                    R{j}.thalamic = thalamic(r);
                    R{j}.hemi     = h;
                    % Mapping of rois is different if cortical or thalamic:
                    if cortical(r)==1 && thalamic(r)==0 % cortical surface rois
                        eval(['D=' corticalFile{r} ';']);
                        R{j}.type     = 'surf_nodes_wb';
                        R{j}.location = find(D.cdata(:,1)==corticalMarker(r));
                        R{j}.white    = [wbFilePrefix '.white.164k.surf.gii'];
                        R{j}.pial     = [wbFilePrefix '.pial.164k.surf.gii'];
                        R{j}.linedef  = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
                        R{j}.image    = mask;    % functional mask
                        if strcmp(corticalFile{r},'D1')
                            R{j}.origFile = regFile1;
                        elseif strcmp(corticalFile{r},'D2')
                            R{j}.origFile = regFile2;
                        end
                    elseif thalamic(r)==1 && cortical(r)==0 % thalamic rois
                        % make temporary functional mask of thalamic nuclei
                        R{j}.type   = 'roi_image';
                        R{j}.value  = thalamicMarker(t);
                        R{j}.file   = fullfile(anatomicalDir,subjName,sprintf('%s_thalamicNuclei.nii',subjName));
                        t = t+1;
                    end
                    % update roi ticker
                    j = j+1;
                end    
            end
            exculdePairs = [1,5; 1,6; 2,5; 2,6; 3,5; 3,6; 8,9]; % exclude voxels that span across CS
            excludePairs = [exculdePairs; exculdePairs+numregions]; % do this exclusion for rois in both hemispsheres
            R = region_calcregions(R,'exclude',excludePairs,'exclude_thres',0.75);
            cd(regDir);
            save(['regions_' subjName '.mat'],'R');
            fprintf('..done\n');
            clear R
        end
    case 'ROI_define_OP1'                                               
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 
        I = pp1_imana('LIST_subjs');
        sn = pp1_imana('getSubjs');
        vararginoptions(varargin,{'sn'});
        
        for s = sn % for each subject
            R = {};
            j = 1; % overall region ticker
            t = 1; % thalamic ticker
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            mask      = fullfile(glmDir{1},subjName,'mask.nii,1');  % load mask file now 
            for h = 1:2 % per hemisphere
                wbFilePrefix = fullfile(wbDir,subjName,[subjName '.' hemLetter{h}]);
                regFile = fullfile(atlasDir,['ROI_pp1_OP1.164k.' hemLetter{h} '.label.gii']);
                D       = gifti(regFile);
                % make R region structure for participant
                    R{j}.name     = [subjName '_OP1_' hem{h}];
                    R{j}.regNum   = j;
                    R{j}.hemi     = h;
                    R{j}.type     = 'surf_nodes_wb';
                    R{j}.location = find(D.cdata(:,1)==7);
                    R{j}.white    = [wbFilePrefix '.white.164k.surf.gii'];
                    R{j}.pial     = [wbFilePrefix '.pial.164k.surf.gii'];
                    R{j}.linedef  = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
                    R{j}.image    = mask;    % functional mask
                    R{j}.origFile = regFile;
                    % update roi ticker
                    j = j+1;
            end
            R = region_calcregions(R,'exclude_thres',0.75);
            cd(regDir);
            save(['regions_OP1' subjName '.mat'],'R');
            fprintf('..done\n');
            clear R
        end
        
    case 'ROI_make_nii'                                                     % OPTIONAL   :  Convert ROI def (.mat) into multiple .nii files (to check!)
        I = pp1_imana('LIST_subjs');
        sn = pp1_imana('getSubjs');
        glm = 4;
        vararginoptions(varargin,{'sn'});        
        
        for s=sn
            subjName = I.origSN{s};
            glmSubjDir = fullfile(glmDir{glm},subjName);

            cd(glmSubjDir);
            % load ROI definition
            load(fullfile(regDir,['regions_' subjName '.mat']));

            % loop over rois
            for rr = [1:28]
                % mask volume
                mask = fullfile(glmSubjDir,'mask.nii');           
                % Save region file as nifti
                cd(regDir);
                region_saveasimg(R{rr},mask);      
            end
            
        end        
    case 'ROI_fitHRF'
        % p    - parameters of the response function (two Gamma functions)
        %                                                          defaults  (seconds)  
        %        p(1) - delay of response (relative to onset)          6
        %        p(2) - delay of undershoot (relative to onset)       16
        %        p(3) - dispersion of response                         1
        %        p(4) - dispersion of undershoot                       1
        %        p(5) - ratio of response to undershoot                6
        %        p(6) - onset (seconds)                                0
        %        p(7) - length of kernel (seconds)                    32

        sn  = 2;
        glm = 4;
        roi = 8; % optimize for S1
        Y   = [];
        vararginoptions(varargin,{'sn','glm','Y'});
        S = pp1_imana('LIST_subjs');
        
        fit = []'; % hrf parameter(s) to be fitted
        if glm==3
            subj_hrfParams = pp1_imana('GLM_hrfParams');
            P0 = subj_hrfParams{sn};
        else
            subj_hrfParams = pp1_imana('GLM_hrfParamsGLM4');
            P0 = subj_hrfParams{sn};
            %P0(1) = 2;
        end
        P0         = P0(1:6)';
        LB         = [2 5 0 0 0 -2]';    
        UB         = [7 16 10 10 10 5]'; 
        %duration   = 4; % duration of 1== unscaled duration (this value scales the duration of the boxcar model by proportion)
        onsetshift = 0;
        eCriteria  = 0.95;
        numIter    = 10;
        
        pre  = 2;
        post = 20;

        warning off
        % display to user which subject we are fitting
        subjName = S.origSN{S.sn==sn};
        fprintf('%s\n',subjName);
        % load appropriate subject data
        cd(fullfile(glmDir{glm},subjName));
        load(fullfile(glmDir{glm},subjName,'SPM.mat'));
        load(fullfile(regDir,sprintf('regions_%s',subjName)));
        % Get data
%         T = load(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)));
%         T = getrow(T,ismember(T.roi,roi) & T.sn==s);
        if isempty(Y)
            Y = region_getts(SPM,{R{roi}});
            Y = mean(Y,2);  % avg. across rois and hemis
        end
        Ypre = spm_filter(SPM.xX.K,SPM.xX.W*Y);    
        Yres = spm_sp('r',SPM.xX.xKXs,Ypre);
        Epre = sum(sum(Yres.^2))/numel(Yres(:));

        % Fit a common hrf
        eRatio = 1;
        P0_    = P0;
        iter   = 1;
        eRatio = 1;
        if ~isempty(fit)
            fprintf('Iteration...');
            while ((eRatio >= eCriteria)&&(iter<numIter))
                fprintf('%d.',iter);
                % fit hrf
                [P,SPM,Yhat,Yres] = spmj_fit_hrf(SPM,Y,...
                    'fit',fit,'LB',LB,'UB',UB,'P0',P0_);
                % update initial value
                P0_(fit) = P0(fit)+0.1*rand*(UB(fit)-LB(fit));
                iter  = iter+1;
                % Check Error after
                Epost  = sum(sum(Yres.^2))/numel(Yres(:));
                eRatio = Epost/Epre;
            end
        elseif isempty(fit)
            % not fitting, just checking current hrf param fit 
            if (length(P0)<7)
                 P0(7)=SPM.Sess(1).U(1).dur(1); 
            end
            P = P0;
            SPM.xBF.bf = spmj_hrf(SPM.xBF.dt,P(1:7));
            SPM=spmj_fMRI_design_changeBF(SPM);
            % return predicted timeseries and residuals 
            % Filter and prepare the data 
            Yy = spm_filter(SPM.xX.K,SPM.xX.W*Y);
            beta  = SPM.xX.pKX*Yy;                    %-Parameter estimates
            Yres  = spm_sp('r',SPM.xX.xKXs,Yy);                        % get the 
            reg_interest=[SPM.xX.iH SPM.xX.iC]; 
            Yhat   = SPM.xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 
        end
        % Parameter values
        fprintf('Epost/Epre: %1.5f\n',eRatio);
        display(P)
        
        % get timeseries
        D     = spmj_get_ons_struct(SPM); % get onsets (onsets are in img #)
        D.sn  = ones(size(D.event,1),1)*sn;
        y_hat = mean(Yhat,2);
        y_res = mean(Yres,2);
        for i=1:size(D.block,1)
            D.y_hat(i,:)=cut(y_hat,pre,round(D.ons(i)),post,'padding','nan')';
            D.y_res(i,:)=cut(y_res,pre,round(D.ons(i)),post,'padding','nan')';
            D.y_adj(i,:)=D.y_hat(i,:)+D.y_res(i,:);
        end
        D=getrow(D,D.event<32); % only chord regressors
        % plot fits
        hold off;
        traceplot([-pre:post],D.y_adj,'errorfcn','stderr');
        hold on;
        traceplot([-pre:post],D.y_res,'linecolor',[0 1 0],'linewidth',3);
        traceplot([-pre:post],D.y_hat,'linecolor',[1 0 0],'linewidth',3);
        xlabel('seconds');
        ylabel('activation');
        xt = get(gca,'xtick');
        set(gca,'xticklabel',xt.*TR_length);
        xlim([-pre post]);
        % draw lines denoting trial events:
        drawline(0,'dir','vert'); % go cue onset
        drawline(6/TR_length,'dir','vert');
        hold off;
        
        keyboard
    case 'ROI_getTimeseries'                                                % (optional) :  Harvest ROI timeseries for specified region.
        % Use this and 'ROI_plot_timeseries' to ensure good GLM fits with
        % measured BOLD in rois.
        
        % Defaults
        sn  = 2:11;
        glm = 4;
        roi = [2:6];
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
                [y_raw, y_adj, y_hat, y_res] = region_getts(SPM,R{reg});          % get SPM info for voxels contained in specified region
                D = spmj_get_ons_struct(SPM);                                     % get trial onsets in TRs- because model was in secs, spmj converts onsets to TR #s by dividing time/TR length (320 trials to 4872 TRs)
                for r=1:size(y_raw,2)
                    for i=1:size(D.block,1)                                    % extract the timeseries of each trial from y_adj, y_hat, & y_res
                        D.y_adj(i,:) = cut(y_adj(:,r),pre,round(D.ons(i)) -1,post,'padding','nan')';
                        D.y_hat(i,:) = cut(y_hat(:,r),pre,round(D.ons(i)) -1,post,'padding','nan')';
                        D.y_res(i,:) = cut(y_res(:,r),pre,round(D.ons(i)) -1,post,'padding','nan')';
                        D.y_raw(i,:) = cut(y_raw(:,r),pre,round(D.ons(i)) -1,post,'padding','nan')';
                    end
                    D.glm = ones(size(D.event,1),1)*glm;
                    D.roi = ones(size(D.event,1),1)*reg;
                    D.sn  = ones(size(D.event,1),1)*s;
                    T=addstruct(T,D);
                end
                fprintf('s%02d roi %d glm %d \n',s,reg,glm);
            end
        end
        save(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)),'-struct','T');
        
        %__________________________________________________________________
    case 'ROI_plotTimeseriesAvg'                                            % (optional) :  Plots timeseries for specified region (avg. across digits, separated by pressing speed)
        glm = 4;
        sn  = 2:11;
        roi = 2;
        vararginoptions(varargin,{'sn','glm','roi'});
        
        D=load(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)));
        D=getrow(D,ismember(D.sn,sn) & D.roi==roi & ismember(D.event,1:31));
        D=tapply(D,{'glm','roi','sn'},{'y_adj','nanmean(x,1)'},{'y_hat','nanmean(x,1)'});
        ii=1;
%         for s=sn
%             T=getrow(D,D.sn==s & D.roi==roi); % only plot chord regressors
%             
%             % plot timeseries
%             plt.trace([-4:20],T.y_adj,'errorfcn','stderr');
%             hold on;
%             traceplot([-4:20],T.y_hat,'linewidth',2,'linecolor',[1 0 0]);
%             hold off;
%             xlabel('TR');
%             ylabel('adjusted activation');
%             xlim([-4 11]);
%             title(sprintf('subj: %d  .', s));
%             legend off
%             ii = ii+1;
% 
% 
%             drawline(0,'dir','vert');
%             drawline(10,'dir','vert');
%             drawline(15,'dir','vert');
%             drawline(0,'dir','horz');
% 
%             
%         end
        % plot-friendly structure
        T=[];
        v=ones(2,1);
        for jj=1:size(D.sn,1)
            %d.bin = v.*[-4:20];
            d.bin = v.*[-3:20];
            d.y  = [D.y_adj(jj,2:end); D.y_hat(jj,1:end-1)];
            d.type = [1;2]; % 1=real, 2=predicted
            d.sn = v.*D.sn(jj);
            d.roi = v.*D.roi(jj);
            d.glm = v.*D.glm(jj);
            T=addstruct(T,d);
        end
        % plot
        sty = style.custom({'black','red'});
        plt.trace([-3:20],T.y,'split',T.type,'errorfcn','stderr','style',sty);
        % event lines
        drawline(0,'dir','vert');
        drawline(10,'dir','vert');
        drawline(15,'dir','vert');
        drawline(0,'dir','horz');
        
        varargout = {T};
        
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
        conds = 1:10;
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
                title(sprintf('s%02d tt:%d', s,d));
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
        % case to extract first-level regression coefficients (betas/
        % activity patterns) from rois per subject
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = [1:9];%[1:28];
        append = 1; % add betas to currently existing datastructure?
        vararginoptions(varargin,{'sn','glm','roi','append'});
        
        if glm==2
            numRegressors = 31; 
        elseif glm==3
            numRegressors = 32; % 31 chords + thumb response
        elseif glm==4
            numRegressors = 32; % 31 chords + thumb response
        end
        fprintf('extracting betas\n');
        T=[];
        if append
            fprintf('...and appending them to already existing file\n');
            T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        end
        % harvest
        for s=sn % for each subj
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            % load files
            load(fullfile(glmDir{glm}, subjName,'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            D = load(fullfile(glmDir{glm}, subjName,'SPM_info.mat'));
            load(fullfile(regDir,sprintf('regions_%s.mat',subj_name{s})));          % load subject's region parcellation & depth structure (R)
            % add percent signal change imgs for subject
            Q = {}; 
            for q = 1:numRegressors
                Q{q} = (fullfile(glmDir{glm}, subjName, sprintf('psc_%02d.nii',q))); 
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
                [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','runwise');
                % toss stuff into output structure
                S.sn                 = s;
                S.roi                = r;
                S.regSide            = regSide(r);
                S.regType            = regType(r);
                S.cortical           = cortical(S.regType);
                S.thalamic           = thalamic(S.regType);
                S.tt                 = {D.tt};
                S.run                = {D.run};
                S.numDigits          = {D.numDigits};
                S.sess               = {D.sess};
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
            fprintf('..done\n');
        end
        % save T
        save(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)),'-struct','T'); 
        fprintf('\n')
    case 'ROI_getBetas_OP1'                                                     % STEP 5.3   :  Harvest activity patterns from specified rois
        % case to extract first-level regression coefficients (betas/
        % activity patterns) from OP1 & save to data file per subject
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = 4;
        vararginoptions(varargin,{'sn','glm'});
        
        if glm==2
            numRegressors = 31; 
        elseif glm==3
            numRegressors = 32; % 31 chords + thumb response
        elseif glm==4
            numRegressors = 32; % 31 chords + thumb response
        end
        fprintf('extracting betas\n');
        T=[];
        
        % harvest
        for s=sn % for each subj
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            % load files
            load(fullfile(glmDir{glm}, subjName,'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            D = load(fullfile(glmDir{glm}, subjName,'SPM_info.mat'));
            load(fullfile(regDir,sprintf('regions_OP1%s.mat',subj_name{s})));          % load subject's region parcellation & depth structure (R)
            % add percent signal change imgs for subject
            Q = {}; 
            for q = 1:numRegressors
                Q{q} = (fullfile(glmDir{glm}, subjName, sprintf('psc_%02d.nii',q))); 
            end
            Q = spm_vol(char(Q));
            % TR img info
            V = SPM.xY.VY; 
            % remove run means from patterns
            C0         = indicatorMatrix('identity',D.run); 
            ofInterest = 1:size(C0,1); % indicies for regressors of interest
            
            for r = 1:2 % for the OP1 region in both hemis
                % get raw data/psc for voxels in region
                Y = region_getdata(V,R{r});  % Data Y is N x P (P is in order of transpose of R{r}.depth)
                P = region_getdata(Q,R{r});
                % estimate region betas
                [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','runwise');
                % toss stuff into output structure
                S.sn                 = s;
                S.roi                = r;
                S.regSide            = regSide(r);
                S.regType            = regType(r);
                S.cortical           = 1;
                S.thalamic           = 0;
                S.tt                 = {D.tt};
                S.run                = {D.run};
                S.numDigits          = {D.numDigits};
                S.sess               = {D.sess};
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
            fprintf('..done\n');
        end
        % save T
        for r=1:2
            t=getrow(T,T.roi==r);
            t.roi = t.roi+28;
            save(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,28+r)),'-struct','t');
        end
         
        fprintf('\n')
    
        
    case 'ROI_stats'                                                        % STEP 5.4   :  Calculate stats/distances on activity patterns
        glm = 4;
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        roi = [1:28];%unique(T.roi)';
        vararginoptions(varargin,{'glm'});
        % housekeeping
        numDigits = stimulationChords;
        numDigits = sum(numDigits,2)';
        if glm==2 % 32 chords
            numConds = 31;
        elseif glm==3 % 31 chords + thumb response regressor
            numConds = 32;
        elseif glm==4
            numConds = 32;
        end
        % output structures
        To = [];
        Td = [];
        % get data
        T  = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        
        % do stats
        for s = sn % for each subject
            D  = load(fullfile(glmDir{glm}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
            C0 = indicatorMatrix('identity',D.run);
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
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
                % cosine angles
                So.corr_dist = corr_crossval(rsa_squareIPM(So.G_wmean),'reg','abs');
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
                elseif glm>2 % regressors for chords and thumb response
                    So.psc_chord = [1:32]; 
                    So.psc_numD  = [numDigits,99];
                end
                % Calculate avg. betas for each condition
                Q            = [];
                Q.raw_beta   = S.raw_beta{1};
                Q.tt         = D.tt;
                Q            = tapply(Q,{'tt'},{'raw_beta','mean'});
                So.avg_betas = mean(Q.raw_beta,2)';
                So.avg_tt    = Q.tt';
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
            fprintf('..done\n');
        end % each subject

        % % save
        save(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)),'-struct','To');
        save(fullfile(regDir,sprintf('glm%d_reg_TnumDigits.mat',glm)),'-struct','Td');
        fprintf('done.\n')
    case 'ROI_statsSession'
        % calculate RDMs per session per participant
        glm = 3;
        I = pp1_imana('LIST_subjs');
        I = getrow(I,I.fmri_sessions==2);
        vararginoptions(varargin,{'glm'});
        % housekeeping
        if glm==2 % 32 chords
            numConds = 31;
        elseif glm==3 % 31 chords + thumb response regressor
            numConds = 32;
        end
        % output structure
        To = [];
        % get data
        T   = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        roi = unique(T.roi)';
        % do stats
        for s = unique(I.sn)' % for each subject
            D  = load(fullfile(glmDir{glm}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
            fprintf('\nSubject: %d\n',s)
            
            for sess=1:2 % for each session (2 per subj)
                d = getrow(D,D.sess==sess);
                d.run = d.run-(min(d.run)-1);
                C0         = indicatorMatrix('identity',D.run(D.sess==sess));
                num_run    = length(unique(D.run(D.sess==sess)));
                ofInterest = 1:(numConds*num_run); % indicies for regressors of interest
                for r = roi % for each region
                    S = getrow(T,(T.sn==s & T.roi==r)); % subject's region data
                    betaW       = S.betaW{1}(D.sess==sess,:);
                    betaW_nmean = betaW(ofInterest,:)-C0*pinv(C0)*betaW(ofInterest,:); % run mean subtraction  
                    % % Toverall structure stats
                    % crossval second moment matrix
                    [G,Sig]      = pcm_estGCrossval(betaW_nmean(ofInterest,:),d.run,d.tt);
                    So.sig       = rsa_vectorizeIPM(Sig);
                    So.G         = rsa_vectorizeIPM(G);
                    So.G_wmean   = rsa_vectorizeIPM(pcm_estGCrossval(betaW(ofInterest,:),d.run,d.tt));
                    % calculate empirical G (not cv)
                    tmp.betaW    = betaW_nmean;
                    tmp.tt       = d.tt;
                    tmp          = tapply(tmp,{'tt'},{'betaW','mean'});
                    So.G_emp     = rsa_vectorizeIPM(cov(tmp.betaW')); 
                    % squared dissimilarities
                    So.ldc_wmean = rsa.distanceLDC(betaW,d.run,d.tt);        % rdm crossvalidated, on patterns without run mean patterns removed
                    So.ldc       = rsa.distanceLDC(betaW_nmean,d.run,d.tt);  % rdm crossvalidated, patterns with run means removed
                    % do condition classification for all chords (NOTE: this pools data across sessions)
                    H = d;
                    H.beta = betaW;
                    H  = getrow(H,H.chord<32);
                    [So.nnEstAcc,So.nnChanceAcc,So.nnPval] = pp1_imana('ROI_nnClassify',H.beta,H.chord,H.run);
                    % do condition classification for single fingers 
                    H  = getrow(H,H.chord<6);
                    [So.nnEstAccSF,So.nnChanceAccSF,So.nnPvalSF] = pp1_imana('ROI_nnClassify',H.beta,H.chord,H.run);
                    % Calculate avg. betas for each condition
                    Q            = [];
                    Q.raw_beta   = S.raw_beta{1}(D.sess==sess,:);
                    Q.tt         = d.tt;
                    Q            = tapply(Q,{'tt'},{'raw_beta','mean'});
                    So.avg_betas = mean(Q.raw_beta,2)';
                    So.avg_tt    = Q.tt';
                    % indexing fields
                    So.sn       = s;
                    So.roi      = r;
                    So.sess     = sess;
                    So.numVox   = size(betaW,2);
                    So.regSide  = regSide(r);
                    So.regType  = regType(r);
                    To          = addstruct(To,So);
                    fprintf('%d.',r)
                end % each region
            end % each session
        end % each subject

        % % save
        save(fullfile(regDir,sprintf('glm%d_reg_ToverallSession.mat',glm)),'-struct','To');
        fprintf('done.\n')
    case 'ROI_nnClassify'
        % do nearest neighbour classifications. Also estimate chance level
        numIters = 100; % number of chance shuffles
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
    case 'ROI_pattConsist'                                                
        % Crossvalidated Pattern consistency is a measure of the proportion of explained
        % beta variance across runs within conditions. 
        % 
        % This stat is useful for determining which GLM model yeilds least
        % variable beta estimates. That GLM should be the one you use for
        % future analysis cases.

        glm = 4;
        sn = pp1_imana('getSubjs');
        roi = 1:4;
        removeMean = 1;
        conds = 1:31;
        vararginoptions(varargin,{'sn','glm','roi','removeMean','conds'});
        
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
    case 'ROI_pattConsist_perNumDigits'
        glm = 3;
        sn = pp1_imana('getSubjs');
        roi = 1:4;
        vararginoptions(varargin,{'sn','glm','roi'});

        numDigits = sum(pp1_imana('chords'),2);
        R = [];
        for d = 1:4 % don't do 5-finger chord (mean-removal kills this pattern)
            r = pp1_imana('ROI_pattConsist','glm',glm,'sn',sn,'roi',roi,'removeMean',1,'conds',find(numDigits==d));
            r.numDigits = ones(size(r.sn)).*d;
            R = addstruct(R,r);
        end
        varargout = {R};
    
    case 'ROI_calcPattCorr'                                           % plot w/in subj, w/in speed rdm reliability (Across two partitions), compare with across-speed correlations. Insights into RSA stability    
        % Calculates correlation coefficients witihin conditions.
        % Default setup includes subtraction of each run's mean
        % activity pattern (across conditions).
        Y       = varargin{1}; % activity patterns
        condVec = varargin{2};
        partVec = varargin{3};
        sessVec = varargin{4};
        splitBy = varargin{5}; % how to split the patterns ('oddEven','sess1','sess2','sess_cv_odd','sess_cv_even','sess_cv_evenOdd','sess_cv_oddEven')

        % Prep some variables (output structure, indexing variables, counter)
        Q = [];
        % loop across and correlate

        % cluster run numbers into partition splits
        switch splitBy
            case 'sess_cv_odd' % compare odd splits across sessions
                isOdd      = logical(mod(partVec,2));
                partitions = [{unique(partVec(isOdd & sessVec==1))},{unique(partVec(isOdd & sessVec==2))}];
            case 'sess_cv_even' % compare even splits across sessions
                isOdd      = logical(mod(partVec,2));
                partitions = [{unique(partVec(~isOdd & sessVec==1))},{unique(partVec(~isOdd & sessVec==2))}];
            case 'sess_cv_evenOdd' % compare even splits sess 1 to odd splits sess 2
                isOdd      = logical(mod(partVec,2));
                partitions = [{unique(partVec(~isOdd & sessVec==1))},{unique(partVec(isOdd & sessVec==2))}];
            case 'sess_cv_oddEven' % compare odd splits sess 1 to Even splits sess 2
                isOdd      = logical(mod(partVec,2));
                partitions = [{unique(partVec(isOdd & sessVec==1))},{unique(partVec(~isOdd & sessVec==2))}];
            case 'sess1'   % odd-even split within session 1
                isOdd      = logical(mod(partVec,2));
                partitions = [{unique(partVec(isOdd & sessVec==1))},{unique(partVec(~isOdd & sessVec==1))}];
            case 'sess2'   % odd-even split within session 2
                isOdd      = logical(mod(partVec,2));
                partitions = [{unique(partVec(isOdd & sessVec==2))},{unique(partVec(~isOdd & sessVec==2))}];
        end
        
        % avg. patterns across runs within partition split
        d.Y     = Y;
        d.cond  = condVec;
        d.split = zeros(size(condVec));
        d.split(ismember(partVec,partitions{1})) = 1;
        d.split(ismember(partVec,partitions{2})) = 2;
        d = getrow(d,d.split>0);
        
%         C0  = indicatorMatrix('identity',d.split);
%         d.Y = d.Y - C0*pinv(C0)* d.Y; % remove mean pattern within each partition split
        
        d = tapply(d,{'split','cond'},{'Y','mean'});
        % correlate patterns and harvest correlations accordingly
        sameCond = tril(bsxfun(@eq,d.cond,d.cond'),-1);
        diffPart = tril(bsxfun(@(x,y) x~=y,d.split,d.split'),-1);
        take     = sameCond & diffPart;
        R        = corr(d.Y');
        r        = mean(R(take));

        varargout = {r};  
    case 'ROI_compareReliability'
        % compares within-session pattern reliabilities to across-session
        % reliabilities. 
        % Reliabilities are odd-even splits.
        % corr(A_s1,B_s1) = r_s1
        % corr(A_s2,B_s2) = r_s2
        % corr(A_s1,B_s2) = r_AB
        % corr(B_s1,A_s2) = r_BA
        % ssqrt(r_s1 * r_s2) vs. ssqrt(r_A * r_B)
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = 4;
        vararginoptions(varargin,{'roi','glm','sn'});
        
        % Load betas
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        T = getrow(T,ismember(T.roi,roi) & ismember(T.sn,sn));
        % do calculations
        R = [];
        v = [1;1];
        for ii = 1:length(T.sn)
            condVec = T.tt{ii}(T.tt{ii}<32);
            partVec = T.run{ii}(T.tt{ii}<32);
            sessVec = T.sess{ii}(T.tt{ii}<32);
            if sn(ii)==11 % subj 11 was modeled with 4 glms b/c of breaks during both scanning sessions
                % sessions 1 and 2 were from first scanning session
                % sessions 3 and 4 were from second scanning session
                sessVec(sessVec==2) = 1;
                sessVec(sessVec>2) = 2;
            end
            % remove mean pattern (across conditions) from each run:
            Y  = T.raw_beta{ii}(T.tt{ii}<32,:); % take only chord regressors - Univariately whitened! (drop the thumb press and intercepts)
            C0 = indicatorMatrix('identity',partVec);
            Y  = Y - C0*pinv(C0)* Y; % run mean subtraction  
            % calc reliabilities
            r_sess1 = pp1_imana('ROI_calcPattCorr',Y,condVec,partVec,sessVec,'sess1');
            r_sess2 = pp1_imana('ROI_calcPattCorr',Y,condVec,partVec,sessVec,'sess2');
            r_E     = pp1_imana('ROI_calcPattCorr',Y,condVec,partVec,sessVec,'sess_cv_odd');
            r_O     = pp1_imana('ROI_calcPattCorr',Y,condVec,partVec,sessVec,'sess_cv_even');
            r_OE    = pp1_imana('ROI_calcPattCorr',Y,condVec,partVec,sessVec,'sess_cv_oddEven');
            r_EO    = pp1_imana('ROI_calcPattCorr',Y,condVec,partVec,sessVec,'sess_cv_evenOdd');

            within  = ssqrt(r_sess1 * r_sess2);
            %between = ssqrt(r_OE * r_EO);
            signR = sign(r_O*r_E*r_OE*r_EO);
            between = nthroot( abs(r_O*r_E*r_OE*r_EO) ,4)*signR; %mean([r_O,r_E,r_OE,r_EO]);
            % add to output structure
            r.corr  = [within; between];
            r.corrMore{1,1} = [r_sess1,r_sess2];
            r.corrMore{2,1} = [r_O,r_E,r_OE,r_EO];
            r.type  = [1;2];
            r.sn    = v.*T.sn(ii);
            r.roi   = v.*T.roi(ii);
            r.glm   = v.*glm;
            R = addstruct(R,r);
        end
        %plt.box(R.type,R.corr,'plotall',2);
        fprintf('-----within > 0-----\n');
        ttest(R.corr(R.type==1),[],1,'onesample');
        fprintf('-----between > 0-----\n');
        ttest(R.corr(R.type==2),[],1,'onesample');
        fprintf('-----bt vs wn paired-----\n');
        ttest(R.corr(R.type==1),R.corr(R.type==2),2,'paired')
        save(fullfile(regDir,sprintf('patternReliability_roi%d_glm%d.mat',roi,glm)),'-struct','R');
        varargout = {R};    
    case 'ROI_rdmStability'
        glm = 4;
        roi = [1:7];
        conds = [1:31]; % which conditions to compared distances from
        sn = pp1_imana('getSubjs');
        vararginoptions(varargin,{'glm','roi','sn'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        % housekeeping
        T = getrow(T,ismember(T.sn,sn));
        D = []; % output structure
        take = logical(tril(ones(numel(sn)),-1));
        for r = roi
            t = getrow(T,T.roi==r);
            if ~isempty(conds)
                ldc = [];
                for s=unique(t.sn)'
                    rdm = rsa_squareRDM(t.ldc(t.sn==s,:));
                    rdm = rdm(conds,conds);
                    ldc(end+1,:) = rsa_vectorizeRDM(rdm);
                end
            else
                ldc = t.ldc;
            end
            R = corr(ldc');
            d = [];
            d.numSN = numel(sn); 
            d.roi   = r;
            d.corrs = R(take)';
            d.corr  = mean(R(take));
            % calc confidence bounds
            rz        = fisherz(d.corrs)'; 
            d.is_mean = fisherinv(mean(rz));
            d.is_LB   = fisherinv(mean(rz) - 1.96*stderr(rz));
            d.is_UB   = fisherinv(mean(rz) + 1.96*stderr(rz));
            D = addstruct(D,d);
        end
        varargout = {D};  
    case 'ROI_MDS_overall'                                                  % (optional) :  Plots the scaled representational structure. 
        % enter region, glm #, sn (if desired)
        cplot = 'one';
        glm   = 4;
        roi   = 2; % default primary sensory cortex   
        sn    = pp1_imana('getSubjs');
        clrCode = 0; % if >0, color all chords with digit X red, and all other chords black.
        vararginoptions(varargin,{'roi','glm','cplot','clrCode','sn'});
        % cplot = 'all' to plot all 4 MDS figures (i.e. no contrast and 3 contrasts)- default
        % cplot = 'one'  to plot only no contrast MDS figure        

        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        T = getrow(T,T.roi==roi & ismember(T.sn,sn));
        IPM = nanmean(T.G_wmean,1); % take G_wmean- here, run means are retained
        G = rsa_squareIPM(IPM);
        G = G([1:31],[1:31]); % drop any non-chord regressors
        IPM = rsa_vectorizeIPM(G);
        
        [chords,labels]  = pp1_imana('chords');
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
                CnumDigits = indicatorMatrix('identity',[numDigits]);
                CnumDigits = bsxfun(@minus,CnumDigits,mean(CnumDigits,2));
                Cdigit     = indicatorMatrix('identity',[1:5]);
                Cdigit     = bsxfun(@minus,Cdigit,mean(Cdigit,2));
                Call       = eye(31)-ones(31)/31;

                G = rsa_squareIPM(IPM);
                G = G([1:5],[1:5]);
                IPMdigit = rsa_vectorizeIPM(G);
                
                Y{1} = rsa_classicalMDS(IPMdigit,'mode','IPM','contrast',Cdigit);
                Y{2} = rsa_classicalMDS(IPM,'mode','IPM');
                Y{3} = rsa_classicalMDS(IPM,'mode','IPM','contrast',Call);
                Y{4} = rsa_classicalMDS(IPM,'mode','IPM','contrast',CnumDigits);
                
                subplot(1,4,1);
                pp1_imana('MISC_scatterplotMDS',Y{1}(1:end,1:3),'split',split,'label',[1:5]');
                title('contrast: digit');
                subplot(1,4,2);
                pp1_imana('MISC_scatterplotMDS',Y{2}(1:end,1:3),'split',split,'label',labels');
                title('no contrast');
                subplot(1,4,3);
                pp1_imana('MISC_scatterplotMDS',Y{3}(1:end,1:3),'split',split,'label',labels');
                title('contrast: diff b/t all conds');
                subplot(1,4,4);
                pp1_imana('MISC_scatterplotMDS',Y{4}(1:end,1:3),'split',split,'label',labels');
                title('contrast: num digits');
                
            case 'one' % only do and plot no contrast MDS
                Y{1} = rsa_classicalMDS(IPM,'mode','IPM');
                
                % plot variance explained
                subplot(1,2,1);
                totvar = sum(sum(Y{1}.*Y{1}));
                contribution = sum(Y{1}.*Y{1},1)/totvar*100;
                bar(contribution);xlabel('Dimensions');ylabel('%variance explained');
                yyaxis right; plot(cumsum(contribution),'linewidth',2,'Color','k'); ylabel('cumulative %var explained');
                drawline(95,'dir','horz','linestyle',':');
                
                % plot mds projection
                subplot(1,2,2);
                pp1_imana('MISC_scatterplotMDS',Y{1}(1:end,1:3),'split',split,'label',labels');
        end
%          keyboard
        varargout = {contribution};
    case 'ROI_getChordPSC'
        % makes a nearly-universal plotting structure from ROI_stats output
        glm = [];
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
        D=getrow(D,D.numDigits<=5); % drop thumb-press responses
        varargout = {D};
    case 'ROI_getSingleFingerRatio'
        glm = [];
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
    
    case 'ROI_pscScaling'
        % calculate linear and log-linear fits of overall scaling as a
        % function of # fingers stimulated.
        glm = [];
        roi = [];
        sn  = 2:11;
        vararginoptions(varargin,{'roi','glm','sn'});
        R=pp1_imana('ROI_getChordPSC','glm',glm); % get psc per chord
        R=getrow(R,ismember(R.roi,roi) & ismember(R.sn,sn)); % restrict to roi and subjs and remove thumb press activity
        T=[]; % output
        for rr=roi
            for ss=sn
                Rs=getrow(R,R.sn==ss & R.roi==rr);
                Rs=tapply(Rs,{'numDigits'},{'psc','mean'}); % avg. across chords per # digits
                tss = sum((Rs.psc - mean(Rs.psc)).^2);
                b0 = ones(size(Rs.numDigits)); % intercept regressor
                for mm=1:2
                    switch mm
                        case 1 % linear scaling
                            X = [b0 Rs.numDigits];
                        case 2 % log-linear scaling
                            X = [b0 log(Rs.numDigits)];
                    end
                    b = pinv(X)*Rs.psc;
                    res = Rs.psc - X*b;
                    rss = sum(res.^2);
                    t.model=mm;
                    t.r2 = 1-(rss/tss);
                    t.sn = ss;
                    t.roi = rr;
                    t.glm = glm;
                    t.pred = [X*b]';
                    t.true = Rs.psc';
                    T=addstruct(T,t);
                end
            end
        end
        varargout = {T};
        
    case 'toyMFmodel'
        % use subject dataset to do some toy encoding models:
        Y = [];
        chordVec = [];
        partVec = [];
        chords = [];
        load(fullfile(baseDir,'pp1_s09_glm4_roi2.mat'));
        % move into structure:
        T.chordNum = chordVec;
        T.run      = partVec;
        T.betas    = Y;
        T.chords   = chords(T.chordNum,:);
        T.numDigits= sum(T.chords,2);

        % not crossvalidated:
        t = tapply(T,{'chordNum'},{'betas','mean'});
        t.chords = chords;
        t.numDigits = sum(chords,2);
        t.beta_pred_linear = nan(size(t.betas));
        for cc=6:31
            digitsInChords = find(chords(cc,:));
            t.beta_pred_linear(cc,:) = sum(t.betas(digitsInChords,:));
        end
        t = getrow(t,t.numDigits>1); % test fit on multi finger chords only
        t.err = t.betas-t.beta_pred_linear;
        r_overall = corr(t.betas(:),t.beta_pred_linear(:));
        r2_overall = 1-(sum(t.err(:).^2) / sum(t.betas(:).^2));
        % make nonlinear model on non-crossval data:
        opts = optimset('fminsearch');
        opts.MaxIter = 10000;
        theta = [0.5,0.5];
        [theta,err,gEst] = fminsearch(@(theta) fitPowerScaling(theta,t.beta_pred_linear,t.betas),theta,opts);
        signB = sign(t.beta_pred_linear);
        absB  = abs(t.beta_pred_linear);
        t.beta_pred_nl = ((absB.^theta(1)).*signB).*theta(2);
        t.errNL = t.betas-t.beta_pred_nl;
        r_overall_nl = corr(t.betas(:),t.beta_pred_nl(:));
        r2_overall_nl = 1-(sum(t.errNL(:).^2) / sum(t.betas(:).^2));
        
        % do leave-one-out crossvalidation:
        D=[];
        for rr = unique(T.run)'
            % train with single finger data
            trainData = getrow(T,T.run~=rr);
            trainData = tapply(trainData,{'chordNum','numDigits'},{'betas','mean'});
            % test against multi-finger data
            testData  = getrow(T,T.run==rr & T.numDigits>1);
            testData.beta_pred = nan(size(testData.betas));
            % get measure of best-possible fits by calculating r2 across
            % training and test sets for same multi-finger chords:
            testData.beta_pred_Lceil = trainData.betas(t.chordNum>5,:);
            trainData = getrow(trainData,trainData.numDigits==1);
            % do linear model
            for cc=6:31
                digitsInChords = find(chords(cc,:));
                testData.beta_pred_linear(testData.chordNum==cc,:) = sum(trainData.betas(ismember(trainData.chordNum,digitsInChords),:));
            end
            % do nonlienar model:
            [theta,err,gEst] = fminsearch(@(theta) fitPowerScaling(theta,testData.beta_pred_linear,testData.betas),theta,opts);
            signB = sign(testData.beta_pred_linear);
            absB  = abs(testData.beta_pred_linear);
            testData.beta_pred_nl = ((absB.^theta(1)).*signB).*theta(2);
            D = addstruct(D,testData);
        end
        D.err = D.betas - D.beta_pred_linear;
        r_cv = corr(D.betas(:),D.beta_pred_linear(:));
        r2_cv = 1-(sum(D.err(:).^2) / sum(D.betas(:).^2));
        D.errNL = D.betas - D.beta_pred_nl;
        r_cv_nl = corr(D.betas(:),D.beta_pred_nl(:));
        r2_cv_nl = 1-(sum(D.errNL(:).^2) / sum(D.betas(:).^2));
        r_Lceil = corr(D.betas(:),D.beta_pred_Lceil(:));
        D.errLceil = D.betas - D.beta_pred_Lceil;
        r2_Lceil =  1-(sum(D.errCeil(:).^2) / sum(D.betas(:).^2));
        
        keyboard
    
    case 'getBetaHist'
        % makes histogram plot of beta values for chords split by #digits
        sn = pp1_imana('getSubjs');
        roi = 2;
        vararginoptions(varargin,{'roi'});
        glm = 4;
        rmvMean = 0;
        posJustify = 0;
        [Y,~,condVec] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm,...
            'rmvMean',rmvMean,'posJustify',posJustify);
        
        chords = pp1_imana('chords');
        numDigits = sum(chords,2);
        clrs = {[0.30196 0.30196 0.30196] [0.43922 0 0.33333] [0.4549 0.15686 0.6549] [0.81569 0 0.96078] [0.89412 0.56471 1]};
        sty = style.custom(clrs);
        sty.hist.facealpha = 0.2;
        sty.hist.edgealpha = 1;
        D = []; % plotting structure
        warning off
        for ii=1:numel(sn)
            % (1). calculate avg. beta per num digits, using run avg. betas
            % (mlt whitened)
            C0 = pcm_indicatorMatrix('identity',condVec{ii}); % avg. across runs
            CD = pcm_indicatorMatrix('identity',numDigits);  % avg. per #dgts
            U  = [pinv(CD)*pinv(C0)*Y{ii}]';
            numVox = size(U,1);
            v = ones(numVox*5,1);
            d.numDigits = kron([1:5]',ones(numVox,1));
            d.beta = U(:);
            d.sn = v.*sn(ii);
            d.roi = v.*roi;
            d.glm = v.*glm;
            D = addstruct(D,d);
            % plot
            subplot(2,5,ii);
            plt.hist(D.beta,'split',D.numDigits,'subset',D.sn==sn(ii),'style',sty); % plot so that foremost histogram is for single finger betas
            drawline(0,'dir','vert');
            title([subj_name{sn(ii)} ' ' regname{roi}]);
            for jj=1:5
                drawline(mean(D.beta(D.numDigits==jj & D.sn==sn(ii))),'dir','vert','color',clrs{jj},'linewidth',3);
            end
            legend off
        end
        % arrage into friendlier scatterplot structure:
        T=[];
        for ii=sn
            v=ones(length(D.beta(D.sn==ii & D.numDigits==1)),1);
            t.sn=v.*ii;
            t.roi=v.*roi;
            t.glm=v.*glm;
            t.beta1 = D.beta(D.numDigits==1 & D.sn==ii);
            t.beta2 = D.beta(D.numDigits==2 & D.sn==ii);
            t.beta3 = D.beta(D.numDigits==3 & D.sn==ii);
            t.beta4 = D.beta(D.numDigits==4 & D.sn==ii);
            t.beta5 = D.beta(D.numDigits==5 & D.sn==ii);
            T=addstruct(T,t);
        end
        warning on
        varargout = {D,T};    
    case 'ROI_calcBetaVar'
        sn  = pp1_imana('getSubjs');
        roi = [1:4,7,9];
        glm = 4;
        vararginoptions(varargin,{'roi'});
        chords = pp1_imana('chords');
        numDigits = sum(chords,2);
        D=[];
        for rr=roi
            % get region data
            [Y,~,condVec] = pp1_imana('PCM_getData','sn',sn,'roi',rr,'glm',glm,'rmvMean',0,'posJustify',0);
            for ii=1:numel(sn) % loop through subjects
                % (1). calculate avg. beta per num digits, using run avg. betas
                % (mlt whitened)
                C0 = pcm_indicatorMatrix('identity',condVec{ii}); % avg. across runs
                CD = pcm_indicatorMatrix('identity',numDigits);  % avg. per #dgts
                U  = pinv(CD)*pinv(C0)*Y{ii};
                numVox = size(U,2);
                % calcualte variance & mean of betas:
                d.betaVar = diag(U*U'./numVox);
                d.numDigits = [1:5]';
                d.betaMean = sum(U,2)./numVox;
                v = ones(5,1);
                d.numVox = v.*numVox;
                d.sn = v.*sn(ii);
                d.roi = v.*rr;
                d.glm = v.*glm;
                D = addstruct(D,d);
            end
        end
        varargout = {D};
    case 'scatterBetas'
        % scatterplot betas for 1 and 
        roi = [1:4];%,7,9];
        sty = style.custom(plt.helper.get_shades(10,'jet'));
        sty.scatter.sizedata = 5;
        sty.general.linewidth = 3;
        T=[];
        D=[];
        for ii=1:numel(roi)
            % get data
            [~,t]=pp1_imana('getBetaHist','roi',roi(ii));
            T=addstruct(T,t);
            % plot
            subplot(1,numel(roi),ii);
            title(regname{roi(ii)});
            [d.r2,b,t,p] = plt.scatter(t.beta1,t.beta4,'split',t.sn,'style',sty);
            if ii<numel(roi)
                legend off
            end
            xlims = xlim;
            ylims = ylim;
            minL  = min([min(xlims) min(ylims)]);
            maxL  = max([max(xlims) max(ylims)]);
            xlim([minL maxL]);
            ylim([minL maxL]);
            axis square
            drawline(0,'dir','vert');
            drawline(0,'dir','horz');
            xlabel('betas (1 finger)');
            ylabel('betas (4 fingers)');
            % save regression fits
            d.int   = b(1,:);
            d.b     = b(2,:);
            d.int_t = t(1,:);
            d.b_t   = t(2,:);
            d.int_p = p(1,:);
            d.b_p   = p(2,:);
            d.roi   = roi(ii);
            D=addstruct(D,d);
        end
        varargout = {T,D};
    
    case 'ROI_estSNR'
        % estimates SNR per subj for specified roi
        % signal var = mean covariance across runs (pattern matrix is
        % vectorized per run)
        % error var = estimated err var accounting for signal var and
        % correlation across runs
        I = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = [1:6];
        vararginoptions(varargin,{'glm','roi','sn'});
        conds = 1:5;
        D=[];
        for rr = roi
            % load data
            T=load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,rr));
            for ii=1:size(T.sn,1)
                % get condition specific betas
                b          = [];
                b.beta     = T.betaW{ii};%T.raw_beta{ii}; %
                b.tt       = T.tt{ii};
                b.run      = T.run{ii};
                b          = getrow(b,ismember(b.tt,conds));
                % remove avg. pattern per run (avg. across only the single
                % finger conditions)
                C0 = indicatorMatrix('identity',b.run);
                b.beta = b.beta -C0*pinv(C0)*b.beta;

                % estimate signal and noise variances (after run mean removal
                % to ensure signal and variance components scale similarly)
                [d.evar,d.svar] = pp1_imana('SFT:estimateVariances',b.beta,b.tt,b.run);
                if d.evar<0; d.evar=0; end
                if d.svar<0; d.svar=0; end
                d.roi = T.roi(ii);
                d.conds = conds;
                d.sn = T.sn(ii);
                d.glm = glm;
                D=addstruct(D,d);
            end
        end
        D.snr = D.svar./D.evar;
        varargout = {D};
        
    case '0' % ------------ SFT: single-finger tuning analyses. -----------
    case 'SFT:calcROI'
        % wrapper to perform SFT analysis, using all voxels
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = 2;
        vararginoptions(varargin,{'glm','roi','sn'});
        
        numSim = 10; % # simulated datasets per model per participant (mvnrnd and sparse tuning models)
        conds  = 1:5; % conditions to analyze (single finger conditions)
        numConds = numel(conds);
        % load data
        T=[];
        for rr=roi
            t = load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T=addstruct(T,t);
        end
        
        fprintf('\nsubj\troi\tsvar\t\tevar\t\ts/e ratio');
        fprintf('\n----\t---\t----\t\t----\t\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.betaUW{ii};%T.raw_beta{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            
            % remove avg. pattern per run (avg. across only the single
            % finger conditions)
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            
            % some simulation params needed:
            numVox  = ceil(size(b.beta,2)/numConds)*numConds; % round up so equal # of voxels per condition (for sparse patterns)
            numRun  = numel(unique(T.run{ii}));
            
            % estimate signal and noise variances (after run mean removal
            % to ensure signal and variance components scale similarly)
            % These estimates will be used to simulate mvnrnd and sparsely
            % tuned datasets as comparison models.
            [evar,svar] = pp1_imana('SFT:estimateVariances',b.beta,b.tt,b.run);
            
            % for noisy data the signal variance can become negative. 
            % It's a bummer. When we encounter this, we set signal to 0.
            if svar<0; svar = 0; end
            
            b    = tapply(b,{'tt'},{'beta','mean'}); % avg. betas across runs (mean beta per voxel = 0)
           % Gtt  = cov(b.beta'); % single-finger G
           % H = eye(numConds)-ones(numConds)./numConds;
            Gtt  = eye(numConds)-ones(numConds)./numConds;
            % 1. calc tuning of actual voxel data
            sftBeta = pp1_imana('SFT:estimateSFT',b.beta);
            sftBeta = mean(sftBeta);
            
            % 2. calc expected tuning of voxels with ~N(0,G)
            [sftGaussEV,evarEst_gauss,svarEst_gauss] = pp1_imana('SFT:expectedValue_G',evar,svar,Gtt,numVox,numRun,numSim); % expected value of the null
            
            % 3. calc expected tuning of voxels with ~sparse tuning, but where signal-to-noise equals that voxel data
            [sftSparseEV,evarEst_sp,svarEst_sp] = pp1_imana('SFT:expectedValue_Sparse',evar,svar,numConds,numVox,numRun,numSim); % 1= tuned to one condition (perfectly sparse)
            %[sftSp2,sftDistSp2] = pp1_imana('SFT:calcExpectedValueSparse',evar,svar,2,numVox,numRun,numSim); % 2= tuned to two conditions
            
            % 4. do prob test for each subject on their gauss sft distribution
            pEV = sum(sftBeta<=sftGaussEV)/length(sftGaussEV);
            if isempty(pEV); pEV = realmin; end
            
            % 5. do prob test for each subject on their sparse sft
            % distributions
            pSp = sum(sftBeta<=sftSparseEV)/length(sftSparseEV);
            if isempty(pSp); pSp = realmin; end
            
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.numVoxTot = v.*size(T.raw_beta{ii},2);
            d.glm       = v.*glm;
            d.sft       = [sftBeta; mean(sftGaussEV); mean(sftSparseEV)];
            d.sftProb   = [0; pEV; pSp];
            d.svar_est  = [svar; mean(svarEst_gauss); mean(svarEst_sp)];
            d.evar_est  = [evar; mean(evarEst_gauss); mean(evarEst_sp)];
            d.snr       = d.svar_est ./ d.evar_est;
            d.isEV      = [0; 1; 2]; 
            D = addstruct(D,d);
            % display to user:
            fprintf('%s\t%02d\t%2.3f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),svar,evar,svar/evar);
        end
        save(fullfile(regDir,sprintf('sft_glm%d',glm)),'-struct','D');
        varargout = {D};   
    case 'SFT:calcROI_fthres'
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = [1:4];
        vararginoptions(varargin,{'glm','roi','sn'});
        
        numSim = 1000; % # simulated datasets per model per participant (random and sparse tuning models)
        conds  = 1:5; % conditions to analyze (single finger conditions)
        numConds = numel(conds);
        % load data
        T=[];
        for rr=roi
            t = load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T=addstruct(T,t);
        end
        T=getrow(T,ismember(T.sn,sn));
        
        fprintf('\nsubj\troi\tsig voxels (%%)\t\tsvar\t\tevar\t\ts/e ratio');
        fprintf('\n----\t---\t--------------\t\t----\t\t----\t\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.betaUW{ii};%T.raw_beta{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            
            % remove avg. pattern per run (avg. across only the single
            % finger conditions)
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            
            % restrict analyses to voxels with significant F-stats:
            [F,Fcrit]  = pp1_imana('SFT:calcFstat',b.beta,b.tt,b.run);
            numVoxOrig = size(b.beta,2);
            numVoxSig  = sum(F>=Fcrit);
            if numVoxSig==0 % if no voxels meet significance threshold, boot
                d.passive   = v;
                d.sn        = v.*T.sn(ii);
                d.roi       = v.*T.roi(ii);
                d.numVoxF   = v.*numVoxSig;
                d.numVoxTot = v.*size(T.raw_beta{ii},2);
                d.avgF      = v.*nan;
                d.glm       = v.*glm;
                d.sft       = [nan;nan;nan];
                d.sftProb   = [nan;nan;nan];
                d.isEV      = [0;1;2];
                D = addstruct(D,d);
                % display to user:
                fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,nan,nan,nan);
                continue
            end
            % some simulation params needed:
            numVoxSim  = ceil(numVoxSig/numConds)*numConds; % round up so equal # of voxels per condition (for sparse patterns)
            %numVoxSim = 100;
            numRun  = numel(unique(T.run{ii}));
            
            % remove avg. pattern per run (avg. across only the single
            % finger conditions) from the significant voxels
            b.beta      = b.beta(:,F>=Fcrit); % drop non-sig. voxels
            C0          = indicatorMatrix('identity',b.run);
            b.beta      = b.beta -C0*pinv(C0)*b.beta; % remove run means
            [evar,svar] = pp1_imana('SFT:estimateVariances',b.beta,b.tt,b.run);
            % estimate signal and noise variances (after run mean removal
            % to ensure signal and variance components scale similarly)
            % These estimates will be used to simulate mvnrnd and sparsely
            % tuned datasets as comparison models.
            
            % for noisy data the signal variance can become negative. 
            % It's a bummer. When we encounter this, we set signal to 0.
            if svar<0; svar = 0; end
            
            b    = tapply(b,{'tt'},{'beta','mean'}); % avg. betas across runs (mean beta per voxel = 0)
           % Gtt  = cov(b.beta'); % single-finger G
           % H = eye(numConds)-ones(numConds)./numConds;
            Gtt  = eye(numConds);
            % 1. calc tuning of actual voxel data
            sftBeta = pp1_imana('SFT:estimateSFT',b.beta);
            sftBeta = mean(sftBeta);
            
            % 2. calc expected tuning of voxels with ~N(0,G)
            [sftGaussEV,evarEst_gauss,svarEst_gauss] = pp1_imana('SFT:expectedValue_G',evar,svar,Gtt,numVoxSim,numRun,numSim); % expected value under random tuning (null)
            
            % 3. calc expected tuning of voxels with ~sparse tuning, but where signal-to-noise equals that voxel data
            [sftSparseEV,evarEst_sp,svarEst_sp] = pp1_imana('SFT:expectedValue_Sparse',evar,svar,numConds,numVoxSim,numRun,numSim); % 1= tuned to one condition (perfectly sparse)
            %[sftSp2,sftDistSp2] = pp1_imana('SFT:calcExpectedValueSparse',evar,svar,2,numVox,numRun,numSim); % 2= tuned to two conditions
            
            % 4. do prob test for each subject on their gauss sft distribution
            pEV = sum(sftBeta<=sftGaussEV)/length(sftGaussEV);
            if isempty(pEV); pEV = realmin; end
            
            % 5. do prob test for each subject on their sparse sft
            % distributions
            pSp = sum(sftBeta<=sftSparseEV)/length(sftSparseEV);
            if isempty(pSp); pSp = realmin; end
            
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.numVoxF   = v.*numVoxSig;
            d.numVoxTot = v.*size(T.raw_beta{ii},2);
            d.avgF      = v.*mean(F(F>=Fcrit));
            d.glm       = v.*glm;
            d.sft       = [sftBeta; mean(sftGaussEV); mean(sftSparseEV)];
            d.sftProb   = [0; pEV; pSp];
            d.svar_est  = [svar; mean(svarEst_gauss); mean(svarEst_sp)];
            d.evar_est  = [evar; mean(evarEst_gauss); mean(evarEst_sp)];
            d.snr       = d.svar_est ./ d.evar_est;
            d.isEV      = [0; 1; 2]; 
            D = addstruct(D,d);
            % display to user:
            fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,svar,evar,svar/evar);
        end
        %save(fullfile(regDir,sprintf('sft_glm%d_fthres',glm)),'-struct','D');
        varargout = {D};   
    case 'SFT:new'
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = [1:4];
        vararginoptions(varargin,{'glm','roi','sn'});
        
        numSim = 1000; % # simulated datasets per model per participant (mvnrnd and sparse tuning models)
        conds  = 1:5; % conditions to analyze (single finger conditions)
        numConds = numel(conds);
        % load data
        T=[];
        for rr=roi
            t = load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T=addstruct(T,t);
        end
        T=getrow(T,ismember(T.sn,sn));
        
        fprintf('\nsubj\troi\tsig voxels (%%)\t\tsvar\t\tevar\t\ts/e ratio');
        fprintf('\n----\t---\t--------------\t\t----\t\t----\t\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.betaUW{ii};%T.raw_beta{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            
            % estimate signal and noise strengths
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            [evar,svar] = pp1_imana('SFT:estimateVariances',b.beta,b.tt,b.run);
            if svar<0; svar = 0; end
            
            % restrict analyses to voxels with significant F-stats:
            [F,Fcrit]  = pp1_imana('SFT:calcFstat',b.beta,b.tt,b.run);
            numVoxOrig = size(b.beta,2);
            numVoxSig  = sum(F>=Fcrit);
            b.beta      = b.beta(:,F>=Fcrit); % drop non-sig. voxels
            
            if numVoxSig==0 % if no voxels meet significance threshold, boot
                d.passive   = v;
                d.sn        = v.*T.sn(ii);
                d.roi       = v.*T.roi(ii);
                d.numVoxF   = v.*numVoxSig;
                d.numVoxTot = v.*size(T.raw_beta{ii},2);
                d.avgF      = v.*nan;
                d.glm       = v.*glm;
                d.sft       = [nan;nan;nan];
                d.sftProb   = [nan;nan;nan];
                d.isEV      = [0;1;2];
                D = addstruct(D,d);
                % display to user:
                fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,nan,nan,nan);
                continue
            end
            
            b    = tapply(b,{'tt'},{'beta','mean'}); % avg. betas across runs (mean beta per voxel = 0)
            %Gtt  = eye(numConds);
            G = pp1_encoding('getRegionG','roi',roi,'glm',glm);
            Gtt = G(1:5,1:5);
            % 1. calc tuning of actual voxel data
            sftBeta = pp1_imana('SFT:estimateSFT',b.beta);
            sftBeta = mean(sftBeta);
            % random tuning
            numRun = numel(unique(T.run{ii}));
            [sftGaussEV,evarEst_gauss,svarEst_gauss] = pp1_imana('SFT:expectedValue_G_fthres',evar,svar,Gtt,numVoxOrig,numRun,numSim);
            % sparse tuning
            [sftSparseEV,evarEst_sp,svarEst_sp] = pp1_imana('SFT:expectedValue_Sparse_fthres',evar,svar,numConds,numVoxOrig,numRun,numSim);
            
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.numVoxF   = v.*numVoxSig;
            d.numVoxTot = v.*size(T.raw_beta{ii},2);
            d.avgF      = v.*mean(F(F>=Fcrit));
            d.glm       = v.*glm;
            d.sft       = [sftBeta; mean(sftGaussEV); mean(sftSparseEV)];
            d.svar_est  = [svar; mean(svarEst_gauss); mean(svarEst_sp)];
            d.evar_est  = [evar; mean(evarEst_gauss); mean(evarEst_sp)];
            d.snr       = d.svar_est ./ d.evar_est;
            d.isEV      = [0; 1; 2]; 
            d.model     = [1;2;3];
            D = addstruct(D,d);
            % display to user:
            fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,svar,evar,svar/evar);
        end
        %save(fullfile(regDir,sprintf('sft_glm%d_fthres_new',glm)),'-struct','D');
        varargout = {D};   
    case 'SFT:new2'
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        % here, simulated data is also F-thresholded.
        
        I   = pp1_imana('LIST_subjs');
        sn  = [2,6];%pp1_imana('getSubjs');
        glm = 4;
        roi = [1];
        vararginoptions(varargin,{'glm','roi','sn'});
        
        numSim = 1000; % # simulated datasets per model per participant (mvnrnd and sparse tuning models)
        conds  = 1:5; % conditions to analyze (single finger conditions)
        numConds = numel(conds);
        % load data
        T=[];
        for rr=roi
            t = load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T=addstruct(T,t);
        end
        T=getrow(T,ismember(T.sn,sn));
        
        fprintf('\nsubj\troi\tsig voxels (%%)\t\tsvar\t\tevar\t\ts/e ratio');
        fprintf('\n----\t---\t--------------\t\t----\t\t----\t\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.betaUW{ii};%T.raw_beta{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            
            % remove run means
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
                        
            % restrict analyses to voxels with significant F-stats:
            [F,Fcrit]  = pp1_imana('SFT:calcFstat',b.beta,b.tt,b.run);
            numVoxOrig = size(b.beta,2);
            numVoxSig  = sum(F>=Fcrit);
            b.beta     = b.beta(:,F>=Fcrit); % drop non-sig. voxels
            
            if numVoxSig==0 % if no voxels meet significance threshold, boot
                d.passive   = v;
                d.sn        = v.*T.sn(ii);
                d.roi       = v.*T.roi(ii);
                d.numVoxF   = v.*numVoxSig;
                d.numVoxTot = v.*size(T.raw_beta{ii},2);
                d.avgF      = v.*nan;
                d.glm       = v.*glm;
                d.sft       = [nan;nan;nan];
                d.sftProb   = [nan;nan;nan];
                d.isEV      = [0;1;2];
                D = addstruct(D,d);
                % display to user:
                fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,nan,nan,nan);
                continue
            end
            
            % estimate signal and noise strengths
            [evar,svar] = pp1_imana('SFT:estimateVariances',b.beta,b.tt,b.run);
            if svar<0; svar = 0; end % rarely, if ever, happens (do this to avoid complex values)
            if evar<0; svar = 0; end % never happens
            
            b    = tapply(b,{'tt'},{'beta','mean'}); % avg. betas across runs (mean beta per voxel = 0)
            %Gtt  = eye(numConds);
            G = pp1_encoding('getRegionG','roi',T.roi(ii),'glm',glm);
            Gtt = G(1:5,1:5);
            % 1. calc tuning of actual voxel data
            sftBeta = pp1_imana('SFT:estimateSFT',b.beta);
            sftBeta = mean(sftBeta);
            % random tuning
            numRun = numel(unique(T.run{ii}));
            [sftGaussEV,evarEst_gauss,svarEst_gauss,numSig_gauss] = pp1_imana('SFT:expectedValue_G_fthres',evar,svar,Gtt,numVoxSig,numRun,numSim);
            % sparse tuning
            [sftSparseEV,evarEst_sp,svarEst_sp,numSig_sparse] = pp1_imana('SFT:expectedValue_Sparse_fthres',evar,svar,numConds,numVoxSig,numRun,numSim);
            
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.numVoxF   = [numVoxSig;numSig_gauss;numSig_sparse];
            d.numVoxTot = v.*size(T.raw_beta{ii},2);
            d.avgF      = v.*mean(F(F>=Fcrit));
            d.glm       = v.*glm;
            d.sft       = [sftBeta; nanmean(sftGaussEV); nanmean(sftSparseEV)];
            d.svar_est  = [svar; nanmean(svarEst_gauss); nanmean(svarEst_sp)];
            d.evar_est  = [evar; nanmean(evarEst_gauss); nanmean(evarEst_sp)];
            d.snr       = d.svar_est ./ d.evar_est;
            d.isEV      = [0; 1; 2]; 
            d.model     = [1;2;3];
            D = addstruct(D,d);
            % display to user:
            fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,svar,evar,svar/evar);
        end
        % normalize sfts:
        D.sft_norm = D.sft - kron(D.sft(D.model==2),ones(3,1)); % 0=random tuning
        D.sft_norm = D.sft_norm./kron(D.sft_norm(D.model==3),ones(3,1)); % 1=selective tuning
        
        %save(fullfile(regDir,sprintf('sft_glm%d_fthres_NEW2.mat',glm)),'-struct','D');
        varargout = {D};   
    case 'SFT:perVox'
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        I   = pp1_imana('LIST_subjs');
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = [1:4];
        vararginoptions(varargin,{'glm','roi','sn'});
        
        numSim = 100; % # simulated datasets per model per participant (mvnrnd and sparse tuning models)
        conds  = 1:5; % conditions to analyze (single finger conditions)
        numConds = numel(conds);
        % load data
        T=[];
        for rr=roi
            t = load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T=addstruct(T,t);
        end
        T=getrow(T,ismember(T.sn,sn));
        
        fprintf('\nsubj\troi\tsig voxels (%%)\t\tsvar\t\tevar\t\ts/e ratio');
        fprintf('\n----\t---\t--------------\t\t----\t\t----\t\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.betaUW{ii};%T.raw_beta{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            
            % estimate signal and noise strengths
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            [evar,svar,voxr] = pp1_imana('SFT:estimateVariances_perVox',b.beta,b.tt,b.run);
            
            % restrict analyses to voxels with significant F-stats:
            [F,Fcrit]  = pp1_imana('SFT:calcFstat',b.beta,b.tt,b.run);
            sigIdx = F>=Fcrit;
            q = ones(size(svar));
            q(sigIdx) = 2;
%             evar(evar<0)=0;
%             svar(svar<0)=0;
%             phat_e = gamfit(evar);
%             phat_s = gamfit(svar);
            
            keyboard
            
            % ************ not done
            
            numVoxOrig = size(b.beta,2);
            numVoxSig  = sum(F>=Fcrit);
            b.beta     = b.beta(:,F>=Fcrit); % drop non-sig. voxels
            
            if numVoxSig==0 % if no voxels meet significance threshold, boot
                d.passive   = v;
                d.sn        = v.*T.sn(ii);
                d.roi       = v.*T.roi(ii);
                d.numVoxF   = v.*numVoxSig;
                d.numVoxTot = v.*size(T.raw_beta{ii},2);
                d.avgF      = v.*nan;
                d.glm       = v.*glm;
                d.sft       = [nan;nan;nan];
                d.sftProb   = [nan;nan;nan];
                d.isEV      = [0;1;2];
                D = addstruct(D,d);
                % display to user:
                fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,nan,nan,nan);
                continue
            end
            
            b    = tapply(b,{'tt'},{'beta','mean'}); % avg. betas across runs (mean beta per voxel = 0)
            %Gtt  = eye(numConds);
            G = pp1_encoding('getRegionG','roi',roi,'glm',glm);
            Gtt = G(1:5,1:5);
            % 1. calc tuning of actual voxel data
            sftBeta = pp1_imana('SFT:estimateSFT',b.beta);
            sftBeta = mean(sftBeta);
            % random tuning
            numRun = numel(unique(T.run{ii}));
            [sftGaussEV,evarEst_gauss,svarEst_gauss] = pp1_imana('SFT:expectedValue_G_fthres',evar,svar,Gtt,numVoxOrig,numRun,numSim);
            % sparse tuning
            [sftSparseEV,evarEst_sp,svarEst_sp] = pp1_imana('SFT:expectedValue_Sparse_fthres',evar,svar,numConds,numVoxOrig,numRun,numSim);
            
            % add to output structure
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.numVoxF   = v.*numVoxSig;
            d.numVoxTot = v.*size(T.raw_beta{ii},2);
            d.avgF      = v.*mean(F(F>=Fcrit));
            d.glm       = v.*glm;
            d.sft       = [sftBeta; mean(sftGaussEV); mean(sftSparseEV)];
            d.svar_est  = [svar; mean(svarEst_gauss); mean(svarEst_sp)];
            d.evar_est  = [evar; mean(evarEst_gauss); mean(evarEst_sp)];
            d.snr       = d.svar_est ./ d.evar_est;
            d.isEV      = [0; 1; 2]; 
            d.model     = [1;2;3];
            D = addstruct(D,d);
            % display to user:
            fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,svar,evar,svar/evar);
        end
        %save(fullfile(regDir,sprintf('sft_glm%d_fthres_new',glm)),'-struct','D');
        varargout = {D};   
   
        
    case 'SFT:estimateSFT'
        % calculate single finger tuning using normalzied distance approach
        X         = varargin{1}; % CxN matrix of data.
        numC      = size(X,1);
        maxX      = max(X,[],1);
        avgDistsN = (sum(maxX-X)./(numC-1)) ./ (maxX-min(X,[],1));
        varargout = {avgDistsN};
    case 'SFT:estimateVariances'
        % empirically estimates error variance in activity patterns across
        % runs. 
        % Estimate is more accurate when run means have been removed.
        
        % % NOTE: we integrate across conditions & voxels
        
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
        
        varargout = {var_E,var_S,r};
    case 'SFT:estimateVariances_perVox'
        % empirically estimates error variance in activity patterns across
        % runs. 
        % Estimate is more accurate when run means have been removed.
        
        % % NOTE: we integrate across conditions for each voxel
        
        Y = varargin{1};    % patterns   [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        p = varargin{3};    % partitions [regressors x 1]
        
        nV = size(Y,2);         % # voxels
        nP = numel(unique(p));  % # partitions
        take = logical(tril(ones(nP),-1)); % lower-triangular index
        
        % remove run means:
        C0 = indicatorMatrix('identity',p);
        Y  = Y - C0*pinv(C0)*Y;
        
        for vv = 1:nV
            ytc = pivottable(c,p,Y(:,vv),'mean'); % [condition x partition]
            G   = cov(ytc);  % covariances between runs (each row = one run, has zero mean)
            var_s(vv)= sum(G(take))/sum(sum(take)); % avg. crossvalidated signal co-variance
            
            R = corr(ytc); % correlations between runs
            r(vv) = sum(R(take))/sum(sum(take)); % avg. crossvalidated correlation
            var_e(vv) = var_s(vv)/r(vv) - var_s(vv);
        end
        
        varargout = {var_e,var_s,r};
    case 'getTuningCurves'
        % wrapper to perform SFT analysis, specifically focused on voxels
        % that show significant finger resposnes (voxel F-stat > F-crit)
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = [1:6];
        vararginoptions(varargin,{'glm','roi','sn'});
        
        conds  = 1:5; % conditions to analyze (single finger conditions)
        % load data
        T=[];
        for rr=roi
            t = load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,rr)));
            T=addstruct(T,t);
        end
        T=getrow(T,ismember(T.sn,sn));
        
        D = []; % output structure
        Q = [];
        v = ones(5,1); % helper vector
        for ii = 1:size(T.sn,1)
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.betaUW{ii};%T.raw_beta{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            
            % estimate signal and noise strengths
            C0 = indicatorMatrix('identity',b.run);
            b.beta = b.beta -C0*pinv(C0)*b.beta;
            [evar,svar] = pp1_imana('SFT:estimateVariances',b.beta,b.tt,b.run);
            if svar<0; svar = 0; end
            if evar<0; evar = 0; end
            % restrict analyses to voxels with significant F-stats:
            [F,Fcrit]  = pp1_imana('SFT:calcFstat',b.beta,b.tt,b.run);
            sigIdx  = F>=Fcrit;
            if sum(sigIdx)==0
                keyboard
            end
            b.beta  = b.beta(:,sigIdx); % drop non-sig. voxels
            b = tapply(b,{'tt'},{'beta','mean'}); % avg. betas across runs (mean beta per voxel = 0)
            
            sftBeta = pp1_imana('SFT:estimateSFT',b.beta);
            
            % normalize values by lowest and max response
            [maxVal,maxIdx] = max(b.beta);
            [minVal,minIdx] = min(b.beta);
            b.betaNorm = b.beta-minVal;
            maxNormVal = max(b.betaNorm);
            b.betaNorm = b.betaNorm./maxNormVal;
            Cm = indicatorMatrix('identity',maxIdx');
            
            % check if we are missing any fingers:
            isMissing = ~ismember([1:5],unique(maxIdx));
            if sum(isMissing)>0
                Cm_old = Cm;
                Cm = zeros(size(Cm_old,1),5);
                Cm(:,~isMissing) = Cm_old;
            end
            d=[];
            d.tuningNorm = pinv(Cm)*b.betaNorm'; % each row is one finger
            d.tuningRaw  = pinv(Cm)*b.beta'; % each row is one finger
            d.avgOtherFingers = pinv(Cm) * ((sum(b.betaNorm)-1)./3)';
            d.avgSelectivity = pinv(Cm)*sftBeta';
            d.tuningNorm(isMissing,:) = nan;
            d.tuningRaw(isMissing,:) = nan;
            d.avgOtherFingers(isMissing,:) = nan;
            d.numVox     = sum(Cm)'; % # voxels averaged into each finger
            d.numVoxRoi  = v.*size(b.beta,2);
            d.svar  = v.*svar;
            d.evar  = v.*evar;
            d.digitMax = [1:5]';
            d.roi   = v.*T.roi(ii);
            d.sn    = v.*T.sn(ii);
            d.glm   = v.*glm;
            D = addstruct(D,d);
            
            % simpler structure:
            q.roi   = T.roi(ii);
            q.sn    = T.sn(ii);
            q.glm   = glm;
            q.numVoxF = sum(sigIdx);
            q.numVoxTot = size(b.beta,2);
            q.svar  = svar;
            q.evar  = evar;
            q.avgSFT = nanmean(pp1_imana('SFT:estimateSFT',b.betaNorm)); % avg. response to "less-preferred" fingers
            Q=addstruct(Q,q);
            
            
        end
        
        % make into more universal plotting structure:
        T=[];
        for ii=1:size(D.sn)
            t=[];
            t.tuningNorm = D.tuningNorm(ii,:)';
            t.tuningRaw  = D.tuningRaw(ii,:)';
            t.avgSelectivity = v.*D.avgSelectivity(ii);
            t.numVox     = v.*D.numVox(ii);
            t.numVoxRoi  = v.*D.numVoxRoi(ii);
            t.svar  = v.*D.svar(ii);
            t.evar  = v.*D.evar(ii);
            t.digitMax = v.*D.digitMax(ii);
            t.digit    = [1:5]';
            t.roi      = v.*D.roi(ii);
            t.glm      = v.*D.glm(ii);
            t.sn       = v.*D.sn(ii);
            T=addstruct(T,t);
        end

        varargout = {T,D,Q};      
    
    case 'SFT:expectedValue_G'
        % calculates expected value (avg. sft across voxels) for voxels generated under specfified G
        
        % simulation params
        errvar  = varargin{1};
        sigvar  = varargin{2};
        G       = varargin{3};
        numVox  = varargin{4};
        numRun  = varargin{5};
        numSim  = varargin{6};

        % generate data
        D = pp1_imana('SFT:model_Gnoise',G,numVox,numRun,numSim,sigvar,errvar);
        % calc expected tuning on simulated datasets
        v = nan(1,numSim);
        sft = v;
        evar_est = v;
        svar_est = v; % preallocate
        for s = 1:numSim
            d      = getrow(D,D.sn==s);
            % remove run means and estimate simulated err and sig vars
            C0     = indicatorMatrix('identity',d.run); 
            d.Y    = d.Y -C0*pinv(C0)*d.Y; % remove run means
            [evar_est(s),svar_est(s)] = pp1_imana('SFT:estimateVariances',d.Y,d.cond,d.run);
            d = tapply(d,{'cond'},{'Y','mean'}); % avg. betas across runs for this simulated subject
            tmp_sft = pp1_imana('SFT:estimateSFT',d.Y);
            sft(s) = mean(tmp_sft); % mean tuning across voxels for this simulated dataset
        end
        varargout = {sft,mean(evar_est),mean(svar_est)};
    case 'SFT:expectedValue_Sparse'
        % calculates expected value (avg. sft across voxels) for voxels
        % with sparse tuning but corrupted by some degree of noise
        
        % simulation params
        errvar  = varargin{1};
        sigvar  = varargin{2};
        numCond = varargin{3};
        numVox  = varargin{4};
        numRun  = varargin{5};
        numSim  = varargin{6};

        % generate data
%         sparsity = 1;
%         D = pp1_imana('SFT:model_Sparse_OLD',sparsity,numVox,numRun,numSim,sigvar,errvar); % sparse patterns with noise
        D = pp1_imana('SFT:model_Sparse',numCond,numVox,numRun,numSim,sigvar,errvar); % sparse patterns with noise
        % calc expected tuning on simulated datasets
        v = nan(1,numSim);
        sft = v;
        evar_est = v;
        svar_est = v; % preallocate
        for s = 1:numSim
            d      = getrow(D,D.sn==s);
            % remove run means and estimate simulated err and sig vars
            C0     = indicatorMatrix('identity',d.run); 
            d.Y    = d.Y -C0*pinv(C0)*d.Y; % remove run means
            [evar_est(s),svar_est(s)] = pp1_imana('SFT:estimateVariances',d.Y,d.cond,d.run);
            d = tapply(d,{'cond'},{'Y','mean'}); % avg. betas across runs for this simulated subject
            tmp_sft = pp1_imana('SFT:estimateSFT',d.Y);
            sft(s) = mean(tmp_sft); % mean tuning across voxels for this simulated dataset
        end
        varargout = {sft,mean(evar_est),mean(svar_est)};
    case 'SFT:expectedValue_G_fthres'
        % calculates expected value (avg. sft across voxels) for voxels generated under specfified G
        
        % simulation params
        errvar  = varargin{1};
        sigvar  = varargin{2};
        G       = varargin{3};
        numVox  = varargin{4};
        numRun  = varargin{5};
        numSim  = varargin{6};

        % generate data
        D = pp1_imana('SFT:model_Gnoise',G,numVox,numRun,numSim,sigvar,errvar);
        % calc expected tuning on simulated datasets
        v = nan(1,numSim);
        sft = v;
        evar_est = v;
        svar_est = v; % preallocate
        numSig   = v;
        for s = 1:numSim
            d      = getrow(D,D.sn==s);
            Cd = indicatorMatrix('identity',d.cond);
            % remove run means and estimate simulated err and sig vars
            C0     = indicatorMatrix('identity',d.run); 
            d.Y    = d.Y -C0*pinv(C0)*d.Y; % remove run means
            [evar_est(s),svar_est(s)] = pp1_imana('SFT:estimateVariances',d.Y,d.cond,d.run);
            
            [F,Fcrit] = pp1_imana('SFT:calcFstat',d.Y,d.cond,d.run);
            sigIdx = F>=Fcrit;
            Ysig = pinv(Cd)*d.Y(:,sigIdx); % avg. betas for significant voxels across runs for this simulated subject
            numSig(s) = sum(sigIdx);
            
            tmp_sft = pp1_imana('SFT:estimateSFT',Ysig);
            sft(s) = mean(tmp_sft); % mean tuning across voxels for this simulated dataset
        end
        varargout = {sft,nanmean(evar_est),nanmean(svar_est),nanmean(numSig)};
    case 'SFT:expectedValue_Sparse_fthres'
        % calculates expected value (avg. sft across voxels) for voxels
        % with sparse tuning but corrupted by some degree of noise
        
        % simulation params
        errvar  = varargin{1};
        sigvar  = varargin{2};
        numCond = varargin{3};
        numVox  = varargin{4};
        numRun  = varargin{5};
        numSim  = varargin{6};

        % generate data
        D = pp1_imana('SFT:model_Sparse',numCond,numVox,numRun,numSim,sigvar,errvar); % sparse patterns with noise
        % calc expected tuning on simulated datasets
        v = nan(1,numSim);
        sft = v;
        evar_est = v;
        svar_est = v; % preallocate
        for s = 1:numSim
            d      = getrow(D,D.sn==s);
            Cd = indicatorMatrix('identity',d.cond);
            % remove run means and estimate simulated err and sig vars
            C0     = indicatorMatrix('identity',d.run); 
            d.Y    = d.Y -C0*pinv(C0)*d.Y; % remove run means
            [evar_est(s),svar_est(s)] = pp1_imana('SFT:estimateVariances',d.Y,d.cond,d.run);
            
            [F,Fcrit] = pp1_imana('SFT:calcFstat',d.Y,d.cond,d.run);
            sigIdx = F>=Fcrit;
            Ysig = pinv(Cd)*d.Y(:,sigIdx); % avg. betas for significant voxels across runs for this simulated subject
            numSig(s) = sum(sigIdx);
            
            tmp_sft = pp1_imana('SFT:estimateSFT',Ysig);
            sft(s)  = mean(tmp_sft); % mean tuning across significant voxels for this simulated dataset
        end
        varargout = {sft,mean(evar_est),mean(svar_est),mean(numSig)};
    
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
        G       = varargin{1}; % G is already double centred b/c run means removed
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5};
        noise   = varargin{6};
        
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation

        numCond = size(G,1);
        % need to scale all elements of G by mean of the diagonal elements
        % (variances) to ensure appropriate signal scaling:
        G = G./ (sum(diag(G)) / (numCond-1));
        
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
    case 'SFT:model_Gnoise_test'
        % true betas drawn from normal distribution, with added i.i.d. noise
        G       = varargin{1}; % G is already double centred b/c run means removed
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5};
        noise   = varargin{6};
        
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation

        numCond = size(G,1);
        % need to scale all elements of G by mean of the diagonal elements
        % (variances) to ensure appropriate signal scaling:
        G = G./ (sum(diag(G)) / (numCond-1));
        
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
            
            d.Y = X*trueU;
            % indexing fields
            d.run   = kron([1:numRun],ones(1,numCond))';
            d.cond  = kron(ones(1,numRun),[1:numCond])';
            d.sn    = ones(size(d.run)).*s;
            d.noise = ones(size(d.run)).*0;
            D = addstruct(D,d);
            
            % Now add the random noise 
            Noise  = noiseDist(pNoise); 
            Noise  = bsxfun(@times,Noise,sqrt(noise)); % Multiply by noise scaling factor
            d.Y    = X*trueU + Noise; % pull through condition specific patterns and add i.i.d. noise
            % indexing fields
            d.run   = kron([1:numRun],ones(1,numCond))';
            d.cond  = kron(ones(1,numRun),[1:numCond])';
            d.sn    = ones(size(d.run)).*s;
            d.noise = ones(size(d.run)).*1;
            D = addstruct(D,d);
        end
        varargout = {D};
    case 'SFT:model_Sparse'
        % true patterns are sparse (0's and 1's), with iid noise
        numCond = varargin{1};% sparsity level (1=totally sparse, 2-two conditions tuned, etc.)
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5}; % signal variance
        noise   = varargin{6}; % noise variance
        
        % get # conditions based on sparsity level
        numVoxPerCond = ceil(numVox/numCond); % voxels tuned per chord
        numVox = numVoxPerCond*numCond;
        signal = signal*(numCond/1); % rescale the signal by # conditions (each condition contributes independent amount to signal)

        % define signal generation
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        D = []; % output structure
        for s = 1:numSubj % per simulated dataset
            % draw uniform values for signal and noise (to be inverted
            % through any arbitrary function later)
            pNoise = unifrnd(0,1,numCond*numRun,numVox); 
            % Generate true sparse patterns
            U = kron(eye(numCond),ones(1,numVoxPerCond)); % true patterns are 1s for voxels tuned to fingers, 0s for non-tuned
            % scale patterns to match specified signal strength
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
    case 'SFT:model_Sparse_test'
        % true patterns are sparse (0's and 1's), with iid noise
        numCond = varargin{1};% sparsity level (1=totally sparse, 2-two conditions tuned, etc.)
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5}; % signal variance
        noise   = varargin{6}; % noise variance
        
        % get # conditions based on sparsity level
        numVoxPerCond = ceil(numVox/numCond); % voxels tuned per chord
        numVox = numVoxPerCond*numCond;
        signal = signal*(numCond/1); % rescale the signal by # conditions (each condition contributes independent amount to signal)

        % define signal generation
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        D = []; % output structure
        for s = 1:numSubj % per simulated dataset
            % draw uniform values for signal and noise (to be inverted
            % through any arbitrary function later)
            pNoise = unifrnd(0,1,numCond*numRun,numVox); 
            % Generate true sparse patterns
            U = kron(eye(numCond),ones(1,numVoxPerCond)); % true patterns are 1s for voxels tuned to fingers, 0s for non-tuned
            % scale patterns to match specified signal strength
            U = bsxfun(@times,U,sqrt(signal));
            %U = U-mean(U);
            X = kron(ones(numRun,1),eye(numCond)); % design matrix
            trueU = X*U;
            
            d.Y = trueU;
            % indexing fields
            d.run   = kron([1:numRun],ones(1,numCond))';
            d.cond  = kron(ones(1,numRun),[1:numCond])';
            d.sn    = ones(size(d.run)).*s;
            d.noise = ones(size(d.run)).*0;
            D = addstruct(D,d);
            
            % Now add the random noise 
            Noise  = noiseDist(pNoise); 
            Noise  = bsxfun(@times,Noise,sqrt(noise)); % sacle noise
            d.Y    = trueU + Noise;
            % indexing fields
            d.run   = kron([1:numRun],ones(1,numCond))';
            d.cond  = kron(ones(1,numRun),[1:numCond])';
            d.sn    = ones(size(d.run)).*s;
            d.noise = ones(size(d.run)).*1;
            D = addstruct(D,d);
        end
        varargout = {D};
    case 'SFT:model_Sparse_OLD'
        % true patterns are sparse (0's and 1's), with iid noise
        sparsity = varargin{1};% sparsity level (1=totally sparse, 2-two conditions tuned, etc.)
        numVox  = varargin{2};
        numRun  = varargin{3};
        numSubj = varargin{4};
        signal  = varargin{5}; % signal variance
        noise   = varargin{6}; % noise variance
        
        % get # conditions based on sparsity level
        chords = pp1_imana('chords');
        numDigits = sum(chords,2);
        chords = chords(numDigits==sparsity,:);
        [numChord,numCond] = size(chords);
        numVoxPerChord = ceil(numVox/numChord); % voxels tuned per chord
        numVox = numVoxPerChord*numChord;
        signal = signal*(numCond/sparsity); % rescale the signal by # conditions (each condition contributes independent amount to signal)

        % define signal generation
        noiseDist  = @(x) norminv(x,0,1);  % Standard normal inverse for Noise generation 
        D = []; % output structure
        for s = 1:numSubj % per simulated dataset
            % draw uniform values for signal and noise (to be inverted
            % through any arbitrary function later)
            pNoise = unifrnd(0,1,numCond*numRun,numVox); 
            % Generate true sparse patterns
            U = kron(chords',ones(1,numVoxPerChord)); % true patterns are 1s for voxels tuned to fingers, 0s for non-tuned
            % scale patterns to match specified signal strength
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

        
    case 'SFT:calcFstat'
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
            res = bsxfun(@minus,Y(idx,:),muK(ii,:)); % residuals from the grand mean across all observations of this condition
            SSR = SSR + sum(res.^2,1) ./ n(ii); % SSR (across observations) scaled by number of observations (in case # obs differ per condition)
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
    case 'SFT:calcFstatCV'
        
        % * * * REVIEW THIS! NOT SURE IF CORRECT! * * *
        
        error('please review ''calcFstatCV'' before use');
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
        df2 = numCond*numRun - numCond - numRun +1; % account for leave-one-out cv
        % remove run means 
        C0  = indicatorMatrix('identity',runVec);
        Y   = Y - C0*pinv(C0)*Y; 
        
        % calculate CrossValidated F
        A = zeros(numCond,numVox,numRun); % zero-pad
        ApredCV = A;
        for ii = 1:numRun
            A(:,:,ii) = Y(runVec==ii,:);
        end
        muK = mean(A,3);
        res = bsxfun(@minus,A,muK);
        SSR = sum(sum(res.^2,1),3); % common residual covariance
        SSB = sum((muK.^2).*numRun,1);
        F   = (SSB./df1) ./ (SSR./df2);
        for ii = runs
            testRuns = runs~=ii;
            ApredCV(:,:,ii) = mean(A(:,:,testRuns),3);
        end
        % crossval f-stat
        %Apred = repmat(mean(A,3),1,1,numRun); % predicted voxel tunings
        TSS   = sum(sum((A).^2,3),1);         % total SS 
        %RSS   = sum(sum((A-Apred).^2,3),1);   % finger model (5 params)
        RSScv = sum(sum((A-ApredCV).^2,3),1); % unrestricted SSR
        SSBcv = TSS-RSScv;
        FstatCV = (SSBcv./df1) ./ (RSScv./df2);
        FcritCV = finv(0.95,df1,df2); % 95% cutoff for F-stat
        varargout = {FstatCV,FcritCV};
        %keyboard
    
    case 'PLOT_sft'
         % plot single-finger preference analysis results
        sn  = pp1_imana('getSubjs');
        roi = [1:4];
        glm = 4;
        vararginoptions(varargin,{'roi','glm','sn'});
        % load data
        T = load(fullfile(regDir,sprintf('sft_glm%d.mat',glm)));
        T = getrow(T,ismember(T.roi,roi) & ismember(T.sn,sn) & T.isEV<3);
        % get roi labels
        labels1 = {};
        labels2 = {};
        for r=roi
            labels1{end+1} = '';
            labels1{end+1} = regname{r};
            labels1{end+1} = '';
            labels2{end+1} = regname{r};
        end
        % plot styles
        sty1 = style.custom({'blue','darkgray','red','orange'}); 
        sty2 = style.custom({[0 0.49412 0.49804]});
        %sty2 = style.custom({'black'});           
        sty1.general.markersize = 4;
        sty2.general.markersize = 4;
        T.sft = real(T.sft);
        % plot sft and mean (of simulated subject samples) expected value
        subplot(3,4,[1,2,5,6]);
        plt.box(T.roi,T.sft,'split',T.isEV,'style',sty1,'plotall',1);
        set(gca,'xticklabel',labels1,'xticklabelrotation',45);
        ylabel('selectivity index');
        title('single-finger selectivity');
        % plot normalized sft (0=expected noise, 1=fully sparse given
        % estimated noise level)
        Tn = getrow(T,T.isEV==0);
        Tn.sftNorm  = (T.sft(T.isEV==0)-T.sft(T.isEV==1))./(T.sft(T.isEV==2)-T.sft(T.isEV==1));
        Tn.sftRatio = T.sft(T.isEV==0)./T.sft(T.isEV==2);
        
        subplot(3,4,[3,4,7,8]);
        %plt.line(Tn.roi,Tn.sftNorm,'style',sty2,'split',Tn.sn);
        plt.box(Tn.roi,Tn.sftNorm,'style',sty2,'plotall',0);
        set(gca,'xticklabel',labels2,'xticklabelrotation',45);
        ylabel('normalized selectivity');
        title('normalized selectivity index');
        drawline(1,'dir','horz','linestyle','-','linewidth',1);
        drawline(0,'dir','horz','linestyle','-.','linewidth',1);
        
        % plot estimated signal-to-noise variance ratios:
        subplot(3,4,[9,10]);
        plt.box(T.roi,T.svar_est./T.evar_est,'split',T.isEV,'style',sty1,'plotall',1);
        ylabel(sprintf('signal / noise\nvariance estimates'));

        % anova for normalized selectivity index
        anovaMixed(Tn.sftNorm,Tn.sn,'within',Tn.roi,{'roi'});
        % anova to check no differences between model SE ratios
        anovaMixed(T.svar_est./T.evar_est, T.sn,'within',[T.roi T.isEV+1],{'roi','model'});
        % ttests for 3b>1, 3b>2, 1>2
        ttest(Tn.sftNorm(Tn.roi==2),Tn.sftNorm(Tn.roi==3),1,'paired');
        ttest(Tn.sftNorm(Tn.roi==2),Tn.sftNorm(Tn.roi==4),1,'paired');
        ttest(Tn.sftNorm(Tn.roi==3),Tn.sftNorm(Tn.roi==4),1,'paired');
        
        varargout = {T,Tn};
    case 'PLOT_sft_fthres'
        % plot single-finger preference analysis results
        sn  = pp1_imana('getSubjs');
        roi = [1:7];
        glm = 4;
        vararginoptions(varargin,{'roi','glm','sn'});
        % load data
        T = load(fullfile(regDir,sprintf('sft_glm%d_fthres.mat',glm)));
       % T = load(fullfile(regDir,sprintf('sft_glm%d_fthres_allConds.mat',glm)));
        T = getrow(T,ismember(T.roi,roi) & ismember(T.sn,sn) & T.isEV<3);
        % get roi labels
        labels1 = {};
        labels2 = {};
        for r=roi
            labels1{end+1} = '';
            labels1{end+1} = regname{r};
            labels1{end+1} = '';
            labels2{end+1} = regname{r};
        end
        % plot styles
        sty1 = style.custom({'blue','darkgray','red','orange'}); 
        sty2 = style.custom({[0 0.49412 0.49804]});
        %sty2 = style.custom({'black'});           
        sty1.general.markersize = 4;
        sty2.general.markersize = 4;
        T.sft = real(T.sft);
        % plot sft and mean (of simulated subject samples) expected value
        subplot(3,4,[1,2,5,6]);
        plt.dot(T.roi,T.sft,'split',T.isEV,'style',sty1,'plotall',1);
        set(gca,'xticklabel',labels1,'xticklabelrotation',45);
        ylabel('selectivity index');
        title('single-finger selectivity');
        % plot normalized sft (0=expected noise, 1=fully sparse given
        % estimated noise level)
        Tn = getrow(T,T.isEV==0);
        Tn.sftNorm  = (T.sft(T.isEV==0)-T.sft(T.isEV==1))./(T.sft(T.isEV==2)-T.sft(T.isEV==1));
        Tn.sftRatio = T.sft(T.isEV==0)./T.sft(T.isEV==2);
        
        subplot(3,4,[3,4,7,8]);
        %plt.line(Tn.roi,Tn.sftNorm,'style',sty2,'split',Tn.sn);
        plt.box(Tn.roi,Tn.sftNorm,'style',sty2,'plotall',1);
        set(gca,'xticklabel',labels2,'xticklabelrotation',45);
        ylabel('normalized selectivity');
        title('normalized selectivity index');
        drawline(1,'dir','horz','linestyle','-','linewidth',1);
        drawline(0,'dir','horz','linestyle','-.','linewidth',1);
        
        % plot estimated signal-to-noise variance ratios:
        subplot(3,4,[9,10]);
        plt.box(T.roi,T.svar_est./T.evar_est,'split',T.isEV,'style',sty1,'plotall',1);
        ylabel(sprintf('signal / noise\nvariance estimates'));
        
        % plot number of voxels included (i.e. those significant after F-test):
        subplot(3,4,11);
        plt.dot(T.roi,T.numVoxF,'subset',T.isEV==0,'plotall',1);
        set(gca,'xticklabel',labels2,'xticklabelrotation',45);
        %title('# of included voxels');
        ylabel(sprintf('# voxels\nincluded'));
        subplot(3,4,12);
        plt.dot(T.roi,(T.numVoxF./T.numVoxTot).*100,'subset',T.isEV==0,'plotall',1);
        set(gca,'xticklabel',labels2,'xticklabelrotation',45);
        %title('proportion of included voxels');
        ylabel(sprintf('%% voxels\nincluded'));

        % anova for normalized selectivity index
        anovaMixed(Tn.sftNorm,Tn.sn,'within',Tn.roi,{'roi'});
        % anova to check no differences between model SE ratios
%         anovaMixed(T.svar_est./T.evar_est, T.sn,'within',[T.roi T.isEV+1],{'roi','model'});
        % ttests for 3b>3a, 3b>1, 3b>2, 1>2
        fprintf('\n3b >3a\n'); ttest(Tn.sftNorm(Tn.roi==2),Tn.sftNorm(Tn.roi==1),1,'paired');
        fprintf('\n3b >1\n'); ttest(Tn.sftNorm(Tn.roi==2),Tn.sftNorm(Tn.roi==3),1,'paired');
        fprintf('\n3b >2\n'); ttest(Tn.sftNorm(Tn.roi==2),Tn.sftNorm(Tn.roi==4),1,'paired');
        fprintf('\n2 >1\n'); ttest(Tn.sftNorm(Tn.roi==3),Tn.sftNorm(Tn.roi==4),1,'paired');
        
        varargout = {T,Tn};
    case 'PLOT_sft_fthres_lines'
        % plot single-finger preference analysis results
        sn  = pp1_imana('getSubjs');
        roi = [1:4];
        glm = 4;
        vararginoptions(varargin,{'roi','glm','sn'});
        % load data
        T = load(fullfile(regDir,sprintf('sft_glm%d_fthres.mat',glm)));
        T = getrow(T,ismember(T.roi,roi) & ismember(T.sn,sn) & T.isEV<3);
        % get roi labels
        labels1 = {};
        labels2 = {};
        for r=roi
            labels1{end+1} = '';
            labels1{end+1} = regname{r};
            labels1{end+1} = '';
            labels2{end+1} = regname{r};
        end
        % plot styles
        sty1 = style.custom({'blue','darkgray','red','orange'}); 
        sty2 = style.custom({[0 0.49412 0.49804]});
        %sty2 = style.custom({'black'});           
        sty1.general.markersize = 4;
        sty2.general.markersize = 4;
        T.sft = real(T.sft);
        % plot sft and mean (of simulated subject samples) expected value
        subplot(3,4,[1,2,5,6]);
        plt.box(T.roi,T.sft,'split',T.isEV,'style',sty1,'plotall',1);
        set(gca,'xticklabel',labels1,'xticklabelrotation',45);
        ylabel('selectivity index');
        title('single-finger selectivity');
        % plot normalized sft (0=expected noise, 1=fully sparse given
        % estimated noise level)
        Tn = getrow(T,T.isEV==0);
        Tn.sftNorm  = (T.sft(T.isEV==0)-T.sft(T.isEV==1))./(T.sft(T.isEV==2)-T.sft(T.isEV==1));
        Tn.sftRatio = T.sft(T.isEV==0)./T.sft(T.isEV==2);
        
        subplot(3,4,[3,4,7,8]);
        %plt.line(Tn.roi,Tn.sftNorm,'style',sty2,'split',Tn.sn);
        plt.box(Tn.roi,Tn.sftNorm,'style',sty2,'plotall',1);
        set(gca,'xticklabel',labels2,'xticklabelrotation',45);
        ylabel('normalized selectivity');
        title('normalized selectivity index');
        drawline(1,'dir','horz','linestyle','-','linewidth',1);
        drawline(0,'dir','horz','linestyle','-.','linewidth',1);
        
        % plot estimated signal-to-noise variance ratios:
        subplot(3,4,[9,10]);
        plt.box(T.roi,T.svar_est./T.evar_est,'split',T.isEV,'style',sty1,'plotall',1);
        ylabel(sprintf('signal / noise\nvariance estimates'));
        
        % plot number of voxels included (i.e. those significant after F-test):
        subplot(3,4,11);
        plt.box(T.roi,T.numVoxF,'subset',T.isEV==0,'plotall',1);
        set(gca,'xticklabel',labels2,'xticklabelrotation',45);
        %title('# of included voxels');
        ylabel(sprintf('# voxels\nincluded'));
        subplot(3,4,12);
        plt.box(T.roi,(T.numVoxF./T.numVoxTot).*100,'subset',T.isEV==0,'plotall',1);
        set(gca,'xticklabel',labels2,'xticklabelrotation',45);
        %title('proportion of included voxels');
        ylabel(sprintf('%% voxels\nincluded'));

        % anova for normalized selectivity index
        anovaMixed(Tn.sftNorm,Tn.sn,'within',Tn.roi,{'roi'});
        % anova to check no differences between model SE ratios
%         anovaMixed(T.svar_est./T.evar_est, T.sn,'within',[T.roi T.isEV+1],{'roi','model'});
        % ttests for 3b>3a, 3b>1, 3b>2, 1>2
        fprintf('\n3b >3a\n'); ttest(Tn.sftNorm(Tn.roi==2),Tn.sftNorm(Tn.roi==1),1,'paired');
        fprintf('\n3b >1\n'); ttest(Tn.sftNorm(Tn.roi==2),Tn.sftNorm(Tn.roi==3),1,'paired');
        fprintf('\n3b >2\n'); ttest(Tn.sftNorm(Tn.roi==2),Tn.sftNorm(Tn.roi==4),1,'paired');
        fprintf('\n2 >1\n'); ttest(Tn.sftNorm(Tn.roi==3),Tn.sftNorm(Tn.roi==4),1,'paired');
        
        varargout = {T,Tn};   
        
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
        elseif glm>2
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
            % find rows where the contrast is applied to any two chord conditions
            if(sum(any(C(ii,numD~=0 & ismember(numD,numDigits)),1))==2)
                P = [P;ii];
            end
        end
        varargout = {P};
    case 'GET_LDCperNumDigits'    
        sn  = pp1_imana('getSubjs');
        glm = [];
        roi = [];
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
                d.ldc_wmean = T.ldc_wmean(i,P')';
                d.ldc       = T.ldc(i,P')';
                d.corr_dist = T.corr_dist(i,P')';
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
       
    
    case '0' % ------------ TESSEL: do analyses with tessel rois. ----
    case 'TSL:setThres'
        % case to ensure consistent thresholding across TSL cases:
        nTessel = 1442; % options: 162, 362, 642, 1002, 1442
        prop=0.2;
        distThres=0.005;
        varargout = {nTessel,prop,distThres};
    case 'TSL:check'
        % threshold tessels that would be selected based on distance
        % threshold and save as surface label gifti for visual inspection.
        [nTessel,prop,distThres] = pp1_imana('TSL:setThres');
        % get tessels per hemi based on distance threshold criteria:
        for h=1:2
            T = pp1_imana('TSL:select','nTessel',nTessel,'hemi',h,'prop',prop,'distThres',distThres);
            I = gifti(fullfile(atlasDir,sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemLetter{h})));
            I.cdata((~ismember(I.cdata(:,1),T)),1) = nan; % drop tessels not included
            save(I,fullfile(wbDir,sprintf('tessel_check%d.164k.%s.label.gii',nTessel,hemLetter{h})));
        end
        
    case 'TSL:select'
        % here select based on the distance mask
        hemi = [];
        [nTessel,prop,distThres] = pp1_imana('TSL:setThres');
        vararginoptions(varargin,{'hemi','nTessel','distThres','prop'});
        % load in distances and icosahedron
        I = gifti(fullfile(atlasDir,sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemLetter{hemi})));
        G = gifti(fullfile(wbDir,'group_164k',sprintf('summary.LDC.%s.glm4.164k.func.gii',hemLetter{hemi})));
        maskDist = G.cdata(:,1)>distThres;
        tessels = unique(I.cdata)';
        tessels = tessels(tessels~=0); % exclude 0 - medial wall
        choice = [];
        for ii=tessels
          numAll = sum(I.cdata==ii);
          distPres = sum(maskDist(I.cdata==ii)); % how many times significant distance (dist presence)
          if distPres>=(numAll*prop)
             choice = [choice,ii];
          end
        end
        varargout={double(choice)};
    case 'TSL:define'
        % define rois using tesselated surface map
        % restrict tessels to be those where prop% of vertices within tessel
        % have ldc>=distThres
        sn = 2:11;
        vararginoptions(varargin,{'sn','nTessel'});
        nTessel = pp1_imana('TSL:setThres');
        
%         % get tessels per hemi based on distance threshold criteria:
%         for h=1:2
%             T{h} = pp1_imana('TSL:select','hemi',h);
%         end
        % map all tessels, collect beta from all tessels, then threshold
        % later
        
        for s=sn
            fprintf('\n%s.',subj_name{s});
            R=[];
            idx=1;
            mask    = fullfile(glmDir{1},subj_name{s},'mask.nii,1');  % load mask file now 
            subjDir = fullfile(wbDir,subj_name{s});
            for h=1%:2
                I = gifti(fullfile(atlasDir,sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemLetter{h})));
                tessels = unique(I.cdata)';
                tessels = tessels(tessels~=0); % exclude 0 - medial wall                
                for r=tessels%1:numel(T{h})
                    R{idx}.type     = 'surf_nodes'; 
                    R{idx}.location = find(I.cdata(:,1)==r);
                    R{idx}.white    = fullfile(subjDir,sprintf('%s.%s.white.164k.surf.gii',subj_name{s},hemLetter{h}));
                    R{idx}.pial     = fullfile(subjDir,sprintf('%s.%s.pial.164k.surf.gii',subj_name{s},hemLetter{h}));
                    R{idx}.linedef  = [5,0,1];
                    R{idx}.image    = mask;
                    R{idx}.hemi     = h;
                    R{idx}.tesselNum= r;
                    R{idx}.name     = sprintf('%s.tessel-%d',hemLetter{h},r);
                    idx = idx+1;
                end
            end
            R = region_calcregions(R);
            save(fullfile(regDir,sprintf('regions_tessels%d_%s.mat',nTessel,subj_name{s})),'R');
            
        end
    case 'TSL:getBetas'
        % harvest betas from tesselated rois
        I   = pp1_imana('LIST_subjs');
        sn  = 2:11;
        glm = 4;
        nTessel = pp1_imana('TSL:setThres');
        vararginoptions(varargin,{'sn','glm'});
        
        if glm==2
            numRegressors = 31; 
        elseif glm>2
            numRegressors = 32; % 31 chords + thumb response
        end
        fprintf('extracting betas\n');
        % harvest
        for s=sn % for each subj
            T=[];
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            % load files
            load(fullfile(glmDir{glm}, subjName,'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            D = load(fullfile(glmDir{glm}, subjName,'SPM_info.mat'));
            load(fullfile(regDir,sprintf('regions_tessels%d_%s.mat',nTessel,subj_name{s})));          % load subject's region parcellation & depth structure (R)
            % TR img info
            V = SPM.xY.VY; 
            
            C0         = indicatorMatrix('identity',D.run); 
            ofInterest = 1:size(C0,1); % indicies for regressors of interest
            
            for r = 1:numel(R) % for each region on LEFT hemi
                % toss stuff into output structure
                S.sn         = s;
                S.roi        = r;
                S.tesselNum  = R{r}.tesselNum;
                S.hemi       = R{r}.hemi;
                S.tt         = {D.tt};
                S.run        = {D.run};
                S.numDigits  = {D.numDigits};
                S.sess       = {D.sess};
                S.xyzcoord   = {R{r}.data'}; 
                % get raw data/psc for voxels in region
                Y = region_getdata(V,R{r});  % Data Y is N x P (P is in order of transpose of R{r}.depth)
                if ~isempty(Y) % where there betas for this roi?
                    % estimate region betas
                    [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','runwise');
                    % remove nuisance regressor betas
                    betaUW       = bsxfun(@rdivide,beta,sqrt(resMS));
                    betaUW       = betaUW(ofInterest,:);
                    betaW        = betaW(ofInterest,:);
                    raw_beta     = beta(ofInterest,:);
                    % add data to output structure
                    S.betaW              = {betaW};        
                    S.betaUW             = {betaUW};  
                    S.raw_beta           = {raw_beta};
                    S.resMS              = {resMS};
                else % no data for this subj in this tessel
                    S.betaW              = {nan};        
                    S.betaUW             = {nan};  
                    S.raw_beta           = {nan};
                    S.resMS              = {nan};
                end
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
            fprintf('..done\n');
            % save T
            save(fullfile(regDir,sprintf('glm%d_%s_tessels%d_betas.mat',glm,subj_name{s}, nTessel)),'-struct','T'); 
        end
        % save T
       % save(fullfile(regDir,sprintf('glm%d_tessels%d_betas.mat',glm,nTessel)),'-struct','T'); 
        fprintf('\n')
    
    case 'TSL:selectivity'
        % wrapper to perform SFT analysis with data from tessels
        % restricts analysis to voxels that show significant finger resposnes (voxel F-stat > F-crit)
        I   = pp1_imana('LIST_subjs');
        sn  = 2:11;
        glm = 4;
        nTessel = pp1_imana('TSL:setThres');
        vararginoptions(varargin,{'glm','sn'});
        
        numSim = 10; % # simulated datasets per model per participant (mvnrnd and sparse tuning models)
        conds  = 1:5; % conditions to analyze (single finger conditions)
        
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_tessels%d_betas.mat',glm,nTessel)));
        T = getrow(T,ismember(T.sn,sn));

        fprintf('\nsubj\troi\tsig voxels (%%)\t\tsvar\t\tevar\t\ts/e ratio');
        fprintf('\n----\t---\t--------------\t\t----\t\t----\t\t---------\n');
        D = []; % output structure
        v = ones(3,1); % helper vector
        for ii = 1:size(T.sn,1)
            % output structure:
            d.passive   = v;
            d.sn        = v.*T.sn(ii);
            d.roi       = v.*T.roi(ii);
            d.tesselNum = v.*T.tesselNum(ii);
            d.hemi      = v.*T.hemi(ii);
            d.numVoxF   = v.*nan;
            d.numVoxTot = v.*nan;
            d.avgF      = v.*nan;
            d.glm       = v.*glm;
            d.sft       = v.*nan;
            d.sftProb   = v.*nan;
            d.svar_est  = v.*nan;
            d.evar_est  = v.*nan;
            d.snr       = v.*nan;
            d.isEV      = [0; 1; 2];
            if isnan(T.raw_beta{ii})
                D = addstruct(D,d);
                continue; % some participants don't have data from a few tessels
            end
            % for each voxel in subj-roi, get avg. betas for single fingers
            b          = [];
            b.beta     = T.raw_beta{ii};
            b.tt       = T.tt{ii};
            b.run      = T.run{ii};
            b          = getrow(b,ismember(b.tt,conds));
            % restrict analyses to voxels with significant F-stats:
            [F,Fcrit]  = pp1_imana('SFT:calcFstat',b.beta,b.tt,b.run);
            numVoxOrig = size(b.beta,2);
            numVoxSig  = sum(F>=Fcrit);
            if numVoxSig==0 % if no voxels meet significance threshold, boot
                d.numVoxTot = v.*size(T.raw_beta{ii},2);
                D = addstruct(D,d);
                % display to user:
                fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.roi(ii),(numVoxSig/numVoxOrig)*100,nan,nan,nan);
                continue
            end
            % some simulation params needed:
            numVox  = ceil(numVoxSig/5)*5; % round up so equal # of voxels per condition (for sparse patterns)
            numRun  = numel(unique(T.run{ii}));
            
            % remove avg. pattern per run (avg. across only the single
            % finger conditions) from the significant voxels
            b.beta      = b.beta(:,F>=Fcrit); % drop non-sig. voxels
            C0          = indicatorMatrix('identity',b.run);
            b.beta      = b.beta -C0*pinv(C0)*b.beta; % remove run means
            [evar,svar] = pp1_imana('SFT:estimateVariances',b.beta,b.tt,b.run);
            % estimate signal and noise variances (after run mean removal
            % to ensure signal and variance components scale similarly)
            % These estimates will be used to simulate mvnrnd and sparsely
            % tuned datasets as comparison models.
            
            % for noisy data the signal variance can become negative. 
            % It's a bummer. When we encounter this, we set signal to 0.
            if svar<0; svar = 0; end
            
            b    = tapply(b,{'tt'},{'beta','mean'}); % avg. betas across runs (mean beta per voxel = 0)
            Gtt  = cov(b.beta'); % single-finger G
            
            % 1. calc tuning of actual voxel data
            sftBeta = pp1_imana('SFT:estimateSFT',b.beta);
            sftBeta = mean(sftBeta);
            
            % 2. calc expected tuning of voxels with ~N(0,G)
            [sftGaussEV,evarEst_gauss,svarEst_gauss] = pp1_imana('SFT:expectedValue_G',evar,svar,Gtt,numVox,numRun,numSim); % expected value of the null
            
            % 3. calc expected tuning of voxels with ~sparse tuning, but where signal-to-noise equals that voxel data
            [sftSparseEV,evarEst_sp,svarEst_sp] = pp1_imana('SFT:expectedValue_Sparse',evar,svar,1,numVox,numRun,numSim); % 1= tuned to one condition (perfectly sparse)
            %[sftSp2,sftDistSp2] = pp1_imana('SFT:calcExpectedValueSparse',evar,svar,2,numVox,numRun,numSim); % 2= tuned to two conditions
            
            % 4. do prob test for each subject on their gauss sft distribution
            pEV = sum(sftBeta<=sftGaussEV)/length(sftGaussEV);
            if isempty(pEV); pEV = realmin; end
            
            % 5. do prob test for each subject on their sparse sft
            % distributions
            pSp = sum(sftBeta<=sftSparseEV)/length(sftSparseEV);
            if isempty(pSp); pSp = realmin; end
            
            % add to output structure
            d.numVoxF   = v.*numVoxSig;
            d.numVoxTot = v.*size(T.raw_beta{ii},2);
            d.avgF      = v.*mean(F(F>=Fcrit));
            d.sft       = [sftBeta; mean(sftGaussEV); mean(sftSparseEV)];
            d.sftProb   = [0; pEV; pSp];
            d.svar_est  = [svar; mean(svarEst_gauss); mean(svarEst_sp)];
            d.evar_est  = [evar; mean(evarEst_gauss); mean(evarEst_sp)];
            d.snr       = d.svar_est ./ d.evar_est;
            d.isEV      = [0; 1; 2]; 
            D = addstruct(D,d);
            % display to user:
            fprintf('%s\t%02d\t%2.3f\t\t\t%2.4f\t\t%2.4f\t\t%1.5f\n',I.origSN{T.sn(ii)},T.tesselNum(ii),(numVoxSig/numVoxOrig)*100,svar,evar,svar/evar);
        end
        save(fullfile(regDir,sprintf('sft_glm%d_fthres_tessel%d.mat',glm,nTessel)),'-struct','D');
        varargout = {D};   
    case 'TSL:pcmFit'
        % case to fit models to participant data
        sn  = [2:11]; % subj 1 is excluded from analyses (no fieldmaps)
        glm = 4;
        nTessel = pp1_imana('TSL:setThres');
        vararginoptions(varargin,{'roi'});
        
        % get data (load in once for all subjs to save time):
        B = load(fullfile(regDir,sprintf('glm%d_tessels%d_betas.mat',glm,nTessel)));
        D = [];
        for ii=unique(B.roi)'
            fprintf('\nroi %d.',ii);
            tesselNum = unique(B.tesselNum(B.roi==ii));
            % package the data from this roi across subjects
            Y = {};
            partVec = {};
            condVec = {};
            for jj = 1:length(sn)
                % get subject data
                s = sn(jj);
                b = getrow(B,B.sn==s & B.roi==ii);
                bb = [];
                bb.run   = cell2mat(b.run);
                bb.chord = cell2mat(b.tt);
                bb.betas = cell2mat(b.betaW);
                bb = getrow(bb,ismember(bb.chord,1:31)); % restrict to passive stimulation conditions (tt 32 is thumb press)
                % put subj data into pcm variables
                Y{jj}         = bb.betas;
                partVec{jj}   = bb.run;
                condVec{jj}   = bb.chord;
            end
            
            % do fitting:
            d = pp1_encoding('modelFittingWrapper',Y,partVec,condVec);
            
            % add appropriate indexing labels to each row of output data
            % structure
            v = ones(size(d.sn));
            d.sn  = pcm_indicatorMatrix('identity',d.sn)*sn';
            d.roi = v.*ii;
            d.glm = v.*glm;
            d.tesselNum = v.*tesselNum;
            d = rmfield(d,{'g_pred'});
            D = addstruct(D,d);
        end
        save(fullfile(regDir,sprintf('sft_glm%d_PCM_tessel%d.mat',glm,nTessel)),'-struct','D');
        varargout = {D};
    case 'TSL:mapTessels_selectivity'
        hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'
        glm = 4;
        nTessel = pp1_imana('TSL:setThres');
        % calculate normalized selectivity indexes:
        To = load(fullfile(regDir,sprintf('sft_glm%d_fthres_tessel%d.mat',glm,nTessel))); % selectivity values
        T = getrow(To,To.isEV==0);
        T1 = getrow(To,To.isEV==1);
        T2 = getrow(To,To.isEV==2);
        T.sftNorm = (T.sft - T1.sft)./T2.sft;
        clear To T1 T2
        for h=1%:2 % hemisphere
            Th = getrow(T,T.hemi==h);
            G = gifti(fullfile(atlasDir,sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemLetter{h})));
            data = nan(size(G.cdata));
            for r=unique(Th.tesselNum)' % loop through tessels we have analyzed data from
              t = getrow(Th,Th.tesselNum==r);
              idx = find(G.cdata(:,1)==r);
              if (sum(isnan(t.sftNorm))/numel(t.sftNorm))<0.3% ensure more than 70% have data
                data(idx,1) = nanmean(t.sftNorm);
              else
                data(idx,1) = nan;
              end
            end
            outfile = fullfile(wbDir, 'group_164k', sprintf('sft_fthres.tessels%d.%s.glm%d.164k.func.gii', nTessel, hemLetter{h}, glm));
            G       = surf_makeFuncGifti(data,'anatomicalStruct',hname{h},'columnNames',{sprintf('sft_glm%d_fthres',glm)});
            save(G, outfile);
            fprintf('Hemisphere: %d\n',h);
        end
        
    case 'TSL:makeTesselG'
        % wrapper to perform SFT analysis with data from tessels
        % restricts analysis to voxels that show significant finger resposnes (voxel F-stat > F-crit)
        I   = pp1_imana('LIST_subjs');
        sn  = 2:11;
        glm = 4;
        vararginoptions(varargin,{'glm','sn'});
        
        % define which tessels to include in left-hemi:
        hemi = 1;
        [nTessel,prop,distThres] = pp1_imana('TSL:setThres');
        takeTessels = pp1_imana('TSL:select','nTessel',nTessel,'hemi',hemi,'prop',prop,'distThres',distThres);
        
        D = []; % participant-average G structure
        for ii=sn
            % load subject's tessel betas
            T = load(fullfile(regDir,sprintf('glm%d_%s_tessels%d_betas.mat',glm,subj_name{ii},nTessel)));
            T = getrow(T,T.hemi==hemi & ismember(T.tesselNum,takeTessels));
            for jj=takeTessels
                t=getrow(T,T.tesselNum==jj); % data from this tessel region
                b.tt = cell2mat(t.tt);
                b.run = cell2mat(t.run);
                b.betaW = cell2mat(t.betaW);
                try
                    b = getrow(b,b.tt<32); % drop all non-chord regressors
                    Gsubj = pcm_estGCrossval(b.betaW,b.run,b.tt);
                    Gsubj = pcm_makePD(Gsubj);
                    d.g = rsa_vectorizeIPM(Gsubj);
                    d.hasData = true;
                catch % no data from this tessel for this participant
                    d.g = nan(1,496);
                    d.hasData = false;
                end
                d.sn = t.sn;
                d.roi = t.roi;
                d.hemi = t.hemi;
                D=addstruct(D,d);
            end
        end
        D=tapply(D,{'roi'},{'g','nanmean'});
        save(fullfile(regDir,sprintf('glm%d_tessels%d_G.mat',glm,nTessel)),'-struct','D');
        varargout = {D};
    case 'TSL:encodingModels'
        sn  = 2:11;
        glm = 4;
        
        % define which tessels to include in left-hemi:
        hemi = 1;
        [nTessel,prop,distThres] = pp1_imana('TSL:setThres');
        takeTessels = pp1_imana('TSL:select','nTessel',nTessel,'hemi',hemi,'prop',prop,'distThres',distThres);
        
        G = load(fullfile(regDir,sprintf('glm%d_tessels%d_G.mat',glm,nTessel)));
        D=[];
        for ii=sn
            Y = {}; partVec = {}; condVec = {}; d=[];
            % load subject's tessel betas
            T = load(fullfile(regDir,sprintf('glm%d_%s_tessels%d_betas.mat',glm,subj_name{ii},nTessel)));
            T = getrow(T,T.hemi==hemi & ismember(T.tesselNum,takeTessels));
            % arrange into data structure for encoding model fitting:
            for jj=takeTessels
                t=getrow(T,T.tesselNum==jj); % data from this tessel region
                b.tt = cell2mat(t.tt);
                b.run = cell2mat(t.run);
                b.betaW = cell2mat(t.betaW);
                try
                    b = getrow(b,b.tt<32); % drop all non-chord regressors
                    Y{end+1}       = b.betaW;
                    partVec{end+1} = b.run;
                    condVec{end+1} = b.tt;
                catch % no data from this tessel for this participant
                    Y{end+1}       = [];
                    partVec{end+1} = [];
                    condVec{end+1} = [];
                    continue
                end
            end
            % fit encoding models to data from these tessels in this
            % subject:
            d = pp1_encoding('encoding_Passive_tessels',Y,partVec,condVec,G,takeTessels);
            v = ones(size(d.tessel));
            d.sn = v.*ii;
            save(fullfile(regDir,sprintf('glm%d_%s_tessels%d_encodingFits.mat',glm,subj_name{ii},nTessel)),'-struct','d');
            D = addstruct(D,d);
        end
        save(fullfile(regDir,sprintf('glm%d_tessels%d_repModelFits.mat',glm,nTessel)),'-struct','D');
        varargout = {D};
    case 'TSL:mapTessels_encoding_medianFits'
        hname = {'CortexLeft','CortexRight'}; 
        glm = 4;
        nTessel = pp1_imana('TSL:setThres');
        % load model fits for thresholded tessels:
        T = load(fullfile(regDir,sprintf('glm%d_tessels%d_repModelFits.mat',glm,nTessel)));
        % calculate normalized model fits:
        nCeil = 9;
        nNull = 1;
        numModels=9;
        T.r_norm = T.r_test - kron(T.r_test(T.model==nNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1));
        
        % map left hemi
        h = 1;
        G = gifti(fullfile(atlasDir,sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemLetter{h})));
        data = nan(size(G.cdata,1),10);
        for r=unique(T.tessel)' % loop through tessels we have analyzed data from
          t = getrow(T,T.tessel==r);
          idx = G.cdata(:,1)==r;
          for ii=2:8
             data(idx,ii-1) = median(t.r_norm(t.model==ii));
          end
          data(idx,ii) = median(t.r_norm(t.model==3)-t.r_norm(t.model==2)); % 2f int - linear model fits
          data(idx,ii+1) = median(t.r_norm(t.model==3)-t.r_norm(t.model==6)); % 2f int - flexible model fits
          data(idx,ii+2) = median(1-t.r_norm(t.model==3)); % 1 - 2f int fits
        end
        colnames = {'linear_fit','2f_fit','3f_fit','4f_fit','flex_fit','2f_noAdj','2f_noNonAdj','2f-lin_fits','2f-flex_fits','nc-2f fits'};
        outfile = fullfile(wbDir, 'group_164k', sprintf('repModelFits_median.tessels%d.%s.glm%d.164k.func.gii', nTessel, hemLetter{h}, glm));
        G       = surf_makeFuncGifti(data,'anatomicalStruct',hname{h},'columnNames',colnames);
        save(G, outfile);
        fprintf('Hemisphere: %d\n',h);
    case 'TSL:mapTessels_encoding_signtest'
        hname = {'CortexLeft','CortexRight'}; 
        glm = 4;
        nTessel = pp1_imana('TSL:setThres');
        % load model fits for thresholded tessels:
        T = load(fullfile(regDir,sprintf('glm%d_tessels%d_repModelFits.mat',glm,nTessel)));
        % calculate normalized model fits:
        nCeil = 9;
        nNull = 1;
        numModels=9;
        T.r_norm = T.r_test - kron(T.r_test(T.model==nNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1));
        
        % map left hemi
        h = 1;
        G = gifti(fullfile(atlasDir,sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemLetter{h})));
        data = nan(size(G.cdata,1),10);
        for r=unique(T.tessel)' % loop through tessels we have analyzed data from
          t = getrow(T,T.tessel==r);
          idx = G.cdata(:,1)==r;
          % model fits:
          for ii=2:8
             p=signtest(t.r_norm(t.model==ii)); %2-sided signtest against 0 (null fit)
             if p<0.05
                 data(idx,ii-1) = median(t.r_norm(t.model==ii));
             end
          end
          % difference in model fits:
          x = t.r_norm(t.model==3) - t.r_norm(t.model==2); % 2f int - linear model fits
          p=signtest(x);
          if p<0.05
              data(idx,ii) = median(x);
          end
          x = t.r_norm(t.model==3) - t.r_norm(t.model==6); % 2f int - flexible model fits
          p=signtest(x);
          if p<0.05
              data(idx,ii+1) = median(x);
          end
          x = 1 - t.r_norm(t.model==3); % noise ceiling - 2f int fits
          p=signtest(x);
          if p<0.05
              data(idx,ii+2) = median(x);
          end
        end
        colnames = {'linear_fit','2f_fit','3f_fit','4f_fit','flex_fit','2f_noAdj','2f_noNonAdj','2f-lin_fits','2f-flex_fits','nc-2f fits'};
        outfile = fullfile(wbDir, 'group_164k', sprintf('repModelFits_signtest.tessels%d.%s.glm%d.164k.func.gii', nTessel, hemLetter{h}, glm));
        G       = surf_makeFuncGifti(data,'anatomicalStruct',hname{h},'columnNames',colnames);
        save(G, outfile);
        fprintf('Hemisphere: %d\n',h);
         

    case 'TSL:varExpMap'
        % map best fitting models using % variance explained criterion
        % Map is purely for visualization purposes!
        hname = {'CortexLeft','CortexRight'}; 
        glm = 4;
        nTessel = pp1_imana('TSL:setThres');
        varExpThres = 0.95;
        % load model fits for thresholded tessels:
        T = load(fullfile(regDir,sprintf('glm%d_tessels%d_repModelFits.mat',glm,nTessel)));
        % calculate normalized model fits:
        nCeil = 9;
        nNull = 1;
        numModels=9;
        T.r_norm = T.r_test - kron(T.r_test(T.model==nNull),ones(numModels,1));
        T.r_norm = T.r_norm./kron(T.r_norm(T.model==nCeil),ones(numModels,1));
        T = getrow(T,ismember(T.model,[3,4,5,6])); % linear, 2f, 3f, 4f, flexible
        T.model(T.model==2)=1; 
        T.model(T.model==6)=2; % linear, flexible, 2f, 3f, 4f
        % map left hemi
        h = 1;
        G = gifti(fullfile(atlasDir,sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemLetter{h})));
        data = zeros(size(G.cdata,1),1);
        for r=unique(T.tessel)' % loop through tessels we have analyzed data from
          t = getrow(T,T.tessel==r);
          idx = G.cdata(:,1)==r;
          r_norm = pivottable([],t.model,t.r_norm,'median');
          best_model = find(r_norm >= varExpThres,1);
          if isempty(best_model)
              best_model=0;
          end
          data(idx,1) = int32(best_model); % best simplest model
        end
        modelNames = {'linear','flexible','2finger','3finger','4finger'};
        outfile = fullfile(wbDir, 'group_164k', sprintf('repModelFits_WTA_median.tessels%d.%s.glm%d.164k.label.gii', nTessel, hemLetter{h}, glm));
        
        rgba = [0 0 0 0;...
                69 114 180 255;...
                116 173 209 255;...
                254 224 144 255;...
                253 174 97 255;...
                244 109 67 255;...
                215 48 39 255;...
                255 95 41 255]; 
        rgba=rgba./255;
        
        G = surf_makeLabelGifti(data,...
                    'columnNames',{'WTA_model_median'},...
                    'labelNames',modelNames,...
                    'labelRGBA',rgba,...
                    'anatomicalStruct',hname{h});
        
        
        save(G, outfile);
        fprintf('Hemisphere: %d\n',h);
        
    case '0' % ------------ PLOT: figures. --------------------------------
        % The FIG cases usually harvest some data from ROI stats structures.
        % Sometimes they also do stats (since the data is harvested
        % accordingly).
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'PLOT_avgBetas'
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

    case 'PLOT_RDMline_LDC'
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
        sn    = pp1_imana('getSubjs');
        glm   = 4;
        roi   = [1:7];
        vararginoptions(varargin,{'sn','roi','glm','conds'});
        D=[];
        xtick = {};
        for rr=roi
            d = load(fullfile(regDir,sprintf('patternReliability_roi%d_glm%d.mat',rr,glm)));
            D=addstruct(D,d);
            roiName = [regname{rr-(regSide(rr)-1)*length(regname)} ' (' hem{regSide(rr)} ')'];
            xtick = {xtick{:},roiName,''};
        end
        D = getrow(D,ismember(D.sn,sn));
        plt.bar(D.roi,D.corr,'split',D.type);
        plt.labels('','Pearson''s r','raw pattern reliability');
        plt.set('xticklabel',xtick,'xticklabelrotation',45);
        legend off
        drawline(0,'dir','horz');
%         ylims = ylim;
%         ylim([0 ylims(2)]);
        
        varargout = {D};
    case 'PLOT_PSCperNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn  = 2:11;%pp1_imana('getSubjs');
        glm = 4;
        roi = [2:4]; %[1:7,12:14,15:21,26:28]; % 
        vararginoptions(varargin,{'sn','roi','glm'});
        
        D  = pp1_imana('ROI_getChordPSC','glm',glm);
        D  = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        D  = getrow(D,D.numDigits<6);
        Dr = tapply(D,{'sn','roi','numDigits'},{'psc','mean'});
        % plot
        style.use('numDigitsPurple');
        %sty = style.custom({[0.75 0.75 0.75],[0.5 0.5 0.5],[0.25 0.25 0.25],[0 0 0]});
        %sty = style.custom(plt.helper.get_shades(numel(roi),'jet'));
        %sty =style.custom({[0.2549 0.25882 0.5451],[0 0.47451 0.56471],[0 0.67059 0.50588],[0.33 0.83529 0.20784]});
        plt.box([Dr.roi Dr.numDigits],Dr.psc,'split',Dr.numDigits);
        %plt.line([Dr.numDigits],Dr.psc,'split',Dr.roi,'style',sty);
%         xtick = {};
%         for r = roi
%             roiName = [regname{r-(regSide(r)-1)*length(regname)} ' (' hem{regSide(r)} ')'];
%             xtick = {xtick{:} '','',roiName,'',''};
%         end
%         plt.set('xticklabel',xtick,'xticklabelrotation',45);
        plt.legend('east',{'1 digit','2 digits','3 digits','4 digits','5 digits'});
        plt.labels('# digits stimulated','% signal change','task-modulated scaling of activity');
        drawline(0,'dir','horz');
        
        varargout = {Dr,D};
    case 'PLOT_PSCperNumDigits_lines'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn  = 2:11;%pp1_imana('getSubjs');
        glm = 4;
        roi = [1:7];%,15:19]; %[1:7,12:14,15:21,26:28]; % 
        vararginoptions(varargin,{'sn','roi','glm'});
        
        D  = pp1_imana('ROI_getChordPSC','glm',glm);
        D  = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi));
        D  = getrow(D,D.numDigits<6); % drop thumb press responses
        Dr = tapply(D,{'sn','roi','numDigits'},{'psc','mean'}); % avg. across chords per # digit in each subject
        % plot
%        CAT_s1.linecolor = {[0.39216 0.56078 1],[ 0.47059 0.36863 0.94118],[0.86275 0.14902 0.49804],[0.99608 0.38039 0],[1 0.6902 0]};
      %  CAT_s1.linecolor = {[0.43922 0 0.33333],[0.4549 0.15686 0.6549],[0.81569 0 0.96078],[0.89412 0.56471 1]};
        CAT_s1.linecolor = {[0.6 0.6 0.6] [0.5 0 0] [0.9 0 0] [1 0.6 0] [0 0 0.8] [0.3 0.6 0.9],[0.6 0 0.9]};  
        CAT_s1.linestyle = {'-'};%{'-','-','-','-','-',':',':',':',':',':'};
        CAT_s1.markersize = 8;
        CAT_s1.errorcolor = CAT_s1.linecolor;
        CAT_s1.patchcolor = CAT_s1.linecolor;
        CAT_s1.transp = 0.15;
        CAT_s1.markertype = 'none';
        CAT_s1.markerfill = CAT_s1.linecolor;
        CAT_s1.linewidth = 2.5;
        CAT_s1.errorwidth = 1.5;
        
        
        % [0 0 0],[0.9 0 0],[1 0.6 0]
        
%         dat = pivottable([Dr.roi Dr.sn],Dr.numDigits,Dr.psc,'mean');
%         roiVec = kron(roi',ones(numel(sn),1));
%         traceplot(1:5,dat,'split',roiVec,'CAT',CAT_s1,'errorfcn','stderr');
        lineplot(Dr.numDigits,Dr.psc,'split',Dr.roi,'CAT',CAT_s1,'errorfcn','stderr');
        drawline(0,'dir','horz','linestyle','-');
        title('task-modulated BOLD activity');
        ylabel('% signal change');
        xlabel('# digit stimulated');

        varargout = {Dr,D};
        
    case 'PLOT_RDMsquare_avg'
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = 1;
        chords = 1:31;
        vararginoptions(varargin,{'sn','roi','glm','chords'});
        % load data
        T = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat']));
        T = getrow(T,ismember(T.sn,sn) & T.roi==roi);
        % get rdm
        rdmVec = T.corr_dist;
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
        %patchimg(rdm);
        imagesc(rdm);
        title(regname{roi});
        set(gca,'xtick',[0.5:length(chords)],'xticklabel',chords);
        set(gca,'ytick',[0.5:length(chords)],'yticklabel',chords);
%         ylim([0 length(chords)]);
%         xlim([0 length(chords)]);
        
        varargout = {rdm,interSubjCorr};
    case 'PLOT_G_avg'
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = 1;
        chords = 1:31;
        vararginoptions(varargin,{'sn','roi','glm','chords'});
        % load data
        T = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat']));
        T = getrow(T,ismember(T.sn,sn) & T.roi==roi);
        % get second moment
        ipmVec = T.G_wmean; %ssqrt(T.ldc);
        ipm = [];
        for i = 1:size(ipmVec,1)
            tmp = ipmVec(i,:);
            tmp = rsa_squareIPM(tmp);
            tmp = tmp(chords,chords);
            ipm(i,:) = rsa_vectorizeIPM(tmp);
        end
        % quickly calculate inter-subject rdm correlations
        r = corr(ipm');
        interSubjCorr = 1-squareform(1-r)';
        % plot
        ipm = rsa_squareIPM(mean(ipm));
        %patchimg(ipm);
        imagesc(ipm);
        title(regname{roi});
        set(gca,'xtick',[0.5:length(chords)],'xticklabel',chords);
        set(gca,'ytick',[0.5:length(chords)],'yticklabel',chords);
%         ylim([0 length(chords)]);
%         xlim([0 length(chords)]);
        
        varargout = {ipm,interSubjCorr};
        
    case 'PLOT_avgDist'
        % plots avg. distances b/t all chord conditions per roi
        dist = 'ldc'; % or corr
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = [1:7,12]; % [3a, 3b, 1, 2, M1]
        vararginoptions(varargin,{'sn','roi','glm','dist'});
        % load data
        T = load(fullfile(regDir,['glm' num2str(glm) '_reg_Toverall.mat']));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        % drop the distance b/t thumb press and chord stimulations (glm==3, 4)
        if glm>2
            for ii=1:size(T.sn,1)
                rdm = rsa_squareRDM(T.ldc_wmean(ii,:));
                rdm = rdm(1:31,1:31);
                T.LDC(ii,:) = rsa_vectorizeRDM(rdm);
                rdm = rsa_squareRDM(T.corr_dist(ii,:));
                rdm = rdm(1:31,1:31);
                T.COS(ii,:) = rsa_vectorizeRDM(rdm);
            end
        end
        switch dist
            case 'ldc'
                T.dist = mean(T.LDC,2);
            case 'cos'
                T.dist = mean(T.COS,2);
        end
        
        % plot
        xtick = {};
        for r = roi
            roiName = [regname{r-(regSide(r)-1)*length(regname)} ' (' hem{regSide(r)} ')'];
            xtick = {xtick{:},roiName};
        end
        plt.line(T.roi,T.dist);
        plt.set('xticklabel',xtick);
        plt.labels('',sprintf('avg. %s b/t chords',dist),'avg. distance/region');
        drawline(0,'dir','horz');
        
        varargout = {T};
        
        
    case 'PLOT_LDCperNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn  = pp1_imana('getSubjs');
        glm = 4;
        roi = 1:4; % [3a, 3b, 1, 2, M1]
        vararginoptions(varargin,{'sn','roi','glm'});
        
        D = pp1_imana('GET_LDCperNumDigits','sn',sn,'glm',glm,'roi',roi);
        %D.ldc = ssqrt(D.ldc);
        D = tapply(D,{'sn','roi','numDigits'},{'ldc_wmean','mean'},{'ldc','mean'},{'corr_dist','mean'});
        % plot
        style.use('numDigits');
        plt.box(D.roi,D.ldc,'split',D.numDigits);
        plt.labels('roi','ldc^2 between chords with same # digits')
        drawline(0,'dir','horz');
        varargout = {D};
    case 'PLOT_LDCacrossNumDigits'
        % plot distance between chords with differing number of digits.
        sn    = pp1_imana('getSubjs');
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
    case 'PLOT_LDCacrossAvgNumDigits'
        % plot distance between chords with differing number of digits.
        % distances are calculated between avg. chord patterns (done in
        % ROI_stats);
        sn    = pp1_imana('getSubjs');
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
    case 'PLOT_plotLDCS1'
        sn = pp1_imana('getSubjs');
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

    case 'PLOT_sft_MF'
        % plot tuning across chords analysis results
        roi  = 2:4;
        glm = 4;
        vararginoptions(varargin,{'roi','glm'});
        % load data
        T = load(fullfile(regDir,sprintf('sft_mf_glm%d.mat',glm)));
        T = getrow(T,ismember(T.roi,roi) & ismember(T.sn,[1:8]));
        % normalize fits by expected[mvnrnd] and expected[sparse]:
        numD = 1:4;
        D = [];
        for r = roi
            for f=numD
               t = getrow(T,T.roi==r & T.numDigits==f); 
               d.sftNorm = (t.sft(t.isEV==0) - t.sft(t.isEV==1)) ./ (t.sft(t.isEV==2) - t.sft(t.isEV==1));
               t = getrow(t,t.isEV==0);
               d.sn = t.sn;
               d.roi = t.roi;
               d.glm = t.glm;
               d.r2  = t.r2;
               d.numDigits = t.numDigits;
               d.passive = t.passive;
               D = addstruct(D,d);
            end
        end
        
        % plot normalized sft
        sty = style.custom(plt.helper.get_shades(numel(roi),'jet','descend'));
        subplot(1,2,1);
        plt.line(D.numDigits,D.sftNorm,'split',D.roi,'style',sty);
        ylabel('single-finger tuning index');
        xlabel('# fingers stimulated');
        title('PASSIVE');
        ylim([-0.2 1]);
        % plot r2 of estimated single-finger patterns
        subplot(1,2,2);
        plt.line(D.numDigits,D.r2,'split',D.roi,'style',sty);
        title('single-finger pattern estimates')
        ylabel('r2');
        xlabel('# fingers stimulated');
        ylim([0 1]);

        varargout = {D,T};
    
    case 'PLOT_nnClassification'
        % plots nearest neighbour classification accuracies for regions
        % (accuracies for all chords and single fingers in separate
        % subplots)
        glm = 4;
        roi = [1:12];
        sn  = pp1_imana('getSubjs');
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
        plt.bar(D.roi,D.acc.*100,'split',D.isEV,'style',sty,'subset',D.sf==0);
        title('all conditions');
        ylabel('classification accuracy (%)');
        xlabel('region');
        drawline((1/31)*100,'dir','horz','linestyle',':');
        % plot single finger accuracy
        subplot(1,2,2);
        plt.bar(D.roi,D.acc.*100,'split',D.isEV,'style',sty,'subset',D.sf==1);
        title('single fingers');
        ylabel('classification accuracy (%)');
        xlabel('region');
        drawline(20,'dir','horz','linestyle',':');

        varargout = {D};

    case '0' % ------------ encoding models -------------------------------
    case 'doEncoding'
        % wrapper to do encoding models for different ROIs:
        roi = [1:4,7,9];
        D=[];
        for rr=roi
           d = pp1_imana('encoding_noCV_oneROI','roi',rr);
           D=addstruct(D,d);
        end
        % plot
        labels = {'','','BA 3a','','','','','BA 3b','','','','','BA 1','','','','','BA 2','','','','','SII','','','','','M1','',''};
        sty = style.custom({[0 1 0],[1 0 0],[0 0 1],[0 0 0]});
        subplot(2,1,1);
        plt.bar([D.roi D.model],D.r2,'style',sty,'split',D.model);
        set(gca,'xticklabel',labels);
        ylabel('R2');
        drawline(0,'dir','horz');
        title('subject-level encoding models (not cv)')
        legend off
        
        subplot(2,1,2);
        plt.bar([D.roi D.model],D.r,'style',sty,'split',D.model);
        set(gca,'xticklabel',labels);
        ylabel('Pearson''s R');
        %drawline(0,'dir','horz');
        
        varargout = {D};
    case 'doEncodingCV'
        % wrapper to do encoding models for different ROIs:
        roi = [1:4,7,9];
        D=[];
        for rr=roi
           d = pp1_imana('encoding_CV_oneROI','roi',rr);
           % avg. across CV folds:
           d = tapply(d,{'roi','sn','model'},{'act_pred','mean'},{'act_true','mean'},...
               {'r','mean'},{'r2','mean'},{'rpool','mean'},{'r2pool','mean'});
           D=addstruct(D,d);
        end
        % plot
        labels = {'','','BA 3a','','','','','BA 3b','','','','','BA 1','','','','','BA 2','','','','','SII','','','','','M1','',''};
        sty = style.custom({[0 1 0],[1 0 0],[0 0 1],[0.5 0.5 0.5],[0 0 0]});
        subset = D.model~=4;
        subplot(2,1,1);
        plt.bar([D.roi D.model],D.r2,'style',sty,'split',D.model,'subset',subset);
        set(gca,'xticklabel',labels);
        ylabel(sprintf('R2 \n(avg. over CV folds)'));
        drawline(0,'dir','horz');
        title('subject-level CV encoding models');
        legend off
        
%         subplot(2,2,2);
%         plt.bar([D.roi D.model],D.r2pool,'style',sty,'split',D.model,'subset',subset);
%         set(gca,'xticklabel',labels);
%         ylabel('R2 (pooled)');
%         drawline(0,'dir','horz');
        
        subplot(2,1,2);
        plt.bar([D.roi D.model],D.r,'style',sty,'split',D.model,'subset',subset);
        set(gca,'xticklabel',labels);
        ylabel(sprintf('Pearson''s R \n(avg. over CV folds)'));
        %drawline(0,'dir','horz');
        
%         subplot(2,2,4);
%         plt.bar([D.roi D.model],D.rpool,'style',sty,'split',D.model,'subset',subset);
%         set(gca,'xticklabel',labels);
%         ylabel('R (pooled)');
%         %drawline(0,'dir','horz');
        
        varargout = {D};
    case 'encoding_noCV_oneROI'
        % encoding model- no crossvalidation
        % Here, we predict multifinger patterns (U_mf_hat) as the linear
        % summation of constituent single finger patterns (U_sf).
        % This is done using patterns averaged across runs, no run mean
        % removal.
        % We evalute the fit of U_mf_hat to the true multifinger patterns
        % (U_mf) using R (pearson's corr) and R2. 
        % This is done for each subject.
        sn  = pp1_imana('getSubjs');
        roi = [];
        glm = 4;
        vararginoptions(varargin,{'roi'});
        model = {'null','summation','summationNL','noiseceiling'};
        chords = pp1_imana('chords');
        [Y,partVec,condVec] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm,'rmvMean',0,'posJustify',0);
        numDigits = sum(chords(6:31,:),2);
        CD = pcm_indicatorMatrix('identity',numDigits);
        D=[];
        for ii=1:numel(Y) % for each subject
            d.sn=sn(ii);
            d.roi=roi;
            d.glm=glm;
            % avg. betas across runs:
            C0 = pcm_indicatorMatrix('identity',condVec{ii});
            U  = pinv(C0)*Y{ii};
            % split into single and multi-finger chords:
            U_sf = U(1:5,:); % single finger
            U_mf = U(6:31,:); % multi finger
            for mm=1:numel(model) % for each model
                % predict multi finger chords according to model type
                switch model{mm}
                    case 'null' % fit an intercept (the multi-finger mean)
                        theta = mean(U_mf(:));
                        U_mf_hat = ones(size(U_mf)).*theta;
                    case 'summation'
                        X = chords(6:31,:);
                        U_mf_hat = X*U_sf;
                        theta=nan;
                    case 'summationNL'
                        X = chords(6:31,:);
                        U_mf_hat = X*U_sf; % do summation
                        % estimate theta
                        [theta,err,gEst] = fminsearch(@(x) fitPower(x,U_mf_hat(:),U_mf(:)),0.5);
                        % squish patterns with theta param:
                        sign_Uhat= sign(U_mf_hat);
                        abs_Uhat = abs(U_mf_hat);
                        U_mf_hat = sign_Uhat.*(abs_Uhat.^theta);
                    case 'noiseceiling' % as a test- evaluation metrics should be perfect
                        U_mf_hat = U_mf;
                        theta=nan;
                end
                % evaluate model predictions
                % R2:
                rss = sum((U_mf(:)-U_mf_hat(:)).^2);
                tss = sum(U_mf(:).*U_mf(:));
                d.r2 = 1-(rss/tss);
                % R:
                ssC = sum(U_mf(:).*U_mf_hat(:));
                ss1 = sum(U_mf(:).*U_mf(:));
                ss2 = sum(U_mf_hat(:).*U_mf_hat(:));
                d.r = ssC/sqrt(ss1.*ss2);
                %d.rtest = corr(U_mf(:),U_mf_hat(:)); 
                d.model = mm;
                d.theta = theta;
                % calculate overall activity estimates per numdigits:
                d.act_pred = mean(pinv(CD)*U_mf_hat,2)';
                d.act_true = mean(pinv(CD)*U_mf,2)';
                D=addstruct(D,d);
            end
        end
        varargout = {D};
    case 'encoding_CV_oneROI'
        % simplest encoding model- WITH crossvalidation
        % Here, we predict multifinger patterns (U_mf_hat) as the linear
        % summation of constituent single finger patterns (U_sf).
        % This is done in a leave-one-out crossvalidated fashion, where we:
        % - avg. patterns across runs in training set
        % - use avg. single finger patterns to predict multi-finger patterns in test set
        % We evalute the fits on each fold using R (pearson's corr) and R2. 
        % This is done for each subject.
        sn  = pp1_imana('getSubjs');
        roi = [];
        glm = 4;
        vararginoptions(varargin,{'roi'});
        model = {'null','summation','summationNL','summationNLScaled','noiseCeilingCV','noiseCeiling'};
        crossvalScheme='leaveOneOut';
        chords = pp1_imana('chords');
        [Y,partVec,condVec] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm,'rmvMean',0,'posJustify',0);
        numDigits = sum(chords(6:31,:),2);
        CD = pcm_indicatorMatrix('identity',numDigits);
        D=[];
        for ii=1:numel(Y) % for each subject
            d.sn=sn(ii);
            d.roi=roi;
            d.glm=glm;
            % define crossvalidation folds:
            cV      = condVec{ii};
            pV      = partVec{ii};
            parts   = unique(pV);
            numPart = numel(parts);
            partI   = {};
            switch crossvalScheme
                case 'leaveOneOut'
                    for p=1:numPart
                        partI{p}=parts(p);
                    end
                case 'evenOdd'
                    for p=1:2
                        partI{p}=parts(mod(parts+p,2)==0);
                    end
            end
            % avg. betas across runs (for overall noise ceiling):
            C0 = pcm_indicatorMatrix('identity',condVec{ii});
            U  = pinv(C0)*Y{ii};
            % loop through crossval folds and test models:
            numFolds = numel(partI);
            Dp=[];
            for p=1:numFolds
                trainIdx = ~ismember(pV,partI{p}); 
                testIdx  = ismember(pV,partI{p}); 
                % make design matrices for training and test splits
                Xtrain = pcm_indicatorMatrix('identity',cV(trainIdx));
                Xtest  = pcm_indicatorMatrix('identity',cV(testIdx));
                % get data for training and test, avg. across conditions
                % within each split:
                Ytrain = pinv(Xtrain) * Y{ii}(trainIdx,:);
                Ytest  = pinv(Xtest) * Y{ii}(testIdx,:);
                % split into single finger and multi finger data
                Y_sf_train = Ytrain(1:5,:);
                Y_mf_train = Ytrain(6:31,:);
                Y_sf_test  = Ytest(1:5,:);
                Y_mf_test  = Ytest(6:31,:); % multi finger chords test set
                % loop through models and predict mf chords accordingly:
                for mm=1:numel(model) % for each model
                    % predict multi finger chords according to model type
                    switch model{mm}
                        case 'null' % fit an intercept (the multi-finger mean)
                            theta = mean(Y_mf_train(:));
                            Y_mf_hat = ones(size(Y_mf_train)).*theta;
                        case 'summation'
                            X        = chords(6:31,:);
                            Y_mf_hat = X * Y_sf_train;
                            theta=nan;
                        case 'summationNL'
                            X        = chords(6:31,:);
                            Y_mf_hat = X * Y_sf_train;
                            % estimate theta using training data
                            [theta,err,gEst] = fminsearch(@(x) fitPower(x,Y_mf_hat(:),Y_mf_train(:)),0.5);
                            % squish patterns with theta param:
                            sign_Yhat= sign(Y_mf_hat);
                            abs_Yhat = abs(Y_mf_hat);
                            Y_mf_hat = sign_Yhat.*(abs_Yhat.^theta);
                        case 'summationNLScaled'
                            X        = chords(6:31,:);
                            Y_mf_hat = X * Y_sf_train;
                            % estimate theta using training data
                            [theta,err,gEst] = fminsearch(@(x) fitPowerScaling(x,Y_mf_hat(:),Y_mf_train(:)),[0.5,1]);
                            % squish patterns with theta param:
                            sign_Yhat= sign(Y_mf_hat);
                            abs_Yhat = abs(Y_mf_hat);
                            Y_mf_hat = theta(2)*sign_Yhat.*(abs_Yhat.^theta(1));
                        case 'noiseCeilingCV' % crossvalidated fits
                            Y_mf_hat = Y_mf_train; % actual mf chords from training data
                            theta=nan;
                        case 'noiseCeiling' % full data fit
                            Y_mf_hat = U(6:31,:); % overall avg. data (includes all runs)
                            theta=nan;
                    end
                    % evaluate model predictions
                    % R2:
                    rss = sum((Y_mf_test(:)-Y_mf_hat(:)).^2);
                    tss = sum(Y_mf_test(:).*Y_mf_test(:));
                    d.r2 = 1-(rss/tss);
                    d.rss=rss;
                    d.tss=tss;
                    % R:
                    ssC = sum(Y_mf_test(:).*Y_mf_hat(:));
                    ss1 = sum(Y_mf_test(:).*Y_mf_test(:));
                    ss2 = sum(Y_mf_hat(:).*Y_mf_hat(:));
                    d.r = ssC/sqrt(ss1.*ss2);
                    d.ssC=ssC;
                    d.ss1=ss1;
                    d.ss2=ss2;
                    % info fields
                    d.model = mm;
                    d.modelName = {model{mm}};
                    d.theta = {theta};
                    d.fold  = p;
                    % calculate overall activity estimates per numdigits:
                    d.act_pred = mean(pinv(CD)*Y_mf_hat,2)';
                    d.act_true = mean(pinv(CD)*Y_mf_test,2)';
                    Dp=addstruct(Dp,d);
                end
            end
            mX = pcm_indicatorMatrix('identity',Dp.model);
            % calculate pooled R2 metric:
            RSS = pivottable([],Dp.model,Dp.rss,'sum');
            TSS = pivottable([],Dp.model,Dp.tss,'sum');
            R2pool = 1-(RSS./TSS);
            Dp.r2pool = mX*R2pool';
            % calculate pooled R metric (sum vars and covars across folds):
            SSC   = pivottable([],Dp.model,Dp.ssC,'sum');
            SS1   = pivottable([],Dp.model,Dp.ss1,'sum');
            SS2   = pivottable([],Dp.model,Dp.ss2,'sum');
            Rpool = SSC./sqrt(SS1.*SS2);
            Dp.rpool = mX*Rpool';
            D=addstruct(D,Dp);
        end
        varargout = {D};
    case 'encoding_CV_oneROI_scaling'
        % simplest encoding model- WITH crossvalidation
        % Here, we predict multifinger patterns (U_mf_hat) as the linear
        % summation of constituent single finger patterns (U_sf).
        % This is done in a leave-one-out crossvalidated fashion, where we:
        % - avg. patterns across runs in training set
        % - use avg. single finger patterns to predict multi-finger patterns in test set
        % We evalute the fits on each fold using R (pearson's corr) and R2. 
        % This is done for each subject.
        sn  = pp1_imana('getSubjs');
        roi = [];
        glm = 4;
        vararginoptions(varargin,{'roi'});
        model = {'null','summation','summationNL','summationNLScaled','noiseCeilingCV','noiseCeiling'};
        fitScale = 1; % include scaling param for each cv fold that accounts for changes in overall signal strength
        crossvalScheme='leaveOneOut';
        % get data
        chords = pp1_imana('chords');
        [Y,partVec,condVec] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm,'rmvMean',0,'posJustify',0);
        D=[];
        for ii=1:numel(Y) % for each subject
            d.sn=sn(ii);
            d.roi=roi;
            d.glm=glm;
            d.scaleParam=0;
            % define crossvalidation folds:
            cV      = condVec{ii};
            pV      = partVec{ii};
            parts   = unique(pV);
            numPart = numel(parts);
            partI   = {};
            switch crossvalScheme
                case 'leaveOneOut'
                    for p=1:numPart
                        partI{p}=parts(p);
                    end
                case 'evenOdd'
                    for p=1:2
                        partI{p}=parts(mod(parts+p,2)==0);
                    end
            end
            % avg. betas across runs (for overall noise ceiling):
            C0 = pcm_indicatorMatrix('identity',condVec{ii});
            U  = pinv(C0)*Y{ii};
            % loop through crossval folds and test models:
            numFolds = numel(partI);
            for p=1:numFolds
                trainIdx = ~ismember(pV,partI{p}); 
                testIdx  = ismember(pV,partI{p}); 
                % make design matrices for training and test splits
                Xtrain = pcm_indicatorMatrix('identity',cV(trainIdx));
                Xtest  = pcm_indicatorMatrix('identity',cV(testIdx));
                % get data for training and test, avg. across conditions
                % within each split:
                Ytrain = pinv(Xtrain) * Y{ii}(trainIdx,:);
                Ytest  = pinv(Xtest) * Y{ii}(testIdx,:);
                % split into single finger and multi finger data
                Y_sf_train = Ytrain(1:5,:);
                Y_mf_train = Ytrain(6:31,:);
                Y_sf_test  = Ytest(1:5,:);
                Y_mf_test  = Ytest(6:31,:); % multi finger chords test set
                % fit scalar to account for different SNRs across CV folds:
                % (use data from all conditions for this estimate)
                if fitScale
                    scaleParam = mean(Ytest(:)) / mean(Ytrain(:));
                end
                % loop through models and predict mf chords accordingly:
                for mm=1:numel(model) % for each model
                    % predict multi finger chords according to model type
                    switch model{mm}
                        case 'null' % fit an intercept (the multi-finger mean)
                            theta = mean(Y_mf_train(:));
                            Y_mf_hat = ones(size(Y_mf_train)).*theta;
                        case 'summation'
                            X        = chords(6:31,:);
                            Y_mf_hat = X * Y_sf_train;
                            theta=nan;
                        case 'summationNL'
                            X        = chords(6:31,:);
                            Y_mf_hat = X * Y_sf_train;
                            % estimate theta using training data
                            [theta,err,gEst] = fminsearch(@(x) fitPower(x,Y_mf_hat(:),Y_mf_train(:)),0.5);
                            % squish patterns with theta param:
                            sign_Yhat= sign(Y_mf_hat);
                            abs_Yhat = abs(Y_mf_hat);
                            Y_mf_hat = sign_Yhat.*(abs_Yhat.^theta);
                        case 'summationNLScaled'
                            X        = chords(6:31,:);
                            Y_mf_hat = X * Y_sf_train;
                            % estimate theta using training data
                            [theta,err,gEst] = fminsearch(@(x) fitPowerScaling(x,Y_mf_hat(:),Y_mf_train(:)),[0.5,1]);
                            % squish patterns with theta param:
                            sign_Yhat= sign(Y_mf_hat);
                            abs_Yhat = abs(Y_mf_hat);
                            Y_mf_hat = theta(2)*sign_Yhat.*(abs_Yhat.^theta(1));
                        case 'noiseCeilingCV' % crossvalidated fits
                            Y_mf_hat = Y_mf_train; % actual mf chords from training data
                            theta=nan;
                        case 'noiseCeiling' % full data fit
                            Y_mf_hat = U(6:31,:); % overall avg. data (includes all runs)
                            theta=nan;
                    end
                    if fitScale % same scalar applied to each model
                        Y_mf_hat   = Y_mf_hat.*scaleParam;
                        d.scaleParam = scaleParam;
                    end
                    % evaluate model predictions
                    % R2:
                    rss = sum((Y_mf_test(:)-Y_mf_hat(:)).^2);
                    tss = sum(Y_mf_test(:).*Y_mf_test(:));
                    d.r2 = 1-(rss/tss);
                    d.rss=rss;
                    d.tss=tss;
                    % R:
                    ssC = sum(Y_mf_test(:).*Y_mf_hat(:));
                    ss1 = sum(Y_mf_test(:).*Y_mf_test(:));
                    ss2 = sum(Y_mf_hat(:).*Y_mf_hat(:));
                    d.r = ssC/sqrt(ss1.*ss2);
                    d.ssC=ssC;
                    d.ss1=ss1;
                    d.ss2=ss2;
                    % info fields
                    d.model = mm;
                    d.modelName = {model{mm}};
                    d.theta = {theta};
                    d.fold  = p;
                    D=addstruct(D,d);
                end
            end
        end
        varargout = {D};
        
    case '0' % ------------ PCM: pcm analyses. ----------------------------
    case 'PCM_getData'
        % Get betas for roi from subjects in PCM-friendly format.
        % Betas do not have run means removed.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        rmvMean=0; % remove run means?
        posJustify=0; % make all beta values positive?
        vararginoptions(varargin,{'sn','glm','roi','rmvMean','posJustify'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        % load betas
        betaType = 'betaW'; % multivariately prewhitened (or betaUW, raw_beta)
        B = load(fullfile(regDir,sprintf('glm%d_roi%d_betas.mat',glm,roi)));
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
            % remove run means from each run?
            if rmvMean
                C0 = indicatorMatrix('identity',bb.run); % run-mean centring matrix
                bb.betas = bb.betas - C0*pinv(C0)*bb.betas;
            end
            % positively justify patterns? (set mininum value to zero):
            if posJustify
                bb.betas = bb.betas - min(min(bb.betas));
            end
            % put subj data into pcm variables
            Y{ii}         = bb.betas;
            partVec{ii}   = bb.run;
            condVec{ii}   = bb.chord;
            G_hat(:,:,ii) = pcm_estGCrossval(Y{ii},partVec{ii},condVec{ii});
        end
        varargout = {Y,partVec,condVec,G_hat}; 
    
    case 'PCM_fitGroup'
        % Does group-level pcm fitting for multiple rois
        sn  = pp1_imana('getSubjs');
        roi = [1:9];%[1:12,21:23];
        glm = 4;
        rmvMean=0; % remove run means? If so, run effect becomes fixed
        posJustify=0; % make all beta values positive?
        vararginoptions(varargin,{'sn','roi','glm','rmvMean','posJustify'});
        fprintf('PCM GROUP fitting.\nsubjs included in model fitting : ');
        for s=sn
            fprintf('%d.',s);
        end
        for r = roi
            fprintf('\nhemi: %s  |  roi: %s\n',hem{regSide(r)},regname{r-(regSide(r)-1)*length(regname)});
            pp1_imana('PCM_fitModels_oneROI','sn',sn,'roi',r,'glm',glm,'saveit',1,'rmvMean',rmvMean,'posJustify',posJustify);
        end
    case 'PCM_fitModels_oneROI'
        % fits pcm models to data from one region
        sn     = [];
        glm    = [];
        roi    = [];
        saveit = [];        % save fits
        rmvMean=0; % remove run means? If so, run effect becomes fixed
        posJustify=0; % make all beta values positive?
        vararginoptions(varargin,{'sn','glm','roi','saveit','rmvMean','posJustify'});

        [Y,partVec,condVec,G_hat_all] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm,...
            'rmvMean',rmvMean,'posJustify',posJustify); % get data
        G_hat = mean(G_hat_all,3);
        
        % get model structure
        M = pp1_imana('PCM_defineModels',G_hat); 
        
        % choose proper way to deal with run effects
        runEffect = 'random';
        if rmvMean
            % if removing run means, model run effects as fixed (diffs
            % should be removed)
            runEffect = 'fixed';
        end
        % fit all models
        [T,theta,G_pred]        = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'isCheckDeriv',0,'verbose',1); % perform derivative checks (to ensure I've done correct deriv implementation for nonlinear models)
        [Tcv,theta_cv,Gcv_pred] = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta,'fitScale',1,'verbose',1,'isCheckDeriv',0);
        
        % do group CV gainExp model fittings:
        [M,T,theta,G_pred,Tcv,Gcv_pred] = pp1_imana('PCM_cvPowerScalingFits',Y,partVec,condVec,G_hat(1:5,1:5),M,T,theta,G_pred,Tcv,Gcv_pred);
        
        % save fits?
        if saveit
           outfile = fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d_%d%d',glm,roi,rmvMean,posJustify));
           save(outfile,'M','T','theta','Tcv','theta_cv','G_hat_all','G_pred','Gcv_pred'); 
        end
        
        %keyboard
        varargout = {Tcv,T,M,theta_cv,G_pred,Y,partVec,condVec};
    case 'PCM_defineModels'
        % case to define models fitted with PCM
        G_hat = varargin{1}; % [31,31] second moment matrix
        G_hat = mean(G_hat,3); % ensure it is avg. across subjects
        % define model params and starting values
        G_sf        = G_hat(1:5,1:5); % group avg. single-finger G
        Asf         = pcm_diagonalize(G_sf); 
        chords      = pp1_imana('chords');
        scaleParams = log([0.8 0.6 0.4 0.2])';
        
        % 1. null model
        M{1} = pp1_imana('pcm_null');
        M{end}.name = 'null I';
        % - - - - -
        % 2. nonlinear null model
        M{end+1} = pp1_imana('pcm_nullNonlinear');
        M{end}.name = 'null NL';
        % - - - - -
        % 3. linear scaling of avg. single-finger G
        M{end+1}  = pp1_imana('pcm_linearScale_FixedG');
        M{end}.Gc = chords*(Asf*Asf')*chords';
        M{end}.name   = 'L scale Gsf';
        % - - - - -
        % 4. non-linear scaling of avg. single-finger G (scalar per #fingers)
        M{end+1}      = pp1_imana('pcm_nonlinearScale');
        M{end}.theta0 = scaleParams;
        M{end}.Ac     = Asf;
        M{end}.name   = 'NL scale Gsf';
        % - - - - -
        % 5-n. nonlinear scaling null model
        M = pp1_imana('pcm_powerScaling',M,G_sf,linspace(0.02,1,50));
        % - - - - -
        % n+1. noiseceiling (freedirect)
        M{end+1} = pp1_imana('pcm_freedirect'); 
        % results are nearly identical to freechol, but faster estimation
        
        varargout = {M};
        
%         % 6. chord usage G (predicted from 'ef1_gloveData.m')
%         M{end+1} = pp1_imana('pcm_componentUsageChords');
%         % - - - - -
%         % 7. chord usage G + non-linear linear model
%         M{end+1} = pp1_imana('pcm_nonlinearScale_plusUsageChords');
%         M{end}.theta0 = [scaleParams;0];
%         M{end}.Gc(1:5,1:5,1)=G_hat(1:5,1:5);
    case 'PCM_cvPowerScalingFits'
        % case to make and do CV fits of the PowerScaling model.
        % Specifically, we will calculate the avg. group gain value (leaving each
        % subject out).
        % Then, we fit these models, and append the results to the current
        % model structures.
        % Finally, we return these model structures so that they can be
        % saved.
        % This case is a work in progress and so is pretty ugly.
        
        % deal with inputs:
        Y       = varargin{1}; % subject beta structure
        partVec = varargin{2}; % partition vector structure
        condVec = varargin{3}; % condition vector structure
        G_sf    = varargin{4}; % group-avg. single-finger second moment matrix
        M       = varargin{5}; % pcm model structure
        T       = varargin{6}; % overall pcm fits
        theta   = varargin{7}; % overall theta params
        G_pred  = varargin{8}; % overall predicted Gs
        Tcv     = varargin{9}; % crossvalidated pcm fits
        Gcv_pred = varargin{10};% crossvalidated predicted Gs
        
        % find the gain exponent models (loop through models and check
        % names):
        midx = cellfun(@(x) startsWith(x.name,'power'),M,'uni',1);
        % determine best fitting gain model (based on cv loglikelihoods):
        logLike = Tcv.likelihood(:,midx);
        [~,idx] = max(logLike,[],2);
        Mtmp = {};
        for ii=find(midx)
            Mtmp{end+1} = M{ii};
        end
        % use the best gains from each subject to create CV gain models:
        bestParams = cellfun(@(x) x.powerParam,Mtmp,'uni',1);
        bestParams = bestParams(idx);
        Mg = pp1_imana('pcm_powerScaling',{},G_sf,bestParams(1));
        for ii=2:numel(Y)
            tmp = pp1_imana('pcm_powerScaling',{},G_sf,bestParams(ii));
            Mg{1}.powerParam(end+1) = bestParams(ii);
            Mg{1}.Gc(:,:,end+1) = tmp{1}.Gc;
        end
        Mg{1}.name = 'power CV';
        Mg{1}.fitAlgorithm = 'minimize';
        % now fit this model:
        [Tge,theta_ge,Gge_pred] = pcm_fitModelGroup(Y,Mg,partVec,condVec,'runEffect','random','fitScale',1,'isCheckDeriv',0,'verbose',1); % perform derivative checks (to ensure I've done correct deriv implementation for nonlinear models)
        [Tge_cv,~,Gge_cv_pred]  = pcm_fitModelGroupCrossval(Y,Mg,partVec,condVec,'runEffect','random','fitScale',1,'verbose',1,'isCheckDeriv',0);
        % finally, append this model info accordingly (will be second last model)
        M{end+1} = M{end};
        M{end-1} = Mg{1};
        theta{end+1} = theta{end};
        theta{end-1} = theta_ge{1};
        G_pred{end+1} = G_pred{end};
        G_pred{end-1} = Gge_pred{1};
        Gcv_pred{end+1} = Gcv_pred{end};
        Gcv_pred{end-1} = Gge_cv_pred{1};
        fields = {'iterations','time','likelihood','noise','scale','run'};
        for ii = 1:numel(fields)
            T.(fields{ii})(:,end+1) = T.(fields{ii})(:,end);
            T.(fields{ii})(:,end-1) = Tge.(fields{ii});
            Tcv.(fields{ii})(:,end+1) = Tcv.(fields{ii})(:,end);
            Tcv.(fields{ii})(:,end-1) = Tge_cv.(fields{ii});
        end
        T.reg(:,end+1) = T.reg(:,end);
        T.reg(:,end-1) = 0;
        varargout = {M,T,theta,G_pred,Tcv,Gcv_pred};
    
    case 'PCM_fitIndivid'  
        % Does individual-level pcm fitting for multiple rois
        sn  = pp1_imana('getSubjs');
        roi = [1:4,9,7];
        glm = 4;
        rmvMean=0; % remove run means? If so, run effect becomes fixed
        posJustify=0; % make all beta values positive?
        vararginoptions(varargin,{'sn','roi','glm','rmvMean','posJustify'});
        fprintf('PCM INDIVID CROSSVAL fitting.\nsubjs included in model fitting : ');
        for s=sn
            fprintf('%d.',s);
        end
        for r = roi
            fprintf('\nhemi: %s  |  roi: %s\n',hem{regSide(r)},regname{r-(regSide(r)-1)*length(regname)});
            pp1_imana('PCM_fitModelsIndivid_oneROI','sn',sn,'roi',r,'glm',glm,'saveit',1,'rmvMean',rmvMean,'posJustify',posJustify);
        end    
    case 'PCM_fitModelsIndivid_oneROI'
        % fits pcm models to data from one region
        sn     = [];
        glm    = [];
        roi    = [];
        saveit = [];        % save fits
        rmvMean=0; % remove run means? If so, run effect becomes fixed
        posJustify=0; % make all beta values positive?
        vararginoptions(varargin,{'sn','glm','roi','saveit','rmvMean','posJustify'});

        [Y,partVec,condVec,G_hat_all] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm,...
            'rmvMean',rmvMean,'posJustify',posJustify); % get data
        
        % choose proper way to deal with run effects
        runEffect = 'random';
        if rmvMean
            % if removing run means, model run effects as fixed (diffs
            % should be removed)
            runEffect = 'fixed';
        end
        % define output structures:
        T      = [];
        Tcv    = [];
        DD     = [];
        Params = [];
        G_pred = cell(1,55);
        % fit all models
        % since fitting multiple fixed models, we need to define models for each subject, and fit them separately:
        for ii=1:numel(Y) % for each subj
            fprintf('s%02d\n',sn(ii));
            % get run-averaged single finger patterns for this subject:
            C0 = pcm_indicatorMatrix('identity',condVec{ii});
            Yavg = pinv(C0)*Y{ii};
            Ysf  = Yavg(1:5,:);
            % test model fits using only multi-finger conditions:
            testIdx = condVec{ii}>5;
            Ytest = Y{ii}(testIdx,:);
            pV    = partVec{ii}(testIdx);
            cV    = condVec{ii}(testIdx);
            % get model structure (models defined for each participant separately because use their own G)
            M = pp1_imana('PCM_defineModelsIndivid',Ysf); 
            % overall model fitting:
            [t,theta_subj,subj_G_pred] = pcm_fitModelIndivid({Ytest},M,{pV},{cV},'runEffect',runEffect,'verbose',0);
            % crossval model fitting:
            [tcv,dd,theta_cv] = pcm_fitModelIndividCrossval({Ytest},M,{pV},{cV},...
                'evaluation',{'likelihood','R2','R'},'crossvalScheme','oddEven',...
                'runEffect',runEffect,'verbose',1);
            % note correct subect number:
            t.SN   = sn(ii);
            tcv.SN = sn(ii);
            dd.SN = ones(size(dd.SN)).*tcv.SN;
            % add fits and info to output variables
            T   = addstruct(T,t);
            Tcv = addstruct(Tcv,tcv);
            DD  = addstruct(DD,dd);
            for mm=1:numel(M) % per model
                params.SN    = sn(ii);
                params.model = mm;
                params.theta = {theta_subj{mm}'}; % make cells to deal with models with diff numParams
                params.thetaCV = {theta_cv{mm}'};
                Params = addstruct(Params,params);
                G_pred{mm}(:,:,ii) = subj_G_pred{mm};
            end
        end
        
        % save fits?
        if saveit
           outfile = fullfile(pcmDir,sprintf('pcmFitsIndivid_glm%d_roi%d_%d%d',glm,roi,rmvMean,posJustify));
           %save(outfile,'M','T','theta','G_hat_all','G_pred','INFO','DD','Tcv','theta_cv'); 
           save(outfile,'M','G_hat_all','G_pred','DD','T','Tcv','Params'); 
        end
        
        %keyboard
        %varargout = {T,M,G_pred,INFO,Y,partVec,condVec};
        varargout = {Tcv,M,DD,Y,partVec,condVec};
    case 'PCM_defineModelsIndivid'
        % case to define models fitted with PCM
        Ysf   = varargin{1}; % subject-specific single finger activity patterns (avg. across runs) [5xP]
        
        % 1. null model
        M{1} = pp1_imana('pcm_null');
        M{1}.Gc = eye(26);
        M{end}.name = 'null I';
        % - - - - -
        % 2. nonlinear null model
        M{end+1} = pp1_imana('pcm_nullNonlinear');
        M{end}.numCond = 26;
        M{end}.name = 'null NL';
        % - - - - -
        % 3. linear scaling of avg. single-finger G
        M{end+1}  = pp1_imana('pcm_linearScale_FixedG');
        M{end}.Gc = pp1_imana('pcm_subjectModel_Glinear',Ysf);
        M{end}.name   = 'L scale Gsf';
        M{end}.fitAlgorithm = 'minimize';
        % - - - - -
        % 4-n. nonlinear scaling null model on patterns
        M = pp1_imana('pcm_subjectModel_PowerScaling',M,Ysf,[0.02:0.02:0.98]);
        % - - - - -
        % n+1. noiseceiling (freedirect)
        M{end+1} = pp1_imana('pcm_freedirect'); 
        M{end}.name = 'noiseceiling_cv';
        M{end+1} = M{end};
        M{end}.name = 'noiseceiling_overall';
        % - - - -
        % fancy noise ceiling:
%         M{end+1}.type   = 'freechol';
%         M{end}.numCond  = 31;
%         M{end}.name     = 'noiseceiling_freechol';
%         M{end}          = pcm_prepFreeModel(M{end});
% - - - - -
%         % 4. non-linear scaling of avg. single-finger G (scalar per #fingers)
%         M{end+1}      = pp1_imana('pcm_nonlinearScale');
%         M{end}.theta0 = scaleParams;
%         M{end}.Ac     = Asf;
%         M{end}.name   = 'NL scale Gsf';
%         M{end}.fitAlgorithm = 'minimize';
        
        % results are nearly identical to freechol, but faster estimation
        
        varargout = {M};    
    
    case 'usageGsf'
        U = load('/Users/sarbuckle/DATA/passivePatterns1/fmri/PCM_models/usageGsf.mat');
        varargout = {U.G,U.G_cent};
    case 'usageGchord'
        % made from ef1_gloveData('calcDistances_chords')
        d = [];
        usageGchords = [];
        load('/Users/sarbuckle/DATA/passivePatterns1/fmri/PCM_models/usageGchord.mat');
        varargout = {usageGchords};
    case 'pcm_null'
        % Model Null: chords patterns are fully independent
        M.type       = 'component';
        M.numGparams = 1;
        %M.theta0     = [];
        M.name       = 'null';
        M.Gc(:,:,1)  = eye(31);
        varargout = {M};
    case 'pcm_nullNonlinear'
        % nonlinear null model where all distances are equal, but
        % conditions are not fully independent
        M.type         = 'nonlinear'; 
        M.modelpred    = @pp1_modelpred_null;
        M.fitAlgorithm = 'minimize'; 
        M.numGparams   = 1;
        M.theta0       = 0.3;
        M.name         = 'null nonlinear';
        M.numCond      = 31;
        varargout = {M};
    case 'pcm_null_nonlinearScale'
        % nonlinear null model where all chords of the same # fingers are
        % equally distinct, and chords with more fingers are more distinct.
        % This model attempts to fit overall activity scaling.
        M.type         = 'nonlinear'; 
        M.modelpred    = @pp1_modelpred_nonlinearScale_fixedFinger;
        M.fitAlgorithm = 'minimize'; 
        M.numGparams = 4;
        M.Ac(:,:,1)  = [];
        M.name         = 'null nonlinear scale';
        varargout = {M};
    case 'pcm_linearScale_FixedG'
        % Model Linear: chord patterns are linear summations of single
        % fingers (independent)
        M.type       = 'fixed';
        M.numGparams = 0;
        M.theta0     = [];
        M.name       = 'L scale single finger G';
%         M.Ac(:,:,1)  = pp1_simulations('chords');
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
    case 'pcm_nonlinearScale_usageG'
        % Model Linear: chord patterns are linear summations of single
        % finger patterns scaled by number of fingers pressed.
        % Estimates both scaling params with usage finger model.
        M.type       = 'nonlinear';
        M.name       = 'NL scale Usage fingers';
        M.modelpred  = @pp1_modelpred_nonlinearScale_fixedFinger;
        M.numGparams = 4;
        %M.theta0     = [];
        M.Gc(:,:,1)  = pp1_imana('usageGsf');
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
    case 'pcm_componentUsageChords'
        % Model component: second moment is estimated from glove data.
        M.type       = 'fixed'; % allow for scaling factor
        M.name       = 'chord usage';
        M.numGparams = 0;
        M.theta0     = [];
        M.Gc(:,:,1)  = pp1_imana('usageGchord');
        varargout = {M};
    case 'pcm_nonlinearScale_plusUsageChords'
        % Model nonLinear plus Chord usage: chord patterns are function of:
        % 1. linear summations of single fingers scaled by number of fingers pressed.
        % 2. plus some weighted amount of the chord usage second moment
        % (from 'ef1_gloveData.m')
        % Estimates ONLY scaling and component params
        M.type       = 'nonlinear';
        M.name       = 'NL scale + Usage';
        M.modelpred  = @pp1_modelpred_nonlinearScale_plusUsage;
        M.numGparams = 5;
        M.Gc(:,:,1)  = zeros(31); % leave empty for user to enter single-finger G
        M.Gc(:,:,2)  = pp1_imana('usageGchord');
        M.theta0     = [];
        varargout = {M};
    case 'pcm_powerScaling'
        % case to define a series of fixed pos def models, where G is
        % transformed by theta*G.^b
        M     = varargin{1}; % pcm model structure (we append the fixed models to current structure)
        Gsf   = varargin{2}; % single finger G (we make semi-PD in this case)
        b_vec = varargin{3}; % assumes b is row vector
        chords= pp1_imana('chords'); % indicator matrix for chords
        Asf   = pcm_diagonalize(Gsf); % make Gsf semi-pd
        Gmf   = chords*(Asf*Asf')*chords';
        signGmf = sign(Gmf);
        Gmf_abs = abs(Gmf); % hacky solution to scaling negative covariances (remove sign, then add sign back in later)
        numModels = numel(M);
        for ii = 1:numel(b_vec)
            midx = numModels+ii;
            M{midx}.type       = 'fixed';
            M{midx}.numGparams = 0;
            %M{midx}.theta0     = [];
            M{midx}.Gc         = (Gmf_abs.^b_vec(ii)) .* signGmf; % this model is not crossvalidated (it's group avg. single finger G)
            M{midx}.name       = sprintf('power %1.2f',b_vec(ii));
            M{midx}.powerParam = b_vec(ii);
        end
        varargout = {M};
    case 'pcm_freedirect'
        % Naive averaring model- noise ceiling method 1- totall free model
        M.type       = 'freedirect';  
        M.numGparams = 0;
        M.theta0     = [];
        M.name       = 'fd noiseceiling';
        varargout = {M};
    case 'pcm_freechol'
        M.type       = 'freechol';  
        M.numCond    = 31;
        M.name       = 'fc noiseceiling';
        M            = pcm_prepFreeModel(M);
        varargout = {M};
 
    case 'pcm_subjectModel_PowerScaling'
        % case to define a series of fixed pos def models, where model G is
        % second moment of M*U.^b * U'.^b * M'
        M     = varargin{1}; % pcm model structure (we append the fixed models to current structure)
        Ysf   = varargin{2}; % single finger patterns (avg. across runs) for one subj
        b_vec = varargin{3}; % assumes b is row vector
        chords= pp1_imana('chords'); % indicator matrix for chords
        
        if iscell(Ysf)
            error('can only define fixed model for one subject at a time!')
        end
        numModels = numel(M);
        for ii = 1:numel(b_vec)
            midx = numModels+ii;
            M{midx}.type       = 'fixed';
            M{midx}.fitAlgorithm = 'minimize';
            M{midx}.numGparams = 0;%1;
           % M{midx}.theta0     = [];
            % first do linear model predictions:
            X       = chords(6:31,:);
            Ylinear = X*Ysf; %linear chord prediction per run
            % now do nonlinear scaling:
            % for negative values, they stay negative values!
            signY = sign(Ylinear);
            absY  = abs(Ylinear);
            Ymf   = (absY.^b_vec(ii)).*signY; % do nonlinearity AFTER summation of single finger patterns
            M{midx}.Gc = Ymf*Ymf' ./ size(Ysf,2);
            M{midx}.name       = sprintf('power %1.2f',b_vec(ii));
            M{midx}.powerParam = b_vec(ii);
        end
        varargout = {M};
    case 'pcm_subjectModel_Glinear'
        Ysf = varargin{1}; % [5 x P]
        chords = pp1_imana('chords');
        X      = chords(6:31,:); 
        Ymf    = X*Ysf; % predict multi-finger conditions
        Gmf    = Ymf*Ymf' ./ size(Ysf,2);
        varargout = {Gmf,Ymf};
   
    case 'PCM_plotGpredGroup'  
        % loads fit results per roi and plots them.
        glm = [];
        roi = [];
        model = [];
        vararginoptions(varargin,{'glm','roi','model'});
        if length(roi)>1
           error('can only call case with one roi'); 
        end
        % load fits
        load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d_00.mat',glm,roi)));
        % plot
        numPlots = numel(model);
        j = 1;
        for ii = model
            % plot fits
            %subplot(2,ceil(numPlots/2),j);
            subplot(1,numPlots,j);
            modelGpred = G_pred{ii}.*mean(T.scale(:,ii)); % rescale G accordingly to account for shifts in SnR;
            imagesc(modelGpred); 
            %title(sprintf('%s : %d params',M{i}.name,M{i}.numGparams));
            title(sprintf('GROUP %s : %s',[regname{roi-(regSide(roi)-1)*length(regname)} ' (' hem{regSide(roi)} ')'],M{ii}.name));
            axis square,
            colorbar
            % drawlines for chord distinctions
            drawline(5.5,'dir','horz');
            drawline(5.5,'dir','vert');
            drawline(15.5,'dir','horz');
            drawline(15.5,'dir','vert');
            drawline(25.5,'dir','horz');
            drawline(25.5,'dir','vert');
            drawline(30.5,'dir','horz');
            drawline(30.5,'dir','vert');
            j = j+1;
        end 
        varargout = {G_pred};     
    case 'PCM_plotFitsGroup'
        type = varargin{1}; % bar, box, or r2
        glm = 4;
        roi = [1:6];
        nNull = 1;
        nCeil = 56;
        modelsToPlot = [1,3:4,56];%[nNull:4,nCeil];
        
        
        inputs = {'glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil,'modelsToPlot',modelsToPlot};
        switch type
            case 'bar'
                D=pp1_imana('PCM_plotFitsBar',inputs{:});
            case 'box'
                D=pp1_imana('PCM_plotFitsBox',inputs{:});
            case 'r2'
                D=pp1_imana('PCM_plotPseudoR2',inputs{:});
        end
        varargout={D};
    case 'PCM_getFits'
        % Gets models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        glm   = [];
        roi   = [];
        nNull = []; % which model is null?
        nCeil = []; % which model is noise ceiling?
        vararginoptions(varargin,{'glm','roi','nNull','nCeil'});
        D   = []; % output structure
        
        for r = roi
            Tcv = [];
            % if exists, load pcm fits for region (otherwise, skip region)
            fprintf('\nroi %d...',r);
            try
                load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d_00.mat',glm,r)));
            catch
                fprintf('no file.');
                continue
            end
            % scale crossvalidated likelihoods
            Tcv.likelihood_norm = bsxfun(@minus,Tcv.likelihood,Tcv.likelihood(:,nNull));
            % arrange into plotting structure
            numSubjs   = size(Tcv.SN,1);
            numModels  = numel(M);
            nameModels = {};
            Q = [];
            for m = 1:numModels
                % get model names
                nameModels{end+1,1} = M{m}.name;
                % get thetas (if model has any)
                if M{m}.numGparams==0
                    q.thetaCV = num2cell(nan(numSubjs,1),2);
                else
                    q.thetaCV = num2cell(theta_cv{m}',2);
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
        glm = [];
        roi = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        vararginoptions(varargin,{'glm','roi','nNull','nCeil','modelsToPlot'});
        % get pcm fits
        D = pp1_imana('PCM_getFits','glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.model,modelsToPlot));
        % calc pseudo R2
        upperNoiseCeil = kron(D.likeNorm(D.model==nCeil),ones(numel(unique(D.model)),1));
        D.pseudoR2 = D.likeNormCV./upperNoiseCeil;
        % check: [D.roi D.sn upperNoiseCeil]
        % % each subject should be scaled by a different upper noise ceiling model per roi
        
        % plot pcm fits
        numModels  = numel(modelsToPlot);
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        sty = style.custom(plt.helper.get_shades(numModels-2,'hot','descend'));
        sty.general.markersize=5;%3;
        for ii = 1:numPlots
            r = roi(ii);
            roiName = [regname{r-(regSide(r)-1)*length(regname)} ' (' hem{regSide(r)} ')'];
            subplot(1,numPlots,ii);
            plt.box(D.model,D.pseudoR2,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty,'plotall',1);
            %plt.box(D.model,D.pseudoR2,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
%             xticks=get(gca,'xtick');
%             xlabels=linspace(0.02,1,50);
%             plt.set('xtick',xticks(1:3:end),'xticklabel',xlabels(1:3:end),'xticklabelrotation',45);
            plt.labels('','pseudo model R2',roiName);
            % plot noise ceilings
            drawline(median(D.pseudoR2(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-.'); % lower noise ceiling
            drawline(1,'dir','horz','linestyle','-');
            drawline(0,'dir','horz','linestyle','-','color',[0.7 0.7 0.7]);
            legend off
            ylim([0 1]);
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotFitsBox'  
        % loads fit results per roi and plots them.
        glm = [];
        roi = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        vararginoptions(varargin,{'glm','roi','nNull','nCeil','modelsToPlot'});
        % get pcm fits
        D = pp1_imana('PCM_getFits','glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.model,modelsToPlot));
        % plot pcm fits
        numModels  = length(unique(D.model));
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        for i = 1:numPlots
            r = roi(i);
            subplot(1,numPlots,i);
            sty = style.custom(plt.helper.get_shades(numModels-2,'jet','descend'));
            plt.box(D.model,D.likeNormCV,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
            plt.labels('','relative log likelihood',[regname{r-(regSide(r)-1)*length(regname)} ' (' hem{regSide(r)} ')']);
            % plot noise ceilings
            drawline(median(D.likeNormCV(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-.');
            drawline(median(D.likeNorm(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-');
            legend off
            ylims = ylim;
            %ylim([0 ylims(2)]);
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotFitsBar'  
        % loads fit results per roi and plots them.
        glm = [];
        roi = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        vararginoptions(varargin,{'glm','roi','nNull','nCeil','modelsToPlot'});
        numPlots = numel(roi);
        figure('Color',[1 1 1]);
        for i = 1:numPlots
            r = roi(i);
            load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d',glm,r)));
            %T   = getrow(T,ismember(T.SN,sn));
            %Tcv = getrow(Tcv,ismember(Tcv.SN,sn));
            
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
            title([regname{r-(regSide(r)-1)*length(regname)} ' (' hem{regSide(r)} ')']);
            box off;
        end
        plt.match('y');
        varargout = {Tpcv};
    
    case 'PCM_plotGpredIndivid'  
        % loads fit results per roi and plots them.
        glm = [];
        roi = [];
        model = [];
        vararginoptions(varargin,{'glm','roi','model'});
        if length(roi)>1
           error('can only call case with one roi'); 
        end
        % load fits
        load(fullfile(pcmDir,sprintf('pcmFitsIndivid_glm%d_roi%d_00.mat',glm,roi)));
        % plot
        numPlots = numel(model);
        j = 1;
        for ii = model
            % plot fits
            %subplot(2,ceil(numPlots/2),j);
            subplot(1,numPlots,j);
            modelGpred = mean(G_pred{ii},3); % rescale G accordingly to account for shifts in SnR;
            imagesc(modelGpred); 
            %title(sprintf('%s : %d params',M{i}.name,M{i}.numGparams));
            title(sprintf('INDIVID %s : %s',[regname{roi-(regSide(roi)-1)*length(regname)} ' (' hem{regSide(roi)} ')'],M{ii}.name));
            axis square,
            colorbar
            % drawlines for chord distinctions
            drawline(5.5,'dir','horz');
            drawline(5.5,'dir','vert');
            drawline(15.5,'dir','horz');
            drawline(15.5,'dir','vert');
            drawline(25.5,'dir','horz');
            drawline(25.5,'dir','vert');
            drawline(30.5,'dir','horz');
            drawline(30.5,'dir','vert');
            j = j+1;
        end 
        varargout = {G_pred};      
    case 'PCM_plotFitsIndivid'    
        glm = 4;
        roi = [2:4,7,9];
        nNull = 2;
        nCeil = 54;
        modelsToPlot = [2,3,8,13,18,23,28,33,38,43,48,53,54];%[nNull:nCeil];
        metric = 'R'; % which evaluation metric do we use (R, R2, like)
        type = 'rel'; % plot normalized ('norm'), relative ('rel'), or raw ('raw') metrics
        vararginoptions(varargin,{'glm','roi','metric','type'});
        inputs = {'glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil,'modelsToPlot',modelsToPlot,'metric',metric};
        switch type
            case 'raw'
                D=pp1_imana('PCM_plotIndividRawFits',inputs{:});
            case 'rel'
                D=pp1_imana('PCM_plotIndividRelativeFits',inputs{:});
            case 'norm'
                D=pp1_imana('PCM_plotIndividNormFits',inputs{:});
        end
        varargout={D};
    case 'PCM_getFitsIndivid'
        % Gets individual models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        glm   = [];
        roi   = [];
        nNull = []; % which model is null?
        nCeil = []; % which model is noise ceiling?
        vararginoptions(varargin,{'glm','roi','nNull','nCeil'});
        D   = []; % output structure
        
        for r = roi
            % if exists, load pcm fits for region (otherwise, skip region)
            fprintf('\nroi %d...',r);
            try
                load(fullfile(pcmDir,sprintf('pcmFitsIndivid_glm%d_roi%d_00.mat',glm,r)));
            catch
                fprintf('no file.');
                continue
            end
%             % scale crossvalidated evaluation metrics
%             Tcv.R_norm = bsxfun(@minus,Tcv.R,Tcv.R(:,nNull)); % set null model to 0
%             Tcv.R_norm = bsxfun(@rdivide,Tcv.R_norm,Tcv.R_norm(:,nCeil)); % set upper noise ceiling to 1
%             Tcv.R2_norm = bsxfun(@minus,Tcv.R2,Tcv.R2(:,nNull));
%             Tcv.R2_norm = bsxfun(@rdivide,Tcv.R2_norm,Tcv.R2_norm(:,nCeil));
%             Tcv.likeNorm = bsxfun(@minus,Tcv.likelihood,Tcv.likelihood(:,nNull));
%             Tcv.likeNorm = bsxfun(@rdivide,Tcv.likeNorm,Tcv.likeNorm(:,nCeil));
            % arrange into plotting structure
            numSubjs   = size(Tcv.SN,1);
            numModels  = numel(M);
            nameModels = {};
            for m = 1:numModels
                % get model names
                nameModels{end+1,1} = M{m}.name;
            end
            v = ones(numModels,1);
            for j = 1:numSubjs
                d.sn  = v.*Tcv.SN(j);
                d.roi = v.*r;
                d.model = [1:numModels]';
                d.modelName  = nameModels;
                % raw evaluation values:
                d.R2CV    = Tcv.R2(j,:)';
                d.RCV     = Tcv.R(j,:)';
                d.likeCV  = Tcv.likelihood(j,:)';
                % relative evaluation values:
%                 d.R2NormCV = Tcv.R2_norm(j,:)';
%                 d.RNormCV  = Tcv.R_norm(j,:)';
%                 d.likeNormCV = Tcv.likeNorm(j,:)';
                D = addstruct(D,d);
            end
            fprintf('done.');
        end
        fprintf('\n');
        varargout = {D};
    case 'PCM_plotIndividNormFits'
        % plots pseudo R2 value (0 = null, 1 = upper noise ceiling)
        glm = [];
        roi = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        metric = []; % R, R2, or like
        vararginoptions(varargin,{'glm','roi','nNull','nCeil','modelsToPlot','metric'});
        % get pcm fits
        D = pp1_imana('PCM_getFitsIndivid','glm',glm,'roi',roi);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.model,modelsToPlot));
        % calc normalized fit (lower bound = nNull, upper = nCeil)
        switch metric
            case 'R'
                D.normFit = D.RCV - kron(D.RCV(D.model==nNull),ones(numel(unique(D.model)),1)); % subtract null
            case 'R2'
                D.normFit = D.R2CV - kron(D.R2CV(D.model==nNull),ones(numel(unique(D.model)),1));
            case 'like'
                D.normFit = D.likeCV - kron(D.likeCV(D.model==nNull),ones(numel(unique(D.model)),1));
        end
        upperNoiseCeil = kron(D.normFit(D.model==nCeil),ones(numel(unique(D.model)),1));
        D.normFit      = D.normFit./upperNoiseCeil; % normalize to upper noise ceiling
        % check: [D.roi D.sn upperNoiseCeil]
        % % each subject should be scaled by a different upper noise ceiling model per roi
        
        % plot pcm fits
        numModels  = numel(modelsToPlot);
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        sty = style.custom(plt.helper.get_shades(numModels-2,'jet','descend'));
        sty.general.markersize=5;%3;
        for ii = 1:numPlots
            r = roi(ii);
            roiName = [regname{r-(regSide(r)-1)*length(regname)} ' (' hem{regSide(r)} ')'];
            subplot(1,numPlots,ii);
            plt.box(D.model,D.normFit,'subset',D.model~=nNull & D.model~=nCeil & D.model~=nCeil-1 & D.roi==r,'split',D.model,'style',sty,'plotall',1);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
%             xticks=get(gca,'xtick');
%             xlabels=linspace(0.02,1,50);
%             plt.set('xtick',xticks(1:3:end),'xticklabel',xlabels(1:3:end),'xticklabelrotation',45);
            plt.labels('',sprintf('%s (noralized)',metric),roiName);
            % plot noise ceilings
            lowerNoiseCeil = mean(D.normFit(D.model==nCeil-1 & D.roi==r));
            drawline(lowerNoiseCeil,'dir','horz','linestyle','-.'); % lower noise ceiling
            drawline(1,'dir','horz','linestyle','-'); % upper bound
            drawline(0,'dir','horz','linestyle','-','color',[0.7 0.7 0.7]); % null bound
            legend off
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotIndividRelativeFits'
        % plots evaulation criterions relative to null model (0 = null)
        glm = [];
        roi = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        metric = []; % R, R2, or like
        vararginoptions(varargin,{'glm','roi','nNull','nCeil','modelsToPlot','metric'});
        % get pcm fits
        D = pp1_imana('PCM_getFitsIndivid','glm',glm,'roi',roi);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.model,modelsToPlot));
        % calc normalized fit (lower bound = nNull, upper = nCeil)
        switch metric
            case 'R'
                D.relFit = D.RCV - kron(D.RCV(D.model==nNull),ones(numel(unique(D.model)),1));
            case 'R2'
                D.relFit = D.R2CV - kron(D.R2CV(D.model==nNull),ones(numel(unique(D.model)),1));
            case 'like'
                D.relFit = D.likeCV - kron(D.likeCV(D.model==nNull),ones(numel(unique(D.model)),1));
        end
        % plot pcm fits
        numModels  = numel(modelsToPlot);
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        sty = style.custom(plt.helper.get_shades(numModels-1,'jet','descend'));
        sty.general.markersize=5;%3;
        for ii = 1:numPlots
            r = roi(ii);
            roiName = [regname{r-(regSide(r)-1)*length(regname)} ' (' hem{regSide(r)} ')'];
            subplot(1,numPlots,ii);
            plt.box(D.model,D.relFit,'subset',D.model~=nNull & D.roi==r,'split',D.model,'style',sty,'plotall',1);
            plt.set('xticklabel',{nameModels{2:numModels}},'xticklabelrotation',45);
%             xticks=get(gca,'xtick');
%             xlabels=linspace(0.02,1,50);
%             plt.set('xtick',xticks(1:3:end),'xticklabel',xlabels(1:3:end),'xticklabelrotation',45);
            plt.labels('',sprintf('%s (relative)',metric),roiName);
            % plot noise ceilings
            upperNoiseCeil = mean(D.relFit(D.model==nCeil & D.roi==r));
            lowerNoiseCeil = mean(D.relFit(D.model==nCeil-1 & D.roi==r)); % assume lower bound is model prior to upper noise ceiling
            drawline(lowerNoiseCeil,'dir','horz','linestyle','-.'); % lower noise ceiling
            drawline(upperNoiseCeil,'dir','horz','linestyle','-'); % upper bound
            legend off
        end
        plt.match('y');
        varargout = {D};
    case 'PCM_plotIndividRawFits'
        % plots raw evaulation criterions
        glm = [];
        roi = [];
        nNull = [];
        nCeil = [];
        modelsToPlot = [];
        metric = []; % R, R2, or like
        vararginoptions(varargin,{'glm','roi','nNull','nCeil','modelsToPlot','metric'});
        % get pcm fits
        D = pp1_imana('PCM_getFitsIndivid','glm',glm,'roi',roi);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.model,modelsToPlot));
        % calc normalized fit (lower bound = nNull, upper = nCeil)
        switch metric
            case 'R'
                D.rawFit = D.RCV;
            case 'R2'
                D.rawFit = D.R2CV;
            case 'like'
                D.rawFit = D.likeCV;
        end
        
        % plot pcm fits
        numModels  = numel(modelsToPlot);
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        sty = style.custom(plt.helper.get_shades(numModels,'jet','descend'));
        sty.general.markersize=5;%3;
        for ii = 1:numPlots
            r = roi(ii);
            roiName = [regname{r-(regSide(r)-1)*length(regname)} ' (' hem{regSide(r)} ')'];
            subplot(1,numPlots,ii);
            plt.box(D.model,D.rawFit,'split',D.model,'style',sty,'plotall',1,'subset',D.roi==r);
            plt.set('xticklabel',{nameModels{1:numModels}},'xticklabelrotation',45);
%             xticks=get(gca,'xtick');
%             xlabels=linspace(0.02,1,50);
%             plt.set('xtick',xticks(1:3:end),'xticklabel',xlabels(1:3:end),'xticklabelrotation',45);
            plt.labels('',sprintf('%s (raw)',metric),roiName);
            % plot noise ceilings
            upperNoiseCeil = mean(D.rawFit(D.model==nCeil & D.roi==r));
            lowerNoiseCeil = mean(D.rawFit(D.model==nCeil-1 & D.roi==r)); % assume lower bound is model prior to upper noise ceiling
            nullFit        = mean(D.rawFit(D.model==nNull & D.roi==r));
            drawline(lowerNoiseCeil,'dir','horz','linestyle','-.'); % lower noise ceiling
            drawline(upperNoiseCeil,'dir','horz','linestyle','-'); % upper bound
            drawline(nullFit,'dir','horz','linestyle','-');
            legend off
        end
        plt.match('y');
        varargout = {D};
        
        
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
                   
    case 'pcm_null_rmvCM'
        % Model Null: chords patterns are fully independent
        M.type       = 'component';
        M.numGparams = 1;
        M.name       = 'null';
        M.Gc(:,:,1)  = eye(30);
        varargout = {M};
    case 'pcm_nonlinearScale_rmvCM'
        % Model nonLinear: chord patterns are linear summations of single
        % fingers (independent) scaled by number of fingers pressed.
        % Estimates ONLY scaling params- no finger params
        M.type       = 'nonlinear';
        M.name       = 'NL scale I fingers';
        M.modelpred  = @pp1_modelpred_nonlinearScale_fixedFinger_rmvCM;
        M.numGparams = 3;
        M.Gc(:,:,1)  = eye(5);
        %M.theta0     = [];
        varargout = {M};
    case 'pcm_freedirect_rmvCM'
        % Naive averaring model- noise ceiling method 1- totall free model
        M.type       = 'freedirect';  
        M.numGparams = 0;
        M.name       = 'fd noiseceiling';
        varargout = {M};
    
    case 'PCM_fit_rmvCM'
        % Does pcm fitting for multiple rois
        I   = pp1_imana('LIST_subjs');
        I   = getrow(I,I.fmri_sessions==2);   % get subjects
        sn  = unique(I.sn);
        roi = 3;
        glm = 3;
        vararginoptions(varargin,{'sn','roi','glm'});
        fprintf('chord means removed\n');
        for r = roi
            fprintf('\nroi: %d...',r);
            pp1_imana('PCM_fitModels_oneROI_rmvCM','sn',sn,'roi',r,'glm',glm,'saveit',1);
        end
    case 'PCM_fitModels_oneROI_rmvCM'
        % fits pcm models to data from one region
        sn     = [];
        glm    = [];
        roi    = [];
        saveit = 1;        % save fits
        vararginoptions(varargin,{'sn','glm','roi','saveit'});
        
        %[Y,partVec,condVec,G_hat] = pp1_imana('PCM_getData_rmvCM','sn',sn,'roi',roi,'glm',glm); % get data
        [Y,partVec,condVec,G_hat] = pp1_imana('PCM_getDataSimulated_rmvCM','sn',sn,'roi',roi,'glm',glm); % get simulated data
        G_hat = mean(G_hat,3);
        M = pp1_imana('PCM_defineModels_rmvCM',G_hat); % get models
        outfile = fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d_noChordMeans',glm,roi));

        % fit models
        runEffect = 'random';
        [T,theta_hat,G_pred] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
        [Tcv,theta_cv]       = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'verbose',1);

        % save fits?
        if saveit
           save(outfile,'M','T','theta_hat','G_pred','Tcv','theta_cv'); 
        end
        
        %keyboard
        varargout = {Tcv,T,M,theta_cv,G_pred,Y,partVec};
    case 'PCM_defineModels_rmvCM'
        % case to define models fitted with PCM
        G_hat = varargin{1}; % 30x30 group avg. second moment matrix- no 5 finger chords
        
        % define model params and starting values
        scaleParams = log([0.8 0.6 0.4])';
        
        % 1. null model
        M{1} = pp1_imana('pcm_null_rmvCM');
        % - - - - -
        % 2. non-linear scaling of avg. single-finger G
        M{end+1}      = pp1_imana('pcm_nonlinearScale_rmvCM');
        M{end}.theta0 = scaleParams;
        M{end}.Gc     = G_hat(1:5,1:5);
        M{end}.name   = 'NL scale single finger G';
        % - - - - -
        % 3. noiseceiling (freedirect)
        M{end+1} = pp1_imana('pcm_freedirect_rmvCM'); 
        % results are nearly identical to freechol, but faster estimation
        
        varargout = {M};
    case 'PCM_getData_rmvCM'
        % Get betas for roi from subjects in PCM-friendly format.
        % Betas do not have run means removed.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        numDigits = sum(pp1_imana('chords'),2);
        numDigits = numDigits(1:30);
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
            bb = getrow(bb,ismember(bb.chord,1:30)); % restrict to passive patterns only and exclude the five-finger chord (rmv numDigits mean will kill this pattern entirely)
            % remove chord means (separately per run)
            C0 = indicatorMatrix('identity',numDigits);
            C0 = kron(eye(numel(unique(bb.run))),C0);
            bb.betas = bb.betas -C0*pinv(C0)*bb.betas;
            % put subj data into pcm variables
            Y{i}         = bb.betas;
            partVec{i}   = bb.run;
            condVec{i}   = bb.chord;
            G_hat(:,:,i) = pcm_estGCrossval(Y{i},partVec{i},condVec{i});
        end
        varargout = {Y,partVec,condVec,G_hat};
    case 'PCM_getDataSimulated_rmvCM'
        % Get simulated for roi from subjects in PCM-friendly format.
        % Simulated from non-linear linear combination model.
        % This is used for sanity checks.
        sn  = [];
        glm = [];
        roi = []; % only one roi supported
        vararginoptions(varargin,{'sn','glm','roi'});
        if length(roi)>1
            error('only 1 roi supported per call to case');
        end
        numDigits = sum(pp1_imana('chords'),2);
        numDigits = numDigits(1:30);
        % load model fits (model 4)
        M = load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d.mat',glm,roi)));
        modelNum = 3;
        signal   = M.Tcv.scale(:,modelNum);
        noise    = M.Tcv.noise(:,modelNum);
        thetas   = M.theta_cv{modelNum};
        model    = M.M{modelNum};
        clear M
        numSubj  = size(thetas,2);
        numRun   = 11;
        % get num voxels for realistic simulations
        V = pp1_imana('numVoxels','sn',1:numSubj,'roi',roi);
        % outputs
        Y = {};
        partVec = {};
        condVec = {};
        for s = 1:numSubj
            % simulate subject data
            b = [];
            [b.Y,b.pV,b.cV] = pcm_generateData(model,thetas(:,s),...
                'numPart',numRun,'numVox',V.numVox(V.sn==s),'signal',signal(s),'noise',noise(s));
            % drop the five finger chord condition, and remove mean pattern
            % across chords with the same number of fingers stimulated
            b.Y = b.Y{1};
            b = getrow(b,b.cV<31);
            % remove chord means (separately per run)
            C0 = indicatorMatrix('identity',numDigits);
            C0 = kron(eye(numRun),C0);
            Y{s} = b.Y -C0*pinv(C0)*b.Y;
%             Y{s} = b.Y;
            partVec{s} = b.pV;
            condVec{s} = b.cV;
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},partVec{s},condVec{s});
        end
        varargout = {Y,partVec,condVec,G_hat};
    case 'PCM_getFits_rmvCM'
        % Gets models fits across regions & arranges into plotting
        % structure.
        % Assumes null model is model 1 and noiseceiling is last model.
        glm   = [];
        roi   = [];
        nNull = []; % which model is null?
        nCeil = []; % which model is noise ceiling?
        vararginoptions(varargin,{'glm','roi','nNull','nCeil'});
        D   = []; % output structure
        
        for r = roi
            Tcv = [];
            % if exists, load pcm fits for region (otherwise, skip region)
            fprintf('\nroi %d...',r);
            try
                load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d_noChordMeans.mat',glm,r)));
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
    case 'PCM_plotPseudoR2_rmvCM'
        % plots pseudo R2 value (0 = null, 1 = upper noise ceiling)
        I   = pp1_imana('LIST_subjs');
        I   = getrow(I,I.fmri_sessions==2);   % get subjects
        sn  = unique(I.sn);
        glm = 3;
        roi = 3;
        nNull = 1;
        nCeil = 3;
        modelsToPlot = [nNull:nCeil];
        vararginoptions(varargin,{'glm','roi','sn'});
        % get pcm fits
        D = pp1_imana('PCM_getFits_rmvCM','glm',glm,'roi',roi,'nNull',nNull,'nCeil',nCeil);
        D = getrow(D,ismember(D.roi,roi) & ismember(D.sn,sn) & ismember(D.model,modelsToPlot));
        % calc pseudo R2
        D.pseudoR2 = D.likeNormCV./kron(D.likeNorm(D.model==nCeil),ones(numel(unique(D.model)),1));
        % plot pcm fits
        numModels  = numel(modelsToPlot);
        nameModels = D.modelName(1:numModels);
        numPlots   = numel(roi);
        for i = 1:numPlots
            r = roi(i);
            subplot(1,numPlots,i);
            sty = style.custom(plt.helper.get_shades(6,'jet','descend'));
            plt.box(D.model,D.pseudoR2,'subset',D.model~=nNull & D.model~=nCeil & D.roi==r,'split',D.model,'style',sty);
            plt.set('xticklabel',{nameModels{2:numModels-1}},'xticklabelrotation',45);
            plt.labels('','pseudo model R2',regname{r});
            % plot noise ceilings
            drawline(mean(D.pseudoR2(D.model==nCeil & D.roi==r)),'dir','horz','linestyle','-.'); % lower noise ceiling
            drawline(1,'dir','horz','linestyle','-');
            legend off
            ylims = ylim;
            ylim([0 ylims(2)]);
        end
        plt.match('y');
        varargout = {D};
    
    case 'sim_divNorm'
        % simulate single finger data
        % scale data linearly
        % also scale data linearly, then apply divisive normalization
        
        % note: the params for this model are arbitrary and I fiddled with
        % them to produce an okay simulation. The simulation is not about
        % these parameters- instead, I am using this case to gain an
        % intuition about applying divisive normalization to the voxels vs.
        % the second moment.
        
        Gsf = varargin{1};
        Gsf = Gsf(1:5,1:5);
        theta = varargin{2};
        Gsf = pp1_imana('usageGchord');
        Gsf = Gsf(1:5,1:5).*0.5; % use scaled nat stats as single finger model
        % generate single finger voxel patterns
        p = 100; % #voxels
        U = mvnrnd_exact(Gsf,p);
        Z = pp1_imana('chords');
        Y = Z*U; % linearly scaled patterns, positive mean
        % do divisive norm scaling of voxel patterns, BUT DON'T APPLY TO
        % SINGLE FINGER PATTERNS
        Ymf = Y(6:end,:);
        fN  = Ymf.^theta(2);
        fD  = fN + theta(3)^theta(2);
        Yd  = [ Y(1:5,:);  real(theta(1) * (fN./fD))];
        % calculate second moments of patterns:
        Gl  = (Y*Y')./p;   % linear scaling
        Gd  = (Yd*Yd')./p; % div norm @ voxel level
        % do div norm on second moment:
        Gdg = pp1_imana('divNorm_noLatInhibition',Gl,theta);
        % plot results:
        subplot(1,2,1);
        scatter(rsa_vectorizeIPM(Gl),rsa_vectorizeIPM(Gd),12,'k');   % normalization of voxel patterns
        hold on
        scatter(rsa_vectorizeIPM(Gl),rsa_vectorizeIPM(Gdg),12,'r');  % normalization of second moment
        hold off
        xl = xlim;
        yl = ylim;
        limMin = min([min(xl) min(yl)]);
        limMax = max([max(xl) max(yl)]);
        ylim([limMin limMax]);
        xlim([limMin limMax]);
        axis square
        
        subplot(1,2,2);
        scatter(rsa_vectorizeIPM(Gd),rsa_vectorizeIPM(Gdg),12,'b','filled');
        refline(1,0);
        xlabel('voxel div norm');
        ylabel('covariance div norm');
        axis equal
        axis square
        
        %keyboard
    case 'plot_modelScaling'
        % makes scatter plot of linear vs. nonlinear scaling model
        % (co)variances.
        roi = [];
        glm = [];
        model = [];
        vararginoptions(varargin,{'roi','glm','model'});
        m = 4; % model # for linear scaling
        % get predicted gs
        load(fullfile(pcmDir,sprintf('pcmFits_glm%d_roi%d_old.mat',glm,roi)));
        g_true = rsa_vectorizeIPM(mean(G_hat_all,3));               % actual data covariances
        g_lin  = rsa_vectorizeIPM(G_pred{m} .* mean(T.scale(:,m))); % linear scaling model covariances

        % do scatterplot
        scatter(g_lin,g_true,12,'k','filled');
        ylabel('(co)variances (group avg. Gcv)');
        xlabel('(co)variances (group avg. linear model)');
        xl = xlim;
        yl = ylim;
        limMin = min([min(xl) min(yl)]);
        limMax = max([max(xl) max(yl)]);
        ylim([limMin limMax]);
        xlim([limMin limMax]);
        axis square
        
        % fit models:
        opts = optimset('fminsearch');
        opts.MaxIter = 10000;
        xp = linspace(min(g_lin),max(g_lin),20);
        switch model
            case 1 % fitCompress
                theta = [0.5,0.5];
                [theta,err,gEst] = fminsearch(@(theta) fitCompress(theta,g_lin,g_true),theta,opts);
                yp = real((theta(1)*xp).^(1/theta(2))); % make fitted line to overlay on scatterplot
            case 2 % fit exponent
                theta = 2;
                [theta,err,gEst] = fminsearch(@(theta) fitExponent(theta,g_lin,g_true),theta,opts);
                yp = real(xp.^theta(1));
            case 3 % fit gain exponent (ax^b)
                theta = [2,2];
                [theta,err,gEst] = fminsearch(@(theta) fitGainExponent(theta,g_lin,g_true),theta,opts);
                yp = real(theta(1)* (xp.^theta(2)));
            case 4 % fit divisive norm
                % eq. 1 from 10.1016/j.neuron.2010.04.009
                theta = [0.1,0.1,0.1];
                [theta,err,gEst] = fminsearch(@(theta) fitDN_noLatInhibition(theta,g_lin,g_true),theta,opts);
                fN = xp.^theta(2);
                fD = fN + theta(3)^(2*theta(2));
                yp = real(theta(1)^2 * (fN./fD));
        end
        if gEst==0
            error('nonlinear fit could not be achieved.');
        end
        g_nomean = g_true-mean(g_true);
        r2 = 1-(err/sum(((g_nomean).^2)));
        hold on
        plot(xp,yp,'linestyle','-','color','r','marker','none','linewidth',2);
        h = refline(1,0);
        h.Color = 'k';
        h.LineStyle = ':';
        hold off
        legend({'data','fit'});
        legend boxoff
        title(sprintf('%s | glm%d | r2=%0.2f',regname{roi},glm,r2));
        
        varargout = {r2,theta,g_lin,g_true};
    
    case 'pcm_simIndivid'
        % simulate data under linear model. Fit linear and nonlinear
        % scaling both CV and not CV at individual level.
        % Fit of nonlinear CV should fall.
        rmvMean    = 0;
        posJustify = 0;
        % get subject data from a region:
        [~,partVec,condVec,G_hat_all] = pp1_imana('PCM_getData','sn',2:11,'roi',2,'glm',4,...
            'rmvMean',rmvMean,'posJustify',posJustify);
        G_hat = mean(G_hat_all,3); % group-avg second moment
        G_sf        = G_hat(1:5,1:5); % group avg. single-finger G
        U           = mvnrnd_exact(G_sf,100); % single-finger patterns from group-avg single finger G
        Asf         = pcm_diagonalize(G_sf); 
        chords      = pp1_imana('chords');
        scaleParams = log([0.8 0.6 0.4 0.2])';
        powerParam  = 0.5;
        % linearly scaled G:
        G_linear = chords*(Asf*Asf')*chords';
        % nonlinear scaled G (scaling at pattern level):
        Uchord = (chords * U);
        signU = sign(Uchord);
        absU = abs(Uchord);
        OM = (absU.^powerParam).*signU; % do nonlinearity AFTER summation of single finger patterns
        G_nl_power = OM*OM'./ size(U,2); % normalize G to # voxels
        
        
        % define models:
        M{1}.type       = 'component';
        M{1}.numGparams = 1;
        M{1}.Gc         = eye(31);
        M{1}.name       = 'null';
        
        M{2}.type       = 'component';
        M{2}.numGparams = 1;
        M{2}.Gc         = G_linear;
        M{2}.name       = 'linear';
        
        M{3}            = pp1_imana('pcm_nonlinearScale');
        M{3}.theta0     = scaleParams;
        M{3}.Ac         = Asf;
        M{3}.name       = 'nonlinear scalar';
        M{3}.Gc         = M{3}.modelpred(M{3}.theta0,M{3});
        M{3}.fitAlgorithm = 'minimize';
        
        M{4}.type       = 'component';
        M{4}.numGparams = 1;
        M{4}.Gc         = G_nl_power;
        M{4}.name       = 'nonlinear power';
        M{4}.fitAlgorithm = 'minimize';
        
        M{end+1} = pp1_imana('pcm_freedirect'); 
        M{end}.name = 'noiseceiling_cv';
        M{end}.fitAlgorithm = 'minimize';
        M{end+1} = M{end};
        M{end}.name = 'noiseceiling_overall';
        
%         M{5}.type       = 'freechol';
%         M{5}.numCond    = 31;
%         M{5}.name       = 'noiseceiling';
%         M{5}            = pcm_prepFreeModel(M{5});

        trueModel = 2; % 2-linear, 3-nonlinearscalar, 4-nonlinear power
        
        % simulate data under one of these models:
        signal = 1;
        numVox = 100;
        numSim = 10;
        D = [];
        switch trueModel
            case 2 % linear model
                theta0 = [-1];
            case 3 % nonlinear scalar
                theta0 = scaleParams;
            case 4 % nonlinear power
                theta0 = [-1];
        end
        Y = pcm_makeDataset(M{trueModel},theta0,'numVox',numVox,'signal',signal,'noise',1,...
            'numSim',numSim,'samesignal',false,'exact',false, ...
            'design',condVec{1});
        % remove run means from each run?
        if rmvMean
            C0 = indicatorMatrix('identity',partVec{1}); % run-mean centring matrix
            Y = cellfun(@(X) X-C0*pinv(C0)*X,Y,'uni',0);
        end
        % positively justify patterns? (set mininum value to zero):
        if posJustify
            Y = cellfun(@(X) X- min(min(X)),Y,'uni',0);
        end
   
        % fit all models and assess likelihoods
        [Tall,th1] = pcm_fitModelIndivid(Y,M,partVec{1},condVec{1},'runEffect','random');
        [Tcross,D] = pcm_fitModelIndividCrossval(Y,M,partVec{1},condVec{1},'runEffect','random','crossvalScheme','leaveTwoOut');
        
        % plot
        subplot(4,1,1);
        pcm_plotModelLikelihood(Tall,M,'normalize',1,'upperceil',Tall.likelihood(:,end)); ylabel('likelihood');
        subplot(4,1,2);
        pcm_plotModelLikelihood(Tcross,M,'normalize',1,'upperceil',Tcross.likelihood(:,end)); ylabel('crossval likelihood');
        subplot(4,1,3);
        barplot([],Tcross.R2); ylabel('crossval R2');
        subplot(4,1,4);
        barplot([],Tcross.R); ylabel('crossval R');
        
        % ttests:
        ttest(Tcross.likelihood(:,trueModel),Tcross.likelihood(:,3),2,'paired')
        ttest(Tcross.R2(:,trueModel),Tcross.R2(:,3),2,'paired')
        ttest(Tcross.R(:,trueModel),Tcross.R(:,3),2,'paired')
        
        keyboard
    case 'pcm_toySubj'
        % make lineplots of avg. betas per num digits.
        % plots real avg-betas, and model-predicted betas (at the pattern
        % level)
        roi = 2;
        vararginoptions(varargin,{'roi'});
        sn = 2:11;
        glm = 4;
        rmvMean = 0;
        posJustify = 0;
        [Y,partVec,condVec,G_hat] = pp1_imana('PCM_getData','sn',sn,'roi',roi,'glm',glm,...
            'rmvMean',rmvMean,'posJustify',posJustify);
        
        chords = pp1_imana('chords');
        d.numDigits = sum(chords,2);
        d.chordNum = [1:31]';
        d.patternLvl = ones(31,1);
        v=ones(31,1);
        D=[];
        G=[];
        %chords = chords./sum(chords,2);
        % do modelling at pattern level
        for ii=1:numel(sn)
            d.sn = v.*sn(ii);
            g.sn = sn(ii);
            % (1). calculate avg. beta per num digits, using true betas
            C0 = pcm_indicatorMatrix('identity',condVec{ii});
            U  = pinv(C0)*Y{ii};
            d.act = mean(U,2);
            d.model = v;
            D = addstruct(D,d);
            g.ipm = rsa_vectorizeIPM(U*U'./size(U,2)); % (co)var normed by #vox
            g.model = 1;
            g.patternLvl = 1;
            G=addstruct(G,g);
            % (2). calc for linear model prediction (do this per run)
            for jj=unique(partVec{ii})'
                
            end
            U_linear = chords*U(1:5,:);
            d.act = mean(U_linear,2);
            d.ipm = rsa_vectorizeIPM(U*U'./size(U,2));
            d.model = v.*2;
            D = addstruct(D,d);
            g.ipm = rsa_vectorizeIPM(U*U'./size(U,2));
            g.model = 2;
            g.patternLvl = 1;
            G=addstruct(G,g);
            % (3). calc for nonlinear model prediction
            signU = sign(U_linear); % apply scaling to linear model predictions
            absU  = abs(U_linear);
            U_power = (absU.^0.05).*signU;
            d.act = mean(U_power,2);
            d.ipm = rsa_vectorizeIPM(U*U'./size(U,2));
            d.model = v.*3;
            D = addstruct(D,d);
            g.ipm = rsa_vectorizeIPM(U*U'./size(U,2));
            g.model = 3;
            g.patternLvl = 1;
            G=addstruct(G,g);
        end
        % do modelling on G:
        d.patternLvl = zeros(31,1);
        for ii=1:numel(sn)
            % (1) G from betas
            G_subj  = G_hat(:,:,ii);
            d.sn    = v.*sn(ii);
            d.act   = diag(G_subj);
            d.model = v;
            D = addstruct(D,d);
            g.ipm = rsa_vectorizeIPM(G_subj);
            g.model = 1;
            g.patternLvl = 0;
            G=addstruct(G,g);
            % (2) linear model G
            Asf = pcm_diagonalize(G_subj(1:5,1:5));
            Gsf = Asf*Asf';
            G_linear = chords*Gsf*chords';
            d.act = diag(G_linear);
            d.model = v.*2;
            D = addstruct(D,d);
            g.ipm = rsa_vectorizeIPM(G_linear);
            g.model = 2;
            g.patternLvl = 0;
            G=addstruct(G,g);
            % (3) power model G
            signG = sign(G_linear); % apply scaling to linear model predictions
            absG  = abs(G_linear);
            G_power = (absG.^0.05).*signG;
            d.act = diag(G_power);
            d.model = v.*3;
            D = addstruct(D,d);
            g.ipm = rsa_vectorizeIPM(G_power);
            g.model = 3;
            g.patternLvl = 0;
            G=addstruct(G,g);
        end
        
        
        % avg. across chords, grouped by num digits
        D = tapply(D,{'sn','patternLvl','model','numDigits'},{'act','mean'});
        D.scaling = ones(size(D.sn));
        % calculate subject-specific model scaling factors
        for ll = 1:-1:0
            for mm=2:3
                C0 = pcm_indicatorMatrix('identity',D.sn(D.model==mm & D.patternLvl==ll)); % subject avg. scaling factors
                D.scaling(D.model==mm & D.patternLvl==ll) = C0*pinv(C0) * [D.act(D.model==1 & D.patternLvl==ll)./D.act(D.model==mm & D.patternLvl==ll)];
            end
        end
        % plot model G from pattern models:
        modelName = {'avg','linear','power'};
        for mm=1:3
            modelG = mean(G.ipm(G.model==mm & G.patternLvl==1,:),1);
            subplot(2,4,mm);
            imagesc(rsa_squareIPM(modelG));
            title(modelName{mm});
        end
        subplot(2,4,4);
        title('pattern models');
        plt.line(D.numDigits,D.act.*D.scaling,'split',D.model,'style',style.custom({[1 0 0],[0 1 0],[0 0 1]}),'subset',D.patternLvl==1);
        ylabel('mean beta');
        xlabel('# digits');
        % plot model G from second moment models:
        G_hat = rsa_vectorizeIPM(mean(G_hat,3));
        for mm=1:3
            modelG = mean(G.ipm(G.model==mm & G.patternLvl==0,:),1);
            scaling = G_hat/modelG;
            subplot(2,4,mm+4);
            imagesc(rsa_squareIPM(modelG.*scaling));
            title(modelName{mm});
%             if mm>1 % apply same scaling across models
%                 caxis(cs)
%             else
%                 cs=caxis;
%             end
        end
        subplot(2,4,8);
        title('G models');
        plt.line(D.numDigits,D.act.*D.scaling,'split',D.model,'style',style.custom({[1 0 0],[0 1 0],[0 0 1]}),'subset',D.patternLvl==0);
        ylabel('mean variance');
        xlabel('# digits');
        
        varargout = {D,G};  
        
        
        
        
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
        sn  = pp1_imana('getSubjs');
        glm = 4;
        vararginoptions(varargin,{'sn','glm'});
        hemisphere = [1:2];   % left hemi

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
                
                [M,vox2Node]  = caret_vol2surf_own(C1.data,C2.data,images,'topo',topo,'exclude_thres',0.75,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,sprintf('%s_glm%d_hemi%d_finger.metric',subj_name{s},glm,h)),M);
                save(fullfile(caretSDir,sprintf('%s_glm%d_hemi%d_finger_vox2Node.mat',subj_name{s},glm,h)),'vox2Node');
            end
        end;   
    case 'surf_fingerpatterns'             % Make finger pattern jpegs
        %close all;
        sn  = 1;
        glm = 3;
        numDigits = [1,2];
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
    case 'surf_fingerpatterns2Digits'             
        %close all;
        sn  = 1;
        glm = 4;
        vararginoptions(varargin,{'sn','glm'});
        
        h = 1; % left hemi
        groupDir = [gpCaretDir filesep hemName{h} ];
        cd(groupDir);

        switch(h)
            case 1 % left hemi
                coord  = 'lh.FLAT.coord';
                topo   = 'lh.CUT.topo';
                xlims=[-4 18]; % may need to adjust locations for pics
                ylims=[-9 20];
        end
        
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
        close gcf;
        
        % plot patterns for single finger stimulation
        chords        = pp1_imana('chords'); % stimulation chords
        digitsInChord = sum(chords,2);       % how many digits stimulated per chord
        figure('Color',[1 1 1]); % make figure 
        alreadyPlotted = [];
        jj=[1:5];
        titles = {'thumb','index','middle','fourth','little'};
        for dd=1:5 % for each digit:
            % (1) plot this digit's pattern
            subplot(5,5,dd);
            M = caret_plotflatmap('M',M,'col',dd,'data',data,...
                        'border',B.Border,'bordersize',10,'topo',topo,'coord',coord,'bordercolor',{'k.'});
            title(titles{dd});
            % (2) find which two-finger patterns to plot (check they have not already been plotted!)
            toPlot = find(chords(:,dd)==1 & digitsInChord==2);
            toPlot = toPlot(~ismember(toPlot,alreadyPlotted)); % check the patterns have not already been plotted
            alreadyPlotted = [alreadyPlotted; toPlot]; % keep track of what we have plotted
            % (3) plot the appropriate patterns
            for ii=toPlot'
                jj(end+1) = ii+sum(1:dd); % record subplot #s
                subplot(5,5,jj(end));
                M = caret_plotflatmap('M',M,'col',ii,'data',data,...
                        'border',B.Border,'bordersize',10,'topo',topo,'coord',coord,'bordercolor',{'k.'});
            end
        end
        colormap('parula');
        
        mm = 3; % force colour scaling on patterns
        % loop through both figures and all subplots to: force colour
        % scaling and label conditions
        for ii = jj
            subplot(5,5,ii);
            caxis([-mm/2 mm]);   % scale color across plots
            set(gca,'XTick',[]); % remove X and Y axis ticks
            set(gca,'YTick',[]);
            %axis equal;
            box on
            ax = get(gca);
            ax.XAxis.LineWidth = 1;
            ax.YAxis.LineWidth = 1;
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
    
    
        
        
    case '0' % ------------------------------------------------------------
    case '0' % DEPRECIATED / OLD / UNUSED cases
    case '0' % * note that these cases likely won't run without edits
    case '0' % ------------------------------------------------------------
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
    case 'SURF_xhemireg'     % Cross-Register surfaces left / right hem  
        vararginoptions(varargin,{'sn'});
        freesurfer_registerXhem({subj_name{sn}},freesurferDir,'hemisphere',[1 2]); % For debug... [1 2] orig
    case 'SURF_map_ico'      % Align subj surface to atlas surface 
        vararginoptions(varargin,{'sn'});
        freesurfer_mapicosahedron_xhem(subj_name{sn},freesurferDir,'smoothing',1,'hemisphere',[1:2]);
    case 'SURF_make_caret'   % convert freesurfer to caret  
        vararginoptions(varargin,{'sn'});
        caret_importfreesurfer(['x' subj_name{sn}],freesurferDir,caretDir);
    case 'ROI_makeBApaint_UNUSED'   % Make paint file for brodmann area ROIs cut to hand area (saves as ROI_3.paint)
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
    end        
    fprintf('Done.\n');
    case 'ROI_define_OLD'                                                       % Define rois: BA rois, M1/S1 cut to hand area, BA rois cut to hand area
        % Define the ROIs of the group fsaverage atlas for each subject's
        % surface reconstruction. 
        % Output saved for each subject ('s#_regions.mat').
        % The save variable for each subject is a cell array of size
        % {1,#rois}. 
        I = pp1_imana('LIST_subjs');
        sn = pp1_imana('getSubjs');
        vararginoptions(varargin,{'sn'});

        linedef = [5,0,1]; % take 5 steps along node between white (0) and pial (1) surfaces
        
        thalamicMarker = [8133,8115,8109; 8233,8215,8209]; % values of voxels assigned to each of the 3 thalamic regions (vpl,mgn,lgn) split by hemisphere
        corticalMarker = [1,2,3,4,5,6,1,2,7,8]; % values of surface nodes assigned to the cortical rois
        corticalFile   = {'D1','D1','D1','D1','D1','D1','D2','D2','D2','D2'}; % different files for coritcal rois
        
        for s = sn % for each subject
            R = {};
            j = 1; % overall region ticker
            subjName = I.origSN{s};
            fprintf('%s...',subjName);
            caretSubjDir   = fullfile(caretDir,[atlasA subjName]);
            for h = 1:2 % per hemisphere
                t = 1; % thalamic ticker
                D1 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI_pp1.paint'])); % ba3A, ba3B, ba1, ba2, rM1, cM1
                D2 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI.paint']));     % premade ROIs from fsaverage_sym: M1 and S1
                for r = 1:numregions % make regions
                    % make R region structure for participant
                    R{j}.name  = [subjName '_' regname{r} '_' hem{h}];
                    % Mapping of rois is different if cortical or thalamic:
                    if cortical(r)==1 && thalamic(r)==0 % cortical rois
                        eval(['D=' corticalFile{r} ';']);
                        R{j}.type     = 'surf_nodes';
                        R{j}.location = find(D.data(:,1)==corticalMarker(r));
                        R{j}.white    = fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                        R{j}.pial     = fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                        R{j}.topo     = fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                        R{j}.linedef  = linedef;
                        R{j}.flat     = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                        R{j}.image    = fullfile(glmDir{1},subjName,'mask.nii,1'); % load glm mask from first glm
                        if strcmp(corticalFile{r},'D1')
                            R{j}.origFile = 'ROI_pp1.paint';
                        elseif strcmp(corticalFile{r},'D2')
                            R{j}.origFile = 'ROI.paint';
                        end
                    elseif thalamic(r)==1 && cortical(r)==0 % thalamic rois
                        % make temporary functional mask of thalamic nuclei
                        R{j}.type     = 'roi_image';
                        R{j}.value    = thalamicMarker(h,t);
                        R{j}.file     = fullfile(anatomicalDir,subjName,sprintf('%s_thalamicNuclei.nii',subjName));
                        t = t+1;
                    end
                    % add indexing to structure R
                    R{j}.regNum   = j;
                    R{j}.regType  = r;
                    R{j}.cortical = cortical(r);
                    R{j}.thalamic = thalamic(r);
                    R{j}.hemi     = h;
                    j = j+1;
                end    
            end
            % threshold across CS, both hemispheres
            exculdePairs = [1,5; 1,6; 2,5; 2,6; 3,5; 3,6; 7,8];
            excludePairs = [exculdePairs; exculdePairs+numregions]; % do same exclusion for rois in both hemispsheres
            R = region_calcregions(R,'exclude',excludePairs,'exclude_thres',0.75);
            cd(regDir);
            save(['regions_' subjName '.mat'],'R');
            fprintf('..done\n');
            clear R
        end
    case 'ROI_fitHRF_OLD'                                                       % (optional) :  Fit the hrf to the data of a region.
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
        roi = 7;
        eCriteria = 0.964;
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
                fprintf('Epost/Epre: %1.5f\n',e);
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
    case 'ROI_compareReliability_OLD'
        % compares within-session pattern reliabilities to across-session
        % reliabilities. 
        % Reliabilities are odd-even splits.
        % corr(A_s1,B_s1) = r_s1
        % corr(A_s2,B_s2) = r_s2
        % corr(A_s1,A_s2) = r_A
        % corr(B_s1,B_s2) = r_B
        % ssqrt(r_s1 * r_s2) vs. ssqrt(r_A * r_B)
        sn  = 1:7;
        glm = 3;
        roi = [1:6];
        vararginoptions(varargin,{'roi','glm','sn'});
        
        conds = 1:31;
        mean_sub = 1;
        betaType = 'raw';
        
        R = [];
        for region=roi
            for s = sn
                r_s1 = pp1_imana('ROI_patternReliability','sn',s,'roi',region,'glm',glm,'split','sess1','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
                r_s2 = pp1_imana('ROI_patternReliability','sn',s,'roi',region,'glm',glm,'split','sess2','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
                r_A  = pp1_imana('ROI_patternReliability','sn',s,'roi',region,'glm',glm,'split','sess_cv_odd','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
                r_B  = pp1_imana('ROI_patternReliability','sn',s,'roi',region,'glm',glm,'split','sess_cv_even','mean_subtract',mean_sub,'betaType',betaType,'conds',conds);
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
                r.roi   = [region;region];
                r.glm   = [glm;glm];
                R = addstruct(R,r);
            end
        end
        % plt.dot(R.sn,R.corr,'split',R.type);
        %save(fullfile(regDir,sprintf('patternReliability_roi%d_glm%d.mat',roi,glm)),'-struct','R');
        varargout = {R};
    case 'ROI_patternReliability_UNUSED'                                              
        % Splits data for each session into two partitions (even and odd runs).
        % Calculates correlation coefficients between each condition pair 
        % between all partitions.
        % Default setup includes subtraction of each run's mean
        % activity pattern (across conditions).
        glm           = 3;
        roi           = 7; % default roi
        sn            = pp1_imana('getSubjs');
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
    case 'ROI_getSingleFingerTuning_multiFinger_UNUSED'
        % wrapper to perform SFT analysis
        %here, we test conditions with multiple fingers pressed.
        sn  = pp1_imana('getSubjs');
        glm = 3;
        roi = 1:4;
        vararginoptions(varargin,{'glm','roi','sn'});
        
        numSim = 10; % # simulated datasets per model per participant (mvnrnd and sparse tuning models)
        chords = pp1_imana('chords');
        numDigits = sum(chords,2);
        numCond = 5;
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        T = getrow(T,ismember(T.sn,sn) & ismember(T.roi,roi));
        
        D = []; % output structure
        v = ones(3,1); % helper vector
        for numD = 2:4 % can't do all 5 fingers chord
            fprintf('%d-finger configurations:\n',numD);
            conds = find(numDigits==numD); % conditions to analyze
            numChords = numel(unique(conds));
            for ii = 1:size(T.sn,1)
                % some simulation params needed:
                numVox  = ceil(size(T.raw_beta{ii},2)/numCond)*numCond; % round up so equal # of voxels per condition (for sparse patterns)
                numRun  = numel(unique(T.run{ii}));
                
                % for each voxel in subj-roi, get avg. betas for single fingers
                b           = [];
                b.beta      = T.raw_beta{ii}; %T.betaUW{ii};
                b.tt        = T.tt{ii};
                b.run       = T.run{ii};
                b           = getrow(b,ismember(b.tt,conds));
                %C0          = indicatorMatrix('identity',b.run);
                %b.beta      = b.beta - C0*pinv(C0)*b.beta; % remove run means to avoid potential for scaling across runs to bias averaging

                % don't need to mean centre patterns because the metric is
                % insensitive to where the mean is (and is calculated per voxel)
                
                % for multi-finger analysis, we need to convert from chord
                % conditions to responses to single fingers. 
                % A simple option is to use OLS to estimate single-finger 
                % patterns per run (split per # fingers in chords)
                r2 = nan;
                if numD>1
                    % ols
                    X = kron(eye(numRun),chords(conds,:));
                    X = [X,kron(eye(numRun),ones(numChords,1))]; % add run intercepts
                    B = pinv(X)*b.beta;
                    Y_hat = X*B;
                    res   = b.beta - Y_hat;
                    r2    = 1- (sum(res.^2,1)./sum(b.beta.^2,1)); % voxel-specific R2s
                    % remake b structure with single finger patterns
                    b = [];
                    b.beta = B(1:numRun*numCond,:);
                    b.tt   = kron(ones(numRun,1),[1:5]');
                    b.run  = kron([1:numRun]',ones(5,1));
                end

                % get estimate of signal, estimate of noise, and
                % single-finger G
                [evar,svar] = pp1_imana('estErrVar',b.beta,b.tt,b.run);
                b       = tapply(b,{'tt'},{'beta','mean'});
                Gtt     = cov(b.beta');

                % 1. calc tuning of actual voxel data
                sftBeta = pp1_imana('estSingleFingerTuning',b.beta);
                sftBeta = mean(sftBeta);

                % 2. calc expected tuning of voxels with ~N(0,G)
                [sftEV,sftDistEV] = pp1_imana('SFT:calcExpectedValueG',evar,svar,Gtt,numVox,numRun,numSim); % expected value of the null

                % 3. calc expected tuning of voxels with ~sparse tuning, but where signal-to-noise equals that voxel data
                [sftSp,sftDistSp] = pp1_imana('SFT:calcExpectedValueSparse',evar,svar,numCond,numVox,numRun,numSim);

                % 4. do prob test for each subject on their gauss sft distribution
                pEV = sum(sftBeta<=sftDistEV)/length(sftDistEV);
                if isempty(pEV)
                    pEV = realmin;
                end

                % 5. do prob test for each subject on their sparse sft distribution
                pSp = sum(sftBeta<=sftDistSp)/length(sftDistSp);
                if isempty(pSp)
                    pSp = realmin;
                end

                % add to output structure
                d = [];
                d.passive   = v;
                d.r2        = v.*mean(r2);
                d.numDigits = v.*numD;
                d.sn        = v.*T.sn(ii);
                d.roi       = v.*T.roi(ii);
                d.glm       = v.*glm;
                d.sft       = [sftBeta;sftEV;sftSp];
                d.sftProb   = [0;pEV;pSp];
                d.isEV      = [0;1;2];
                D = addstruct(D,d);
                fprintf('s%02d roi%02d done\n',T.sn(ii),T.roi(ii));
            end
        end
        save(fullfile(regDir,sprintf('sft_mf_glm%d',glm)),'-struct','D');
        varargout = {D};
    case 'SFT:estErrVarMeans_DEPRECIATED'
        % WARNING: REMOVE RUN MEANS PRIOR TO ESTIMATION! 
        error('this case is depreciated')
        
        % empirically estimates error variance in activity patterns across
        % runs. 
        
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
    case 'SFT:getSingleFingerG_UNUSED'
        glm = 3;
        roi = 1:8;
        sn = pp1_imana('getSubjs');
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
    case 'calcFstat_OLD'
        % calculates F-statistic per voxel to determine if voxel is
        % significantly modulated by finger(s)
        Y        = varargin{1}; % N x P matrix of data. (N=numCond*numRun, P=numVox)
        condVec  = varargin{2}; % N x 1 vector of condition assignments
        runVec   = varargin{3}; % N x 1 vector of run assignments
        % housekeeping
        numVox   = size(Y,2);
        numCond  = length(unique(condVec));
        numRun   = length(unique(runVec));
        % remove run means 
        C0  = indicatorMatrix('identity',runVec);
        Y   = Y - C0*pinv(C0)*Y; 
        % calc F-stat per voxel
        A = zeros(numCond,numVox,numRun);
        ApredCV = A;
        for ii = 1:numRun
            A(:,:,ii) = Y(runVec==ii,:);
        end
        runs = 1:numRun;
        for ii = 1:numRun
            testRuns = runs~=ii;
            ApredCV(:,:,ii) = mean(A(:,:,testRuns),3);
        end

        % non-crossval f-stat
        Apred = repmat(mean(A,3),1,1,numRun); % predicted voxel tunings
        TSS   = sum(sum((A).^2,3),1);         % total SS (null model- intercept)
        RSS   = sum(sum((A-Apred).^2,3),1);   % finger model (5 params)
        dfN   = numCond - 1;                  % numerator DF is num params (5 fingers) - 1 (intercept)    
        dfD   = numCond*numRun - numCond -numel(unique(runVec));
        Fstat = ((TSS-RSS)./dfN) ./ (RSS./dfD);
        % crossval f-stat
        RSScv = sum(sum((A-ApredCV).^2,3),1); % unrestricted SSR
        FstatCV = ((TSS-RSScv)./dfN) ./ (RSS./dfD);
        Fcrit   = finv(0.95,dfN,dfD); % 95% cutoff for F-stat
        varargout = {Fstat,FstatCV,Fcrit};
        %keyboard
        
        %         % calculate Non-CrossValidated F
%         A = zeros(numCond,numVox,numRun); % zero-pad
%         ApredCV = A;
%         for ii = 1:numRun
%             A(:,:,ii) = Y(runVec==ii,:);
%         end
%         muK = mean(A,3);
%         res = bsxfun(@minus,A,muK);
%         SSR = sum(sum(res.^2,1),3); % common residual covariance
%         SSB = sum((muK.^2).*numRun,1);
%         F   = (SSB./df1) ./ (SSR./df2);
%         
%         % calculate CrossValidated F
%         df2 = numCond*numRun - numCond - numRun - 1;
%         for ii = runs
%             testRuns = runs~=ii;
%             ApredCV(:,:,ii) = mean(A(:,:,testRuns),3);
%         end
%         muK = mean(ApredCV,3);
%         rescv = A-ApredCV;
%         SSRcv = sum(sum(res.^2,1),3); % common residual covariance
%         SSBcv = sum((ApredCV.^2).*numRun,1);
%         F   = (SSB./df1) ./ (SSR./df2);
%         
% 
% 
%         % non-crossval f-stat
%         Apred = repmat(mean(A,3),1,1,numRun); % predicted voxel tunings
%         TSS   = sum(sum((A).^2,3),1);         % total SS (null model- intercept)
%         RSS   = sum(sum((A-Apred).^2,3),1);   % finger model (5 params)
%         dfN   = numCond - 1;                  % numerator DF is num params (5 fingers) - 1 (intercept)    
%         dfD   = numCond*numRun - numCond -numel(unique(runVec));
%         Fstat = ((TSS-RSS)./dfN) ./ (RSS./dfD);
%         % crossval f-stat
%         RSScv = sum(sum((A-ApredCV).^2,3),1); % unrestricted SSR
%         FstatCV = ((TSS-RSScv)./dfN) ./ (RSS./dfD);
%         Fcrit   = finv(0.95,dfN,dfD); % 95% cutoff for F-stat
    case 'ROI_snr' % Calculate intercept and slope for each activation value
        
        sn  = pp1_imana('getSubjs');
        glm = 3;
        roi = [1:4];
        vararginoptions(varargin,{'roi','glm','sn'});
        
        % Load betas
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)));
        T = getrow(T,ismember(T.roi,roi) & ismember(T.sn,sn));
        % do calculations
        R = [];
        for ii = 1:length(T.sn)
            Y   = T.raw_beta{ii}(T.tt{ii}<6,:); % take only single fingers
            c   = T.tt{ii}(T.tt{ii}<6);
            snr = pp1_imana('ROI_calcSnr',Y,c);
            % add to output structure
            v = ones(length(snr),1);
            r.snr   = snr';
            r.vox   = [1:length(snr)]';
            r.sn    = v.*T.sn(ii);
            r.roi   = v.*T.roi(ii);
            r.glm   = v.*glm;
            R = addstruct(R,r);
        end
        varargout = {R};
    case 'ROI_calcSnr'
        % signal^2 ./ variance per voxel
        Y = varargin{1};    % patterns [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        
        P   = size(Y,2);     % # of voxels
        cc  = unique(c);     % condition labels
        nC  = length(cc);    % # of conditions
        muK = zeros([nC P]); % zero-pad for condition means
        Sw  = muK;           % zero-pad for w/in class (& voxel) variabilities
        
        for i=1:nC
            j        = find(c==cc(i));                 % select rows of Y in this class
            n        = length(j);                      % number of trials per class
            muK(i,:) = sum(Y(j,:),1)./n;               % condition means
            res      = bsxfun(@minus,Y(j,:),muK(i,:)); % residuals of trials in this class
            Sw(i,:)  = sum(res.^2,1)./n;               % variance of this class
        end 
        
        snr = (muK.^2)./Sw; % snr per voxel per condition
        snr = mean(snr,1);  % avg. snr per voxel across conditions 
        
        varargout ={snr};
    case 'ROI_calcSnr_cv'
        % signal^2 ./ variance per voxel, crossvalidated fashion
        Y = varargin{1};    % patterns [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        r = varargin{3};    % partitions
        
        %Y = Y'; % reorient to [PxN]
        nV  = size(Y,2);     % # of voxels
        cc  = unique(c);     % condition labels
        nC  = length(cc);    % # of conditions
        rr  = unique(r);
        nR  = length(rr);
        muK = zeros([nC nV nR]); % zero-pad for condition means
        Sw  = muK;           % zero-pad for w/in class (& voxel) variabilities
        
        for q=1:nR
            for i=1:nC
                test       = find(c==cc(i) & r==rr(q));      % select rows of Y in this class from this partition
                train      = find(c==cc(i) & r~=rr(q));      % select rows of Y in this class from this partition
                n          = length(train);                  % number of trials per class
                muK(i,:,q) = sum(Y(train,:),1);        % condition mean of training class data
                res        = Y(test,:)-muK(i,:,q);      % residuals of test trial in this class
                Sw(i,:,q)  = sum(res.^2,1)./length(test);  % variance of this class
            end 
        end 
       % Sw  = sum(Sw,3)./(nR);
       % muK = sum(muK,3)./(nR);
        
        snr = ((muK.^2)-Sw)./Sw; % snr per voxel per condition
        snr = mean(mean(snr,3),1);  % avg. snr per voxel across conditions 
        
        varargout ={snr};
    case 'ROI_calcSnr_r2cv'
        % signal^2 ./ variance per voxel
        Y = varargin{1};    % patterns [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        r = varargin{3};    % partitions [regressors x 1]
        
        nV  = size(Y,2);        % # of voxels
        cc  = unique(c);        % condition labels
        nC  = length(cc);       % # of conditions
        rr  = unique(r);        % partition labels
        nR  = length(rr);       % # of partitions
        X   = indicatorMatrix('identity_p',c);
        muKa = zeros([nC nV nR]); % zero-pad for condition means
        muKb = muKa;
        
        % remove run means:
        C0 = indicatorMatrix('identity',r);
        Y0 = Y - C0*pinv(C0)*Y;
        
        % Estimate condition means within each run
        r2cv_cond = zeros(nC,nR);
        r2cv_vox  = zeros(nV,nR);
        for i=1:nC
            % test data
            Xa = X(r==rr(i),:);
            Ya = Y0(r==rr(i),:);
            muKa(:,:,i) = pinv(Xa)*Ya;
            % training data
            Xb = X(r~=rr(i),:);%./(nR-1);
            Yb = Y0(r~=rr(i),:);
            muKb(:,:,i) = pinv(Xb)*Yb;
        end 

        res = muKa-muKb;
        SSR = sum(sum(res.^2,3),1);
        SST = sum(sum(muKa.^2,3),1);
        R2  = 1-SSR./SST;
        
        varargout ={R2};
    case 'ROI_calcSnr_rcv'
        % calculates snr as crossvalidated pearson's correlation (per voxel)
        Y = varargin{1};    % patterns   [regressors x voxels]
        c = varargin{2};    % conditions [regressors x 1]
        r = varargin{3};    % partitions [regressors x 1]
        
        nV  = size(Y,2);        % # of voxels
        cc  = unique(c);        % condition labels
        nC  = length(cc);       % # of conditions
        rr  = unique(r);        % partition labels
        nR  = length(rr);       % # of partitions
        X   = indicatorMatrix('identity_p',c);
        muKa = zeros([nC nV nR]); % zero-pad for condition means
        muKb = muKa;
        
        % Estimate condition means within each run
        rcv_cond = zeros(nC,nR);
        rcv_vox  = zeros(nV,nR);
        for i=1:nR
            % test data
            Xa = X(r==rr(i),:);
            Ya = Y(r==rr(i),:);
            muKa(:,:,i) = pinv(Xa)*Ya;
            % training data
            Xb = X(r~=rr(i),:);%./(nR-1);
            Yb = Y(r~=rr(i),:);
            muKb(:,:,i) = pinv(Xb)*Yb;
            % correlate across test and training
            rcv_cond(:,i) = diag(corr(muKa(:,:,i)',muKb(:,:,i)'));
            rcv_vox(:,i)  = diag(corr(muKa(:,:,i),muKb(:,:,i)));
            
%             if (removeMean) % per partition
%                 A(:,:,i)=bsxfun(@minus,A(:,:,i),sum(A(:,:,i),1)/numCond);
%             end
        end
        
        % avg. correlations across testing folds
        rcv_cond = mean(rcv_cond,2); % not invariant to scaling across voxels
        rcv_vox  = mean(rcv_vox,2);  % invariant to scaling across voxels
        
        varargout ={rcv_vox,rcv_cond};
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
%     case 'depreciated_CM_getGrayBetas'  % resample patterns from all gray matter voxels (irrespective of region)
%         sn  = 1;
%         glm = 2;
%         % load in beta vol data
%         HCPDir=fullfile(atlasDir,'HCP/');
%         cd(HCPDir);
%         HCPFiles=dir('*HCP_*');
%         
%         % get V and volIndx
%         load(fullfile(studyDir{2},'encoding','glm4','cereb_avrgDataStruct.mat'));
%         
%         % get HCP contrasts
%         for i=1:length(HCPFiles),
%             VA{i}=spm_vol(fullfile(HCPDir,HCPFiles(i).name));
%         end
%         
%         % now sample the contrasts into the same space
%         [i,j,k]=ind2sub(V.dim,volIndx);
%         [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
%         [i1,j1,k1]=spmj_affine_transform(x,y,z,inv(VA{1}.mat));
%         for i=1:length(VA),
%             map(i,:) = spm_sample_vol(VA{i},i1,j1,k1,0);
%             colNames{i,1}=HCPFiles(i).name(5:end-9);
%         end
%         
%         % normalise data
%         X_C=bsxfun(@minus,map,nanmean(map));
%         
%         C{1}.dim=V.dim;
%         C{1}.mat=V.mat;
%         
%         % make volume
%         Vi=zeros(size(map,1), [V.dim(1)*V.dim(2)*V.dim(3)]);
%         Vi(:,volIndx)=map;
%         
%         % map vol2surf
%         for i=1:size(map,1),
%             data=reshape(Vi(i,:,:,:),[C{1}.dim]);
%             C{i}.dat=data;
%         end
%         S=caret_suit_map2surf(C,'space','SUIT','stats','nanmean','column_names',colNames);  % MK created caret_suit_map2surf to allow for output to be used as input to caret_save
%         
%         caret_save(fullfile(HCPDir,'HCPContrasts.metric'),S);
%         varargout={X_C,colNames};    
%     case 'depreciated_CM_logicalGrayMatterMask'
%         % define t/f gray matter mask for cerebellum
%         sn = 1;
%         vararginoptions(varargin,{'sn'});
%         G                       = spm_vol(fullfile(cerebAnatDir,subj_name{sn},[subj_name{sn} '_anatomical_seg1.nii']));
%         G.data                  = spm_read_vols(G);
%         G.data(G.data>0.25)     = 1; % mask out gray matter voxels
%         G.data(G.data<1)        = 0;
%         G.fname                 = fullfile(cerebAnatDir,subj_name{sn},'gray_mask.nii');
%         G                       = rmfield(G,{'pinfo'});
%         spm_write_vol(G,G.data);
%         fprintf('Done %s cerebellar gray matter mask.\n',subj_name{sn});
%     case 'depreciated_CM_defineFunctionalMask'
%         % make functional mask for volumetric searchlight of cerebellar
%         % data.
%         glm = 2;
%         sn  = 1;
%         vararginoptions(varargin,{'sn','glm'});
%         
%         % load mask files to combine (functional and anatomical gray)
%         F      = spm_vol(fullfile(glmDir{glm},subj_name{sn},'mask.nii'));
%         F.data = spm_read_vols(F);
%         G      = spm_vol(fullfile(cerebAnatDir,subj_name{sn},'gray_mask.nii'));
%         G.data = spm_read_vols(G);
%         Vin(1) = F;
%         Vin(2) = G;
%         Vo       = Vin(1); 
%         Vo.fname = 'func_grayMask.nii'; 
%         % find gray matter voxels of cerebellum for which we have
%         % functional imaging data (make mask img in functional coords)
%         spm_imcalc(Vin,Vo,'i1.*i2');
%         fprintf('Done.\n')
%     case 'CM_defineSearchVol'           % define volumetric searchlight in cerebellum 
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
            J.eoptions.tpm      = {'/Users/sarbuckle/MATLAB/spm12/toolbox/Anatomy/wgrey.nii'};%{'/Users/sarbuckle/Documents/MotorControl/matlab/spm12/tpm/TPM.nii'};
            J.eoptions.affreg   = 'mni';
            J.eoptions.reg      = [0 0.001 0.5 0.05 0.2];
            J.eoptions.fwhm     = 0;
            J.eoptions.samp     = 3;
            J.woptions.bb       = [-78 -112 -70; 100 76 85];%[-78 -112 -70; 78 76 85];
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
    
        
        
    otherwise
        fprintf('%s: no such case.\n',what)

end
cd(cwd);
%% ------------------------- Local funcs ----------------------------------
function dircheck(dir)
% Checks existance of specified directory. Makes it if it does not exist.
% SA 01/2016
if ~exist(dir,'dir');
    %warning('%s didn''t exist, so this directory was created.\n',dir);
    mkdir(dir);
end
end

function err = fitPowerScaling(theta,x,y)
% theta : params
% x : data used to fit y
% y : data we are fitting
signX = sign(x);
absX  = abs(x);
yh = ((absX.^theta(1)).*signX).*theta(2);
err = sum(((y(:)-yh(:)).^2));
end
function err = fitPower(theta,x,y)
% theta : params
% x : data used to fit y
% y : data we are fitting
signX = sign(x);
absX  = abs(x);
yh = ((absX.^theta(1)).*signX);
% ssCov = sum(y(:).*yh(:));
% ssYh  = sum(yh(:).*yh(:));
% ssY   = sum(y(:).*y(:));
% err   = ssCov/sqrt(ssYh.*ssY); %corr is cost
err = sum(((y(:)-yh(:)).^2)); %SSR is cost
end

function err = fitQuad(theta,x,y)
% theta : [1x3]
% x : data used to fit y
% y : data we are fitting
yh  = theta(1) + theta(2)*x + theta(3)*(x.^2);
err = sum(((y-yh).^2));
end
function err = fitCompress(theta,x,y)
% case to fit a compressive nonlinearity
% theta : power of nthroot
% x : data used to fit y
% y : data we are fitting
yh  = real((theta(1) * (x.^theta(2))));
err = sum(((y-yh).^2));
end
function err = fitExponent(theta,x,y)
% case to fit an exponent scaling (for decimals, this makes them smaller,
% so works like sqrt scaling)
% theta : power of nthroot
% x : data used to fit y
% y : data we are fitting
yh  = real(x.^(theta(1)));
% add back in single finger elements:
G_est  = rsa_squareIPM(yh);
G_true = rsa_squareIPM(y);
G_est(1:5,1:5) = G_true(1:5,1:5);
yh = rsa_vectorizeIPM(G_est);
% calc error:
err = sum(((y-yh).^2));

end
function err = fitGainExponent(theta,x,y)
% case to fit a gain & squareroot
% theta : power of nthroot
% x : data used to fit y
% y : data we are fitting
yh  = real(theta(1) * (x.^theta(2)));
% add back in single finger elements:
G_est  = rsa_squareIPM(yh);
G_true = rsa_squareIPM(y);
G_est(1:5,1:5) = G_true(1:5,1:5);
yh = rsa_vectorizeIPM(G_est);
% calc error:
err = sum(((y-yh).^2));

end
function err = fitDN_noLatInhibition(theta,x,y)
% fits population level divisive normalization
% (eq1. from 10.1016/j.neuron.2010.04.009)
% theta : [1x4]
% x : data used to fit y    -> covariance elements of linear model
% y : data we are fitting   -> covariance elements of data

% numerator   = x.^theta(2);
% denominator = numerator + theta(3)^theta(2);
% yh  = real(theta(1) * (numerator./denominator));

% divisive normalization at single voxel level can be treated as f(ZU).
% Therefore, G = f(ZU) * (f(ZU))'
% where f = q1( (ZU)^q2 / ((ZU)^q2 + q3^q2))
% We expand this, and it yeilds: G = q1^2( (ZUU'Z')^q2 / ((ZUU'Z')^q2 + q3^2q2))
% numerator = x.^theta(2);
% denominator = numerator.^2 + 2*numerator*theta(3).^theta(2) + theta(3)^(2*theta(2));
% yh  = real(theta(1)^2 * (numerator./denominator));

% we only fit multi-finger elements, as the single finger elements might
% already be suppressed. So we just add the single finger element back into
% the final G.

numerator = x.^theta(2);
denominator = numerator + theta(3)^(2*theta(2));
yh  = real(theta(1)^2 * (numerator./denominator));
% add back in single finger elements:
G_est  = rsa_squareIPM(yh);
G_true = rsa_squareIPM(y);
G_est(1:5,1:5) = G_true(1:5,1:5);
yh = rsa_vectorizeIPM(G_est);
% calc error:
err = sum(((y-yh).^2));
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