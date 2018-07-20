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

% SArbuckle, Motor Control Group, 2018, UWO
% saarbuckle@gmail.com


% -------------------------------------------------------------------------
cwd = cd; % get current directory when called, and return at end of script
% ------------------------- Directories -----------------------------------
codeDir         ='/Users/sarbuckle/Dropbox (Diedrichsenlab)/Arbuckle_code/projects/project_passivePatterns';
baseDir         ='/Users/sarbuckle/Documents/MotorControl/data/passivePatterns1';   % base directory for analysis
% directories for "raw" data (i.e. data during/from preprocessing)
behavDir        =[baseDir '/data'];          % behavioural data directory
dicomDir        =[baseDir '/data_dicom'];    % imgs hot off the scanner
imagingDirRaw   =[baseDir '/imaging_data_raw']; 
phaseDirRaw     =[baseDir '/phase_data_raw'];
% working data directories (this data has been preprocessed and is ready to go)
fieldmapDir     =[baseDir '/fieldmaps/'];    % this explicit path is not required here, but is helpful to know it is required by realign_unwarp
imagingDir      =[baseDir '/imaging_data'];               
phaseDir        =[baseDir '/phase_data'];
anatomicalDir   =[baseDir '/anatomicals'];                
freesurferDir   =[baseDir '/surfaceFreesurfer'];          
caretDir        =[baseDir '/surfaceCaret'];     
gpCaretDir      =[caretDir '/fsaverage_sym'];
regDir          =[baseDir '/RegionOfInterest/'];   % where most of the important analysis data structures will be saved       
% update glmDir when adding new glms
glmDir          ={[baseDir '/glm1'],[baseDir '/glm2'],[baseDir '/glm3'],[baseDir '/glm4'],[baseDir '/glm5']};              
    
filePrefix = 'pp1';
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
numDummys  = 2;  % dummy images at the start of each run (these are discarded)
numTRs     = {[360,360,360,370,370,370,370]};  % total # of images per run (including dummies)
TR_length  = 1.5;  % seconds
run        = {[1:7]};

% ------------------------- ROI things ------------------------------------
hem        = {'lh','rh'};                                                   % left & right hemi folder names/prefixes
regname    = {'sB3a','sB3b','sB1','sB2','S1','M1','B1','B2','B3','B4','B6'}; % Cortical ROIs, 5 = S1, 6 = M1};                                             % roi names, independent of hemisphere    
regSide    = [ones(size(regname)),...                                       % Hemisphere of the roi
                ones(size(regname)).*2];                                    % [1 = left hemi (contra), 2 = right hemi (ipsi)]
regType    = [1:length(regname),...                                         % roi # (within hemisphere)
                1:length(regname)];
numregions = max(regType);                                                  % total number of regions 
% title of regions ordered according to numerical call id (eg. 2 = Lh M1)
reg_title  = {'LHsB3a','LHsB3b','LHsB1','LHsB2','LHS1','LHM1','LHB1','LHB2','LHB3','LHB4','LHB6'}; % Cortical ROIs, 5 = S1, 6 = M1}; % for example..-> these included hemisphere in the names

% ------------------------- Freesurfer things -----------------------------         
atlasA    = 'x';                                                            % freesurfer filename prefix
atlasname = 'fsaverage_sym';                                                % freesurfer average atlas
hemName   = {'LeftHem','RightHem'};                                         % freesurfer hemisphere folder names    

% ------------------------- Voxel Depth/Layer things ----------------------
% Although not true cortical 'layers', these 'layers' facilitate 
% harvesting of data from voxels at specified depths along the grey matter 
% sheet. 
%
% The layers are:
%       1  :  'all' voxels
%       2  :  'superficial' voxels
%       3  :  'deep' voxels
%
% Voxels with a depth of zero (or negative) have their centroid located on
% (or above) the grey matter surface (constructed with freesurfer tools).
% Voxels with a depth of 1 (or greater) have their centroid located at (or
% in) the grey & pial matter juncture.
% Voxels may have centroids that don't fall within the grey matter sheet
% because they are still included in the grey matter mask, given portions
% of the voxel (and thus portions of their signal) originate from grey
% matter.
layers     = { [-Inf Inf], [-Inf 0.5], [0.5 Inf] };
layer_name = {'all','superficial','deep'};

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
subj_name  = {'p01'};
loc_AC     = {[-107 -164 -173]};
DicomName  = {'2018_07_12_P01.MR.Diedrichsen_PassivePatterns'};
NiiRawName = {'2018_07_12_P01'};
fscanNum   = {[18,21,25,28,31,34,43]};                   
pscanNum   = {[19,22,26,29,32,35,44]};            
anatNum    = {[]};                         
fieldNum   = {[45,46]};                                               

% ------------------------- Analysis Cases --------------------------------
switch(what)
    case '0' % ------------ MISC: some aux. things ------------------------
    case 'MISC_checkTime'                                                  % Check alignment of scanner and recorded time (sanity check): enter sn
        vararginoptions(varargin,{'sn'});

        D=dload(fullfile(behavDir,sprintf('%s_%s.dat',filePrefix,subj_name{sn})));
        figure('Name',sprintf('Timing of Task Onsets vs. TR onsets for Subj %d',sn),'NumberTitle','off')
        
        % plot alignment of TR time and trial onset time
        for b = unique(D.BN)'
            subplot(2,length(run{sn}),b);
            d = getrow(D,D.BN==b);
            %subplot(2,1,1); plot(D.realStartTime/1000,(D.realStartTR-1)*0.7+D.realStartTRTime/1000)
            plot(d.StartTimeMeas,(d.mStartTR-1)*TR_length*1000 + d.mStartTRTime,'LineWidth',1.5,'Color','k')
            title(sprintf('run %d',run{sn}(b)));
            xlabel('trial start time (ms)');
            ylabel('tr start time (ms)');
            grid on
            axis equal
        end
        
        % plot difference of TR time and trial onset time
        subplot(2,length(run{sn}),[length(run{sn})+1 : length(run{sn})*2]); 
        plot(D.StartTimeMeas - ((D.mStartTR-1)*TR_length*1000 + D.mStartTRTime),'LineWidth',1.5,'Color','k')
        ylabel('trial onset - tr time (ms)');
        xlabel('trial number');
        title('Difference of Trial onset time and TR time')
        xlim([0 length(D.BN)]);
        hold on
        % draw line marking each new run
        for r = 2:length(run{sn})
            drawline((r-1)*62,'dir','vert','linestyle',':');
        end
        hold off
        %keyboard
        
        %__________________________________________________________________
    case 'MISC_checkMovement'                                              % Check movement of subject. Requires GLM 3 for specified subject.
        vararginoptions(varargin,{'sn'});
        glm = 1;
        load(fullfile(glmDir{glm},subj_name{sn},'SPM.mat'));
        spm_rwls_resstats(SPM)          
    case 'MISC_getSNR'                                                      % UNDER DEVELOPMENT
		% Harvests signal-to-noise ratio of subjects.
		% SNR calculated by dividing variance from voxel within hand knob (m1/s1) by that from some white matter voxel.
		% It's crude but also helpful.
		
		sn = 1;
		vararginoptions(varargin,{'sn'});
		
		gm_coord = {[31,36,18]};
		wm_coord = {[31,49,14]};
		D = [];
		for r = 1:numel(run)
			fname    = fullfile(imagingDir,subj_name{sn},sprintf('rs%02d_run_%02d.nii',sn,r));
			ts       = spm_vol(fname);
			ts       = spm_read_vols(ts);
            numimgs  = [1:size(ts,4)];
			d.gm_var = var(squeeze(ts([gm_coord{sn},numimgs])));
			d.wm_var = var(squeeze(ts([wm_coord{sn},numimgs])));
			d.SN     = sn;
			D.run    = r;
			D = addstruct(D,d);
		end
		keyboard
	case 'MISC_CorticalWidth'                                               % Harvest cortical width (distance in mm between white and grey surfaces) in specified roi.
        % Harvest cortical widths along specified roi. 
        % Widths are in mm.
        sn   = 1;
        hemi = 1;     % left hemi default (2 is right hemisphere)
        roi  = 10;    % set to 0 for no specified roi
        vararginoptions(varargin,{'sn','numbins','hemi','roi'})
        
        if roi>12
            hemi=2;
        end
        D = [];
        for s=sn;
            % load pial and white-gray surfaces
            a = caret_load(fullfile(caretDir,['x',subj_name{s}],hemName{hemi},sprintf('%s.PIAL.coord',hem{hemi})));
            b = caret_load(fullfile(caretDir,['x',subj_name{s}],hemName{hemi},sprintf('%s.WHITE.coord',hem{hemi})));
            % Both structures contain .data subfield, which has [x,y,z]
            % coords in mm of each vertex location.
            % Calc distance between surfaces:  d = sqrt((x2 - x1)^2)
            d = sqrt(sum((a.data-b.data).^2,2));
            % load region data- want R.location arrays  
            load(fullfile(regDir,[subj_name{s},'_regions.mat'])); 
            for r=roi;
                if r==0 % roi = 0 if we take all rois together
                    dindx = [1:length(d)];
                    dat   = d(dindx);
                else
                    dindx = R{r}.location;
                    dat   = d(dindx);
                end
                t.depth = mean(dat);
                t.SN    = s;
                t.reg   = r;
                t.hemi  = hemi;
                D       = addstruct(D,t);
            end
        end
        varargout = {D};
    case 'MISC_SEARCH_calculate_contrast'                                   % Called by 'SEARCH_map_contrast': calculates distances from searchlight results at each node using submitted contrast matrix.
        vararginoptions(varargin,{'sn','glm','C','file'});
        % C = contrast matrix
        % file = filename of searchlight results (nifti extension)
        %   - file can be string array if calling for multiple subjs
        % % (1). Create variables that remain identical across subjs
        K = 20;                     % num conditions
        H = eye(K) - ones(K)/K;     % centering matrix (sets baseline to avg. for G)
        Y = [];                     % output structure
        
        % % Loop through subjects
        for s = sn
            % % (2). Load subject surface searchlight results (1 vol per paired conds)
            if length(sn)>1
                [subjDir,fname,ext] = fileparts(file{s-(s-1)});             % if looping over many subjs
            else
                [subjDir,fname,ext] = fileparts(file);                      % if doing for only one subj
            end
            cd(subjDir);
            vol  = spm_vol([fname ext]);
            vdat = spm_read_vols(vol);                                      % searchlight data
            % % (3). Get distances at each voxel
            [xVox,yVox,nslices,lRDM] = size(vdat);
            V.all = zeros((xVox*yVox*nslices),lRDM);                        % preallocate voxel full RDM field
            for i = 1:lRDM
                V.all(:,i) = reshape(vdat(:,:,:,i),[],1);                   % string out voxels to one dimension (easier indexing) & get distances described at each voxel
            end
            clear vdat
            
            % % (3). Calc G for each voxel, est. new distances
            LDC{s} = zeros((xVox*yVox*nslices),1);                          % preallocate new voxel RDM field
            % loop through voxels
            for v = 1:size(V.all,1)                     
                RDM         = rsa_squareRDM(V.all(v,:));                    % get RDM (all K conds)
                G           = -0.5*(H*RDM*H);                               % make G w/ centering matrix
                avgDist     = sum((C*G).*C,2)';                             % est. new distances
                LDC{s}(v,1) = nansum(avgDist)/size(C,1);                    % average across new distances for each voxel
            end
            clear V.all
            
            % % (4). Make output structure
            Ys.LDC = reshape(LDC{s},xVox,yVox,nslices);                     % re-arrange voxels to original space
            Ys.SN  = s;
            Ys.C   = C;
            Ys.dim = vol(1).dim;
            Ys.dt  = vol(1).dt;
            Ys.mat = vol(1).mat;
            Y = addstruct(Y,Ys);
            clear vol                                                      
        end
        varargout = {Y};
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

        if (size(Y,1)==32) % foot regressor/error regressor
            hold on;
            plot3(Y(32,1),Y(32,2),Y(32,3),'o','MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]);
            hold off;
        end
         % rest crosshairs
        hold on;
        plot3(0,0,0,'+','MarkerFaceColor',[0.75, 0, 0.75],'MarkerEdgeColor',[0.75, 0, 0.75],'MarkerSize',8);
        hold off;
        axis equal;
        
        %__________________________________________________________________
          
    case '0' % ------------ BEHA: behavioural force data cases. ----------- % These functions are from fdf2/3 so need editing for your paradigm
    case 'BEHA_get_forces'                                                  % Harvest figner force data from fingerbox force traces
        % harvest pressing force data for each trial
        sn       = 1:9;
        vararginoptions(varargin,{'sn'});
        cwd = pwd;
        A = []; % A = all subject force data
        for s = sn;
            cd(behavDir);
            [D,ND] = fivedigitFreq2_subj(sprintf('%d',s),0,1,1,'prefix','fdf');
            if s==2
                %ND = getrow(ND,ND.BN~=7); % doesn't work this way- 
                D  = getrow(D,D.BN~=7);
            end
            % D  = force data of cued fingers
            % ND = force data of non-cued fingers
            D.SN = ones(length(D.BN),1).*s;
            save(fullfile(behavDir,sprintf('fdf2_forces_s%d.mat',s)),'D');%,'ND');
            A = addstruct(A,D);
        end
        save(fullfile(behavDir,'fdf2_allSubjForces.mat'),'A');
        cd(cwd);
    
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
    case 'WRAPPER_dicom_import'                                             % imports dicoms of various data types for subjects specified w/ option 'sn'
        % Converts dicom to nifti files w/ spm_dicom_convert.
        % Comment out series_type not required.
        vararginoptions(varargin,{'sn'});
        for s = sn
            % Did we acquire this subject's anatomical during this scan
            % session? If not, don't process.
            if ~isempty(anatNum{s})
                fprintf('\nImporting ANATOMICAL runs- subj %d\n',s);
                pp1_imana('PREP_dicom_import','sn',s,'series_type','anatomical');
            end
            if ~isempty(fscanNum{s})
                fprintf('\nImporting FUNCTIONAL runs- subj %d\n',s);
                pp1_imana('PREP_dicom_import','sn',s,'series_type','functional');
            end
%             if ~isempty(pscanNum{s})
%                 fprintf('\nImporting PHASE runs- subj %d\n',s);
%                 pp1_imana('PREP_dicom_import','sn',s,'series_type','phase');
%             end
            if ~isempty(fieldNum{s})
                fprintf('\nImporting FIELDMAP runs- subj %d\n',s);
                pp1_imana('PREP_dicom_import','sn',s,'series_type','fieldmap');
            end
        end
    case 'PREP_dicom_import'                                                % STEP 1.1   :  Import functional/anatomical dicom series: enter sn
        % converts dicom to nifti files w/ spm_dicom_convert
        series_type = 'functional';
        vararginoptions(varargin,{'sn','series_type'});
        cwd = pwd;
        switch series_type
            case 'functional'
                seriesNum = fscanNum;
            case 'anatomical'
                seriesNum = anatNum;  
            case 'fieldmap'
                seriesNum = fieldNum;
        end
        
        cd(fullfile(dicomDir,subj_name{sn}));
        % For each series number of this subject (in 'Subject Things')
        for i=1:length(seriesNum{sn})
            r     = seriesNum{sn}(i);
            % Get DICOM FILE NAMES
            folder = fullfile(dicomDir,subj_name{sn},[sprintf('%4.4d',r)]);
            cd(folder)
            DIR   = dir(sprintf('%s.%4.4d.*.dcm',DicomName{sn},r));  
            Names = vertcat(DIR.name);
            % Convert the dicom files with these names.
            if (~isempty(Names))
                % Load dicom headers
                HDR=spm_dicom_headers(Names,1);  
                % Make a directory for series{r} for this subject.
                % The nifti files will be saved here.
                dirname = fullfile(dicomDir,subj_name{sn},sprintf('series%2.2d',r));
                dircheck(dirname);
                % Go to the dicom directory of this subject
                cd(dirname);
                % Convert the data to nifti
                spm_dicom_convert(HDR,'all','flat','nii');                  
                cd ..
            end
            display(sprintf('Series %d done \n',seriesNum{sn}(i)))
        end
        % Display verbose messages to user. 
        % Lazy and won't include none-verbose version here.
        switch series_type
            case 'functional'
                fprintf('Subject %02d functional runs imported. Copy the unique .nii name for subj files and place into ''Subject Things''.\n',sn)
            case 'anatomical'
                fprintf('Anatomical runs have been imported for subject %d.\n',sn); 
                fprintf('Please locate the T1 weighted anatomical img. Copy it to the anatomical folder.\n')
                fprintf('Rename this file to ''s%02d_anatomical_raw.nii'' in the anatomical folder.\n',sn); 
            case 'fieldmap'
                fprintf('Subject %02d fieldmaps imported. Copy the unique .nii name for subj files and place into ''Subject Things''.\n',sn)
        end
        cd(cwd); 
    case 'WRAPPER_preprocess1'                                              % Wrapper for steps 1.2-1.7 preprocessing for a subject     
        % need to have check_time, dicom_import, anat_dicom_import done prior to this step.
        fieldmap_correct = 1; % do you apply fieldmap corrections?
        vararginoptions(varargin,{'sn','fieldmap_correct'});
        
        for s = sn
            % concatinate scans from the same run into one file, and remove
            % dummy images at the start
            pp1_imana('PREP_make_4dNifti','sn',s);
            % check if fieldmap correcting.
            if fieldmap_correct
                % If yes, unwarping, distortion, and realignment is done
                pp1_imana('PREP_fieldmap_make','sn',s);            % make fieldmap correction imgs
                pp1_imana('PREP_fieldmap_RealignUnwarp','sn',s);   % realign with fieldmap correction
            else
                % If no, only realignment is done
                pp1_imana('PREP_realign','sn',s);                  % realign without fieldmap correction
            end
            % At this point, files for each run are saved in 'imaging_dir_raw'
            % Move the files to a working directory and operate on those,
            % keeping the 'raw' versions untouched as backup.
            pp1_imana('PREP_move_data','sn',s);
            % Now, reslice the anatomical into LPI coordinates and centre
            % the origin at the anterior commissure (note you need to
            % provide these coordinates at the top of this script).
            if ~isempty(anatNum{s})
                % These anatomical processing cases are only applied if the
                % anatomical scan was acquired during this session. If you
                % already have this subject's processed anatomical, we
                % should skip this b/c it's already done.
                pp1_imana('PREP_reslice_LPI','sn',s);
                pp1_imana('PREP_centre_AC','sn',s);
            end
        end
    case 'PREP_make_4dNifti'                                                % STEP 1.2   :  Converts dicoms to 4D niftis out of your raw data files
        vararginoptions(varargin,{'sn'});
        % For each functional run
        for i = 1:length(fscanNum{sn})                                      
            outfilename = fullfile(imagingDirRaw,subj_name{sn},sprintf('%s_run_%2.2d.nii',subj_name{sn},i));
            % Create a 4d nifti of all functional imgs in this run.
            % Don't include the first few dummy scans in this 4d nifti.
            for j = 1:(numTRs{sn}(i)-numDummys)                                        
                P{j}=fullfile(dicomDir,subj_name{sn},sprintf('series%2.2d',fscanNum{sn}(i)),...
                    sprintf('f%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{sn},fscanNum{sn}(i),j+numDummys,j+numDummys));
            end;
            dircheck(fullfile(imagingDirRaw,subj_name{sn}))
            spm_file_merge(char(P),outfilename);
            fprintf('Run %d done\n',i);
        end
    case 'PREP_fieldmap_make'
        vararginoptions(varargin,{'sn'});

        prefix = '';
        runNames = {};
        for i = 1:numel(run{sn})
            runNames{end+1} = sprintf('_%02d',i);
        end
        
        vararginoptions(varargin,{'sn'});
        spmj_makefieldmap(baseDir, subj_name{sn}, runNames,'prefix',prefix);
    case 'PREP_fieldmap_RealignUnwarp'
        error('needs editing to account for different numTRs per run')
        vararginoptions(varargin,{'sn'});
        prefix  ='';
        runNames = {};
        for i = 1:numel(run{sn})
            runNames{end+1} = sprintf('_%02d',i);
        end
        
        vararginoptions(varargin,{'sn'});
        spmj_realign_unwarp_sess(baseDir,subj_name{sn},runNames,numTRs{sn},'prefix',prefix);
    case 'PREP_realign'                                                     % STEP 1.3   :  Realign functinoal runs
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
        
        % Appends prefix 'r' to realigned imgs.
        vararginoptions(varargin,{'sn'});

        cd(fullfile(imagingDirRaw,subj_name{sn}));
        data={};
        for i=1:length(fscanNum{sn});
            for j=1:(numTRs{sn}(i)-numDummys);
                data{i}{j,1}=sprintf('%s_run_%2.2d.nii,%d',subj_name{sn},i,j);
            end;
        end;
        spmj_realign(data);
        fprintf('Subj %d realigned\n',sn);


    %__________________________________________________________________
    case 'PREP_move_data'                                                   % STEP 1.4   :  Moves subject data from raw directories to working directories
        % Moves image data from imaging_dicom_raw into a "working dir":
        % imaging_dicom.
        vararginoptions(varargin,{'sn'});

        prefix='r';
        dircheck(fullfile(baseDir, 'imaging_data',subj_name{sn}))
        for r=1:length(run{sn});

            source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, sprintf('%s%s_run_%02d.nii',prefix,subj_name{sn},run{sn}(r)));
            dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, sprintf('%s%s_run_%02d.nii',prefix,subj_name{sn},run{sn}(r)));

            copyfile(source,dest);
            source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, sprintf('rp_%s_run_%02d.txt',subj_name{sn},run{sn}(r)));
            dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, sprintf('rp_%s_run_%02d.txt',subj_name{sn},run{sn}(r)));

            copyfile(source,dest);
        end;
        % naming convention changes between fieldmap corrected and
        % not-corrected mean epi img.
        if strcmp(prefix,'u')
            source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, sprintf('mean%s%s_run_%02d.nii',prefix,subj_name{sn},run{sn}(r)));
        else
            source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, sprintf('mean%s_run_%02d.nii',subj_name{sn},run{sn}(1)));
        end
        dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, sprintf('%smeanepi_%s.nii',prefix,subj_name{sn}));

        copyfile(source,dest);
        fprintf('Moved niftis to working dir.')

    %__________________________________________________________________
    case 'PREP_reslice_LPI'                                                 % STEP 1.5   :  Reslice anatomical image within LPI coordinate systems
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
        display 'Done'


    %___________
    case 'PREP_centre_AC'                                                   % STEP 1.6   :  Re-centre AC in anatomical image
        % Set origin of anatomical to anterior commissure (must provide
        % coordinates in section (4)).
        vararginoptions(varargin,{'sn'});

        img    = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii']);
        V               = spm_vol(img);
        dat             = spm_read_vols(V);
        V.mat(1:3,4)    = loc_AC{sn};
        spm_write_vol(V,dat);
        display 'Done'


    %_____
    case 'PREP_meanimage_bias_correction'                                   % STEP 1.7   :  Bias correct mean image prior to coregistration
        prefix = 'u';
        vararginoptions(varargin,{'sn'});
        prefix = getCorrectPrefix(sn,prefix,what);
        
        % make copy of original mean epi, and work on that
        source  = fullfile(baseDir, 'imaging_data',subj_name{sn},[char(prefix) 'meanepi_' subj_name{sn} '.nii']);
        dest    = fullfile(baseDir, 'imaging_data',subj_name{sn},['b' char(prefix) 'meanepi_' subj_name{sn} '.nii']);
        copyfile(source,dest);
        
        % bias correct mean image for grey/white signal intensities 
        P{1}    = dest;
        spmj_bias_correct(P);
    case 'PREP_coreg'                                                       % STEP 1.8   :  Coregister meanepi to anatomical image
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi
        %   image
        vararginoptions(varargin,{'sn'});
        prefix = 'r';
        cd(fullfile(anatomicalDir,subj_name{sn}));
        %coregtool;
        %keyboard();
        
        % (2) Automatically co-register functional and anatomical images
        
        J.ref = {fullfile(anatomicalDir,subj_name{sn},[ subj_name{sn}, '_anatomical','.nii'])};
        J.source = {fullfile(imagingDir,subj_name{sn},['r' char(prefix) 'meanepi_' subj_name{sn} '.nii'])}; 
        J.other = {''};
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
    case 'WRAPPER_preprocess2'                                              % Wrapper for steps 1.9-1.11 preprocessing for a subject     
        vararginoptions(varargin,{'sn'});
        
        for s=sn
            pp1_imana('PREP_make_samealign','sn',s);
            pp1_imana('PREP_segmentation','sn',s);
            pp1_imana('PREP_make_maskImage','sn',s);
            %display('Run spmj_checksamealign to check alignment of run_epi to rmean_epi')
            %spmj_checksamealign
        end
    case 'PREP_make_samealign'                                              % STEP 1.9   :  Align to first image (rbmeanepi_* of first session)
        prefix  = 'r';
        vararginoptions(varargin,{'sn'});

        cd(fullfile(imagingDir,subj_name{sn}));

        % Select image for reference
        P{1} = fullfile(imagingDir,subj_name{sn},sprintf('%srmeanepi_%s.nii',prefix,subj_name{sn}));

        % Select images to be realigned
        Q={};
        for r=1:numel(run{sn})
            for i=1:(numTRs{sn}(r)-numDummys);
                Q{end+1}    = fullfile(imagingDir,subj_name{sn},...
                    sprintf('%s%s_run_%02d.nii,%d',prefix, subj_name{sn},r,i));
            end;
        end;

        % Run spmj_makesamealign_nifti to bring all functional runs into
        % same space as realigned mean epis
        spmj_makesamealign_nifti(char(P),char(Q));
        fprintf('Done. Run spmj_checksamealign to check alignment.\n')
        spmj_checksamealign
    case 'PREP_segmentation'                                                % STEP 1.10  :  Segmentation & normalization
        vararginoptions(varargin,{'sn'});

        SPMhome=fileparts(which('spm.m'));

        J.channel.vols     = {fullfile(anatomicalDir,subj_name{sn},[subj_name{sn},'_anatomical.nii,1'])};
        J.channel.biasreg  = 0.001;
        J.channel.biasfwhm = 60;
        J.channel.write    = [0 0];

        J.tissue(1).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,1')};
        J.tissue(1).ngaus  = 1;
        J.tissue(1).native = [1 0];
        J.tissue(1).warped = [0 0];

        J.tissue(2).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,2')};
        J.tissue(2).ngaus  = 1;
        J.tissue(2).native = [1 0];
        J.tissue(2).warped = [0 0];

        J.tissue(3).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,3')};
        J.tissue(3).ngaus  = 2;
        J.tissue(3).native = [1 0];
        J.tissue(3).warped = [0 0];

        J.tissue(4).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,4')};
        J.tissue(4).ngaus  = 3;
        J.tissue(4).native = [1 0];
        J.tissue(4).warped = [0 0];

        J.tissue(5).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,5')};
        J.tissue(5).ngaus  = 4;
        J.tissue(5).native = [1 0];
        J.tissue(5).warped = [0 0];

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


    %__________________________________________________________________
    case 'PREP_make_maskImage'                                              % STEP 1.11  :  Make mask images (noskull and gray_only)
        vararginoptions(varargin,{'sn'});
        prefix = 'r';
        cd(fullfile(imagingDir,subj_name{sn}));

        nam{1}  = fullfile(imagingDir,subj_name{sn}, [prefix 'rmeanepi_' subj_name{sn} '.nii']);
        nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
        nam{3}  = fullfile(anatomicalDir, subj_name{sn}, ['c2' subj_name{sn}, '_anatomical.nii']);
        nam{4}  = fullfile(anatomicalDir, subj_name{sn}, ['c3' subj_name{sn}, '_anatomical.nii']);
        spm_imcalc_ui(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2')

        nam={};
        nam{1}  = fullfile(imagingDir,subj_name{sn}, [prefix 'rmeanepi_' subj_name{sn} '.nii']);
        nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
        spm_imcalc_ui(nam, 'rmask_gray.nii', 'i1>1 & i2>0.4')

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
    case 'SURF_freesurfer'                                                  % STEP 2.1
        vararginoptions(varargin,{'sn'});
        freesurfer_reconall(freesurferDir,subj_name{sn},fullfile(anatomicalDir,subj_name{sn},[subj_name{sn} '_anatomical.nii']));
    case 'SURF_xhemireg'                                                    % STEP 2.2   :  Cross-Register surfaces left / right hem
        vararginoptions(varargin,{'sn'});
        freesurfer_registerXhem({subj_name{sn}},freesurferDir,'hemisphere',[1 2]); % For debug... [1 2] orig
    case 'SURF_map_ico'                                                     % STEP 2.3   :  Align to the new atlas surface (map icosahedron)
        vararginoptions(varargin,{'sn'});
        freesurfer_mapicosahedron_xhem(subj_name{sn},freesurferDir,'smoothing',1,'hemisphere',[1:2]);
    case 'SURF_make_caret'                                                  % STEP 2.4   :  Translate into caret format
        vararginoptions(varargin,{'sn'});
        caret_importfreesurfer(['x' subj_name{sn}],freesurferDir,caretDir);
            
    case '0' % ------------ GLM: SPM GLM fitting. Expand for more info. ---
        % The GLM cases fit general linear models to subject data with 
        % SPM functionality.
        %
        % All functions can be called with ('GLM_processAll','sn',[Subj#s]).
        %
        % You can view reconstructed surfaces with Caret software.
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'WRAPPER_GLM'                                                  
        vararginoptions(varargin,{'sn','glm'});
        % You can call this case to do all the GLM estimation and contrasts.
        for s = sn
            for g = glm
                pp1_imana('GLM_make','sn',s,'glm',g);
                pp1_imana('GLM_estimate','sn',s,'glm',g);
                if glm==1
                    pp1_imana('GLM_contrastglm1','sn',s);
                elseif glm==2
                    pp1_imana('GLM_contrastglm2','sn',s);
                    pp1_imana('PSC_calcImgsglm2','sn',s);
                elseif glm==3
                    pp1_imana('GLM_contrastglm3','sn',s);
                    pp1_imana('PSC_calcImgsglm3','sn',s);
                elseif glm==4
                    pp1_imana('GLM_contrastglm4','sn',s);
                    pp1_imana('PSC_calcImgsglm4','sn',s);
                end
            end
        end
    case 'GLM_make'                                                         % STEP 3.1   :  Make the SPM.mat and SPM_info.mat files (prep the GLM)
        % makes the GLM file for each subject, and a corresponding aux.
        % SPM_info.mat file. The latter file contains nice summary
        % information of the model regressors, condition names, etc.
        glm = 1;
        vararginoptions(varargin,{'sn','glm'});
        % Set some constants.
        prefix		 = 'r';
        T			 = [];
        announcetime = 0.5;                                                 % length of task announce time- % THIS NEEDS ADJUSTING TO YOUR STUDY
        dur          = 0.5; % task duration (seconds)
        % Gather appropriate GLM presets.
        subj_hrfParams = {[2.7046 10.348 1 1 0.16446]};
        
        
        switch glm                                                          % THIS NEEDS ADJUSTING TO YOUR STUDY
            case 1
                hrf_params = []; % use defaults
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
                % set hrf params to those optimized in glm 2
                if sn==1
                    hrf_params = [];
                end
            case 3
                % model all chords, model error trials as separate regressor
                % define the 7 parameters of the HRF
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 32; % 31 chords and 1 error regressor
                % set hrf params to those optimized in glm 2
                if sn==1
                    hrf_params = [];
                end    
            case 4
                % model all chords, model error trials as separate regressor
                % define the 7 parameters of the HRF
                hrf_params = subj_hrfParams{sn};
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 32; % 31 chords and 1 foot regressor
                % set hrf params to those optimized in glm 2
                if sn==1
                    hrf_params = [];
                end  
            case 5
                % model all chords, model error trials as separate regressor
                % define the 7 parameters of the HRF
                hrf_params = [];
                hrf_cutoff = inf;
                cvi_type   = 'fast';
                numConds   = 33; % 31 chords, 1 error regressor (subsequent trial following error), & 1 foot regressor
                % set hrf params to those optimized in glm 2
                if sn==1
                    hrf_params = [];
                end  
        end
        % Load subject's .dat file (has info on each trial)
        D	 = dload(fullfile(behavDir,sprintf('pp1_%s.dat',subj_name{sn})));
        
        % Do some subject structure fields.
        J.dir 			 = {fullfile(glmDir{glm}, subj_name{sn})};
        J.timing.units   = 'secs';                                          % timing unit that all timing in model will be
        J.timing.RT 	 = TR_length;                                       % TR (in seconds, as per 'J.timing.units')
        J.timing.fmri_t  = 16;
        J.timing.fmri_t0 = 1;
        % Loop through runs. 
        for r = 1:numel(run{sn})                                            
            R = getrow(D,D.BN==r);
            for i = 1:(numTRs{sn}(r)-numDummys)                                    % get nifti filenames, correcting for dummy scancs
                N{i} = [fullfile(baseDir, 'imaging_data',subj_name{sn}, ...
                        sprintf('%s%s_run_%02d.nii,%d',prefix,subj_name{sn},run{sn}(r),i))];
                        %[prefix subj_name{sn},'_run_0',run{r},'.nii,',num2str(i)])];
            end;
            J.sess(r).scans = N;                                            % number of scans in run
            % Loop through conditions.
            for c = 1:numConds
                switch glm
                    case 1
                        idx	= logical(R.chordNum>0);     
                        J.sess(r).cond(c).name = 'movement_on';
                        % Correct start time for numDummys removed & announce time, & convert to seconds
                        J.sess(r).cond(c).onset    = [R.StartTimeMeas(idx)/1000 - J.timing.RT*numDummys+announcetime];    
                        J.sess(r).cond(c).duration = dur;                           % durations of task we are modeling (not length of entire trial)
                        J.sess(r).cond(c).tmod     = 0;
                        J.sess(r).cond(c).orth     = 0;
                        J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                    case 2
                        % include all trials regardless of judgement
                        idx	= find(R.chordNum==c); % find indx of all trials in run of that condition 
                        J.sess(r).cond(c).name = sprintf('D%d',R.chordNum(idx(1)));  % make condition name (for user readability)
                        S.isError = 0;

                        % Correct start time for numDummys removed & announce time, & convert to seconds
                        J.sess(r).cond(c).onset    = [R.StartTimeMeas(idx)/1000 - J.timing.RT*numDummys+announcetime];    
                        J.sess(r).cond(c).duration = dur;                           % durations of task we are modeling (not length of entire trial)
                        J.sess(r).cond(c).tmod     = 0;
                        J.sess(r).cond(c).orth     = 0;
                        J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                    case 3
                        % model error trials as separate regressor
                        if c<numConds
                            idx	= find(R.chordNum==c & R.isError==0); % find indx of all trials in run of that condition 
                            if ~isempty(idx)
                                % Correct start time for numDummys removed & announce time, & convert to seconds
                                J.sess(r).cond(c).onset    = [R.StartTimeMeas(idx)/1000 - J.timing.RT*numDummys+announcetime];    
                                J.sess(r).cond(c).duration = dur;                           % durations of task we are modeling (not length of entire trial)
                                J.sess(r).cond(c).tmod     = 0;
                                J.sess(r).cond(c).orth     = 0;
                                J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                                S.incl = 1;
                            else
                                % All trials of this cond type were errors in
                                % this block. Zero-pad regressor column.
                                % Correct start time for numDummys removed & announce time, & convert to seconds
                                J.sess(r).cond(c).onset    = 0;    
                                J.sess(r).cond(c).duration = 0;                          
                                J.sess(r).cond(c).tmod     = 0;
                                J.sess(r).cond(c).orth     = 0;
                                J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                                S.incl = 0;
                            end
                            S.isError = 0;
                        else
                            % these are error trials to exclude
                            idx	= find(R.isError==1);
                            % Correct start time for numDummys removed & announce time, & convert to seconds
                            J.sess(r).cond(c).onset    = [R.StartTimeMeas(idx)/1000 - J.timing.RT*numDummys+announcetime];    
                            J.sess(r).cond(c).duration = dur;                           % durations of task we are modeling (not length of entire trial)
                            J.sess(r).cond(c).tmod     = 0;
                            J.sess(r).cond(c).orth     = 0;
                            J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                            J.sess(r).cond(c).name = 'error_trials'; 
                            S.isError = 1;
                            S.incl    = 1;
                        end
                        J.sess(r).cond(c).name = sprintf('D%d',R.chordNum(find(R.chordNum==c,1)));  % make condition name (for user readability)
                    case 4
                        % model error trials as separate regressor
                        if c<numConds
                            idx	= find(R.chordNum==c); % find indx of all trials in run of that condition 
                            % Correct start time for numDummys removed & announce time, & convert to seconds
                            J.sess(r).cond(c).onset    = [R.StartTimeMeas(idx)/1000 - J.timing.RT*numDummys+announcetime];    
                            J.sess(r).cond(c).duration = dur;                           % durations of task we are modeling (not length of entire trial)
                            J.sess(r).cond(c).tmod     = 0;
                            J.sess(r).cond(c).orth     = 0;
                            J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                            
                            S.incl = 1;
                            S.isFoot  = 0;
                        else
                            % these are trials with foot movement
                            idx	= find(R.isError==1 & R.RT>0);
                            if ~isempty(idx)
                                % Correct start time for numDummys removed & announce time, & convert to seconds
                                J.sess(r).cond(c).onset    = [R.StartTimeMeas(idx)/1000 - J.timing.RT*numDummys + announcetime + R.RT(idx)/1000];    
                                J.sess(r).cond(c).duration = 1;                           % duration of foot movement- let's say 1 second
                                J.sess(r).cond(c).tmod     = 0;
                                J.sess(r).cond(c).orth     = 0;
                                J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                                J.sess(r).cond(c).name = 'footmovement_trials'; 
                                S.incl   = 1;
                            else
                                % No foot movement this block
                                J.sess(r).cond(c).onset    = 0;    
                                J.sess(r).cond(c).duration = 0;                          
                                J.sess(r).cond(c).tmod     = 0;
                                J.sess(r).cond(c).orth     = 0;
                                J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                                S.incl = 0;
                            end
                            S.isFoot = 1;
                        end
                        J.sess(r).cond(c).name = sprintf('D%d',R.chordNum(find(R.chordNum==c,1)));  % make condition name (for user readability)
                    case 5
                        error('not yet implemented')
                        % model error trials as separate regressor
                        if c<32
                            idx	= find(R.chordNum==c & R.isError==0); % find indx of all trials in run of that condition 
                            if ~isempty(idx)
                                % Correct start time for numDummys removed & announce time, & convert to seconds
                                J.sess(r).cond(c).onset    = [R.StartTimeMeas(idx)/1000 - J.timing.RT*numDummys+announcetime];    
                                J.sess(r).cond(c).duration = dur;                           % durations of task we are modeling (not length of entire trial)
                                J.sess(r).cond(c).tmod     = 0;
                                J.sess(r).cond(c).orth     = 0;
                                J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                                S.incl = 1;
                            else
                                % All trials of this cond type were errors in
                                % this block. Zero-pad regressor column.
                                % Correct start time for numDummys removed & announce time, & convert to seconds
                                J.sess(r).cond(c).onset    = 0;    
                                J.sess(r).cond(c).duration = 0;                          
                                J.sess(r).cond(c).tmod     = 0;
                                J.sess(r).cond(c).orth     = 0;
                                J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                                S.incl = 0;
                            end
                            S.isError = 0;
                        else
                            % these are error trials to exclude
                            idx	= find(R.isError==1);
                            % Correct start time for numDummys removed & announce time, & convert to seconds
                            J.sess(r).cond(c).onset    = [R.StartTimeMeas(idx)/1000 - J.timing.RT*numDummys+announcetime];    
                            J.sess(r).cond(c).duration = dur;                           % durations of task we are modeling (not length of entire trial)
                            J.sess(r).cond(c).tmod     = 0;
                            J.sess(r).cond(c).orth     = 0;
                            J.sess(r).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                            J.sess(r).cond(c).name = 'error_trials'; 
                            S.isError = 1;
                            S.incl    = 1;
                        end
                    J.sess(r).cond(c).name = sprintf('D%d',R.chordNum(find(R.chordNum==c,1)));  % make condition name (for user readability)
                
                end
                % Do some subject info for fields in SPM_info.mat.
                S.sn  = sn;
                S.run = run{sn}(r);
                if c<32
                    S.chord 		= R.chordNum(find(R.chordNum==c,1));                       
                    S.numDigits     = R.numDigits(find(R.chordNum==c,1)); 
                    S.targetForce   = R.forceNtarget(find(R.chordNum==c,1));
                    S.numStim       = R.numStim(find(R.chordNum==c,1));
                else
                    S.chord 		= c;                          
                    S.numDigits     = 0; 
                    S.targetForce   = 3;
                    S.numStim       = 5;
                end  
                S.tt			= c;
                S.regtype		= 'Task';
                T				= addstruct(T,S);
            end;
            % Add any additional aux. regressors here.
            J.sess(r).multi 	= {''};
            J.sess(r).regress 	= struct('name', {}, 'val', {});
            J.sess(r).multi_reg = {''};                                
            % Define high pass filter cutoff (in seconds): see glm cases.
            J.sess(r).hpf 		= hrf_cutoff;
        end;
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
    case 'GLM_contrastglm4'                                                 % STEP 3.3   :  Make t-stat contrasts for specified GLM estimates.
        % Make t-stat contrasts for specified GLM estimates.
        % enter sn, glm #
        % models each chord and also error trials
        vararginoptions(varargin,{'sn'});
        glm = 4;
        cwd = pwd;
        % Go to subject's directory
        cd(fullfile(glmDir{glm}, subj_name{sn}));
        load SPM;
        SPM = rmfield(SPM,'xCon');
        T   = load('SPM_info.mat');

        %_____t contrast for chords
        for d = 1:32
            con                = zeros(1,size(SPM.xX.X,2));
            con(:,T.tt==d)     = 1;
            con                = con/sum(con);
            SPM.xCon(d)        = spm_FcUtil('Set',sprintf('chord_%d',d), 'T', 'c',con',SPM.xX.xKXs);
        end;

        %_____t contrast overall chords
        con                    = zeros(1,size(SPM.xX.X,2));
        con(:,T.tt<32)         = 1;
        con                    = con/sum(con);
        SPM.xCon(33)           = spm_FcUtil('Set',sprintf('overall'), 'T', 'c',con',SPM.xX.xKXs);

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
        cwd = pwd;
        % Go to subject's directory
        cd(fullfile(glmDir{glm}, subj_name{sn}));
        load SPM;
        SPM = rmfield(SPM,'xCon');
        T   = load('SPM_info.mat');

        %_____t contrast for chords
        for d = 1:32
            con                = zeros(1,size(SPM.xX.X,2));
            con(:,T.tt==d & T.incl==1) = 1;
            con                = con/sum(con);
            SPM.xCon(d)        = spm_FcUtil('Set',sprintf('chord_%d',d), 'T', 'c',con',SPM.xX.xKXs);
        end;

        %_____t contrast overall chords
        con                    = zeros(1,size(SPM.xX.X,2));
        con(:,T.tt>0 & T.incl==1) = 1;
        con                    = con/sum(con);
        SPM.xCon(33)           = spm_FcUtil('Set',sprintf('overall'), 'T', 'c',con',SPM.xX.xKXs);

        %____do the constrasts
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
    case 'PSC_calcImgsglm4'                                                     % calculates % signal change for digits of both trial types vs. rest (based on betas/baseline in each run)
        % calculate psc for all digits vs. rest - based on betas from glm 1    
        sn  = 1;
        vararginoptions(varargin,{'sn'});
        
        glm = 4;
        
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
            end;
            for con=1:31   % 31 chords
                P{numB+1}=sprintf('con_%04d.nii',con);
                outname=sprintf('psc_%02d.nii',con); % ,subj_name{s}
                formula = '100.*%f.*i%1.0f./((';
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
            end;
            fprintf('Subject %d: %3.3f\n',s,h);
        end;
    case 'PSC_calcImgsglm2'                                                     % calculates % signal change for digits of both trial types vs. rest (based on betas/baseline in each run)
        % calculate psc for all digits vs. rest - based on betas from glm 1    
        sn  = 1;
        vararginoptions(varargin,{'sn'});
        
        glm = 2;
        
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
            end;
            for con=1:31   % 31 chords
                P{numB+1}=sprintf('con_%04d.nii',con);
                outname=sprintf('psc_%02d.nii',con); % ,subj_name{s}
                formula = '100.*%f.*i%1.0f./((';
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
            end;
            fprintf('Subject %d: %3.3f\n',s,h);
        end;
    case 'PSC_calcImgsglm3'                                                     % calculates % signal change for digits of both trial types vs. rest (based on betas/baseline in each run)
        % calculate psc for all digits vs. rest - based on betas from glm 3
        % Skips over data from runs where all trials of a certain task type
        % were modeled as errors.
        sn  = 1;
        vararginoptions(varargin,{'sn'});
        
        glm = 3;
        
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
            end;
            for con=1:32   % 31 chords + error regressor
                ii = 1;
                numMissing = 0;
                P{numB+1}=sprintf('con_%04d.nii',con);
                outname=sprintf('psc_%02d.nii',con); % ,subj_name{s}
                formula = '100.*%f.*i%1.0f./((';
                for i = 1:numB
%                     if ~T.incl(logical(T.run==i & T.tt==con))
%                         keyboard
%                     end
                    
                    if i==1 && T.incl(logical(T.run==i & T.tt==con)) % first run and modeled as task
                        fadd = sprintf('i%1.0f',i);
                        ii = ii+1;
                    elseif i~=numB && T.incl(logical(T.run==i & T.tt==con)) % subsequent runs and modeled as task
                        fadd = sprintf('+i%1.0f',i);
                        ii = ii+1;
                    elseif i==numB && T.incl(logical(T.run==i & T.tt==con)) % last run and modeled as task
                        fadd = sprintf('+i%1.0f)/',i);
                        ii = ii+1;
                    elseif i==numB && ~T.incl(logical(T.run==i & T.tt==con)) % last run and not modeled as task
                        fadd = '/';
                    else
                        fadd = '';
                        ii = ii+1;
                        numMissing = numMissing +1;
                    end
                    formula = [formula fadd];
                end
                formula = [formula num2str(ii-1-numMissing) ')'];
                formula = sprintf(formula,h,ii);

                spm_imcalc_ui(P,outname,formula,{0,[],spm_type(16),[]});        % Calculate percent signal change
            end;
            fprintf('Subject %d: %3.3f\n',s,h);
        end;
        
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
            end
            %pp1_imana('SEARCH_map_contrast','sn',s,'glm',glm);
            %pp1_imana('SEARCH_map_LDC','sn',s,'glm',glm);
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
        glm = 2;
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
        % make index vectors
        conditionVec  = kron(ones(numel(run{sn}),1),[1:numConds]');
        partitionVec  = kron(run{sn}',ones(numConds,1));
        % go to subject's glm directory 
        cd(fullfile(glmDir{glm},subj_name{sn}));
        % load their searchlight definitions and SPM file
        L = load(fullfile(anatomicalDir,subj_name{sn},sprintf('%s_searchlight_120.mat',subj_name{sn})));
        load SPM;
        SPM  = spmj_move_rawdata(SPM,fullfile(imagingDir,subj_name{sn}));

        name = sprintf('%s_glm%d',subj_name{sn},glm);
        % run the searchlight
        rsa.runSearchlightLDC(L,'conditionVec',conditionVec,'partition',partitionVec,'analysisName',name,'idealBlock',block);
        cd(cwd);
    case 'SEARCH_map_contrast'                                              % STEP 4.3   :  Averaged LDC values for specified contrasts
        % Calls 'MISC_SEARCH_calculate_contrast'
        sn  = 1;
        glm = 1;
        con = {'avg','1digitChords','2digitChords','3digitChords','4digitChords','numDigits'};
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
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm);
                    % average across all paired dists (excluding error/foot conds)
                    Y.LDC   = vdat(:,:,:,gidx);
                    Y.LDC   = nanmean(Y.LDC,4);
                    % prep output file
                    Y.dim   = vol(1).dim;
                    Y.dt    = vol(1).dt;
                    Y.mat   = vol(1).mat;    
                case '1digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',1);
                    % average across all paired dists (excluding error/foot conds)
                    Y.LDC   = vdat(:,:,:,gidx);
                    Y.LDC   = nanmean(Y.LDC,4);
                    % prep output file
                    Y.dim   = vol(1).dim;
                    Y.dt    = vol(1).dt;
                    Y.mat   = vol(1).mat;    
                case '2digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',2);
                    % average across all paired dists (excluding error/foot conds)
                    Y.LDC   = vdat(:,:,:,gidx);
                    Y.LDC   = nanmean(Y.LDC,4);
                    % prep output file
                    Y.dim   = vol(1).dim;
                    Y.dt    = vol(1).dt;
                    Y.mat   = vol(1).mat;    
                case '3digitChords'
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',3);
                    % average across all paired dists (excluding error/foot conds)
                    Y.LDC   = vdat(:,:,:,gidx);
                    Y.LDC   = nanmean(Y.LDC,4);
                    % prep output file
                    Y.dim   = vol(1).dim;
                    Y.dt    = vol(1).dt;
                    Y.mat   = vol(1).mat;    
                case '4digitChords' 
                    gidx    = pp1_imana('GET_idxPerNumDigits','glm',glm,'numDigits',4);
                    % average across all paired dists (excluding error/foot conds)
                    Y.LDC   = vdat(:,:,:,gidx);
                    Y.LDC   = nanmean(Y.LDC,4);
                    % prep output file
                    Y.dim   = vol(1).dim;
                    Y.dt    = vol(1).dt;
                    Y.mat   = vol(1).mat;    
                case 'numDigits'
                    gidx    = pp1_imana('GET_idxAcrossNumDigits','glm',glm);
                    % average across all paired dists (excluding error/foot conds)
                    Y.LDC   = vdat(:,:,:,gidx);
                    Y.LDC   = nanmean(Y.LDC,4);
                    % prep output file
                    Y.dim   = vol(1).dim;
                    Y.dt    = vol(1).dt;
                    Y.mat   = vol(1).mat; 
                case {'ttt'}
                    switch con{c}
                        case 'numDigits'
                            C = kron(ones([1,4]),rsa.util.pairMatrix(5));  % digit contrast  
                            C = C/4; % scale contrast vectors across speeds
                        case 'speed'
                            C = kron(rsa.util.pairMatrix(4),ones(1,5));    % speed contrast
                            C = C/5; % scale contrast vectors across digits
                    end
                    % get new distances from searchlight results
                    LDC_file  = fullfile(glmDir{glm},subj_name{sn},sprintf('s0%d_glm%d_LDC.nii',sn,glm));
                    Y         = fivedigitFreq2_imana('MISC_SEARCH_calculate_contrast','sn',sn,'glm',glm,'C',C,'file',LDC_file);
            end

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

        for s=sn
            caretSubjDir = fullfile(caretDir,[atlasA subj_name{s}]);
            file         = fullfile(glmDir,subj_name{s},'mask.nii,1');
                
            for h=1:2 % ma
                D1 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI_3.paint']));   % sBa1 sBa2 sB3a sB3b
                D2 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROI.paint']));     % premade ROIs from fsaverage_sym: M1 and S1
                D3 = caret_load(fullfile(caretDir,atlasname,hemName{h},['ROIsm.paint']));   % premade ROIs from fsaverage_sym: brodmann areas
                for r = 1:numregions
                    if r<5; D = D1; rr = r;         % premade ROIs from fsaverage_sym: brodmann areas
                    elseif r<7 D = D2; rr = r-4;    % premade ROIs from fsaverage_sym: S1 and M1
                    else D = D3; rr = r-6; end      % sBa1 sBa2 sB3a sB3b: brodmann areas cut with S1 boundary
                    idx = r+(h-1)*numregions;
                    R{idx}.type     = 'surf_nodes';
                    R{idx}.location = find(D.data(:,1)==rr);
                    R{idx}.white    = fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                    R{idx}.pial     = fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                    R{idx}.topo     = fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                    R{idx}.linedef  = linedef;
                    R{idx}.image    = file;
                    R{idx}.name     = [subj_name{s} '_' regname{r} '_' hem{h}];
                    R{idx}.flat     = fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.FLAT.coord']);
                end    
            end;
            %R = region_calcregions(R,'exclude',[1,6; 1,2; 2,3; 3,4; 5,6; 6,11; 11,10; 10,9; 9,7; 7,8],'exclude_thres',0.75);
            %R = region_calcregions(R,'exclude',[1,6; 2,6; 3,6; 5,6],'exclude_thres',0.75);
            R = region_calcregions(R,'exclude',[2,6; 5,6],'exclude_thres',0.75);
            cd(regDir);
            save(['regions_' subj_name{s} '.mat'],'R');
            fprintf('\n %s done\n',subj_name{s})
            clear R
        end;
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
        P0  = [1  6  1   1    6]';               % Default parameters for the SPM hrf
        LB  = [0 0 0 0 0]';%[0  9   0.2 0.2  3  -2 0.2]';     
        UB  = [7  16  10   10    10]'; 
        duration = 1;
        onsetshift = 0;
        fit = [1,2,5]'; % hrf parameter to be fitted
        roi = 1;
        eCriteria = 0.97;
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
                end;
                SPM.Sess(r).U=spm_get_ons(SPM,r);
            end;
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
                for i=1:size(D.block,1);
                    D.y_hat(i,:)=cut(y_hat,pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_res(i,:)=cut(y_res,pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_adj(i,:)=D.y_hat(i,:)+D.y_res(i,:);
                end;
                D.regType   = ones(size(D.event,1),1)*r;
                D.regSide   = ones(size(D.event,1),1)*hemi;
                D.sn        = ones(size(D.event,1),1)*s;
                Ts          = addstruct(Ts,D);
            end
            warning on
        end
        
        
        % plot
%         titles = ['Delay','','','','','Onset','Duration'];
%         figure('name',[what,'_parameters']);
%         subplot(2,2,1)
%         lineplot(T.regSide,T.P(:,1),'split',T.SN); title(titles(fit(1)))
%         subplot(2,2,2)
%         try;lineplot(T.regSide,T.P(:,2),'split',T.SN); title(titles(fit(2)));catch;end;
%         subplot(2,2,3)
%         lineplot(T.regSide,T.R2,'split',T.SN); title('R^2')
%         subplot(2,2,4)
%         lineplot(T.regSide,T.Eratio,'split',T.SN); title('ErrorAfter/ErrorBefore')
        
        figure('name',[what,'_timeseries'])
        title(hemName{hemi});
        for reg = roi
            subset = (Ts.regType==reg & Ts.regSide==hemi);
            traceplot([-pre:post],Ts.y_adj,'errorfcn','stderr',...
                'split',[Ts.regType],'leg',regname(reg),'subset',subset); % ,
            hold on;
            traceplot([-pre:post],Ts.y_hat,'linestyle','--',...
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
    case 'ROI_getTimeseries'                                                % (optional) :  Harvest ROI timeseries for specified region.
        % Use this and 'ROI_plot_timeseries' to ensure good GLM fits with
        % measured BOLD in rois.
        
        % Defaults
        sn  = 1;
        glm = 1;
        roi = 1:11;
        vararginoptions(varargin,{'sn','glm','roi'});
        pre  = 4;                                                                  % how many TRs before trial onset (2.8 secs)
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
                end
        end
        save(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)),'-struct','T');
        display(sprintf('Done %s (region # %d) for glm %d \n',reg_title{reg},reg,glm));
        
        %__________________________________________________________________
    case 'ROI_plotTimeseriesglm1'                                               % (optional) :  Plots timeseries for specified region (avg. across digits, separated by pressing speed)
        glm = 1;
        sn  = 1;
        roi = 5;
        vararginoptions(varargin,{'sn','glm','roi'});
        
        numDigits = stimulationChords;
        numDigits = sum(numDigits,2)';
        
        figure('Color',[1 1 1]);
        D=load(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)));
        i = 1;
        for s=sn
            T=getrow(D,D.sn==s & D.roi==roi); % only plot chord regressors
            
                % plot timeseries
                plt.trace([-4:20],T.y_adj,'errorfcn','stderr');
                hold on;
                traceplot([-4:20],T.y_hat,'linestyle',':','linewidth',2);
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
        glm = 1;
        sn  = 1;
        roi = 5;
        vararginoptions(varargin,{'sn','glm','roi'});
        
        numDigits = stimulationChords;
        numDigits = sum(numDigits,2)';
        
        figure('Color',[1 1 1]);
        D=load(fullfile(regDir,sprintf('glm%d_reg_timeseries.mat',glm)));
        i = 1;
        for s=sn
            T=getrow(D,D.sn==s & D.roi==roi); % only plot chord regressors
            for d = 1:5 % plot separately for numDigits
                subplot(length(sn),6,i)
                % plot timeseries
                plt.trace([-4:20],T.y_adj,'errorfcn','stderr','subset',ismember(T.event,find(numDigits==d)));
                hold on;
                traceplot([-4:20],T.y_hat,'linestyle',':','linewidth',2,'subset',ismember(T.event,find(numDigits==d)));
                hold off;
                xlabel('TR');
                ylabel('adjusted activation');
                xlim([-4 11]);
                title(sprintf('subj: %d numD: %d  .', s,d));
                legend off
                i = i+1;
                
            end
            if glm>2
                subplot(length(sn),6,i)
                plt.trace([-4:20],T.y_adj,'errorfcn','stderr','subset',T.event==32);
                hold on;
                traceplot([-4:20],T.y_hat,'linestyle',':','linewidth',2,'subset',T.event==32);
                hold off;
                xlabel('TR');
                ylabel('adjusted activation');
                xlim([-4 11]);
                title('foot reg  .');
                legend off
                i = i+1;
            end
            plt.match('y');
            for j = i-6:i-1
                subplot(length(sn),6,j);
                drawline(0,'dir','vert');
                drawline(10,'dir','vert');
                drawline(15,'dir','vert');
                drawline(0,'dir','horz');
            end
            
        end
        
        
        %__________________________________________________________________
    case 'ROI_plotTimeseriesCond'                                               % (optional) :  Plots timeseries for specified region (avg. across digits, separated by pressing speed)
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
        roi = [1:11];
        vararginoptions(varargin,{'sn','glm','roi'});
        
        T=[];
            
        % harvest
        for s=sn % for each subj
            fprintf('\nSubject: %d\n',s) % output to user
            
            % load files
            load(fullfile(glmDir{glm}, subj_name{s},'SPM.mat'));  % load subject's SPM data structure (SPM struct)
            load(fullfile(regDir,sprintf('regions_%s.mat',subj_name{s})));          % load subject's region parcellation & depth structure (R)
            
            % add percent signal change imgs for subject
            Q = {}; 
            for q = 1:31; Q{q} = (fullfile(glmDir{glm}, subj_name{s}, sprintf('psc_%02d.nii',q))); end
            Q = spm_vol(char(Q));
            
            % TR img info
            V = SPM.xY.VY; 
            
            for r = roi % for each region
                % get raw data/psc for voxels in region
                Y = region_getdata(V,R{r});  % Data Y is N x P (P is in order of transpose of R{r}.depth)
                P = region_getdata(Q,R{r});
                % estimate region betas
                [betaW,resMS,SW_raw,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM,'normmode','runwise');
                S.betaW                   = {betaW};        % cells for voxel data b/c diff numVoxels
                S.betaUW                  = {bsxfun(@rdivide,beta,sqrt(resMS))};  
                S.raw_beta                = {beta};
                S.psc                     = {P};
                S.resMS                   = {resMS};
                S.depth                   = {R{r}.depth(R{r}.excl==0)'};
                S.sn                      = s;
                S.roi                     = r;
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
        end
        % save T
        save(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm)),'-struct','T'); 
        fprintf('\n')
    case 'ROI_stats'                                                        % STEP 5.4   :  Calculate stats/distances on activity patterns
        glm = 2;
        sn  = 1;
        roi = [1:22];
        vararginoptions(varargin,{'sn','glm','roi'});
        
        T = load(fullfile(regDir,sprintf('glm%d_reg_betas.mat',glm))); % loads region data (T)
        
        numDigits = stimulationChords;
        numDigits = sum(numDigits,2)';
        
        % output structures
        To = [];
        Td = [];
        
        if glm==2
            numConds = 31;
        elseif glm==3;
            numConds = 32;
        elseif glm==4 % 31 chords and 1 error regressor
            numConds = 32; %31 chords and 1 foot regressor
        end
        
        % do stats
        for s = sn % for each subject
            D = load(fullfile(glmDir{glm}, subj_name{s}, 'SPM_info.mat'));   % load subject's trial structure
            C0  = indicatorMatrix('identity',D.run);
            fprintf('\nSubject: %d\n',s)
            num_run = length(run{s});
            
            for r = roi % for each region
                S = getrow(T,(T.sn==s & T.roi==r)); % subject's region data
                fprintf('%d.',r)
                for L = 1:length(layers) % for each layer defined in 'layers'
                    L_indx = (S.depth{1} > layers{L}(1)) & (S.depth{1} < layers{L}(2)); % index of voxels for layer depth
                    betaW  = S.betaW{1}(:,L_indx); 
                    betaW_nmean = betaW(1:(numConds*num_run),:)-C0*pinv(C0)*betaW(1:(numConds*num_run),:); % run mean subtraction  
                    %beta   = S.beta{1}(:,L_indx);
                    % % Toverall structure stats
                    % crossval second moment matrix
                    [G,Sig]     = pcm_estGCrossval(betaW_nmean(1:(numConds*num_run),:),D.run,D.tt);
                    So.G        = rsa_vectorizeIPM(G);
                    So.G_wmean  = rsa_vectorizeIPM(pcm_estGCrossval(betaW(1:(numConds*num_run),:),D.run,D.tt));
                    So.sig      = rsa_vectorizeIPM(Sig);
                    % squared dissimilarities
                    %So.RDM_wmean= rsa.distanceLDC(betaW,D.run,D.tt);        % rdm crossvalidated, on patterns without run mean patterns removed
                    So.ldc      = rsa.distanceLDC(betaW_nmean,D.run,D.tt);  % rdm crossvalidated, patterns with run means removed
                    So.corr_dist= corr_crossval(pcm_makePD(G)); 
                    %So.rdm_nocv = distance_euclidean(betaW_nmean',D.tt)';   % rdm no CV, on patterns with run means removed
                    % PSC
                    So.psc       = mean(S.psc{1},2)';
                    So.psc_chord = [1:31]; % no error/foot imgs calculated for psc
                    So.psc_numD  = numDigits;
                    % Calculate avg. betas for each condition + intercepts
                    Q = [];
                    Q.raw_beta = S.raw_beta{1};
                    Q.tt = D.tt;
                    Q.tt(end+1:end+num_run) = numConds +1; % intercept betas
                    Q = tapply(Q,{'tt'},{'raw_beta','mean'});
                    So.avg_betas = mean(Q.raw_beta,2)';
                    So.avg_tt    = Q.tt';
                    % indexing fields
                    So.sn       = s;
                    So.roi      = r;
                    So.layer    = L;
                    So.numVox   = size(betaW,2);
                    So.regSide  = regSide(r);
                    So.regType  = regType(r);
                    To          = addstruct(To,So);
                    
                    % calc avg. chord patterns for each number of digits 
                    d = [];
                    d.betaW = betaW_nmean;
                    d.numDigits = D.numDigits;
                    d.run = D.run;
                    d.roi = ones(size(d.run)).*r;
                    d.chord = D.chord;
                    d = getrow(d,d.chord<32);
                    d5 = getrow(d,d.numDigits==5);
                    d = getrow(d,d.numDigits<5 & d.numDigits>0);
                    d = tapply(d,{'numDigits','run','roi'},{'betaW','mean'});
                    d = addstruct(d,d5);
                    d = rmfield(d,{'chord'});
                    dg = pcm_estGCrossval(d.betaW,d.run,d.numDigits);
                    % calc distance between avg patterns for 1 finger up to
                    % 5 finger chords:
                    td.ldc = rsa.distanceLDC(d.betaW,d.run,d.numDigits)';
                    td.corr_dist = corr_crossval(pcm_makePD(dg))';
                    td.distPair = [1:10]';
                    td.digitDiff = [1;2;3;4;1;2;3;1;2;1];
                    td.roi   = ones(10,1).*r;
                    td.sn    = ones(10,1).*s;
                    td.layer = ones(10,1).*L;
                    Td = addstruct(Td,td);
                    
                end; % each layer
            end; % each region
        end; % each subject

        % % save
        save(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)),'-struct','To');
        save(fullfile(regDir,sprintf('glm%d_reg_TnumDigits.mat',glm)),'-struct','Td');
        fprintf('\nDone.\n')
    
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
        layer = 1; % all voxels
        removeMean = 1;
        betaType = 3; % 1: raw, 2: uni-whitened, 3: multi-whitened
        vararginoptions(varargin,{'sn','glm','roi','removeMean','betaType','layer'});
        
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
                    % - restrict patterns from sepecific layer
                    % - index patterns for conditions and runs
                    % - exclude any foot/error/intercept regressors
                    % - calculate non-CV pattern consistency
                    
                    S = getrow(T,(T.sn==s & T.roi==r));
                    % get voxels for this layer
                    L_indx = (S.depth{1} > layers{layer}(1)) & (S.depth{1} < layers{layer}(2)); 
                    switch betaType 
                        case 1 % raw betas
                            betas = S.raw_beta{1}(:,L_indx);
                        case 2 % univariately prewhitened betas
                            betas = S.betaUW{1}(:,L_indx);
                        case 3 % multivariately prewhitened betas
                            betas = S.betaW{1}(:,L_indx);
                    end
                    % make vectors for pattern consistency func
                    condVec = kron(ones(numel(run{s}),1),conds);
                    partVec = kron(run{s}',ones(length(conds),1));
                    partVec(condVec==0) = 0;
                    % calculate the pattern consistency
                    rs.r2  = rsa_patternConsistency(betas,partVec,condVec,'removeMean',removeMean);
                    rs.sn  = s;
                    rs.roi = r;
                    rs.glm = g;
                    rs.layer = layer;
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
        layer = 1; % all voxels
        removeMean = 1;
        betaType = 3; % 1: raw, 2: uni-whitened, 3: multi-whitened
        vararginoptions(varargin,{'sn','glm','roi','removeMean','betaType','layer'});
        
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
                    % - restrict patterns from sepecific layer
                    % - index patterns for conditions and runs
                    % - exclude any foot/error/intercept regressors
                    % - calculate non-CV pattern consistency
                    
                    S = getrow(T,(T.sn==s & T.roi==r));
                    % get voxels for this layer
                    L_indx = (S.depth{1} > layers{layer}(1)) & (S.depth{1} < layers{layer}(2)); 
                    switch betaType 
                        case 1 % raw betas
                            betas = S.raw_beta{1}(:,L_indx);
                        case 2 % univariately prewhitened betas
                            betas = S.betaUW{1}(:,L_indx);
                        case 3 % multivariately prewhitened betas
                            betas = S.betaW{1}(:,L_indx);
                    end
                    % make vectors for pattern consistency func
                    condVec = kron(ones(numel(run{s}),1),conds);
                    partVec = kron(run{s}',ones(length(conds),1));
                    partVec(condVec==0) = 0;
                    % calculate the pattern consistency
                    [rs.r2,rs.r] = rsa_patternConsistency_crossval(betas,partVec,condVec,'removeMean',removeMean);
                    rs.sn  = s;
                    rs.roi = r;
                    rs.glm = g;
                    rs.layer = layer;
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
        layer = 'all';
        % Correlate distances of conditions in roi across subjs.
        % Regression line is not forced through origin so ignores distance 
        % scaling across subjects.
        vararginoptions(varargin,{'roi','glm','layer'});

        switch layer
            case 'all'
                L=1;
            case 'superficial'
                L=2;
            case 'deep'
                L=3;
        end
        
        D   = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));      
        D   = getrow(D,D.region==roi);
        D   = getrow(D,D.layer==L);
        Cs  = corr(D.RDM');
        
        varargout = {Cs};
    case 'ROI_MDS_overall'                                                  % (optional) :  Plots the scaled representational structure. 
        % enter region, glm #, sn (if desired)
        cplot = 'one';
        glm   = 1;
        layer = 1;
        roi   = 2; % default primary motor cortex    
        clrCode = 1; % if >0, color all chords with digit X red, and all other chords black.
        vararginoptions(varargin,{'roi','glm','cplot','layer','clrCode'});
        % cplot = 'all' to plot all 4 MDS figures (i.e. no contrast and 3 contrasts)- default
        % cplot = 'one'  to plot only no contrast MDS figure        

        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        T = getrow(T,T.layer==layer & T.roi==roi);
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
                    % account for error/foot regressors
                    CnumDigits = indicatorMatrix('identity',[numDigits;0]);
                    CnumDigits = bsxfun(@minus,CnumDigits,mean(CnumDigits,2));
                    Call   = eye(32)-ones(32)/32;
                else
                    CnumDigits = indicatorMatrix('identity',[numDigits]);
                    CnumDigits = bsxfun(@minus,CnumDigits,mean(CnumDigits,2));
                    Call   = eye(31)-ones(31)/31;
                end

                Y{1} = rsa_classicalMDS(IPM,'mode','IPM');
                Y{2} = rsa_classicalMDS(IPM,'mode','IPM','contrast',Call);
                Y{3} = rsa_classicalMDS(IPM,'mode','IPM','contrast',CnumDigits);
                
                figure('Color',[1 1 1]);
                subplot(1,3,1);
                pp1_imana('MISC_scatterplotMDS',Y{1}(1:32,1:3),'split',split,'label',[1:31]','fig',gca);
                title('no contrast');
                subplot(1,3,2);
                pp1_imana('MISC_scatterplotMDS',Y{2}(1:32,1:3),'split',split,'label',[1:31]','fig',gca);
                title('contrast: diff b/t all conds');
                subplot(1,3,3);
                pp1_imana('MISC_scatterplotMDS',Y{3}(1:32,1:3),'split',split,'label',[1:31]','fig',gca);
                title('contrast: num digits');
            case 'one' % only do and plot no contrast MDS
                Y{1} = rsa_classicalMDS(IPM,'mode','IPM');
                pp1_imana('MISC_scatterplotMDS',Y{1}(1:32,1:3),'split',split,'label',[1:31]');
        end
%         keyboard
    
    case 'ROI_splithalfPattReliability'                                          % plot w/in subj, w/in speed rdm reliability (Across two partitions), compare with across-speed correlations. Insights into RSA stability    
        % Splits data for each session into two partitions (even and odd runs).
        % Calculates correlation coefficients between each condition pair 
        % between all partitions.
        % Default setup includes subtraction of each run's mean
        % activity pattern (across conditions).
        glm = 4;
        roi = 5; % default roi
        sn  = 1;
        mean_subtract = 1; % subtract run means
        betaType = 'raw'; % use raw, not normalized betas
        conds    = 1:31; % restrict to specific conditions?
        % Correlate patterns across even-odd run splits within subjects.
        % Does correlation across all depths.
        vararginoptions(varargin,{'roi','glm','sn','mean_subtract','betaType','conds'});
        
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
            % cluster run numbers into partition splits. Here, is even-odd splits.
            oddRuns    = logical(mod(run{sn},2));
            partitions = [{run{s}(oddRuns)},{run{s}(~oddRuns)}];
            
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
                for i = 1:2
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
                end;
                
            end
        end;
        varargout = {Q};  
    

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
        layer = 1;
        vararginoptions(varargin,{'sn','glm','roi','layer'});
        
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
                d.layer     = v.*T.layer(i);
                d.numDigits = v.*dd;
                d.sn        = v.*T.sn(i);
                D = addstruct(D,d);
                d = [];
            end
        end
    
        D = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi) & ismember(D.layer,layer));
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
        layer = 1;
        vararginoptions(varargin,{'sn','glm','roi','layer'});
        
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
            d.layer     = v.*T.layer(i);
            d.sn        = v.*T.sn(i);
            D = addstruct(D,d);
            d = [];
        end
    
        D = getrow(D,ismember(D.sn,sn) & ismember(D.roi,roi) & ismember(D.layer,layer));
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
        layer = 1;
        fig = [];
        vararginoptions(varargin,{'sn','glm','roi','layer','fig'});
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
       
        D = [];
        d = [];
        v = ones(size(T.avg_tt,2),1);
        for i = 1:size(T.sn,1);
            d.avg_beta = T.avg_betas(i,:)';
            d.tt       = T.avg_tt(i,:)';
            d.layer = v.*T.layer(i);
            d.roi   = v.*T.roi(i);
            d.sn    = v.*T.sn(i);
            D = addstruct(D,d);
            d = [];
        end

        D = getrow(D,D.layer==layer & ismember(D.roi,roi) & ismember(D.sn,sn));
        % plot
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.tt,D.avg_beta,'split',D.roi);
        plt.labels('roi','avg raw beta')
        drawline(0,'dir','horz');
        
        varargout = {D};
    case 'FIG_RDMline_LDC'
        sn = 1;
        glm = 4;
        roi = 5;
        numDigits = 1;
        layer = 1;
        fig = [];
        vararginoptions(varargin,{'sn','glm','roi','numDigits','layer','fig'});
        
        % % LDC lineplot for specified roi(s).
        % - load distances form ROI_stats
        % - determine indicies of distance pairs from overall RDM
        % - harvest specified distances from indicies
        % - arrange into output structure & plot
        
        D = pp1_imana('GET_LDCperNumDigits','sn',sn,'glm',glm,'roi',roi,'layer',layer);
        D = getrow(D,D.numDigits==numDigits);
        % plot
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.distpair,D.ldc,'split',D.roi);
        plt.labels('distance pair','ldc^2')
        drawline(0,'dir','horz');
        
        varargout = {D};
    case 'FIG_pattReliabilitySingleFinger'
        % plot split-half reliabilities for patterns evoked by single finger
        % chords.
        % no layer-specific analysis ability in this case.
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
    case 'FIG_pattReliabilityAllChords'
        % plot split-half reliabilities for patterns evoked by single finger
        % chords.
        % no layer-specific analysis ability in this case.
        sn  = 1;
        glm = 4;
        roi = 5;
        fig = [];
        vararginoptions(varargin,{'sn','glm','roi','fig'});
        
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
        
    case 'FIG_LDCperNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn    = 1;
        glm   = 2;
        layer = 1; % all voxels
        roi = [1:4,6]; % [3a, 3b, 1, 2, M1]
        fig = [];
        vararginoptions(varargin,{'sn','roi','layer','glm','fig'});
        
        D = pp1_imana('GET_LDCperNumDigits','sn',sn,'glm',glm,'roi',roi,'layer',layer);
        D = tapply(D,{'sn','roi','layer','numDigits'},{'ldc','mean'});
        % plot
        %style.use('5shades')
        sty = style.custom(plt.helper.get_shades(5,'gray','decrease',10));
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.roi,D.ldc,'split',D.numDigits,'style',sty);
        plt.labels('roi','ldc^2 between chords with same # digits')
        drawline(0,'dir','horz');
        
        varargout = {D};
    case 'FIG_PSCperNumDigits'
        % plot distance between chords (split by number of fingers
        % stimulated in chord) across M1 and S1 subregions.
        sn    = 1;
        glm   = 2;
        layer = 1; % all voxels
        roi = [1:4,6]; % [3a, 3b, 1, 2, M1]
        fig = [];
        vararginoptions(varargin,{'sn','roi','layer','glm','fig'});
        
        % load data
        T = load(fullfile(regDir,sprintf('glm%d_reg_Toverall.mat',glm)));
        
        D = [];
        d = [];
        v = ones(31,1);
        for i = 1:size(T.sn,1);
            d.psc = T.psc(i,:)';
            d.chord = T.psc_chord(i,:)';
            d.numDigits = T.psc_numD(i,:)';
            d.layer = v.*T.layer(i);
            d.roi   = v.*T.roi(i);
            d.sn    = v.*T.sn(i);
            D = addstruct(D,d);
            d = [];
        end
        
        %style.use('5shades');
        sty = style.custom(plt.helper.get_shades(5,'gray','decrease',10));
        Dr = tapply(D,{'sn','roi','layer','numDigits'},{'psc','mean'});
        Dr = getrow(Dr,Dr.layer==layer & ismember(Dr.roi,roi) & ismember(Dr.sn,sn));
        % plot
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(Dr.roi,Dr.psc,'split',Dr.numDigits,'style',sty);
        plt.labels('roi','psc')
        drawline(0,'dir','horz');
        
        varargout = {Dr,D};
    case 'FIG_LDCacrossNumDigits'
        % plot distance between chords with differing number of digits.
        sn    = 1;
        glm   = 2;
        layer = 1; % all voxels
        roi = [1:4,6]; % [3a, 3b, 1, 2, M1]
        fig = [];
        vararginoptions(varargin,{'sn','roi','layer','glm','fig'});
        
        
        
        
        D = pp1_imana('GET_LDCacrossNumDigits','sn',sn,'glm',glm,'roi',roi,'layer',layer);
        D = tapply(D,{'sn','roi','layer'},{'ldc','mean'});
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
        glm   = 2;
        layer = 1; % all voxels
        roi = [1:4,6]; % [3a, 3b, 1, 2, M1]
        fig = [];
        vararginoptions(varargin,{'sn','roi','layer','glm','fig'});
        
        D = load(fullfile(regDir,sprintf('glm%d_reg_TnumDigits.mat',glm)));
        D = getrow(D,ismember(D.sn,sn) & D.layer==layer & ismember(D.roi,roi));
        % plot
        style.use('default');
        if isempty(fig); figure('Color',[1 1 1]); else fig; end
        plt.line(D.roi,D.corr_dist);
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
        
    case '0' % ------------ fingerpics: project patterns onto surface & take pictures.
        % You will absolutely need to edit these cases. 
        % These versions of the cases are frome ef1_imana
        % (extensionflexion).
        % 1 values of coordinates in xlims and ylims correspond to 1/10mm
        % distances.
        % So a range of 20 = 2mm distance on surface projection.
    case 'fingerpics'                                                       % Makes jpegs of finger activity patterns on cortical surface M1/S1
        sn  = 1;
        glm = 4;
        vararginoptions(varargin,{'sn','glm'});
        
        for s=sn;
            for g=glm;
                pp1_imana('surf_mapFingers','sn',s,'glm',g)
                pp1_imana('surf_fingerpatterns','sn',s,'glm',g)
            end
        end
    case 'surf_mapFingers'                                                % Map locations of finger patterns- run after glm estimation
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
                caret_save(fullfile(caretSDir,sprintf('%s_glm%d_hemi%d_finger.metric',subj_name{s},glm,h)),M);
            end
        end;   
    case 'surf_fingerpatterns'             % Make finger pattern jpegs
        %close all;
        sn  = 1;
        glm = 4;
        xlims = [7 -12]; % optimized for left-hemi (contralateral to movement)
        ylims = [-5 14];
        vararginoptions(varargin,{'sn','glm','xlims','ylims'});
        
        h = 1; % left hemi
        groupDir = [gpCaretDir filesep hemName{h} ];
        cd(groupDir);
        switch(h)
            case 1 % left hemi
                coord  = 'lh.FLAT.coord';
                topo   = 'lh.CUT.topo';
                %data  = 'lh.surface_shape';  
                xlims=[-4 15]; % may need to adjust locations for pics
                ylims=[-2 17];
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
        for i = 1:5
            subplot(1,5,i);
            [M,d]   = caret_plotflatmap('M',M,'col',i,'data',data,...
                        'border',B.Border,'bordersize',10,'topo',topo,'coord',coord);
            colormap('jet');
        end;
        
        mm = 4; % force colour scaling on patterns
        % loop through both figures and all subplots to: force colour
        % scaling and label conditions
        for i = 1:5
            subplot(1,5,i);
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
function prefix = getCorrectPrefix(sn,prefix,step)
    % Store subject-specific prefixes here for realigned images ('r'),
    % fieldmap corrected relaigned images ('u'), etc.
    % These are important for specific preprocessing stages.
if isempty(prefix)
    switch step
        case 'PREP_meanimage_bias_correction'
            p = {'','','r','','','','r','','','u'};
            prefix = p(sn);
        case {'PREP_move_data','PREP_coreg','PREP_make_samealign','PREP_make_maskImage','GLM_make'}
            p = {'','','r','','','','r','','','u'};
            prefix = p(sn);
    end
end    
end

function dircheck(dir)
% Checks existance of specified directory. Makes it if it does not exist.
% SA 01/2016
if ~exist(dir,'dir');
    %warning('%s didn''t exist, so this directory was created.\n',dir);
    mkdir(dir);
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

cd(cwd); % return user to their directory when function was called
end