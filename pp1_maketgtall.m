function varargout=pp1_maketgtall(varargin)
% maketgtall specific params
numRuns     = 15;        % number of target files to be made
sn          = {'s03'};  % subject number
outDir      = cd;       % output directory for target files
% target file params
numStim     = 2;        % number of stimulations per trial
numReps     = 2;        % repetitions per chord
numFalse    = 3;        % number of false response trials per block
chordNum    = [1:31];   % 1:5 are flexion, 6:10 are extension
forceN      = 3;        % force applied to stimulated finger
cueTime     = 500;      % ms
stimTime    = 4000;     % time for finger stimulation (in ms)
respTime    = 1500;     % time to wait for chord response (in ms)
fbTime      = 500;      % feedback time per trial (ms)
TR          = 1500;     % ms
dummyscans  = 2;        % number of dummy scans @ start of block
vararginoptions(varargin,{'numStim','numReps','forceN','numFalse',...
    'chordNum','cueTime','stimTime','trialTime',...
    'TR','dummyscans','numRuns','sn','outDir'});
if ~iscell(sn); sn = {sn}; end

D=[]; 
for s=1:length(sn)
    for i=1:numRuns 
        T=pp1_maketgt('numStim',numStim,...
                     'numReps',numReps,...
                     'numFalse',numFalse,...
                     'cueTime',cueTime,...
                     'stimTime',stimTime,...
                     'respTime',respTime,...
                     'fbTime',fbTime,...    
                     'chordNum',chordNum,...
                     'forceN',forceN,...
                     'TR',TR,...
                     'dummyscans',dummyscans); 
        dsave(fullfile(outDir,sprintf('fmri_%s_%d.tgt',sn{s},i)),T);
        D=addstruct(D,T); 
    end
end
varargout = {D};

