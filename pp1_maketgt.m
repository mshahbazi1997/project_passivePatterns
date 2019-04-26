function varargout=pp1_maketgt(varargin) 
numStim     = 2;        % number of stimulations per trial
numReps     = 2;        % repetitions per chord
chordNum    = 1:31;     % 1:5 are flexion, 6:10 are extension
forceN      = 2;        % force applied to stimulated finger
numFalse    = 3;        % number of incorrect chord presentations
cueTime     = 500;      % ms
stimTime    = 4000;     % time for finger stimulation (in ms)
respTime    = 1500;     % time to wait for chord response (in ms)
fbTime      = 500;      % feedback time per trial (ms)
TR          = 1500;     % ms
dummyscans  = 2;        % number of dummy scans @ start of block
vararginoptions(varargin,{'numStim','numReps','forceN','numFalse',...
    'chordNum','cueTime','stimTime','respTime','fbTime',...
    'TR','dummyscans'});

numConds  = length(chordNum);
numTrials = numReps*numConds;

%- - - - - - -
% determine random rest itis
%- - - - - - -
iti         = [1:1.5:10].*1000;                 % variable inter-trial-interval duration (in ms)
niti        = numel(iti);                       % how many iti types
prb         = geopdf(1:niti, .35);               % geometric probability distribution function with p = 0.3
pct         = prb./sum(prb);                    % proportion of trials per ITI type 
nt_iti      = round(numTrials * pct);  % how many trials per ITI type
% correcting because rest itis are between trials, so num itis = numTrials-1
% While we could correct earlier, this way we ensure that the shortest iti
% is the one affected, and not longer itis.

% Make randomly shuffled rest itis
nt_iti(1)   = nt_iti(1)-1;                      
lengthRests = [];
for i = 1:niti
    lengthRests = [lengthRests;ones(nt_iti(i),1).*iti(i)];
end
lengthRests = sample_wor(lengthRests,numTrials-1,1); 
lengthRests = [0;lengthRests];

%- - - - - - -
% make trial vector
%- - - - - - -
trials = repmat(chordNum,1,numReps);                   % sequential trial list
trials = sample_wor(trials',numTrials,1);       % pseudorandomize trial order
% determine which trials will have false chord presentations presented at
% response window. Ensure that only one trial for chord conditions are choosen. 
sameChord = 1;
while sameChord 
    % get randomly selected trials
    falseIdx = sample_wor(1:numTrials,numFalse,1);
    % check no condition appears twice in false trials
    sameChord = pdist([trials(falseIdx)],@(x,y) x-y) == 0;
end
falseCue = zeros(numTrials,1);
falseCue(falseIdx) = 1;

%- - - - - - -
% set stimulation chords
%- - - - - - -
chords = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;...             % singles            5
          1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;...                        % doubles (thumb)    4
          0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;...                                   % doubles            3
          0 0 1 1 0; 0 0 1 0 1;...                                              % doubles            2
          0 0 0 1 1;...                                                         % doubles            1
          1 1 1 0 0; 1 1 0 1 0; 1 1 0 0 1; 1 0 1 1 0; 1 0 1 0 1; 1 0 0 1 1;...  % triples (thumb)    6
          0 1 1 1 0; 0 1 1 0 1; 0 1 0 1 1; 0 0 1 1 1;...                        % triples            4
          1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 0 1 1 1; 0 1 1 1 1;                % quadruples         5
          1 1 1 1 1]';                                                          % all five           1


%- - - - - - -
% add timing to trial list (calculates timing based on TR counter: x)
%- - - - - - -
T = [];            % output target structure
x = dummyscans+1;  % set the TR counter to the first TR we start at (used to calculate timing)
trialTime = cueTime + stimTime + respTime + fbTime;
for n=1:numTrials
    % (1) determine trial timing
    % determine start time
    x = x + lengthRests(n)/TR; % compensate for rest with TR counter (first trial is zero iti before)
    t.startTime = (x-1)*TR; 
    t.startTR   = x;
    t.endTime   = t.startTime + trialTime;
    t.cueTime   = cueTime;
    t.stimTime  = stimTime;
    t.respTime  = respTime;
    t.fbTime    = fbTime;
    t.trialTime = trialTime;
    % (2) add condition information
    t.chordNum  = trials(n);
    t.numStim   = numStim;
    t.forceNtarget = forceN;
    t.falseResponse = falseCue(n);
    chordCon = chords(:,t.chordNum);
    % digits to stimulate
    t.d1 = chordCon(1);
    t.d2 = chordCon(2);
    t.d3 = chordCon(3);
    t.d4 = chordCon(4);
    t.d5 = chordCon(5);
    % digits to present response cue to
    if t.falseResponse
        while isequal(chordCon,[t.d1;t.d2;t.d3;t.d4;t.d5])
            if sum(chordCon)==5
                chordCon(sample_wor(1:5,1,1)) = 0;
            end
            chordCon = sample_wor(chordCon,5,1);
        end
    end
    t.c1 = chordCon(1);
    t.c2 = chordCon(2);
    t.c3 = chordCon(3);
    t.c4 = chordCon(4);
    t.c5 = chordCon(5);
    
    T = addstruct(T,t); 
    % update TR counter
    x = x + (trialTime)/TR;
end
T.endTime(end) = T.endTime(end) + 3000; % add some buffer seconds 
varargout={T}; 