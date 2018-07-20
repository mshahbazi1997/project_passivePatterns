function varargout=pp1_maketgt(varargin); 
numStim     = 5;        % number of stimulations per trial
numReps     = 2;        % repetitions per chord
chordNum    = [1:31];   % 1:5 are flexion, 6:10 are extension
forceN      = 3;        % force applied to stimulated finger
numRests    = [6,1];    % number of rest periods during block
numFalse    = 5;        % number of incorrect chord presentations
lengthRests = [9000,15000]; % length of one rest period (ms)
cueTime     = 500;      % ms
stimTime    = 4000;     % time for finger stimulation (in ms)
respTime    = 1500;     % time to wait for chord response (in ms)
fbTime      = 500;      % feedback time per trial (ms)
ITI         = 1000;     % ms 
TR          = 1500;     % ms
dummyscans  = 2;        % number of dummy scans @ start of block
vararginoptions(varargin,{'numStim','numReps','forceN','numFalse',...
    'numRests','lengthRests','chordNum','cueTime','stimTime','respTime','fbTime',...
    'ITI','TR','dummyscans'});

% check rests
if sum((size(numRests)~=size(lengthRests)))>0
    error('pp1_maketgt: Please check numRests and lengthRests- should be same size')
end

numConds = length(chordNum);

%- - - - - - -
% make rest vecotr
%- - - - - - -
if (length(numRests)==1)
    lengthRests=ones(numRests,1)*lengthRests; 
else
    lrests      = lengthRests;
    lengthRests = [];
    for i = 1:length(lrests)
        lengthRests = [lengthRests, ones([1,numRests(i)]).*lrests(:,i)];
    end        
    numRests = sum(numRests);
end; 
% randomly order rests
lengthRests=lengthRests(randperm(numRests)); 
% determine where rests will be in tgt file (ensure rest is not before
% first trial)
atStart = 1;
while atStart
    % assign rest before these trials:
    rests = sort(sample_wor([1:(numReps*numConds)]',numRests,1)); 
    % check that rest is not assigned before the first trial:
    if isempty(find(rests==1, 1)); atStart=0; end
end

%- - - - - - -
% make trial vector
%- - - - - - -
trials = repmat(chordNum,1,numReps);                   % sequential trial list
trials = sample_wor(trials',numReps*numConds,1);       % pseudorandomize trial order
% determine which trials will have false chord presentations presented at
% response window. Ensure that only one trial for chord conditions are choosen. 
sameChord = 1;
while sameChord 
    % get randomly selected trials
    falseIdx = sample_wor(1:numReps*numConds,numFalse,1);
    % check no condition appears twice in false trials
    sameChord = pdist([trials(falseIdx)],@(x,y) x-y) == 0;
end
falseCue = zeros(numReps*numConds,1);
falseCue(falseIdx) = 1;

%-
% set stimulation chords
%-
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
r = 1;             % rest counter (used to index values in rest vector)
trialTime = cueTime + stimTime + respTime + fbTime;
for n=1:length(trials) 
    % (1) determine trial timing
    % adding rest time before this collection of movements? (not before first trial)
    if (isincluded(rests,n))
       if n==1; error('rest before first trial'); end
       x = x + lengthRests(r)/TR; % compensate for rest with TR counter
       r = r + 1; 
    end;     
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
    x = x + (trialTime + ITI)/TR;
end; 
varargout={T}; 