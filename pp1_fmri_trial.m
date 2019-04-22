function T = pp1_fmri_trial(mov,d,fig)
% Script to process force traces for fmri portions of passivePatterns1.
% mov = .mov data trace for trial
% d   = .dat file row for trial
% fig = boolean (t/f) to plot individual trial plot traces to figure

% 1. Extract data:
state   = mov(:,1);        % trial state (4 = stimulation, 6 = response)
time    = mov(:,4);        % time (ms)
Force   = mov(:,5:9);      % 5 = thumb, 9 = pinky
T = [];
if length(unique(state))==1
    fprintf('trial %02d run %02d is missing data\n',d.TN,d.BN);
    return;
end

% 2. Indices for important experimental variables
% - state timing events  
Ia1 = find(state==4,1,'first'); % stimulation start
Ia2 = find(state==4,1,'last');  % stimulation end

% 3. Estimating peak forces
%digits  = [d.d1 d.d2 d.d3 d.d4 d.d5].*[1:5];
%digits  = digits(digits>0);
fthres  = 0.3;
Force   = smooth_kernel(Force,4);   % Smoothing with Gaussian kernel
%vF  = velocity_discr(Force); % first derivative of force (velocity)

for i = 1:5
    
    % prep output structure for trial
    t = [];
    t.run            = d.BN;
    t.trial          = d.TN;
    t.chordNum       = d.chordNum;
    t.numDigits      = d.numDigits;
    t.numStim        = d.numStim;
    t.falseResp      = d.falseResponse;
    t.isError        = d.isError;
    t.RT             = d.RT;
    t.forceStim      = d.targetForceStim;
    t.finger         = i;
    t.stimulated     = eval(sprintf('[d.d%d]',i));
    t.time_stimOnset = nan;
    t.peakF_raw      = eval(sprintf('[d.peakF_d%d]',i));
    t.peakF_filt     = nan;
    t.peakF_stims    = [nan nan];
    t.peakF_times    = [nan nan];
    
    % get times of stimulation state
    idx    = find(Force(:,i)>fthres);
    idx    = idx(idx>=Ia1 & idx<=Ia2); % cut stimulation times to be within range of stimulation state
    ftime  = time(idx);
    if ~isempty(idx) && t.stimulated
        t.time_stimOnset = time(idx(1)) - time(Ia1);
        t.peakF_filt     = max(Force(idx,i),[],1);
        % find times for peak forces of fingers (split into windows to find
        % peaks b/c peak estimate is more stable than toying with
        % minpeakdistance)
        [n,e,b] = histcounts(idx,d.numStim);
        pf = [];
        pt = [];
        for j = 1:d.numStim
            jidx = idx(b==j);
            try
                [pf(j),pt(j)] = findpeaks(Force(jidx,i),'minpeakheight',fthres,'npeaks',1,'SortStr','descend');
            catch ME
                keyboard
            end
            pt(j) = jidx(pt(j));
        end
        %[pf,pt] = findpeaks(Force(idx,i),'minpeakdistance',80,'minpeakheight',fthres,'npeaks',d.numStim,'SortStr','descend');
        t.peakF_stims = pf;
        t.peakF_times = time(pt)' - time(Ia1); % correct peak stim times to start of stimulation epoch
    end
    T = addstruct(T,t);
end



% -------------------------------------------------------------------------
% Display trial
if (fig) %&& sum(ismember(d.TN,[2,7,9,11,6,30,43,20,27,15]))==1)
    time = time(Ia1:Ia2)./1000;
    time = time-ones(size(time)).*time(1);
    plot(time,Force(Ia1:Ia2,1),'Color',[0.2 0.14 0.53],'LineWidth',3); hold on
    plot(time,Force(Ia1:Ia2,2),'Color',[0.29 0.49 0.36],'LineWidth',3);
    plot(time,Force(Ia1:Ia2,3),'Color',[0.97 0 0],'LineWidth',3);
    plot(time,Force(Ia1:Ia2,4),'Color',[0 0.58 0.77],'LineWidth',3);
    plot(time,Force(Ia1:Ia2,5),'Color',[0.8 0.27 0.5],'LineWidth',3);
    leg = {'thumb','index','middle','fourth','little'};
    legend(leg,'Location','NorthWest');   
    legend boxoff
    box off
    xlabel('stimulation time (seconds)');
    ylabel('force (N)');
    drawline(0,'dir','horz');
    xlim([time(1) time(end)]);
    hold off
    keyboard;
end;