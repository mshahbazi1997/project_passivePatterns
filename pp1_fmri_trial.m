function t = pp1_fmri_trial(mov,d,fig)
% Script to process force traces for fmri portions of passivePatterns1.
% mov = .mov data trace for trial
% d   = .dat file row for trial
% fig = boolean (t/f) to plot individual trial plot traces to figure

% 1. Extract data:
state   = mov(:,1);        % trial state (4 = stimulation, 6 = response)
time    = mov(:,4);        % time (ms)
Force   = mov(:,5:9);      % 5 = thumb, 9 = pinky

t.run   = d.BN;
t.trial = d.TN;
t.chordNum   = d.chordNum;
t.numDigits  = d.numDigits;
t.numStim    = d.numStim;
t.falseResp  = d.falseResponse;
t.isError    = d.isError;
t.RT         = d.RT;
t.forceStim  = d.targetForceStim;
t.stimulated = [d.d1 d.d2 d.d3 d.d4 d.d5];
t.time_stimOnset = [0 0 0 0 0];
t.peakF_raw  = [d.peakF_d1 d.peakF_d2 d.peakF_d3 d.peakF_d4 d.peakF_d5];
t.peakF_filt = [0 0 0 0 0];
t.meanF_filt = [0 0 0 0 0];

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
fthres  = 0.33;
Force   = smooth_kernel(Force,4);   % Smoothing with Gaussian kernel
%vF  = velocity_discr(Force); % first derivative of force (velocity)
for i = 1:5
    idx    = find(Force(:,i)>fthres);
    idx    = idx(idx>=Ia1 & idx<=Ia2);
    if ~isempty(idx)
        t.time_stimOnset(i) = time(idx(1)) - time(Ia1);
        t.peakF_filt(i)     = max(Force(idx,i),[],1);
        t.meanF_filt(i)     = mean(Force(idx,i));
    end
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