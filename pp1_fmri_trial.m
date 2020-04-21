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
movGood = (length(unique(state))~=1);
if movGood
    % 2. Indices for important experimental variables
    % - state timing events  
    Ia1 = find(state==4,1,'first'); % stimulation start
    Ia2 = find(state==4,1,'last');  % stimulation end
    idx = Ia1:Ia2;
    % 3. Estimating peak forces
    %digits  = [d.d1 d.d2 d.d3 d.d4 d.d5].*[1:5];
    %digits  = digits(digits>0);
    fthres  = 0.2;
    Force   = smooth_kernel(Force,4);   % Smoothing with Gaussian kernel
    %vF  = velocity_discr(Force); % first derivative of force (velocity)
else
    %fprintf('trial %02d run %02d is missing data\n',d.TN,d.BN);
end



for i = 1:5 % for each finger:
    
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
    t.peakF_filt     = t.peakF_raw;
    t.peakF_stims    = [nan nan];
    t.peakF_times    = [nan nan];
    
    if movGood
        t.peakF_filt = max(Force(idx,i),[],1);
        
        if t.stimulated
            t.time_stimOnset = time(find(Force(idx,i)>fthres,1) + Ia1-1) - time(Ia1);
            if isempty(t.time_stimOnset)
                t.time_stimOnset = nan;
            end
            t.peakF_filt     = max(Force(idx,i),[],1);
            % find times for peak forces of fingers (split into windows to find
            % peaks b/c peak estimate is more stable than toying with
            % minpeakdistance)
            %[n,~,b] = histcounts(1:size(idx,2),d.numStim);
            b = floor(length(idx)/d.numStim);
            pf = [];
            pt = [];
            for j = 1:d.numStim
                b1 = (1+b*(j-1));
                b2 = b*j;
                jidx = idx(b1:b2);
                try
                    [pf(j),pt(j)] = findpeaks(Force(jidx,i),'minpeakheight',fthres,'npeaks',1,'SortStr','descend');
                    t.peakF_stims(j) = pf(j);
                    t.peakF_times(j) = time(pt(j))' - time(Ia1); % correct peak stim times to start of stimulation epoch
                catch
                    keyboard
                end
            end
        end
    end
    
    
    T = addstruct(T,t);
end



% -------------------------------------------------------------------------
% Display trial
if (fig) %&& sum(ismember(d.TN,[2,7,9,11,6,30,43,20,27,15]))==1)
    time = time(Ia1:Ia2)./1000;
    time = time-ones(size(time)).*time(1);
    plot(time,Force(Ia1:Ia2,1),'Color',[0.2 0.14 0.53],'LineWidth',2); hold on
    plot(time,Force(Ia1:Ia2,2),'Color',[0.29 0.49 0.36],'LineWidth',2);
    plot(time,Force(Ia1:Ia2,3),'Color',[0.97 0 0],'LineWidth',2);
    plot(time,Force(Ia1:Ia2,4),'Color',[0 0.58 0.77],'LineWidth',2);
    plot(time,Force(Ia1:Ia2,5),'Color',[0.8 0.27 0.5],'LineWidth',2);
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