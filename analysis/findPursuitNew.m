% FUNCTION to find pursuit following a combination of the method described
%   in Blohm & % Lefèvre, 2010, JNP (also some old papers they cited), and
%   our old piece-wise linear fit method instead of the regression they used.
% *Requires changeDetectSlope.m, evalPWLslope.m, and ms2frames.m
% The original findPursuit.m works fine for "classical" conditions;
%   however, if your data is more tricky, could try this function instead.
% In general, the first step is to find a proper window which ideally contains
%   the fixation right before pursuit onset and pursuit initiation, then
%   do piecewise linear fit to the window. To find the proper window, first
%   find a saccade-free fixation baseline before target onset of a certain
%   baselineIntervalLength, then starting from pursuitSearchStart, looking
%   for when pursuit velocity consistently exceeds 4sd (or 3) of the baseline
%   mean (pursuitIntervalStart). The fitting window is then from target onset
%   to pursuitIntervalLength after pursuitIntervalStart; you could adjust
%   based on your needs.
%
% Could more accurately find pursuit onset by using a serials of restrictions:
% --pursuitIntervalLength: how long after the starting point of the pursuit
%   interval to include for the piecewise linear fit.
% --restrictions for the piecewise linear fit (see explanation in
%   changeDetectSlope.m): slopeRange (restrict slope1--slope for fixation period,
%   and slope2--slope for pursuit initiation),
%   cxRange (restrict time of the break point, aka pursuit onset)
% Restrict possible window for pursuit onset according to your specific
%   condition by defining cxRange
%
% General suggestions for adjusting the parameter values to accurately
%   detect pursuit onset (could be different depending on the specific data):
% If onset detected too early, try shorter pursuitIntervalLength (such as
%   50ms), or larger slope2 (larger lower bound)
% If onset detected too late, try smaller magnitude of slope1 (more strict restriction of its range)
%
% History
% Jan-04-2021, created by XW; xiuyunwu5@gmail.com
%
% input: trial --> structure containing relevant current trial information
% output: pursuit --> structure containing info about pursuit onset

function [pursuit] = findPursuitNew(trial)

if trial.log.eyeType==0 % fixation
    pursuit.onset = NaN;
else
    %% set parameters
    slopeRange = [-0.01 0.01; 0.01 inf]; % if you are using changeDetectRange.m, see explanations there
    % set the interval lengths in frames
    baselineIntervalLength = ms2frames(100);
    pursuitIntervalLength = ms2frames(50);
    
    % define the interval for fixation baseline, time in relation to target onset
    fixationIntervalStart = ms2frames(-150);
    fixationIntervalEnd = fixationIntervalStart + baselineIntervalLength-1;
    
    % define the search window for the starting point of the pursuit interval
    pursuitSearchStart = ms2frames(50); % when should we start looking for pursuit onset, x ms after stimulus onset
    pursuitSearchEnd = ms2frames(250);
    
    %% First, look for a pursuit interval that consistently exceeds a velocity larger than 4SD of the baseline
    % get the time points for the baseline
    baselineStart = trial.stim_onset + fixationIntervalStart;
    baselineEnd = trial.stim_onset + fixationIntervalEnd;
    % if a saccade overlaps with the baseline interval, search for earlier saccade-free phases
    while any(trial.saccades.onsets>=baselineStart & trial.saccades.onsets<=baselineEnd) || ...
            any(trial.saccades.offsets>=baselineStart & trial.saccades.offsets<=baselineEnd) ||...
            any(baselineStart>=trial.saccades.onsets & baselineEnd<=trial.saccades.offsets) % the last account for the interval being within a saccade
        baselineStart = baselineStart - 1;
        baselineEnd = baselineEnd - 1;
    end
    if baselineStart < 1 % try to search towards target onset if earlier saccade-free phases not available
        baselineStart = trial.stim_onset + fixationIntervalStart;
        baselineEnd = trial.stim_onset + fixationIntervalEnd;
        while any(trial.saccades.onsets>=baselineStart & trial.saccades.onsets<=baselineEnd) || ...
                any(trial.saccades.offsets>=baselineStart & trial.saccades.offsets<=baselineEnd) ||...
                any(baselineStart>=trial.saccades.onsets & baselineEnd<=trial.saccades.offsets)
            baselineStart = baselineStart + 1;
            baselineEnd = baselineEnd + 1;
            if baselineEnd > trial.stim_onset + ms2frames(50) % not later than 50ms after target onset; if still not working...
                baselineStart = -999;
                break
            end
        end
    end
    
    % confirm if we have a valid baseline interval
    if baselineStart == -999
        disp('Error of fixation interval in findPursuitNew.m!')
        pursuit.onset = NaN; % no pursuit onset
        pursuit.baselineStart = NaN;
        pursuit.baselineEnd = NaN;
        pursuit.meanBaseline = NaN;
        pursuit.sdBaseline = NaN;
        pursuit.pursuitIntervalStart = NaN;
        pursuit.pursuitIntervalEnd = NaN;
        return
    end
    
    % get the baseline mean and SD, in 2D
    baselineVelX = trial.eyeDX_filt(baselineStart:baselineEnd);
    baselineVelY = trial.eyeDY_filt(baselineStart:baselineEnd);
    baselineVel2D = sqrt(baselineVelX.^2 + baselineVelY.^2);
    meanBaseline = nanmean(baselineVel2D);
    sdBaseline = nanstd(baselineVel2D);
    
    % find the pursuit interval that is consistently 4SD away from baseline
    % mean; use the interpol velocity to ignore saccades for now
    pursuitWindowStart = trial.stim_onset + pursuitSearchStart;
    pursuitWindowEnd = trial.stim_onset + pursuitSearchEnd;
    pursuitVelX = trial.DX_interpolSac;
    pursuitVelY = trial.DY_interpolSac;
    pursuitVel2D = sqrt(pursuitVelX.^2 + pursuitVelY.^2);
    
    diffAbsVel = abs(pursuitVel2D(pursuitWindowStart:pursuitWindowEnd)-meanBaseline)>4*sdBaseline;
    changeDiff = diff(diffAbsVel);
    
    if all(changeDiff~=1) && ~any(diffAbsVel)
        % usually the case would be that the baseline interval is too noisy, or
        % the pursuit magnitude is small; try directly use the piecewise linear
        % fit in the window
        time = trial.stim_onset:trial.stim_onset+pursuitSearchEnd-1;
        % first use the saccade interpolated trace, deal with saccades later
        dataxy_tmp = sqrt( (trial.DX_interpolSac-nanmean(baselineVelX)).^2 + (trial.DY_interpolSac-nanmean(baselineVelY)).^2 );
        XY = dataxy_tmp(time);
        % run changeDetect.m to find the time of pursuit onset
        % Basically we assume the velocity curve from fixation to steady-state
        % pursuit to be a two-breakpoint piecewise linear model, and this
        % function can be used to find either one breakpoint (should only include
        % one breakpoint in XY to be accurate).
        %                      ---------------- steady-state
        %                     / initiation phase end
        %                    /
        % fixation__________/ pursuit onset
        if any(isnan(XY))
            pursuit.onset = NaN;
        else
            cxRange = [trial.stim_onset; time(end)]; % not really restricting the onsets now...
            [cx,cy,slope1,slope2] = changeDetectSlope(time, XY, cxRange, slopeRange);
            pursuit.onset = round(cx);
        end
        
        pursuit.pursuitIntervalStart = NaN;
        pursuit.pursuitIntervalEnd = NaN;
    else
        if all(changeDiff~=1) % if already much higher than fixation, start from 50ms and eventually use 0-200 ms
            potentialStartingPoints = 1;
        else
            % to make it more stable if necessary, try to find the first
            % starting point of an interval of at least 50 ms that
            % consistently exceeds the threshold
            idxAll = find(changeDiff==1);
            if length(idxAll)>1
                while idxAll(2)-idxAll(1)<100
                    idxAll(1) = [];
                    if length(idxAll)==1
                        break
                    end
                end
            end
            potentialStartingPoints = idxAll(1)+1;
        end
        pursuitIntervalStart = pursuitWindowStart + potentialStartingPoints-1;
        % this is the first time point when velocity consistently exceeds the
        % threshold of mean + 4*sd
        
        %% use piece-wise linear fit to find the exact break point
        % now we choose the interval to look for pursuit onset around this
        % point, 150 before to 150 after
        time = trial.stim_onset:(pursuitIntervalStart+pursuitIntervalLength);
        % first use the saccade interpolated trace, deal with saccades later
        dataxy_tmp = sqrt( (trial.DX_interpolSac-nanmean(baselineVelX)).^2 + (trial.DY_interpolSac-nanmean(baselineVelY)).^2 );
        XY = dataxy_tmp(time);
        % run changeDetect.m to find the time of pursuit onset
        % Basically we assume the velocity curve from fixation to steady-state
        % pursuit to be a two-breakpoint piecewise linear model, and this
        % function can be used to find either one breakpoint (should only include
        % one breakpoint in XY to be accurate).
        %                      ---------------- steady-state
        %                     / initiation phase end
        %                    /
        % fixation__________/ pursuit onset
        if any(isnan(XY))
            pursuit.onset = NaN;
        else
            cxRange = [trial.stim_onset; pursuitIntervalStart]; % restrict the onset to be before the starting point of the pursuit interval
            [cx,cy,slope1,slope2] = changeDetectSlope(time, XY, cxRange, slopeRange); % restrict the range of slopes
            pursuit.onset = round(cx);
        end
        
        pursuit.pursuitIntervalStart = pursuitIntervalStart; % for debugging...
        pursuit.pursuitIntervalEnd = pursuitIntervalStart+pursuitIntervalLength-1;
        
        %% using linear regression, could result in super early pursuit onsets, unstable
        %   % could try restrictions if you want, but I still recommend the
        %   piecewise linear fit method which is more robust
        %     pursuitIntervalEnd = pursuitIntervalStart + pursuitIntervalLength-1;
        %
        %     if pursuitIntervalStart
        %         % now fit a simple regression line to find the exact onset, which is the time of the intersection
        %         x = (pursuitIntervalStart:pursuitIntervalEnd)';
        %         x = [ones(size(x)) x];
        %         y = pursuitVel2D(pursuitIntervalStart:pursuitIntervalEnd);
        %         b = x\y;
        %         pursuit.onset = round((meanBaseline-b(1))/b(2));
        %         % restrict pursuit onset to be after 50ms after target onset
        %         if pursuit.onset < trial.log.targetOnset + ms2frames(50)
        %             pursuit.onset = trial.log.targetOnset + ms2frames(50);
        %         end
        %
        %         % for debugging
        %         pursuit.pursuitIntervalStart = pursuitIntervalStart;
        %         pursuit.pursuitIntervalEnd = pursuitIntervalEnd;
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
%     % for debugging
%     pursuit.baselineStart = baselineStart;
%     pursuit.baselineEnd = baselineEnd;
%     pursuit.meanBaseline = meanBaseline;
%     pursuit.sdBaseline = sdBaseline;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end