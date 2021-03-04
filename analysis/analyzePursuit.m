% FUNCTION to analyze open and closed loop pursuit when viewing moving
% stimuli; requires ms2frames.m
% history
% 07-2012       JE created analyzePursuit.m
% 2012-2018     JF added stuff to and edited analyzePursuit.m
% 14-07-2018    JF commented to make the script more accecable for future
%               VPOM students; for questions email jolande.fooken@rwth-aachen.de
% 13-Jan-2021   XW modified based on using findPursuitNew.m; xiuyunwu5@gmail.com
%
% input: trial --> structure containing relevant current trial information
%        pursuit --> structure containing pursuit onset
% output: pursuit --> structure containing relevant all open and closed
%                     loop pursuit measures

function [trial] = analyzePursuit(trial, pursuit)

%% set parameters
% defint the length of pursuit initiation phase
openLoopDuration = ms2frames(140);

pursuitOff = trial.log.targetOffset-ms2frames(100); % the end of the closed loop phase

%% Open loop
% first, check if any saccades overlap with the pursuit initiation phase
pursuit.openLoopEnd = pursuit.onset + openLoopDuration - 1;
sacIdx = find((trial.saccades.onsets>=pursuit.onset & trial.saccades.onsets<=pursuit.openLoopEnd) | ...
    (trial.saccades.offsets>=pursuit.onset & trial.saccades.offsets<=pursuit.openLoopEnd) | ...
    (pursuit.onset>=trial.saccades.onsets & pursuit.openLoopEnd<=trial.saccades.offsets));
if isempty(sacIdx)
    % if no saccades, mark as type 0 and go ahead to analyze pursuit initiation
    pursuit.onsetType = 0;
else % in general, just do not analyze pursuit initiation contaminated by saccades
    % check if the onset is within a saccade
    onsetSacIdx = find(trial.saccades.onsets<=pursuit.onset & trial.saccades.offsets>=pursuit.onset);
    if ~isempty(onsetSacIdx) % if yes, move the onset to the end of the saccade and mark it as type 1
        pursuit.onset = trial.saccades.offsets(onsetSacIdx)+1;
        pursuit.onsetType = 1;
    else % for other situations, mark as type 2
        pursuit.onsetType = 2;
    end
end

if pursuit.onsetType==0
    startFrame = nanmax([trial.log.targetOnset, pursuit.onset]); % if there is no pursuit onset we still want to analyze eye movement quaility
    endFrame = nanmax([trial.log.targetOnset + openLoopDuration, pursuit.openLoopEnd]);
    
    % first analyze initial pursuit in X
    velocityX = trial.eyeDX_filt(startFrame:endFrame);
    pursuit.initialMeanVelocityX = nanmean(velocityX); % why we used abs in the earlier version?... anyway
    peakIdx = find(abs(velocityX)==nanmax(abs(velocityX)));
    pursuit.initialPeakVelocityX = velocityX(peakIdx(1));
    
    accelerationX = trial.eyeDDX_filt(startFrame:endFrame);
    pursuit.initialMeanAccelerationX = nanmean(accelerationX);
    peakIdx = find(abs(accelerationX)==nanmax(abs(accelerationX)));
    pursuit.initialPeakAccelerationX = accelerationX(peakIdx(1));
    
    % then, fit a regression line to get an estimate of the overall
    % acceleration during the open loop phase
    xFit = (startFrame:endFrame)';
    xFit = [ones(size(xFit)) xFit];
    yFit = velocityX;
    bFit = xFit\yFit;
    pursuit.initialAccelerationFitX = bFit(2)*trial.log.eyeSampleRate;
    
    % next analyze initial pursuit in Y
    velocityY = trial.eyeDY_filt(startFrame:endFrame);
    pursuit.initialMeanVelocityY = nanmean(velocityY); % why we used abs in the earlier version?... anyway
    peakIdx = find(abs(velocityY)==nanmax(abs(velocityY)));
    pursuit.initialPeakVelocityY = velocityY(peakIdx(1));
    
    accelerationY = trial.eyeDDY_filt(startFrame:endFrame);
    pursuit.initialMeanAccelerationY = nanmean(accelerationY);
    peakIdx = find(abs(accelerationY)==nanmax(abs(accelerationY)));
    pursuit.initialPeakAccelerationY = accelerationY(peakIdx(1));
    
    yFit = velocityY;
    bFit = xFit\yFit;
    pursuit.initialAccelerationFitY = bFit(2)*trial.log.eyeSampleRate;
    
    % combine X and Y, 2D
    velocity2D = sqrt(trial.eyeDX_filt(startFrame:endFrame).^2 + trial.eyeDY_filt(startFrame:endFrame).^2);
    pursuit.initialMeanVelocity2D = nanmean(velocity2D); % why we used abs in the earlier version?... anyway
    peakIdx = find(abs(velocity2D)==nanmax(abs(velocity2D)));
    pursuit.initialPeakVelocity2D = velocity2D(peakIdx(1));
    
    acceleration2D = sqrt(trial.eyeDDX_filt(startFrame:endFrame).^2 + trial.eyeDDY_filt(startFrame:endFrame).^2);
    pursuit.initialMeanAcceleration2D = nanmean(acceleration2D);
    peakIdx = find(acceleration2D==nanmax(abs(acceleration2D)));
    pursuit.initialPeakAcceleration2D = acceleration2D(peakIdx);
    
    yFit = velocity2D;
    bFit = xFit\yFit;
    pursuit.initialAccelerationFit2D = bFit(2)*trial.log.eyeSampleRate;
else
    pursuit.initialMeanVelocityX = NaN;
    pursuit.initialPeakVelocityX = NaN;
    pursuit.initialMeanAccelerationX = NaN;
    pursuit.initialPeakAccelerationX = NaN;
    pursuit.initialMeanVelocityY = NaN;
    pursuit.initialPeakVelocityY = NaN;
    pursuit.initialMeanAccelerationY = NaN;
    pursuit.initialPeakAccelerationY = NaN;
    pursuit.initialMeanVelocity2D = NaN;
    pursuit.initialPeakVelocity2D = NaN;
    pursuit.initialMeanAcceleration2D = NaN;
    pursuit.initialPeakAcceleration2D = NaN;
    
    pursuit.initialAccelerationFitX = NaN;
    pursuit.initialAccelerationFitY = NaN;
    pursuit.initialAccelerationFit2D = NaN;
end

%% now analyze closed loop
% if there is no pursuit onset, use stimulus onset as onset
closedLoopFrames = [nanmax([trial.log.targetOnset + openLoopDuration, pursuit.openLoopEnd]); pursuitOff];
decisionFrames = [pursuitOff-150; pursuitOff-50];
if pursuitOff-100 < closedLoopFrames(1)
    decisionFrames = closedLoopFrames;
end

pursuit.closedLoopDuration = closedLoopFrames(2) - closedLoopFrames(1);
pursuit.decisionWindowDuration = decisionFrames(2) - decisionFrames(1);

if pursuit.closedLoopDuration < 0 % closed loop phase too short
    pursuit.gain = NaN;
    pursuit.gainDecision = NaN;
    pursuit.dirGain = NaN;
    pursuit.dirGainDecision = NaN;
    pursuit.dirError = NaN;
    pursuit.dirErrorDecision = NaN;
%     pursuit.velocityError = NaN;
else
    % calculate gain first,  the ratio of magnitude only
    speedXY_noSac = sqrt(trial.DX_noSac.^2 + trial.DY_noSac.^2);
    absoluteTargetVel = sqrt(trial.target.velocityX.^2 + trial.target.velocityY.^2);
    absoluteTargetVel(absoluteTargetVel < 0.05) = NaN;        
    pursuitGain = (speedXY_noSac)./absoluteTargetVel; 
    zScore = zscore(pursuitGain(~isnan(pursuitGain)));
    pursuitGain((zScore > 3 | zScore < -3)) = NaN;
    
    pursuit.targetVelCLP = nanmean(absoluteTargetVel(closedLoopFrames));
    pursuit.clpVel2D = nanmean(speedXY_noSac(closedLoopFrames));
    pursuit.decisionVel2D = nanmean(speedXY_noSac(decisionFrames));
    pursuit.gain = nanmean(pursuitGain(closedLoopFrames));
    pursuit.gain(pursuit.gain>2.5) = NaN;
    pursuit.gainDecision = nanmean(pursuitGain(decisionFrames));
    pursuit.gainDecision(pursuit.gainDecision>2.5) = NaN;
    
    % consider the direction as well, the ratio of (magnitude of projection on the target velocity direction)/(magnitude of target velocity) 
    eyeDot = [trial.DX_noSac'; trial.DY_noSac'];
    targetDot = [trial.target.velocityX'; trial.target.velocityY'];
    dotProduct = dot(eyeDot, targetDot)';
    dirGainAll = dotProduct./(trial.target.velocityX.^2 + trial.target.velocityY.^2);
    pursuit.dirGain = nanmean(dirGainAll(closedLoopFrames));
    pursuit.dirGainDecision = nanmean(dirGainAll(decisionFrames));
    
    % the mean direction error of eye velocity compared to target velocity
    targetDir = (atan2(trial.DY_noSac, trial.DX_noSac))/pi*180;
    diffDir = (atan2(trial.DY_noSac, trial.DX_noSac)-atan2(trial.target.velocityY, trial.target.velocityX))/pi*180;
    pursuit.targetMeanDirCLP = nanmean(targetDir(closedLoopFrames));
    pursuit.dirError = nanmean(diffDir(closedLoopFrames)); % in degs, how ccw the pursuit velocity rotated from the target velocity
    pursuit.dirErrorDecision = nanmean(diffDir(decisionFrames));
    
    % % calculate position error--not applicable to RDKs
    % horizontalError = trial.X_noSac(closedLoopFrames)-trial.target.X(closedLoopFrames);
    % verticalError = trial.target.Y(closedLoopFrames)-trial.Y_noSac(closedLoopFrames);
    % pursuit.positionError = nanmean(sqrt(horizontalError.^2+ verticalError.^2));
    
    %     % calculate velocity error
    %     pursuit.velocityError = nanmean(sqrt((trial.target.velocityX(closedLoopFrames) - trial.DX_noSac(closedLoopFrames)).^2 + ...
    %         (trial.target.velocityY(closedLoopFrames) - trial.DY_noSac(closedLoopFrames)).^2)); %auch 2D
end

trial.pursuit = pursuit;
end