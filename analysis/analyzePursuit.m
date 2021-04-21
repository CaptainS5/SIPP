% FUNCTION to analyze open and closed loop pursuit when viewing moving
% stimuli; requires ms2frames.m
% history
% 07-2012       JE created analyzePursuit.m
% 2012-2018     JF added stuff to and edited analyzePursuit.m
% 14-07-2018    JF commented to make the script more accecable for future
%               VPOM students; for questions email jolande.fooken@rwth-aachen.de
% 11-Mar-2021   XW modified based on using findPursuitNew.m; xiuyunwu5@gmail.com
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
if isempty(sacIdx) && trial.log.eyeType==1 % pursuit onset without saccades in the pursuit condition
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
    
    pursuit.latency = startFrame-trial.log.targetOnset;
    
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
closedLoopFrames = [nanmax([trial.log.targetOnset + openLoopDuration, pursuit.openLoopEnd]):pursuitOff];

pursuit.closedLoopDuration = closedLoopFrames(2) - closedLoopFrames(1);

% if pursuit.closedLoopDuration < 0 % closed loop phase too short
%     pursuit.gainXexternal = NaN;
%     pursuit.gainYinternal = NaN;
%     pursuit.gainXaverage = NaN;
%     pursuit.gainYaverage = NaN;
%     pursuit.gain2D = NaN;
% else
%     % the averaged target vector (external+internal)
%     externalVecX = trial.target.velocityX;
%     externalVecY = 0;
%     if trial.log.rdkCoh==0
%         internalVecX = 0;
%         internalVecY = 0;
%     else
%         if trial.log.rdkApertureDir==0 % rightward
%             internalVecX = trial.log.rdkInternalSpeed.*cos(trial.log.rdkInternalDir/180*pi);
%         else
%             internalVecX = trial.log.rdkInternalSpeed.*cos((180-trial.log.rdkInternalDir)/180*pi);
%         end
%         internalVecY = trial.log.rdkInternalSpeed.*sin(trial.log.rdkInternalDir/180*pi);
%     end
%     averageVecX = externalVecX+internalVecX;
%     averageVecY = externalVecY+internalVecY;
%     % make them into vectors the same length as eye velocity
%     averageVecX = repmat(averageVecX, size(trial.DX_noSac));
%     averageVecY = repmat(averageVecY, size(trial.DX_noSac));
%     
%     % calculate gain first,  the ratio of magnitude only
%     % horizontal pursuit gain, if only account for external aperture motion
%     targetVelX = trial.target.velocityX;
%     targetVelX(abs(targetVelX) < 0.05) = NaN;        
%     pursuitGainX = trial.DX_noSac./targetVelX; 
%     zScore = zscore(pursuitGainX(~isnan(pursuitGainX)));
%     pursuitGainX((zScore > 3 | zScore < -3)) = NaN;
%     pursuit.gainXexternal = nanmean(pursuitGainX(closedLoopFrames));
%     pursuit.gainXexternal(pursuit.gainXexternal>2.5) = NaN;
%     % horizontal gain in relation to average target motion
%     targetVelX = averageVecX;
%     targetVelX(abs(targetVelX) < 0.05) = NaN;        
%     pursuitGainX = trial.DX_noSac./targetVelX; 
%     zScore = zscore(pursuitGainX(~isnan(pursuitGainX)));
%     pursuitGainX((zScore > 3 | zScore < -3)) = NaN;
%     pursuit.gainXaverage = nanmean(pursuitGainX(closedLoopFrames));
%     pursuit.gainXaverage(pursuit.gainXaverage>2.5) = NaN;
%     
%     % vertical pursuit gain, if only account for relative internal dot motion
%     targetVelY = trial.log.rdkInternalSpeed.*sin(trial.log.rdkInternalDir/180*pi);
%     targetVelY(abs(targetVelY) < 0.05) = NaN;        
%     pursuitGainY = trial.DY_noSac./targetVelY; 
%     zScore = zscore(pursuitGainY(~isnan(pursuitGainY)));
%     pursuitGainY((zScore > 3 | zScore < -3)) = NaN;
%     pursuit.gainYinternal = nanmean(pursuitGainY(closedLoopFrames));
%     pursuit.gainYinternal(pursuit.gainYinternal>2.5) = NaN;
%     % vertical pursuit gain, in relation to average target motion
%     targetVelY = averageVecY;
%     targetVelY(abs(targetVelY) < 0.05) = NaN;        
%     pursuitGainY = trial.DY_noSac./targetVelY; 
%     zScore = zscore(pursuitGainY(~isnan(pursuitGainY)));
%     pursuitGainY((zScore > 3 | zScore < -3)) = NaN;
%     pursuit.gainYaverage = nanmean(pursuitGainY(closedLoopFrames));
%     pursuit.gainYaverage(pursuit.gainYaverage>2.5) = NaN;
%     
%     % pursuit gain, if account for external+internal vector average
%     speedXY_noSac = sqrt(trial.DX_noSac.^2 + trial.DY_noSac.^2); 
%     % calculate the averaged vector velocity in 2D
%     absoluteTargetVel = sqrt(averageVecX.^2 + averageVecY.^2);
%     absoluteTargetVel(absoluteTargetVel < 0.05) = NaN;        
%     pursuitGain = (speedXY_noSac)./absoluteTargetVel; 
%     zScore = zscore(pursuitGain(~isnan(pursuitGain)));
%     pursuitGain((zScore > 3 | zScore < -3)) = NaN;
%     pursuit.gain2D = nanmean(pursuitGain(closedLoopFrames));
%     pursuit.gain2D(pursuit.gain2D>2.5) = NaN;
%     
%     % consider the direction as well, the ratio of (magnitude of projection on the target velocity direction)/(magnitude of target velocity) 
%     eyeDot = [trial.DX_noSac'; trial.DY_noSac'];
%     targetDot = [averageVecX'; averageVecY'];
%     dotProduct = dot(eyeDot, targetDot)';
%     dirGainAll = dotProduct./(averageVecX.^2 + averageVecY.^2);
%     pursuit.dirGain = nanmean(dirGainAll(closedLoopFrames));
%     
%     % the mean direction of eye velocity
%     % first, normalize all velocity vectors (at each time point)
%     eyeDir = [trial.DX_noSac, trial.DY_noSac]; % each row is one velocity vector
%     normEye = vecnorm(eyeDir, 2, 2);
%     eyeDir = eyeDir./repmat(normEye, 1, 2); % now each vector is of unit length
%     % then, calculate the average velocity vector; just record the vector coordinate values;  
%     % if need to calculate for example the average direction of a participant, still need to use these
%     % velocity vectors
%     xTemp = nansum(eyeDir(closedLoopFrames, 1));
%     yTemp = nansum(eyeDir(closedLoopFrames, 2));
%     % again, normalize to unit length vector for later calculations
%     pursuit.dirClpX = xTemp/sqrt(xTemp^2+yTemp^2);
%     pursuit.dirClpY = yTemp/sqrt(xTemp^2+yTemp^2);
%     
%     %%%%%%%%%%%%%%%% old method using atan2... might be tricky
%     %     dir = atan2(trial.DY_noSac, trial.DX_noSac)/pi*180;
%     %     % need to make the angle "averagable"...
%     %     if trial.log.rdkApertureDir==180 % leftward, need to average ~180 degs
%     %         dir(dir<0) = 360+dir(dir<0);
%     %     end
%     %     pursuit.dirClp = nanmean(dir(closedLoopFrames)); % in degs, horizontal right is 0, ccw is positive
%     %%%%%%%%%%%%%%%%
%     
%     % the mean direction error of eye velocity compared to averaged target
%     % velocity (external+internal)
%     % normalize target velocity vectors
%     targetDir = [averageVecX, averageVecY];
%     normTarget = vecnorm(targetDir, 2, 2);
%     targetDir = targetDir./repmat(normTarget, 1, 2);
%     % now calculate the angle between each pairs of the vectors
%     dotProduct = dot(eyeDir, targetDir, 2);
%     dirErrors = acos(dotProduct)/pi*180; % both are already unit length vectors, norms are all 1
%     %%%%%%%%%%%%%%%% old method using atan2... might be tricky
%     %     targetDir = (atan2(vecY, vecX))/pi*180;
%     %     if trial.log.rdkApertureDir==180
%     %         targetDir(targetDir<0) = 360+targetDir(targetDir<0);
%     %     end
%     %     dirErrors = abs(dir-targetDir); % errors, absolute value
%     %     pursuit.targetMeanDirCLP = nanmean(targetDir(closedLoopFrames));
%     %%%%%%%%%%%%%%%%
%     pursuit.dirError = nanmean(dirErrors(closedLoopFrames)); % in degs
%     
%     % velocity variation
%     pursuit.velCovX = nanstd(trial.DX_noSac(closedLoopFrames))/nanmean(trial.DX_noSac(closedLoopFrames));
%     pursuit.velCovY = nanstd(trial.DY_noSac(closedLoopFrames))/nanmean(trial.DY_noSac(closedLoopFrames));
%     pursuit.velCov2D = nanstd(speedXY_noSac(closedLoopFrames))/nanmean(speedXY_noSac(closedLoopFrames));
%     
%     %     % calculate velocity error
%     %     pursuit.velocityError = nanmean(sqrt((trial.target.velocityX(closedLoopFrames) - trial.DX_noSac(closedLoopFrames)).^2 + ...
%     %         (trial.target.velocityY(closedLoopFrames) - trial.DY_noSac(closedLoopFrames)).^2)); %auch 2D
% end

trial.pursuit = pursuit;
end