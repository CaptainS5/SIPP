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
% for the perturbation study, we mainly want to compare the pursuit velocity &
% direction after perturbation

%% set parameters
% define the analysis window
windowStart = trial.log.perturbationOnset + ms2frames(100); % 100ms after perturbation on
windowEnd = trial.log.targetOffset; % perturbation end

pursuitOnset = pursuit.summary.onset; 
openLoopDuration = ms2frames(140);
openLoopEnd = pursuitOnset+openLoopDuration-1;

% label the perturbation window situation
% 0: during initiation-before pursuit onset
% 1: during initiaiton-after pursuit onset
% -1: during initiation-after openloop... should definitely exclude these
% trials
% 2: during steady-state
if windowStart < pursuitOnset
    pursuit.summary.perturbationLabel = 0;
elseif windowStart>=pursuitOnset && windowStart<openLoopEnd
    pursuit.summary.perturbationLabel = 1;
elseif windowStart>=openLoopEnd && windowStart<trial.log.targetOnset+ms2frames(500)
    pursuit.summary.perturbationLabel = -1;
else
    pursuit.summary.perturbationLabel = 2;
end

% normalize all velocity vectors (at each time point)
eyeDir = [trial.DX_noSac, trial.DY_noSac]; % each row is one velocity vector
normEye = vecnorm(eyeDir, 2, 2);
pursuit.eyeDir = eyeDir./repmat(normEye, 1, 2); % now each vector is of unit length

%% Open loop type
% first, check if any saccades overlap with the pursuit initiation phase
pursuit.summary.openLoopEnd = openLoopEnd;
sacIdx = find((trial.saccades.onsets>=pursuit.summary.onset & trial.saccades.onsets<=pursuit.summary.openLoopEnd) | ...
    (trial.saccades.offsets>=pursuit.summary.onset & trial.saccades.offsets<=pursuit.summary.openLoopEnd) | ...
    (pursuit.summary.onset>=trial.saccades.onsets & pursuit.summary.openLoopEnd<=trial.saccades.offsets));
if isempty(sacIdx) && trial.log.eyeType==1 % pursuit onset without saccades in the pursuit condition
    % if no saccades, mark as type 0 and go ahead to analyze pursuit initiation
    pursuit.summary.onsetType = 0;
else % in general, just do not analyze pursuit initiation contaminated by saccades
    % check if the onset is within a saccade
    onsetSacIdx = find(trial.saccades.onsets<=pursuit.summary.onset & trial.saccades.offsets>=pursuit.summary.onset);
    if ~isempty(onsetSacIdx) % if yes, mark it as type 1
%         pursuit.summary.onset = trial.saccades.offsets(onsetSacIdx)+1; % move the onset to the end of the saccade 
        pursuit.summary.onsetType = 1;
    else % for other situations, mark as type 2
        pursuit.summary.onsetType = 2;
    end
end

%% Perturbation phase analysis
% for all perturbation situations, analyze mean + peak velocity,
% acceleration(?), and direction
% if windowEnd-windowStart <= 50 % shouldn't happen... 
%     pursuit.summary.gainXexternal = NaN;
%     pursuit.summary.gainYexternal = NaN;
%     pursuit.summary.gainYaverage = NaN;
%     pursuit.summary.gain2Dexternal = NaN;
%     pursuit.summary.gain2Daverage = NaN;
%     pursuit.summary.dirGainExternal = NaN;
%     pursuit.summary.dirPertubX = NaN;
%     pursuit.summary.dirPertubY = NaN;
% %     pursuit.summary.dirClpEarlyX = NaN;
% %     pursuit.summary.dirClpEarlyY = NaN;
% %     pursuit.summary.dirClpLateX = NaN;
% %     pursuit.summary.dirClpLateY = NaN;
%     pursuit.summary.dirRange = NaN;
%     pursuit.summary.dirError = NaN;
%     pursuit.summary.disCenterMean = NaN;
% %     pursuit.summary.disCenterEarly = NaN;
% %     pursuit.summary.disCenterMeanLate = NaN;
%     pursuit.summary.travelPertubDis = NaN;
%     pursuit.summary.targetPertubDis = NaN;
% else
    % calculate gain first,  the ratio of magnitude only
    % horizontal pursuit gain, if only account for aperture motion
    targetVelX = trial.target.velX;
    targetVelX(abs(targetVelX) < 0.05) = NaN;        
    pursuitGainX = trial.DX_noSac./targetVelX; 
    zScore = zscore(pursuitGainX(~isnan(pursuitGainX)));
    pursuitGainX((zScore > 3 | zScore < -3)) = NaN;
    pursuit.summary.gainXexternal = nanmean(pursuitGainX(windowStart:windowEnd));
    pursuit.summary.gainXexternal(pursuit.summary.gainXexternal>2.5) = NaN;
%     % horizontal gain in relation to average target motion
%     targetVelX = averageVecX;
%     targetVelX(abs(targetVelX) < 0.05) = NaN;        
%     pursuitGainX = trial.DX_noSac./targetVelX; 
%     zScore = zscore(pursuitGainX(~isnan(pursuitGainX)));
%     pursuitGainX((zScore > 3 | zScore < -3)) = NaN;
%     pursuit.summary.gainXaverage = nanmean(pursuitGainX(closedLoopFrames));
%     pursuit.summary.gainXaverage(pursuit.summary.gainXaverage>2.5) = NaN;
    
    % vertical pursuit gain, if only account for aperture motion
    targetVelY = trial.target.velY;
    targetVelY(abs(targetVelY) < 0.05) = NaN;        
    pursuitGainY = trial.DY_noSac./targetVelY; 
    zScore = zscore(pursuitGainY(~isnan(pursuitGainY)));
    pursuitGainY((zScore > 3 | zScore < -3)) = NaN;
    pursuit.summary.gainYexternal = nanmean(pursuitGainY(windowStart:windowEnd));
    pursuit.summary.gainYexternal(pursuit.summary.gainYexternal>2.5) = NaN;
    % vertical pursuit gain, in relation to average target motion (aperture+internal)
    targetVelY = trial.target.averageVelY;
    targetVelY(abs(targetVelY) < 0.05) = NaN;        
    pursuitGainY = trial.DY_noSac./targetVelY; 
    zScore = zscore(pursuitGainY(~isnan(pursuitGainY)));
    pursuitGainY((zScore > 3 | zScore < -3)) = NaN;
    pursuit.summary.gainYaverage = nanmean(pursuitGainY(windowStart:windowEnd));
    pursuit.summary.gainYaverage(pursuit.summary.gainYaverage>2.5) = NaN;
     
    % pursuit gain 2D, if only account for external motion
    speedXY_noSac = sqrt(trial.DX_noSac.^2 + trial.DY_noSac.^2); 
    % calculate the averaged vector velocity in 2D
    absoluteTargetVel = sqrt(trial.target.velX.^2 + trial.target.velY.^2);
    absoluteTargetVel(absoluteTargetVel < 0.05) = NaN;        
    pursuitGain = (speedXY_noSac)./absoluteTargetVel; 
    zScore = zscore(pursuitGain(~isnan(pursuitGain)));
    pursuitGain((zScore > 3 | zScore < -3)) = NaN;
    pursuit.summary.gain2Dexternal = nanmean(pursuitGain(windowStart:windowEnd));
    pursuit.summary.gain2Dexternal(pursuit.summary.gain2Dexternal>2.5) = NaN;
    % if account for external+internal vector average
    absoluteTargetVel = sqrt(trial.target.averageVelX.^2 + trial.target.averageVelY.^2);
    absoluteTargetVel(absoluteTargetVel < 0.05) = NaN;
    pursuitGain = (speedXY_noSac)./absoluteTargetVel;
    zScore = zscore(pursuitGain(~isnan(pursuitGain)));
    pursuitGain((zScore > 3 | zScore < -3)) = NaN;
    pursuit.summary.gain2Daverage = nanmean(pursuitGain(windowStart:windowEnd));
    pursuit.summary.gain2Daverage(pursuit.summary.gain2Daverage>2.5) = NaN;
    
    % consider the direction as well, the ratio of (magnitude of projection on the target velocity direction)/(magnitude of target velocity) 
    eyeDot = [trial.DX_noSac'; trial.DY_noSac'];
    targetDot = [trial.target.velX'; trial.target.velY'];
    dotProduct = dot(eyeDot, targetDot)';
    dirGainAll = dotProduct./(trial.target.velX.^2 + trial.target.velY.^2);
    pursuit.summary.dirGainExternal = nanmean(dirGainAll(windowStart:windowEnd));
    
    % the mean direction of eye velocity (smooth phase)
    % calculate the average velocity vector; just record the vector coordinate values;  
    % if need to calculate for example the average direction of a participant, still need to use these
    % velocity vectors
    xTemp = nansum(pursuit.eyeDir(windowStart:windowEnd, 1));
    yTemp = nansum(pursuit.eyeDir(windowStart:windowEnd, 2));
    % again, normalize to unit length vector for later calculations
    pursuit.summary.dirPerturbX = xTemp/sqrt(xTemp^2+yTemp^2);
    pursuit.summary.dirPerturbY = yTemp/sqrt(xTemp^2+yTemp^2);
    
%     % separate into early and late phases (half & half), then calculate mean
%     % direction again
%     % the first half, early phase
%     xTemp = nansum(eyeDir(closedLoopFrames(1):closedLoopFrames(round(pursuit.summary.closedLoopDuration/2)), 1));
%     yTemp = nansum(eyeDir(closedLoopFrames(1):closedLoopFrames(round(pursuit.summary.closedLoopDuration/2)), 2));
%     % again, normalize to unit length vector for later calculations
%     pursuit.summary.dirClpEarlyX = xTemp/sqrt(xTemp^2+yTemp^2);
%     pursuit.summary.dirClpEarlyY = yTemp/sqrt(xTemp^2+yTemp^2);
%     
%     % the second half, late phase
%     xTemp = nansum(eyeDir(closedLoopFrames(round(pursuit.summary.closedLoopDuration/2)+1):closedLoopFrames(end), 1));
%     yTemp = nansum(eyeDir(closedLoopFrames(round(pursuit.summary.closedLoopDuration/2)+1):closedLoopFrames(end), 2));
%     % again, normalize to unit length vector for later calculations
%     pursuit.summary.dirClpLateX = xTemp/sqrt(xTemp^2+yTemp^2);
%     pursuit.summary.dirClpLateY = yTemp/sqrt(xTemp^2+yTemp^2);
%     
%     dirClpEarly = atan2(pursuit.summary.dirClpEarlyY, pursuit.summary.dirClpEarlyX)/pi*180;
%     dirClpLate = atan2(pursuit.summary.dirClpLateY, pursuit.summary.dirClpLateX)/pi*180;
%     pursuit.summary.dirClpChange = dirClpLate-dirClpEarly;
    
    dir = atan2(eyeDir(windowStart:windowEnd, 2), eyeDir(windowStart:windowEnd, 1))/pi*180;
    pursuit.summary.dirRange = nanmax(dir)-nanmin(dir); % the difference between the min and max direction; to identify trials with dramatic direction changes
    
%     %%%%%%%%%%%%%%%% old method using atan2... might be tricky
%     %     dir = atan2(trial.DY_noSac, trial.DX_noSac)/pi*180;
%     %     % need to make the angle "averagable"...
%     %     if trial.log.rdkApertureDir==180 % leftward, need to average ~180 degs
%     %         dir(dir<0) = 360+dir(dir<0);
%     %     end
%     %     pursuit.summary.dirClp = nanmean(dir(closedLoopFrames)); % in degs, horizontal right is 0, ccw is positive
%     %%%%%%%%%%%%%%%%
%     
    % the mean direction error of eye velocity compared to target
    % velocity (external only)
    % normalize target velocity vectors
    targetDir = [trial.target.velX, trial.target.velY];
    normTarget = vecnorm(targetDir, 2, 2);
    targetDir = targetDir./repmat(normTarget, 1, 2);
    % now calculate the angle between each pairs of the vectors
    dotProduct = dot(eyeDir, targetDir, 2);
    dirErrors = acos(dotProduct)/pi*180; % both are already unit length vectors, norms are all 1
    %%%%%%%%%%%%%%%% old method using atan2... might be tricky
    %     targetDir = (atan2(vecY, vecX))/pi*180;
    %     if trial.log.rdkApertureDir==180
    %         targetDir(targetDir<0) = 360+targetDir(targetDir<0);
    %     end
    %     dirErrors = abs(dir-targetDir); % errors, absolute value
    %     pursuit.summary.targetMeanDirCLP = nanmean(targetDir(closedLoopFrames));
    %%%%%%%%%%%%%%%%
    pursuit.summary.dirError = nanmean(dirErrors(windowStart:windowEnd)); % in degs
    
    % calculate the average distance of eye position to aperture center
    disCenter = sqrt((trial.eyeX_filt-trial.target.posX).^2 + ...
        (trial.eyeY_filt-trial.target.posY).^2);
    pursuit.summary.disCenterMean = nanmean(disCenter(windowStart:windowEnd));
%     pursuit.summary.disCenterMeanEarly = nanmean(disCenter(closedLoopFrames(1):closedLoopFrames(round(pursuit.summary.closedLoopDuration/2))));
%     pursuit.summary.disCenterMeanLate = nanmean(disCenter(closedLoopFrames(round(pursuit.summary.closedLoopDuration/2)+1):closedLoopFrames(end)));
    
    % mainly used for exclusion criteria, the total 2D distance traveled by
    % the eye during the whole clp duration
    pursuit.summary.travelClpDis = sqrt((trial.eyeX_filt(windowEnd)-trial.eyeX_filt(windowStart)).^2 ...
        + (trial.eyeY_filt(windowEnd)-trial.eyeY_filt(windowStart)).^2);
    pursuit.summary.targetClpDis = sqrt((trial.target.posX(windowEnd)-trial.target.posX(windowStart)).^2 ...
        + (trial.target.posY(windowEnd)-trial.target.posY(windowStart)).^2);
    
%     % velocity variation
%     pursuit.summary.velCovX = nanstd(trial.DX_noSac(closedLoopFrames))/nanmean(trial.DX_noSac(closedLoopFrames));
%     pursuit.summary.velCovY = nanstd(trial.DY_noSac(closedLoopFrames))/nanmean(trial.DY_noSac(closedLoopFrames));
%     pursuit.summary.velCov2D = nanstd(speedXY_noSac(closedLoopFrames))/nanmean(speedXY_noSac(closedLoopFrames));
% end

trial.pursuit = pursuit;
end