% FUNCTION to analyze saccade parameters

% history
% 07-2012       JE created analyzeSaccades.m
% 2012-2018     JF added stuff to and edited analyzeSaccades.m
% 13-07-2018    JF commented to make the script more accecable for future
%               VPOM students
% for questions email jolande.fooken@rwth-aachen.de
% Jan-08-2021    XW modified (again) the way to find peak velocity... it
% doesn't make sense to consider the saccade velocity baseline around
% target velocity, in my opinion...; and other edits for the translagint RDK exp,
% only using the combined saccade onsets/offsets. xiuyunwu5@gmail.com
%
% input: trial --> structure containing relevant current trial information
%        saccades --> output from findSaccades.m; contains on- & offsets
% output: trial --> structure containing relevant current trial information
%                   with saccades added
%         saccades --> edited saccade structure

function [trial, saccades] = analyzeSaccades(trial)
% define the window you want to analyze saccades in
% all saccade properties are within this window
startFrame = trial.log.perturbationOnset;
endFrame = trial.log.targetOffset;

% then find the proper onsets and offsets
idx = find(trial.saccades.onsets>=startFrame & trial.saccades.onsets<=endFrame);
trial.saccades.onsets_pursuit = trial.saccades.onsets(idx);
trial.saccades.offsets_pursuit = trial.saccades.offsets(idx);
trial.saccades.summary.num = length(trial.saccades.onsets_pursuit); % this is the number of all saccades during pursuit, including both olp and clp...
% however, if we simply exclude trials with saccades during opl, then the
% number would be clp only in other trials...

% xIdxL = find(trial.saccades.X_left.onsets>=startFrame & trial.saccades.X_left.onsets<=endFrame);
% xIdxR = find(trial.saccades.X_right.onsets>=startFrame & trial.saccades.X_right.onsets<=endFrame);
% trial.saccades.X_left.onsets_pursuit = trial.saccades.X_left.onsets(xIdxL);
% trial.saccades.X_left.offsets_pursuit = trial.saccades.X_left.offsets(xIdxL);
% trial.saccades.X_right.onsets_pursuit = trial.saccades.X_right.onsets(xIdxR);
% trial.saccades.X_right.offsets_pursuit = trial.saccades.X_right.offsets(xIdxR);

% calculate saccade amplitudes
% use the combined trial.saccades_onset, and calculate x/y/2D accordingly
if numel(trial.saccades.onsets_pursuit) == 0
    trial.saccades.amplitudes2D = NaN;
    trial.saccades.amplitudesX = NaN;
    trial.saccades.amplitudesY = NaN;
    trial.saccades.dirXY = NaN;
    trial.saccades.dirAngle = NaN;
    trial.saccades.eyePosOnset = NaN;
    trial.saccades.eyePosOffset = NaN;
    trial.saccades.targetPosOnset = NaN;
    trial.saccades.targetPosOffset = NaN; 
    
    trial.saccades.summary.meanAmp2D = NaN;
    trial.saccades.summary.meanAmpXLeft = NaN;
    trial.saccades.summary.meanAmpXRight = NaN;
    trial.saccades.summary.meanAmpYUp = NaN;
    trial.saccades.summary.meanAmpYDown = NaN;
    
    trial.saccades.summary.sumAmp2D = NaN;
    trial.saccades.summary.sumAmpXLeft = NaN;
    trial.saccades.summary.sumAmpXRight = NaN;
    trial.saccades.summary.sumAmpYUp = NaN;
    trial.saccades.summary.sumAmpYDown = NaN;
    
    trial.saccades.summary.numXLeft = NaN;
    trial.saccades.summary.numXRight = NaN;
    trial.saccades.summary.numYUp = NaN;
    trial.saccades.summary.numYDown = NaN;
else
    trial.saccades.amplitudes2D = sqrt((trial.eyeX_filt(trial.saccades.offsets_pursuit) - trial.eyeX_filt(trial.saccades.onsets_pursuit)).^2 ...
        + (trial.eyeY_filt(trial.saccades.offsets_pursuit) - trial.eyeY_filt(trial.saccades.onsets_pursuit)).^2);
    trial.saccades.amplitudesX = trial.eyeX_filt(trial.saccades.offsets_pursuit) - trial.eyeX_filt(trial.saccades.onsets_pursuit); % with direction
    trial.saccades.amplitudesY = trial.eyeY_filt(trial.saccades.offsets_pursuit) - trial.eyeY_filt(trial.saccades.onsets_pursuit); % with direction
    
    % the starting and ending location of each saccade and the target
    % location at that moment
    trial.saccades.eyePosOnset = [trial.eyeX_filt(trial.saccades.onsets_pursuit) trial.eyeY_filt(trial.saccades.onsets_pursuit)];
    trial.saccades.eyePosOffset = [trial.eyeX_filt(trial.saccades.offsets_pursuit) trial.eyeY_filt(trial.saccades.offsets_pursuit)];
    trial.saccades.targetPosOnset = [trial.target.posX(trial.saccades.onsets_pursuit) trial.target.posY(trial.saccades.onsets_pursuit)];
    trial.saccades.targetPosOffset = [trial.target.posX(trial.saccades.offsets_pursuit) trial.target.posY(trial.saccades.offsets_pursuit)];    
    
    % direction vector (unit length) of each saccade
    eyeVec = [trial.eyeX_filt(trial.saccades.offsets_pursuit) - trial.eyeX_filt(trial.saccades.onsets_pursuit), ...
        trial.eyeY_filt(trial.saccades.offsets_pursuit) - trial.eyeY_filt(trial.saccades.onsets_pursuit)]; % each row is one vector of the saccade
    normEye = vecnorm(eyeVec, 2, 2);
    eyeVec = eyeVec./repmat(normEye, 1, 2); % now each vector is of unit length
    trial.saccades.dirXY = eyeVec; % direction vectors, the first column is x, the second column is y
    trial.saccades.dirAngle = atan2(eyeVec(:,2), eyeVec(:, 1))/pi*180;
    
    trial.saccades.summary.meanAmp2D = nanmean(trial.saccades.amplitudes2D);
    trial.saccades.summary.meanAmpXLeft = -nanmean(trial.saccades.amplitudesX(trial.saccades.amplitudesX<=0));
    trial.saccades.summary.meanAmpXRight = nanmean(trial.saccades.amplitudesX(trial.saccades.amplitudesX>0));
    trial.saccades.summary.meanAmpYUp = nanmean(trial.saccades.amplitudesY(trial.saccades.amplitudesY>0));
    trial.saccades.summary.meanAmpYDown = -nanmean(trial.saccades.amplitudesY(trial.saccades.amplitudesY<=0));
    
    trial.saccades.summary.sumAmp2D = nansum(trial.saccades.amplitudes2D);
    trial.saccades.summary.sumAmpXLeft = -nansum(trial.saccades.amplitudesX(trial.saccades.amplitudesX<=0));
    trial.saccades.summary.sumAmpXRight = nansum(trial.saccades.amplitudesX(trial.saccades.amplitudesX>0));
    trial.saccades.summary.sumAmpYUp = nansum(trial.saccades.amplitudesY(trial.saccades.amplitudesY>0));
    trial.saccades.summary.sumAmpYDown = -nansum(trial.saccades.amplitudesY(trial.saccades.amplitudesY<=0));
    
    trial.saccades.summary.numXLeft = length(find((trial.saccades.amplitudesX<=0)));
    trial.saccades.summary.numXRight = length(find((trial.saccades.amplitudesX>0)));
    trial.saccades.summary.numYUp = length(find((trial.saccades.amplitudesY>0)));
    trial.saccades.summary.numYDown = length(find((trial.saccades.amplitudesY<=0)));
end

% xSacL = length(trial.saccades.X_left.onsets_pursuit);
% xSacR = length(trial.saccades.X_right.onsets_pursuit);
% if ~isempty(xSacL)
%     trial.saccades.X_left.amplitudes = abs(trial.eyeX_filt(trial.saccades.X_left.offsets_pursuit) - trial.eyeX_filt(trial.saccades.X_left.onsets_pursuit));
%     % trial.saccades.X_left.amplitudes = sqrt((trial.eyeX_filt(trial.saccades.X_left.offsets_pursuit) - trial.eyeX_filt(trial.saccades.X_left.onsets_pursuit)).^2 ...
%     %         + (trial.eyeY_filt(trial.saccades.X_left.offsets_pursuit) - trial.eyeY_filt(trial.saccades.X_left.onsets_pursuit)).^2);
% else
%     trial.saccades.X_left.amplitudes = NaN;
% end
% if ~isempty(xSacR)
%     trial.saccades.X_right.amplitudes = abs(trial.eyeX_filt(trial.saccades.X_right.offsets_pursuit) - trial.eyeX_filt(trial.saccades.X_right.onsets_pursuit));
% else
%     trial.saccades.X_right.amplitudes = NaN;
% end

% caluclate mean and max amplitude, mean duration, total number, &
% cumulative saccade amplitude (saccadic sum)
% if isempty(trial.saccades.onsets_pursuit)
%     trial.saccades.meanAmplitude = NaN;
%     trial.saccades.maxAmplitude = NaN;   
%     trial.saccades.X.meanDuration = NaN;
%     trial.saccades.Y.meanDuration = NaN;
%     trial.saccades.meanDuration = NaN;
%     trial.saccades.number = NaN;
%     trial.saccades.sacSum = NaN;
% else
%     trial.saccades.meanAmplitude = nanmean(trial.saccades.amplitudes);
%     trial.saccades.maxAmplitude = max(trial.saccades.amplitudes);
%     trial.saccades.X.meanDuration = mean(trial.saccades.X.offsets_pursuit - trial.saccades.X.onsets_pursuit);
%     trial.saccades.Y.meanDuration = mean(trial.saccades.Y.offsets_pursuit - trial.saccades.Y.onsets_pursuit);
%     trial.saccades.meanDuration = nanmean(sqrt(trial.saccades.X.meanDuration.^2 + ...
%                                                trial.saccades.Y.meanDuration.^2));
%     trial.saccades.number = length(trial.saccades.onsets_pursuit);
%     trial.saccades.sacSum = sum(trial.saccades.amplitudes);
% end

% % calculate mean and peak velocity for each saccade; then find average
% trial.saccades.X.peakVelocity = NaN;
% trial.saccades.Y.peakVelocity = NaN;
% trial.saccades.X.meanVelocity = NaN;
% trial.saccades.Y.meanVelocity = NaN;
% saccadesXXpeakVelocity = NaN(length(trial.saccades.X.onsets_pursuit),1);
% saccadesXYpeakVelocity = NaN(length(trial.saccades.X.onsets_pursuit),1);
% saccadesXXmeanVelocity = NaN(length(trial.saccades.X.onsets_pursuit),1);
% saccadesXYmeanVelocity = NaN(length(trial.saccades.X.onsets_pursuit),1);
% if ~isempty(trial.saccades.X.onsets_pursuit)
%     for i = 1:length(trial.saccades.X.onsets_pursuit)
%         saccadesXXpeakVelocity(i) = nanmax(abs(trial.eyeDX_filt(trial.saccades.X.onsets_pursuit(i):trial.saccades.X.offsets_pursuit(i))));
%         saccadesXYpeakVelocity(i) = nanmax(abs(trial.eyeDY_filt(trial.saccades.X.onsets_pursuit(i):trial.saccades.X.offsets_pursuit(i))));
%         saccadesXXmeanVelocity(i) = nanmean(abs(trial.eyeDX_filt(trial.saccades.X.onsets_pursuit(i):trial.saccades.X.offsets_pursuit(i))));
%         saccadesXYmeanVelocity(i) = nanmean(abs(trial.eyeDY_filt(trial.saccades.X.onsets_pursuit(i):trial.saccades.X.offsets_pursuit(i))));
%     end
% end
% saccadesYYpeakVelocity = NaN(length(trial.saccades.Y.onsets_pursuit),1);
% saccadesYXpeakVelocity = NaN(length(trial.saccades.Y.onsets_pursuit),1);
% saccadesYYmeanVelocity = NaN(length(trial.saccades.Y.onsets_pursuit),1);
% saccadesYXmeanVelocity = NaN(length(trial.saccades.Y.onsets_pursuit),1);
% if ~isempty(trial.saccades.Y.onsets_pursuit)
%     for i = 1:length(trial.saccades.Y.onsets_pursuit)
%         saccadesYYpeakVelocity(i) = nanmax(abs(trial.eyeDY_filt(trial.saccades.Y.onsets_pursuit(i):trial.saccades.Y.offsets_pursuit(i))));
%         saccadesYXpeakVelocity(i) = nanmax(abs(trial.eyeDX_filt(trial.saccades.Y.onsets_pursuit(i):trial.saccades.Y.offsets_pursuit(i))));
%         saccadesYYmeanVelocity(i) = nanmean(abs(trial.eyeDY_filt(trial.saccades.Y.onsets_pursuit(i):trial.saccades.Y.offsets_pursuit(i))));
%         saccadesYXmeanVelocity(i) = nanmean(abs(trial.eyeDX_filt(trial.saccades.Y.onsets_pursuit(i):trial.saccades.Y.offsets_pursuit(i))));
%     end
% end
% if ~isempty(trial.saccades.X.onsets_pursuit) || ~isempty(trial.saccades.Y.onsets_pursuit)
%     trial.saccades.X.peakVelocity = nanmax([saccadesXXpeakVelocity; saccadesYXpeakVelocity]);
%     trial.saccades.Y.peakVelocity = nanmax([saccadesXYpeakVelocity; saccadesYYpeakVelocity]);
%     trial.saccades.X.meanVelocity = nanmean([saccadesXXmeanVelocity; saccadesYXmeanVelocity]);
%     trial.saccades.Y.meanVelocity = nanmean([saccadesXYmeanVelocity; saccadesYYmeanVelocity]);
% end
% trial.saccades.peakVelocity = nanmean(sqrt(trial.saccades.X.peakVelocity.^2 + trial.saccades.Y.peakVelocity.^2));
% trial.saccades.meanVelocity = nanmean(sqrt(trial.saccades.X.meanVelocity.^2 + trial.saccades.Y.meanVelocity.^2));

end
