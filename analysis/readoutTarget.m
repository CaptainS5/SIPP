% FUNCTION to read in information about the target in each trial; align it
% with the eye data samples
%
% history
% 13-Jan-2021   created by XW for the Gaussian distributed direction RDK in 
% the speed-accuracy task; xiuyunwu5@gmail.com

function [target] = readoutTarget(eyeData, dotSpeed, currentSubjectPath, currentTrial, eventLog)%, rdkFrameLog)
% based on the direction of each dot recorded, treat each as a velocity
% vector centered at zero, then calculate the mean target velocity vector

trialIdxInData = eventLog.trialIdxInData(currentTrial, 1);
fileName = dir([currentSubjectPath, '\rdkSeeds\rdkseed_t', num2str(trialIdxInData, '%03d'), '*.mat']);
load([currentSubjectPath, '\rdkSeeds\', fileName.name]);

% align the target frames to eye data samples
rdkDuration = eventLog.rdkOff(currentTrial)-eventLog.rdkOn(currentTrial);
frameN = floor(rdkDuration/12);
if frameN>50
    frameN = 50;
end

targetFrameIdx = [1:frameN];
targetFrameStamps = [eventLog.rdkOn(currentTrial)+round((targetFrameIdx-1)*1000/85)-1, eventLog.rdkOff(currentTrial)]; 
if targetFrameStamps(end)<targetFrameStamps(end-1)
    targetFrameStamps = [eventLog.rdkOn(currentTrial)+ceil((targetFrameIdx-1)*1000/85)-1, eventLog.rdkOff(currentTrial)]; 
end
% making them up now, assuming each target frame lasted 12 data frame...; will have the records in later versions

targetDirAll = NaN(size(rdkDirFrames, 2), length(eyeData.DX)); % each column is one frame, each row is one dot

for ii = 1:length(targetFrameIdx)-1 % fill in the target velocities into corresponding data sample
    % the last RDK frame has no motion 
    startF = targetFrameStamps(ii)-eyeData.timeStamp(1)+1;
    endF = targetFrameStamps(ii+1)-1-eyeData.timeStamp(1)+1;
    targetDirAll(:, startF:endF) = -repmat(rdkDirFrames(targetFrameIdx(ii), :)', 1, endF-startF+1); % CHANGE for later analysis; only flip for earlier versions
end
% then, calculate the mean target velocity vector in x and y
% previous pilot version--rdkFrameDir, in degs; 0/2pi equals to the horizontal right, rotates CW
% new versions--rdkFrameDir, in radians; direction already flipped, positive is up, negative is down
% to get the average, average the x coordinates (cosines) and y coordinates (sines) first
velXAll = dotSpeed.*cos(targetDirAll/180*pi);
velYAll = dotSpeed.*sin(targetDirAll/180*pi);

target.velocityX = nanmean(velXAll)';
target.velocityY = nanmean(velYAll)';
end