% FUNCTION to read in information about the target in each trial; align it
% with the eye data samples
%
% history
% 11-Mar-2021   created by XW for recovering the translating RDK in the micropursuit exp; xiuyunwu5@gmail.com

function [target] = readoutTarget(eyeData, Experiment, trial, currentSubjectPath, currentTrial, eventLog, rdkFrameLog)
% initialize target velocity and location vectors
target.velX = NaN(trial.log.trialEnd, 1);
target.velY = NaN(trial.log.trialEnd, 1);
target.internalVelX = NaN(trial.log.trialEnd, 1);
target.internalVelY = NaN(trial.log.trialEnd, 1);
target.averageVelX = NaN(trial.log.trialEnd, 1);
target.averageVelY = NaN(trial.log.trialEnd, 1);

target.posX = NaN(trial.log.trialEnd, 1);
target.posY = NaN(trial.log.trialEnd, 1);

% position values
% load the rdk center position for the current trial
trialIdxInData = eventLog.trialIdxInData(currentTrial, 1);
fileName = dir([currentSubjectPath, '\rdkSeeds\rdkseed_t', num2str(trialIdxInData, '%03d'), '*.mat']);
load([currentSubjectPath, '\rdkSeeds\', fileName.name]);

eyeFrames = eyeData.timeStamp(trial.log.targetOnset:trial.log.targetOffset);
rdkFrames = rdkFrameLog{currentTrial}.eyeLinkTimeStamp;
% to align the range of the frames (might be one frame difference...)
eyeFrames(1) = min(eyeFrames(1), rdkFrames(1));
rdkFrames(1) = eyeFrames(1);
eyeFrames(end) = max(eyeFrames(end), rdkFrames(end));
rdkFrames(end) = eyeFrames(end);
% displayed center positions, transferred into deg (screen center is (0,0))
rdkPosX = NaN(size(rdkFrames));
rdkPosY = NaN(size(rdkFrames));
for ii = 1:length(rdkAperturePos)
    rdkPosX(ii) = (rdkAperturePos{ii}(1)-(Experiment.screen.widthPX/2))*Experiment.screen.dpp;
    rdkPosY(ii) = ((Experiment.screen.heightPX/2)-rdkAperturePos{ii}(2))*Experiment.screen.dpp;
end
% interpolate target center positions
target.posX(trial.log.targetOnset:trial.log.targetOffset, 1) = interp1(rdkFrames, rdkPosX, eyeFrames);
target.posY(trial.log.targetOnset:trial.log.targetOffset, 1) = interp1(rdkFrames, rdkPosY, eyeFrames);

% velocity values
% calculate target (aperture) velocity based on position values
target.velX(trial.log.targetOnset:trial.log.targetOffset, 1) = [diff(target.posX(trial.log.targetOnset:trial.log.targetOffset, 1))*eyeData.sampleRate; NaN];
target.velY(trial.log.targetOnset:trial.log.targetOffset, 1) = [diff(target.posY(trial.log.targetOnset:trial.log.targetOffset, 1))*eyeData.sampleRate; NaN];

if trial.log.rdkCoh==0 % static pattern
    target.internalVelX(trial.log.targetOnset:trial.log.targetOffset, 1) ...
        = 0; % relative internal motion
    target.internalVelY(trial.log.targetOnset:trial.log.targetOffset, 1) ...
        = 0; % relative internal motion
else
    target.internalVelX(trial.log.targetOnset:trial.log.targetOffset, 1) ...
        = cos(trial.log.rdkInternalDir/180*pi)*trial.log.rdkInternalSpeed; % relative internal motion
    target.internalVelY(trial.log.targetOnset:trial.log.targetOffset, 1) ...
        = sin(trial.log.rdkInternalDir/180*pi)*trial.log.rdkInternalSpeed; % relative internal motion
end
target.averageVelX(trial.log.targetOnset:trial.log.targetOffset, 1) ...
    = target.velX(trial.log.targetOnset:trial.log.targetOffset, 1)+target.internalVelX(trial.log.targetOnset:trial.log.targetOffset, 1);
target.averageVelY(trial.log.targetOnset:trial.log.targetOffset, 1) ...
    = target.velY(trial.log.targetOnset:trial.log.targetOffset, 1)+target.internalVelY(trial.log.targetOnset:trial.log.targetOffset, 1);
end