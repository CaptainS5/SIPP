% FUNCTION to read in information about the target in each trial; align it
% with the eye data samples
%
% history
% 11-Mar-2021   created by XW for recovering the translating RDK in the micropursuit exp; xiuyunwu5@gmail.com

function [target] = readoutTarget(eyeData, Experiment, trial, currentSubjectPath, currentTrial, eventLog, rdkFrameLog)
% initialize target velocity and location vectors
target.velocityX = zeros(trial.log.trialEnd, 1);
target.velocityY = zeros(trial.log.trialEnd, 1);
target.internalVelX = zeros(trial.log.trialEnd, 1);
target.internalVelY = zeros(trial.log.trialEnd, 1);
target.averageVelX = zeros(trial.log.trialEnd, 1);
target.averageVelY = zeros(trial.log.trialEnd, 1);

target.posX = NaN(trial.log.trialEnd, 1);
target.posY = NaN(trial.log.trialEnd, 1);

% fill in the velocity values
% baseline
target.velocityX(trial.log.targetOnset:trial.log.targetOffset, 1) ...
    = Experiment.const.rdk.apertureSpeed*cos(trial.log.rdkApertureDirBefore/180*pi);
target.averageVelX(trial.log.targetOnset:trial.log.targetOffset, 1) ...
    = Experiment.const.rdk.apertureSpeed*cos(trial.log.rdkApertureDirBefore/180*pi);
% % during perturbation... --to be implemented...
% if trial.log.rdkApertureDirBefore==0
% target.velocityX(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
%     = Experiment.const.rdk.apertureSpeed*cos((trial.log.rdkApertureDirBefore+trial.log.rdkApertureDirPerturbation)/180*pi); % apertureDirBefore + apertureDirPerturbation
% target.velocityY(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
%     = Experiment.const.rdk.apertureSpeed*sin((trial.log.rdkApertureDirBefore+trial.log.rdkApertureDirPerturbation)/180*pi); % apertureDirBefore + apertureDirPerturbation
% target.averageVelX(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
%     = cos((trial.log.rdkApertureDirBefore+trial.log.rdkApertureDirPerturbation+trial.log.rdkInternalDirPerturbation)/180*pi); % aDirBefore+aDirPerturbation+internalDirPerturbation
% target.averageVelY(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
%     = sin((trial.log.rdkApertureDirBefore+trial.log.rdkApertureDirPerturbation+trial.log.rdkInternalDirPerturbation)/180*pi); % aDirBefore+aDirPerturbation+internalDirPerturbation
% elseif trial.log.rdkApertureDirBefore==180
% target.velocityX(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
%     = Experiment.const.rdk.apertureSpeed*cos((trial.log.rdkApertureDirBefore-trial.log.rdkApertureDirPerturbation)/180*pi); % apertureDirBefore + apertureDirPerturbation
% target.velocityY(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
%     = Experiment.const.rdk.apertureSpeed*sin((trial.log.rdkApertureDirBefore-trial.log.rdkApertureDirPerturbation)/180*pi); % apertureDirBefore + apertureDirPerturbation
% target.averageVelX(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
%     = cos((trial.log.rdkApertureDirBefore-trial.log.rdkApertureDirPerturbation-trial.log.rdkInternalDirPerturbation)/180*pi); % aDirBefore+aDirPerturbation+internalDirPerturbation
% target.averageVelY(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
%     = sin((trial.log.rdkApertureDirBefore-trial.log.rdkApertureDirPerturbation-trial.log.rdkInternalDirPerturbation)/180*pi); % aDirBefore+aDirPerturbation+internalDirPerturbation
% end
if trial.log.rdkCohPerturbation~=0
    target.internalVelX(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
        = cos(trial.log.rdkInternalDirPerturbation/180*pi); % relative internalDirPerturbation
    target.internalVelY(trial.log.perturbationOnset+1:trial.log.perturbationOffset, 1) ...
        = sin(trial.log.rdkInternalDirPerturbation/180*pi); % relative internalDirPerturbation
end

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
end