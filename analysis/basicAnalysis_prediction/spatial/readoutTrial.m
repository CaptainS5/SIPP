% FUNCTION to set up data structure for reading out and later saving all
% relevant experimental info, data, and parameters
% history
% 07-2012       JE created readoutTrial.m
% 05-2014       JF edited readoutTrial.m
% 13-07-2018    JF commented to make the script more accecable for future 
%               VPOM students
% for questions email jolande.fooken@rwth-aachen.de
% 
% input: eyeData --> structure containing filtered eye movements
%        currentSubject --> selected subject using selectSubject.m
%        parameters --> all other relevant experimental information
%        currentTrial --> current trial
% output: trial --> a structure containing all relevant information for
%                   each trial

function [trial] = readoutTrial(eyeData, currentSubject, analysisPath, Experiment, currentTrial, eventLog)
% get eyeData for this trial
trial.eyeX_filt = eyeData.X_filt;
trial.eyeY_filt = eyeData.Y_filt;

trial.eyeDX_filt = eyeData.DX_filt;
trial.eyeDY_filt = eyeData.DY_filt;

trial.eyeDDX_filt = eyeData.DDX_filt;
trial.eyeDDY_filt = eyeData.DDY_filt;

trial.eyeDDDX = eyeData.DDDX;
trial.eyeDDDY = eyeData.DDDY;

% save some info for each trial and store it in trial.log
% for example stimulus speed, fixation duration and other events
trial.log.subject = currentSubject;
trial.log.trialNumber = currentTrial;
trialIdxInData = eventLog.trialIdxInData(currentTrial, 1);
if Experiment.const.startExp==-1
    trial.log.eyeType = 0; % fixation condition
elseif Experiment.const.startExp==1
    trial.log.eyeType = 1; % pursuit condition
end
trial.log.blockN = Experiment.trialData.blockN(trialIdxInData, 1);
trial.log.rdkApertureDir = Experiment.trialData.rdkApertureDir(trialIdxInData, 1); % either left or right, the "base" direction
trial.log.rdkApertureSpeed = Experiment.trialData.rdkApertureSpeed(trialIdxInData, 1); % positive is up, negative is down, relative to the aperture direction
trial.log.rdkApertureAngle = Experiment.trialData.rdkApertureAngle(trialIdxInData, 1);

trial.log.rdkInternalSpeed = Experiment.const.rdk.internalSpeed;
trial.log.rdkInternalCons = Experiment.trialData.rdkInternalCons(trialIdxInData, 1);
if trial.log.rdkInternalCons==-1
    trial.log.rdkInternalDir = 0;
    trial.log.rdkCoh = 0;
else
    trial.log.rdkInternalDir = Experiment.trialData.rdkInternalCons(trialIdxInData, 1); % relative direction within the RDK
    trial.log.rdkCoh = 1;
end

trial.log.shiftDis = Experiment.trialData.shiftDis(trialIdxInData, 1);
trial.log.shiftDir = Experiment.trialData.shiftDir(trialIdxInData, 1);
trial.log.shiftTime = Experiment.trialData.shiftTime(trialIdxInData, 1);

trial.log.response = Experiment.trialData.response(trialIdxInData, 1);
trial.log.eyeSampleRate = eyeData.sampleRate;

% frame indices of all events; after comparing eventLog with eyeData.frameIdx
trial.log.trialStart = 1; % the first frame, fixation onset, decided in readEyeData
trial.log.targetOnset = find(eyeData.timeStamp==eventLog.rdkOn(currentTrial, 1)); % rdk onset
trial.log.occlusionOnset = find(eyeData.timeStamp==eventLog.occlusionOn(currentTrial, 1)); % rdk onset
trial.log.occlusionOffset = find(eyeData.timeStamp==eventLog.occlusionOff(currentTrial, 1)); % rdk offset
trial.log.targetOffset = find(eyeData.timeStamp==eventLog.rdkOff(currentTrial, 1)); % rdk offset
trial.log.trialEnd = find(eyeData.timeStamp==eventLog.trialEnd(currentTrial, 1)); % response given
end