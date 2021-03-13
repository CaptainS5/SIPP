% FUNCTION to read in information about the target in each trial; align it
% with the eye data samples
%
% history
% 11-Mar-2021   created by XW for recovering the translating RDK in the micropursuit exp; xiuyunwu5@gmail.com

function [target] = readoutTarget(eyeData, dotSpeed, currentSubjectPath, currentTrial, eventLog, rdkFrameLog)
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

end