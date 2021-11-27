% Exp2, Xiuyun Wu, 22/Nov/2021
% getting the raw processed eye data... will be much more convenient for
% later analysis; run this after getting the errorfiles

clear all; close all; clc

names = {'w01'};
subStartI = 1;

cd ..
cd ..
analysisPath = pwd; % folder for the eye movement preprocessing codes
dataPath = ['..\data\']; % still need to go into specific folders

%% All trials
if subStartI>1 % if starting halway, load the current eyeTrialDataAll
    load([analysisPath '\furtherAnalysis\exp2\eyeTrialData_all.mat'])
else
    clear eyeTrialData
end

for subN = subStartI:length(names)
    % define the sub folder
    currentSubject = names{subN}
    currentSubjectPath = [dataPath, currentSubject];
    % load files
    load([currentSubjectPath, '\eventLog.mat']);
    load([currentSubjectPath, '\rdkFrameLog.mat']);
    load([currentSubjectPath, '\info_Experiment.mat'])
    
    % error file
    errorFilePath = [analysisPath, '\ErrorFiles\'];
    load([errorFilePath, 'Sub_', currentSubject, '_errorFile.mat']);
    
    clear eyeTrialDataSub
    
    % saccade algorithm threshold
    saccadeThreshold = 400; % acceleration
    microSaccadeThreshold = 200;
    
    for currentTrial = 1:size(errorStatus, 1)
        trialIdxInData = eventLog.trialIdxInData(currentTrial, 1);
        eyeTrialData.sub{subN, currentTrial} = currentSubject;
        eyeTrialData.trialIdx(subN, currentTrial) = currentTrial;
        eyeTrialData.trialIdxInData(subN, currentTrial) = trialIdxInData;
        eyeTrialData.blockN(subN, currentTrial) = Experiment.trialData.blockN(trialIdxInData, 1);
        eyeTrialData.errorStatus(subN, currentTrial) = errorStatus(currentTrial, 1);
        
        analyzeTrial;
        eyeTrialData.rdkApertureDir(subN, currentTrial) = trial.log.rdkApertureDir; % direction before perturbation, right (0) or left (180)
        eyeTrialData.rdkApertureSpeed(subN, currentTrial) = trial.log.rdkApertureSpeed;
        eyeTrialData.rdkApertureAngle(subN, currentTrial) = trial.log.rdkApertureAngle; % angle during perturbation
        eyeTrialData.rdkInternalSpeed(subN, currentTrial) = trial.log.rdkInternalSpeed; % during perturbation
        eyeTrialData.rdkInternalCon(subN, currentTrial) = trial.log.rdkInternalCon;
        eyeTrialData.rdkInternalDir(subN, currentTrial) = trial.log.rdkInternalDir; % during perturbation
        eyeTrialData.rdkCoh(subN, currentTrial) = trial.log.rdkCoh;
        eyeTrialData.perturbPhase(subN, currentTrial) = trial.log.perturbTime;
        eyeTrialData.response(subN, currentTrial) = trial.log.response;
        
        eyeTrialData.target{subN, currentTrial} = trial.target;
        
        if errorStatus(currentTrial, 1)==0
            eyeTrialData.frameLog.fixationOn(subN, currentTrial) = trial.log.trialStart;
            eyeTrialData.frameLog.rdkOn(subN, currentTrial) = trial.log.targetOnset;
            eyeTrialData.frameLog.perturbationOn(subN, currentTrial) = trial.log.perturbationOnset;
            eyeTrialData.frameLog.rdkOff(subN, currentTrial) = trial.log.targetOffset;
            
            % all pursuit data
            fields = fieldnames(trial.pursuit.summary);
            for ii = 1:length(fields)
                eyeTrialData.pursuit.(fields{ii})(subN, currentTrial) = trial.pursuit.summary.(fields{ii});
            end
            % all saccade summary data
            fields = fieldnames(trial.saccades.summary);
            for ii = 1:length(fields)
                eyeTrialData.saccades.(fields{ii})(subN, currentTrial) = trial.saccades.summary.(fields{ii});
            end
            eyeTrialData.saccades.dirXY{subN, currentTrial} = trial.saccades.dirXY;
            eyeTrialData.saccades.dirAngle{subN, currentTrial} = trial.saccades.dirAngle;
            eyeTrialData.saccades.eyePosOnset{subN, currentTrial} = trial.saccades.eyePosOnset;
            eyeTrialData.saccades.eyePosOffset{subN, currentTrial} = trial.saccades.eyePosOffset;
            eyeTrialData.saccades.targetPosOnset{subN, currentTrial} = trial.saccades.targetPosOnset;
            eyeTrialData.saccades.targetPosOffset{subN, currentTrial} = trial.saccades.targetPosOffset;
            
            eyeTrialDataSub.trial{1, currentTrial}.eyeX_filt = trial.eyeX_filt; % for velocity traces
            eyeTrialDataSub.trial{1, currentTrial}.eyeY_filt = trial.eyeY_filt;
            eyeTrialDataSub.trial{1, currentTrial}.eyeDX_filt = trial.eyeDX_filt;
            eyeTrialDataSub.trial{1, currentTrial}.eyeDY_filt = trial.eyeDY_filt;
            eyeTrialDataSub.trial{1, currentTrial}.X_noSac = trial.X_noSac;
            eyeTrialDataSub.trial{1, currentTrial}.Y_noSac = trial.Y_noSac;
            eyeTrialDataSub.trial{1, currentTrial}.DX_noSac = trial.DX_noSac;
            eyeTrialDataSub.trial{1, currentTrial}.DY_noSac = trial.DY_noSac;
            eyeTrialDataSub.trial{1, currentTrial}.X_interpolSac = trial.X_interpolSac;
            eyeTrialDataSub.trial{1, currentTrial}.Y_interpolSac = trial.Y_interpolSac;
            eyeTrialDataSub.trial{1, currentTrial}.DX_interpolSac = trial.DX_interpolSac;
            eyeTrialDataSub.trial{1, currentTrial}.DY_interpolSac = trial.DY_interpolSac;
            eyeTrialDataSub.trial{1, currentTrial}.targetVelX = trial.target.velX;
            eyeTrialDataSub.trial{1, currentTrial}.targetVelY = trial.target.velY;
            eyeTrialDataSub.trial{1, currentTrial}.pursuit.dirVec = trial.pursuit.eyeDir;
        else
            eyeTrialData.frameLog.fixationOn(subN, currentTrial) = NaN;
            eyeTrialData.frameLog.rdkOn(subN, currentTrial) = NaN;
            eyeTrialData.frameLog.perturbationOn(subN, currentTrial) = NaN;
            eyeTrialData.frameLog.rdkOff(subN, currentTrial) = NaN;

            fields = fieldnames(trial.pursuit.summary);
            for ii = 1:length(fields)
                eyeTrialData.pursuit.(fields{ii})(subN, currentTrial) = NaN;
            end
            fields = fieldnames(trial.saccades.summary);
            for ii = 1:length(fields)
                eyeTrialData.saccades.(fields{ii})(subN, currentTrial) = NaN;
            end
            eyeTrialData.saccades.dirXY{subN, currentTrial} = NaN;
            eyeTrialData.saccades.dirAngle{subN, currentTrial} = NaN;
            eyeTrialData.saccades.eyePosOnset{subN, currentTrial} = NaN;
            eyeTrialData.saccades.eyePosOffset{subN, currentTrial} = NaN;
            eyeTrialData.saccades.targetPosOnset{subN, currentTrial} = NaN;
            eyeTrialData.saccades.targetPosOffset{subN, currentTrial} = NaN;
            
            eyeTrialDataSub.trial{1, currentTrial} = NaN; % for velocity traces
        end
    end
    save([analysisPath '\furtherAnalysis\exp2\eyeTrialDataSub_' names{subN} '.mat'], 'eyeTrialDataSub');
end
save([analysisPath '\furtherAnalysis\exp2\eyeTrialData_all.mat'], 'eyeTrialData');
