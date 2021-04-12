% Xiuyun Wu, 13/Jan/2021
% getting the raw processed eye data... will be much more convenient for
% later analysis; run this after getting the errorfiles

clear all; close all; clc

names = {'w00' 'w01' 'w02' 'w03' 'w04' 'w05' 'w06' 'w07' 'w08' 'w09' 'w10'};
subStartI = 1;

cd ..
analysisPath = pwd; % folder for the eye movement preprocessing codes
dataPath = ['..\data\']; % still need to go into specific folders

%% All trials
if subStartI>1 % if starting halway, load the current eyeTrialDataAll
    load([analysisPath '\furtherAnalysis\eyeTrialData_all.mat'])
else
    clear eyeTrialData
end

for subN = subStartI:length(names)
    % define the sub folder
    currentSubject = names{subN}
    currentSubjectPath = [dataPath, currentSubject];
    % load files
    load([currentSubjectPath, '\eventLog.mat']);
    load([currentSubjectPath, '\info_Experiment.mat'])
    
    % error file
    errorFilePath = [analysisPath, '\ErrorFiles\'];
    load([errorFilePath, 'Sub_', currentSubject, '_errorFile.mat']);
    
    clear eyeTrialDataSub
    
    if strcmp(currentSubject, 'w07') || strcmp(currentSubject, 'w10') % fixation condition...
        % saccade algorithm threshold
        saccadeThreshold = 300; % acceleration
        microSaccadeThreshold = 200;
    else
        % saccade algorithm threshold
        saccadeThreshold = 500; % acceleration
        microSaccadeThreshold = 200;
    end
    
    for currentTrial = 1:size(errorStatus, 1)
        trialIdxInData = eventLog.trialIdxInData(currentTrial, 1);
        eyeTrialData.sub{subN, currentTrial} = currentSubject;
        eyeTrialData.trialIdx(subN, currentTrial) = currentTrial;
        eyeTrialData.trialIdxInData(subN, currentTrial) = trialIdxInData;
        eyeTrialData.blockN(subN, currentTrial) = Experiment.trialData.blockN(trialIdxInData, 1);
        eyeTrialData.errorStatus(subN, currentTrial) = errorStatus(currentTrial, 1);
        
        if errorStatus(currentTrial, 1)==0
            analyzeTrial;
            % to get target info
            eyeTrialData.rdkApertureDir(subN, currentTrial) = trial.log.rdkApertureDir; % positive is up, negative is down
            eyeTrialData.rdkInternalDir(subN, currentTrial) = trial.log.rdkInternalDir; % direction std
            eyeTrialData.rdkInternalSpeed(subN, currentTrial) = trial.log.rdkInternalSpeed; %
            eyeTrialData.rdkCoh(subN, currentTrial) = trial.log.rdkCoh;
            
            eyeTrialData.frameLog.fixationOn(subN, currentTrial) = trial.log.trialStart;
            eyeTrialData.frameLog.rdkOn(subN, currentTrial) = trial.log.targetOnset;
            eyeTrialData.frameLog.rdkOff(subN, currentTrial) = trial.log.targetOffset;
            %             eyeTrialData.frameLog.respond(subN, currentTrial) = trial.log.trialEnd;
            %             eyeTrialData.target{subN, currentTrial} = trial.target;
            
            % all pursuit data
            fields = fieldnames(trial.pursuit);
            for ii = 1:length(fields)
                eyeTrialData.pursuit.(fields{ii})(subN, currentTrial) = trial.pursuit.(fields{ii});
            end
            
            eyeTrialData.saccades.number(subN, currentTrial) = trial.saccades.num;
            eyeTrialData.saccades.meanAmp2D(subN, currentTrial) = trial.saccades.meanAmp2D;
            eyeTrialData.saccades.meanAmpXLeft(subN, currentTrial) = trial.saccades.meanAmpXLeft;
            eyeTrialData.saccades.meanAmpXRight(subN, currentTrial) = trial.saccades.meanAmpXRight;
            eyeTrialData.saccades.meanAmpYUp(subN, currentTrial) = trial.saccades.meanAmpYUp;
            eyeTrialData.saccades.meanAmpYDown(subN, currentTrial) = trial.saccades.meanAmpYDown;
            eyeTrialData.saccades.sumAmp2D(subN, currentTrial) = trial.saccades.sumAmp2D;

            %                 % record saccades in both directions...
            % %                 if trial. log.rdkDir>0 || (trial.log.rdkDir==0 && trial.pursuit.closedLoopMeanVelX>=0)% first use rdk dir to judge, then see pursuit; trial.pursuit.closedLoopMeanVelX>=0 % right ward pursuit
            %                     eyeTrialData.saccades.X_right.number(subN, currentTrial) = trial.saccades.X_right.number;
            %                     eyeTrialData.saccades.X_right.meanAmplitude(subN, currentTrial) = trial.saccades.X_right.meanAmplitude;
            %                     %                 eyeTrialData.saccades.X.maxAmplitude(subN, currentTrial) = trial.saccades.X_right.maxAmplitude;
            %                     eyeTrialData.saccades.X_right.meanDuration(subN, currentTrial) = trial.saccades.X_right.meanDuration;
            %                     eyeTrialData.saccades.X_right.sumAmplitude(subN, currentTrial) = trial.saccades.X_right.sumAmplitude;
            %                     %                 eyeTrialData.saccades.X.peakVelocity(subN, currentTrial) = trial.saccades.X_right.peakVelocity;
            %                     %                 eyeTrialData.saccades.X.meanVelocity(subN, currentTrial) = trial.saccades.X_right.meanVelocity;
            %                     eyeTrialData.saccades.X_right.onsets_pursuit{subN, currentTrial} = trial.saccades.X_right.onsets_pursuit;
            %                     eyeTrialData.saccades.X_right.offsets_pursuit{subN, currentTrial} = trial.saccades.X_right.offsets_pursuit;
            % %                 elseif trial. log.rdkDir<0 || (trial.log.rdkDir==0 && trial.pursuit.closedLoopMeanVelX<0)
            %                     eyeTrialData.saccades.X_left.number(subN, currentTrial) = trial.saccades.X_left.number;
            %                     eyeTrialData.saccades.X_left.meanAmplitude(subN, currentTrial) = trial.saccades.X_left.meanAmplitude;
            %                     %                 eyeTrialData.saccades.X.maxAmplitude(subN, currentTrial) = trial.saccades.X_left.maxAmplitude;
            %                     eyeTrialData.saccades.X_left.meanDuration(subN, currentTrial) = trial.saccades.X_left.meanDuration;
            %                     eyeTrialData.saccades.X_left.sumAmplitude(subN, currentTrial) = trial.saccades.X_left.sumAmplitude;
            %                     %                 eyeTrialData.saccades.X.peakVelocity(subN, currentTrial) = trial.saccades.X_left.peakVelocity;
            %                     %                 eyeTrialData.saccades.X.meanVelocity(subN, currentTrial) = trial.saccades.X_left.meanVelocity;
            %                     eyeTrialData.saccades.X_left.onsets_pursuit{subN, currentTrial} = trial.saccades.X_left.onsets_pursuit;
            %                     eyeTrialData.saccades.X_left.offsets_pursuit{subN, currentTrial} = trial.saccades.X_left.offsets_pursuit;
            % %                 end
            % %                 if ~isnan(eyeTrialData.saccades.X.sumAmplitude(subN, currentTrial))
            % %                     eyeTrialData.pursuit.gainSacSumAmpX(subN, currentTrial) = eyeTrialData.pursuit.gainX(subN, currentTrial)+eyeTrialData.saccades.X.sumAmplitude(subN, currentTrial)/10;
            % %                 else
            % %                     eyeTrialData.pursuit.gainSacSumAmpX(subN, currentTrial) = eyeTrialData.pursuit.gainX(subN, currentTrial);
            % %                 end
            
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
            eyeTrialDataSub.trial{1, currentTrial}.targetVelX = trial.target.velocityX;
            eyeTrialDataSub.trial{1, currentTrial}.targetVelY = trial.target.velocityY;
        else
            eyeTrialData.frameLog.fixationOn(subN, currentTrial) = NaN;
            eyeTrialData.frameLog.rdkOn(subN, currentTrial) = NaN;
            eyeTrialData.frameLog.rdkOff(subN, currentTrial) = NaN;
            %             eyeTrialData.frameLog.respond(subN, currentTrial) = NaN;
            %             eyeTrialData.target{subN, currentTrial} = trial.target;
            
            eyeTrialData.rdkApertureDir(subN, currentTrial) = NaN; % positive is up, negative is down
            eyeTrialData.rdkInternalDir(subN, currentTrial) = NaN; % direction std
            eyeTrialData.rdkInternalSpeed(subN, currentTrial) = NaN; %
            eyeTrialData.rdkCoh(subN, currentTrial) = NaN;
            
            fields = fieldnames(trial.pursuit);
            for ii = 1:length(fields)
                eyeTrialData.pursuit.(fields{ii})(subN, currentTrial) = NaN;
            end
            
            eyeTrialData.saccades.number(subN, currentTrial) = NaN;
            eyeTrialData.saccades.meanAmp2D(subN, currentTrial) = NaN;
            eyeTrialData.saccades.meanAmpXLeft(subN, currentTrial) = NaN;
            eyeTrialData.saccades.meanAmpXRight(subN, currentTrial) = NaN;
            eyeTrialData.saccades.meanAmpYUp(subN, currentTrial) = NaN;
            eyeTrialData.saccades.meanAmpYDown(subN, currentTrial) = NaN;
            eyeTrialData.saccades.sumAmp2D(subN, currentTrial) = NaN;
            
            eyeTrialDataSub.trial{1, currentTrial} = NaN; % for velocity traces
        end
    end
    save([analysisPath '\furtherAnalysis\eyeTrialDataSub_' names{subN} '.mat'], 'eyeTrialDataSub');
end
save([analysisPath '\furtherAnalysis\eyeTrialData_all.mat'], 'eyeTrialData');
