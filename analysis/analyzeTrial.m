%% THIS m-file is the magic script that analyzes all your eye data
% this script requries several functions:
% readEyeData.m, processEyeData.m, readoutTrial.m, findSaccades.m,
% analyzeSaccades.m
% optional: findPursuit.m, analyzePursuit.m, removeSaccades.m,
% findMicroSaccades.m
% optional: we also have scripts to read and analyze target and finger data
% for the EyeCatch/EyeStrike paradigm that can be added here --> ask if you
% need more info

% history
% 07-2012       JE created analyzeTrial.m
% 2012-2018     JF added stuff to and edited analyzeTrial.m
% 13-07-2018    JF commented to make the script more accecable for future
%               VPOM students
% for questions email jolande.fooken@rwth-aachen.de
% 13-Jan-2021   XW modified based on using findSaccadeAcc.m and findPursuitNew.m; xiuyunwu5@gmail.com

%% Eye Data
%  eye data need to have been converted using convert_edf2asc.m
eyeFile = [currentSubject 't' num2str(currentTrial, '%03d') '.mat']; % mat file, eye data transformed from edf
% make sure they are included in the experiment code
eyeData = readEyeData(eyeFile, dataPath, currentSubject, currentTrial, analysisPath, eventLog, Experiment);
if length(eyeData.timeStamp)<=10
    trial.signalLoss = 1;
else
    eyeData = processEyeData(eyeData);
    
    %% extract all relevant experimental data and store it in trial variable
    trial = readoutTrial(eyeData, currentSubject, analysisPath, Experiment, currentTrial, eventLog);
%     trial.target = readoutTarget(eyeData, Experiment.const.rdk.apertureSpeed, currentSubjectPath, currentTrial, eventLog, rdkFrameLog);
    trial.stim_onset = trial.log.targetOnset;
    trial.stim_offset = trial.log.targetOffset;
    trial.length = trial.log.trialEnd;
    
    if isempty(trial.stim_offset)
        trial.signalLoss = 1;
    else
        trial.signalLoss = 0;
        %% find saccades
        saccadeThresholds = [evalin('base', 'saccadeThreshold'), evalin('base', 'microSaccadeThreshold')];
        onset = 1;
        offset = trial.length; %min(trial.stim_offset, size(trial.eyeDX_filt, 1)); % to be able to detect saccades at the end of display
        
        % use acceleration to find saccades...
        [saccades.X.onsets, saccades.X.offsets] = findSaccadesAcc(onset, offset, trial.eyeDX_filt, trial.eyeDDX_filt, trial.eyeDDDX, saccadeThresholds, trial.log);
        [saccades.Y.onsets, saccades.Y.offsets] = findSaccadesAcc(onset, offset, trial.eyeDY_filt, trial.eyeDDY_filt, trial.eyeDDDY, saccadeThresholds, trial.log);
        
        % remove saccades
        trial = removeSaccades(trial, saccades);
        clear saccades;
        
%         if trial.log.eyeCondition==1 % pursuit trials
            %% find and analyze pursuit
            pursuit = findPursuitNew(trial);
            % analyze pursuit
            trial = analyzePursuit(trial, pursuit);
            
            %% analyze saccades
            trial = analyzeSaccades(trial);
%         else
            %% OPTIONAL: find micro saccades
            % % remove saccades
            % trial = removeSaccades(trial);
            % m_threshold = evalin('base', 'microSaccadeThreshold');
            % [saccades.X.onsets, saccades.X.offsets] = findSaccades(onset, offset, trial.DX_noSac, trial.DDX_noSac, m_threshold, 0);
            % [saccades.Y.onsets, saccades.Y.offsets] = findSaccades(onset, offset, trial.DY_noSac, trial.DDY_noSac, m_threshold, 0);
            % % analyze micro-saccades
            % [trial] = analyzeMicroSaccades(trial, saccades);
%         end
    end
end

if trial.signalLoss
    % for updateText
    trial.log.subject = currentSubject;
    trial.log.trialNumber = currentTrial;
    trialIdxInData = eventLog.trialIdxInData(currentTrial, 1);
    trial.log.blockN = Experiment.trialData.blockN(trialIdxInData, 1);
    trial.log.rdkApertureDir = Experiment.trialData.rdkApertureDir(trialIdxInData, 1); % positive is up, negative is down
    trial.log.rdkInternalDir = Experiment.trialData.rdkInternalDir(trialIdxInData, 1); % direction std
    trial.log.rdkInternalSpeed = Experiment.const.rdk.internalSpeed; %Experiment.trialData.rdkInternalSpeed(trialIdxInData, 1);
    trial.log.rdkCoh = Experiment.trialData.rdkCoh(trialIdxInData, 1);
    trial.stimulus.absoluteVelocity = Experiment.const.rdk.apertureSpeed;
end
