% trial analysis regarding pursuit direction change
initializeParas;

% % choose which plot to look at now
% individualPlots = 1;
% averagedPlots = 0;
% subStart = 12;
% subEnd = 14;

%% first, descriptive plots of the direction change situation within each participant
for subN = 1:size(names, 2)
%     % if you want to sort all direction vectors in each trial... edit codes
%     % below
%     cd(analysisFolder)
%     load(['eyeTrialDataSub_' names{subN} '.mat']);
%     
%     % find all valid trials
%     idxT = find(eyeTrialData.errorStatus(subN, :)==0);
%     
%     for trialN = 1:length(idx)
%         ySubSummary{subN}.rdkApertureAngle(trialN) = eyeTrialData.rdkApertureAngle(subN, idxT(trialN));
%         ySub{subN}.internalCon(trialN) = eyeTrialData.internalCon(subN, idxT(trialN));
%         ySub{subN}.response(trialN) = eyeTrialData.response(subN, idxT(trialN));
%         
%     end


% plot the direction change distribution in clp for each participant
figure
for angleN = 1:7
    % calculate the direction change in baseline as the baseline
idxT = find(eyeTrialData.errorStatus(subN, :)==0 & ...
    eyeTrialData.rdkInternalCon(subN, :)==0 & ...
    eyeTrialData.rdkApertureAngle(subN, :)==apertureAngles(angleN));
baselineDirClpChange(subN, angleN) = nanmean(eyeTrialData.pursuit.dirClpChange(subN, idxT));

% calculate difference in direction changes in each


    subplot(3, 3, angleN)
    hist()
    title(['angle ', num2str(apertureAngles(angleN))])
    xlabel()
end
end