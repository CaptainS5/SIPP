% prepare trial by trial data for correlation analysis
initializeParas;

%%
idxV = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkInternalCon~=0);
trialT = table;
for ii = 1:length(idxV)
    trialT.sub(ii, 1) = str2num(eyeTrialData.sub{idxV(ii)}(2:end));
    trialT.rdkInternalDir(ii, 1) = eyeTrialData.rdkInternalCon(idxV(ii));
    trialT.rdkApertureAngle(ii, 1) = eyeTrialData.rdkApertureAngle(idxV(ii));
    
    % find baseline
    idx = find(summaryData.sub==trialT.sub(ii, 1) & summaryData.rdkInternalDir==0 & summaryData.rdkApertureAngle==trialT.rdkApertureAngle(ii, 1));
    baseClp = summaryData.dirClp(idx);
    basePerceptual = summaryData.response(idx);
    
    % calculate bias of the current trial
    trialClp = atan2(eyeTrialData.pursuit.dirClpY(idxV(ii)), eyeTrialData.pursuit.dirClpX(idxV(ii)))/pi*180;
    trialT.dirClp(ii, 1) = trialClp-baseClp;
    trialT.response(ii, 1) = eyeTrialData.response(idxV(ii))-basePerceptual;
end

writetable(trialT, [RFolder, 'trialDataCorr.csv'])