initializeParas;

checkVariables = {'response', 'gainXexternal', 'gainYexternal', 'gainYaverage', 'gain2Dexternal', 'gain2Daverage', ...
    'dirGainExternal', 'dirClp', 'dirError', 'disCenterMean'}; % always put all saccade parameters at last
saccadeVarStart = 100; % if there is no saccade variables, just use a number larger than the number of all variables to be plotted

% delete obvious error trials
idxT = find(eyeTrialData.errorStatus==0 & ...
    abs(eyeTrialData.rdkApertureAngle-eyeTrialData.response)>20 & ...
    eyeTrialData.rdkApertureAngle.*eyeTrialData.response<=0); % leftward valid trials;
eyeTrialData.errorStatus(idxT) = -2;

% delete trials with little eye movements
idxT = find(eyeTrialData.errorStatus==0 & ...
    eyeTrialData.pursuit.travelClpDis<eyeTrialData.pursuit.targetClpDis/2); 
eyeTrialData.errorStatus(idxT) = -3;

% initialization
summaryData = table;
count = 1;
%%
for subN = 1:size(names, 2)
    for internalConN = 1:size(internalCons, 2) % each column is one internal condition
        for angleN = 1:length(apertureAngles)
            summaryData.sub(count, 1) = subN;
            if internalCons(internalConN)==0
                summaryData.rdkCoh(count, 1) = 0;
                summaryData.rdkInternalDir(count, 1) = 0;
            else
                summaryData.rdkCoh(count, 1) = 1;
                summaryData.rdkInternalDir(count, 1) = internalCons(internalConN);
            end
            summaryData.rdkApertureAngle(count, 1) = apertureAngles(angleN);
            
            idxT = find(eyeTrialData.rdkCoh(subN, :)==summaryData.rdkCoh(count, 1) & ...
                eyeTrialData.rdkApertureAngle(subN, :)==summaryData.rdkApertureAngle(count, 1) & ...
                eyeTrialData.rdkInternalDir(subN, :)==summaryData.rdkInternalDir(count, 1) & ...
                eyeTrialData.errorStatus(subN, :)==0);
            
            for varN = 1:length(checkVariables)
                if strcmp(checkVariables{varN}, 'latency') % needs to calculate from onset
                    onsetT = eyeTrialData.pursuit.onset(subN, idxT);
                    rdkOnT = eyeTrialData.frameLog.rdkOn(subN, idxT);
                    
                    yMeanSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanmean(onsetT-rdkOnT);
                    yStdSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanstd(onsetT-rdkOnT);
                    if internalConN>1 % calculate difference from baseline
                        yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(onsetT-rdkOnT-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(onsetT-rdkOnT-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                    end
                elseif strcmp(checkVariables{varN}, 'dirClp') % needs to calculate from vectors
                    dir = atan2(eyeTrialData.pursuit.dirClpY(subN, idxT), eyeTrialData.pursuit.dirClpX(subN, idxT))/pi*180;
                    
                    yMeanSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanmean(dir);
                    yStdSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanstd(dir);
                    if internalConN>1 % calculate difference from baseline
                        yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(dir-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(dir-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                    end
                else
                    if strcmp(checkVariables{varN}, 'response')
                        yMeanSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanmean(eyeTrialData.(checkVariables{varN})(subN, idxT));
                        yStdSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanstd(eyeTrialData.(checkVariables{varN})(subN, idxT));
                        if internalConN>1 % calculate difference from baseline
                            yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(eyeTrialData.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(eyeTrialData.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        end
                    elseif varN>=saccadeVarStart
                        % saccade parameters
                        yMeanSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanmean(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                        yStdSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanstd(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                        if internalConN>1 % calculate difference from baseline
                            yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        end
                    else
                        % pursuit parameters
                        yMeanSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                        yStdSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                        if internalConN>1 % calculate difference from baseline
                            yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        end
                    end
                end
                summaryData.(checkVariables{varN})(count, 1) = yMeanSub.(checkVariables{varN}){subN}(angleN, internalConN);
            end
            count = count+1;
            
        end
    end
end
save('summaryDataSub.mat', 'yMeanSub', 'yStdSub', 'yMeanDiffSub', 'yStdDiffSub')
save('summaryData.mat', 'summaryData')