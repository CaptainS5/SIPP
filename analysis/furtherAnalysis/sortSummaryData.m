clear all; close all; clc

names = {'w00'};
apertureAngles = [-9 -6 -3 0 3 6 9];
internalCons = [-90, 90];
perturbPhase = [1, 2]; % 1: perturbation during initiation; 2: perturbation during steady-state
load('eyeTrialData_all2.mat');

checkVariables = {'response', 'latency', ...
    'gainXexternal', 'gainYexternal', 'gainYaverage', 'gain2Dexternal', 'gain2Daverage', 'dirGainExternal', ...
    'dirPerturb', 'dirError', 'disCenterMean', ...
    'num', 'numXLeft', 'numXRight', 'numYUp', 'numYDown', 'meanAmp2D', 'meanAmpXLeft', 'meanAmpXRight', 'meanAmpYUp', 'meanAmpYDown',...
    'sumAmp2D', 'sumAmpXLeft', 'sumAmpXRight', 'sumAmpYUp', 'sumAmpYDown'}; % always put all saccade parameters at last
% openloopVarEnd = 6;
saccadeVarStart = 12; % if there is no saccade variables, just use a number larger than the number of all variables to be plotted

% delete trials with little eye movements
% idxT = find(eyeTrialData.errorStatus==0 & ...
%     eyeTrialData.pursuit.gainXexternal<0.5);
% eyeTrialData.errorStatus(idxT) = -3;

idxT = find(eyeTrialData.errorStatus==0 & ...
    eyeTrialData.pursuit.travelClpDis<eyeTrialData.pursuit.targetClpDis/2);
eyeTrialData.errorStatus(idxT) = -3;

% initialization
summaryData = table;
% summaryDataDiff = table;
count = 1;
countDiff = 1;
%%
for subN = 1:size(names, 2)
    for perturbN = 1:length(perturbPhase)
        for internalConN = 1:size(internalCons, 2) % each column is one internal condition
            for angleN = 1:length(apertureAngles)
                summaryData.sub(count, 1) = subN;
                summaryData.perturbPhase(count, 1) = perturbPhase(perturbN);
                if internalCons(internalConN)==0
                    summaryData.rdkCoh(count, 1) = 0;
                    summaryData.rdkInternalDir(count, 1) = 0;
                else
                    summaryData.rdkCoh(count, 1) = 1;
                    summaryData.rdkInternalDir(count, 1) = internalCons(internalConN);
                end
                summaryData.rdkApertureAngle(count, 1) = apertureAngles(angleN);
                
                %             if internalConN>1 % calculate the difference from baseline
                %                 summaryDataDiff.sub(countDiff, 1) = subN;
                %                 summaryDataDiff.rdkCoh(countDiff, 1) = 1;
                %                 summaryDataDiff.rdkInternalDir(countDiff, 1) = internalCons(internalConN);
                %                 summaryDataDiff.rdkApertureAngle(countDiff, 1) = apertureAngles(angleN);
                %             end
                
                for varN = 1:length(checkVariables)
                    % if open-loop parameters, choose trials that have a valid
                    % open-loop
                    %                 if varN>=2 && varN <=openloopVarEnd % open-loop parameters
                    %                     idxT = find(eyeTrialData.rdkCoh(subN, :)==summaryData.rdkCoh(count, 1) & ...
                    %                         eyeTrialData.rdkApertureAngle(subN, :)==summaryData.rdkApertureAngle(count, 1) & ...
                    %                         eyeTrialData.rdkInternalDir(subN, :)==summaryData.rdkInternalDir(count, 1) & ...
                    %                         eyeTrialData.pursuit.onsetType(subN, :)==0 & ...
                    %                         eyeTrialData.errorStatus(subN, :)==0);
                    %                 else % perception, closed-loop, or saccades
                    if perturbPhase(perturbN)==1 % during initiation, delete those started after openloop
                        idxT = find(eyeTrialData.rdkCoh(subN, :)==summaryData.rdkCoh(count, 1) & ...
                            eyeTrialData.rdkApertureAngle(subN, :)==summaryData.rdkApertureAngle(count, 1) & ...
                            eyeTrialData.rdkInternalDir(subN, :)==summaryData.rdkInternalDir(count, 1) & ...
                            eyeTrialData.perturbPhase(subN, :)==summaryData.perturbPhase(count, 1) & ...
                            eyeTrialData.pursuit.perturbationLabel(subN, :)>=0 & ...
                            eyeTrialData.errorStatus(subN, :)==0);
                    else
                        idxT = find(eyeTrialData.rdkCoh(subN, :)==summaryData.rdkCoh(count, 1) & ...
                            eyeTrialData.rdkApertureAngle(subN, :)==summaryData.rdkApertureAngle(count, 1) & ...
                            eyeTrialData.rdkInternalDir(subN, :)==summaryData.rdkInternalDir(count, 1) & ...
                            eyeTrialData.perturbPhase(subN, :)==summaryData.perturbPhase(count, 1) & ...
                            eyeTrialData.errorStatus(subN, :)==0);
                    end
                    %                 end
                    %
                    if strcmp(checkVariables{varN}, 'latency') % needs to calculate from onset
                        onsetT = eyeTrialData.pursuit.onset(subN, idxT);
                        rdkOnT = eyeTrialData.frameLog.rdkOn(subN, idxT);
                        
                        yMeanSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = nanmean(onsetT-rdkOnT);
                        yStdSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = nanstd(onsetT-rdkOnT);
                        %                     if internalConN>1 % calculate difference from baseline
                        %                         yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(onsetT-rdkOnT-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        %                         yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(onsetT-rdkOnT-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        %                     end
                    elseif strcmp(checkVariables{varN}, 'dirPerturb') % || strcmp(checkVariables{varN}, 'dirOlp') ...
                        %                         || strcmp(checkVariables{varN}, 'dirClpEarly') || strcmp(checkVariables{varN}, 'dirClpLate') % needs to calculate from vectors
                        dirXName = [checkVariables{varN}, 'X'];
                        dirYName = [checkVariables{varN}, 'Y'];
                        dir = atan2(eyeTrialData.pursuit.(dirYName)(subN, idxT), eyeTrialData.pursuit.(dirXName)(subN, idxT))/pi*180;
                        
                        %                     % delete extreme values for Olp pursuit direction... or
                        %                     % better, just do not analyze olp pursuit or really
                        %                     % figure out the problem with finding pursuit onset
                        %                     if strcmp(checkVariables{varN}, 'dirOlp')
                        %                         dir(abs(dir)>90)=[]; % has to be wrong, shouldn't be moving to the left...
                        %                     end
                        
                        yMeanSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = nanmean(dir);
                        yStdSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = nanstd(dir);
                        %                     if internalConN>1 % calculate difference from baseline
                        %                         yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(dir-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        %                         yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(dir-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                        %                     end
                    else
                        if strcmp(checkVariables{varN}, 'response')
                            % rate of reporting up
                            yMeanSub.choiceUpNum{subN, perturbN}(angleN, internalConN) = length(find(eyeTrialData.(checkVariables{varN})(subN, idxT)==1));
                            yMeanSub.totalTrialN{subN, perturbN}(angleN, internalConN) = length(idxT);
                            yMeanSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = length(find(eyeTrialData.(checkVariables{varN})(subN, idxT)==1))/length(idxT);
                            yStdSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = NaN;
                            %                         % Mean values of angles
                            %                         yMeanSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanmean(eyeTrialData.(checkVariables{varN})(subN, idxT));
                            %                         yStdSub.(checkVariables{varN}){subN}(angleN, internalConN) = nanstd(eyeTrialData.(checkVariables{varN})(subN, idxT));
                            %                         if internalConN>1 % calculate difference from baseline
                            %                             yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(eyeTrialData.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            %                             yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(eyeTrialData.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            %                         end
                        elseif varN>=saccadeVarStart
                            % saccade parameters
                            yMeanSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = nanmean(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                            yStdSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = nanstd(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                            %                         if internalConN>1 % calculate difference from baseline
                            %                             yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            %                             yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            %                         end
                        else
                            % pursuit parameters
                            yMeanSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                            yStdSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                            %                         if internalConN>1 % calculate difference from baseline
                            %                             yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            %                             yStdDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT)-yMeanSub.(checkVariables{varN}){subN}(angleN, 1));
                            %                         end
                        end
                    end
                    summaryData.(checkVariables{varN})(count, 1) = yMeanSub.(checkVariables{varN}){subN, perturbN}(angleN, internalConN);
                    %                 if internalConN>1 % calculate the difference from baseline
                    %                     summaryDataDiff.(checkVariables{varN})(countDiff, 1) = yMeanDiffSub.(checkVariables{varN}){subN}(angleN, internalConN-1);
                    %                 end
                end
                count = count+1;
                %             if internalConN>1 % calculate the difference from baseline
                %                 countDiff = countDiff+1;
                %             end
            end
        end
    end
end
save('summaryDataSub.mat', 'yMeanSub', 'yStdSub'); %, 'yMeanDiffSub', 'yStdDiffSub')
save('summaryData.mat', 'summaryData')
% save('summaryDataDiff.mat', 'summaryDataDiff')

%% generate csv for plotting in R
% writetable(summaryData, [RFolder, 'summaryData.csv'])
% writetable(summaryDataDiff, [RFolder, 'summaryDataDiff.csv'])
