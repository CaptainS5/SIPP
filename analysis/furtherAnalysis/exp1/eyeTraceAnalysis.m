% Exp1
% analyze the eye traces further to see the temporal dynamics in pursuit
% direction...
initializeParas;
%%
load eyeTracesAll.mat
% the organization in indiMean and allMean are:
% indiMean{groupN}.vel/pos(Diff){angleN, internalConN}{dimN}(subN, :)
% allMean{groupN}.vel/pos(Diff){angleN, internalConN}{dimN}
% for diff, thr first internalConN column is empty (baseline condition)
% time point is from -199 ms before RDK onset, one frame per second 

subStart = 1;
subEnd = 20;

groupName = {'visualDir'};
% naming by trial type (could include grouping rules) + group based on which direction (visual or perceived)
groupN = [1];
timePoints = -199:1:787;

%% calculate the individual "turning point" in y vel diff
for ii = 1:length(groupN)
    for subN = subStart:subEnd
        for internalConN = 2:length(internalCons)
            dataTemp = [];
            for angleN = 1:length(apertureAngles)
                dataTemp = [dataTemp; indiMean{ii}.velDiff{angleN, internalConN}{2}(subN, :)]; % y vel diff traces
%                 dataTemp = indiMean{ii}.velDiff{angleN, internalConN}{2}(subN, :);
%                 if internalCons(internalConN)<0
%                     tp = find(dataTemp == min(dataTemp))-200;
%                 else
%                     tp = find(dataTemp == max(dataTemp))-200;
%                 end
%                 %             figure
%                 %             hold on
%                 %             plot(timePoints, dataTempM)
%                 %             line([tp, tp], [-0.5, 0.5])
%                 %             close all
%                 
%                 % feed into summaryDataDiff
%                 idxS = find(summaryDataDiff.sub==subN & summaryDataDiff.rdkInternalDir==internalCons(internalConN) & summaryDataDiff.rdkApertureAngle==apertureAngles(angleN));
%                 summaryDataDiff.turningPoint(idxS, 1) = tp;
            end
            dataTempM = nanmean(dataTemp);
            tp = find(abs(dataTempM) == max(abs(dataTempM)))-200;
            %             figure
            %             hold on
            %             plot(timePoints, dataTempM)
            %             line([tp, tp], [-0.5, 0.5])
            %             close all
            turningPoint(subN, internalConN) = tp;
            
            % feed into summaryDataDiff
            idxS = find(summaryDataDiff.sub==subN & summaryDataDiff.rdkInternalDir==internalCons(internalConN));
            summaryDataDiff.turningPoint(idxS, 1) = tp;
        end
    end
end
save('summaryDataDiff.mat', 'summaryDataDiff')

% now... well just rerun sortEyeDataRaw to get the updated pursuit
% parameters
