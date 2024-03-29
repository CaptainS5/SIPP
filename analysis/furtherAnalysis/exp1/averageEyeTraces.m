% Exp 1: plot position&velocity traces for translating RDK
% also generate csv file for further plotting in R

% currently each figure is one internal motion condition; plot the difference in
% velocity traces between the exp conditions and the baseline (static)
% condition
initializeParas;

% choose which plot to look at now
individualPlots = 0;
averagedPlots = 1;
subStart = 1;
subEnd = 20;

groupName = {'visualDir'};
% naming by trial type (could include grouping rules) + group based on which direction (visual or perceived)
groupN = [1]; % corresponds to the listed rules... can choose multiple, just list as a vector
% when choosing multiple groupN, will plot each group rule in one figure

% prepare for the grouping...
allCons.internalCons = internalCons;
allCons.apertureAngles = apertureAngles;

% plotting settings
dimNames = {'Horizontal', 'Vertical'};
yRangeX = [6 12];
yRangeY = [-3 3];

%% align RDK onset, frame data for all trials
frameLength = NaN(size(names, 2), 1);
latency = NaN(size(names));
for subN = subStart:subEnd
    cd(analysisFolder)
    load(['eyeTrialDataSub_' names{subN} '.mat']);
    
    idxT = find(eyeTrialData.errorStatus(subN, :)==0);
    tempL = eyeTrialData.frameLog.rdkOff(subN, idxT)-eyeTrialData.frameLog.rdkOn(subN, idxT)+200; % starting from 200 ms before RDK onset, to the end of RDK
    tempL(tempL==0) = [];
    frameLength(subN, 1) = min(tempL);
    
    for internalConN = 1:length(allCons.internalCons)
        if allCons.internalCons(internalConN)==0
            rdkCoh = 0;
            rdkInternalDir = 0;
        else
            rdkCoh = 1;
            rdkInternalDir = allCons.internalCons(internalConN);
        end
        idxT = find(eyeTrialData.errorStatus(subN, :)==0 & ...
            eyeTrialData.rdkCoh(subN, :)==rdkCoh & ...
            eyeTrialData.rdkInternalDir(subN, :)==rdkInternalDir);
        %                 eyeTrialData.rdkApertureAngle(subN, :)==allCons.apertureAngles(angleN));
        
        % latency
        latency(subN, internalConN) = nanmean(eyeTrialData.pursuit.onset(subN, idxT)-eyeTrialData.frameLog.rdkOn(subN, idxT));
        
        lengthT = length(idxT);
        frames.onset{subN, internalConN}.posX = NaN(lengthT, frameLength(subN, 1));
        frames.onset{subN, internalConN}.posY = NaN(lengthT, frameLength(subN, 1));
        frames.onset{subN, internalConN}.velX = NaN(lengthT, frameLength(subN, 1));
        frames.onset{subN, internalConN}.velY = NaN(lengthT, frameLength(subN, 1));
        
        for trialN = 1:lengthT
            % align at 200 ms before RDK onset
            startI = eyeTrialData.frameLog.rdkOn(subN, idxT(trialN))-200+1;
            endI = startI+frameLength(subN, 1)-1;
            frames.onset{subN, internalConN}.posX(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.X_interpolSac(startI:endI);
            frames.onset{subN, internalConN}.posX(trialN, :) = frames.onset{subN, internalConN}.posX(trialN, :)-frames.onset{subN, internalConN}.posX(trialN, 1); % all start from 0
            frames.onset{subN, internalConN}.posY(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.Y_interpolSac(startI:endI);
            frames.onset{subN, internalConN}.posY(trialN, :) = frames.onset{subN, internalConN}.posY(trialN, :)-frames.onset{subN, internalConN}.posY(trialN, 1);
            frames.onset{subN, internalConN}.velX(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.DX_interpolSac(startI:endI);
            frames.onset{subN, internalConN}.velY(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.DY_interpolSac(startI:endI);
        end
        %         end
    end
end
maxFrameLength = max(frameLength);

% plotting parameters
minFrameLength = min(frameLength);
framePerSec = 1/sampleRate;
timePointsOnset = [(-199:1:(minFrameLength-200))]*framePerSec*1000; % rdk onset is 0

%% calculate mean traces
for ii = 1:length(groupN)
    [indiMean{ii}, allMean{ii}, trialNumber{ii}] = getMeanTraces(eyeTrialData, frames, frameLength, names, allCons, groupN(ii), subStart, subEnd);
end

%% Draw velocity trace plots
for ii = 1:length(groupN)
    % plot mean traces for each participant
    if individualPlots
        for subN = subStart:subEnd
            latencySI = round(latency(subN)+200);
                        % position trace
                        figure
                        for internalConN = 1:length(allCons.internalCons)
                            subplot(3, 1, internalConN)
                            hold on
                            for angleN = 1:length(allCons.apertureAngles)
                                if allCons.apertureAngles(angleN)>0 % upward internal dir
                                    lineStyle = '-';
                                else % downward internal dir
                                    lineStyle = '--';
                                end
                                p{angleN} = plot(indiMean{ii}.pos{angleN, internalConN}{1}(subN, :), indiMean{ii}.pos{angleN, internalConN}{2}(subN, :), ...
                                    'LineStyle', lineStyle, 'color', colorObjAngles(angleN, :)); %, 'LineWidth', 1);
            
                                % average pursuit onset
                                plot(indiMean{ii}.pos{angleN, internalConN}{1}(subN, latencySI), indiMean{ii}.pos{angleN, internalConN}{2}(subN, latencySI), '+k')
                                % average olp end
                                plot(indiMean{ii}.pos{angleN, internalConN}{1}(subN, latencySI+140), indiMean{ii}.pos{angleN, internalConN}{2}(subN, latencySI+140), '+b')
                                % average early/late clp phase border
                                plot(indiMean{ii}.pos{angleN, internalConN}{1}(subN, latencySI+140+round((900-latencySI-140)/2)), ...
                                    indiMean{ii}.pos{angleN, internalConN}{2}(subN, latencySI+140+round((900-latencySI-140)/2)), '+b')
                                % end of clp analysis window
                                plot(indiMean{ii}.pos{angleN, internalConN}{1}(subN, 900), indiMean{ii}.pos{angleN, internalConN}{2}(subN, 900), '+k')
                            end
                            title([names{subN}, ', ', internalConNames{internalConN}])
                            xlabel('Horizontal position (deg)')
                            ylabel('Vertical position (deg)')
                            %                 xlim([-200 800])
                            if internalConN==1
                                legend([p{:}], apertureAngleNames, 'Location', 'best')
                            end
                            box off
                        end
                        saveas(gcf, [eyeTracesFolder, 'individuals\posTrace_' names{subN} '.pdf'])
            
%             % position; difference from baseline
%             figure
%             for internalConN = 2:length(allCons.internalCons)
%                 subplot(2, 1, internalConN-1)
%                 hold on
%                 for angleN = 1:length(allCons.apertureAngles)
%                     if allCons.apertureAngles(angleN)>0 % upward internal dir
%                         lineStyle = '-';
%                     else % downward internal dir
%                         lineStyle = '--';
%                     end
%                     p{angleN} = plot(timePointsOnset, indiMean{ii}.posDiff{angleN, internalConN}{2}(subN, :), 'LineStyle', lineStyle, 'color', colorPlot(angleN, :)); %, 'LineWidth', 1);
%                     % average pursuit onset
%                     line([latency(subN) latency(subN)], [-0.5, 0.5])
%                     % average olp end
%                     line([latency(subN)+140 latency(subN)+140], [-0.5, 0.5])
%                     % average early/late clp phase border
%                     line([latency(subN)+140+(700-latency(subN)-140)/2 latency(subN)+140+(700-latency(subN)-140)/2], [-0.5, 0.5])
%                     % end of clp analysis window
%                     line([700 700], [-0.5, 0.5])
%                     
%                     %                     % average pursuit onset
%                     %                     plot(indiMean{ii}.posDiff{angleN, internalConN}{1}(subN, latencySI), indiMean{ii}.posDiff{angleN, internalConN}{2}(subN, latencySI), 'ok')
%                     %                     % average olp end
%                     %                     plot(indiMean{ii}.posDiff{angleN, internalConN}{1}(subN, latencySI+140), indiMean{ii}.posDiff{angleN, internalConN}{2}(subN, latencySI+140), 'ob')
%                     %                     % average early/late clp phase border
%                     %                     plot(indiMean{ii}.posDiff{angleN, internalConN}{1}(subN, latencySI+140+round((900-latency(subN)-140)/2)), ...
%                     %                         indiMean{ii}.posDiff{angleN, internalConN}{2}(subN, latencySI+140+round((900-latency(subN)-140)/2)), 'ob')
%                     %                     % end of clp analysis window
%                     %                     plot(indiMean{ii}.posDiff{angleN, internalConN}{1}(subN, 900), indiMean{ii}.posDiff{angleN, internalConN}{2}(subN, 900), 'ok')
%                 end
%                 title([names{subN}, ', ', internalConNames{internalConN}])
%                 xlabel('Time (ms)')
%                 ylabel('Vertical position diff (deg)')
%                 
%                 if internalConN==2
%                     legend([p{:}], apertureAngleNames, 'Location', 'best')
%                 end
%                 box off
%             end
%             saveas(gcf, [eyeTracesFolder, 'individuals\posTraceDiff_' names{subN} '.pdf'])
            
            %             % velocity; each internal cons is one figure; difference from
            %             % baseline
            %             for internalConN = 2:length(allCons.internalCons)
            %                 figure
            %                 for dimN = 1:2
            %                     subplot(2, 1, dimN)
            %                     hold on
            %                     for angleN = 1:length(allCons.apertureAngles)
            %                         if allCons.apertureAngles(angleN)>0 % upward internal dir
            %                             lineStyle = '-';
            %                         else % downward internal dir
            %                             lineStyle = '--';
            %                         end
            %                         p{angleN} = plot(timePointsOnset, indiMean{ii}.velDiff{angleN, internalConN}{dimN}(subN, :), 'LineStyle', lineStyle, 'color', colorPlot(angleN, :)); %, 'LineWidth', 1);
            %                     end
            %                     title([names{subN}, ', ', internalConNames{internalConN}])
            %                     xlabel('Time from RDK onset (ms)')
            %                     ylabel([dimNames{dimN}, ' eye velocity difference (deg/s)'])
            %                     xlim([-200 800])
            %
            %                     if dimN==1
            %                         legend([p{:}], apertureAngleNames, 'Location', 'best')
            %                     end
            %                     box off
            %                 end
            %                 saveas(gcf, [eyeTracesFolder, 'individuals\velTraceDiff_internalCon' num2str(allCons.internalCons(internalConN)) '_' names{subN} '.pdf'])
            %             end
            %             close all
        end
    end
    
    % plot mean traces of all participants in all probabilities
    if averagedPlots
        meanLatencyI = round(mean(latency))+200;
        % align at rdk onset
        
%         % position, averaged across all conditions, one trace from
%         % each--NOT GONNA WORK
%         % subgroup
%         load('subGroupListLPA.mat')
%         
%         % averaging in each observer, need to minus individual baseline
%         basePosAllX = []; % each row is one observer
%         biasPosAllX = []; 
%         basePosAllY = []; % each row is one observer
%         biasPosAllY = []; 
%         
%         for subN = 1:size(names, 2)
%             biasPosAllSubX = [];
%             biasPosAllSubY = [];
%             
%             for internalConN = 1:length(allCons.internalCons)
%                 if allCons.internalCons(internalConN)>0 % upward internal dir
%                     signI = 1;
%                 else % downward internal dir
%                     signI = -1;
%                 end
%                 
%                 subTempX = [];
%                 subTempY = [];
%                 for angleN = 1:length(allCons.apertureAngles)
%                     subTempX(angleN, :) = indiMean{1}.pos{angleN, internalConN}{1}(subN, :);
%                     subTempY(angleN, :) = indiMean{1}.pos{angleN, internalConN}{2}(subN, :);
%                 end
%                 
%                 if allCons.internalCons(internalConN)==0
%                     basePosAllX(subN, :) = nanmean(subTempX);
%                     basePosAllY(subN, :) = nanmean(subTempY);
%                 else
%                     biasPosAllSubX = [biasPosAllSubX; nanmean(subTempX)];
%                     biasPosAllSubY = [biasPosAllSubY; signI.*nanmean(subTempY)];
%                 end
%             end
%             biasPosAllX(subN, :) = nanmean(biasPosAllSubX);
%             biasPosAllY(subN, :) = nanmean(biasPosAllSubY);
%         end
%         
%         % averaging in each subgroup
%         biasPosX = []; % first rwo is assimilation group, second row is contrast group
%         biasPosY = [];        
%         basePosX = []; % first rwo is assimilation group, second row is contrast group
%         basePosY = [];
%         
%         for subGroupN = 1:size(subGroup, 2)
%             idx = subGroup{subGroupN};
%             basePosX(subGroupN, :) = nanmean(basePosAllX(idx, :));
%             basePosY(subGroupN, :) = nanmean(basePosAllY(idx, :));
%             
%             biasPosX(subGroupN, :) = nanmean(biasPosAllX(idx, :));
%             biasPosY(subGroupN, :) = nanmean(biasPosAllY(idx, :));
%         end
%         
%         figure
%         hold on
%         for subGroupN = 1:size(subGroup, 2)
%             for conN = 1:2
%                 if conN==1
%                     xPlot = basePosX(subGroupN, :);
%                     yPlot = basePosY(subGroupN, :);
%                     lineStyle = '--';
%                 else
%                     xPlot = biasPosX(subGroupN, :);
%                     yPlot = biasPosY(subGroupN, :);
%                     lineStyle = '-';
%                 end
%             
%             p{subGroupN} = plot(xPlot, yPlot, 'lineStyle', lineStyle, 'color', colorGroup(subGroupN, :)); %, 'LineWidth', 1);
%             
%             % average pursuit onset
%             plot(xPlot(meanLatencyI), yPlot(meanLatencyI), '+k')
%             % average olp end
%             plot(xPlot(meanLatencyI+140), yPlot(meanLatencyI+140), '+k')
%             % average early/late clp phase border
%             plot(xPlot(meanLatencyI+140+round((900-meanLatencyI-140)/2)), ...
%                 yPlot(meanLatencyI+140+round((900-meanLatencyI-140)/2)), '+k')
%             % end of clp analysis window
%             plot(xPlot(900), yPlot(900), '+k')
%         end
%         xlim([0, 6])
%         ylim([-1, 1])
%         xlabel('Horizontal eye position (deg)')
%         ylabel(['Vertical eye\n position (deg)'])
%     
%         legend([p{:}], {'assimilation-base', 'assimilation-bias', 'contrast-base', 'contrast-bias'}, 'Location', 'best')
%         box off
%         saveas(gcf, [eyeTracesFolder, 'posTrace_subgroup.pdf'])
        
        % position, all object motion conditions
        
        for internalConN = 1:length(allCons.internalCons)
            %             subplot(3, 1, internalConN)
            figure
            hold on
            for angleN = 1:length(allCons.apertureAngles)
                if allCons.apertureAngles(angleN)>0 % upward internal dir
                    lineStyle = '-';
                else % downward internal dir
                    lineStyle = '--';
                end
                p{angleN} = plot(allMean{ii}.pos{angleN, internalConN}{1}, allMean{ii}.pos{angleN, internalConN}{2}, 'LineStyle', lineStyle, 'color', colorObjAngles(angleN, :)); %, 'LineWidth', 1);
                
                % average pursuit onset
                plot(allMean{ii}.pos{angleN, internalConN}{1}(meanLatencyI(internalConN)), allMean{ii}.pos{angleN, internalConN}{2}(meanLatencyI(internalConN)), '+k', 'MarkerSize', 15)
                % average olp end
                plot(allMean{ii}.pos{angleN, internalConN}{1}(meanLatencyI(internalConN)+140), allMean{ii}.pos{angleN, internalConN}{2}(meanLatencyI(internalConN)+140), '+k', 'MarkerSize', 15)
                % average early/late clp phase border
%                 plot(allMean{ii}.pos{angleN, internalConN}{1}(meanLatencyI+140+round((900-meanLatencyI-140)/2)), ...
%                     allMean{ii}.pos{angleN, internalConN}{2}(meanLatencyI+140+round((900-meanLatencyI-140)/2)), '+k', 'MarkerSize', 5)
                % end of clp analysis window
                plot(allMean{ii}.pos{angleN, internalConN}{1}(900), allMean{ii}.pos{angleN, internalConN}{2}(900), '+k', 'MarkerSize', 15)
            end
            
            xlim([0, 6])
            ylim([-1, 1])
%             axis square
            legend([p{:}], apertureAngleNames, 'Location', 'best', 'box', 'off')
%             title(['all, ', internalConNames{internalConN}])
            xlabel('Horizontal eye position (deg)')
            ylabel(['Vertical eye\n position (deg)'])
            saveas(gcf, [eyeTracesFolder, 'posTrace_all_', internalConNames{internalConN}, '.pdf'])
        end
        
%         if internalConN==1
%             legend([p{:}], apertureAngleNames, 'Location', 'best', 'box', 'off')
%         end
        
%         % position difference from baseline
%         figure
%         for internalConN = 2:length(allCons.internalCons)
%             subplot(2, 1, internalConN-1)
%             hold on
%             for angleN = 1:length(allCons.apertureAngles)
%                 if allCons.apertureAngles(angleN)>0 % upward internal dir
%                     lineStyle = '-';
%                 else % downward internal dir
%                     lineStyle = '--';
%                 end
%                 p{angleN} = plot(timePointsOnset, allMean{ii}.posDiff{angleN, internalConN}{2}, 'LineStyle', lineStyle, 'color', colorPlot(angleN, :)); %, 'LineWidth', 1);
%                 % average pursuit onset
%                 line([mean(latency) mean(latency)], [-0.5, 0.5])
%                 % average olp end
%                 line([mean(latency)+140 mean(latency)+140], [-0.5, 0.5])
%                 % average early/late clp phase border
%                 line([mean(latency)+140+(700-mean(latency)-140)/2 mean(latency)+140+(700-mean(latency)-140)/2], [-0.5, 0.5])
%                 % end of clp analysis window
%                 line([700 700], [-0.5, 0.5])
%                 
%                 %                 % average pursuit onset
%                 %                 plot(allMean{ii}.posDiff{angleN, internalConN}{1}(subN, meanLatencyI), allMean{ii}.posDiff{angleN, internalConN}{2}(subN, meanLatencyI), 'ok')
%                 %                 % average olp end
%                 %                 plot(allMean{ii}.posDiff{angleN, internalConN}{1}(subN, meanLatencyI+140), allMean{ii}.posDiff{angleN, internalConN}{2}(subN, meanLatencyI+140), 'ob')
%                 %                 % average early/late clp phase border
%                 %                 plot(allMean{ii}.posDiff{angleN, internalConN}{1}(subN, meanLatencyI+140+round((900-meanLatencyI-140)/2)), ...
%                 %                     allMean{ii}.posDiff{angleN, internalConN}{2}(subN, meanLatencyI+140+round((900-meanLatencyI-140)/2)), 'ob')
%                 %                 % end of clp analysis window
%                 %                 plot(allMean{ii}.posDiff{angleN, internalConN}{1}(subN, 900), allMean{ii}.posDiff{angleN, internalConN}{2}(subN, 900), 'ok')
%                 title(['all, ', internalConNames{internalConN}])
%                 xlabel('Horizontal position (deg)')
%                 ylabel(['Vertical position (deg)'])
%             end
%             
%             if internalConN==2
%                 legend([p{:}], apertureAngleNames, 'Location', 'best')
%             end
%             box off
%             saveas(gcf, [eyeTracesFolder, 'posTraceDiff_all.pdf'])
%         end
        
%         % velocity, all conditions
%         for internalConN = 1:length(allCons.internalCons)
%             figure
%             for dimN = 1:2
%                 subplot(2, 1, dimN)
%                 hold on
%                 for angleN = 1:length(allCons.apertureAngles)
%                     if allCons.apertureAngles(angleN)>0 % upward internal dir
%                         lineStyle = '-';
%                     else % downward internal dir
%                         lineStyle = '--';
%                     end
%                     p{angleN} = plot(timePointsOnset, allMean{ii}.vel{angleN, internalConN}{dimN}, 'LineStyle', lineStyle, 'color', colorPlot(angleN, :)); %, 'LineWidth', 1);
%                 end
%                 % average pursuit onset
%                 line([mean(latency) mean(latency)], [-0.5, 0.5])
%                 % average olp end
%                 line([mean(latency)+140 mean(latency)+140], [-0.5, 0.5])
%                 % average early/late clp phase border
%                 line([mean(latency)+140+(700-mean(latency)-140)/2 mean(latency)+140+(700-mean(latency)-140)/2], [-0.5, 0.5])
%                 % end of clp analysis window
%                 line([700 700], [-0.5, 0.5])
%                 
%                 title(['all, ', internalConNames{internalConN}])
%                 xlabel('Time from RDK onset (ms)')
%                 ylabel([dimNames{dimN}, ' eye velocity (deg/s)'])
%                 xlim([-200 800])
%                 
%                 if dimN==1
%                     legend([p{:}], apertureAngleNames, 'Location', 'best')
%                 end
%                 box off
%             end
%             saveas(gcf, [eyeTracesFolder, 'velTrace_internalCon' num2str(allCons.internalCons(internalConN)) '_all.pdf'])
%         end
%         
%         % difference from baseline
%         for internalConN = 2:length(allCons.internalCons)
%             figure
%             for dimN = 1:2
%                 subplot(2, 1, dimN)
%                 hold on
%                 for angleN = 1:length(allCons.apertureAngles)
%                     if allCons.apertureAngles(angleN)>0 % upward internal dir
%                         lineStyle = '-';
%                     else % downward internal dir
%                         lineStyle = '--';
%                     end
%                     p{angleN} = plot(timePointsOnset, allMean{ii}.velDiff{angleN, internalConN}{dimN}, 'LineStyle', lineStyle, 'color', colorPlot(angleN, :)); %, 'LineWidth', 1);
%                 end
%                 % average pursuit onset
%                 line([mean(latency) mean(latency)], [-0.5, 0.5])
%                 % average olp end
%                 line([mean(latency)+140 mean(latency)+140], [-0.5, 0.5])
%                 % average early/late clp phase border
%                 line([mean(latency)+140+(700-mean(latency)-140)/2 mean(latency)+140+(700-mean(latency)-140)/2], [-0.5, 0.5])
%                 % end of clp analysis window
%                 line([700 700], [-0.5, 0.5])
%                 
%                 title(['all, ', internalConNames{internalConN}])
%                 xlabel('Time from RDK onset (ms)')
%                 ylabel([dimNames{dimN}, ' eye velocity difference (deg/s)'])
%                 xlim([-200 800])
%                 
%                 if dimN==1
%                     legend([p{:}], apertureAngleNames, 'Location', 'best')
%                 end
%                 box off
%             end
%             saveas(gcf, [eyeTracesFolder, 'velTraceDiff_internalCon' num2str(allCons.internalCons(internalConN)) '_all.pdf'])
%         end
    end
end
% save('eyeTracesAll.mat', 'indiMean', 'allMean')

%% generate csv files, each file for one probability condition
% % each row is the mean velocity trace of one participant
% % use the min frame length--the lengeth where all participants have
% % valid data points
% % cd(RsaveFolder)
% % averaged traces
% for ii = 1:length(groupN)
%     % raw velocity values
%     for internalConN = 1:length(allCons.internalCons)
%         velTAverageSub = [];
%         for binN = 1:2
%             if binN==1
%                 dataTemp = indiMean{ii}{probNmerged}.left(:, (maxFrameLength-minFrameLength+1):end);
%             else
%                 dataTemp = indiMean{ii}{probNmerged}.right(:, (maxFrameLength-minFrameLength+1):end);
%             end
%             for subN = 1:size(names, 2)
%                 velTAverageSub((binN-1)*length(names)+subN, :) = dataTemp(subN, :);
%             end
%         end
%         csvwrite(['velocityTrace_' groupName{groupN(ii)} '_internalCon' num2str(), '.csv'], velTAverageSub)
%     end
% end

%%
function [indiMean, allMean, trialNumber] = getMeanTraces(eyeTrialData, frames, frameLength, names, allCons, groupN, subStart, subEnd)
% calculate mean traces
% indiMean: each row is one participant
% allMean: averaged across participants
% trialNumber: corresponds to indiMean, the trial number for each element

minFrameLength = nanmin(frameLength);
for internalConN = 1:length(allCons.internalCons)
    if allCons.internalCons(internalConN)==0
        rdkCoh = 0;
        rdkInternalDir = 0;
    else
        rdkCoh = 1;
        rdkInternalDir = allCons.internalCons(internalConN);
    end
    
    for angleN = 1:length(allCons.apertureAngles)
        
        indiMean.vel{angleN, internalConN}{1} = NaN(length(names), minFrameLength); % horizontal
        indiMean.vel{angleN, internalConN}{2} = NaN(length(names), minFrameLength); % vertical
        
        indiMean.pos{angleN, internalConN}{1} = NaN(length(names), minFrameLength); % horizontal
        indiMean.pos{angleN, internalConN}{2} = NaN(length(names), minFrameLength); % vertical
        
        for subN = subStart:subEnd
            idxT = find(eyeTrialData.errorStatus(subN, :)==0 & ...
                eyeTrialData.rdkCoh(subN, :)==rdkCoh & ...
                eyeTrialData.rdkInternalDir(subN, :)==rdkInternalDir);
            
            switch groupN
                case 1 % trials by visual motion
                    idx = find(eyeTrialData.rdkApertureAngle(subN, idxT)==allCons.apertureAngles(angleN));
                    
                    %                     leftIdx = find(eyeTrialData.rdkApertureDir(subN, idxT)==180);
                    %                     rightIdx = find(eyeTrialData.rdkApertureDir(subN, idxT)==0);
                    
                    %                 upIdx = find(eyeTrialData.rdkInternalDir(subN, idxT)>0 & eyeTrialData.rdkCoh(subN, idxT)==allCons(conN, 2));
                    %                 downIdx = find(eyeTrialData.rdkInternalDir(subN, idxT)<0 & eyeTrialData.rdkCoh(subN, idxT)==allCons(conN, 2));
                    %             case 2 % trials by perceived motion
                    %                 upIdx = find(eyeTrialData.choice(subN, idxT)>0 & eyeTrialData.rdkDirSD(subN, idxT)==sdValue);
                    %                 downIdx = find(eyeTrialData.choice(subN, idxT)<0 & eyeTrialData.rdkDirSD(subN, idxT)==sdValue);
            end
            
            % individual mean traces
            indiMean.vel{angleN, internalConN}{1}(subN, :) = ...
                nanmean(frames.onset{subN, internalConN}.velX(idx, 1:minFrameLength), 1);
            indiMean.vel{angleN, internalConN}{2}(subN, :) = ...
                nanmean(frames.onset{subN, internalConN}.velY(idx, 1:minFrameLength), 1);
            indiMean.pos{angleN, internalConN}{1}(subN, :) = nanmean(frames.onset{subN, internalConN}.posX(idx, 1:minFrameLength), 1);
            indiMean.pos{angleN, internalConN}{2}(subN, :) = nanmean(frames.onset{subN, internalConN}.posY(idx, 1:minFrameLength), 1);
            
            if internalConN>=2
                indiMean.velDiff{angleN, internalConN}{1}(subN, :) = ...
                    indiMean.vel{angleN, internalConN}{1}(subN, :)-indiMean.vel{angleN, 1}{1}(subN, :);
                indiMean.velDiff{angleN, internalConN}{2}(subN, :) = ...
                    indiMean.vel{angleN, internalConN}{2}(subN, :)-indiMean.vel{angleN, 1}{2}(subN, :);
                indiMean.posDiff{angleN, internalConN}{1}(subN, :) = ...
                    indiMean.pos{angleN, internalConN}{1}(subN, :)-indiMean.pos{angleN, 1}{1}(subN, :);
                indiMean.posDiff{angleN, internalConN}{2}(subN, :) = ...
                    indiMean.pos{angleN, internalConN}{2}(subN, :)-indiMean.pos{angleN, 1}{2}(subN, :);
            end
            
            %             indiMean.pos{angleN, internalConN}{1}(subN, :) = nanmean([-frames.onset{subN, internalConN}.posX(leftIdx, 1:minFrameLength); frames.onset{subN, internalConN}.posX(rightIdx, 1:minFrameLength)], 1);
            %             indiMean.pos{angleN, internalConN}{2}(subN, :) = nanmean([frames.onset{subN, internalConN}.posY(leftIdx, 1:minFrameLength); frames.onset{subN, internalConN}.posY(rightIdx, 1:minFrameLength)], 1);
            %
            trialNumber{angleN, internalConN}(subN, 1) = length(idx);
        end
        
        % collapsed all participants
        allMean.vel{angleN, internalConN}{1} = nanmean(indiMean.vel{angleN, internalConN}{1}, 1);
        allMean.vel{angleN, internalConN}{2} = nanmean(indiMean.vel{angleN, internalConN}{2}, 1);
        allMean.pos{angleN, internalConN}{1} = nanmean(indiMean.pos{angleN, internalConN}{1}, 1);
        allMean.pos{angleN, internalConN}{2} = nanmean(indiMean.pos{angleN, internalConN}{2}, 1);
        if internalConN>=2 % hmmmm... shouldn't this be the mean of the individual differences?
            allMean.velDiff{angleN, internalConN}{1} = nanmean(indiMean.velDiff{angleN, internalConN}{1}, 1);
            allMean.velDiff{angleN, internalConN}{2} = nanmean(indiMean.velDiff{angleN, internalConN}{2}, 1);
            allMean.posDiff{angleN, internalConN}{1} = nanmean(indiMean.posDiff{angleN, internalConN}{1}, 1);
            allMean.posDiff{angleN, internalConN}{2} = nanmean(indiMean.posDiff{angleN, internalConN}{2}, 1);
%             allMean.velDiff{angleN, internalConN}{1} = ...
%                 allMean.vel{angleN, internalConN}{1}-allMean.vel{angleN, 1}{1};
%             allMean.velDiff{angleN, internalConN}{2} = ...
%                 allMean.vel{angleN, internalConN}{2}-allMean.vel{angleN, 1}{2};
%             allMean.posDiff{angleN, internalConN}{1} = ...
%                 allMean.pos{angleN, internalConN}{1}-allMean.pos{angleN, 1}{1};
%             allMean.posDiff{angleN, internalConN}{2} = ...
%                 allMean.pos{angleN, internalConN}{2}-allMean.pos{angleN, 1}{2};
        end
    end
end
end