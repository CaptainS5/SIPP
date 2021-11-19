% plot position&velocity traces for translating RDK
% also generate csv file for further plotting in R

% currently each figure is one internal motion condition; plot the difference in
% velocity traces between the exp conditions and the baseline (static)
% condition
initializeParas;

% choose which plot to look at now
individualPlots = 1;
averagedPlots = 0;
subStart = 1;
subEnd = 1;

groupName = {'visualDir'};
% naming by trial type (could include grouping rules) + group based on which direction (visual or perceived)
groupN = [1]; % corresponds to the listed rules... can choose multiple, just list as a vector
% when choosing multiple groupN, will plot each group rule in one figure

% prepare for the grouping...
allCons.perturbPhase = perturbPhase;
allCons.internalCons = internalCons; 
allCons.apertureAngles = apertureAngles;

% plotting settings
dimNames = {'Horizontal', 'Vertical'};
yRangeX = [6 12];
yRangeY = [-3 3];

%% align Perturbation onset, frame data for all trials
frameLength = NaN(size(names, 2), 1);
latency = NaN(size(names));
for subN = subStart:subEnd
    cd(analysisFolder)
    load(['eyeTrialDataSub_' names{subN} '.mat']);
    
    idxT = find(eyeTrialData.errorStatus(subN, :)==0);
    tempL = eyeTrialData.frameLog.rdkOff(subN, idxT)-eyeTrialData.frameLog.perturbationOn(subN, idxT)+100; % starting from 100 ms before perturbation onset, to the end of RDK
    tempL(tempL==0) = [];
    frameLength(subN, 1) = min(tempL);
    %     % latency
    %     latency(subN) = nanmean(eyeTrialData.pursuit.onset(subN, idxT)-eyeTrialData.frameLog.rdkOn(subN, idxT));
    
    for perturbN = 1:length(perturbPhase)
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
                eyeTrialData.rdkInternalDir(subN, :)==rdkInternalDir & ...
                eyeTrialData.perturbPhase(subN, :)==perturbPhase(perturbN));
            
            lengthT = length(idxT);
            %         frames.onset{subN, internalDirN}.posX = NaN(lengthT, frameLength(subN, 1));
            %         frames.onset{subN, internalDirN}.posY = NaN(lengthT, frameLength(subN, 1));
            frames.onset{subN}{perturbN, internalConN}.velX = NaN(lengthT, frameLength(subN, 1));
            frames.onset{subN}{perturbN, internalConN}.velY = NaN(lengthT, frameLength(subN, 1));
            
            for trialN = 1:lengthT
                % align at 100 ms before perturbation onset
                startI = eyeTrialData.frameLog.perturbationOn(subN, idxT(trialN))-100+1;
                endI = startI+frameLength(subN, 1)-1;
                %             frames.onset{subN, internalDirN}.posX(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.X_interpolSac(startI:endI);
                %             frames.onset{subN, internalDirN}.posY(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.Y_interpolSac(startI:endI);
                frames.onset{subN}{perturbN, internalConN}.velX(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.DX_interpolSac(startI:endI);
                frames.onset{subN}{perturbN, internalConN}.velY(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.DY_interpolSac(startI:endI);
            end
            %         end
        end
    end
end
maxFrameLength = max(frameLength);

% plotting parameters
minFrameLength = min(frameLength);
framePerSec = 1/sampleRate;
timePointsOnset = [(-99:1:(minFrameLength-100))]*framePerSec*1000; % perturbation onset is 0

%% calculate mean traces
for ii = 1:length(groupN)
    [indiMean{ii}, allMean{ii}, trialNumber{ii}] = getMeanTraces(eyeTrialData, frames, frameLength, names, allCons, groupN(ii), subStart, subEnd);
end

%% Draw velocity trace plots
for ii = 1:length(groupN)
    % plot mean traces for each participant
    if individualPlots
        for subN = subStart:subEnd
            
            for perturbN = 1:length(allCons.perturbPhase)
                figure
                for internalConN = 1:length(allCons.internalCons)
                    for dimN = 1:2
                        subplot(2, 2, (internalConN-1)*2+dimN)
                        hold on
                        for angleN = 1:length(allCons.apertureAngles)
                            if allCons.apertureAngles(angleN)>0 % upward internal dir
                                lineStyle = '-';
                            else % downward internal dir
                                lineStyle = '--';
                            end
                            p{angleN} = plot(timePointsOnset, indiMean{ii}.vel{perturbN}{angleN, internalConN}{dimN}(subN, :), 'LineStyle', lineStyle, 'color', colorPlot(angleN, :)); %, 'LineWidth', 1);
                        end
                        title([names{subN}, ', ', internalConNames{internalConN}])
                        xlabel('Time from perturbation onset (ms)')
                        ylabel([dimNames{dimN}, ' eye velocity (deg/s)'])
                        xlim([-100 300])
                        
                        if internalConN==1 && dimN==1
                            legend([p{:}], apertureAngleNames, 'Location', 'best')
                        end
                        box off
                    end
                    saveas(gcf, [eyeTracesFolder, 'individuals\exp2\velTraceDiff_perturbPhase' num2str(allCons.perturbPhase(perturbN)) '_' names{subN} '.pdf'])
                end
            end
            
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
        end
    end
    
%     % plot mean traces of all participants in all probabilities
%     if averagedPlots
%         % align at rdk onset
%         % all conditions
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
%     end
end

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
for perturbN = 1:length(allCons.perturbPhase)
    for internalConN = 1:length(allCons.internalCons)
        if allCons.internalCons(internalConN)==0
            rdkCoh = 0;
            rdkInternalDir = 0;
        else
            rdkCoh = 1;
            rdkInternalDir = allCons.internalCons(internalConN);
        end
        
        for angleN = 1:length(allCons.apertureAngles)
            
            indiMean.vel{perturbN}{angleN, internalConN}{1} = NaN(length(names), minFrameLength); % horizontal
            indiMean.vel{perturbN}{angleN, internalConN}{2} = NaN(length(names), minFrameLength); % vertical
            
            %     indiMean.pos{conN, 1} = NaN(length(names), minFrameLength); % horizontal
            %     indiMean.pos{conN, 2} = NaN(length(names), minFrameLength); % vertical
            
            for subN = subStart:subEnd
                idxT = find(eyeTrialData.errorStatus(subN, :)==0 & ...
                    eyeTrialData.rdkCoh(subN, :)==rdkCoh & ...
                    eyeTrialData.rdkInternalDir(subN, :)==rdkInternalDir & ...
                    eyeTrialData.perturbPhase(subN, :)==allCons.perturbPhase(perturbN));
                
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
                indiMean.vel{perturbN}{angleN, internalConN}{1}(subN, :) = ...
                    nanmean(frames.onset{subN}{perturbN, internalConN}.velX(idx, 1:minFrameLength), 1);
                indiMean.vel{perturbN}{angleN, internalConN}{2}(subN, :) = ...
                    nanmean(frames.onset{subN}{perturbN, internalConN}.velY(idx, 1:minFrameLength), 1);
%                 if internalConN>=2
%                     indiMean.velDiff{angleN, internalConN}{1}(subN, :) = ...
%                         indiMean.vel{angleN, internalConN}{1}(subN, :)-indiMean.vel{angleN, 1}{1}(subN, :);
%                     indiMean.velDiff{angleN, internalConN}{2}(subN, :) = ...
%                         indiMean.vel{angleN, internalConN}{2}(subN, :)-indiMean.vel{angleN, 1}{2}(subN, :);
%                 end
                %         indiMean.pos{conN, 1}(subN, :) = nanmean([-frames.onset{subN, internalDirN}.posX(leftIdx, 1:minFrameLength); frames.onset{subN, internalDirN}.posX(rightIdx, 1:minFrameLength)], 1);
                %         indiMean.pos{conN, 2}(subN, :) = nanmean([frames.onset{subN, internalDirN}.posY(leftIdx, 1:minFrameLength); frames.onset{subN, internalDirN}.posY(rightIdx, 1:minFrameLength)], 1);
                
                trialNumber{perturbN}{angleN, internalConN}(subN, 1) = length(idx);
            end
            
            % collapsed all participants
            allMean.vel{perturbN}{angleN, internalConN}{1} = nanmean(indiMean.vel{perturbN}{angleN, internalConN}{1}, 1);
            allMean.vel{perturbN}{angleN, internalConN}{2} = nanmean(indiMean.vel{perturbN}{angleN, internalConN}{2}, 1);
%             if internalConN>=2
%                 allMean.velDiff{angleN, internalConN}{1} = ...
%                     allMean.vel{angleN, internalConN}{1}-allMean.vel{angleN, 1}{1};
%                 allMean.velDiff{angleN, internalConN}{2} = ...
%                     allMean.vel{angleN, internalConN}{2}-allMean.vel{angleN, 1}{2};
%             end
        end
    end
end
end