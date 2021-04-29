% plot position&velocity traces for translating RDK
% also generate csv file for further plotting in R
initializeParas;

% choose which plot to look at now
individualPlots = 1;
averagedPlots = 0;
textFontSize = 8;
subStart = 3;
subEnd = 3;

groupName = {'visualDir'};
% naming by trial type (could include grouping rules) + group based on which direction (visual or perceived)
groupN = [1]; % corresponds to the listed rules... can choose multiple, just list as a vector
% when choosing multiple groupN, will plot each group rule in one figure

% prepare for the grouping...
allCons = internalDirCons'; % first column is internal dir
conNames = internalDirNames;

dimNames = {'Horizontal', 'Vertical'};

%% align rdk onset, frame data for all trials
frameLength = NaN(size(names, 2), 1);
for subN = subStart:subEnd
    cd(analysisFolder)
    load(['eyeTrialDataSub_' names{subN} '.mat']);

    idxT = find(eyeTrialData.errorStatus(subN, :)==0);
    tempL = eyeTrialData.frameLog.rdkOff(subN, idxT)-eyeTrialData.frameLog.rdkOn(subN, idxT);
    tempL(tempL==0) = [];
    frameLength(subN, 1) = min(tempL);
    
    for internalDirN = 1:length(allCons)
        idxT = find(eyeTrialData.errorStatus(subN, :)==0 & eyeTrialData.rdkInternalDir(subN, :)==internalDirCons(internalDirN));

        lengthT = length(idxT);
%         frames.onset{subN, internalDirN}.posX = NaN(lengthT, frameLength(subN, 1));
%         frames.onset{subN, internalDirN}.posY = NaN(lengthT, frameLength(subN, 1));
        frames.onset{subN, internalDirN}.velX = NaN(lengthT, frameLength(subN, 1));
        frames.onset{subN, internalDirN}.velY = NaN(lengthT, frameLength(subN, 1));
        
        for trialN = 1:lengthT
            % align at target onset
            startI = eyeTrialData.frameLog.rdkOn(subN, idxT(trialN));
            endI = startI+frameLength(subN, 1)-1;
%             frames.onset{subN, internalDirN}.posX(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.X_interpolSac(startI:endI);
%             frames.onset{subN, internalDirN}.posY(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.Y_interpolSac(startI:endI);
            frames.onset{subN, internalDirN}.velX(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.DX_interpolSac(startI:endI);
            frames.onset{subN, internalDirN}.velY(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.DY_interpolSac(startI:endI);
        end
    end

end
maxFrameLength = max(frameLength);

% plotting parameters
minFrameLength = min(frameLength);
framePerSec = 1/sampleRate;
timePointsOnset = [(1:minFrameLength)-1]*framePerSec*1000; % rdk onset is 0

%% calculate mean traces
for ii = 1:length(groupN)
    [indiMean{ii}, allMean{ii}, trialNumber{ii}] = getMeanTraces(eyeTrialData, frames, frameLength, names, allCons, groupN(ii), subStart, subEnd);
end

%% Draw velocity trace plots
for ii = 1:length(groupN)
    % plot mean traces for each participant
    if individualPlots
        for subN = subStart:subEnd
            % velocity
            figure
            for dimN = 1:2
                subplot(2, 1, dimN)
                hold on
                for conN = 1:size(allCons, 1)
                    if allCons(conN, 1)>0 % upward internal dir
                        lineStyle = '-';
                    else % downward internal dir
                        lineStyle = '--';
                    end
                    p{conN} = plot(timePointsOnset, indiMean{ii}.vel{conN, dimN}(subN, :), 'LineStyle', lineStyle, 'color', colorPlot(mod(conN, 3)+1, :)); %, 'LineWidth', 1);
                end
                if dimN==1
                legend([p{:}], conNames, 'Location', 'best')
                end
                title(names{subN})
                xlabel('Time from RDK onset (ms)')
                ylabel([dimNames{dimN}, ' eye velocity (deg/s)'])
                %             xlim([-500 700])
%                             ylim([-0.5, 2.5])
                box off
            end
            saveas(gcf, [eyeTracesFolder, 'velTrace_' groupName{groupN(ii)} '_coh1_' names{subN} '.pdf'])
        end
    end
    
%     % plot mean traces of all participants in all probabilities
%     if averagedPlots
%         % align at rdk onset
%         figure
%         for dimN = 1:2
%             subplot(2, 1, dimN)
%             hold on
%             for conN = 1:size(allCons, 1)
%                 internalDirN = allCons(conN, 1)+1;
%                 if allCons(conN, 2)==30 % low sd, solid line
%                     lineStyle = '-';
%                 else % high sd
%                     lineStyle = '--';
%                 end
%                 p{conN} = plot(timePointsOnset{internalDirN}, allMean{ii}.onset{conN, dimN}, 'LineStyle', lineStyle, 'color', colorPlot(conN, :)); %, 'LineWidth', 1);
%             end
%             legend([p{:}], conNames, 'Location', 'best')
%             title(['up&down merged, vel in the ', groupName{groupN(ii)}])
%             xlabel('Time from RDK onset (ms)')
%             ylabel([dimNames{dimN}, ' eye velocity (deg/s)'])
%             %             xlim([-500 700])
%             %             ylim(yRange)
%             box off
%         end
%         saveas(gcf, [eyeTracesFolder, '\velTrace_' groupName{groupN(ii)} '_onsetAlign_all.pdf'])
%     end
end

%% generate csv files, each file for one probability condition
% % each row is the mean velocity trace of one participant
% % use the min frame length--the lengeth where all participants have
% % valid data points
% % cd(RsaveFolder)
% % averaged traces
% for ii = 1:length(groupN)
%     for probNmerged = 1:probTotalN
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
%         csvwrite(['velocityTrace_' groupName{groupN(ii)} '_exp' num2str(expN) '_prob' num2str(probCons(probNmerged+probTotalN-1)), '.csv'], velTAverageSub)
%     end
% end

%%
function [indiMean, allMean, trialNumber] = getMeanTraces(eyeTrialData, frames, frameLength, names, allCons, groupN, subStart, subEnd)
% calculate mean traces
% indiMean: each row is one participant
% allMean: averaged across participants
% trialNumber: corresponds to indiMean, the trial number for each element

minFrameLength = nanmin(frameLength);
for conN = 1:size(allCons, 1) 
    % initialize
    internalDirN = conN;

    indiMean.vel{conN, 1} = NaN(length(names), minFrameLength); % horizontal
    indiMean.vel{conN, 2} = NaN(length(names), minFrameLength); % vertical
    
%     indiMean.pos{conN, 1} = NaN(length(names), minFrameLength); % horizontal
%     indiMean.pos{conN, 2} = NaN(length(names), minFrameLength); % vertical
    
    for subN = subStart:subEnd
        idxT = find(eyeTrialData.errorStatus(subN, :)==0 & eyeTrialData.rdkInternalDir(subN, :)==allCons(conN, 1));
        
        switch groupN
            case 1 % trials by visual motion
                leftIdx = find(eyeTrialData.rdkApertureDir(subN, idxT)==180);
                rightIdx = find(eyeTrialData.rdkApertureDir(subN, idxT)==0);
                
                %                 upIdx = find(eyeTrialData.rdkInternalDir(subN, idxT)>0 & eyeTrialData.rdkCoh(subN, idxT)==allCons(conN, 2));
                %                 downIdx = find(eyeTrialData.rdkInternalDir(subN, idxT)<0 & eyeTrialData.rdkCoh(subN, idxT)==allCons(conN, 2));
                %             case 2 % trials by perceived motion
                %                 upIdx = find(eyeTrialData.choice(subN, idxT)>0 & eyeTrialData.rdkDirSD(subN, idxT)==sdValue);
                %                 downIdx = find(eyeTrialData.choice(subN, idxT)<0 & eyeTrialData.rdkDirSD(subN, idxT)==sdValue);
        end
        
        % individual mean traces
        indiMean.vel{conN, 1}(subN, :) = nanmean([-frames.onset{subN, internalDirN}.velX(leftIdx, 1:minFrameLength); frames.onset{subN, internalDirN}.velX(rightIdx, 1:minFrameLength)], 1);
        indiMean.vel{conN, 2}(subN, :) = sign(allCons(conN, 1))*nanmean([frames.onset{subN, internalDirN}.velY(leftIdx, 1:minFrameLength); frames.onset{subN, internalDirN}.velY(rightIdx, 1:minFrameLength)], 1);
%         indiMean.pos{conN, 1}(subN, :) = nanmean([-frames.onset{subN, internalDirN}.posX(leftIdx, 1:minFrameLength); frames.onset{subN, internalDirN}.posX(rightIdx, 1:minFrameLength)], 1);
%         indiMean.pos{conN, 2}(subN, :) = nanmean([frames.onset{subN, internalDirN}.posY(leftIdx, 1:minFrameLength); frames.onset{subN, internalDirN}.posY(rightIdx, 1:minFrameLength)], 1);
        
        trialNumber{conN}(subN, 1) = length(leftIdx)+length(rightIdx);
    end
    
    % collapsed all participants
    allMean.vel{conN, 1} = nanmean(indiMean.vel{conN, 1}, 1);
    allMean.vel{conN, 2} = nanmean(indiMean.vel{conN, 2}, 1);
end
end