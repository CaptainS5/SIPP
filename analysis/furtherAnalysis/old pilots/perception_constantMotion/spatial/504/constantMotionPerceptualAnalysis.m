% plot perceptual results for the constant internal motion exp
initializeParas;

% choose which plot to look at now
individualPlots = 1;
averagedPlots = 0;

% plot settings
textFontSize = 8;
plotSub = {'504'};

% make rdkApertureAngle left/right flipped--only the up/down values
% relative to the aperture direction
idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkApertureDir==180); % leftward valid trials
eyeTrialData.rdkApertureAngle(idxT) = 180-eyeTrialData.rdkApertureAngle(idxT);

% % flip the down aperture directions
% idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkApertureAngle<0); 
% eyeTrialData.response(idxT) = -eyeTrialData.response(idxT);

%%
% initialization
summaryData = table;
count = 1;
for subPlotN = 1:size(plotSub, 2)
    subN = find(strcmp(names, plotSub{subPlotN}));
    for internalConN = 1:size(internalCons, 2) % each column is one internal condition
        summaryData.sub(count, 1) = subPlotN;
        if internalCons(internalConN)==0
            summaryData.rdkCoh(count, 1) = 0;
            summaryData.rdkInternalDir(count, 1) = 0;
        else
            summaryData.rdkCoh(count, 1) = 1;
            summaryData.rdkInternalDir(count, 1) = internalCons(internalConN);
        end
        summaryData.rdkInternalSpeed(count, 1) = 5;
        
        for angleN = 1:length(apertureAngles)
            summaryData.rdkApertureAngle(count, 1) = apertureAngles(angleN);
            idxT = find(eyeTrialData.rdkCoh(subN, :)==summaryData.rdkCoh(count, 1) & ...
                eyeTrialData.rdkApertureAngle(subN, :)==summaryData.rdkApertureAngle(count, 1) & ...
                eyeTrialData.rdkInternalDir(subN, :)==summaryData.rdkInternalDir(count, 1) & ...
                eyeTrialData.errorStatus(subN, :)==0);
            
            yMeanSub.responseAngle{subPlotN}(angleN, internalConN) = nanmean(eyeTrialData.response(subN, idxT));
            yStdSub.responseAngle{subPlotN}(angleN, internalConN) = nanstd(eyeTrialData.response(subN, idxT));
            
            summaryData.responseAngle(count, 1) = yMeanSub.responseAngle{subPlotN}(angleN, internalConN);
        end
        count = count+1;
    end
end
% save('summaryData.mat', 'summaryData')
%%
% load('summaryData')
if individualPlots
    for subPlotN = 1:length(plotSub)
        figure
        hold on
        
        for internalConN = 1:size(internalCons, 2)
            errorbar(apertureAngles, yMeanSub.responseAngle{subPlotN}(:, internalConN), yStdSub.responseAngle{subPlotN}(:, internalConN), 'color', colorCons(internalConN, :))
        end
        %         xticks(1:length(plotSub))
        %         xticklabels(plotSub)
        axis square
        xlim([-10 10])
        ylim([-20 20])
        legend(internalConNames, 'box', 'on', 'location', 'northwest', 'color', 'w')
        xlabel('Aperture trajectory angle')
        ylabel('Reported angle')
        saveas(gcf, [perceptFolder, '\responseAngle_', plotSub{subPlotN}, '.pdf'])
    end
    %     close all
end