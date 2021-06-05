% plot perceptual results for the constant internal motion exp
initializeParas;

% choose which plots to look at
individualPlots = 1;
averagePlots = 1;

% plot settings
textFontSize = 8;

load('summaryData')
load('summaryDataSub')

%%
if individualPlots
    for subN = 1:length(names)
        figure
        hold on
        
        for internalConN = 1:size(internalCons, 2)
            errorbar(apertureAngles, yMeanSub.response{subN}(:, internalConN), yStdSub.response{subN}(:, internalConN), 'color', colorCons(internalConN, :))
        end
        %         xticks(1:length(names))
        %         xticklabels(names)
        title(names{subN})
        axis square
        xlim([-10 10])
        ylim([-30 30])
        legend(internalConNames, 'box', 'on', 'location', 'northwest', 'color', 'w')
        xlabel('Aperture trajectory angle')
        ylabel('Reported angle')
        saveas(gcf, [perceptFolder, '\individuals\responseAngle_', names{subN}, '.pdf'])
    end
    %     close all
end

%%
if averagePlots
    figure
    hold on
    
    for internalConN = 1:size(internalCons, 2)
        meanAngles = NaN(size(apertureAngles));
        stdAngles =NaN(size(apertureAngles));
        if internalCons(internalConN)==0
            rdkCoh = 0;
            rdkInternalDir = 0;
        else
            rdkCoh = 1;
            rdkInternalDir = internalCons(internalConN);
        end
        
        for angleN = 1:length(apertureAngles)
            idxT = find(summaryData.rdkCoh(:, 1)==rdkCoh & ...
                summaryData.rdkInternalDir(:, 1)==rdkInternalDir & ...
                summaryData.rdkApertureAngle(:, 1)==apertureAngles(angleN));
            meanAngles(angleN) = nanmean(summaryData.response(idxT, 1));
            stdAngles(angleN) = nanstd(summaryData.response(idxT, 1));
        end
        errorbar(apertureAngles, meanAngles, stdAngles, 'color', colorCons(internalConN, :))
    end
    title('all')
    axis square
    xlim([-10 10])
    ylim([-30 30])
    legend(internalConNames, 'box', 'on', 'location', 'northwest', 'color', 'w')
    xlabel('Aperture trajectory angle')
    ylabel('Reported angle')
    saveas(gcf, [perceptFolder, '\responseAngle_all.pdf'])
end