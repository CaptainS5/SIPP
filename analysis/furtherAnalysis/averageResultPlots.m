% grouped bars plots for any applicable results in different conditions
initializeParas;

% choose which plot to look at
individualPlots = 0;
averagePlots = 1;

% input the parameters to plot
plotVariables = {'response', 'gainXexternal', 'gainYexternal', 'gainYaverage', 'gain2Dexternal', 'gain2Daverage', ...
    'dirGainExternal', 'dirClp', 'dirError', 'disCenterMean'};
saccadeVarStart = 100; % if there is no saccade variables, just use a super large number (larger than the number of all variables to be plotted)
% % for saccades
% checkVariables = {'number', 'meanAmp2D', 'sumAmp2D'};

% plot settings
textFontSize = 8;

load('summaryData')
load('summaryDataSub')

%%
if individualPlots
    for varN = 2:length(plotVariables)
        for subN = 1:size(names, 2)
            % plot difference from baseline
            figure
            hold on
            % line plot
            for internalConN = 1:size(yMeanDiffSub.(plotVariables{varN}){subN}, 2)
                errorbar(apertureAngles, yMeanDiffSub.(plotVariables{varN}){subN}(:, internalConN), yStdDiffSub.(plotVariables{varN}){subN}(:, internalConN), 'color', colorCons(internalConN+1, :))
            end
            
%             % bar plot
%             b = bar(yMeanDiffSub.(plotVariables{varN}){subPlotN});
%             for ii = 1:size(yMeanDiffSub.(plotVariables{varN}){subPlotN}, 2)
%                 xtips{ii} = b(ii).XEndPoints;
%                 ytips{ii} = b(ii).YEndPoints;
%                 for jj = 1:length(b(ii).YData)
%                     labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
%                 end
%                 text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
%                     'VerticalAlignment','bottom')
%                 errorbar(xtips{ii},ytips{ii},yStdDiffSub.(plotVariables{varN}){subPlotN}(:, ii), 'lineStyle', 'none', 'color', 'k')
%             end
%             xticks(1:length(apertureAngles))
%             xticklabels(apertureAngleNames)
            xlabel('Aperture angle (deg)')
            legend(internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
            
            ylabel(['Diff ', plotVariables{varN}])
            title(names{subN})
            if varN>=saccadeVarStart
                saveas(gcf, [saccadeFolder, '\individuals\sacDiff_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            elseif varN==1
                saveas(gcf, [perceptFolder, 'individuals\diff_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            else
                saveas(gcf, [pursuitFolder, '\individuals\pursuitDiff_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            end
            
            % plot all conditions
            figure
            hold on
            % line plot
            for internalConN = 1:size(internalCons, 2)
                errorbar(apertureAngles, yMeanSub.(plotVariables{varN}){subN}(:, internalConN), yStdSub.(plotVariables{varN}){subN}(:, internalConN), 'color', colorCons(internalConN, :))
            end
            
%             % bar plot
%             b = bar(yMeanSub.(plotVariables{varN}){subN});
%             for ii = 1:size(yMeanSub.(plotVariables{varN}){subN}, 2)
%                 xtips{ii} = b(ii).XEndPoints;
%                 ytips{ii} = b(ii).YEndPoints;
%                 for jj = 1:length(b(ii).YData)
%                     labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
%                 end
%                 text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
%                     'VerticalAlignment','bottom')
%                 errorbar(xtips{ii},ytips{ii},yStdSub.(plotVariables{varN}){subN}(:, ii), 'lineStyle', 'none', 'color', 'k')
%             end
%             xticks(1:length(apertureAngles))
%             xticklabels(apertureAngleNames)

            xlabel('Aperture angle (deg)')
            legend(internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
            
            ylabel(plotVariables{varN})
            title(names{subN})
            if varN>=saccadeVarStart
                saveas(gcf, [saccadeFolder, 'individuals\sac_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            elseif varN==1
                saveas(gcf, [perceptFolder, 'individuals\', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            else
                saveas(gcf, [pursuitFolder, 'individuals\pursuit_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            end
            
            close all
        end
    end
end

%%
if averagePlots
    for varN = 1:length(plotVariables)
        % plot difference from baseline
        figure
        hold on
        for internalConN = 2:size(internalCons, 2)
            meanDiffAll = NaN(size(apertureAngles));
            stdDiffAll =NaN(size(apertureAngles));
            
            for angleN = 1:length(apertureAngles)
                baselineIdx = find(summaryData.rdkCoh(:, 1)==0 & ...
                    summaryData.rdkInternalDir(:, 1)==0 & ...
                    summaryData.rdkApertureAngle(:, 1)==apertureAngles(angleN));
                baselineMean = nanmean(summaryData.(plotVariables{varN})(baselineIdx, 1));
                
                idxT = find(summaryData.rdkCoh(:, 1)==1 & ...
                    summaryData.rdkInternalDir(:, 1)==internalCons(internalConN) & ...
                    summaryData.rdkApertureAngle(:, 1)==apertureAngles(angleN));
                meanDiffAll(angleN) = nanmean(summaryData.(plotVariables{varN})(idxT, 1)-baselineMean);
                stdDiffAll(angleN) = nanstd(summaryData.(plotVariables{varN})(idxT, 1)-baselineMean);
            end
            errorbar(apertureAngles, meanDiffAll, stdDiffAll, 'color', colorCons(internalConN, :))
        end
        title('all')
        xlim([-10 10])
        legend(internalConNames(2:3), 'box', 'on', 'location', 'northwest', 'color', 'w')
        xlabel('Aperture trajectory angle (deg)')
        ylabel(['Diff ', plotVariables{varN}])
        if varN>=saccadeVarStart
            saveas(gcf, [saccadeFolder, 'sacDiff_', plotVariables{varN}, '_lineplot_all.pdf'])
        elseif varN==1
            saveas(gcf, [perceptFolder, 'diff_', plotVariables{varN}, '_lineplot_all.pdf'])
        else
            saveas(gcf, [pursuitFolder, 'pursuitDiff_', plotVariables{varN}, '_lineplot_all.pdf'])
        end
        
        % plot all conditions
        figure
        hold on
        for internalConN = 1:size(internalCons, 2)
            meanDiffAll = NaN(size(apertureAngles));
            stdDiffAll =NaN(size(apertureAngles));
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
                meanDiffAll(angleN) = nanmean(summaryData.(plotVariables{varN})(idxT, 1));
                stdDiffAll(angleN) = nanstd(summaryData.(plotVariables{varN})(idxT, 1));
            end
            errorbar(apertureAngles, meanDiffAll, stdDiffAll, 'color', colorCons(internalConN, :))
        end
        title('all')
%         axis square
        xlim([-10 10])
%         ylim([-30 30])
        legend(internalConNames, 'box', 'on', 'location', 'northwest', 'color', 'w')
        xlabel('Aperture trajectory angle (deg)')
        ylabel(plotVariables{varN})
        if varN>=saccadeVarStart
            saveas(gcf, [saccadeFolder, 'sac_', plotVariables{varN}, '_lineplot_all.pdf'])
        elseif varN==1
            saveas(gcf, [perceptFolder, plotVariables{varN}, '_lineplot_all.pdf'])
        else
            saveas(gcf, [pursuitFolder, 'pursuit_', plotVariables{varN}, '_lineplot_all.pdf'])
        end
    end
end