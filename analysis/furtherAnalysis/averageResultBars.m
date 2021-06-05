% grouped bars plots for any applicable results in different conditions
initializeParas;

% choose which plot to look at
individualPlots = 1;
averagePlots = 0;
% calculate the difference from baseline, then compare two internal
% conditions in different aperture angles

% input the parameters to plot
% checkVariables = {'initialAccelerationFit2D', 'latency', 'gainXexternal', 'gainXaverage', 'gainYinternal', 'gainYaverage', 'velCovX', 'velCovY', 'velCov2D', 'dirClp', 'dirError', 'dirGain'}; % for generating summaryData and save the mat file
plotVariables = {'gainXexternal', 'gainYexternal', 'gainYaverage', 'gain2Dexternal', 'gain2Daverage', ...
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
    for varN = 1:length(plotVariables)
        for subN = 1:size(names, 2)
            figure
            hold on
            
            %             % plot difference from baseline
            %             b = bar(yMeanDiffSub.(checkVariables{varN}){subPlotN});
            %             for ii = 1:size(yMeanDiffSub.(checkVariables{varN}){subPlotN}, 2)
            %                 xtips{ii} = b(ii).XEndPoints;
            %                 ytips{ii} = b(ii).YEndPoints;
            %                 for jj = 1:length(b(ii).YData)
            %                     labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
            %                 end
            %                 text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
            %                     'VerticalAlignment','bottom')
            %                 errorbar(xtips{ii},ytips{ii},yStdDiffSub.(checkVariables{varN}){subPlotN}(:, ii), 'lineStyle', 'none', 'color', 'k')
            %             end
            %             xticks(1:length(apertureAngles))
            %             xticklabels(apertureAngleNames)
            %             xlabel('Aperture angle (deg)')
            %             legend(internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
            %
            %             ylabel(['Diff ', checkVariables{varN}])
            %             title(names{subPlotN})
            %             if varN>=saccadeVarStart
            %                 saveas(gcf, [saccadeFolder, '\individuals\sacDiff_', checkVariables{varN}, '_barplot_', names{subPlotN}, '.pdf'])
            %             else
            %                 saveas(gcf, [pursuitFolder, '\individuals\pursuitDiff_', checkVariables{varN}, '_barplot_', names{subPlotN}, '.pdf'])
            %             end
            
            % plot all conditions
            b = bar(yMeanSub.(plotVariables{varN}){subN});
            for ii = 1:size(yMeanSub.(plotVariables{varN}){subN}, 2)
                xtips{ii} = b(ii).XEndPoints;
                ytips{ii} = b(ii).YEndPoints;
                for jj = 1:length(b(ii).YData)
                    labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
                end
                text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
                    'VerticalAlignment','bottom')
                errorbar(xtips{ii},ytips{ii},yStdSub.(plotVariables{varN}){subN}(:, ii), 'lineStyle', 'none', 'color', 'k')
            end
            xticks(1:length(apertureAngles))
            xticklabels(apertureAngleNames)
            xlabel('Aperture angle (deg)')
            legend(internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
            
            ylabel(plotVariables{varN})
            title(names{subN})
            if varN>=saccadeVarStart
                saveas(gcf, [saccadeFolder, '\individuals\sac_', plotVariables{varN}, '_barplot_', names{subN}, '.pdf'])
            else
                saveas(gcf, [pursuitFolder, '\individuals\pursuit_', plotVariables{varN}, '_barplot_', names{subN}, '.pdf'])
            end
            
            %             figure
            %             hold on
            %             if strcmp(checkVariables{varN}, 'accuracy')
            %                 % fixation vs. pursuit in one plot
            %                 b = bar(eyeConditions, yMeanSub.(checkVariables{varN}){subN});
            %                 legend(barNames, 'box', 'on', 'location', 'southeast', 'color', 'w')
            %                 for ii = 1:size(yMeanSub.(checkVariables{varN}){subN}, 2)
            %                     xtips{ii} = b(ii).XEndPoints;
            %                     ytips{ii} = b(ii).YEndPoints;
            %                     for jj = 1:length(b(ii).YData)
            %                         labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
            %                     end
            %                     text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
            %                         'VerticalAlignment','bottom')
            %                 end
            %                 xticks([0, 1])
            %                 xticklabels({'fixation', 'pursuit'})
            %                 ylabel(checkVariables{varN})
            %                 title(names{subN})
            %                 saveas(gcf, [perceptFolder, 'individuals\', checkVariables{varN}, '_barplot_', names{subN}, '.pdf'])
            %             else
            %                 % 2 x 2, instruction on x axis, dirSD in colour
            %                 b = bar(instructionCons, yMeanSub.(checkVariables{varN}){subN});
            %                 for ii = 1:size(yMeanSub.(checkVariables{varN}){subN}, 2)
            %                     xtips{ii} = b(ii).XEndPoints;
            %                     ytips{ii} = b(ii).YEndPoints;
            %                     for jj = 1:length(b(ii).YData)
            %                         labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
            %                     end
            %                     text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
            %                         'VerticalAlignment','bottom')
            %                     % adding errorbars
            %                     errorbar(xtips{ii}, ytips{ii}, yStdSub.(checkVariables{varN}){subN}(:, ii)', 'LineStyle', 'none')
            %                 end
            %                 legend(sdNames, 'box', 'on', 'location', 'southeast', 'color', 'w')
            %                 xticks([0, 1])
            %                 xticklabels({'fast', 'accurate'})
            %                 ylabel(checkVariables{varN})
            %                 title(names{subN})
            %                 saveas(gcf, [pursuitFolder, 'individuals\', checkVariables{varN}, '_barplot_', names{subN}, '.pdf'])
            %             end
        end
    end
    %     close all
end

% % calculate average values
% for internalDirN = 1:2
%     for cohN = 1:2
%         for eyeN = 1:2
%             idxT = find(summaryData.instruction(:, 1)==instructionCons(internalDirN) ...
%                 & summaryData.rdkDirSD(:, 1)==sdCons(cohN) & summaryData.eyeCondition(:, 1)==eyeConditions(eyeN));
%             for varN = 1:length(checkVariables)
%                 if strcmp(checkVariables{varN}, 'accuracy')
%                     meanY.(checkVariables{varN})(eyeN, (internalDirN-1)*2+cohN) = mean(summaryData.accuracy(idxT, 1));
%                     stdY.(checkVariables{varN})(eyeN, (internalDirN-1)*2+cohN) = std(summaryData.accuracy(idxT, 1));
%                 elseif eyeConditions(eyeN)==1
%                     meanY.(checkVariables{varN})(internalDirN, cohN) = nanmean(summaryData.(checkVariables{varN})(idxT, 1));
%                     stdY.(checkVariables{varN})(internalDirN, cohN) = nanstd(summaryData.(checkVariables{varN})(idxT, 1));
%                 end
%             end
%         end
%     end
% end
%
% if averagePlots
%     for varN = 1:length(checkVariables)
%         figure
%         hold on
%
%         if strcmp(checkVariables{varN}, 'accuracy')
%             % fixation vs. pursuit in one plot
%             b = bar(eyeConditions, meanY.(checkVariables{varN}));
%         else
%             % 2 x 2, instruction on x axis, dirSD in colour
%             b = bar(instructionCons, meanY.(checkVariables{varN}));
%         end
%
%         for ii = 1:size(meanY.(checkVariables{varN}), 2)
%             xtips{ii} = b(ii).XEndPoints;
%             ytips{ii} = b(ii).YEndPoints;
%             for jj = 1:length(b(ii).YData)
%                 labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
%             end
%             text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
%                 'VerticalAlignment','bottom')
%             % adding errorbars
%             errorbar(xtips{ii}, ytips{ii}, stdY.(checkVariables{varN})(:, ii)', 'LineStyle', 'none')
%         end
%         ylabel(checkVariables{varN})
%         title('all')
%         if strcmp(checkVariables{varN}, 'accuracy')
%             legend(barNames, 'box', 'on', 'location', 'southeast', 'color', 'w')
%             xticks([0, 1])
%             xticklabels({'fixation', 'pursuit'})
%             saveas(gcf, [perceptFolder, checkVariables{varN}, '_barplot_all.pdf'])
%         else
%             legend(sdNames, 'box', 'on', 'location', 'southeast', 'color', 'w')
%             xticks([0, 1])
%             xticklabels({'fast', 'accurate'})
%             saveas(gcf, [pursuitFolder, checkVariables{varN}, '_barplot_all.pdf'])
%         end
%     end
% end