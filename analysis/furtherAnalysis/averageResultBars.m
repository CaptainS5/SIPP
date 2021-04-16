% grouped bars plots for any applicable results in different conditions

initializeParas;

% choose which plot to look at
individualPlots = 1;
averagePlots = 0;
% currently we want to compare between the participants(conditions)
% so maybe one giant plot with each condition as one group of bars

% input the parameters to plot
checkVariables = {'initialAccelerationFit2D', 'latency', 'gainXexternal', 'gainXaverage', 'gainYinternal', 'gainYaverage', 'velCovX', 'velCovY', 'velCov2D', 'dirClp', 'dirError', 'dirGain'}; % for generating summaryData and save the mat file
% checkVariables = {'gainXexternal', 'gainYinternal', 'velCovX', 'velCovY', 'velCov2D', 'dirClp', 'dirError', 'dirGain'};
% % for saccades
% checkVariables = {'number', 'meanAmp2D', 'sumAmp2D'};

% plot settings
textFontSize = 8;
plotSub = {'xw0' 'dc0' 'ib0'};
cons = [45; 90; 135]; % absolute internal direction
legendNames = {'45' '90' '135'};

% flip left directions
idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkApertureDir==180); % leftward valid trials
eyeTrialData.pursuit.dirClpX(idxT) = -eyeTrialData.pursuit.dirClpX(idxT);
% then also flip the internal up/down directions
idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkInternalDir<0); % leftward valid trials
eyeTrialData.pursuit.dirClpY(idxT) = -eyeTrialData.pursuit.dirClpY(idxT);

% initialization
summaryData = table;
count = 1;
%%
for subPlotN = 1:size(plotSub, 2)
    subN = find(strcmp(names, plotSub{subPlotN}));
    for conN = 1:size(cons, 1) 
        summaryData.sub(count, 1) = subPlotN;
        summaryData.rdkCoh(count, 1) = 1;
        summaryData.rdkInternalSpeed(count, 1) = 5;
        % internal direction merged
        idxT = find(eyeTrialData.rdkCoh(subN, :)==1 & ...
            eyeTrialData.rdkInternalSpeed(subN, :)==5 & ...
            abs(eyeTrialData.rdkInternalDir(subN, :))==cons(conN) & ...
            eyeTrialData.pursuit.onsetType(subN, :)==0 & ...
            eyeTrialData.errorStatus(subN, :)==0);
        summaryData.rdkInternalDir(count, 1) = cons(conN);
                
        for varN = 1:length(checkVariables)
            if strcmp(checkVariables{varN}, 'latency') % needs to calculate from onset
                onsetT = eyeTrialData.pursuit.onset(subN, idxT);
                rdkOnT = eyeTrialData.frameLog.rdkOn(subN, idxT);
                
                yMeanSub.(checkVariables{varN})(subPlotN, conN) = nanmean(onsetT-rdkOnT);
                yStdSub.(checkVariables{varN})(subPlotN, conN) = nanstd(onsetT-rdkOnT);
            elseif strcmp(checkVariables{varN}, 'dirClp') % needs to calculate from vectors
                dir = atan2(eyeTrialData.pursuit.dirClpY(subN, idxT), eyeTrialData.pursuit.dirClpX(subN, idxT))/pi*180;
                
                yMeanSub.(checkVariables{varN})(subPlotN, conN) = nanmean(dir);
                yStdSub.(checkVariables{varN})(subPlotN, conN) = nanstd(dir);
            else
                %                 % saccade parameters
                %                 yMeanSub.(checkVariables{varN})(subPlotN, conN) = nanmean(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                %                 yStdSub.(checkVariables{varN})(subPlotN, conN) = nanstd(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                % pursuit parameters
                yMeanSub.(checkVariables{varN})(subPlotN, conN) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                yStdSub.(checkVariables{varN})(subPlotN, conN) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
            end
            summaryData.(checkVariables{varN})(count, 1) = yMeanSub.(checkVariables{varN})(subPlotN, conN);
            
        end
        count = count+1;
    end
end
% save('summaryData.mat', 'summaryData')
%%
% load('summaryData')
if individualPlots
    for varN = 1:length(checkVariables)
        % for the pilot data, do a giant plot including all conditions...
        figure
        hold on
        b = bar(yMeanSub.(checkVariables{varN}));
        
        for ii = 1:size(yMeanSub.(checkVariables{varN}), 2)
            xtips{ii} = b(ii).XEndPoints;
            ytips{ii} = b(ii).YEndPoints;
            for jj = 1:length(b(ii).YData)
                labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
            end
            text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
            errorbar(xtips{ii},ytips{ii},yStdSub.(checkVariables{varN})(:, ii), 'lineStyle', 'none', 'color', 'k')
        end
        xticks(1:length(plotSub))
        xticklabels(plotSub)
        legend(legendNames, 'box', 'on', 'location', 'best', 'color', 'w') 
        
        ylabel(checkVariables{varN})
        %                 title(names{subN})
%         saveas(gcf, [saccadeFolder, '\pursuit_sac_', checkVariables{varN}, '_barplot_xw0dc0.pdf'])
                saveas(gcf, [pursuitFolder, '\fixation_', checkVariables{varN}, '_barplot_xw0dc0.pdf'])
        
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