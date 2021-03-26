% grouped bars plots for any applicable results in different conditions

initializeParas;

% choose which plot to look at
individualPlots = 1;
averagePlots = 0;
% currently we want to compare between the participants(conditions)
% so maybe one giant plot with each condition as one group of bars

% input the parameters to plot
checkVariables = {'initialAccelerationFit2D', 'latency', 'gainX', 'gainY', 'gain2D', 'velCovX', 'velCovY', 'velCov2D', 'dirClp', 'dirError'}; % for generating summaryData and save the mat file
% checkVariables = {'initialAccelerationFit2D'};

% plot settings
textFontSize = 8;
barNames = {'down-0' 'down-0.5' 'down-1' 'up-0' 'up-0.5' 'up-1'};

% flip left directions, only analyze trials with no saccadic initiation
idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.pursuit.onsetType==0 & eyeTrialData.rdkApertureDir==180); % leftward valid trials
% eyeTrialData.pursuit.initialAccelerationFit2D(idxT) = -eyeTrialData.pursuit.initialAccelerationFit2D(idxT);
eyeTrialData.pursuit.dirClp(idxT) = 180-eyeTrialData.pursuit.dirClp(idxT);
% then also flip the internal up/down directions
eyeTrialData.pursuit.dirClp = abs(eyeTrialData.pursuit.dirClp);
% absolute error
eyeTrialData.pursuit.dirError = abs(eyeTrialData.pursuit.dirError);

% initialization
summaryData = table;
count = 1;
for subN = 1:size(names, 2)
%     for internalDirN = 1:2
        for cohN = 1:3
            % internal direction merged
            idxT = find(eyeTrialData.rdkCoh(subN, :)==cohCons(cohN) & ...
                eyeTrialData.pursuit.onsetType(subN, :)==0 & ...
                eyeTrialData.errorStatus(subN, :)==0);
            summaryData.sub(count, 1) = subN;
            summaryData.rdkCoh(count, 1) = cohCons(cohN);
            
            for varN = 1:length(checkVariables)
                if strcmp(checkVariables{varN}, 'latency') % needs to calculate from onset
                    onsetT = eyeTrialData.pursuit.onset(subN, idxT);
                    rdkOnT = eyeTrialData.frameLog.rdkOn(subN, idxT);
                    yMeanSub.(checkVariables{varN})(subN, cohN) = nanmean(onsetT-rdkOnT);
                    yStdSub.(checkVariables{varN})(subN, cohN) = nanstd(onsetT-rdkOnT);
                else
                    yMeanSub.(checkVariables{varN})(subN, cohN) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                    yStdSub.(checkVariables{varN})(subN, cohN) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                end
                summaryData.(checkVariables{varN})(count, 1) = yMeanSub.(checkVariables{varN})(subN, cohN);
            end
            
%             % get the variable for plotting, internal direction not
%             merged
%             idxT = find(eyeTrialData.rdkInternalDir(subN, :)*internalDirCons(internalDirN)>1 & ...
%                 eyeTrialData.rdkCoh(subN, :)==cohCons(cohN) & ...
%                 eyeTrialData.pursuit.onsetType(subN, :)==0 & ...
%                 eyeTrialData.errorStatus(subN, :)==0);
%             summaryData.sub(count, 1) = subN;
%             summaryData.internalDir(count, 1) = internalDirCons(internalDirN);
%             summaryData.rdkCoh(count, 1) = cohCons(cohN);
%             
%             for varN = 1:length(checkVariables)
%                 if strcmp(checkVariables{varN}, 'latency') % needs to calculate from onset
%                     onsetT = eyeTrialData.pursuit.onset(subN, idxT);
%                     rdkOnT = eyeTrialData.frameLog.rdkOn(subN, idxT);
%                     yMeanSub.(checkVariables{varN}){subN}(internalDirN, cohN) = nanmean(onsetT-rdkOnT);
%                     yStdSub.(checkVariables{varN}){subN}(internalDirN, cohN) = nanstd(onsetT-rdkOnT);
%                 else
%                     yMeanSub.(checkVariables{varN}){subN}(internalDirN, cohN) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
%                     yStdSub.(checkVariables{varN}){subN}(internalDirN, cohN) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
%                 end
%                 summaryData.(checkVariables{varN})(count, 1) = yMeanSub.(checkVariables{varN}){subN}(internalDirN, cohN);
%             end
            count = count+1;
        end
%     end
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
            legend(cohNames, 'box', 'on', 'location', 'best', 'color', 'w')
            for ii = 1:size(yMeanSub.(checkVariables{varN}), 2)
                xtips{ii} = b(ii).XEndPoints;
                ytips{ii} = b(ii).YEndPoints;
                for jj = 1:length(b(ii).YData)
                    labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
                end
                text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
                    'VerticalAlignment','bottom')
            end
            xticks(1:length(names))
            xticklabels(names)
            ylabel(checkVariables{varN})
            %                 title(names{subN})
            saveas(gcf, [pursuitFolder, '\', checkVariables{varN}, '_barplot_allPilots.pdf'])
            
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