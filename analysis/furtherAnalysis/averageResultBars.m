% grouped bars plots for any applicable results in different conditions

initializeParas;

% choose which plot to look at
individualPlots = 1;
averagePlots = 0;
% currently we want to compare between the participants(conditions)
% so maybe one giant plot with each condition as one group of bars

% input the parameters to plot
checkVariables = {'initialAccelerationFit2D', 'latency', 'gainXexternal', 'gainXaverage', 'gainYaverage', 'velCovX', 'velCovY', 'velCov2D', 'dirClp', 'dirError', 'dirGain'}; % for generating summaryData and save the mat file
% checkVariables = {'dirError'};
% % for saccades
% checkVariables = {'number', 'meanAmp2D', 'sumAmp2D'};

% plot settings
textFontSize = 8;
% plotSub = {'w00' 'w01' 'w02' 'w03' 'w04' 'w05' 'w06' 'w08'}; % which participants to plot
plotSub = {'w01' 'w08' 'w03' 'w09'};

% % w00-w06, w08
% cons = cohCons;
% barNames = {'down-0' 'down-0.5' 'down-1' 'up-0' 'up-0.5' 'up-1'};

% % w07
% cons = [5, 0.5; 5, 1; 10, 0.5; 10, 1]; % first column is internal speed, second column is coh
% barNames = {'5-coh 0.5' '5-coh 1' '10-coh 0.5' '10-coh 1'}; % w07

% % w07 & w10
% barNames = {'unlimited dot lifetime' '200 ms dot lifetime'};

% w01, w08, w03, w09
cons = [45; 90; 135]; % absolute internal direction
legendNames = {'unlimited lifetime' '200ms lifetime'};

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
    for conN = 1:size(cons, 1) % w00-w06, w08; w07; w01, w08, w03, w09
        % for conN = 1:1 % w07&w10
        
        summaryData.sub(count, 1) = subPlotN;
        
        % internal direction merged
        
        %         % w00-w06, w08:
        %         idxT = find(eyeTrialData.rdkCoh(subN, :)==cons(conN) & ...
        %             eyeTrialData.pursuit.onsetType(subN, :)==0 & ...
        %             eyeTrialData.errorStatus(subN, :)==0); % for pursuit conditions
        %         summaryData.rdkCoh(count, 1) = cohCons(cohN);
        
        %                 % w07
        %                 idxT = find(eyeTrialData.rdkCoh(subN, :)==cons(conN, 2) & ...
        %                     eyeTrialData.rdkInternalSpeed(subN, :)==cons(conN, 1) & ...
        %                     eyeTrialData.errorStatus(subN, :)==0);
        %                 summaryData.rdkCoh(count, 1) = cons(conN, 2);
        %                 summaryData.rdkInternalSpeed(count, 1) = cons(conN, 1);
        
        % % w07&w10
        % idxT = find(eyeTrialData.rdkCoh(subN, :)==1 & ...
        %     eyeTrialData.rdkInternalSpeed(subN, :)==5 & ...
        %     eyeTrialData.errorStatus(subN, :)==0);
        % summaryData.rdkCoh(count, 1) = 1;
        % summaryData.rdkInternalSpeed(count, 1) = 5;
        
        % w01, w08, w03, w09
        if strcmp(plotSub(subPlotN), 'w09')
            idxT = find(eyeTrialData.rdkCoh(subN, :)==1 & ...
                eyeTrialData.rdkInternalSpeed(subN, :)==5 & ...
                abs(eyeTrialData.rdkInternalDir(subN, :))==cons(conN) & ...
                eyeTrialData.errorStatus(subN, :)==0);
            summaryData.rdkInternalDir(count, 1) = cons(conN);
        else
            idxT = find(eyeTrialData.rdkCoh(subN, :)==1 & ...
                eyeTrialData.rdkInternalSpeed(subN, :)==5 & ...
                eyeTrialData.errorStatus(subN, :)==0);
            summaryData.rdkInternalDir(count, 1) = cons(subPlotN);
            conN = subPlotN;
        end
        summaryData.rdkCoh(count, 1) = 1;
        summaryData.rdkInternalSpeed(count, 1) = 5;
        
        for varN = 1:length(checkVariables)
            if strcmp(checkVariables{varN}, 'latency') % needs to calculate from onset
                onsetT = eyeTrialData.pursuit.onset(subN, idxT);
                rdkOnT = eyeTrialData.frameLog.rdkOn(subN, idxT);
                
                % w01, w08, w03, w09
                if strcmp(plotSub(subPlotN), 'w09')
                    yMeanSub.(checkVariables{varN})(conN, 2) = nanmean(onsetT-rdkOnT);
                    yStdSub.(checkVariables{varN})(conN, 2) = nanstd(onsetT-rdkOnT);
                else
                    yMeanSub.(checkVariables{varN})(conN, 1) = nanmean(onsetT-rdkOnT);
                    yStdSub.(checkVariables{varN})(conN, 1) = nanstd(onsetT-rdkOnT);
                end
                
                %                 % w00-w06, w08; w07; w07&w10
                %                 yMeanSub.(checkVariables{varN})(subPlotN, conN) = nanmean(onsetT-rdkOnT);
                %                 yStdSub.(checkVariables{varN})(subPlotN, conN) = nanstd(onsetT-rdkOnT);
            elseif strcmp(checkVariables{varN}, 'dirClp') % needs to calculate from vectors
                dir = atan2(eyeTrialData.pursuit.dirClpY(subN, idxT), eyeTrialData.pursuit.dirClpX(subN, idxT))/pi*180;
                
                %                 % w00-w06, w08; w07; w07&w10
                %                 yMeanSub.(checkVariables{varN})(subPlotN, conN) = nanmean(dir);
                %                 yStdSub.(checkVariables{varN})(subPlotN, conN) = nanstd(dir);
                
                % w01, w08, w03, w09
                if strcmp(plotSub(subPlotN), 'w09')
                    yMeanSub.(checkVariables{varN})(conN, 2) = nanmean(dir);
                    yStdSub.(checkVariables{varN})(conN, 2) = nanstd(dir);
                else
                    yMeanSub.(checkVariables{varN})(conN, 1) = nanmean(dir);
                    yStdSub.(checkVariables{varN})(conN, 1) = nanstd(dir);
                end
            else
                % w01, w08, w03, w09
                if strcmp(plotSub(subPlotN), 'w09')
%                     % saccade parameters
%                     yMeanSub.(checkVariables{varN})(2, conN) = nanmean(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
%                     yStdSub.(checkVariables{varN})(2, conN) = nanstd(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                    % pursuit parameters
                    yMeanSub.(checkVariables{varN})(conN, 2) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                    yStdSub.(checkVariables{varN})(conN, 2) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                else
%                     % saccade parameters
%                     yMeanSub.(checkVariables{varN})(1, conN) = nanmean(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
%                     yStdSub.(checkVariables{varN})(1, conN) = nanstd(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                    % pursuit parameters
                    yMeanSub.(checkVariables{varN})(conN, 1) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                    yStdSub.(checkVariables{varN})(conN, 1) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                end
                
                %                 % w00-w06, w08; w07; w07&w10
                %                 % saccade parameters
                %                 yMeanSub.(checkVariables{varN})(subPlotN, conN) = nanmean(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                %                 yStdSub.(checkVariables{varN})(subPlotN, conN) = nanstd(eyeTrialData.saccades.(checkVariables{varN})(subN, idxT));
                %                 %                                 % pursuit parameters
                %                 %                                 yMeanSub.(checkVariables{varN})(subPlotN, conN) = nanmean(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
                %                 %                                 yStdSub.(checkVariables{varN})(subPlotN, conN) = nanstd(eyeTrialData.pursuit.(checkVariables{varN})(subN, idxT));
            end
            %             % w00-w06, w08; w07; w07&w10
            %             summaryData.(checkVariables{varN})(count, 1) = yMeanSub.(checkVariables{varN})(subPlotN, conN);
            
            % w01, w08, w03, w09
            if strcmp(plotSub(subPlotN), 'w09')
                 summaryData.(checkVariables{varN})(count, 1) = yMeanSub.(checkVariables{varN})(conN, 2);
            else
                 summaryData.(checkVariables{varN})(count, 1) = yMeanSub.(checkVariables{varN})(conN, 1);
            end
            
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
        % w01, w08, w03, w09
        if ~strcmp(plotSub(subPlotN), 'w09')
            break
        end
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
        %         % w07
        %         b = bar(yMeanSub.(checkVariables{varN})(end, :));
        
        % w00-w06, w08; w07&w10; w01, w08, w03, w09
        b = bar(yMeanSub.(checkVariables{varN}));
        
        for ii = 1:size(yMeanSub.(checkVariables{varN}), 2)
            xtips{ii} = b(ii).XEndPoints;
            ytips{ii} = b(ii).YEndPoints;
            for jj = 1:length(b(ii).YData)
                labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
            end
            text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
            % w00-w06, w08; w07&w10; w01, w08, w03, w09
            errorbar(xtips{ii},ytips{ii},yStdSub.(checkVariables{varN})(:, ii), 'lineStyle', 'none', 'color', 'k')
            
            %             % w07
            %             errorbar(xtips{ii},ytips{ii},yStdSub.(checkVariables{varN})(end, :), 'lineStyle', 'none', 'color', 'k')
        end
        %         % w00-w06, w08
        %         xticks(1:length(plotSub))
        %         xticklabels(plotSub)
        %         legend(cohNames, 'box', 'on', 'location', 'best', 'color', 'w') % w00-w06
        
        % w01, w08, w03, w09
        xticks(1:3)
        xticklabels({'45', '90', '135'})
        legend(legendNames, 'box', 'on', 'location', 'best', 'color', 'w') % w00-w06
        
%         % w07; w07&w10
%         xticks(1:length(barNames))
%         xticklabels(barNames)
        
        ylabel(checkVariables{varN})
        %                 title(names{subN})
%         saveas(gcf, [saccadeFolder, '\pursuit_sac_', checkVariables{varN}, '_barplot_w09w010308.pdf'])
                saveas(gcf, [pursuitFolder, '\fixation_', checkVariables{varN}, '_barplot_w09w010308.pdf'])
        
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