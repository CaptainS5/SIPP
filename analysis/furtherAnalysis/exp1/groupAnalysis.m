% group people by some perception/eye movements features, then plot
% accordingly
initializeParas;

%% choose which plot to look at
individualPlots = 0;
averagePlots = 1;
plotVarStart = 6;
plotVarEnd = 13;

%% first, generate the list of different groups of participants to separately look at
% calculate the average difference in all aperture angles in each
% internal condition, then separate the participants into assimilation and
% contrast groups
for subN = 1:size(names, 2)
    idx1 = find(summaryDataDiff.sub==subN & ...
        summaryDataDiff.rdkInternalDir==-90);
    idx2 = find(summaryDataDiff.sub==subN & ...
        summaryDataDiff.rdkInternalDir==90);
    tempData(subN, 1) = mean(summaryDataDiff.response(idx2))-mean(summaryDataDiff.response(idx1));
end
subGroup{1} = find(tempData>=1); % assimilation group
subGroup{2} = find(tempData<=-1); % contrast group
subGroup{3} = find(tempData<1 & tempData>-1); % those who doesn't have a large bias

%% saving csv files for plotting in R
% for groupN = 1:2
%     summaryDataSubgroup = table;
%
%     subList = subGroup{groupN};
%     V = intersect(summaryDataDiff.sub, subList);
%     bSummary = ismember(summaryDataDiff.sub, V);
%     idxT = find(bSummary);
%
%     summaryDataDiffSubgroup = summaryDataDiff(idxT, :);
%     writetable(summaryDataDiffSubgroup, [RFolder, 'summaryDataDiffGroup2', num2str(groupN), '.csv'])
% end

%%
if individualPlots
    for varN = plotVarStart:plotVarEnd
        for groupN = 1:3
            subList = subGroup{groupN};
            for subN = 1:length(subList)
                % plot difference from baseline
                figure
                hold on
                % line plot
                for internalConN = 1:size(yMeanDiffSub.(plotVariables{varN}){subList(subN)}, 2)
                    errorbar(apertureAngles, yMeanDiffSub.(plotVariables{varN}){subList(subN)}(:, internalConN), yStdDiffSub.(plotVariables{varN}){subList(subN)}(:, internalConN), 'color', colorCons(internalConN+1, :))
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
                title(names{subList(subN)})
                %                 if varN>=saccadeVarStart
                %                     saveas(gcf, [saccadeFolder, '\individuals\sacDiff_', plotVariables{varN}, '_lineplot_', names{subList(subN)}, '.pdf'])
                %                 elseif varN==1
                %                     saveas(gcf, [perceptFolder, 'individuals\diff_', plotVariables{varN}, '_lineplot_', names{subList(subN)}, '.pdf'])
                %                 else
                %                     saveas(gcf, [pursuitFolder, '\individuals\pursuitDiff_', plotVariables{varN}, '_lineplot_', names{subList(subN)}, '.pdf'])
                %                 end
                %
                %             % plot all conditions
                %             figure
                %             hold on
                %             % line plot
                %             for internalConN = 1:size(internalCons, 2)
                %                 errorbar(apertureAngles, yMeanSub.(plotVariables{varN}){subN}(:, internalConN), yStdSub.(plotVariables{varN}){subN}(:, internalConN), 'color', colorCons(internalConN, :))
                %             end
                %
                % %             % bar plot
                % %             b = bar(yMeanSub.(plotVariables{varN}){subN});
                % %             for ii = 1:size(yMeanSub.(plotVariables{varN}){subN}, 2)
                % %                 xtips{ii} = b(ii).XEndPoints;
                % %                 ytips{ii} = b(ii).YEndPoints;
                % %                 for jj = 1:length(b(ii).YData)
                % %                     labels{ii}{jj} = num2str(b(ii).YData(jj), '%.2f');
                % %                 end
                % %                 text(xtips{ii},ytips{ii},labels{ii},'HorizontalAlignment','center',...
                % %                     'VerticalAlignment','bottom')
                % %                 errorbar(xtips{ii},ytips{ii},yStdSub.(plotVariables{varN}){subN}(:, ii), 'lineStyle', 'none', 'color', 'k')
                % %             end
                % %             xticks(1:length(apertureAngles))
                % %             xticklabels(apertureAngleNames)
                %
                %             xlabel('Aperture angle (deg)')
                %             legend(internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
                %
                %             ylabel(plotVariables{varN})
                %             title(names{subN})
                %             if varN>=saccadeVarStart
                %                 saveas(gcf, [saccadeFolder, 'individuals\sac_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
                %             elseif varN==1
                %                 saveas(gcf, [perceptFolder, 'individuals\', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
                %             else
                %                 saveas(gcf, [pursuitFolder, 'individuals\pursuit_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
                %             end
                
                %             close all
            end
        end
    end
end

%%
if averagePlots
    for varN = plotVarStart:plotVarEnd
        if varN==2 % one value per participant, compare across groups
            for groupN = 1:size(subGroup, 2)
                subList = subGroup{groupN};
                [C, idx, ib] = intersect(summaryDataDiff.sub, subList);
                tLatencyM(groupN) = mean(summaryDataDiff.turningPoint(idx));
            end
            figure
            plot(tLatencyM)
            xlabel('assimilation-contrast-none')
            saveas(gcf, [pursuitFolder, 'pursuitTurningPoint_lineplot_subGroups.', num2str(size(subGroup, 2)), 'pdf'])
        end
        if varN~=2
            for groupN = 1:size(subGroup, 2)
                subList = subGroup{groupN};
                %
                % plot difference from baseline
                figure
                hold on
                for internalConN = 2:size(internalCons, 2)
                    meanDiffAll = NaN(size(apertureAngles));
                    stdDiffAll =NaN(size(apertureAngles));
                    
                    for angleN = 1:length(apertureAngles)
                        idxT = find(summaryDataDiff.rdkCoh(:, 1)==1 & ...
                            summaryDataDiff.rdkInternalDir(:, 1)==internalCons(internalConN) & ...
                            summaryDataDiff.rdkApertureAngle(:, 1)==apertureAngles(angleN));
                        % locate the participants in this subgroup
                        [C, ia, ib] = intersect(summaryDataDiff.sub(idxT, 1), subList);
                        
                        meanDiffAll(angleN) = nanmean(summaryDataDiff.(plotVariables{varN})(idxT(ia), 1));
                        stdDiffAll(angleN) = nanstd(summaryDataDiff.(plotVariables{varN})(idxT(ia), 1));
                    end
                    errorbar(apertureAngles, meanDiffAll, stdDiffAll, 'color', colorCons(internalConN, :))
                end
                title(['subGroup2', num2str(groupN), ', n=', num2str(length(subList))])
                %                 xlim([-10 10])
                %                 ylim([-15 15])
                legend(internalConNames(2:3), 'box', 'on', 'location', 'northwest', 'color', 'w')
                xlabel('Aperture trajectory angle (deg)')
                ylabel(['Diff ', plotVariables{varN}])
                if varN>=saccadeVarStart
                    saveas(gcf, [saccadeFolder, 'sacDiff_', plotVariables{varN}, '_lineplot_subGroup', num2str(size(subGroup, 2)), num2str(groupN), '.pdf'])
                elseif varN==1
                    saveas(gcf, [perceptFolder, 'diff_', plotVariables{varN}, '_lineplot_subGroup', num2str(size(subGroup, 2)), num2str(groupN), '.pdf'])
                else
                    saveas(gcf, [pursuitFolder, 'pursuitDiff_', plotVariables{varN}, '_lineplot_subGroup', num2str(size(subGroup, 2)), num2str(groupN), '.pdf'])
                end
                
%                 % plot all conditions
%                 figure
%                 hold on
%                 for internalConN = 1:size(internalCons, 2)
%                     meanAll = NaN(size(apertureAngles));
%                     stdAll =NaN(size(apertureAngles));
%                     if internalCons(internalConN)==0
%                         rdkCoh = 0;
%                         rdkInternalDir = 0;
%                     else
%                         rdkCoh = 1;
%                         rdkInternalDir = internalCons(internalConN);
%                     end
%                     
%                     for angleN = 1:length(apertureAngles)
%                         idxT = find(summaryData.rdkCoh(:, 1)==rdkCoh & ...
%                             summaryData.rdkInternalDir(:, 1)==rdkInternalDir & ...
%                             summaryData.rdkApertureAngle(:, 1)==apertureAngles(angleN));
%                         % locate the participants in this subgroup
%                         [C, ia, ib] = intersect(summaryData.sub(idxT, 1), subList);
%                         
%                         meanAll(angleN) = nanmean(summaryData.(plotVariables{varN})(idxT(ia), 1));
%                         stdAll(angleN) = nanstd(summaryData.(plotVariables{varN})(idxT(ia), 1));
%                     end
%                     errorbar(apertureAngles, meanAll, stdAll, 'color', colorCons(internalConN, :))
%                 end
%                 title(['subGroup2', num2str(groupN), ', n=', num2str(length(subList))])
% %                         axis square
%                 xlim([-10 10])
%                 %         ylim([-30 30])
%                 legend(internalConNames, 'box', 'on', 'location', 'northwest', 'color', 'w')
%                 xlabel('Aperture trajectory angle (deg)')
%                 ylabel(plotVariables{varN})
%                 if varN>=saccadeVarStart
%                     saveas(gcf, [saccadeFolder, 'sac_', plotVariables{varN}, '_lineplot_subGroup', num2str(size(subGroup, 2)), num2str(groupN), '.pdf'])
%                 elseif varN==1
%                     saveas(gcf, [perceptFolder, plotVariables{varN}, '_lineplot_subGroup', num2str(size(subGroup, 2)), num2str(groupN), '.pdf'])
%                 else
%                     saveas(gcf, [pursuitFolder, 'pursuit_', plotVariables{varN}, '_lineplot_subGroup', num2str(size(subGroup, 2)), num2str(groupN), '.pdf'])
%                 end
            end
        end
    end
end