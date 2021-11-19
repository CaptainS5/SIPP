% grouped bars plots for any applicable results in different conditions
initializeParas;

% choose which plot to look at
individualPlots = 1;
averagePlots = 0;
plotVarStart = 1;
plotVarEnd = 1;

% Psychometric function fitting settings
PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull,
%PAL_Quick, PAL_logQuick, PAL_Logistic
%PAL_CumulativeNormal, PAL_HyperbolicSecant

%Threshold, Slope, and lapse rate are free parameters, guess is fixed
paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter

%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = -1:.01:3;
searchGrid.beta = logspace(0,0.25,100);
searchGrid.gamma = 0:0.01:0.05;  %scalar here (since fixed) but may be vector
searchGrid.lambda = 0:0.001:0.05;  %ditto

%%
if individualPlots
    for varN = plotVarStart:plotVarEnd
        for subN = 1:size(names, 2)
            % plot all conditions
%             figure
            for perturbN = 1:2
%                 subplot(2, 1, perturbN)
figure
                hold on
                if varN==1 % psychometric functions for response
                    for internalConN = 1:size(internalCons, 2)
                        numRight = yMeanSub.choiceUpNum{subN, perturbN}(:, internalConN);
                        outOfNum =yMeanSub.totalTrialN{subN, perturbN}(:, internalConN);
                        % Perform fit
                        [paramsValues LL exitflag output] = PAL_PFML_Fit(apertureAngles, numRight, ...
                            outOfNum, searchGrid, paramsFree, PF, 'lapseLimits',[0 0.1]);
                        
                        % plotting
                        ProportionCorrectObserved=numRight./outOfNum;
                        StimLevelsFineGrain=[min(apertureAngles):max(apertureAngles)./1000:max(apertureAngles)];
                        ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);
                        
                        if internalConN==1
                            f{internalConN} = plot(StimLevelsFineGrain, ProportionCorrectModel,'-','color', colorCons(internalConN, :), 'linewidth', 2);
                        else
                            f{internalConN} = plot(StimLevelsFineGrain, ProportionCorrectModel,'--','color', colorCons(internalConN, :), 'linewidth', 2);
                        end
                        plot(apertureAngles, ProportionCorrectObserved,'.', 'color', colorCons(internalConN, :), 'markersize', 30);
                    end
                    set(gca, 'fontsize',16);
                    set(gca, 'Xtick', apertureAngles);
                    axis([min(apertureAngles) max(apertureAngles) 0 1]);
                    xlabel('Aperture perturbation angle');
                    ylabel('Proportion up');
                    legend([f{:}], internalConNames, 'box', 'off', 'location', 'northwest')
                    title([names{subN}, ', perturbPhase', num2str(perturbPhase(perturbN))])
                else
                    % line plot
                    for internalConN = 1:size(internalCons, 2)
                        errorbar(apertureAngles, yMeanSub.(plotVariables{varN}){subN, perturbN}(:, internalConN), yStdSub.(plotVariables{varN}){subN, perturbN}(:, internalConN), 'color', colorCons(internalConN, :))
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
                    title([names{subN}, ', perturbPhase', num2str(perturbPhase(perturbN))])
                end
            end
            if varN>=saccadeVarStart
                saveas(gcf, [saccadeFolder, 'individuals\exp2\sac_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            elseif varN==1
                saveas(gcf, [perceptFolder, 'individuals\exp2\', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            else
                saveas(gcf, [pursuitFolder, 'individuals\exp2\pursuit_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
            end
            
            % plot difference between the two internal cons, compare
            % between perturbation phases
%             figure
%             hold on
%             % line plot
%             for internalConN = 1:size(yMeanDiffSub.(plotVariables{varN}){subN}, 2)
%                 errorbar(apertureAngles, yMeanDiffSub.(plotVariables{varN}){subN}(:, internalConN), yStdDiffSub.(plotVariables{varN}){subN}(:, internalConN), 'color', colorCons(internalConN+1, :))
%             end
%             
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
%             xlabel('Aperture angle (deg)')
%             legend(internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
%             
%             ylabel(['Diff ', plotVariables{varN}])
%             title(names{subN})
%             if varN>=saccadeVarStart
%                 saveas(gcf, [saccadeFolder, '\individuals\sacDiff_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
%             elseif varN==1
%                 saveas(gcf, [perceptFolder, 'individuals\diff_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
%             else
%                 saveas(gcf, [pursuitFolder, '\individuals\pursuitDiff_', plotVariables{varN}, '_lineplot_', names{subN}, '.pdf'])
%             end
            
            %             close all
        end
    end
end

%%
if averagePlots
    for varN = plotVarStart:plotVarEnd
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
                meanDiffAll(angleN) = nanmean(summaryDataDiff.(plotVariables{varN})(idxT, 1));
                stdDiffAll(angleN) = nanstd(summaryDataDiff.(plotVariables{varN})(idxT, 1));
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
            meanAll = NaN(size(apertureAngles));
            stdAll =NaN(size(apertureAngles));
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
                meanAll(angleN) = nanmean(summaryData.(plotVariables{varN})(idxT, 1));
                stdAll(angleN) = nanstd(summaryData.(plotVariables{varN})(idxT, 1));
            end
            errorbar(apertureAngles, meanAll, stdAll, 'color', colorCons(internalConN, :))
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