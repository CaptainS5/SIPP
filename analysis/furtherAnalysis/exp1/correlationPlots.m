% Exp1:
% correlation plots
initializeParas;
% load('subGroupList.mat')
% groupCons = unique(groupAll);

%% choose which plot to look at
individualPlots = 0; % within-sub
averagePlots = 1; % across-sub
plotVarStart = 6;
plotVarEnd = 6;

%%
close all
if individualPlots
    corrDiffStats = table;
    corrStats = table;
    count = 1;
    for varN = plotVarStart:plotVarEnd
        for subN = 1:size(names, 2)
            % trial by trial
            % bias vs. bias, only for trials with internal dot motion,
            % difference from the average in baseline individually
            corrE = [];
            corrP = [];
            s = {};
%             figure
%             hold on
            for internalConN = 2:size(internalCons, 2)
                dataEye = [];
                dataP = [];
                % find the eye trial data, all trials for the current
                % internal condition...
                if varN>=3 && varN<=openloopVarEnd
                    idxT = find(eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN) ...
                        & eyeTrialData.errorStatus(subN, :)==0 ...
                        & eyeTrialData.pursuit.onsetType(subN, :)==0);
                elseif varN>=3 && varN<saccadeVarStart
                    if strcmp(plotVariables{varN}, 'dirClp') ...
                            || strcmp(plotVariables{varN}, 'dirEarly') || strcmp(plotVariables{varN}, 'dirLate') ...
                            || strcmp(plotVariables{varN}, 'dirClpEarly') || strcmp(plotVariables{varN}, 'dirClpLate')
                        dirXName = [plotVariables{varN}, 'X'];
                        idxT = find(eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN) ...
                            & eyeTrialData.errorStatus(subN, :)==0 ...
                            & ~isnan(eyeTrialData.pursuit.(dirXName)(subN, :)));
                    else
                        idxT = find(eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN) ...
                            & eyeTrialData.errorStatus(subN, :)==0 ...
                            & ~isnan(eyeTrialData.pursuit.(plotVariables{varN})(subN, :)));
                    end
                else
                    idxT = find(eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN) & eyeTrialData.errorStatus(subN, :)==0);
                end
                
                if strcmp(plotVariables{varN}, 'dirClp') || strcmp(plotVariables{varN}, 'dirOlp') ...
                        || strcmp(plotVariables{varN}, 'dirEarly') || strcmp(plotVariables{varN}, 'dirLate') ...
                        || strcmp(plotVariables{varN}, 'dirClpEarly') || strcmp(plotVariables{varN}, 'dirClpLate')
                    dirXName = [plotVariables{varN}, 'X'];
                    dirYName = [plotVariables{varN}, 'Y'];
                    dir = atan2(eyeTrialData.pursuit.(dirYName)(subN, idxT), eyeTrialData.pursuit.(dirXName)(subN, idxT))/pi*180;
                    
                    % delete extreme values for Olp pursuit direction... or
                    % better, just do not analyze olp pursuit or really
                    % figure out the problem with finding pursuit onset
                    if strcmp(plotVariables{varN}, 'dirOlp')
                        dir(abs(dir)>90)=NaN; % has to be wrong, shouldn't be moving to the left...
                    end
                    
                    dataEye = dir;
                    %                 elseif strcmp(plotVariables{varN}, 'latency')
                else
                    if strcmp(plotVariables{varN}, 'response') % should never be this though...
                        dataEye = eyeTrialData.(plotVariables{varN})(subN, idxT);
                    elseif varN>=saccadeVarStart
                        dataEye = eyeTrialData.saccades.(plotVariables{varN})(subN, idxT);
                    else
                        dataEye = eyeTrialData.pursuit.(plotVariables{varN})(subN, idxT);
                    end
                end
                dataP = eyeTrialData.response(subN, idxT);
                
                % need to go trial by trial since we have to identify the
                % baseline and substract...
                for tN = 1:length(idxT)
                    % condition of the current trial
                    angle = eyeTrialData.rdkApertureAngle(subN, idxT(tN));
                    baseI = find(summaryData.sub==subN & summaryData.rdkApertureAngle==angle & summaryData.rdkInternalDir==0);
                    
                    % substract the baseline in summaryData
                    dataEye(tN) = dataEye(tN)-summaryData.(plotVariables{varN})(baseI); % eye data baseline
                    dataP(tN) = dataP(tN)-summaryData.response(baseI); % perception baseline
                end
                
                %                 if internalConN==1
                %                     corrE = y2;
                %                     corrP = eyeTrialData.response(subN, idxT);
                %                 else
                corrE = [corrE dataEye];
                corrP = [corrP dataP];
                %                 end
%                 s{internalConN-1} = scatter(dataEye, dataP, 'MarkerEdgeColor', colorDotCons(internalConN, :));
            end
            % exclude NaN trials...
            idxD = find(isnan(corrE));
            corrE(idxD) = [];
            corrP(idxD) = [];
            
            corrDiffStats.sub(count) = subN;
            corrDiffStats.varEye{count} = plotVariables{varN};
            [corrDiffStats.rho(count), corrDiffStats.p(count)] = corr(corrE', corrP');
            %
%             ylabel('Bias in perceived direction (deg)')
%             xlabel(['Bias in ' plotVariables{varN}])
%             legend([s{:}], internalConNames{2:3}, 'box', 'on', 'location', 'best', 'color', 'w')
%             title([names{subN}, ', r=', num2str(corrDiffStats.rho(count), '%.2f'), ' p=', num2str(corrDiffStats.p(count), '%.2f')])
%             if varN>=saccadeVarStart
%                 saveas(gcf, [correlationFolder, 'individuals\sacTrialBias_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
%             else
%                 saveas(gcf, [correlationFolder, 'individuals\pursuitTrialBias_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
%             end
            %
%             % raw vs. raw...
%             corrD = [];
%             corrR = [];
%             s = {};
%             figure
%             hold on
%             for internalConN = 1:size(internalCons, 2)
%                 if varN>=3 && varN<=openloopVarEnd
%                     idxT = find(eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN) ...
%                         & eyeTrialData.errorStatus(subN, :)==0 ...
%                         & eyeTrialData.pursuit.onsetType(subN, :)==0);
%                 elseif varN>=3 && varN<saccadeVarStart
%                     if strcmp(plotVariables{varN}, 'dirClp') ...
%                             || strcmp(plotVariables{varN}, 'dirEarly') || strcmp(plotVariables{varN}, 'dirLate') ...
%                             || strcmp(plotVariables{varN}, 'dirClpEarly') || strcmp(plotVariables{varN}, 'dirClpLate')
%                         dirXName = [plotVariables{varN}, 'X'];
%                         idxT = find(eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN) ...
%                             & eyeTrialData.errorStatus(subN, :)==0 ...
%                             & ~isnan(eyeTrialData.pursuit.(dirXName)(subN, :)));
%                     else
%                         idxT = find(eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN) ...
%                             & eyeTrialData.errorStatus(subN, :)==0 ...
%                             & ~isnan(eyeTrialData.pursuit.(plotVariables{varN})(subN, :)));
%                     end
%                 else
%                     idxT = find(eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN) & eyeTrialData.errorStatus(subN, :)==0);
%                 end
%                 if strcmp(plotVariables{varN}, 'dirClp') || strcmp(plotVariables{varN}, 'dirOlp') ...
%                         || strcmp(plotVariables{varN}, 'dirEarly') || strcmp(plotVariables{varN}, 'dirLate') ...
%                         || strcmp(plotVariables{varN}, 'dirClpEarly') || strcmp(plotVariables{varN}, 'dirClpLate')
%                     dirXName = [plotVariables{varN}, 'X'];
%                     dirYName = [plotVariables{varN}, 'Y'];
%                     dir = atan2(eyeTrialData.pursuit.(dirYName)(subN, idxT), eyeTrialData.pursuit.(dirXName)(subN, idxT))/pi*180;
%                     
%                     % delete extreme values for Olp pursuit direction... or
%                     % better, just do not analyze olp pursuit or really
%                     % figure out the problem with finding pursuit onset
%                     if strcmp(plotVariables{varN}, 'dirOlp')
%                         dir(abs(dir)>90)=NaN; % has to be wrong, shouldn't be moving to the left...
%                     end
%                     
%                     y2 = dir;
%                     %                 elseif strcmp(plotVariables{varN}, 'latency')
%                 else
%                     if strcmp(plotVariables{varN}, 'response')
%                         y2 = eyeTrialData.(plotVariables{varN})(subN, idxT);
%                     elseif varN>=saccadeVarStart
%                         y2 = eyeTrialData.saccades.(plotVariables{varN})(subN, idxT);
%                     else
%                         y2 = eyeTrialData.pursuit.(plotVariables{varN})(subN, idxT);
%                     end
%                 end
%                 if internalConN==1
%                     corrD = y2;
%                     corrR = eyeTrialData.response(subN, idxT);
%                 else
%                     corrD = [corrD y2];
%                     corrR = [corrR eyeTrialData.response(subN, idxT)];
%                 end
%                 s{internalConN} = scatter(y2, eyeTrialData.response(subN, idxT), 'MarkerEdgeColor', colorDotCons(internalConN, :));
%             end
%             idxD = find(isnan(corrD));
%             corrD(idxD) = [];
%             corrR(idxD) = [];
%             
%             corrStats.sub(count) = subN;
%             corrStats.varEye{count} = plotVariables{varN};
%             [corrStats.rho(count), corrStats.p(count)] = corr(corrD', corrR');
%             %             %
%             ylabel('Perceptual response (deg)')
%             xlabel([plotVariables{varN}])
%             legend([s{:}], internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
%             title([names{subN}, ', r=', num2str(corrStats.rho(count), '%.2f'), ' p=', num2str(corrStats.p(count), '%.2f')])
%             if varN>=saccadeVarStart
%                 saveas(gcf, [correlationFolder, 'individuals\sacTrial_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
%             else
%                 saveas(gcf, [correlationFolder, 'individuals\pursuitTrial_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
%             end
            
            % one point per angle x internalCon?
            %         % difference from baseline values, compared with perception
            %             figure
            %             hold on
            %             for internalConN = 1:size(yMeanDiffSub.(plotVariables{varN}){subN}, 2)
            %                 s{internalConN} = scatter(yMeanDiffSub.(plotVariables{varN}){subN}(:, internalConN), yMeanDiffSub.response{subN}(:, internalConN), 'MarkerEdgeColor', colorCons(internalConN, :));
            %             end
            %             [rho, pval] = corr(yMeanDiffSub.response{subN}(:), yMeanDiffSub.(plotVariables{varN}){subN}(:));
            %
            %             ylabel('Perceptual response difference (deg)')
            %             xlabel([plotVariables{varN}, ' diff'])
            %             legend([s{1:2}], internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
            %             title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
            %             if varN>=saccadeVarStart
            %                 saveas(gcf, [correlationFolder, 'individuals\diff_sac_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            %             else
            %                 saveas(gcf, [correlationFolder, 'individuals\diff_pursuit_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            %             end
            %         end
            
            %         % absolute difference from baseline values, compared with perception
            %         for varN = plotVarStart:plotVarEnd
            %             figure
            %             hold on
            %             % one point per condition
            %             for internalConN = 1:size(yMeanDiffSub.(plotVariables{varN}){subN}, 2)
            %                 s{internalConN} = scatter(abs(yMeanDiffSub.response{subN}(:, internalConN)), abs(yMeanDiffSub.(plotVariables{varN}){subN}(:, internalConN)), 'MarkerEdgeColor', colorCons(internalConN, :));
            %             end
            %             [rho, pval] = corr(abs(yMeanDiffSub.response{subN}(:)), abs(yMeanDiffSub.(plotVariables{varN}){subN}(:)));
            
            %             % one point per trial?
            %             idxT = find(strcmp(eyeTrialData.sub, names{subN}) & eyeTrialData.rdkInternalCon~=0 & eyeTrialData.errorStatus==0);
            %             absPerception = [];
            %             absEye = [];
            %             for ii = 1:idxT
            %                 apertureA = eyeTrialData.rdkApertureAngle(idxT(ii));
            %                 baseI = find(summaryData.sub==subN & summaryData.rdkInternalDir==0 & summaryData.rdkApertureAngle==apertureA);
            %
            %                 respT = ;
            %                 respB = summaryData.response(baseI);
            %                 eyeT = eyeTrialData.(plotVariables{varN})-----need to do a bit more to get here...;
            %                 eyeB = summaryData.(plotVariables{varN})(baseI);
            %
            %                 absPerception = [absPerception; abs(respT-respB)];
            %                 absEye = [absEye; abs(eyeT-eyeB)];
            %             end
            %             scatter(absPerception, absEye);
            %             [rho, pval] = corr(absPerception, absEye);
            
            %             xlabel('Abs perceptual response difference (deg)')
            %             ylabel(['Abs ', plotVariables{varN}, ' diff'])
            % %             legend([s{1:2}], internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
            %             title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
            %             if varN>=saccadeVarStart
            %                 saveas(gcf, [correlationFolder, 'individuals\absDiffTrial_sac_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            %             else
            %                 saveas(gcf, [correlationFolder, 'individuals\absDiffTrial_pursuit_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            %             end
            %         end
            
            % % one raw value, one difference from baseline value
            %         for varN = plotVarStart:plotVarEnd
            %             figure
            %             hold on
            %             for internalConN = 1:size(yMeanDiffSub.(plotVariables{varN}){subN}, 2)
            %                 s{internalConN} = scatter(yMeanSub.(plotVariables{varN}){subN}(:, internalConN+1), yMeanDiffSub.response{subN}(:, internalConN), 'MarkerEdgeColor', colorCons(internalConN, :));
            %             end
            %             rawValues = yMeanSub.(plotVariables{varN}){subN}(:, 2:3);
            %             [rho, pval] = corr(yMeanDiffSub.response{subN}(:), rawValues(:));
            %
            %             ylabel('Perceptual response difference (deg)')
            %             xlabel([plotVariables{varN}])
            %             legend([s{1:2}], internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
            %             title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
            %             if varN>=saccadeVarStart
            %                 saveas(gcf, [correlationFolder, 'individuals\rawVSdiff_sac_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            %             else
            %                 saveas(gcf, [correlationFolder, 'individuals\rawVSdiff_pursuit_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            %             end
            %         end
            
            %         % 'raw' valuse in all conditions, compared with perception
            %         for varN = plotVarStart:plotVarEnd
            %             figure
            %             hold on
            %             for internalConN = 1:size(yMeanSub.(plotVariables{varN}){subN}, 2)
            %                 s{internalConN} = scatter(yMeanSub.response{subN}(:, internalConN)-apertureAngles', yMeanSub.(plotVariables{varN}){subN}(:, internalConN), 'MarkerEdgeColor', colorCons(internalConN, :));
            %             end
            %             [rho, pval] = corr(yMeanSub.response{subN}(:)-[apertureAngles'; apertureAngles'; apertureAngles'], yMeanSub.(plotVariables{varN}){subN}(:));
            %
            %             xlabel('Perceptual response error (deg)')
            %             ylabel(plotVariables{varN})
            %             legend([s{1:3}], internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
            %             title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
            %             if varN>=saccadeVarStart
            %                 saveas(gcf, [correlationFolder, 'individuals\sac_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            %             else
            %                 saveas(gcf, [correlationFolder, 'individuals\pursuit_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            %             end
            %         end
            %
            %         % comparison between eye measures, 'raw' values
            %         for varN = plotVarStart:plotVarEnd
            %             figure
            %             hold on
            %             for internalConN = 1:size(yMeanSub.(plotVariables{varN}){subN}, 2)
            %                 s{internalConN} = scatter(yMeanSub.dirError{subN}(:, internalConN), yMeanSub.(plotVariables{varN}){subN}(:, internalConN), 'MarkerEdgeColor', colorCons(internalConN, :));
            %             end
            %             [rho, pval] = corr(yMeanSub.dirError{subN}(:), yMeanSub.(plotVariables{varN}){subN}(:));
            %
            %             xlabel('Eye direction error (deg)')
            %             ylabel(plotVariables{varN})
            %             legend([s{1:3}], internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
            %             title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
            %             if varN>=saccadeVarStart
            %                 saveas(gcf, [correlationFolder, 'individuals\sac_', plotVariables{varN}, 'VSdirError_', names{subN}, '.pdf'])
            %             else
            %                 saveas(gcf, [correlationFolder, 'individuals\pursuit_', plotVariables{varN}, 'VSdirError_', names{subN}, '.pdf'])
            %             end
            count = count+1;
        end
        %         close all
        
        % plot summary
%         idxCorr = find(strcmp(corrStats.varEye, plotVariables{varN}));
        idxCorrDiff = find(strcmp(corrDiffStats.varEye, plotVariables{varN}));
        
%         dataCorr = corrStats(idxCorr, :);
        dataCorrDiff = corrDiffStats(idxCorrDiff, :);
        
%         % raw vs. raw
%         figure
%         hold on
%         fill([0.5:3.5 fliplr(0.5:3.5)]', [-0.5*ones(4, 1); 0.5*ones(4, 1)], 0.95*[1 1 1], 'EdgeColor', 'none')
%         for groupN = 1:length(groupCons)
%             dataT = dataCorr(groupAll==groupCons(groupN), :);
%             
%             dataSig = dataT(dataT.p<0.05, :);
%             dataNonsig = dataT(dataT.p>=0.05, :);
%             
%             scatter(repmat(groupN, size(dataSig.rho)), dataSig.rho, 'jitter','on', 'jitterAmount',0.1, 'MarkerFaceColor', colorGroup(groupN, :), 'MarkerEdgeColor', colorGroup(groupN, :))
%             scatter(repmat(groupN, size(dataNonsig.rho)), dataNonsig.rho, 'jitter','on', 'jitterAmount',0.1, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colorGroup(groupN, :))
%         end
%         xlim([0.5, 3.5])
%         %     ylim([0, 1])
%         xlabel('Perceptual bias group')
%         ylabel('Correlation coefficient (r)')
%         title(['Raw ' plotVariables{varN} ' vs. perception'])
%         saveas(gcf, [correlationFolder, 'pursuitTrialSummary_' plotVariables{varN} 'vsPerception_all.pdf'])
        
%         % diff vs. diff
%         figure
%         hold on
%         fill([0.5:3.5 fliplr(0.5:3.5)]', [-0.5*ones(4, 1); 0.5*ones(4, 1)], 0.95*[1 1 1], 'EdgeColor', 'none')
%         for groupN = 1:length(groupCons)
%             dataT = dataCorrDiff(groupAll==groupCons(groupN), :);
%             
%             dataSig = dataT(dataT.p<0.05, :);
%             dataNonsig = dataT(dataT.p>=0.05, :);
%             
%             scatter(repmat(groupN, size(dataSig.rho)), dataSig.rho, 'jitter','on', 'jitterAmount',0.1, 'MarkerFaceColor', colorGroup(groupN, :), 'MarkerEdgeColor', colorGroup(groupN, :))
%             scatter(repmat(groupN, size(dataNonsig.rho)), dataNonsig.rho, 'jitter','on', 'jitterAmount',0.1, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colorGroup(groupN, :))
%         end
%         xlim([0.5, 3.5])
%         %     ylim([-0.8, 0.8])
%         xlabel('Perceptual bias group')
%         ylabel('Correlation coefficient (r)')
%         title(['Bias of ' plotVariables{varN} ' vs. perception'])
%         saveas(gcf, [correlationFolder, 'pursuitTrialBiasSummary_' plotVariables{varN} 'vsPerception_all.pdf'])

% output summary stats
disp([dataCorrDiff.varEye{1}, ', r = ', num2str(mean(abs(dataCorrDiff.rho))), ' +- ', num2str(std(abs(dataCorrDiff.rho))), ', ', num2str(length(find(dataCorrDiff.p<0.05))), ' sig p'])
    end
%     save('corrStatsIndividual.mat', 'corrStats', 'corrDiffStats')
end

%% t-test for individual r...
load('corrDiffStatsIndividual.mat')
phaseName = {'dirOlp', 'dirClpEarly', 'dirClpLate'};
for phaseN=1:3
    idx = find(strcmp(corrDiffStats.varEye, phaseName{phaseN}));
    [h(phaseN), p(phaseN), ci(:, phaseN), stats{phaseN}] = ttest(corrDiffStats.rho(idx));
    std(corrDiffStats.rho(idx))
    stats{phaseN}.sd
    cohensD(phaseN) = mean(corrDiffStats.rho(idx))/std(corrDiffStats.rho(idx));
end

%%
summaryDiffCopy = summaryDataDiff;
if averagePlots
    summaryDiffCopy.response(summaryDiffCopy.rdkInternalDir==-90) = -summaryDiffCopy.response(summaryDataDiff.rdkInternalDir==-90); % flip the responses, bias in the direction of the internal motion
    summaryDiffCopy.dirClp(summaryDiffCopy.rdkInternalDir==-90) = -summaryDiffCopy.dirClp(summaryDiffCopy.rdkInternalDir==-90);
    summaryDiffCopy.dirOlp(summaryDiffCopy.rdkInternalDir==-90) = -summaryDiffCopy.dirOlp(summaryDiffCopy.rdkInternalDir==-90);
    summaryDiffCopy.dirClpEarly(summaryDiffCopy.rdkInternalDir==-90) = -summaryDiffCopy.dirClpEarly(summaryDiffCopy.rdkInternalDir==-90);
    summaryDiffCopy.dirClpLate(summaryDiffCopy.rdkInternalDir==-90) = -summaryDiffCopy.dirClpLate(summaryDiffCopy.rdkInternalDir==-90);
    
    for varN = plotVarStart:plotVarEnd
        corrData = table;
        count = 1;
        
%         summaryDiffCopy.(plotVariables{varN})(summaryDiffCopy.rdkInternalDir==-90) = -summaryDiffCopy.(plotVariables{varN})(summaryDiffCopy.rdkInternalDir==-90);
        for subN = 1:size(names, 2) % first, get the summary data... each row is one sub, each column is one variable
            % one point per person...
%                         % saccades vs. magnitude of clp pursuit bias
%                         idxT = find(summaryDiffCopy.sub==subN);
%                         corrT = summaryDiffCopy(idxT, :);
%                         
%                         corrDT = [corrT.sumAmpYUp(corrT.rdkInternalDir==-90, 1); corrT.sumAmpYDown(corrT.rdkInternalDir==90, 1)];
% %                         corrDT = [corrT.numYUp(corrT.rdkInternalDir==-90, 1); corrT.numYDown(corrT.rdkInternalDir==90, 1)];
%                         corrData.(plotVariables{varN})(subN, 1) = nanmean(corrDT);
%                         
% %                         corrData.(plotVariables{varN})(subN, 1) = nanmean(summaryDiffCopy.(plotVariables{varN})(idxT));
%                         corrData.dirClp(subN, 1) = nanmean(summaryDiffCopy.dirClp(idxT));
%                         corrData.response(subN, 1) = nanmean(summaryDiffCopy.response(idxT));
            
            % diff vs. diff
            idxT = find(summaryDiffCopy.sub==subN);
            corrData.response(subN, 1) = nanmean(summaryDiffCopy.response(idxT));
            corrData.(plotVariables{varN})(subN, 1) = nanmean(summaryDiffCopy.(plotVariables{varN})(idxT));
            
%                         % raw vs. raw
%                         idxT = find(summaryData.sub==subN);
%                         corrData.response(subN, 1) = nanmean(summaryDataDiff.response(idxT));
%                         corrData.(plotVariables{varN})(subN, 1) = nanmean(summaryDataDiff.(plotVariables{varN})(idxT));
            %
            
            %             % raw vs. diff in response
            %             idxT = find(summaryData.sub==subN & ...
            %                 summaryData.rdkInternalDir ~= 0); % exclude the baseline condition
            %             corrData.(plotVariables{varN})(subN, 1) = abs(nanmean(summaryData.(plotVariables{varN})(idxT)));
            
            %             % one point per internalCon per person
            %             for internalN = 1:2
            %                 if internalN==1
            %                     idxT = find(summaryDataDiff.sub==subN & summaryDataDiff.rdkInternalDir>0);
            %                 else
            %                     idxT = find(summaryDataDiff.sub==subN & summaryDataDiff.rdkInternalDir<0);
            %                 end
            %                 corrData.response(count, 1) = nanmean(summaryDataDiff.response(idxT));
            %                 % diff vs. diff
            %                 corrData.(plotVariables{varN})(count, 1) = nanmean(summaryDataDiff.(plotVariables{varN})(idxT));
            %                 count = count+1;
            %             end
        end
        %         % one point per internalCon x angle per person
        %         corrData.response = summaryDataDiff.response;
                
                [rho, pval] = corr(corrData.(plotVariables{varN}), corrData.response);
                
                mycorr = @(x1,x2) corr(x1, x2);
                nIterations = 1000;
                [lower, upper] = bootci(nIterations,{mycorr,corrData.(plotVariables{varN}), corrData.response});
%         [rho, pval] = corr(corrData.(plotVariables{varN}), corrData.dirClp);
        
        figure
        hold on
        % all points together
        scatter(corrData.(plotVariables{varN}), corrData.response,...
                        'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
        lsline(gca)
        
        % colour subgroups
%         for groupN = 1:length(subGroup)
%             scatter(corrData.(plotVariables{varN})(subGroup{groupN}), corrData.response(subGroup{groupN}), ...
%                 'MarkerFaceColor', colorGroup(groupN, :), 'MarkerEdgeColor', colorGroup(groupN, :))
%             %             scatter(corrData.(plotVariables{varN})(subGroup{groupN}), corrData.dirClp(subGroup{groupN}), ...
%             %                 'MarkerFaceColor', colorGroup(groupN, :), 'MarkerEdgeColor', colorGroup(groupN, :))
%         end
        
%         xlabel(['Bias in ', plotVariables{varN}, ' opposite to dot motion'])
        %         xlabel(['Absolute ', plotVariables{varN}])
        xlabel(['Bias in ', plotVariables{varN}])
        
        %         ylabel('Bias in pursuit direction')
        ylabel('Bias in perceived direction')
        ylim([-8, 8])
        title(['r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
        
        if varN>=saccadeVarStart
                        saveas(gcf, [correlationFolder, 'sacDiff_', plotVariables{varN}, 'YVSperceptualBias_all.pdf'])
%             saveas(gcf, [correlationFolder, 'sacDiff_', plotVariables{varN}, 'VSpursuitBias_all.pdf'])
        else
            saveas(gcf, [correlationFolder, 'pursuitDiff_', plotVariables{varN}, 'VSperceptualBias_all.pdf'])
        end
        %         if varN>=saccadeVarStart
        %             saveas(gcf, [correlationFolder, 'sac_raw', plotVariables{varN}, 'VSperceptualBias_all.pdf'])
        %         else
        %             saveas(gcf, [correlationFolder, 'pursuit_raw', plotVariables{varN}, 'VSperceptualBias_all.pdf'])
        %         end
    end
end