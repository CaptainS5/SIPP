% correlation plots
initializeParas;

% choose which plot to look at
individualPlots = 1; % within-sub
averagePlots = 0; % across-sub
plotVarStart = 12;
plotVarEnd = 12;

%%
close all
if individualPlots
    for subN = 1:size(names, 2)
%         % difference from baseline values, compared with perception
%         for varN = plotVarStart:plotVarEnd
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
        
        % absolute difference from baseline values, compared with perception
        for varN = plotVarStart:plotVarEnd
            figure
            hold on
%             % one point per condition
%             for internalConN = 1:size(yMeanDiffSub.(plotVariables{varN}){subN}, 2)
%                 s{internalConN} = scatter(abs(yMeanDiffSub.response{subN}(:, internalConN)), abs(yMeanDiffSub.(plotVariables{varN}){subN}(:, internalConN)), 'MarkerEdgeColor', colorCons(internalConN, :));
%             end
%             [rho, pval] = corr(abs(yMeanDiffSub.response{subN}(:)), abs(yMeanDiffSub.(plotVariables{varN}){subN}(:)));
             
%             % one point per trial
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
            
            xlabel('Abs perceptual response difference (deg)')
            ylabel(['Abs ', plotVariables{varN}, ' diff'])
%             legend([s{1:2}], internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
            title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
            if varN>=saccadeVarStart
                saveas(gcf, [correlationFolder, 'individuals\absDiffTrial_sac_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            else
                saveas(gcf, [correlationFolder, 'individuals\absDiffTrial_pursuit_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            end
        end

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
%         end
        
    end
end
%     close all

%%
if averagePlots % one point per person
    corrData = table;
    summaryDataDiff.response(summaryDataDiff.rdkInternalDir==-90) = -summaryDataDiff.response(summaryDataDiff.rdkInternalDir==-90); % flip the responses, bias in the direction of the internal motion
    
    for varN = plotVarStart:plotVarEnd
        summaryDataDiff.(plotVariables{varN})(summaryDataDiff.rdkInternalDir==-90) = -summaryDataDiff.(plotVariables{varN})(summaryDataDiff.rdkInternalDir==-90);
        summaryData.(plotVariables{varN})(summaryDataDiff.rdkInternalDir==-90) = -summaryData.(plotVariables{varN})(summaryDataDiff.rdkInternalDir==-90);
        
        for subN = 1:size(names, 2) % first, get the summary data... each row is one sub, each column is one variable
            corrData.sub(subN, 1) = subN;
            idxT = find(summaryDataDiff.sub==subN);
            corrData.response(subN, 1) = abs(nanmean(summaryDataDiff.response(idxT)));
            
            % diff vs. diff
            corrData.(plotVariables{varN})(subN, 1) = abs(nanmean(summaryDataDiff.(plotVariables{varN})(idxT)));
            
%             % raw vs. diff in response
%             idxT = find(summaryData.sub==subN & ...
%                 summaryData.rdkInternalDir ~= 0); % exclude the baseline condition
%             corrData.(plotVariables{varN})(subN, 1) = abs(nanmean(summaryData.(plotVariables{varN})(idxT)));
        end
        [rho, pval] = corr(corrData.(plotVariables{varN}), corrData.response);
        
        figure
        scatter(corrData.(plotVariables{varN}), corrData.response)
        xlabel(['Absolute bias in ', plotVariables{varN}])
%         xlabel(['Absolute ', plotVariables{varN}])
        ylabel('Absolute bias in perceived direction')
        title(['r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
        if varN>=saccadeVarStart
            saveas(gcf, [correlationFolder, 'sacDiff_', plotVariables{varN}, 'VSperceptualBias_all.pdf'])
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