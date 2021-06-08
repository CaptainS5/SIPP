% correlation plots
initializeParas;

% choose which plot to look at
individualPlots = 1; % within-sub
averagePlots = 0; % across-sub

% input the eye parameters to plot against perception
eyeVariables = {'gainXexternal', 'gainYexternal', 'gainYaverage', 'gain2Dexternal', 'gain2Daverage', ...
    'dirGainExternal', 'dirClp', 'dirError', 'disCenterMean'};
saccadeVarStart = 100; % if there is no saccade variables, just use a super large number (larger than the number of all variables to be plotted)

% plot settings
textFontSize = 8;

load('summaryData')
load('summaryDataSub')

%%
close all
if individualPlots
    for subN = 1:size(names, 2)
%         % difference from baseline values, compared with perception
%         for varN = 7:7%length(eyeVariables)
%             figure
%             hold on
%             for internalConN = 1:size(yMeanDiffSub.(eyeVariables{varN}){subN}, 2)
%                 s{internalConN} = scatter(yMeanDiffSub.response{subN}(:, internalConN), yMeanDiffSub.(eyeVariables{varN}){subN}(:, internalConN), 'MarkerEdgeColor', colorCons(internalConN, :));
%             end
%             [rho, pval] = corr(yMeanDiffSub.response{subN}(:), yMeanDiffSub.(eyeVariables{varN}){subN}(:));
%             
%             xlabel('Perceptual response difference (deg)')
%             ylabel([eyeVariables{varN}, ' diff'])
%             legend([s{1:2}], internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
%             title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
%             if varN>=saccadeVarStart
%                 saveas(gcf, [correlationFolder, 'individuals\diff_sac_', eyeVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
%             else
%                 saveas(gcf, [correlationFolder, 'individuals\diff_pursuit_', eyeVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
%             end
%         end
        
        % absolute difference from baseline values, compared with perception
        for varN = 7:7%length(eyeVariables)
            figure
            hold on
            for internalConN = 1:size(yMeanDiffSub.(eyeVariables{varN}){subN}, 2)
                s{internalConN} = scatter(abs(yMeanDiffSub.response{subN}(:, internalConN)), abs(yMeanDiffSub.(eyeVariables{varN}){subN}(:, internalConN)), 'MarkerEdgeColor', colorCons(internalConN, :));
            end
            [rho, pval] = corr(abs(yMeanDiffSub.response{subN}(:)), abs(yMeanDiffSub.(eyeVariables{varN}){subN}(:)));
            
            xlabel('Abs perceptual response difference (deg)')
            ylabel(['Abs ', eyeVariables{varN}, ' diff'])
            legend([s{1:2}], internalConNames(2:3), 'box', 'on', 'location', 'best', 'color', 'w')
            title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
            if varN>=saccadeVarStart
                saveas(gcf, [correlationFolder, 'individuals\absDiff_sac_', eyeVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            else
                saveas(gcf, [correlationFolder, 'individuals\absDiff_pursuit_', eyeVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            end
        end

%         % 'raw' valuse in all conditions, compared with perception
%         for varN = 8:8%length(eyeVariables)
%             figure
%             hold on
%             for internalConN = 1:size(yMeanSub.(eyeVariables{varN}){subN}, 2)
%                 s{internalConN} = scatter(yMeanSub.response{subN}(:, internalConN)-apertureAngles', yMeanSub.(eyeVariables{varN}){subN}(:, internalConN), 'MarkerEdgeColor', colorCons(internalConN, :));
%             end
%             [rho, pval] = corr(yMeanSub.response{subN}(:)-[apertureAngles'; apertureAngles'; apertureAngles'], yMeanSub.(eyeVariables{varN}){subN}(:));
%             
%             xlabel('Perceptual response error (deg)')
%             ylabel(eyeVariables{varN})
%             legend([s{1:3}], internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
%             title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
%             if varN>=saccadeVarStart
%                 saveas(gcf, [correlationFolder, 'individuals\sac_', eyeVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
%             else
%                 saveas(gcf, [correlationFolder, 'individuals\pursuit_', eyeVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
%             end
%         end
%         
%         % comparison between eye measures, 'raw' values
%         for varN = 9:9%length(eyeVariables)
%             figure
%             hold on
%             for internalConN = 1:size(yMeanSub.(eyeVariables{varN}){subN}, 2)
%                 s{internalConN} = scatter(yMeanSub.dirError{subN}(:, internalConN), yMeanSub.(eyeVariables{varN}){subN}(:, internalConN), 'MarkerEdgeColor', colorCons(internalConN, :));
%             end
%             [rho, pval] = corr(yMeanSub.dirError{subN}(:), yMeanSub.(eyeVariables{varN}){subN}(:));
%             
%             xlabel('Eye direction error (deg)')
%             ylabel(eyeVariables{varN})
%             legend([s{1:3}], internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
%             title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
%             if varN>=saccadeVarStart
%                 saveas(gcf, [correlationFolder, 'individuals\sac_', eyeVariables{varN}, 'VSdirError_', names{subN}, '.pdf'])
%             else
%                 saveas(gcf, [correlationFolder, 'individuals\pursuit_', eyeVariables{varN}, 'VSdirError_', names{subN}, '.pdf'])
%             end
%         end
        
    end
end
%     close all