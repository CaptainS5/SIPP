% Exp1:
% correlation plots
initializeParas;

%% choose which plot to look at
individualPlots = 1; % within-sub
averagePlots = 0; % across-sub
plotVarStart = 11;
plotVarEnd = 13;

%%
close all
if individualPlots
    for varN = plotVarStart:plotVarEnd
        for subN = 1:size(names, 2)
            % trial by trial
            % raw vs. raw...
            corrD = [];
            corrR = [];
            figure
            hold on
            for internalConN = 1:size(internalCons, 2)
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
                    
%                     % delete extreme values for Olp pursuit direction... or
%                     % better, just do not analyze olp pursuit or really
%                     % figure out the problem with finding pursuit onset
%                     if strcmp(plotVariables{varN}, 'dirOlp')
%                         dir(abs(dir)>90)=[]; % has to be wrong, shouldn't be moving to the left...
%                     end
                    
                    y2 = dir;
%                 elseif strcmp(plotVariables{varN}, 'latency')
                else
                    if strcmp(plotVariables{varN}, 'response')
                        y2 = eyeTrialData.(plotVariables{varN})(subN, idxT);
                    elseif varN>=saccadeVarStart
                        y2 = eyeTrialData.saccades.(plotVariables{varN})(subN, idxT);
                    else
                        y2 = eyeTrialData.pursuit.(plotVariables{varN})(subN, idxT);
                    end
                end
                if internalConN==1
                    corrD = y2;
                    corrR = eyeTrialData.response(subN, idxT);
                else
                    corrD = [corrD y2];
                    corrR = [corrR eyeTrialData.response(subN, idxT)];
                end
                s{internalConN} = scatter(y2, eyeTrialData.response(subN, idxT), 'MarkerEdgeColor', colorCons(internalConN, :));
            end
            [rho, pval] = corr(corrD', corrR');
            %
            ylabel('Perceptual response (deg)')
            xlabel([plotVariables{varN}])
            legend([s{:}], internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
            title([names{subN}, ', r=', num2str(rho, '%.2f'), ' p=', num2str(pval, '%.2f')])
            if varN>=saccadeVarStart
                saveas(gcf, [correlationFolder, 'individuals\sacTrial_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            else
                saveas(gcf, [correlationFolder, 'individuals\pursuitTrial_', plotVariables{varN}, 'VSperception_', names{subN}, '.pdf'])
            end
            
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
        end
        close all
    end
end

%%
if averagePlots
    summaryDataDiff.response(summaryDataDiff.rdkInternalDir==-90) = -summaryDataDiff.response(summaryDataDiff.rdkInternalDir==-90); % flip the responses, bias in the direction of the internal motion
    for varN = plotVarStart:plotVarEnd
        corrData = table;
        count = 1;
        
        % one point per person...
        summaryDataDiff.(plotVariables{varN})(summaryDataDiff.rdkInternalDir==-90) = -summaryDataDiff.(plotVariables{varN})(summaryDataDiff.rdkInternalDir==-90);
        %         summaryData.(plotVariables{varN})(summaryDataDiff.rdkInternalDir==-90) = -summaryData.(plotVariables{varN})(summaryDataDiff.rdkInternalDir==-90);
        for subN = 1:size(names, 2) % first, get the summary data... each row is one sub, each column is one variable
            %             corrData.sub(subN, 1) = subN;
            % one point per person...
            idxT = find(summaryDataDiff.sub==subN);
            corrData.response(subN, 1) = nanmean(summaryDataDiff.response(idxT));
            % diff vs. diff
            corrData.(plotVariables{varN})(subN, 1) = nanmean(summaryDataDiff.(plotVariables{varN})(idxT));
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
        %
        %         % diff vs. diff
        %         corrData.(plotVariables{varN}) = summaryDataDiff.(plotVariables{varN});
        
        [rho, pval] = corr(corrData.(plotVariables{varN}), corrData.response);
        
        figure
        scatter(corrData.(plotVariables{varN}), corrData.response)
        xlabel(['Bias in ', plotVariables{varN}])
        %         xlabel(['Absolute ', plotVariables{varN}])
        ylabel('Bias in perceived direction')
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