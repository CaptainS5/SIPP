% Exp1
% trying to find the best window of pursuit that correlates with
% perception...
initializeParas;

%% set up the parameters
% let's start with a simple one... the window always starts from the last
% moment, and increase with steps
% wLengthShortest = 50; % frames before the RDK offset
% wLengthLongest = 550; % frames before the RDK offset--hmmm now I realized that we cut the last 100ms from the analysis... well should probably include those
% % but we do have contaminations of... blinks and saccades in anticipation
% % to the end of the RDK display...
% wStepSize = 10; % in frames, same as ms
% wStart = [wLengthShortest:wStepSize:wLengthLongest]; % how many frames before the RDK offset, all to loop through
wStart = [200:1:700]; % frames from RDK onset
endT = 788; % the latest time point of all possible windows, frames from RDK onset
wLength = [50:10:endT-min(wStart)];
[windowGridStart, windowGridLength] = meshgrid(wStart, wLength);

% the final parameter table to be fill in; each row is a participant, each
% column is a window length
% RMSE = NaN(length(names), length(wStart)); % the root mean square error for each fit
% R2adjusted = NaN(length(names), length(wStart)); % the root mean square error for each fit

%% looping through each window length and do the fitting
for subN = 1:length(names)
    tic
    R2adjusted{subN} = NaN(length(wLength), length(wStart));
    RMSE{subN} = NaN(length(wLength), length(wStart));
    LogLH{subN} = NaN(length(wLength), length(wStart));
    pval{subN} = NaN(length(wLength), length(wStart));
    Fstat{subN} = NaN(length(wLength), length(wStart));
    
    load(['eyeTrialDataSub_', names{subN}, '.mat'])
    idxS = find(eyeTrialData.errorStatus(subN, :)==0);
    responseSub = eyeTrialData.response(subN, idxS);
    
    for startI = 1:length(wStart)
        for lengthI = 1:length(wLength)
            if wLength(lengthI)+wStart(startI)-2<=endT
                pursuitDir = [];
                for trialN = 1:length(idxS)
                    rdkOff = eyeTrialData.frameLog.rdkOff(subN, idxS(trialN));
                    %             window = [rdkOff-wStart(wI)+1:rdkOff];
                    rdkOn = eyeTrialData.frameLog.rdkOn(subN, idxS(trialN));
                    
                    window = [rdkOn+wStart(startI)-1:rdkOn+wStart(startI)+wLength(lengthI)-2];
                    
                    xTemp = nansum(eyeTrialDataSub.trial{idxS(trialN)}.pursuit.dirVec(window, 1));
                    yTemp = nansum(eyeTrialDataSub.trial{idxS(trialN)}.pursuit.dirVec(window, 2));
                    % normalize
                    dirX = xTemp/sqrt(xTemp^2+yTemp^2);
                    dirY = yTemp/sqrt(xTemp^2+yTemp^2);
                    
                    pursuitDir(trialN, 1) = atan2(dirY, dirX)/pi*180;
                end
                mdl = fitlm(pursuitDir, responseSub);
                [pval{subN}(lengthI, startI),Fstat{subN}(lengthI, startI)] = coefTest(mdl);
                RMSE{subN}(lengthI, startI) = mdl.RMSE;
                LogLH{subN}(lengthI, startI) = mdl.LogLikelihood;
                R2adjusted{subN}(lengthI, startI) = mdl.Rsquared.Adjusted;
            else
                break;
            end
        end
%         disp([names{subN}, ', windows starting from ', num2str(wStart(startI)), ' ms finished'])
    end
    toc
    disp([names{subN}, ' finished'])
end
save('temporalCorrFit.mat', 'windowGridStart', 'windowGridLength', 'R2adjusted', 'RMSE', 'LogLH', 'pval', 'Fstat', '-v7.3')

%% plot R2adjusted
% % individual heatmap
% for subN = 1:length(names)
%     figure
%     heatmap(wStart, wLength, R2adjusted{subN});
%     title('R^2 adjusted')
%     xlabel('Window start from RDK onset (ms)')
%     ylabel('Window length (ms)')
%     saveas(gcf, [correlationFolder, 'regressionR2adjustedHeatmap_pursuitDirVSpercept_', names{subN}, '.pdf'])
% end

% for subN = 1:length(names)
%     plot(wStart, R2adjusted(subN, :));
% end

%% "dependency index"
% for now, just take the max value for each participant
% for subN = 1:length(names)
%     temp = R2adjusted(subN, :);
%     dIndex(subN, 1) = find(temp==max(temp));
% end

% figure
% hold on
% for gN = 1:3
%     scatter(gN/2*ones(length(subGroup{gN}), 1), dIndex(subGroup{gN}), 'MarkerEdgeColor', gN/3*[0.9 0.9 0.9])
% end
% xlim([0, 2])
% xlabel(['assimilation--contrast--neutral perceptual group'])
% ylabel('Time from RDK onset with largest R^2 adjusted')
% saveas(gcf, 'dependencyIndexTemp.jpg')