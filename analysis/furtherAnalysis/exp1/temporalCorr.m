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
wStart = [200:1:550]; % frames from RDK onset

% the final parameter table to be fill in; each row is a participant, each
% column is a window length
RMSE = NaN(length(names), length(wStart)); % the root mean square error for each fit
R2adjusted = NaN(length(names), length(wStart)); % the root mean square error for each fit
pursuitDir = [];
responseSub = [];

%% looping through each window length and do the fitting
for subN = 1:length(names)
    load(['eyeTrialDataSub_', names{subN}, '.mat'])
    idxS = find(eyeTrialData.errorStatus(subN, :)==0);
    responseSub{subN} = eyeTrialData.response(subN, idxS);
    
    for wI = 1:length(wStart)
        for trialN = 1:length(idxS)
            rdkOff = eyeTrialData.frameLog.rdkOff(subN, idxS(trialN));
            %             window = [rdkOff-wStart(wI)+1:rdkOff];
            rdkOn = eyeTrialData.frameLog.rdkOn(subN, idxS(trialN));
            window = [rdkOn+wStart(wI)-1:rdkOn+wStart(wI)+199];
            
            xTemp = nansum(eyeTrialDataSub.trial{idxS(trialN)}.pursuit.dirVec(window, 1));
            yTemp = nansum(eyeTrialDataSub.trial{idxS(trialN)}.pursuit.dirVec(window, 2));
            % normalize
            dirX = xTemp/sqrt(xTemp^2+yTemp^2);
            dirY = yTemp/sqrt(xTemp^2+yTemp^2);
            
            pursuitDir{subN, wI}(trialN, 1) = atan2(dirY, dirX)/pi*180;
        end
        mdl{subN, wI} = fitlm(pursuitDir{subN, wI}, responseSub{subN});
        RMSE(subN, wI) = mdl{subN, wI}.RMSE;
        R2adjusted(subN, wI) = mdl{subN, wI}.Rsquared.Adjusted;
    end
end

%% plot RMSE
% figure
% hold on
% for subN = 1:length(names)
%     plot(wStart, R2adjusted(subN, :));
% end
% title('Window length 200ms')
% xlabel('Window start from RDK onset (ms)')
% ylabel('R^2 adjusted')
% saveas(gcf, 'regression_pursuitDirVSpercept_window200DifferentPos.pdf')

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