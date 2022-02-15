% just some general stats for the publication
initializeParas;

% exclusion criteria in addition to blinking
% obvious error trials
idxT = find(eyeTrialData.errorStatus==0 & ...
    abs(eyeTrialData.rdkApertureAngle-eyeTrialData.response)>20 & ...
    eyeTrialData.rdkApertureAngle.*eyeTrialData.response<=0); % leftward valid trials;
eyeTrialData.errorStatus(idxT) = -2;

% trials with little eye movements
idxT = find(eyeTrialData.errorStatus==0 & ...
    eyeTrialData.pursuit.gainXexternal<0.5); 
eyeTrialData.errorStatus(idxT) = -3;

idxT = find(eyeTrialData.errorStatus==0 & ...
    eyeTrialData.pursuit.travelClpDis<eyeTrialData.pursuit.targetClpDis/2); 
eyeTrialData.errorStatus(idxT) = -3;

%% number of excluded eye trials
% for subN = 1:size(names, 2)
%     dataT = eyeTrialData.errorStatus(subN, :);
%     exclude1(subN, 1) = length(find(dataT==1 | isnan(dataT)));
%     excludePercent1(subN, 1) = exclude1(subN, 1)/length(dataT);
%     exclude2(subN, 1) = length(find(dataT==-2));
%     excludePercent2(subN, 1) = exclude2(subN, 1)/length(dataT);
%     exclude3(subN, 1) = length(find(dataT==-1));
%     excludePercent3(subN, 1) = exclude3(subN, 1)/length(dataT);
% end
% 
% disp([num2str(mean(excludePercent1)*100), ' +- ', num2str(std(excludePercent1)*100), '% blink trials excluded'])
% disp([num2str(mean(excludePercent2)*100), ' +- ', num2str(std(excludePercent2)*100), '% wrong perceptual report trials excluded'])
% disp([num2str(mean(excludePercent3)*100), ' +- ', num2str(std(excludePercent3)*100), '% little eye movement trials excluded'])

%% mean pursuit latency
for subN = 1:size(names, 2)
    idxT = find(eyeTrialData.errorStatus(subN, :)==0);    
    latency(subN, 1) = nanmean(eyeTrialData.pursuit.latency(subN, idxT));
end

disp(['mean pursuit latency is ' num2str(mean(latency)), ' +- ', num2str(std(latency)), ' ms'])
