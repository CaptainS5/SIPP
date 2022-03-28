% plot trace of pursuit bias.... hmmmm
initializeParas;

%%
% average direction vector for each observer in each internal condition
% aligned at RDK onset
frameLength = NaN(size(names, 2), 1);
latency = NaN(size(names));

dirXtrial = [];
dirYtrial = [];
baseXSub = [];
baseYSub = [];
biasXCon = [];
biasYCon = [];
biasXSub = [];
biasYSub = [];

baseDirSub = [];
biasDirCon = [];
biasDirSub = [];

for subN = 1:size(names, 2)
    cd(analysisFolder)
    load(['eyeTrialDataSub_' names{subN} '.mat']);
    
    idxT = find(eyeTrialData.errorStatus(subN, :)==0);
    tempL = eyeTrialData.frameLog.rdkOff(subN, idxT)-eyeTrialData.frameLog.rdkOn(subN, idxT)+1; % starting from RDK onset, to the end of RDK
    tempL(tempL==0) = [];
    frameLength(subN, 1) = min(tempL);
    % latency
    latency(subN) = nanmean(eyeTrialData.pursuit.onset(subN, idxT)-eyeTrialData.frameLog.rdkOn(subN, idxT));
    
    for internalConN = 1:length(internalCons)
        if internalCons(internalConN)==0
            rdkCoh = 0;
            rdkInternalDir = 0;
        else
            rdkCoh = 1;
            rdkInternalDir = internalCons(internalConN);
        end
        
        for angleN = 1:length(apertureAngles)
            idxT = find(eyeTrialData.errorStatus(subN, :)==0 & ...
                ~isnan(eyeTrialData.pursuit.onset(subN, :)) & ...
                eyeTrialData.rdkCoh(subN, :)==rdkCoh & ...
                eyeTrialData.rdkInternalDir(subN, :)==rdkInternalDir & ...
                eyeTrialData.rdkApertureAngle(subN, :)==apertureAngles(angleN));
            
            lengthT = length(idxT);
            dirXtrial{subN}{internalConN, angleN} = NaN(lengthT, frameLength(subN, 1));
            dirYtrial{subN}{internalConN, angleN} = NaN(lengthT, frameLength(subN, 1));
            
            for trialN = 1:lengthT
                % align at RDK onset
                startI = eyeTrialData.frameLog.rdkOn(subN, idxT(trialN));
                pursuitOnsetI = eyeTrialData.pursuit.onset(subN, idxT(trialN))-startI;
                endI = startI+frameLength(subN, 1)-1;
                dirXtrial{subN}{internalConN, angleN}(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.pursuit.dirVec(startI:endI, 1);
                dirYtrial{subN}{internalConN, angleN}(trialN, :) = eyeTrialDataSub.trial{1, idxT(trialN)}.pursuit.dirVec(startI:endI, 2);
                
                % exclude before pursuit onset
                dirXtrial{subN}{internalConN, angleN}(trialN, 1:pursuitOnsetI) = NaN;
                dirYtrial{subN}{internalConN, angleN}(trialN, 1:pursuitOnsetI) = NaN;
                
            end
            xTemp = nansum(dirXtrial{subN}{internalConN, angleN});
            yTemp = nansum(dirYtrial{subN}{internalConN, angleN});
            
            if internalCons(internalConN)==0 % baseline condition
                baseXSub{subN}(angleN, :) = xTemp./sqrt(xTemp.^2+yTemp.^2); % normalize to unit length vectors
                baseYSub{subN}(angleN, :) = yTemp./sqrt(xTemp.^2+yTemp.^2);
                
                baseDirSub{subN}(angleN, :) = atan2(baseYSub{subN}(angleN, :), baseXSub{subN}(angleN, :))/pi*180;
                % exclude extreme and likely wrong values
                baseDirSub{subN}(angleN, abs(baseDirSub{subN}(angleN, :))>90) = NaN;
            else
                biasXtemp = xTemp./sqrt(xTemp.^2+yTemp.^2);
                biasYtemp = yTemp./sqrt(xTemp.^2+yTemp.^2);
                
%                 biasXtemp = biasXtemp-baseXSub{subN}(angleN, :);
%                 biasYtemp = biasYtemp-baseYSub{subN}(angleN, :);
%                 
%                 biasXCon{subN, internalConN-1}(angleN, :) = biasXtemp./sqrt(biasXtemp.^2+biasYtemp.^2);
%                 biasYCon{subN, internalConN-1}(angleN, :) = biasYtemp./sqrt(biasXtemp.^2+biasYtemp.^2);
                
                biasDirCon{subN, internalConN-1}(angleN, :) = atan2(biasYtemp, biasXtemp)/pi*180-baseDirSub{subN}(angleN, :);
                % exclude extreme and likely wrong values
                biasDirCon{subN, internalConN-1}(angleN, abs(biasDirCon{subN, internalConN-1}(angleN, :))>90) = NaN;
            end
        end   
    end
    
    % averaged direction bias
%     xTemp = nansum([biasXCon{subN, 1}; biasXCon{subN, 2}]);
%     yTemp = nansum([-biasYCon{subN, 1}; biasYCon{subN, 2}]);
%     
%     biasXSub(subN, 1:length(xTemp)) = xTemp./sqrt(xTemp.^2+yTemp.^2);
%     biasYSub(subN, 1:length(xTemp)) = yTemp./sqrt(xTemp.^2+yTemp.^2);
    
    biasDirSub(subN, 1:length(biasDirCon{subN, 1})) = nanmean([-biasDirCon{subN, 1}; biasDirCon{subN, 2}]);
end

% average bias in each subgroup
load('subGroupListLPA.mat')
biasDirGroup = [];
biasUpper = [];
biasLower = [];
for subGroupN = 1:size(subGroup, 2)
    biasDirGroup(subGroupN, :) = nanmean(biasDirSub(subGroup{subGroupN}, :));
    biasUpper(subGroupN, :) = biasDirGroup(subGroupN, :)+tinv(0.975,length(subGroup{subGroupN})-1)*nanstd(biasDirSub(subGroup{subGroupN}, :))/sqrt(length(subGroup{subGroupN}));
    biasLower(subGroupN, :) = biasDirGroup(subGroupN, :)+tinv(0.025,length(subGroup{subGroupN})-1)*nanstd(biasDirSub(subGroup{subGroupN}, :))/sqrt(length(subGroup{subGroupN}));
end

%%
figure
hold on
for ii = 1:7
    if ii<4
        lineStyle = '--';
    else
        lineStyle = '-';
    end
    plot(biasDirCon{1, 2}(ii, :), 'lineStyle', lineStyle, 'color', colorObjAngles(ii, :))
end
ylim([-20, 20])

%%
% plotting parameters
startI = floor(mean(latency));

figure
hold on
tMS = [1:size(biasDirGroup, 2)]; % time in ms
p = [];
for drawN = 1:2
   groupN = 3-drawN;
   p{groupN} = plot(tMS(startI:end), biasDirGroup(groupN, startI:end), 'color', colorGroup(groupN, :));
   patch('XData', [tMS(startI:end), fliplr(tMS(startI:end))], 'YData', [biasUpper(groupN, startI:end), fliplr(biasLower(groupN, startI:end))], ...
       'faceColor', colorGroup(groupN, :), 'faceAlpha', 0.2, 'edgeColor', 'none', 'lineWidth', 5)
end
% add time stamps
plot([startI+140, startI+140], [-5, 20], '--')
plot([startI+140+(700-startI-140)/2, startI+140+(700-startI-140)/2], [-5, 20], '--')
plot([700, 700], [-5, 20], '--')

legend([p{:}], groupNames, 'box', 'off')
xlim([100, 800])
xlabel('Time (ms)')
ylabel('Bias in pursuit direction (deg)')
saveas(gcf, ['pursuitBiasTraceSubgroup.pdf'])
