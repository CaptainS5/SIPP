% calculate the bootstrap mean and define subgroups
initializeParas;

nBoot = 1000;
%% bootstrapping
ciLow = [];
ciHigh = [];
baseMean = [];
conBiasMean = [];
conBiasStd = [];
angleBiasMean = [];

for subN = 1:length(names)
    for angleN = 1:length(apertureAngles)
%         idxBase = find(summaryData.sub==subN & summaryData.rdkApertureAngle==apertureAngles(angleN) & summaryData.rdkInternalDir==0);
%         baselineValue = summaryData.response(idxBase);
        for internalConN = 1:length(internalCons)
            idxD = find(eyeTrialData.rdkApertureAngle(subN, :)==apertureAngles(angleN) & eyeTrialData.rdkInternalDir(subN, :)==internalCons(internalConN) & eyeTrialData.errorStatus(subN, :)==0);
            if internalConN==1
                dTemp = eyeTrialData.response(subN, idxD);
                [ci, bootstat] = bootci(nBoot, @mean, dTemp);
                baseMean(subN, angleN) = bootstat(1);
            else
                dTemp = eyeTrialData.response(subN, idxD)-baseMean(subN, angleN);
                [ci, bootstat] = bootci(nBoot, @mean, dTemp);
                ciLow{subN}(angleN, internalConN-1) = ci(1);
                ciHigh{subN}(angleN, internalConN-1) = ci(2);
                conBiasMean{subN}(angleN, internalConN-1) = bootstat(1);
                conBiasStd{subN}(angleN, internalConN-1) = bootstat(2);
            end
        end
    end
    angleBiasMean(subN, :) = (conBiasMean{subN}(:, 2)-conBiasMean{subN}(:, 1))/2; % average values, merge internal directions
    
        
%     figure
%     hold on
%     for ii=1:2
%         errorbar(apertureAngles, conBiasMean{subN}(:, ii), conBiasStd{subN}(:, ii))
%     end
%     xlabel('Object motion direction (deg)')
%     ylabel('Perceptual bias (deg)')
%     legend(internalConNames{2:3}, 'box', 'off')
%     saveas(gcf, ['bootStrappedPerceptualBias_', names{subN}, '.pdf'])
end

%% grouping
% biasMean = [];
% biasStd = [];
% biasCI = [];
% for subN = 1:length(names)
%     biasMean(subN) = mean(angleBiasMean(subN, :));
%     biasStd(subN) = std(angleBiasMean(subN, :));
%     biasCI(subN, 1) = biasMean(subN)-1.96*biasStd(subN);
%     biasCI(subN, 2) = biasMean(subN)+1.96*biasStd(subN);
% end
[sortMean, idxS] = sort(biasMean);
% adding groupIdx
groupSortI = [];
for ii = 1:length(idxS)
    if ~isempty(find(subGroup{1}==idxS(ii)))
        groupSortI(ii) = 1;
    else ~isempty(find(subGroup{2}==idxS(ii)))
        groupSortI(ii) = 2;
    end
end

figure
hold on
% errorbar(sortMean, biasStd(idxS))
for groupN = 1:size(subGroup, 2) % color subgroups
    idxG = find(groupSortI==groupN);
    plot(idxG, sortMean(idxG), '+', 'color', colorGroup(groupN, :))
    errorbar(idxG, sortMean(idxG), biasStd(idxS(idxG)), 'lineStyle', 'none', 'color', colorGroup(groupN, :))
end
xlabel('Object motion direction (deg)')
ylabel('Perceptual bias (deg)')
title('meanDiffBias+-std')
line([1, 20], [0, 0])
% legend(internalConNames{2:3}, 'box', 'off')
saveas(gcf, ['bootStrappedPerceptualBiasDiffMean_all_std.pdf'])

%%
save('bootStrapPerceptualBiasMean.mat', 'biasMean', 'biasStd', 'biasCI', 'baseMean', 'conBiasMean', 'conBiasStd', 'angleBiasMean')

rangeLow = biasMean-biasStd;
rangeHigh = biasMean+biasStd;

subGroup{1} = find(rangeLow > 0);% assimilation group
subGroup{2} = find(rangeHigh < 0);% contrast group
subGroup{3} = find(rangeLow <= 0 & rangeHigh>=0);% undefined group
save('subGroupListBootStrap.mat', 'subGroup')

%% output to R for LPA
dataT = table;
dataT.sub(:, 1) = [1:20]';
dataT(:,2:8) = mat2cell(angleBiasMean, ones(size(angleBiasMean, 1), 1), ones(size(angleBiasMean, 2), 1));

writetable(dataT, [RFolder, 'bootStrapPerceptualBias.csv'])