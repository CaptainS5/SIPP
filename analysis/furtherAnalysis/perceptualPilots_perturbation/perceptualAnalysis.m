% plot perceptual results for the perturbation pilot... x02
% calculate psychometric curves for the three internal motion conditions
initializeParas;

% choose which plot to look at now
individualPlots = 1;
averagedPlots = 0;

% plot settings
textFontSize = 8;

% make -1=choice of down
idxT = find(eyeTrialData.choice==-1); % 
eyeTrialData.choice(idxT) = 0;

% % make rdkApertureAngle left/right flipped--only the up/down values
% % relative to the aperture direction
% idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkApertureDir==180); % leftward valid trials
% eyeTrialData.rdkApertureAngle(idxT) = 180-eyeTrialData.rdkApertureAngle(idxT);
%
% % % flip the down aperture directions
% % idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkApertureAngle<0);
% % eyeTrialData.response(idxT) = -eyeTrialData.response(idxT);

% fitting settings
PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull,
%PAL_Quick, PAL_logQuick, PAL_Logistic
%PAL_CumulativeNormal, PAL_HyperbolicSecant

%Threshold, Slope, and lapse rate are free parameters, guess is fixed
paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter

%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = 0.01:.001:.11;
searchGrid.beta = logspace(0,3,101);
searchGrid.gamma = 0;  %scalar here (since fixed) but may be vector
searchGrid.lambda = 0:0.001:0.05;  %ditto

%% fitting data
for subN = 1:size(names, 2)
    data = table();
    idx = find(eyeTrialData.errorStatus(subN, :)==0);
    data.choice = eyeTrialData.choice(subN, idx)';
    data.rdkApertureDirPerturbation = eyeTrialData.rdkApertureDirPerturbation(subN, idx)';
    data.rdkInternalCons = eyeTrialData.rdkInternalPerturbationCons(subN, idx)';
    
    dataPercept.rdkInternalCons(subN, 1:length(internalCons)) = internalCons;
    
    figure
    hold on
    %% fitting for each coherence level
    for internalConN = 1:length(internalCons)
        dataT = data(data.rdkInternalCons==internalCons(internalConN), :);
        apertureDirLevels = unique(dataT.rdkApertureDirPerturbation); % stimulus levels, negative is down
        dataT.apertureDirIdx = zeros(size(dataT.rdkApertureDirPerturbation));
        for ii = 1:length(apertureDirLevels)
            dataT.apertureDirIdx(dataT.rdkApertureDirPerturbation==apertureDirLevels(ii), 1) = ii;
        end
        numUp{internalConN}(subN, :) = accumarray(dataT.apertureDirIdx, dataT.choice, [], @sum); % choice 1=up, 0=down
        outOfNum{internalConN}(subN, :) = accumarray(dataT.apertureDirIdx, dataT.choice, [], @numel); % total trial numbers
        
        %Perform fit
        [paramsValues{subN, internalConN} LL{subN, internalConN} exitflag{subN, internalConN}] = PAL_PFML_Fit(apertureDirLevels, numUp{internalConN}(subN, :)', ...
            outOfNum{internalConN}(subN, :)', searchGrid, paramsFree, PF);%, 'lapseLimits',[0 0.1]);
        
        % plotting
%         subplot(1, 3, internalConN)
        ProportionCorrectObserved=numUp{internalConN}(subN, :)./outOfNum{internalConN}(subN, :);
        StimLevelsFineGrain=[min(apertureDirLevels):max(apertureDirLevels)./1000:max(apertureDirLevels)];
        ProportionCorrectModel = PF(paramsValues{subN, internalConN},StimLevelsFineGrain);
        
        f{internalConN} = plot(StimLevelsFineGrain, ProportionCorrectModel,'-','color', colorCons(internalConN, :), 'linewidth', 2);
        plot(apertureDirLevels, ProportionCorrectObserved,'.', 'color', colorCons(internalConN, :), 'markersize', 30);
        
%         axis square
%         set(gca, 'fontsize', 16);
%         set(gca, 'Xtick', apertureDirLevels);
%         axis([min(apertureDirLevels) max(apertureDirLevels) 0 1]);
%         xlabel('Stimulus Intensity');
%         ylabel('Proportion up');
%         title(internalConNames{internalConN})
        
        % saving parameters
        dataPercept.alpha(subN, internalConN) = paramsValues{subN, internalConN}(1); % threshold, or PSE
        dataPercept.beta(subN, internalConN) = paramsValues{subN, internalConN}(2); % slope
        dataPercept.gamma(subN, internalConN) = paramsValues{subN, internalConN}(3); % guess rate, or baseline
        dataPercept.lambda(subN, internalConN) = paramsValues{subN, internalConN}(4); % lapse rate
    end
    set(gca, 'fontsize', 16);
    set(gca, 'Xtick', apertureDirLevels);
    axis([min(apertureDirLevels) max(apertureDirLevels) 0 1]);
    xlabel('Stimulus Intensity');
    ylabel('Proportion up');
    legend([f{:}], internalConNames, 'box', 'off', 'location', 'northwest')
    
    saveas(gcf, [perceptFolder, '\pf_', names{subN}, '.pdf'])
end

%% save csv for ANOVA
% cd(perceptFolder)
% cd ..
% cd(['R\Exp' num2str(expN)])
% 
% data = table();
% count = 1;
% for subN = 1:length(names)
%     for probNmerged = 1:probTotalN
%         data.sub(count, 1) = subN;
%         data.prob(count, 1) = probCons(probNmerged+(probTotalN-1));
%         data.PSE(count, 1) = dataPercept.alpha(subN, probNmerged);
%         data.slope(count, 1) = dataPercept.beta(subN, probNmerged);
%         count = count+1;
%     end
% end
% writetable(data, ['dataPercept_Exp' num2str(expN) '.csv'])