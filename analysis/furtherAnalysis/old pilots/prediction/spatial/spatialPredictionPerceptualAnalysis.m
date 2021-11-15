% plot perceptual results for the spatial prediction exp
initializeParas;

% choose which plot to look at now
individualPlots = 1;
averagedPlots = 0;

% plot settings
textFontSize = 8;
plotSub = {'510'};

% % make rdkApertureAngle left/right flipped--only the up/down values
% % relative to the aperture direction
% idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkApertureDir==180); % leftward valid trials
% eyeTrialData.rdkApertureAngle(idxT) = 180-eyeTrialData.rdkApertureAngle(idxT);

% % flip the down aperture directions
% idxT = find(eyeTrialData.errorStatus==0 & eyeTrialData.rdkApertureAngle<0); 
% eyeTrialData.response(idxT) = -eyeTrialData.response(idxT);

% make 0=choice of below
idxT = find(eyeTrialData.response==-1); % 
eyeTrialData.response(idxT) = 0;

% % fitting settings
% PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull,
% %PAL_Quick, PAL_logQuick, PAL_Logistic
% %PAL_CumulativeNormal, PAL_HyperbolicSecant
% 
% %Threshold, Slope, and lapse rate are free parameters, guess is fixed
% paramsFree = [1 1 0 1];  %1: free parameter, 0: fixed parameter
% 
% %Parameter grid defining parameter space through which to perform a
% %brute-force search for values to be used as initial guesses in iterative
% %parameter search.
% searchGrid.alpha = 0.01:.001:.11;
% searchGrid.beta = logspace(0,3,101);
% searchGrid.gamma = 0;  %scalar here (since fixed) but may be vector
% searchGrid.lambda = 0:0.001:0.05;  %ditto

%% fitting data
for subN = 1:size(names, 2)
    figure
    for angleN = 1:length(apertureAngle)
        data = table();
        idx = find(eyeTrialData.errorStatus(subN, :)==0 & eyeTrialData.rdkApertureAngle(subN, :)==apertureAngle(angleN));
        data.response = eyeTrialData.response(subN, idx)';
        data.shiftDis = eyeTrialData.shiftDir(subN, idx)'/90*0.25;
        data.rdkInternalCons = eyeTrialData.rdkInternalCons(subN, idx)';
        
        dataPercept.rdkInternalCons(subN, 1:length(internalCons)) = internalCons;
        
        %     figure
        %     hold on
        subplot(2, 2, angleN)
        hold on
        %% fitting for each coherence level
        for internalConN = 1:length(internalCons)
            dataT = data(data.rdkInternalCons==internalCons(internalConN), :);
            stimLevels = unique(dataT.shiftDis); % stimulus levels, negative is down
            dataT.stimIdx = zeros(size(dataT.shiftDis));
            for ii = 1:length(stimLevels)
                dataT.stimIdx(dataT.shiftDis==stimLevels(ii), 1) = ii;
            end
            numUp{internalConN}(subN, :) = accumarray(dataT.stimIdx, dataT.response, [], @sum); % choice 1=up, 0=down
            outOfNum{internalConN}(subN, :) = accumarray(dataT.stimIdx, dataT.response, [], @numel); % total trial numbers
            
            %         %Perform fit
            %         [paramsValues{subN, internalConN} LL{subN, internalConN} exitflag{subN, internalConN}] = PAL_PFML_Fit(apertureSpeedLevels, numUp{internalConN}(subN, :)', ...
            %             outOfNum{internalConN}(subN, :)', searchGrid, paramsFree, PF);%, 'lapseLimits',[0 0.1]);
            %
            %         % plotting
            % %         subplot(1, 3, internalConN)
            ProportionCorrectObserved=numUp{internalConN}(subN, :)./outOfNum{internalConN}(subN, :);
            %         StimLevelsFineGrain=[min(apertureSpeedLevels):max(apertureSpeedLevels)./1000:max(apertureSpeedLevels)];
            %         ProportionCorrectModel = PF(paramsValues{subN, internalConN},StimLevelsFineGrain);
            %
            %         f{internalConN} = plot(StimLevelsFineGrain, ProportionCorrectModel,'-','color', colorCons(internalConN, :), 'linewidth', 2);
            %         plot(apertureSpeedLevels, ProportionCorrectObserved,'.', 'color', colorCons(internalConN, :), 'markersize', 30);
            
            f{internalConN} = plot(stimLevels, ProportionCorrectObserved,'-','color', colorCons(internalConN, :), 'linewidth', 2);
            
            %         axis square
            %         set(gca, 'fontsize', 16);
            %         set(gca, 'Xtick', apertureDirLevels);
            %         axis([min(apertureDirLevels) max(apertureDirLevels) 0 1]);
            %         xlabel('Stimulus Intensity');
            %         ylabel('Proportion up');
            %         title(internalConNames{internalConN})
            
            %         % saving parameters
            %         dataPercept.alpha(subN, internalConN) = paramsValues{subN, internalConN}(1); % threshold, or PSE
            %         dataPercept.beta(subN, internalConN) = paramsValues{subN, internalConN}(2); % slope
            %         dataPercept.gamma(subN, internalConN) = paramsValues{subN, internalConN}(3); % guess rate, or baseline
            %         dataPercept.lambda(subN, internalConN) = paramsValues{subN, internalConN}(4); % lapse rate
        end
        set(gca, 'fontsize', 16);
        set(gca, 'Xtick', stimLevels);
        axis([-.3 .3 0 1]);
        %     axis([min(stimLevels) max(stimLevels) 0 1]);
        xlabel('Shift distance (deg)');
        ylabel('Proportion above');
        title(['Apeture dir ', num2str(apertureAngle(angleN))])
%         legend([f{:}], internalConNames, 'box', 'off', 'location', 'northwest')
        if angleN==1
        legend([f{:}], internalConNames, 'box', 'off', 'location', 'southeast')
        end
    end
    saveas(gcf, [perceptFolder, '\perceptSeparate_', names{subN}, '.pdf'])
    %     saveas(gcf, [perceptFolder, '\pf_', names{subN}, '.pdf'])
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