% just checking perceptual results without eye data
clear all; close all; clc

names = {'w88'};
subStartI = 1;
internalCons = [0, -90, 90];
internalDir = [0, -1, 1];
internalConNames = {'static', 'down', 'up'};

cd ..
analysisPath = pwd; % folder for the eye movement preprocessing codes
dataPath = ['..\data\']; % still need to go into specific folders

% Psychometric function fitting settings
PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull,
%PAL_Quick, PAL_logQuick, PAL_Logistic
%PAL_CumulativeNormal, PAL_HyperbolicSecant

%Threshold, Slope, and lapse rate are free parameters, guess is fixed
paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter

%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = -3:.01:3;
searchGrid.beta = logspace(0,10,200);
searchGrid.gamma = 0:0.01:0.05;  %scalar here (since fixed) but may be vector
searchGrid.lambda = 0:0.001:0.05;  %ditto

colorCons = [0 0 0; 1 0 0; 0 0 1];

%% All trials
for subN = subStartI:length(names)
    % define the sub folder
    currentSubject = names{subN};
    currentSubjectPath = [dataPath, currentSubject];
    % load files
    load([currentSubjectPath, '\info_Experiment.mat'])
    
    clear trialDataSub
    
    for currentTrial = 1:size(Experiment.trialData, 1)
        trialData.sub{subN, currentTrial} = currentSubject;
        trialData.trialIdx(subN, currentTrial) = currentTrial;
        trialIdxInData = currentTrial;
        %         trialData.trialIdxInData(subN, currentTrial) = currentTrial; %trialIdxInData;
        trialData.blockN(subN, currentTrial) = Experiment.trialData.blockN(trialIdxInData, 1);
        
        %         trialData.rdkApertureDir(subN, currentTrial) = Experiment.trialData.rdkApertureDirBefore(trialIdxInData, 1); % direction before perturbation, right (0) or left (180)
        %         %         trialData.rdkApertureSpeed(subN, currentTrial) = trial.log.rdkApertureSpeed;
        %         trialData.rdkApertureAngle(subN, currentTrial) = Experiment.trialData.rdkApertureAnglePerturbation(trialIdxInData, 1); % angle during perturbation
        %         %         trialData.rdkInternalSpeed(subN, currentTrial) = Experiment.const.rdk.internalSpeed; % during perturbation
        %         trialData.rdkInternalCon(subN, currentTrial) = Experiment.trialData.rdkInternalPerturbationCons(trialIdxInData, 1);
        %         if trialData.rdkInternalCon(subN, currentTrial)==0
        %             trialData.rdkInternalDir(subN, currentTrial) = 0;
        %             trialData.rdkCoh(subN, currentTrial) = 0;
        %         else
        %             trialData.rdkInternalDir(subN, currentTrial) = Experiment.trialData.rdkInternalPerturbationCons(trialIdxInData, 1); % relative direction within the RDK
        %             trialData.rdkCoh(subN, currentTrial) = 1;
        %         end
        %         trialData.rdkInternalDir(subN, currentTrial) = trial.log.rdkInternalDir; % during perturbation
        %         trialData.rdkCoh(subN, currentTrial) = trial.log.rdkCoh;
        %         trialData.perturbPhase(subN, currentTrial) = trial.log.perturbTime;
        
        trialData.rdkApertureDir(subN, currentTrial) = Experiment.trialData.rdkApertureDir(trialIdxInData, 1); % direction before perturbation, right (0) or left (180)
        %         trialData.rdkApertureSpeed(subN, currentTrial) = trial.log.rdkApertureSpeed;
        trialData.rdkApertureAngle(subN, currentTrial) = Experiment.trialData.rdkApertureAngle(trialIdxInData, 1); % angle during perturbation
        %         trialData.rdkInternalSpeed(subN, currentTrial) = Experiment.const.rdk.internalSpeed; % during perturbation
        trialData.rdkInternalCon(subN, currentTrial) = Experiment.trialData.rdkInternalCons(trialIdxInData, 1);
        if trialData.rdkInternalCon(subN, currentTrial)==0
            trialData.rdkInternalDir(subN, currentTrial) = 0;
            trialData.rdkCoh(subN, currentTrial) = 0;
        else
            trialData.rdkInternalDir(subN, currentTrial) = Experiment.trialData.rdkInternalCons(trialIdxInData, 1); % relative direction within the RDK
            trialData.rdkCoh(subN, currentTrial) = 1;
        end
        trialData.response(subN, currentTrial) = Experiment.trialData.choice(trialIdxInData, 1);
    end
    
    % organizing data for plot
    apertureAngles = unique(trialData.rdkApertureAngle(subN, :));
    for internalConN = 1:length(internalCons)
        for angleN = 1:length(apertureAngles)
            idx = find(trialData.rdkInternalCon(subN, :)==internalCons(internalConN) & trialData.rdkApertureAngle(subN, :)==apertureAngles(angleN));
            dataUpN(internalConN, angleN) = length(find(trialData.response(subN, idx)==1));
            dataAllN(internalConN, angleN) = length(idx);
            dataS(internalConN, angleN) = nanmean(dataUpN(internalConN, angleN)/dataAllN(internalConN, angleN));
        end
    end
    
    % plotting
    figure
    hold on
    for internalConN = 1:length(internalCons)
        % line plot
        %         plot(apertureAngles, dataS(internalConN, :))
        
        % psychometric functions
        numRight = dataUpN(internalConN, :);
        outOfNum = dataAllN(internalConN, :);
        % Perform fit
        [paramsValues LL exitflag] = PAL_PFML_Fit(apertureAngles, numRight, ...
            outOfNum, searchGrid, paramsFree, PF, 'lapseLimits',[0 0.1]);
        
        % plotting
        ProportionCorrectObserved=numRight./outOfNum;
        StimLevelsFineGrain=[min(apertureAngles):max(apertureAngles)./1000:max(apertureAngles)];
        ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);
        
        if internalConN==1
            f{internalConN} = plot(StimLevelsFineGrain, ProportionCorrectModel,'-','color', colorCons(internalConN, :), 'linewidth', 2);
        else
            f{internalConN} = plot(StimLevelsFineGrain, ProportionCorrectModel,'--','color', colorCons(internalConN, :), 'linewidth', 2);
        end
        plot(apertureAngles, ProportionCorrectObserved,'.', 'color', colorCons(internalConN, :), 'markersize', 30);
    end

    set(gca, 'fontsize',16);
    set(gca, 'Xtick', apertureAngles);
    axis([min(apertureAngles) max(apertureAngles) 0 1]);
    xlabel('Aperture angle');
    ylabel('Proportion up');
    legend([f{:}], internalConNames, 'box', 'off', 'location', 'northwest')
%     title([names{subN}, ', perturbPhase', num2str(perturbPhase(perturbN))])
    
    % xlabel('Aperture angle (deg)')
    % legend(internalConNames, 'box', 'on', 'location', 'best', 'color', 'w')
    % ylabel('Proportion of perceiving up')
    saveas(gcf, [analysisPath, '\perceptualPlots\individuals\exp2\perceptPSE_', names{subN}, '.pdf'])
end