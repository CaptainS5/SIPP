% just checking perceptual results without eye data
clear all; close all; clc

names = {'qz2'};
subStartI = 1;
durationCons = [0.25, 0.5]; % s
internalCons = [0, -90, 90];
internalDir = [0, -1, 1];
plotConNames = {'250ms-static', '250ms-down', '250ms-up', '500ms-static', '500ms-down', '500ms-up'};

cd ..
cd ..
analysisPath = pwd; % folder for the eye movement preprocessing codes
dataPath = ['..\data\']; % still need to go into specific folders

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
        
        trialData.rdkDuration(subN, currentTrial) = Experiment.trialData.rdkDuration(trialIdxInData, 1);
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
        trialData.response(subN, currentTrial) = Experiment.trialData.reportAngle(trialIdxInData, 1);
    end
    
    % organizing data for plot
    apertureAngles = unique(trialData.rdkApertureAngle(subN, :));
    for durationN = 1:length(durationCons)
        for internalConN = 1:length(internalCons)
            for angleN = 1:length(apertureAngles)
                idx = find(trialData.rdkInternalCon(subN, :)==internalCons(internalConN) ...
                    & trialData.rdkApertureAngle(subN, :)==apertureAngles(angleN) ...
                    & trialData.rdkDuration(subN, :)==durationCons(durationN));
                dataS{durationN}(internalConN, angleN) = nanmean(trialData.response(subN, idx));
            end
        end
    end
    
    % plotting
    figure       
    hold on
    for durationN = 1:length(durationCons)
%         subplot(1, 2, durationN)
        for internalConN = 1:length(internalCons)
            % line plot
            if durationN==1
                plot(apertureAngles, dataS{durationN}(internalConN, :), '--', 'color', colorCons(internalConN, :))
            else
                plot(apertureAngles, dataS{durationN}(internalConN, :), '-', 'color', colorCons(internalConN, :))
            end
        end
    end
    set(gca, 'fontsize',16);
    set(gca, 'Xtick', apertureAngles);
%     axis([min(apertureAngles) max(apertureAngles) -20 25]);
    xlabel('Aperture angle (deg)');
    ylabel('Reported angles (deg)');
    legend(plotConNames, 'box', 'off', 'location', 'northwest')
    title([names{subN}])
    saveas(gcf, [analysisPath, '\perceptualPlots\individuals\exp2temporal\percept_', names{subN}, '.pdf'])
end