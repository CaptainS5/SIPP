% generate the y-t retinal motion image for each participant in each trial
% need: rdk dot position, eye position, + time stamps
% generate images in different phases for a better comparison: before
% pursuit onset, initiation, early/late steady-state--would those be too
% short though?... have a look anyway
initializeParas;

% need to go into each folder to get the RDK trial positions
rawNames = {'lw0' 'ib1' 'tk' 'xw1' 'pd' 'cl' 'pw' 'mc' 'pk' 'yp' 'ts' 'cf' 'hl' 'qz' 'dc1' 'ja' 'mg' 'yz' 'lk' 'as'};
dataPath = [analysisFolder '\..\..\..\data\'];
addpath([analysisFolder '\..\..\..\exp1\Conversion\']) % conversion tools for regenerating the RDK

% note that the trial number of rdkseed corresponds to trialIdxInData in eventlog,
% not the trialN in eventlog; trialN in eventlog simply indicates the n_th
% trial in all valid trials (excluded the repeated trials in the exp)

%% loop into each trial to generate the image
tw = 100; % how long should the time window of the filter be, in ms
for subN = 1:size(names, 2)
    % define the sub folder
    currentSubject = rawNames{subN};
    currentSubjectPath = [dataPath, currentSubject];
    % load files
    load([currentSubjectPath, '\eventLog.mat']);
    load([currentSubjectPath, '\rdkFrameLog.mat']);
    load([currentSubjectPath, '\info_Experiment.mat']);
    load(['eyeTrialDataSub_', names{subN}, '.mat'])
    
    % calculate the average net motion energy in each condition
    
    for trialN = 1:size(eventLog, 1)
        tic
        if eyeTrialData.errorStatus(subN, trialN)==0 % only calculate the valid trials
            % load the rdk file
            originalIdx = eventLog.trialIdxInData(trialN);
            file = dir([currentSubjectPath, '\rdkSeeds\rdkseed_t', num2str(originalIdx, '%.3d'),'*.mat']);
            load([currentSubjectPath, '\rdkSeeds\', file.name])
            
            % need to regenerate dot position from seed...
            if Experiment.trialData.rdkInternalCons(originalIdx)==0
                control.rdkCoh          = 0;
                control.rdkInternalDir  = 0;
            else
                control.rdkCoh          = 1;
                control.rdkInternalDir  = Experiment.trialData.rdkInternalCons(originalIdx);
            end
            control.rdkApertureDir  = Experiment.trialData.rdkApertureDir(originalIdx);
            control.rdkApertureAngle  = Experiment.trialData.rdkApertureAngle(originalIdx);
            control.rdkInternalSpeed = Experiment.const.rdk.internalSpeed;
            
            dotPos = regenerateRDK(seed, Experiment.const, Experiment.screen, control);
            % since seed is recorded after randomizing the aperture center
            % position... only dotPos recreated is reliable; for dotPos, the 0,
            % 0 is the aperture center
            % center is the RDK aperture center, but down is positive;
            % the RDK center position should be read from rdkAperturePos from
            % the saved rdk_seed file; all in pixels,
            % and are in the PTB coordinate system, screen upper left corner
            % is (0, 0), right and down are positive
            dotPosX = [];
            dotPosY = [];
            rdkPosX = [];
            rdkPosY = [];
            
            % convert to screen center being (0, 0)/up and right positive
            for ii = 1:length(rdkAperturePos)
                dotPosX(:, ii) = dotPos{ii}(:, 1);
                dotPosY(:, ii) = -dotPos{ii}(:, 2); % still 0, 0 is aperture center, but down is negative now
                
                rdkPosX(1, ii) = rdkAperturePos{ii}(1)-Experiment.screen.x_mid;
                rdkPosY(1, ii) = Experiment.screen.y_mid-rdkAperturePos{ii}(2);
            end
            % these are in pixels, but not rounded... if want to make images
            % with 1 pixel as one cell, just round as eye position would be
            
%             figure
%             for ii = 1:length(rdkAperturePos)
%                 plot(dotPosX(:, ii), dotPosY(:, ii), '.')
% %                 plot(rdkPosX(:, ii), rdkPosY(:, ii), '.')
%                 xlim([-50, 50])
%                 ylim([-50, 50])
%                 pause
%             end
            
            % calculate the retinal position of the RDK...
            % first, align the timeline of eye position and dot position
            % for eye data, first frame is fixationOn, last frame is trialEnd
            rdkFrameStamps = rdkFrameLog{trialN}.eyeLinkTimeStamp-eventLog.fixationOn(trialN, 1)+1;
            eyePosXdeg = eyeTrialDataSub.trial{trialN}.eyeX_filt(rdkFrameStamps)';
            eyePosYdeg = eyeTrialDataSub.trial{trialN}.eyeY_filt(rdkFrameStamps)';
            % eye position is in deg, screen center is (0, 0), up an right is positive
            % need to convert to pixel, and make it the same coordinate system as in rdkPos...
            % screen upper left corner is 0, 0; center is Experiment.screen.x_mid/y_mid
            
            % eye position in pixels
            eyePosX = eyePosXdeg*Experiment.screen.ppd;
            eyePosY = eyePosYdeg*Experiment.screen.ppd;
            
            % retinal position of the rdk center in each frame
            retinalCenterPosX = rdkPosX-eyePosX;
            retinalCenterPosY = rdkPosY-eyePosY;
            
            % retinal position of the center of each dot, in pixel
            retinalDotPosX = dotPosX + repmat(retinalCenterPosX, size(dotPosX, 1), 1);
            retinalDotPosY = dotPosY + repmat(retinalCenterPosY, size(dotPosY, 1), 1);
            
            % pure image position on the screen of the center of each dot, in pixel
            screenDotPosX = dotPosX + repmat(rdkPosX, size(dotPosX, 1), 1);
            screenDotPosY = dotPosY + repmat(rdkPosY, size(dotPosY, 1), 1);
            
            
            % dot radius for plotting the dots on the image
            [rdkDotRadius, ]= dva2pxl(Experiment.const.rdk.dotRadius, Experiment.const.rdk.dotRadius, Experiment.screen);
            
            t = [1:length(retinalCenterPosY)]'; % in frames
            % how large should the x-y image be; column would be x position, row
            % would be y position; aligned on retinal positions
%             xMin = floor(min(retinalDotPosX(:)))-rdkDotRadius;
%             xMax = ceil(max(retinalDotPosX(:)))+rdkDotRadius;
%             yMin = floor(min(retinalDotPosY(:)))-rdkDotRadius;
%             yMax = ceil(max(retinalDotPosY(:)))+rdkDotRadius;
            
            % no eye movements, screen positions
            xMin = floor(min(screenDotPosX(:)))-rdkDotRadius;
            xMax = ceil(max(screenDotPosX(:)))+rdkDotRadius;
            yMin = floor(min(screenDotPosY(:)))-rdkDotRadius;
            yMax = ceil(max(screenDotPosY(:)))+rdkDotRadius;

% % pure internal motion relative to the aperture...
%             xMin = floor(min(dotPosX(:)))-rdkDotRadius;
%             xMax = ceil(max(dotPosX(:)))+rdkDotRadius;
%             yMin = floor(min(dotPosY(:)))-rdkDotRadius;
%             yMax = ceil(max(dotPosY(:)))+rdkDotRadius;

            [X Y] = meshgrid((xMin:xMax), (yMin:yMax));
            
            imageXYT = [];
            % for each frame, generate the x-y image first
            for frameN = 1:length(t)
                imageXYT(:, :, frameN) = ones(size(X))*255; % the x-y image to be filled in the dots
                
                for dotN = 1:size(retinalDotPosX, 1)
%                     % pure internal motion...
%                     idxXcenter = round(dotPosX(dotN, frameN))-xMin+1;
%                     idxYcenter = round(dotPosY(dotN, frameN))-yMin+1;
                    
                    % no eye movements involved, just pure image motion on
                    % the screen
                    idxXcenter = round(screenDotPosX(dotN, frameN))-xMin+1;
                    idxYcenter = round(screenDotPosY(dotN, frameN))-yMin+1;
                    
%                     % retinal position
%                     idxXcenter = round(retinalDotPosX(dotN, frameN))-xMin+1;
%                     idxYcenter = round(retinalDotPosY(dotN, frameN))-yMin+1;
                    
                    imageXYT(idxYcenter-rdkDotRadius:idxYcenter+rdkDotRadius, ...
                        idxXcenter-rdkDotRadius:idxXcenter+rdkDotRadius, frameN) = 0;
                end
%                 if control.rdkInternalDir~=0
%                     if frameN==1
%                         control.rdkInternalDir
%                     end
%                     imshow(imageXYT(:, :, frameN))
%                 end
            end
            
            % now let's do the motion energy analysis
%             if control.rdkInternalDir~=0
%                 trialN
%                 control.rdkInternalDir
                netME{subN}(:, trialN) = getMotionEnergy(imageXYT, tw, Experiment.screen); % net motion energy as a function of time
%                 plot(netME{subN}(:, trialN))
%             end
        else
            netME{subN}(1:59, trialN) = NaN;
        end
        subN
        trialN
        toc
    end
    
%     % have a look at whether it makes sense--net motion energy in different
%     % dot motion conditions
%     figure
%     hold on
%     for internalConN = 1:length(internalCons)
%         idx = find(eyeTrialData.errorStatus(subN, :)==0 ...
%             & eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN));
%         plot(nanmean(netME{subN}(:, idx), 2), 'color', colorDotCons(internalConN, :))
%     end
%     legend(internalConNames, 'box', 'off')
%     xlabel('Time frame')
%     ylabel('Net motion energy (vertical)')
%     title(names{subN})
%     saveas(gcf, ['netME_', names{subN}, '_', num2str(tw), 'ms.pdf'])
    
    % screen positions
    figure
    hold on
    for internalConN = 1:length(internalCons)
        idx = find(eyeTrialData.errorStatus(subN, :)==0 ...
            & eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN));
        plot(nanmean(netME{subN}(:, idx), 2), 'color', colorDotCons(internalConN, :))
    end
    legend(internalConNames, 'box', 'off')
    xlabel('Time frame')
    ylabel('Net motion energy (vertical)')
    title(names{subN})
    saveas(gcf, ['netME_pureMotionOnScreen_', names{subN}, '_', num2str(tw), 'ms.pdf'])
end
% save(['netMotionEnergy_', num2str(tw), 'ms.mat'], 'netME')
save(['netMotionEnergy_pureMotionOnScreen_', num2str(tw), 'ms.mat'], 'netME')

%% Now we can plot and analyze the net motion energy across participants
% load(['netMotionEnergy_', num2str(tw), 'ms.mat'])
% load('subGroupList.mat')
% allME = [];
% groupME = [];
% 
% for subN = 1:length(names)
%     for internalConN = 1:length(internalCons)
%         idx = find(eyeTrialData.errorStatus(subN, :)==0 ...
%             & eyeTrialData.rdkInternalCon(subN, :)==internalCons(internalConN));
%         allME{internalConN}(subN, :) = nanmean(netME{subN}(:, idx), 2);
%     end
% end
% 
% for groupN = 1:3
%     for internalConN = 1:length(internalCons)
%         groupME{groupN}(internalConN, :) = nanmean(allME{internalConN}(subGroup{groupN}, :));
%     end
% end
% 
% figure
% hold on
% for groupN = 1:3
%    plot(groupME{groupN}', 'color', colorGroup(groupN, :))
% end
% legend(groupNames, 'box', 'off')
% xlabel('Time frame')
% ylabel('Net motion energy (vertical)')
% saveas(gcf, ['motionEnergySubgroup_', num2str(tw), 'ms.pdf'])
