function [const] = constConfig(screen, const, sbj)
% =========================================================================
% constConfig(screen, const)
% =========================================================================
% All display lengths are in degree of visual angle (dva), time in s
%
% -------------------------------------------------------------------------
% Input:
% screen:    strucrure containing screen settings
% const:     structure containing different constant settings
% -------------------------------------------------------------------------
% Output:
% const:     structure containing different constant settings
% -------------------------------------------------------------------------
%%
if sbj.block == 1 && sbj.trial==1
    % Some dsign-related things (These will be used in paramConfig):
    if const.startExp==1 || const.startExp==0
%         if const.internalOnsetType==1 % constant internal motion
%             const.numTrialsPerBlock    = 36*ones(1, 10);                                % Each column = number of trials in one block; number of columns = number of blocks
%         elseif const.internalOnsetType==2 % perturbation of internal motion
            const.numTrialsPerBlock    = 42*ones(1, 10);                                % Each column = number of trials in one block; number of columns = number of blocks
%         end
    elseif const.startExp==-1
        const.numTrialsPerBlock    = 32*ones(1, 10);
    end
    if const.makeVideo; const.numTrialsPerBlock = 1; end
    const.numTrials            = sum(const.numTrialsPerBlock);                  % total number of trials
    
    %% Fixation:
    % const.fixationRadiusEyeVA    = 2.5;                                         % This is radius threshold for eye fixation (in pix?)
    % const.fixationRadiusEyePX    = round(screen.ppd*const.fixationRadiusEyeVA);
    
    %% Stimuli and Timing:
    % fixation
    const.fixation.durationMin = 1;
    const.fixation.durationMax = 1.5;
    const.fixation.dotRadius = 0.15; % in dva
    [const.fixation.dotRadiusPxl, ] = dva2pxl(const.fixation.dotRadius, const.fixation.dotRadius, screen); % in pixel
    const.fixation.colour = [255 255 255]; % fixation colour for fast trials
    const.fixation.windowRadius = 1; % tolerance window for fixation check, deg
    [const.fixation.windowRadiusPxl, ] = dva2pxl(const.fixation.windowRadius, const.fixation.windowRadius, screen); % in pixel
    
    % RDK stimulus
    const.rdk.durationBeforeMin = 0.5; % minimum display duration of the whole RDK, s
    const.rdk.durationBeforeMax = 0.7;
    const.rdk.durationWhole = 1.2;
    const.rdk.durationPerturbation = 0.2; % display duration of the perturbation period, s
    
    const.rdk.dotDensity = 10; % dot per dva^2
    const.rdk.lifeTime = 10; % s; longer than the whold display duration equals to unlimited lifetime
    % how long before a dot disappears and reappears
    const.rdk.labelUpdateTime = 0.050; % change labels and assign new directions for all
    % for Transparent motion, label update time >= the whole rdk duration;
    % for Brownian motion, label update time = one frame
    % for RDK with Gaussian distributed directions, label update time equals
    % the time of one direction for each dot; set to .14 for piloting
    % ========================================================================
    % Notice that whenever label changes, new directions will also be assigned,
    % so be careful with the relationship between lifeTime and labelUpdateTime
    % ========================================================================
    const.rdk.dotRadius = 0.05;
    const.rdk.apertureRadius = 1;
    const.rdk.dotFieldRadius = const.rdk.apertureRadius;
    const.rdk.apertureSpeed = 10; % dva per sec
    const.rdk.colour = screen.white;
    const.rdk.dotNumber = round(const.rdk.dotDensity*pi*const.rdk.dotFieldRadius^2);
    const.rdk.apertureDirBefore = [0]; % left (180) and right (0) 
    const.rdk.apertureDirPerturbation = [-21 -14 -7 0 7 14 21]; % relative to the before perturbation aperture direction  
    % directions are defined as the polar angle in degs away (clockwise is negative) from horizontal right; 
    % this is the initial angle of the perturbation; the second half of the
    % perturbation the RDK will return to the original horizontal level
    % before perturbation
    const.rdk.internalSpeed = 5; % speed of each internal dot
    % internal dots during perturbation: (coherence is 0 before perturbation)
    %     const.rdk.cohPerturbation = [1]; % coherence during the perturbation
    %     const.rdk.internalDirPerturbation = [-90 90]; % relative to the aperture perturbation direction
        const.rdk.internalPerturbationCons = [0, -90, 90]; % coh 0, or conherently vertical down/up
    
%     % for debugging...
%     const.rdk.apertureDirPerturbation = [21];
%     const.rdk.internalPerturbationCons = [-90 90];
    
    % warning beep for feedback on fixation maintainance
    const.beep.samplingRate = 44100;
    const.beep.sound = 0.9 * MakeBeep(300, 0.1, const.beep.samplingRate);
    
    % text size
    const.textSize = 25;
    const.textColour = screen.black;
    
    % inter-trial interval
    const.ITI = 0.05;
    
    %% Photodiode
    const.photoStimSizePX      = 50;                                            % in PX, might have to be adjusted (should not be visible for subjects)
    const.photoStimSizePX2      = 100;
    
    %% DPI Setup Stimulus
    const.calibtationOutRadius     = 0.225;                                     % this is in Degree VA
    const.calibtationInRadius      = 0.075;                                     % this is in Degree VA
    
    %% DPI Calibration Stimuli:
    const.calibrationHorFarVA      = 8; % deg
    % const.calibrationHorFarPX      = round(screen.ppcX*tan(const.calibrationHorFarVA/180*pi)*screen.dist);
    const.calibrationHorFarPX      = round(screen.ppd*const.calibrationHorFarVA);
    const.calibrationHorNearVA     = 4;
    % const.calibrationHorNearPX     = round(screen.ppcX*tan(const.calibrationHorNearVA/180*pi)*screen.dist);
    const.calibrationHorNearPX     = round(screen.ppd*const.calibrationHorNearVA);
    const.calibrationVertFarVA     = 8;
    % const.calibrationVertFarPX     = round(screen.ppcY*tan(const.calibrationVertFarVA/180*pi)*screen.dist);
    const.calibrationVertFarPX     = round(screen.ppd*const.calibrationVertFarVA);
    const.calibrationVertNearVA    = 4;
    % const.calibrationVertNearPX    = round(screen.ppcY*tan(const.calibrationVertNearVA/180*pi)*screen.dist);
    const.calibrationVertNearPX    = round(screen.ppd*const.calibrationVertNearVA);
    
    const.calibPositionsLeft       = [screen.x_mid-const.calibrationHorFarPX   screen.y_mid];
    const.calibPositionsHalfLeft   = [screen.x_mid-const.calibrationHorNearPX  screen.y_mid];
    const.calibPositionsCenter     = [screen.x_mid                             screen.y_mid];
    const.calibPositionsHalfRight  = [screen.x_mid+const.calibrationHorNearPX  screen.y_mid];
    const.calibPositionsRight      = [screen.x_mid+const.calibrationHorFarPX   screen.y_mid];
    
    const.calibPositionsTop        = [screen.x_mid  screen.y_mid-const.calibrationVertFarPX];
    const.calibPositionsHalfTop    = [screen.x_mid  screen.y_mid-const.calibrationVertNearPX];
    const.calibPositionsHalfBottom = [screen.x_mid  screen.y_mid+const.calibrationVertNearPX];
    const.calibPositionsBottom     = [screen.x_mid  screen.y_mid+const.calibrationVertFarPX];
    
    const.calibPositions           = {const.calibPositionsCenter, ...
        const.calibPositionsLeft,       const.calibPositionsHalfLeft,  const.calibPositionsCenter, ...
        const.calibPositionsHalfRight,  const.calibPositionsRight,     const.calibPositionsCenter, ...
        const.calibPositionsTop,        const.calibPositionsHalfTop,   const.calibPositionsCenter, ...
        const.calibPositionsHalfBottom, const.calibPositionsBottom,    const.calibPositionsCenter};
else % just load from info_Experiment
    load([sbj.sbjFolder ,'/info_Experiment.mat'])
    const = Experiment.const;    
end
end

