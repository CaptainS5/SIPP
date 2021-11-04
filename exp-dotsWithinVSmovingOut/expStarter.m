%% ----------------------------- EXP STARTER -------------------------- %%
% ----------------------------------------------------------------------
% This script starts the experiment
% ----------------------------------------------------------------------
% written by Philipp KREYENMEIER (philipp.kreyenmeier@gmail.com)
% 07 / 10 / 2020
% adapted by Xiuyun Wu for micropursuit (xiuyunwu5@gmail.com)
% 13 / Apr / 2021
% Project : micropursuit
% Version : pilot, perturbation
% 
% ----------------------------------------------------------------------
% Description of Experiment
% 
% 
% 
% 
% ----------------------------------------------------------------------
% Description of Version
% 
% Translating RDK (moving aperture), with internal motion presented
% 
% ----------------------------------------------------------------------
% Last Changes:
% 
% 07/10/2020 (PK) - adapted for Koerner S257
% 13/Apr/2021 (XW)
% ----------------------------------------------------------------------
%% ---------------------------------------------------------------------



%% Initial Setup:
clear all; clear mex; clear functions; close all; home; sca;
pathConfig;

const.expName              = 'MicroPursuit'; 								% Experiment name
const.internalOnsetType    = 2;                                             % fixed; 1: constant internal motion (use the other version for this, not in these codes); 2: perturbation of internal motion 
const.startExp             = 0;                                             % 1 = experiment mode; 0 = debugging mode
const.expType              = 1;                                             % 1: experiment; --not implemented...-1: practice; 0: baseline
const.checkEyeFix          = 1;                                             % 1 = checks gaze fixation (this needs to be 1 also when in dummy mode)
% const.feedback             = 1;												% 1 = show task feedback (defined in runSingleTrial); 0 = off
const.makeVideo            = 0;                                             % 1 = creates a video of a single trial(set any conditions manually in expMain); 0 = off (normal experiment mode)
const.runScreenCalib       = 0;                                             % 1 = run screen calibration instead of experiment; 0 = experiment mode
const.showGaze             = 0;

% Eyelink Setup:
eyelink.mode               = 0;                                             % 1 = use eyelink; 0 = off
eyelink.dummy              = 0;                                             % 1 = eyelink in dummy mode; 0 = eyelink dummy off
eyelink.recalib            = true;                                          % true = recalibrate between blocks (recommanded); false = no calibration between blocks
eyelink.dummyEye           = [0,0];                                         % dummy start pos

% Photodiode Setup
photo.mode                 = 0;                                             % 1 = use Photodiode; 0 = off

% DPI Setup:
dpi_set.mode               = 0;                                             % 1 = use DPI; 0 = off
dpi_set.dummy              = 0;                                             % 1 = DPI in dummy mode; 0 = DPI dummy off
dpi_set.recalib            = true;                                          % true = recalibrate between blocks (recommanded); false = no calibration between blocks
dpi_set.dummyEye           = [0,0];                                         % dummy start pos
dpi_init(dpi_set, photo)                                                    % run DPI initialization 

%% Do some configurations (keys, screen, constants, trialData, sbj):
[sbj]                      = sbjConfig(const);
[keys] 					   = keyConfig;                                     % unify and define some keys
[screen]                   = screenConfig(const);                           % screen configurations (update if changes on setup or new setup used); opens PTB!
[const]                    = constConfig(screen, const, sbj);                    % set some constants and variables used in experiment
[trialCondition, trialData]                = paramConfig(const,sbj);        % changed from trialData to trialCondition
% now making the original copy of condition assignment as "trialCondition", 
% which remain untouched during the exp, and then record everything in trialData

% if use sound:
% [mysound]                  = soundConfig(const); 


%% Generate Target Trajectory and Other Stimuli: now doing it trial by trial to save memory
% [const]                    = generateTargetTrajectory(const,screen);        % Moving target trajectory
% [const]                    = generateStationaryStimuli(const,screen,...     % Stationary stimuli
%                                 trialData);       
% if use sound:
% [const]                    = generateSound(const,mysound,trialCondition);   % Stationary sound (at the moment)
   
%% Run the Experiment with defined Settings:
if ~const.runScreenCalib
    expMain(const, screen, keys, eyelink, dpi_set, photo, trialCondition, trialData, sbj);% run experiment
    
    % if use sound
    %expMain(const, screen, eyelink, mysound, trialCondition, sbj);% run experiment
    
    % run convert2ascSynch for the current subject, when done:
%     if IsWin
%         convert2ascSynch(sbj)
%     end
%% Alternatively, run Gamma Calibration procedure:
elseif const.runScreenCalib
    [screen]                = gammaCalib(screen, const, keys);              % run Screen Gamma Calibration
end
