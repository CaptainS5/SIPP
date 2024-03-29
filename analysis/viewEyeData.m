%% Script to manually look at each trial of all subjects
% this script requires selectSubject.m, analyzeTrial.m, plotResults.m
% always update experimental conficurations such as sampling rate distance
% to screen etc.
% you can optionally add a function to manually adjust saccades

% history
% 07-2012       JE created viewTrialByTrial.m
% 2012-2016     JF made edits
% 13-07-2018    JF commented to make the script more accecable for future
%               VPOM students
% for questions email jolande.fooken@rwth-aachen.de
% 09/25/2020    XW edited buttons and how error file is generated, based on Janick Edinger's torsion analysis code. xiuyunwu5@gmail.com

clear all; %clc

%% Define these parameters
currentTrial = 21; % chose trial you want to look at here; default = 1;
% as2, qz2
saccadeThreshold = 400; % acceleration
microSaccadeThreshold = 200;
% saccade algorithm threshold --> depends on your stimulus speed and
% expected saccade size
% note that this threshold is hard-coded! If you want to test different
% values this will not update while clicking through and you will have to
% declare the variable eagain in the command window

%% open a new figure
% (size depends on your current screen size)
name = 'click through eye movement data';
screenSize = get(groot,'ScreenSize');
close all;
fig = figure('Position', [25 50 screenSize(3)-100, screenSize(4)-150],'Name',name);

%% Subject selection
analysisPath = pwd;
% enter your data path here
cd ..
dataPath = fullfile(pwd,'data\');
cd(analysisPath);
currentSubjectPath = selectSubject(dataPath);

cd(currentSubjectPath);
load info_Experiment % load mat file containing experimental info
% delete the invalid trials which don't have the corresponding eye data
% from Experiment.trialData

load eventLog % variable matrix has all the event message frame indice
load rdkFrameLog
cd(analysisPath);

% get sub ID
sidx = strfind(currentSubjectPath, 'data\');
% "data\" should be the folder directly contain the sub data folder
currentSubject = currentSubjectPath(sidx+5:end);
% sidx + n, n depends on how many characters are in the string that you were finding

% prepare error file
errorFilePath = fullfile(analysisPath,'\ErrorFiles\');
if exist(errorFilePath, 'dir') == 0
    % Make folder if it does not exist.
    mkdir(errorFilePath);
end
errorFileName = [errorFilePath 'Sub_' currentSubject '_errorFile.mat'];
try
    load(errorFileName);
    disp('Error file loaded');
catch  %#ok<CTCH>
    errorStatus = NaN(size(eventLog, 1), 1);
    disp('No error file found. Created a new one.');
end

%% run analysis for each trial and plot
analyzeTrial;
plotResults;

buttons.discardTrial = uicontrol(fig,'string','!Discard trial!(1) >>','Position',[0,300,100,30],...
    'callback', 'errorStatus(currentTrial, 1)=1;currentTrial = currentTrial+1;analyzeTrial;plotResults;');

buttons.next = uicontrol(fig,'string','Next (0) >>','Position',[0,130,100,30],...
    'callback','errorStatus(currentTrial, 1)=0;currentTrial = min(currentTrial+1, size(eventLog, 1));analyzeTrial;plotResults;');
buttons.previous = uicontrol(fig,'string','<< Previous','Position',[0,100,100,30],...
    'callback','currentTrial = max(currentTrial-1,1);analyzeTrial;plotResults');
buttons.jumpToTrialn = uicontrol(fig,'string','Jump to trial..','Position',[0,70,100,30],...
    'callback','inputTrial = inputdlg(''Go to trial:'');currentTrial = str2num(inputTrial{:});analyzeTrial;plotResults;');

buttons.exitAndSave = uicontrol(fig,'string','Exit & Save','Position',[0,35,100,30],...
    'callback', 'close(fig);save(errorFileName,''errorStatus'');');