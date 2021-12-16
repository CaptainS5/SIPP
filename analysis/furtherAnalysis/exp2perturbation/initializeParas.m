% Exp2
% use eyeTrialData to do analysis, initialize parameters for the predictive
% pursuit task

clear all; close all; clc

names = {'w01'}; % 
apertureAngles = [-9 -6 -3 0 3 6 9];
apertureAngleNames = {'-9', '-6', '-3', '0', '3', '6', '9'};
internalCons = [-90, 90];
internalConNames = {'dir down', 'dir up'}; 
perturbPhase = [1, 2]; % 1: perturbation during initiation; 2: perturbation during steady-state
perturbPhaseNames = {'during initiation', 'during steady-state'};

% parameter settings for plots
plotVariables = {'response', 'latency', ...
    'gainXexternal', 'gainYexternal', 'gainYaverage', 'gain2Dexternal', 'gain2Daverage', 'dirGainExternal', ...
    'dirPerturb', 'dirError', 'disCenterMean', ...
    'num', 'numXLeft', 'numXRight', 'numYUp', 'numYDown', 'meanAmp2D', 'meanAmpXLeft', 'meanAmpXRight', 'meanAmpYUp', 'meanAmpYDown',...
    'sumAmp2D', 'sumAmpXLeft', 'sumAmpXRight', 'sumAmpYUp', 'sumAmpYDown'}; % always put all saccade parameters at last
% openloopVarEnd = 6;
saccadeVarStart = 12; 

sampleRate = 1000;

analysisFolder = pwd;
load('eyeTrialData_all.mat');
load('summaryData')
% load('summaryDataDiff')
load('summaryDataSub')
perceptFolder = ['..\..\perceptPlots\'];
eyeTracesFolder = ['..\..\eyeTraces\'];
pursuitFolder = ['..\..\pursuitPlots\'];
saccadeFolder = ['..\..\saccadePlots\'];
correlationFolder = ['..\..\corrPlots\'];
perceptFolder = ['..\..\perceptualPlots\'];
RFolder = ['..\..\R\'];

% for plotting
textFontSize = 8;
for t = 1:size(names, 2) % individual color for scatter plots, can do 10 people
    if t<=2
        markerC(t, :) = (t+2)/4*[77 255 202]/255;
    elseif t<=4
        markerC(t, :) = (t)/4*[70 95 232]/255;
    elseif t<=6
        markerC(t, :) = (t-2)/4*[232 123 70]/255;
    elseif t<=8
        markerC(t, :) = (t-4)/4*[255 231 108]/255;
    elseif t<=10
        markerC(t, :) = (t-6)/4*[255 90 255]/255;
    end
end
colorCons = [0 0 0; 1 0 0; 0 0 1];
% colorProb = [8,48,107;198,219,239;8,48,107]/255; % all blue hues
% colorProb = [8,48,107;66,146,198;198,219,239;66,146,198;8,48,107]/255; % all blue hues
% colorPlot = [232 113 240; 15 204 255; 255 182 135; 137 126 255; 113 204 100]/255; % each row is one colour for one probability
colorPlot = [232 113 240; 15 204 255; 255 182 135; 0 0 0; 255 182 135; 15 204 255; 232 113 240]/255; % each row is one colour for one probability