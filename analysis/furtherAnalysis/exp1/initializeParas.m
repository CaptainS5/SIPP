% Exp1
% use eyeTrialData to do analysis, initialize parameters for the predictive
% pursuit task

clear all; close all; clc

names = {'lw0' 'ib1' 'tk' 'xw1' 'pd' 'cl' 'pw' 'mc' 'pk' 'yp' 'ts' 'cf' 'hl' 'qz' 'dc1' 'ja' 'mg' 'yz' 'lk' 'as'}; % also for 500, 501, 504, ib1...
apertureAngles = [-9 -6 -3 0 3 6 9];
apertureAngleNames = {'-9', '-6', '-3', '0', '3', '6', '9'};
internalCons = [0, -90, 90];
internalConNames = {'coh 0', 'dir down', 'dir up'}; % 500 and 501
% internalConNames = {'static', 'dir down', 'dir up'}; % d00

% parameter settings for plots
plotVariables = {'response', 'turningPoint', ...
    'initialMeanVelocity2D', 'initialPeakVelocity2D', 'initialAccelerationFit2D', 'dirOlp', ...
    'dirEarly', 'dirLate', 'dirChange', ...
    'dirClp', 'dirClpEarly', 'dirClpLate', 'dirClpChange', 'dirError', 'disCenterMean', 'disCenterMeanEarly', 'disCenterMeanLate', ...
    'gainXexternal', 'gainYexternal', 'gainYaverage', 'gain2Dexternal', 'gain2Daverage', 'dirGainExternal', ...
    'num', 'numXLeft', 'numXRight', 'numYUp', 'numYDown', 'meanAmp2D', 'meanAmpXLeft', 'meanAmpXRight', 'meanAmpYUp', 'meanAmpYDown',...
    'sumAmp2D', 'sumAmpXLeft', 'sumAmpXRight', 'sumAmpYUp', 'sumAmpYDown'}; % always put all saccade parameters at last
openloopVarEnd = 7;
saccadeVarStart = 20; 

% names = {'w00' 'w01' 'w02' 'w03' 'w04' 'w05' 'w06' 'w07' 'w08' 'w09' 'w10'};
% cohCons = [0; 0.5; 1]; % RDK coherence'
% cohNames = {'coh-0', 'coh-50%', 'coh-100%'};
% internalDirCons = [-1 1]; % -1:left, 1-right
% internalDirNames = {'up', 'down'};

% names = {'xw0' 'dc0' 'ib0'};
% internalDirCons = [-45 -90 -135 45 90 135]; % -1:left, 1-right
% internalDirNames = {'-45', '-90', '-135', '45', '90', '135'};

% names = {'x02'};
% allCons.internalCons = [0, -90, 90];
% internalConNames = {'coh 0', 'dir down', 'dir up'};
% allCons.apertureAngles = [-21 -14 -7 0 7 14 21];
% apertureAngleNames = {'-21', '-14', '-7', '0', '7', '14', '21'};

sampleRate = 1000;

analysisFolder = pwd;
load('eyeTrialData_all.mat');
load('summaryData')
load('summaryDataDiff')
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