% Exp1
% use eyeTrialData to do analysis, initialize parameters for the predictive
% pursuit task

clear all; close all; clc

% names = {'lw0' 'ib1' 'tk' 'xw1' 'pd' 'cl' 'pw' 'mc' 'pk' 'yp' 'ts' 'cf' 'hl' 'qz' 'dc1' 'ja' 'mg' 'yz' 'lk' 'as'}; 
names = {'p1' 'p2' 'p3' 'p4' 'p5' 'p6' 'p7' 'p8' 'p9' 'p10' 'p11' 'p12' 'p13' 'p14' 'p15' 'p16' 'p17' 'p18' 'p19' 'p20'}; 
apertureAngles = [-9 -6 -3 0 3 6 9];
apertureAngleNames = {'-9', '-6', '-3', '0', '3', '6', '9'};
internalCons = [0, -90, 90];
% internalConNames = {'coh 0', 'dir down', 'dir up'}; % 500 and 501
internalConNames = {'static', 'dir down', 'dir up'}; % d00
groupNames = {'Assimilation', 'Contrast', 'No bias'};

% parameter settings for plots
plotVariables = {'response', ...
    'initialMeanVelocity2D', 'initialPeakVelocity2D', 'initialAccelerationFit2D', 'dirOlp', ...
    'dirClp', 'dirClpEarly', 'dirClpLate', 'dirClpChange', 'dirError', 'disCenterMean', 'disCenterMeanEarly', 'disCenterMeanLate', ...
    'gainXexternal', 'gainYexternal', 'gainYaverage', 'gain2Dexternal', 'gain2Daverage', 'dirGainExternal', ...
    'meanPosErrOnsetX', 'meanPosErrOnsetY', 'meanPosErrOnset2D', 'meanPosErrOffsetX', 'meanPosErrOffsetY', 'meanPosErrOffset2D', ...
    'num', 'numXLeft', 'numXRight', 'numYUp', 'numYDown', 'meanAmp2D', ...
    'meanAmpXLeft', 'meanAmpXRight', 'meanAmpYUp', 'meanAmpYDown',...
    'sumAmp2D', 'sumAmpXLeft', 'sumAmpXRight', 'sumAmpYUp', 'sumAmpYDown'}; % always put all saccade parameters at last
openloopVarEnd = 5;
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
% delete obvious error trials
idxT = find(eyeTrialData.errorStatus==0 & ...
    abs(eyeTrialData.rdkApertureAngle-eyeTrialData.response)>20 & ...
    eyeTrialData.rdkApertureAngle.*eyeTrialData.response<=0); % leftward valid trials;
eyeTrialData.errorStatus(idxT) = -2;

% delete trials with little eye movements
idxT = find(eyeTrialData.errorStatus==0 & ...
    eyeTrialData.pursuit.gainXexternal<0.5); 
eyeTrialData.errorStatus(idxT) = -3;

idxT = find(eyeTrialData.errorStatus==0 & ...
    eyeTrialData.pursuit.travelClpDis<eyeTrialData.pursuit.targetClpDis/2); 
eyeTrialData.errorStatus(idxT) = -3;

load('summaryData')
load('summaryDataDiff')
load('summaryDataSub')
load('summarySacData')
perceptFolder = ['..\..\perceptPlots\'];
eyeTracesFolder = ['..\..\eyeTraces\'];
pursuitFolder = ['..\..\pursuitPlots\'];
saccadeFolder = ['..\..\saccadePlots\'];
correlationFolder = ['..\..\corrPlots\'];
perceptFolder = ['..\..\perceptualPlots\'];
MEFolder = ['..\..\motionEnergyPlots\'];
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
colorDotCons = [0 0 0; 248 118 109; 0 191 196]/255; % dot motion conditions
colorGroup = [199 124 255; 124 174 0; 0 0 0]/255; % perceptual bias group, assimilation-contrast-no bias
% colorProb = [8,48,107;198,219,239;8,48,107]/255; % all blue hues
% colorProb = [8,48,107;66,146,198;198,219,239;66,146,198;8,48,107]/255; % all blue hues
% colorPlot = [232 113 240; 15 204 255; 255 182 135; 137 126 255; 113 204 100]/255; % each row is one colour for one probability
colorObjAngles = [232 113 240; 15 204 255; 255 182 135; 0 0 0; 255 182 135; 15 204 255; 232 113 240]/255; % each row is one colour for one object motion angle