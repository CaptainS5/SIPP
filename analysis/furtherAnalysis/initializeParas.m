% use eyeTrialData to do analysis, initialize parameters for the predictive
% pursuit task

clear all; close all; clc

names = {'w00' 'w01' 'w02' 'w03' 'w04' 'w05' 'w06' 'w07' 'w08' 'w09' 'w10'};

internalDirCons = [-1 1]; % -1:left, 1-right
cohCons = [0; 0.5; 1]; % RDK coherence
internalDirNames = {'up', 'down'};
cohNames = {'coh-0', 'coh-50%', 'coh-100%'};

sampleRate = 1000;

analysisFolder = pwd;
load('eyeTrialData_all.mat');
perceptFolder = ['..\perceptPlots\'];
eyeTracesFolder = ['..\eyeTraces\'];
pursuitFolder = ['..\pursuitPlots\'];
saccadeFolder = ['..\saccadePlots\'];
correlationFolder = ['..\corrPlots\'];
% RFolder = pwd;

% for plotting
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
% colorProb = [8,48,107;198,219,239;8,48,107]/255; % all blue hues
% colorProb = [8,48,107;66,146,198;198,219,239;66,146,198;8,48,107]/255; % all blue hues
colorPlot = [232 113 240; 15 204 255; 255 182 135; 137 126 255; 113 204 100]/255; % each row is one colour for one probability