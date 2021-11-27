% Exp1 data, trying to build the model...
% ref: Bogadhi et al. 2011, Bogadhi et al. 2017, recurrent Bayesian networks for pursuit 

clear all; close all; clc

% let's start... with a simple version of simulation?
%% set experimental paramters
% actual parameters in the experiment
displayDuration = 800; % ms
samplingRate = 1000; % Hz
frameN = displayDuration/1000*samplingRate; % 
apertureAngle = 3; % in degs, horizontal right is zero, clockwise is positive
apertureSpeed = 10;
dotAngle = -90;
dotSpeed = 5;

% target velocity, columns are x and y
apertureV = ones(frameN, 2).*repmat([cos(apertureAngle/180*pi), sin(apertureAngle/180*pi)]*apertureSpeed, frameN, 1);
dotV = ones(frameN, 2).*repmat([cos(dotAngle/180*pi), sin(dotAngle/180*pi)]*dotSpeed, frameN, 1);

gain = 1; % pursuit gain...

%% initialize model parameters
% temporal delays, in ms
delayRD = 30; % retinal of dot motion
delayRA = delayRD+45; % retinal of aperture motion
delayEx = 97; % extra-retinal signals
delayEye = 50; % delay from motor command to eye movements

% constant values
sigmaD = 0.5; % eye velocity std to pure dot motion
sigmaA = 0.7; % eye velocity std to pure aperture motion

% these will be updated every frame
% dynamic weights (to be fitted in the actual model, now just simulate...)
wRDV = linspace(0.1, 0.05, frameN)'; % retinal dot velocity
wRAV = 1-wRDV; % retinal aperture velocity
wR = [linspace(1, 0.7, 50)'; linspace(0.1, 0, frameN-50)']; % retinal signal
wEx = 1-wR; % extra-retinal signal

% velocity signals
eyeVel = [0, 0]; % eye velocity
exRV = [0, 0];

%% enter the simulation loop
retinalDV = [];
retinalAV = [];
retinalV = [];
internalRep = [];
for f = 1:frameN
    % processing retinal signals
    retinalDV(f, :) = dotV(f, :)-eyeVel(f, :); % dot
    retinalAV(f, :) = apertureV(f, :)-eyeVel(f, :); % aperture

    % combining retinal signals
    if f<delayRD
        retinalV(f, :) = [0, 0]; % no signals available yet
    elseif f<delayRA
        retinalV(f, :) = retinalDV(f-delayRD+1, :);
    else
        retinalV(f, :) = wRDV(f-delayRD+1).*retinalDV(f-delayRD+1, :) + ...
            wRAV(f-delayRA+1).*retinalAV(f-delayRA+1, :);
    end
    
    % combining with extra-retinal signals
    if f<delayEx
        internalRep(f, :) = retinalV(f, :);
    else
        internalRep(f, :) = wR(f-delayEx+1).*retinalV(f, :) + ...
            wEx(f-delayEx+1).*exRV(f, :);
    end
    
    % updating eyeVel
    if f<delayEye
        eyeVel(f+1, :) = [0, 0];
    else
        eyeVel(f+1, :) = gain*internalRep(f-delayEye+1, :);
    end
    
    % generating the new extra-retinal signals 
    exRV(f+1, :) = eyeVel(f+1, :);
end

% check the plots...
time = 1:1:frameN;
figure
subplot(2, 1, 1)
hold on
plot(time, apertureV(:, 1), '-k')
plot(time, eyeVel(1:end-1, 1), '--b')
plot(time, retinalV(:, 1), '--r')
plot(time, internalRep(:, 1), '--g')
xlim([-10, 800])
ylim([0, 11])

subplot(2, 1, 2)
hold on
plot(time, apertureV(:, 2), '-k')
plot(time, eyeVel(1:end-1, 2), '--b')
plot(time, retinalV(:, 2), '--r')
plot(time, internalRep(:, 2), '--g')
xlim([-10, 800])
ylim([-6, 6])