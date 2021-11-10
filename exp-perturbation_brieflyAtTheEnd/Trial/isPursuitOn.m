function [pursuitOn] = isPursuitOn(screen, eyePos, window, threshold, filter)
% for the online detection of pursuit onset

% filter position data and get velocity data... same way as we do for the
% offline analysis

% % for debug, copy the following and run in command window...
eyePos = control.eyePos(:, 1:2);
window =control.detectWindow;
threshold = control.velThres;
filter = control.filter;

%% position
X_filtP = filtfilt(filter.a, filter.b, eyePos(:, 1)); % in pixels
Y_filtP = filtfilt(filter.a, filter.b, eyePos(:, 2)); % in pixels

X_filt = X_filtP*screen.dpp;
Y_filt = Y_filtP*screen.dpp; % in degs

%% velocity
DX_tmp = diff(X_filt)*filter.sampleRate;
DX_filt = filtfilt(filter.c, filter.d, DX_tmp);
DY_tmp = diff(Y_filt)*filter.sampleRate;
DY_filt = filtfilt(filter.c, filter.d, DY_tmp);

vel2D = sqrt(DX_filt.^2 + DY_filt.^2);
pursuitOn = all(vel2D(end-window+1:end)>threshold);

end

