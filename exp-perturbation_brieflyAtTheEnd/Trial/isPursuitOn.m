function [pursuitOn] = isPursuitOn(eyePos, threshold, filter)
% for the online detection of pursuit onset

% filter position data and get velocity data... same way as we do for the
% offline analysis

%% position
X_filt = filtfilt(filter.a, filter.b, eyePos(:, 1));
Y_filt = filtfilt(filter.a, filter.b, eyePos(:, 2));

%% velocity
DX_tmp = diff(X_filt)*sampleRate;
DX_filt = filtfilt(filter.c, filter.d, DX_tmp);
DY_tmp = diff(Y_filt)*sampleRate;
DY_filt = filtfilt(filter.c, filter.d, DY_tmp);

pursuitOn = all(sqrt(DX_filt.^2 + DY_filt.^2)>threshold);

end

