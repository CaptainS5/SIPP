% FUNCTION to remove saccades from eye movement data; this may be necessary
% when e.g. analyzing smooth eye movement phase

% history
% 07-2012       JE created analyzeSaccades.m
% 2012-2018     JF added stuff to and edited analyzeSaccades.m
% 13-07-2018    JF commented to make the script more accecable for future
%               VPOM students
% for questions email jolande.fooken@rwth-aachen.de
%
% input: trial --> structure containing relevant current trial information
% output: trial --> structure containing relevant current trial information
%                   de-saccaded eye movements added

function [trial] = removeSaccades(trial, saccades)
% add saccades to trial information
% for x
trial.saccades.X.onsets = [];
trial.saccades.X.offsets = [];
duringIdx = 1;

trial.saccades.X_left.onsets = [];
trial.saccades.X_left.offsets = [];
trial.saccades.X_right.onsets = [];
trial.saccades.X_right.offsets = [];

for ii = 1:length(saccades.X.onsets)
    trial.saccades.X.onsets(ii,1) = saccades.X.onsets(ii); % ?... why use the loop
    trial.saccades.X.offsets(ii,1) = saccades.X.offsets(ii);
    if trial.saccades.X.onsets(ii,1)>=trial.stim_onset && trial.saccades.X.offsets(ii,1)<=trial.stim_offset
        trial.saccades.X.onsetsDuring(duringIdx, 1) = trial.saccades.X.onsets(ii,1);
        trial.saccades.X.offsetsDuring(duringIdx, 1) = trial.saccades.X.offsets(ii,1);
        duringIdx = duringIdx + 1;
    end
    % figure out saccade direction
    peakV = max(abs(trial.eyeDX_filt(saccades.X.onsets(ii):saccades.X.offsets(ii))));
    peakVIdx = find(abs(trial.eyeDX_filt)==peakV);
    if length(peakVIdx)>1
        for ii = 1:length(peakVIdx)
            if peakVIdx(ii) > saccades.X.onsets(ii) && peakVIdx(ii) < saccades.X.offsets(ii)
                peakVIdx = peakVIdx(ii);
                break
            end
        end
    end
    if trial.eyeDX_filt(peakVIdx) < 0 % whether velocity is positive or negative
        trial.saccades.X_left.onsets = [trial.saccades.X_left.onsets; saccades.X.onsets(ii)];
        trial.saccades.X_left.offsets = [trial.saccades.X_left.offsets; saccades.X.offsets(ii)];
    else
        trial.saccades.X_right.onsets = [trial.saccades.X_right.onsets; saccades.X.onsets(ii)];
        trial.saccades.X_right.offsets = [trial.saccades.X_right.offsets; saccades.X.offsets(ii)];
    end
end
if duringIdx>1 
    trial.saccades.firstSaccadeOnset = trial.saccades.X.onsetsDuring(1, 1);
else
    trial.saccades.firstSaccadeOnset = [];
    trial.saccades.X.onsetsDuring = [];
    trial.saccades.X.offsetsDuring = [];
end

% and for y
trial.saccades.Y.onsets = [];
trial.saccades.Y.offsets = [];
duringIdxY = 1;
% y-left is down, y-right is up, just the negative & positive--need to
% confirm with data!
trial.saccades.Y_left.onsets = [];
trial.saccades.Y_left.offsets = [];
trial.saccades.Y_right.onsets = [];
trial.saccades.Y_right.offsets = [];
for ii = 1:length(saccades.Y.onsets)    
    trial.saccades.Y.onsets(ii,1) = saccades.Y.onsets(ii);
    trial.saccades.Y.offsets(ii,1) = saccades.Y.offsets(ii);
    if trial.saccades.Y.onsets(ii,1)>=trial.stim_onset && trial.saccades.Y.offsets(ii,1)<=trial.stim_offset
        trial.saccades.Y.onsetsDuring(duringIdxY, 1) = trial.saccades.Y.onsets(ii,1);
        trial.saccades.Y.offsetsDuring(duringIdxY, 1) = trial.saccades.Y.offsets(ii,1);
        duringIdxY = duringIdxY + 1;
    end
    % figure out saccade direction
    peakV = max(abs(trial.eyeDY_filt(saccades.Y.onsets(ii):saccades.Y.offsets(ii))));
    peakVIdx = find(abs(trial.eyeDY_filt)==peakV);
    if length(peakVIdx)>1
        for ii = 1:length(peakVIdx)
            if peakVIdx(ii) > saccades.Y.onsets(ii) && peakVIdx(ii) < saccades.Y.offsets(ii)
                peakVIdx = peakVIdx(ii);
                break
            end
        end
    end
    if trial.eyeDY_filt(peakVIdx) < 0 % whether velocity is positive or negative
        trial.saccades.Y_left.onsets = [trial.saccades.Y_left.onsets; saccades.Y.onsets(ii)];
        trial.saccades.Y_left.offsets = [trial.saccades.Y_left.offsets; saccades.Y.offsets(ii)];
    else
        trial.saccades.Y_right.onsets = [trial.saccades.Y_right.onsets; saccades.Y.onsets(ii)];
        trial.saccades.Y_right.offsets = [trial.saccades.Y_right.offsets; saccades.Y.offsets(ii)];
    end
end
if duringIdxY>1
    trial.saccades.firstSaccadeOnset = min(trial.saccades.firstSaccadeOnset, trial.saccades.Y.onsetsDuring(1, 1));
else
    trial.saccades.Y.onsetsDuring = [];
    trial.saccades.Y.offsetsDuring = [];
end

% store all found on and offsets together
trial.saccades.onsets = [trial.saccades.X.onsets; trial.saccades.Y.onsets];
trial.saccades.offsets = [trial.saccades.X.offsets; trial.saccades.Y.offsets];
% trial.saccades.onsetsDuring = [trial.saccades.X.onsetsDuring; trial.saccades.Y.onsetsDuring];
% trial.saccades.offsetsDuring = [trial.saccades.X.offsetsDuring; trial.saccades.Y.offsetsDuring];
% merge saccades on X and Y that are actually the same...
xSac = length(trial.saccades.X.onsets);
ySac = length(trial.saccades.Y.onsets);
if ~isempty(ySac) && ~isempty(xSac) && numel(trial.saccades.onsets) ~= 0
    testOnsets = sort(trial.saccades.onsets);
    testOffsets = sort(trial.saccades.offsets);
    onsets = [];
    offsets = [];
    diffOnsets = [diff(testOnsets); 999];
    diffOffsets = [diff(testOffsets); 999];
    ii = 1;
    while ii <= length(testOnsets)
        if abs(diffOnsets(ii))<=25 || abs(diffOffsets(ii))<=25
            tempOnset = min(testOnsets(ii), testOnsets(ii+1));
            tempOffset = max(testOffsets(ii), testOffsets(ii+1));
            ii = ii+1;
        else
            tempOnset = testOnsets(ii);
            tempOffset = testOffsets(ii);
        end
        onsets = [onsets; tempOnset];
        offsets = [offsets; tempOffset];
        ii = ii+1;
    end
    trial.saccades.onsets = onsets;
    trial.saccades.offsets = offsets;
%     % original codes
%     testOnsets = sort(trial.saccades.onsets);
%     testOffsets = sort(trial.saccades.offsets);
%     count1 = 1;
%     tempOnset1 = [];
%     tempOffset1 = [];
%     count2 = 1;
%     tempOnset2 = [];
%     tempOffset2 = [];   
%     for ii = 1:length(testOnsets)-1
%         if testOnsets(ii+1)-testOnsets(ii) < 20
%             tempOnset1(count1) = testOnsets(ii);
%             tempOffset1(count1) = testOffsets(ii);
%             count1 = length(tempOnset1) +1;
%         else
%             tempOnset2(count2) = testOnsets(ii+1);
%             tempOffset2(count2) = testOffsets(ii+1);
%             count2 = length(tempOnset2) +1;
%         end
%     end
%     onsets = unique([tempOnset1 tempOnset2 testOnsets(1)])';
%     offsets = unique([tempOffset1 tempOffset2 testOffsets(1)])';
%     trial.saccades.onsets = onsets;
%     trial.saccades.offsets = offsets;
end

% xSac = length(trial.saccades.X.onsetsDuring);
% ySac = length(trial.saccades.Y.onsetsDuring);
% if ~isempty(ySac) && ~isempty(xSac) && numel(trial.saccades.onsetsDuring) ~= 0
%     testOnsets = sort(trial.saccades.onsetsDuring);
%     testOffsets = sort(trial.saccades.offsetsDuring);
%     count1 = 1;
%     tempOnset1 = [];
%     tempOffset1 = [];
%     count2 = 1;
%     tempOnset2 = [];
%     tempOffset2 = [];   
%     for i = 1:length(testOnsets)-1
%         if testOnsets(i+1)-testOnsets(i) < 20
%             tempOnset1(count1) = testOnsets(i);
%             tempOffset1(count1) = testOffsets(i);
%             count1 = length(tempOnset1) +1;
%         else
%             tempOnset2(count2) = testOnsets(i+1);
%             tempOffset2(count2) = testOffsets(i+1);
%             count2 = length(tempOnset2) +1;
%         end
%     end
%     onsets = unique([tempOnset1 tempOnset2 testOnsets(1)])';
%     offsets = unique([tempOffset1 tempOffset2 testOffsets(1)])';
%     trial.saccades.onsetsDuring = onsets;
%     trial.saccades.offsetsDuring = offsets;
% end

% open eye movement data structure
trial.X_noSac = trial.eyeX_filt;
trial.Y_noSac = trial.eyeY_filt;
trial.DX_noSac = trial.eyeDX_filt;
trial.DY_noSac = trial.eyeDY_filt;
trial.DDX_noSac = trial.eyeDDX_filt;
trial.DDY_noSac = trial.eyeDDY_filt;
trial.X_interpolSac = trial.eyeX_filt;
trial.Y_interpolSac = trial.eyeY_filt;
trial.DX_interpolSac = trial.eyeDX_filt;
trial.DY_interpolSac = trial.eyeDY_filt;
trial.DDX_interpolSac = trial.eyeDDX_filt;
trial.DDY_interpolSac = trial.eyeDDY_filt;
trial.quickphases = false(trial.length,1);
% now remove saccadic phase
for ii = 1:length(trial.saccades.onsets)
    % first we calculate the slope between the eye position at saccade on-
    % to saccade offset
    lengthSacX = trial.saccades.offsets(ii) - trial.saccades.onsets(ii);
    slopeX = (trial.eyeX_filt(trial.saccades.offsets(ii))-trial.eyeX_filt(trial.saccades.onsets(ii)))./lengthSacX;
    slopeDX = (trial.eyeDX_filt(trial.saccades.offsets(ii))-trial.eyeDX_filt(trial.saccades.onsets(ii)))./lengthSacX;
    slopeDDX = (trial.eyeDDX_filt(trial.saccades.offsets(ii))-trial.eyeDDX_filt(trial.saccades.onsets(ii)))./lengthSacX;
    slopeY = (trial.eyeY_filt(trial.saccades.offsets(ii))-trial.eyeY_filt(trial.saccades.onsets(ii)))./lengthSacX;
    slopeDY = (trial.eyeDY_filt(trial.saccades.offsets(ii))-trial.eyeDY_filt(trial.saccades.onsets(ii)))./lengthSacX;
    slopeDDY = (trial.eyeDDY_filt(trial.saccades.offsets(ii))-trial.eyeDDY_filt(trial.saccades.onsets(ii)))./lengthSacX;
    % now we can add a completely de-saccaded variable in trial
    trial.X_noSac(trial.saccades.onsets(ii):trial.saccades.offsets(ii)) = NaN;
    trial.Y_noSac(trial.saccades.onsets(ii):trial.saccades.offsets(ii)) = NaN;
    trial.DX_noSac(trial.saccades.onsets(ii):trial.saccades.offsets(ii)) = NaN;
    trial.DY_noSac(trial.saccades.onsets(ii):trial.saccades.offsets(ii)) = NaN;
    trial.DDX_noSac(trial.saccades.onsets(ii):trial.saccades.offsets(ii)) = NaN;
    trial.DDY_noSac(trial.saccades.onsets(ii):trial.saccades.offsets(ii)) = NaN;
    % and finally interpolate the eye position if we later want to plot
    % smooth eye movement traces
    for j = 1:lengthSacX+1
        trial.X_interpolSac(trial.saccades.onsets(ii)-1+j) = trial.eyeX_filt(trial.saccades.onsets(ii)) + slopeX*j;
        trial.Y_interpolSac(trial.saccades.onsets(ii)-1+j) = trial.eyeY_filt(trial.saccades.onsets(ii)) + slopeY*j;
        trial.DX_interpolSac(trial.saccades.onsets(ii)-1+j) = trial.eyeDX_filt(trial.saccades.onsets(ii)) + slopeDX*j;
        trial.DY_interpolSac(trial.saccades.onsets(ii)-1+j) = trial.eyeDY_filt(trial.saccades.onsets(ii)) + slopeDY*j;
        trial.DDX_interpolSac(trial.saccades.onsets(ii)-1+j) = trial.eyeDDX_filt(trial.saccades.onsets(ii)) + slopeDDX*j;
        trial.DDY_interpolSac(trial.saccades.onsets(ii)-1+j) = trial.eyeDDY_filt(trial.saccades.onsets(ii)) + slopeDDY*j;
    end   
end
% % do the exact same thing for y
% for i = 1:length(trial.saccades.Y.onsets)
%     
%     lengthSacY = trial.saccades.Y.offsets(i) - trial.saccades.Y.onsets(i);
%     slopeY = (trial.eyeY_filt(trial.saccades.Y.offsets(i))-trial.eyeY_filt(trial.saccades.Y.onsets(i)))./lengthSacY;
%     slopeDY = (trial.eyeDY_filt(trial.saccades.Y.offsets(i))-trial.eyeDY_filt(trial.saccades.Y.onsets(i)))./lengthSacY;
%     slopeX = (trial.eyeX_filt(trial.saccades.Y.offsets(i))-trial.eyeX_filt(trial.saccades.Y.onsets(i)))./lengthSacY;
%     slopeDX = (trial.eyeDX_filt(trial.saccades.Y.offsets(i))-trial.eyeDX_filt(trial.saccades.Y.onsets(i)))./lengthSacY;
%     
%     trial.X_noSac(trial.saccades.Y.onsets(i):trial.saccades.Y.offsets(i)) = NaN;
%     trial.Y_noSac(trial.saccades.Y.onsets(i):trial.saccades.Y.offsets(i)) = NaN;
%     trial.DX_noSac(trial.saccades.Y.onsets(i):trial.saccades.Y.offsets(i)) = NaN;
%     trial.DY_noSac(trial.saccades.Y.onsets(i):trial.saccades.Y.offsets(i)) = NaN;
%     
%     for j = 1:lengthSacY+1
%         trial.Y_interpolSac(trial.saccades.Y.onsets(i)-1+j) = trial.eyeY_filt(trial.saccades.Y.onsets(i)) + slopeY*j;
%         trial.X_interpolSac(trial.saccades.Y.onsets(i)-1+j) = trial.eyeX_filt(trial.saccades.Y.onsets(i)) + slopeX*j;
%         trial.DY_interpolSac(trial.saccades.Y.onsets(i)-1+j) = trial.eyeDY_filt(trial.saccades.Y.onsets(i)) + slopeDY*j;
%         trial.DX_interpolSac(trial.saccades.Y.onsets(i)-1+j) = trial.eyeDX_filt(trial.saccades.Y.onsets(i)) + slopeDX*j;
%     end    
% end
% done
end
