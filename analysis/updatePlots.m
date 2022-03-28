% FUNCTION to plot eye data in GUI created by viewEyeData.m
% Note that this function probably needs lots of changing depending on what
% you want to look at in your analysis e.g. there is a different version
% for EyeStrike/Eyecatch
% history
% 07-2012       JE created updatePlots.m
% 2012-2018     JF added stuff to and edited updatePlots.m
% 16-07-2018    JF commented to make the script more accecable for future
%               VPOM students
% for questions email jolande.fooken@rwth-aachen.de
%
% input: trial --> structure containing relevant current trial information
%                  (remember by now we have stored saccade information in
%                   trial as well)
% output: this will plot in the open figure

function [] = updatePlots(trial)
% define window for which you want to plot your data
startFrame = max(1, trial.log.targetOnset-200);
endFrame = min(trial.log.targetOffset+200, trial.log.trialEnd); % trial.log.trialEnd; % this is all recorded eye movement data
% if the interval looking at micro-saccades differs define it here
% msStart = trial.log.microSaccade.onset;
% msEnd = trial.log.microSaccade.offset;
stimOnset = trial.log.targetOnset; % this may have to be changed depending on terminology
stimOffset = trial.log.targetOffset;
% range of the x axis
minPosX = -10;
maxPosX = 10;
minPosY = -3;
maxPosY = 3;
minVel = -30;
maxVel = 30;
minAcc = -500;
maxAcc = 500;
minJerk = -300000;
maxJerk = 300000;
% in subplot we divide the screen into potential plotting windows; i.e.
% subplot(2,2) means that the figure is divided into a grid of 2x2 and we
% are plotting in the first position; 'replace' allows us to refresh the
% plot every time we "click through"

%% for the manuscript, single trial example
axis([-3 8 -3 8]);
axis([-3 8 -1 1]);
axis square
hold on
xlabel('x-position (deg)', 'fontsize', 12);
ylabel('y-position (deg)', 'fontsize', 12);
% plot eye x- versus y-position
% plot(trial.eyeX_filt(startFrame:stimOnset)-trial.eyeX_filt(stimOnset), trial.eyeY_filt(startFrame:stimOnset)-trial.eyeY_filt(stimOnset), 'k');
plot(trial.eyeX_filt(stimOnset:stimOffset)-trial.eyeX_filt(stimOnset), trial.eyeY_filt(stimOnset:stimOffset)-trial.eyeY_filt(stimOnset), 'b');
% mark pursuit onset and open-loop phase
plot(trial.eyeX_filt(trial.pursuit.summary.onset)-trial.eyeX_filt(stimOnset), trial.eyeY_filt(trial.pursuit.summary.onset)-trial.eyeY_filt(stimOnset), 'k+', 'MarkerSize', 15);
plot(trial.eyeX_filt(trial.pursuit.summary.onset+140)-trial.eyeX_filt(stimOnset), trial.eyeY_filt(trial.pursuit.summary.onset+140)-trial.eyeY_filt(stimOnset), 'k+', 'MarkerSize', 15);
plot(trial.eyeX_filt(stimOffset-100)-trial.eyeX_filt(stimOnset), trial.eyeY_filt(stimOffset-100)-trial.eyeY_filt(stimOnset), 'k+', 'MarkerSize', 15);

% draw the RDK aperture...
th = 0:pi/50:2*pi;
r = 1;
x1 = r * cos(th) + trial.target.posX(stimOnset)-trial.eyeX_filt(stimOnset);
y1 = r * sin(th) + trial.target.posY(stimOnset)-trial.eyeY_filt(stimOnset);
x2 = r * cos(th) + trial.target.posX(stimOffset)-trial.eyeX_filt(stimOnset);
y2 = r * sin(th) + trial.target.posY(stimOffset)-trial.eyeY_filt(stimOnset);
plot(x1, y1, '--g')
plot(x2, y2, '--g')

% plot(trial.eyeX_filt(stimOffset+1:endFrame), trial.eyeY_filt(stimOffset+1:endFrame), 'k');
% plot target center position
plot(trial.target.posX(stimOnset:stimOffset)-trial.eyeX_filt(stimOnset), trial.target.posY(stimOnset:stimOffset)-trial.eyeY_filt(stimOnset), 'g');
saveas(gcf, 'xw1t21_angle6.pdf')

%% the usual plots...
% % 2D position plot
% subplot(2,2,1,'replace');
% axis([minPosX maxPosX minPosY maxPosY]);
% hold on
% xlabel('x-position (deg)', 'fontsize', 12);
% ylabel('y-position (deg)', 'fontsize', 12);
% % plot eye x- versus y-position
% plot(trial.eyeX_filt(startFrame:stimOnset), trial.eyeY_filt(startFrame:stimOnset), 'k');
% plot(trial.eyeX_filt(stimOnset+1:stimOffset), trial.eyeY_filt(stimOnset+1:stimOffset), 'b');
% % mark pursuit onset and open-loop phase
% plot(trial.eyeX_filt(trial.pursuit.summary.onset), trial.eyeY_filt(trial.pursuit.summary.onset), 'bo');
% plot(trial.eyeX_filt(trial.pursuit.summary.onset+140), trial.eyeY_filt(trial.pursuit.summary.onset+140), 'go');
% 
% plot(trial.eyeX_filt(stimOffset+1:endFrame), trial.eyeY_filt(stimOffset+1:endFrame), 'k');
% % plot target center position
% plot(trial.target.posX(startFrame:endFrame), trial.target.posY(startFrame:endFrame), 'g');
% 
% % eye position plot over time
% subplot(2,2,2,'replace');
% % define some plot parameters
% axis([startFrame endFrame minPosX maxPosX]);
% hold on;
% xlabel('Time(ms)', 'fontsize', 12);
% ylabel('Position (degree)', 'fontsize', 12);
% % plot x- and y- eye position over time
% plot(startFrame:endFrame,trial.eyeX_filt(startFrame:endFrame),'k');
% plot(startFrame:endFrame,trial.eyeY_filt(startFrame:endFrame),'b');
% plot(trial.saccades.X_left.onsets,trial.eyeX_filt(trial.saccades.X_left.onsets),'go');
% plot(trial.saccades.X_left.offsets,trial.eyeX_filt(trial.saccades.X_left.offsets),'mo');
% plot(trial.saccades.X_right.onsets,trial.eyeX_filt(trial.saccades.X_right.onsets),'g*');
% plot(trial.saccades.X_right.offsets,trial.eyeX_filt(trial.saccades.X_right.offsets),'m*');
% plot(trial.saccades.Y.onsets,trial.eyeY_filt(trial.saccades.Y.onsets),'y*');
% plot(trial.saccades.Y.offsets,trial.eyeY_filt(trial.saccades.Y.offsets),'c*');
% % legend({'x pos','y pos', 'sacLeftOn', 'sacLeftOff', 'sacRightOn', 'sacRightOff', ...
% %     'sacOn', 'sacOff'},'Location','NorthWest');%, 'AutoUpdate','off');
% % vertical lines indicate events/target onsets
% line([trial.log.targetOnset trial.log.targetOnset], [minPosX maxPosX],'Color','k','LineStyle','--');
% line([trial.stim_offset trial.stim_offset], [minPosX maxPosX],'Color','k','LineStyle','--');
% % if trial.log.eyeCondition==1 % pursuit trials
%     line([trial.pursuit.summary.onset trial.pursuit.summary.onset], [minPosX maxPosX],'Color','b','LineStyle','--');
%     line([trial.pursuit.summary.onset+140 trial.pursuit.summary.onset+140], [minPosX maxPosX],'Color','b','LineStyle','--');
% % end
% line([startFrame endFrame], [0 0],'Color','k','LineStyle','--')
% % if ~isempty(trial.pursuit.onset)
% %     line([trial.pursuit.onset trial.pursuit.onset], [minPosAbs maxPosAbs],'Color','r','LineStyle','--');
% %     line([trial.pursuit.onsetTrue trial.pursuit.onsetTrue], [minPosAbs maxPosAbs],'Color','r','LineStyle','-.');
% % end
% % if ~isempty(trial.pursuit.onsetSteadyState)
% %     line([trial.pursuit.onsetSteadyState trial.pursuit.onsetSteadyState], [minVel maxVel],'Color','r','LineStyle','--');
% % end
% 
% % velocity plot over time
% subplot(2,2,3,'replace');
% axis([startFrame endFrame minVel maxVel]);
% hold on;
% xlabel('Time(ms)', 'fontsize', 12);
% ylabel('Speed (degree/second)', 'fontsize', 12);
% % plot x- and y- eye velocity over time
% plot(startFrame:endFrame,trial.eyeDX_filt(startFrame:endFrame),'k');
% plot(startFrame:endFrame,trial.eyeDY_filt(startFrame:endFrame),'b');
% % also plot 2D velocity...
% vel2D = sqrt(trial.DX_interpolSac.^2 + trial.DY_interpolSac.^2);
% plot(startFrame:endFrame,vel2D(startFrame:endFrame),'r');
% plot(trial.saccades.onsets,vel2D(trial.saccades.onsets),'y*');
% plot(trial.saccades.offsets,vel2D(trial.saccades.offsets),'c*'); % plot combined saccades found
% % plot saccade onsets in x- and y with different colors
% plot(trial.saccades.X_left.onsets,trial.eyeDX_filt(trial.saccades.X_left.onsets),'go');
% plot(trial.saccades.X_left.offsets,trial.eyeDX_filt(trial.saccades.X_left.offsets),'mo');
% plot(trial.saccades.X_right.onsets,trial.eyeDX_filt(trial.saccades.X_right.onsets),'g*');
% plot(trial.saccades.X_right.offsets,trial.eyeDX_filt(trial.saccades.X_right.offsets),'m*');
% plot(trial.saccades.Y.onsets,trial.eyeDY_filt(trial.saccades.Y.onsets),'y*');
% plot(trial.saccades.Y.offsets,trial.eyeDY_filt(trial.saccades.Y.offsets),'c*');
% % vertical lines indicate events/target onsets
% line([trial.log.targetOnset trial.log.targetOnset], [minVel maxVel],'Color','k','LineStyle','--');
% % line([trial.stim_offset-ms2frames(100) trial.stim_offset-ms2frames(100)], [minVel maxVel],'Color','b','LineStyle','--');
% line([trial.stim_offset trial.stim_offset], [minVel maxVel],'Color','k','LineStyle','--');
% % if trial.log.eyeCondition==1 % pursuit trials
%     line([trial.pursuit.summary.onset trial.pursuit.summary.onset], [minVel maxVel],'Color','b','LineStyle','-');
%     line([trial.pursuit.summary.onset+140 trial.pursuit.summary.onset+140], [minVel maxVel],'Color','b','LineStyle','--');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%% debugging for pursuit latency
% %     line([trial.pursuit.summary.baselineStart trial.pursuit.summary.baselineStart], [minVel maxVel],'Color','g','LineStyle','--');
% %     line([trial.pursuit.summary.baselineEnd trial.pursuit.summary.baselineEnd], [minVel maxVel],'Color','g','LineStyle','--');
% %     line([trial.pursuit.summary.pursuitIntervalStart trial.pursuit.summary.pursuitIntervalStart], [minVel maxVel],'Color','r','LineStyle','--');
% %     line([trial.pursuit.summary.pursuitIntervalEnd trial.pursuit.summary.pursuitIntervalEnd], [minVel maxVel],'Color','r','LineStyle','--');
% %     line([startFrame endFrame], [trial.pursuit.summary.meanBaseline+4*trial.pursuit.summary.sdBaseline trial.pursuit.summary.meanBaseline+4*trial.pursuit.summary.sdBaseline],'Color','b','LineStyle','--')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % end
% line([startFrame endFrame], [0 0],'Color','k','LineStyle','--')
% 
% % if ~isempty(trial.pursuit.onset)
% %     line([trial.pursuit.onset trial.pursuit.onset], [minVel maxVel],'Color','r','LineStyle','--');
% %     line([trial.pursuit.onsetTrue trial.pursuit.onsetTrue], [minVel maxVel],'Color','r','LineStyle','-.');
% %     line([trial.pursuit.openLoopEndFrame trial.pursuit.openLoopEndFrame], [minVel maxVel],'Color','m','LineStyle','--');
% % end
% % if ~isempty(trial.pursuit.onsetSteadyState)
% %     line([trial.pursuit.onsetSteadyState trial.pursuit.onsetSteadyState], [minVel maxVel],'Color','r','LineStyle','--');
% % end
% 
% % acceleration plot over time
% subplot(2,2,4,'replace');
% axis([startFrame endFrame minAcc maxAcc]);
% hold on;
% xlabel('Time(ms)', 'fontsize', 12);
% ylabel('Acceleration', 'fontsize', 12);
% % plot x- and y- eye acceleration over time
% plot(startFrame:endFrame,trial.eyeDDX_filt(startFrame:endFrame),'k');
% plot(startFrame:endFrame,trial.eyeDDY_filt(startFrame:endFrame),'b');
% % plot saccade onsets in x- and y with different colors
% plot(trial.saccades.X_left.onsets,trial.eyeDDX_filt(trial.saccades.X_left.onsets),'go');
% plot(trial.saccades.X_left.offsets,trial.eyeDDX_filt(trial.saccades.X_left.offsets),'mo');
% plot(trial.saccades.X_right.onsets,trial.eyeDDX_filt(trial.saccades.X_right.onsets),'g*');
% plot(trial.saccades.X_right.offsets,trial.eyeDDX_filt(trial.saccades.X_right.offsets),'m*');
% 
% % % absolute values
% % plot(startFrame:endFrame,abs(trial.eyeDDX_filt(startFrame:endFrame)),'k');
% % plot(trial.saccades.X_left.onsets,abs(trial.eyeDDX_filt(trial.saccades.X_left.onsets)),'go');
% % plot(trial.saccades.X_left.offsets,abs(trial.eyeDDX_filt(trial.saccades.X_left.offsets)),'mo');
% % plot(trial.saccades.X_right.onsets,abs(trial.eyeDDX_filt(trial.saccades.X_right.onsets)),'g*');
% % plot(trial.saccades.X_right.offsets,abs(trial.eyeDDX_filt(trial.saccades.X_right.offsets)),'m*');
% 
% plot(trial.saccades.Y.onsets,trial.eyeDDY_filt(trial.saccades.Y.onsets),'y*');
% plot(trial.saccades.Y.offsets,trial.eyeDDY_filt(trial.saccades.Y.offsets),'c*');
% % vertical lines indicate events/target onsets
% line([trial.log.targetOnset trial.log.targetOnset], [minAcc maxAcc],'Color','k','LineStyle','--');
% line([trial.stim_offset trial.stim_offset], [minAcc maxAcc],'Color','k','LineStyle','--');
% % if trial.log.eyeCondition==1 % pursuit trials
%     line([trial.pursuit.summary.onset trial.pursuit.summary.onset], [minAcc maxAcc],'Color','b','LineStyle','--');
%     line([trial.pursuit.summary.onset+140 trial.pursuit.summary.onset+140], [minAcc maxAcc],'Color','b','LineStyle','--');
% % end
% line([startFrame endFrame], [0 0],'Color','k','LineStyle','--')
% % if ~isempty(trial.pursuit.onset)
% %     line([trial.pursuit.onset trial.pursuit.onset], [minAcc maxAcc],'Color','r','LineStyle','--');
% %     line([trial.pursuit.onsetTrue trial.pursuit.onsetTrue], [minAcc maxAcc],'Color','r','LineStyle','-.');
% %     line([trial.pursuit.openLoopEndFrame trial.pursuit.openLoopEndFrame], [minAcc maxAcc],'Color','m','LineStyle','--');
% % end
end
