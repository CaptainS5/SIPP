function [const, control] = drawSingleFrameEyelink(const, trialData, control, screen, photo, rdkControl, eyelink)
% =========================================================================
% runSingleTrial(const, trialData, control, screen, photo)
% =========================================================================
% Single Trial Routine to draw all stimuli on the screen, go through
% different phases of the trial. Note, this function is called in the
% while-loop (i.e it's being run every refresh rate cycle!) - anything
% computationally expensive should be programmed outside the while loop!
% 1. Initial Drawings (BG, Stationary stim)
% 3. Different Phases of Experiment
% 3.1. Ready - check initial fixation
% 3.2. Set - check peripheral fixation
% 3.3. Intercept - target motion etc
% 3.4. finish - feedback etc.
%
% Last Changes:
% 23/10/2019 (PK) - cleaned up function
% 09/29/2020 XW - adapt for the motionIntegration SAT exp
% -------------------------------------------------------------------------
% Inputs:
% const:     structure containing different constant settings
% trialData: structure containing all trial relevant information
% control:   structure containing current trial information
% screen:    strucrure containing screen settings
% photo:     structure containting photodiode settings
% -------------------------------------------------------------------------
%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) INITIAL DRAWING:
%     - Set up the trial timer t (counting up from 0 for every trial!)
%     - Draw Background & stationary things that are always visible during
%     the trial (pre-defined in generateStationaryStimuli)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial timer:
timePassed = control.data_time;                                                      % in expMain this is defined as GetSecs-tMainSync (i.e. equivalent to when experiment started
% if const.makeVideo == 1; timePassed = control.frameCounter*screen.refreshRate; end

curTrial = control.currentTrial;                                         % current trial
Screen('FillRect', screen.window, screen.background);                            % Draw background

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) ENTER DIFFERENT STAGES OF THE TRIAL
%     (1) Fixation:
%     (2) RDK display:
%     (3) If not responded yet, response screen: show "too slow" in fast
%     blocks
%     (4) Finish:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch control.mode
%% (2.1) --------------------------------------------------------------
    case 1                                                                  % PHASE 1: Fixation
        % calculate the step-ramp distance, then place fixation relatively
        % closer to the screen center
        stepTime = 0.15; % how many secs it takes for the aperture center to move from the step to the fixation
        [stepDis, ] = dva2pxl(control.rdkApertureSpeed*stepTime, control.rdkApertureSpeed*stepTime, screen);
        
        if control.rdkApertureDir==0 % moving rightward
            fixationCenter = [rdkControl.apertureCenterPos{1}(1)+stepDis, rdkControl.apertureCenterPos{1}(2)];
        else % moving leftward
            fixationCenter = [rdkControl.apertureCenterPos{1}(1)-stepDis, rdkControl.apertureCenterPos{1}(2)];
        end
        
        % draw target
        if eyelink.mode && ~eyelink.dummy
            if Eyelink( 'NewFloatSampleAvailable') > 0
                % get the sample in the form of an event structure
                evt = Eyelink( 'NewestFloatSample');
                xeye = evt.gx(eyelink.eye_used+1); % +1 as we're accessing MATLAB array
                yeye = evt.gy(eyelink.eye_used+1);

                %% do we have valid data and is the pupil visible?
                if xeye~=eyelink.el.MISSING_DATA && yeye~=eyelink.el.MISSING_DATA && evt.pa(eyelink.eye_used+1)>0
                    % if data is valid, compare gaze position with the limits of the tolerance window
                    diffFix = sqrt((xeye-fixationCenter(1))^2+((yeye-fixationCenter(2))/screen.pixelRatioWidthPerHeight)^2);
                    
                    if diffFix <= const.fixation.windowRadiusPxl % fixation ok
                        PTBdraw_circles(screen, fixationCenter, const.fixation.dotRadiusPxl, const.fixation.colour);
                        control.frameFix = control.frameFix+1;
                    elseif diffFix > const.fixation.windowRadiusPxl % fixation out of range, show warning
                        %                         Snd('Play', const.beep.sound, const.beep.samplingRate, 16);
                        % Plays the sound in case of wrong fixation
                        % show white fixation
                        PTBdraw_circles(screen, fixationCenter, const.fixation.dotRadiusPxl, screen.warningColor);
                    end
                else
                    % if data is invalid (e.g. during a blink), show white
                    % fixation
                    PTBdraw_circles(screen, fixationCenter, const.fixation.dotRadiusPxl, screen.warningColor);
                    disp('Eyelink data invalid')
                end
            else
                error=Eyelink('CheckRecording');
%                 if(error~=0)
                    disp(['EyeLink CheckRecording Error: ', num2str(error)])
%                 end
                control.break = 1;
                control.forceRecalibEL = 1;
                control.repeat = 1;
            end
        else
            PTBdraw_circles(screen, fixationCenter, const.fixation.dotRadiusPxl, const.fixation.colour);
            control.frameFix = control.frameFix+1;
        end
        
        % draw square for photodiode:
         if photo.mode                                                       % Photodiode Event 1: Fixation target 1 on
            PTBdraw_photodiodeStimulus(screen, const.photoStimSizePX2, screen.white);
         end
        
        % Check whether Fixation time is passed:
        if control.frameFix == control.fixationFrames
            b_fix         = 1;
        else
            b_fix         = 0;
        end
        
        % if Fixation time passed, move to next phase
        if b_fix
            control.mode   = 2;                                             % leave ready-phase = fixation was successful, move on to next phase
        end
        
%% (2.2) --------------------------------------------------------------
    case 2                                                                  % PHASE 2: RDK display
        control.frameRDK = control.frameRDK + 1;   % start counting frames for RDK display
        
        % draw target:
        if control.frameRDK <= control.rdkFramesBefore || control.frameRDK > control.rdkFramesBefore + control.rdkFramesDuring % only draw target when it is not during occlusion
            PTBdraw_target_RDK(screen, const, rdkControl.dotPos{control.frameRDK}, rdkControl.apertureTexture, ...
                rdkControl.textureCenterPos{control.frameRDK}, rdkControl.textureWindow{control.frameRDK});
        end
        
        if control.occlusionStart~=0
            control.occlusionStart = 0; % reset
        end
        
        % draw square for photodiode:
        if photo.mode                                                       % Photodiode Event 2: Fixation target 2 on
            PTBdraw_photodiodeStimulus(screen, const.photoStimSizePX2, screen.black);
        end
        
%         % for debugging... can comment out
%         if control.frameRDK <= size(rdkControl.dir, 1)
%             meanRDKdir = mean(rdkControl.dir(control.frameRDK, :));
%             lineLength = 300; % in pixels
%             toH = screen.center(1)+cos(meanRDKdir/180*pi)*lineLength;
%             toV = screen.center(2)+sin(meanRDKdir/180*pi)*lineLength;
%             Screen('DrawLine', screen.window, screen.black, screen.center(1), screen.center(2), toH, toV, 5);
%         end
        
        if control.frameRDK == length(rdkControl.dotPos) % if RDK duration passed but participant did not respond, move to next phase
            b_rdk            = 1;
            control.mode    = 3; % skip response for now
            control.frameRDK = -1;
        elseif control.frameRDK == control.rdkFramesBefore + 1
            control.occlusionStart = 1;
        elseif control.frameRDK == control.rdkFramesBefore + control.rdkFramesDuring + 1
            control.occlusionStart = -1;
        else
            b_rdk            = 0;
        end
        
%% (2.3)---------------------------------------------------------------
    case 3                                                                  % PHASE 3: response screen if needed
        control.frameRDK = control.frameRDK - 1; % need to keep updating the frame numbers...
        
        % keyboard response message
%         PTBwrite_msg(screen, 'ahead of behind?', 'center', 'center', screen.msgfontcolour) % coordinate in relation to screen center
        PTBwrite_msg(screen, 'above or below?', 'center', 'center', screen.msgfontcolour) % coordinate in relation to screen center
        
%         % mouse response
%         [ecc, ] = dva2pxl(const.line.length/2, const.line.length/2, screen); % distance of the cursor from center
%         if isempty(control.mouse_x) % the first response frame, show random angle
%             % show the cursor; put it at the start angle everytime
%             %%% this is for only drawing the line, or only having rightward
%             %%% directions
%             control.respAngle = rand*180-90;
%             %%%
% %             %%% this is for drawing the arrow, range from -180 to 180
% %             control.respAngle = rand*360-180;
% %             %%%
%             SetMouse(screen.x_mid + round(cos(control.respAngle/180*pi)*ecc*2), ...
%                 screen.y_mid - round(sin(control.respAngle/180*pi)*ecc*2), ...
%                 screen.window);
%             ShowCursor;
%         else
%             % changing the angle of the next loop according to the cursor position
%             control.respAngle = atan2(rdkControl.randCenterY-control.mouse_y, control.mouse_x-rdkControl.randCenterX)/pi*180;
%         end
%         %%%%%%%%%% show an arrow for the response... currently just show at
%         % calculate line coordinates and width in pixel
%         [lineWidth, ] = round(dva2pxl(const.line.width, const.line.width, screen));
%         [lineX, lineY] = dva2pxl(cos(control.respAngle/180*pi)*const.line.length/2, sin(control.respAngle/180*pi)*const.line.length/2, screen);
%         lineXY = round([-lineX, lineX; lineY, -lineY]); % this is the main line
%         
%         % now calculate coordinates for the two stroke arrow--make line end
%         % as the center (0, 0)
%         [arrow1X, arrow1Y] = dva2pxl(cos((180-control.respAngle-const.arrowAngle)/180*pi)*const.arrowLength, ...
%             sin((180-control.respAngle-const.arrowAngle)/180*pi)*const.arrowLength, screen);
%         arrow1XY = round([0, arrow1X; 0, arrow1Y]); 
%         [arrow2X, arrow2Y] = dva2pxl(cos((180-control.respAngle+const.arrowAngle)/180*pi)*const.arrowLength, ...
%             sin((180-control.respAngle+const.arrowAngle)/180*pi)*const.arrowLength, screen);
%         arrow2XY = round([0, arrow2X; 0, arrow2Y]); 
%         
%         % draw response line
%         % centered at the screen center
%         Screen('DrawLines', screen.window, lineXY, lineWidth, const.line.colour, [screen.x_mid, screen.y_mid]);
%         Screen('DrawLines', screen.window, [arrow1XY arrow2XY], lineWidth, const.line.colour, [lineX+screen.x_mid, screen.y_mid-lineY]);
        
        % centered at the trajectory center
%         Screen('DrawLines', screen.window, lineXY, lineWidth, const.line.colour, [rdkControl.randCenterX, rdkControl.randCenterY]);
%         Screen('DrawLines', screen.window, [arrow1XY arrow2XY], lineWidth, const.line.colour, [lineX+rdkControl.randCenterX, rdkControl.randCenterY-lineY]);
        %%%%%%%%%%
        
        %%%%%%%%%% this is without the arrow, just a line
%         % make the range within [-90, 90]
%         if control.respAngle>90
%             control.respAngle = control.respAngle-180;
%         elseif control.respAngle<-90
%             control.respAngle = control.respAngle+180;
%         end
%         
%         % calculate line coordinates and width in pixel
%         [lineWidth, ] = round(dva2pxl(const.line.width, const.line.width, screen));
%         [lineX, lineY] = dva2pxl(cos(control.respAngle/180*pi)*const.line.length/2, sin(control.respAngle/180*pi)*const.line.length/2, screen);
%         lineXY = round([-lineX, lineX; lineY, -lineY]);
%         
%         % draw response line
%         Screen('DrawLines', screen.window, lineXY, lineWidth, const.line.colour, [rdkControl.randCenterX, rdkControl.randCenterY]);
        %%%%%%%%%%
        
%         % draw square for photodiode:
%         if photo.mode                                                       % Photodiode Event 2: Fixation target 2 on
%             PTBdraw_photodiodeStimulus(screen, const.photoStimSizePX2, screen.black);
%         end
        
%% (2.4)---------------------------------------------------------------
    case 4                                                                  % PHASE 4: FINISH - WAIT BETWEEN TRIALS
        control.frameRDK = control.frameRDK - 1; % need to keep updating so that the rdk off time is not renewed at every frame!
        tElapse = timePassed - trialData.tRDKoff(curTrial, 1);
        % draw square for photodiode:
        if photo.mode                                                               % Photodiode Event 5: Trial End
            PTBdraw_photodiodeStimulus(screen, const.photoStimSizePX2, screen.black);
        end
        
        if tElapse > const.ITI                         % wait and break
            control.break = 1;
        end
end

