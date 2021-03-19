function expMain(const, screen, keys, eyelink, dpi_set, photo, trialCondition, trialData, sbj)
% =========================================================================
% expMain(const, screen, eyelink, trialData, sbj)
% =========================================================================
% Function created to run the framework of experiment:
% 1. Initialize Experiment/ Show Experiment Instructions
% 2. Setup & Calibrate DPI; Enter Block Loop:
% 3. Single Trial Setup
% 4. Run Single Trials (Enter Trial-While-loop, main flip)
% 5. Finish Trial (close files)
% 6. Close Experiment
% + DPI - a global structure initiated during dpi_init; contains analog
%   input data collection
%
% Last Changes:
% 04/11/2019 (PK) - cleaned up function
% 04/Mar/2021 (XW) - save all data in the main folder instead of in block
% folders
% -------------------------------------------------------------------------
% Inputs:
% const:     structure containing different constant settings
% screen:    strucrure containing screen settings
% keys:      mapping of the response keys
% eyelink:   structure containing eyelink settings
% dpi_set:   structure containing DPI settings
% photo:     structure containting photodiode settings
% trialCondition: a lookup table containing all trial relevant information
% trialData: recorded trial data during the experiment
% sbj:       structure containing subject information
% -------------------------------------------------------------------------
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) INITIALIZE EXPERIMENT:
%     1) initialize experiment, 2) show experiment instruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (1.1) Initialize Psychtoolbox; create el-structure, if Eyelink not used
if const.startExp; ListenChar(2); end          %HideCursor;                  % Disable key output to Matlab window and hide cursor:
global DPI

%% (1.2) Initialize and Calibrate Eyelink
% if eyelink.mode
%     [el]                     = eyelink_init(screen, const, eyelink);
% end


%% (1.3) Show General Experiment Instructions
% disp('EXP: Begin experiment: show instructions');
% PTBinstruction_page(1,screen)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) SETUP & CALIBRATE DPI; ENTER BLOCK LOOP:
%     1) show Calibration Instruction, 2) show DPI Setup screen,
%     3) run calibration, 4) setup trial counter for trial-loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    for block = sbj.block:size(const.numTrialsPerBlock,2)                  % begin block FOR-LOOP
        
        %% (2.1) Show Calibration Instructions
        PTBinstruction_page(2,screen,const);
        
        %% (2.2) Initialize and (2.3) Calibrate DPI
        if dpi_set.mode
            dpi_calibrationSetup(screen, const);                                % show Setup-Screen (Crosshair of bulleys)
            
            start(DPI)                                                          % start recording of analog inputs
            dpi_calibrationProcedure(screen,const);                             % run actual DPI calibration
            pause(0.1)
            stop(DPI)                                                           % stop recording in order to collect calibration data
            pause(0.1)
            filename = sprintf('/calibrationPre_%02d.mat',block);
            dpi_saveData(filename,blockFolder)                                  % save Pre-calibration data
        end
        
        %% or, initialize Eyelink and set up
        if eyelink.mode
            if ~EyelinkInit(eyelink.dummy)
                errorEyelink('Problems with Eyelink connection!');
                return;
            end
            
            eyelink.el=EyelinkInitDefaults(screen.window);
            eyelink.el.backgroundcolour = screen.background;
            EyelinkUpdateDefaults(eyelink.el);
            
            % check which eye is recorded
            eyelink.eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
            if eyelink.eye_used == -1
                eyelink.eye_used = eyelink.el.RIGHT_EYE;
            end
            
            % make sure that we get gaze data from the Eyelink
            Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
            Eyelink('Command', 'link_update_interval = 50'); % extends the keep alive
            % time for the network connection to make sure the connection stays active
            % even during long pauses, nly need to set it once at the start of the experiment
            
            
            % Calibrate the eye tracker
            EyelinkDoTrackerSetup(eyelink.el);
            
            % do a final check of calibration using driftcorrection
            EyelinkDoDriftCorrection(eyelink.el);
            
            % prepare eye recording
            eyelink.edfName = [sbj.name '_b', num2str(block, '%02d'), '.edf'];
            if (size(eyelink.edfName, 2)-4>8)
                errorEyelink('edf filename is too long!'); % Security loop against Eyelink
                % Un-registration of data if namefile
            end
            % open file to record data to
            cd(sbj.sbjFolder)
            Eyelink('Openfile', eyelink.edfName);
        end
        
        %% (2.4) Setup trial counter
        % just use whatever you input for the first time
        if block==sbj.block % assign if just started running; newBlock will be renewed whenever finished a block
            if sbj.trial==1 % if starting a new block, initialize trialData for the current block
                newBlock = 1;
                % delete potential earlier trials in the block, if
                % necessary...
                if size(trialData, 1)>0
                    trialData(trialData.blockN==block, :) = [];
                end
            else
                newBlock = 0;
            end
        end
        
        if newBlock==1
            trialDataNew = table;
            trialDataNew = trialCondition(trialCondition.blockN==block, :); % all conditions for the current block
            % add the temporal measures to be recorded
            trialDataNew.tMainSync(:, 1)          = 0;                   % Time (GetSecs) at trial start, also fixation on.
            trialDataNew.tRDKon(:, 1)        = NaN;
            trialDataNew.tRDKoff(:, 1)  = NaN;                     % actual measured time that the target appeared/disappeared
            %             trialDataNew.tResponse(:, 1)          = NaN;               % also the end of the trial
            trialDataNew.t_start_VBL(:, 1:5)        = NaN;
            trialDataNew.t_rdkOn_VBL(:, 1:5)         = NaN;
            trialDataNew.t_rdkOff_VBL(:, 1:5)        = NaN;
            %             trialDataNew.t_response_VBL(:, 1:5)        = NaN;
            trialDataNew.iterations(:, 1)        = NaN;
            %             trialDataNew.choice(:, 1)          = NaN; % -1=up, 1=down
            %             trialDataNew.choiceCorrect(:, 1)          = NaN; % 0=wrong, 1=correct
            trialDataNew.repeat(:, 1)          = 0; % if the trial was repeated, mark repeat as 1
            trialDataNew.trialCounter           = [1:size(trialDataNew, 1)]';
            %             trialDataNew.rdkFrameResponse(:, 1) = NaN;
            trialData = [trialData; trialDataNew]; % fill in trialData with the trial conditions of the current block
        end
        if block==sbj.block
            currentTrial = find(trialData.blockN==block & trialData.trialCounter==sbj.trial);
        else
            currentTrial = find(trialData.blockN==block & trialData.trialCounter==1);
        end
        
        %% (2.5) Show Block Instrustion (before block starts):
        Screen('FillRect', screen.window, screen.background);                % Draw background
        if photo.mode
            PTBdraw_photodiodeStimulus(screen, const.photoStimSizePX2, screen.black);
        end
        
        PTBwrite_msg(screen, ['Block ', num2str(block)], 'center', 15, screen.black)
        PTBwrite_msg(screen, 'press [space] to continue', 'center', -15, screen.black)
        %         Screen('Flip', screen.window,screen.refreshRate);
        Screen('Flip', screen.window, screen.refreshRate, 0, 0, 2);
          keyPressed = 0;
        while keyPressed==0
            keyPressed = PTBcheck_key_press(keys.space);
        end
        %             PTBwait_anykey_press
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (3) SINGLE TRIAL SETUP
        %      1) enter TRIAL-WHILE-loop: trial setup (control-struct),
        %      2) open target file, 3) Show blank screen, wait for participant
        %      4) mark synctime for target display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% (3.1) Trial Setup:
        while currentTrial <= size(trialData, 1)                                % START TRIAL WHILE-LOOP
            
            % initialize trial control structure
            control                 = [];                                       % Setup a Trial Control Structure
            control.forceRecalibEL  = 0;
            control.currentTrial    = currentTrial;                                 % current trial
            control.trialName       = sprintf('%.3d', currentTrial);
            
            control.rdkCoh          = trialData.rdkCoh(currentTrial);
            control.rdkApertureDir  = trialData.rdkApertureDir(currentTrial);
            control.rdkInternalSpeed = trialData.rdkInternalSpeed(currentTrial);
            control.rdkInternalDir  = trialData.rdkInternalDir(currentTrial);
            control.fixationFrames = ceil(sec2frm(trialData.fixationDuration(currentTrial), screen));
            
            control.mode            = 1;                                        % different phases of the trial (changes values in runSingleTrials)
            control.break           = 0;                                        % trial completion
            control.repeat          = 0;
            control.abort           = 0;                                        % trial abortion (ESCAPE key press)
            control.frameFix        = 0;                                        % frame counter, mark current event
            control.frameRDK        = 0;
            iterations              = 0;                                        % counts iterations of while-loop
            control.repeat          = 0;                                        % if the trial should be repeated at the end of the block
            
            % For Trial Videos:
            if const.makeVideo == 1                                             % settings for single-trial video:
                trialData.tMainSync(currentTrial, 1)        = 0;                       % can change params to show here.
                control.frameCounter                 = 1;
            end
            
            %% (3.2) Open Target-File (record displayed target path):
            % Open Target File
            control.targetFile = ['rdkseed_t' control.trialName];
            %         control.targetFID       = fopen(control.targetFile, 'w');
            [rdkControl seed] = generateTrialRDKpositions(const, screen, control); % generate the position of dots in each frame in the whole trial
            % generate aperture for rdk
%             if const.apertureType==1 % aperture translates across the dot field
%                 for frameN = 1:size(rdkControl.apertureCenterPos, 2)
%                     rdkControl.apertureTexture{frameN} = PTBmakeAperture(const, screen, rdkControl.apertureCenterPos{frameN});
%                 end
%             else % dots move together with the aperture
                rdkControl.apertureTexture = PTBmakeAperture(const, screen, 0);
%             end
            
            fprintf('EXP: begin Block %d Trial %d \n', block, trialData.trialCounter(currentTrial, 1));
            
            
            %% (3.3) Start DPI recording for current trial:
            if dpi_set.mode
                start(DPI)                                                      % start DPI data collection for current trial
                %             pause(0.1)
            end
            %% or, start Eyelink recording
            if eyelink.mode
                Eyelink('Message', 'TrialID: %s', num2str(currentTrial));
                WaitSecs(0.05);
                Eyelink('Command', 'set_idle_mode'); %it puts the tracker into offline mode
                WaitSecs(0.05); % it waits for 50ms before calling the startRecording function
                Eyelink('StartRecording');
                Eyelink('Message', 'SYNCTIME');
            end
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (4) RUN SINGLE TRIALS (WHILE-LOOP):
            %     1) start DPI recording (for current trial),
            %     2) update trial timer, 3) run runSingleTrials.m,
            %     4) main flip,  5) save target data, 6) end the trial/stop DPI
            %     recording
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('EXP: Enter while loop \n');
            
            %% (4.1) Initialize Eyelink Data:
            %             eyelink.Data.coord = [0, 0];                                            % initialize eyelink data (since new sample can be missing for first couple of loops)
            % %             if eyelink.mode == 1
            % %                 control.data_time = Eyelink('TrackerTime')...                       % get time before entering while loop (this should be close to 0)
            % %                     - trialData.tMainSync(currentTrial);
            % %             else
            % %                 control.data_time = GetSecs...                       % get time before entering while loop (this should be close to 0)
            % %                     - trialData.tMainSync(currentTrial);
            % %             end
            
            %% ================================================================
            while  1                                                            % begin WITHIN-TRIAL-WHILE-loop (code will run through this part until break
                iterations = iterations + 1;                                    % this runs every refresh Rate cycle (should be free of heavy computations)

                %% (4.2) Check Eyelink Recording:
                if eyelink.mode == 1
                    errorEyelink = Eyelink('CheckRecording');
                    if(errorEyelink ~= 0)
                        fprintf('EXP: CheckRecording error %d\n', errorEyelink);
                        break;
                    end
                end
                
                %% (4.3) Eyelink Data Acquisition
%                 if eyelink.mode == 1
%                     [eyelink] = eyelink_dataAcquisition(eyelink.el,eyelink,trialData,currentTrial);
%                     %                     control.data_time = eyelink.Data.time;
%                     %                 else                                                                % for dummy mode
%                     %                     control.data_time = GetSecs - trialData.tMainSync(currentTrial);
%                 end
                
                %% Update trial timer
                if iterations>1
                    control.data_time = GetSecs - trialData.tMainSync(currentTrial, 1);    % trial timer: how much time has passed since trial started
                else
                    control.data_time = 0;
                end
                
                %% (4.3) Run Single Trial:
                % =============================================================
                [const, control] = drawSingleFrameEyelink(const, trialData, control, screen, photo, rdkControl, eyelink);
                % =============================================================

                % show gaze or finger if it is turned on
                %         if eyelink.mode && const.showGaze
                %             PTBdraw_circles(screen, eyelink.Data.coord, 10, [255 255 255]);
                %         end
                
                %% (4.4) Main Flip and Time Stemps:
                % =============================================================
                [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', screen.window);
                % =============================================================
                %                 % for debugging:
                %                 if  ~const.startExp && control.frameElapse == control.tFlash
                %                     control.PursuitDir
                %                     control.FlashDir
                %                     control.tFlashPosition
                %                     control.eyeTarget
                %                 end
                
                if iterations == 1 % mark start of the trial
                    trialData.t_start_VBL(currentTrial,:) = [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos];
                    trialData.tMainSync(currentTrial, 1)            = GetSecs;                 % Get the current time as time when the trial started (used for the trial timer)
                    if eyelink.mode
                        Eyelink('Message', 'fixationOn');
                    end
                elseif control.frameRDK >= 1
                    if control.frameRDK==1 % RDK onset
                        trialData.tRDKon(currentTrial, 1) = control.data_time;
                        trialData.t_rdkOn_VBL(currentTrial,:)  = [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos];
                        if eyelink.mode
                            Eyelink('Message', 'rdkOn');
                            Eyelink('Message', ['frameRDK ', num2str(control.frameRDK)]);
                        end
                    else % just record frame number
                        if eyelink.mode
                            Eyelink('Message', ['frameRDK ', num2str(control.frameRDK)]);
                        end
                    end
                elseif control.frameRDK == -1 % RDK offset
                    trialData.tRDKoff(currentTrial, 1) = control.data_time;
                    trialData.t_rdkOff_VBL(currentTrial,:) = [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos];
                    if eyelink.mode
                        Eyelink('Message', 'rdkOff');
                        Eyelink('Message', ['frameRDK ', num2str(control.frameRDK)]);
                    end
                end
                
                if control.repeat==1
                    trialData.repeat(currentTrial, 1) = 1;
                    % add the current trial to the end of the block...
                    trialData = [trialData; trialData(currentTrial, :)];
                    % initialize trial information for the added trial
                    trialData.tMainSync(end, 1)          = 0;                   % Time (GetSecs) at trial start, also fixation on.
                    trialData.tRDKon(end, 1)        = NaN;
                    trialData.tRDKoff(end, 1)  = NaN;                     % actual measured time that the target appeared/disappeared
                    %                     trialData.tResponse(end, 1)          = NaN;               % also the end of the trial
                    trialData.t_start_VBL(end, 1:5)        = NaN;
                    trialData.t_rdkOn_VBL(end, 1:5)         = NaN;
                    trialData.t_rdkOff_VBL(end, 1:5)        = NaN;
                    %                     trialData.t_response_VBL(end, 1:5)        = NaN;
                    trialData.iterations(end, 1)        = NaN;
                    %                     trialData.choice(end, 1)          = NaN; % -1=up, 1=down
                    %                     trialData.choiceCorrect(end, 1)          = NaN; % 0=wrong, 1=correct
                    trialData.repeat(end, 1)          = 0; % if the trial was repeated, mark repeat as 1
                    trialData.trialCounter(end, 1) = trialData.trialCounter(end-1, 1)+1;
                    %                     trialData.rdkFrameResponse(end, 1) = NaN;
                    
                    WaitSecs(1) % show the information on the screen
                    break
                end
                
                % Trial Video:
                if const.makeVideo == 1                                         % if making video, here images are taken of each frame
                    imageArray(:,:,:,control.frameCounter) = Screen('GetImage',screen.window);
                    control.frameCounter = control.frameCounter + 1;            % update frameCounter
                end
                
                %% (4.5) End the trial (i.e. break while-loop; stop DPI recording)
                if control.break
                    if dpi_set.mode
                        pause(0.1)
                        stop(DPI)                                               % stop recording in order to collect calibration data
                        pause(0.1)
                        filename = sprintf('/trial_%03d_tFlash_%02d_FlashDir_%02d_PursuitDir_%02d.mat',...
                            currentTrial,control.tFlash, control.FlashDir, control.PursuitDir);
                        dpi_saveData(filename,blockFolder)                      % save eye trial data
                    end
                    trialData.iterations(currentTrial, 1) = iterations;
                    fprintf('EXP: Block %d Trial %d finished \n', block, trialData.trialCounter(currentTrial, 1));
                    break; % BREAK the while loop, if trial was finished!
                end
                
                % check key press
                [bPressed keyPressed] = PTBcheck_key_press([keys.escape, keys.recalibration]);
                if keyPressed==keys.escape
                    fprintf('EXP: Experiment aborted by pressing esc key \n');
                    control.abort = 1;
                    break;                                                      % BREAK the while loop, if trial was aborted by ESC press
                elseif keyPressed==keys.recalibration
                    control.forceRecalibEL = 1;
                end
                
                % check forced Recalibration:
                if control.forceRecalibEL
                    break;
                end
                
            end                                                                 % end of WITHIN-TRIAL-WHILE-loop
            %% =================================================================
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (5) FINISH TRIAL:
            %     1) Save and move target file/ save trialData
            %     2) Trial Abort, 3) create video, 4) update trial counter,
            %     5) run post calibration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if eyelink.mode
                Eyelink('Message', 'TRIALEND');
                Eyelink('Command','clear_screen 0'); % clears the box from the Eyelink-operator screen
                Eyelink('Command', 'set_idle_mode');
                WaitSecs(0.05);
                Eyelink('StopRecording');
            end
            
            %% (5.1) if recalibration is forced (repeats trial) or trial finished properly:
            if control.forceRecalibEL                                               % check if Eyelink calibration was forced
                Eyelink('Message', 'FORCE_RECALIB_EL');                             % send Eyelink Message
                eyelink_recalibration(control,const,eyelink.el);
            end         
            
            %% (5.2) Save and move target file/ save trialData
            % close TARGET file and save/move file
            %         for target = 1:size(targetInfo,1)
            %            fwrite(control.targetFID, targetInfo{target,:});
            %         end
            rdkAperturePos = rdkControl.apertureCenterPos;
            rdkDirFrames = rdkControl.dotDir; % position matrices take too much space
            
            save([sbj.rdkFolder, '/', control.targetFile, '_', sbj.date, '.mat'], 'seed', 'rdkDirFrames', 'rdkAperturePos')
            fprintf('EXP: Target data is saved in %s %s\n\n', sbj.rdkFolder, control.targetFile);
            clear rdkControl seed
            
            % save updated trialData structure:
            save([sbj.sbjFolder, '/trialData.mat'], 'trialData');                    % save data in subjectFolder
            
            %% (5.3) Experiment aborted
            if control.abort == 1
                Eyelink('CloseFile');
                %                 fclose(control.targetFID);
                %                 movefile(control.targetFile, [sbj.sbjFolder]);
                %         fprintf('EXP: Target data is saved in .\\data\\%s\\%s\n', sbj.filename, control.targetFile);
                throw(MException('EXP: MainLoop','Experiment aborted'));
            end
            
            
            %% (5.4) Create Video, if chosen:
            if const.makeVideo == 1                                             % put the single images into a video clip:
                if ~exist('./data/video', 'dir'); mkdir('./data/video'); end
                m = cat(4,imageArray);
                vidName = input(sprintf('\n\tVideo name:\t'),'s');
                fprintf('\n\tProcessing video - please wait...\n')
                writerObj = VideoWriter(sprintf('./data/video/%s',vidName),'MPEG-4');
                writerObj.FrameRate = 120;
                writerObj.Quality = 100;
                open(writerObj);
                writeVideo(writerObj,m);
                close(writerObj);
            end
            
            %% (5.5) Update trial counter
            currentTrial = currentTrial + 1;                                            % Update trial counter
        end                                                                     % end of TRIAL-WHILE-loop
        %% (5.5) Run Post Calibration: for DPI only
        if dpi_set.mode
            start(DPI)                                                          % start recording of analog inputs
            dpi_calibrationProcedure(screen,const);                             % run actual DPI calibration
            pause(0.1)
            stop(DPI)                                                           % stop recording in order to collect calibration data
            pause(0.1)
            filename = sprintf(['/calibrationPost_b%02d_' sbj.date '.mat'],block);
            dpi_saveData(filename,blockFolder)                                  % save Pre-calibration data
        end
        
        %% (6.1) Save complete files
        if eyelink.mode
            % eye recording output
            %             Eyelink('Command','clear_screen 0'); % clears the box from the Eyelink-operator screen
            %             Eyelink('Command', 'set_idle_mode');
            WaitSecs(0.5);
            Eyelink('CloseFile');
            try
                fprintf('Receiving data file ''%s''\n', eyelink.edfName);
                status=Eyelink('ReceiveFile');
                if status > 0
                    fprintf('ReceiveFile status %d\n', status);
                end
                if 2==exist(eyelink.edfName, 'file')
                    fprintf('Data file ''%s'' can be found in ''%s''\n', eyelink.edfName, sbj.sbjFolder);
                end
            catch
                fprintf('Problem receiving data file ''%s''\n', eyelink.edfName, sbj.sbjFolder);
            end
        end
        Experiment.const       = const;
        Experiment.dpi_set     = dpi_set;
        Experiment.eyelink     = eyelink;
        Experiment.trialData   = trialData;
        Experiment.sbj         = sbj;
        Experiment.screen      = screen;
        save([sbj.sbjFolder '/info_Experiment.mat'], 'Experiment');
        
        % initialize for the next block, if running right after this one
        newBlock = 1;
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (6) CLOSE EXPERIMENT:
    %     1) save trialData and other structures, 2) cleanup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% (6.1) Save complete trialData file
    % already did within each block
    
    %% (6.2) Download Eyelink Files:
    % if eyelink.mode == 1
    %     for i = 1:numel(const.numTrialsPerBlock)
    %         try
    %             fprintf('EXP: Receiving data file ''%s''\n', eyelink.edfFile{i});
    %             status=Eyelink('ReceiveFile', eyelink.edfFile{i}, [sbj.sbjFolder], 1);
    %             if status > 0
    %                 fprintf('EXP: ReceiveFile status %d\n', status);
    %             end
    %         catch rdf
    %             fprintf('EXP: Problem receiving data file ''%s''\n', eyelink.edfFile{i});
    %             rdf;
    %         end
    %     end
    % end
    
    %% (6.3) End Experiment
    disp('EXP: Finish experiment');
    cleanup(eyelink);
    
    
    
    
    
    %% MY ERR
catch myerr                                                                 % this "catch" section executes in case of an error in the "try" section
    if eyelink.mode && (trialData.trialCounter(currentTrial, 1)>1)
        % eye recording output
        %             Eyelink('Command','clear_screen 0'); % clears the box from the Eyelink-operator screen
        %             Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.5);
        Eyelink('CloseFile');
        try
            fprintf('Receiving data file ''%s''\n', eyelink.edfName);
            status=Eyelink('ReceiveFile');
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(eyelink.edfName, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', eyelink.edfName, sbj.sbjFolder);
            end
        catch
            fprintf('Problem receiving data file ''%s''\n', eyelink.edfName, sbj.sbjFolder);
        end
        % rename the edf file so it doesn't replace the "normal" one
        oldName = [sbj.sbjFolder '\' eyelink.edfName];
        currentTime = clock;
        currentDate = sprintf('%d-%d-%d_%d%d', currentTime(1:5));
        newName = [sbj.sbjFolder '\' eyelink.edfName(1:end-4) '_' currentDate '.edf'];
        copyfile(oldName, newName)
    end
    
    save([sbj.sbjFolder, '/trialData.mat'], 'trialData');                    % save data in subjectFolder
    
    Experiment.const       = const;
    Experiment.dpi_set     = dpi_set;
    Experiment.eyelink     = eyelink;
    Experiment.trialData   = trialData;
    Experiment.sbj         = sbj;
    Experiment.screen      = screen;
    save([sbj.sbjFolder '/info_Experiment.mat'], 'Experiment');
    
    % above.  Importantly, it closes the onscreen window if its open.
    fclose('all');                                                          % close any open files
    delete(DPI);
    cleanup(eyelink);
    commandwindow;
    myerr;
    myerr.message
    myerr.stack.line
end                                                                         % end of try...catch
end                                                                         % end of function