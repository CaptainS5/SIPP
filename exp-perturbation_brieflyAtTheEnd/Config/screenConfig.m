function [screen] = screenConfig(const)
% ===============================================================
% screenConfig(const)
% ===============================================================
% written by: Philipp KREYENMEIER (philipp.kreyenmeier@gmail.com)
% 
% Last Changes:
% PK, 22/10/2019: Changed Deg2Pix and Pix2Deg conversion
%     (adapted from M.S. Plaidtracking code)
% 
% ---------------------------------------------------------------
% Inputs: 
% const: structure containing different constant settings
% ---------------------------------------------------------------
% 
% Overview setups in the lab:
% DebugScreen:  w = 47.3, h = 29.6, wp = 1680, hp = 1050, d = 50
% EyeStrike:    w = 43.5, h = 34.7, wp = 1280, hp = 1024, d = 44
% Backroom :    w = 39.6, h = 29.7, wp = 1600, hp = 1200, d = XX
% DPI_Koerner:  w = 36.7, h = 27,   wp = 1280, hp = 1024, d = 85; Hertz = 85
% Eyelink_Koerner:  w = 40.7, h = 30.3,   wp = 1600, hp = 1200, d = 55; Hertz = 85
% ---------------------------------------------------------------


%% (1) Enter desired Screen Settings - this specific to each setup
screen.name                = 'Eyelink_ICICSbackroom';
% Screen('Preference', 'SkipSyncTests', 2);
screen.desired_Hertz       = 85; %											% desired Refresh Rate (Hz) of the Screen
screen.desired_widthPX     = 1600;  %										% desired width (pixels) of the Screen
screen.desired_heightPX    = 1200;  %		 								% desired height(pixels) of the Screen

screen.widthCM             = 39.6; % backroom in ICICS                              % measured size (in cm)
screen.heightCM            = 29.7; % backroom in ICICS                              % measured size (in cm)
screen.widthMM             = screen.widthCM*10;                                     % measured size (in mm) --> this is needed for screen calibration
screen.heightMM            = screen.heightCM*10;                                    % measured size (in mm)
screen.dist                = 55; %                                                  % measured distance Subject-Screen (in cm)

% Get Screen from PTB:
% if ~const.startExp && ~const.runScreenCalib
% %     PsychDebugWindowConfiguration(0,0.5); 									% when in debugging, set transparency to 50%
%     %     screen.number = min(Screen('Screens'));
%     screen.number = max(Screen('Screens'));
%     %     sca
% else
%     screen.number = max(Screen('Screens'));
%     %     screen.number = 1;                                                      % needs to be 1 for DPI Setup with extended screens
% end
%  screen.number = max(Screen('Screens')); % for Koerner
screen.number = 1; % for X717, backroom in ICICS

%% (2) Open Screen - initiating PTB:
% Pixel Size:
screen.clr_depth           = Screen('PixelSize', screen.number);

if ~const.runScreenCalib                                                    % if we don't calibrate the screen, open PTB here
    AssertOpenGL;                                                           % We use PTB-3

    disp(['EXP: Open a graphics window on the main screen ' , ...
     'using the PsychToolbox Screen function.']);

    [screen.window, screen.wRect] = Screen('OpenWindow', screen.number, [0 0 0], [], screen.clr_depth, 2);
    Screen(screen.window, 'BlendFunction', GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    priorityLevel = MaxPriority(screen.window);
    Priority(priorityLevel);
else                                                                        % if calibrate screen, PTB will be opened later
    screen.window = screen.number;
end


%% (3) Validate Screen Size and Refresh rate:
% Validate Screen Size:
[screen.widthPX, screen.heightPX] = Screen('WindowSize', screen.window);
if (screen.widthPX ~= screen.desired_widthPX || screen.heightPX ~= screen.desired_heightPX(1)) && const.startExp
    error('Incorrect screen resolution => Please restart the program after changing the resolution to [%i,%i]',screen.desired_widthPX(1),screen.desired_heightPX(1));
end

% Validate Refresh Rate:
screen.refreshRate         = 1/(Screen('FrameRate',screen.window));
if screen.refreshRate == inf
    screen.refreshRate     = 1/60;
elseif screen.refreshRate == 0
    screen.refreshRate     = 1/60;
end
screen.hertz = 1/(screen.refreshRate);
if (screen.hertz >= 1.1*screen.desired_Hertz || screen.hertz <= 0.9*screen.desired_Hertz) && const.startExp
    error('Incorrect refresh rate => Please restart the program after changing the refresh rate to %i Hz',screen.desired_Hertz);
end


%% (4) Get Basic Screen parameters:

% Screen Center:
screen.x_mid               = screen.widthPX/2;
screen.y_mid               = screen.heightPX/2;
screen.center                 = [screen.x_mid screen.y_mid];

% Pix2Deg and Deg2Pix conversions:
% pixels per cm, horizontal and vertical
% screen.ppcX = screen.widthPX/screen.widthCM; % width
% screen.ppcY = screen.heightPX/screen.heightCM; % height
% screen.pixelRatioWidthPerHeight = screen.ppcY/screen.ppcX;
% % the relationship between pixel and cm are fixed everywhere

screen.width_deg          = 2 * (180/pi) * atan2(screen.widthCM  / 2, screen.dist);
screen.height_deg         = 2 * (180/pi) * atan2(screen.heightCM  / 2, screen.dist);
screen.x_ppd              = screen.widthPX / screen.width_deg;
screen.y_ppd              = screen.heightPX / screen.height_deg;
screen.pixelRatioWidthPerHeight         = screen.x_ppd / screen.y_ppd;

screen.ppd                = exp(mean(log([screen.x_ppd screen.y_ppd])));
screen.dpp                = 1 / screen.ppd;

%% (5) Colors and Luminance Settings:
screen.background          = 60;                                            % Background luminance 1
screen.background2         = 60;                                            % Background luminance 2
screen.white               = WhiteIndex(screen.window);
screen.black               = BlackIndex(screen.window);
screen.gray                = round((screen.white+screen.black)./2);
screen.orange              = [255,150,0];
screen.calibBlack          = BlackIndex(screen.window);
screen.calibWhite          = WhiteIndex(screen.window);
screen.calibGray           = round((screen.calibWhite+screen.calibBlack)./2);
screen.foregroundcolour    = WhiteIndex(screen.window);
screen.msgfontcolour       = WhiteIndex(screen.window);
screen.msgfontcolour2      = WhiteIndex(screen.window);
screen.imgtitlefontcolour  = WhiteIndex(screen.window);
screen.imgtitlecolour      = BlackIndex(screen.window); 
screen.warningColor        = [255, 0, 0];
screen.gaussTargetColor    = screen.calibGray;

%% (6) Load in Gamma Calibration, if Screen is not being calibrated:
% eyelink in S257 not calibrated yet... will need the photometer?!
% screen.calibFlag           = 1;                                             % 1 = use gamma Calibration (recommanded); 0 = run experiment w/o screen calibration
% screen.calibType           = 1;                                             % 1 = Gray linearized; 2 = RGB linearized
% screen.dirCalib            = 'GammaCalib';
% screen.desiredValue        = 16;                                            % number of test colors/grey scales
% if ~const.runScreenCalib
%     if screen.calibFlag  == 1                                               % Load gamma calibration
%         if screen.calibType == 1                                            % Gray linearized calibration loaded
%             [screen]           = loadGRAYLinCalib(screen);
%         elseif screen.calibType == 2                                        % RGB linearized calibration loaded
%             [screen]           = loadRGBLinCalib(screen);
%         end
%         loadgammaCalib(screen);
%     end
% end


end