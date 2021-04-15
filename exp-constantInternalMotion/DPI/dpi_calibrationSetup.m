function dpi_calibrationSetup(screen, const)
% DPI CONFIGURTION SETUP Shows crosshair of calibration dots
%   used to setup the DPI
%   PK, 22/10/2019


%Fill background
Screen('FillRect',screen.window,screen.background);
PTBdraw_photodiodeStimulus(screen, const.photoStimSizePX2, screen.black);
% % the following needs to be deleted:
% Screen('FillRect',screen.window,[118.1856  118.1856  118.1856]);
% Screen('Flip',screen.window);
% PTBwait_anykey_press
% Screen('FillRect',screen.window,[0  154.1445  0]);
% Screen('Flip',screen.window);
% PTBwait_anykey_press

for i = 1:size(const.calibPositions,2)
    PTBdraw_bulleye(screen, const.calibPositions{i}, const.calibtationOutRadius, const.calibtationInRadius, screen.calibWhite, screen.background);
end
Screen('DrawText', screen.window, 'Press [ANY KEY] to start calibration.', 1,1,screen.calibWhite);
Screen('Flip',screen.window);
pause(.1);

PTBwait_anykey_press

end

