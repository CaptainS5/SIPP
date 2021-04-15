function dpi_calibrationProcedure(screen,const)
% DPI CALIBRATION PROCEDURE Shows crosshair of calibration dots
%   used to calibrate the DPI. Press Enter for each fixation (ensures that
%   no fixation spot was missed)
%   Procedure: Center, Left, Center, Right, Center, Top, Center, Bottom,
%   Center (each direction has 2 locations - 13 in total).
%   PK, 22/10/2019
global DPI

for i = 1:size(const.calibPositions,2)
    %Fill background
    Screen('FillRect',screen.window,screen.background);
    
    %Bulleye:
    PTBdraw_bulleye(screen, const.calibPositions{i}, const.calibtationOutRadius, const.calibtationInRadius, screen.calibWhite, screen.background);
    
    %Photodiode
    switch i
        case {1, 3, 5, 7, 9, 11, 13}
            PTBdraw_photodiodeStimulus(screen, const.photoStimSizePX2, screen.white);
        otherwise
            PTBdraw_photodiodeStimulus(screen, const.photoStimSizePX2, screen.black);
    end
    
    Screen('Flip',screen.window);
    % Key press (after proper fixation), to continue to next location)
    PTBwait_anykey_press
end

end

