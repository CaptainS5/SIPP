function PTBdraw_target_RDK(screen, const, rdkPosition, apertureCenter, apertureWindow)

[rdkRadiusX, ]= dva2pxl(const.rdk.dotRadius, const.rdk.dotRadius, screen);
rdkDiameterX = rdkRadiusX*2;

Screen('DrawDots', screen.window, transpose(rdkPosition),...
    rdkDiameterX, const.rdk.colour, apertureCenter, 1);  % change 1 to 0 to draw square dots
Screen('DrawTexture', screen.window, const.rdk.aperture, [], apertureWindow);

end