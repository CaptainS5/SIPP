function PTBdraw_target_RDK(screen, const, rdkPosition, apertureTexture, textureCenter, textureWindow)

[rdkDotRadiusX, ]= dva2pxl(const.rdk.dotRadius, const.rdk.dotRadius, screen);
rdkDotDiameterX = rdkDotRadiusX*2;

Screen('DrawDots', screen.window, transpose(rdkPosition),...
    rdkDotDiameterX, const.rdk.colour, textureCenter, 1);  % change 1 to 0 to draw square dots
Screen('DrawTexture', screen.window, apertureTexture, [], textureWindow);

end