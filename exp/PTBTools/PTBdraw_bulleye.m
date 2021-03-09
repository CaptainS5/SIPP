function PTBdraw_bulleye(screen, center, radiusOUT, radiusIN, colorOUT, colorIN)
% DRAW CALIBRATION BULL EYE draws to ovals. Bull eye used for DPI setup and
%   calibration 
%   PK, 22/10/2019
%   INPUTS: screen    - structure containing screen info
%         positions - all position(s) where Bull Eye should be drawn
%         radius    - radius of the outer/inner bulleye - needs 2 values
%         color     - color of outer/inner bulleye circle - needs 2 colors

%% OUTER OVAL
radiusOUT_PX = round(tan(radiusOUT/180*pi) * screen.dist * screen.ppcX);
left         = center(1) - radiusOUT_PX;
right        = center(1) + radiusOUT_PX;
top          = center(2) - (radiusOUT_PX*screen.pixelRatioWidthPerHeight);
bottom       = center(2) + (radiusOUT_PX*screen.pixelRatioWidthPerHeight);
position_out = [left top right bottom];        
Screen('FillOval', screen.window, colorOUT, position_out);


%% INNER OVAL
radiusIN_PX  = tan(radiusIN) * screen.dist * screen.ppcX;
left         = center(1) - radiusIN_PX;
right        = center(1) + radiusIN_PX;
top          = center(2) - (radiusIN_PX*screen.pixelRatioWidthPerHeight);
bottom       = center(2) + (radiusIN_PX*screen.pixelRatioWidthPerHeight);
position_in  = [left top right bottom];        
Screen('FillOval', screen.window, colorIN, position_in);

end

