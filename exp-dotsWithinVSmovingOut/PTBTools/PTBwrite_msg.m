function PTBwrite_msg(el, msg, coordX, coordY, msgColor)
% coords are the disctance to the screen center
    coordX_ = coordX;
    coordY_ = coordY;
    
    if isnumeric(coordX) || isnumeric(coordY)
        
        if ~isnumeric(coordX)
            coordX_ = 0;
        end
        
        if ~isnumeric(coordY)
            coordY_ = 0;
        end
        
        screenCoord = PTBcenter_to_screen([coordX_, coordY_],el);
        
        coordX_ = screenCoord(1);
        coordY_ = screenCoord(2);
        
        if ~isnumeric(coordX)
            coordX_ = coordX;
            coordY_ = screenCoord(2);
        end
        
        if ~isnumeric(coordY)
            coordX_ = screenCoord(1);
            coordY_ = coordY;
        end
    end

    DrawFormattedText(el.window, msg, coordX_, coordY_, msgColor);
end