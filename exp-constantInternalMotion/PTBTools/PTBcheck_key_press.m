function [bPressed keyPressed]= PTBcheck_key_press(keyname)
% if key named <keyname> is pressed, return 1
% keyname could be a vector of all possible keys that you want to check
    [~, ~, keyCode] = KbCheck;
    
    if any(keyCode(keyname))
        KbReleaseWait;
        idx = find(keyCode(keyname));
        if length(idx)==1 % only one key pressed
            disp(['EXP: ' KbName(keyname(idx)) ' pressed']);
            bPressed = 1;
            keyPressed = keyname(idx);
        else
            bPressed = 0;
            keyPressed = NaN;
        end
    else
        bPressed = 0;
        keyPressed = NaN;
    end
end
