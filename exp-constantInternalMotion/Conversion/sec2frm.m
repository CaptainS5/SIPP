function frameN = sec2frm(sec, screen)
% translate time duration in seconds to in frame numbers
% probably needs to round/ceil/floor after this computation

frameN = sec/screen.refreshRate; % refresh rate in Hz

end

